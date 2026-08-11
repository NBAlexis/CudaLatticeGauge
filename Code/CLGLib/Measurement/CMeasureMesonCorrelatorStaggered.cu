//=============================================================================
// FILENAME : CMeasureMesonCorrelatorStaggered.cu
//
// DESCRIPTION:
// Staggered meson correlator measurement.
// Supports 20 meson channels with up to 3 SO3-degenerate sub-channels each.
// Computes both wall-to-wall (W2W) and point-to-point (P2P) contractions.
//
// DSource_{B,c_2}(x,c_1) = D^-1 source _{A,c}(x,t) = sum_y D ^ -1_{ x, 2y + B; c1,c2 }
// DdSource_{B,c_2}(x,c_1) = Dd^-1 source _{A,c}(x,t) = sum_y Dd ^ -1_{ x, 2y + B; c1,c2 }
//
// sink_{A,c_1} \cdot DSource =\sum _x pick_even_{A,c_1}(DSource_{B,c_2}, x) = sum_{x,y} D ^ -1_{2x+A,2y+B; c1,c2}
// pick_even_{A,C} means pick 2x + A site and c color component from the source
//
// Some explanation of the names, may be wrong names:
//
// ==========================
// Wall to wall contraction:
// C(t) = sum_{ x1,x2,y1,y2,A,B,c1,c2 } s(A)s(B) D ^ -1_{2x1+A,2y1+B; c1,c2} (Dd ^ -1_{2x2+A+d,2y2+B+d; c1,c2})^*
//
// Keep the color pair explicit:
//   p(t, A, B, d) = sum_{c1,c2} S1(t, A, B, c1, c2) * conj(S2(t, A+d, B+d, c1, c2))
// where each S is the spatial wall sum over the corresponding color pair.
//
// ==========================
// Point to point contraction:
// C(t) = sum_{ x,y1,y2,A,B,c1,c2 } s(A)s(B) D ^ -1_{ 2x + A,2y1 + B; c1,c2 } (Dd ^ -1 _{ 2x+A+d,2y2+B+d; c1,c2 })^*
// So we calculate p(t, A, B, d) as
// p(t, A, B, d) = \sum _{c_1,c_2,x} pick_even_{A,c_1}(DSource_{B,c_2},x) * [pick_even_{A+d,c_1}(DdSource_{B+d,c_2},x)]^*
//
// ==========================
// After p(t, A, B, d) is calculated,
// C_d(t) = \sum_{A,B} sign(A)sign(B) p(t, A, B, d)
// for both cases
//
// REVISION:
//  [09/28/2020 nbale]
//  [06/08/2026 nbale] Added SO3 sub-channels and P2P contraction
//=============================================================================

#include "CLGLib_Private.h"
#include "Data/Field/Staggered/CFieldFermionKST.h"
#include "CMeasureMesonCorrelatorStaggered.h"

__BEGIN_NAMESPACE

#pragma region kernels

/**
 * Wall-to-wall (W2W) contraction kernel.
 * Computes the single-field spatial reduction for a given (A, B) pair.
 *
 * p(t, A, B, d) = \sum _{c_1,c_2,x1,x2} [pick_even_{A,c_1}(DSource_{B,c_2}, x1)] * [pick_even_{A+d,c_1}(DdSource_{B+d,c_2}, x2)]^*
 *
 * w1[A*3 + c] is DSource_{A, c}
 * w2[A*3 + c] is DdSource_{A, c}
 *
 * so for each A, B, d, c1, c2
 * pick propagator need to be called two times, one with w1, the other with w2
 * first run: tmp1_{A,B,c_1,c_2} = \sum x pick_even_{A,c_1}(DSource_{B,c_2}, x) = \sum_{even (x+A)}(w1[B*3 + c2][c1])
 *            tmp2_{A,B,c_1,c_2} = \sum x pick_even_{A,c_1}(DdSource_{B,c_2}, x) = \sum_{even (x+A)}(w2[B*3 + c2][c1])
 *            note that one should calculate tmp1 and tmp2 at every site (stored in res) and do the reduced sum
 *
 * then, calculate p(t, A, B, d) = \sum _{c_1,c_2} tmp1_{A,B,c_1,c_2} * conj(tmp2_{A+d,B+d,c_1,c_2})
 * note that A+d and B+d will mod the cube (or the 3 components of the 3-vectors of (A+d) and (B+d) is xor instead of add, if we use bit as x,y,z since x,y,z=0 or 1, it was just A xor d or B xor d)
 *
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelPickPropagatorsW2W(
    const deviceSU3Vector* __restrict__ const * __restrict__ w,
    BYTE byA, BYTE byB, BYTE byFieldId,
    cuDoubleComplex* res)
{
    intokernalInt4;

    for (INT i = 0; i < 9; ++i)
    {
        res[i * _DC_Volume + uiSiteIndex] = _zerocd;
    }

    // Sum over even (sub-lattice) sites only
    if (0 != (sSite4.x & 1)
     || 0 != (sSite4.y & 1)
     || 0 != (sSite4.z & 1))
    {
        return;
    }

    // Extract components of A
    const BYTE ax = byA & 1;
    const BYTE ay = (byA >> 1) & 1;
    const BYTE az = (byA >> 2) & 1;

    // Compute site position: 2x + A
    SSmallInt4 siteA = sSite4;
    siteA.x += ax;
    siteA.y += ay;
    siteA.z += az;

    // Look up actual site index with boundary conditions
    const UINT uiSiteA = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siteA)].m_uiSiteIndex;

    // Keep the 3x3 color-pair channels separate.
    #pragma unroll
    for (BYTE c1 = 0; c1 < 3; ++c1)
    {
        #pragma unroll
        for (BYTE c2 = 0; c2 < 3; ++c2)
        {
            res[(c1 * 3 + c2) * _DC_Volume + uiSiteIndex] = _cToDouble(w[byB * 3 + c2][uiSiteA].m_ve[c1]);
        }
    }
}

/**
 * Point-to-point (P2P) contraction kernel.
 * Computes p(A,B,delta,t) per site for a given (A,B) pair.
 * No shift loop; no sign application.
 * Sum over 9 color combinations (c1, c2) internally.
 *
 * w1[A*3 + c] is DSource_{A, c}
 * w2[A*3 + c] is DdSource_{A, c}
 * p(t, A, B, d) = \sum _{c_1,c_2,x} pick_even_{A,c_1}(DSource_{B,c_2},x) * [pick_even_{A+d,c_1}(DdSource_{B+d,c_2},x)]^*
 *
 * so this will only called once for each A, B, d
 * note that one should calculate tmp(t, A, B, d, x) at every site (stored in res) and do the reduced sum
 * tmp(t, A, B, d, x) = \sum _{c_1,c_2} pick_even_{A,c_1}(DSource_{B,c_2},x) * [pick_even_{A+d,c_1}(DdSource_{B+d,c_2},x)]^*
 *                    = \sum _{c_1,c_2} w1[B*3 + c2][2x + A].m_ve[c1] * conj(w2[(B+d)*3 + c2][2x + (A+d)].m_ve[c1])
 * note that A+d and B+d will mod the cube (or the 3 components of the 3-vectors of (A+d) and (B+d) is xor instead of add, if we use bit as x,y,z since x,y,z=0 or 1, it was just A xor d or B xor d)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelPickPropagatorsP2P(
    const deviceSU3Vector* __restrict__ const * __restrict__ w1,
    const deviceSU3Vector* __restrict__ const * __restrict__ w2,
    BYTE byA, BYTE byB, BYTE byDelta, BYTE byFieldId,
    cuDoubleComplex* res)
{
    intokernalInt4;

    res[uiSiteIndex] = _zerocd;
    if (0 != (sSite4.x & 1)
     || 0 != (sSite4.y & 1)
     || 0 != (sSite4.z & 1))
    {
        return;
    }

    const BYTE byA_p_delta = byA ^ byDelta;
    const BYTE ax = byA & 1;
    const BYTE ay = (byA >> 1) & 1;
    const BYTE az = (byA >> 2) & 1;
    const BYTE apdx = byA_p_delta & 1;
    const BYTE apdy = (byA_p_delta >> 1) & 1;
    const BYTE apdz = (byA_p_delta >> 2) & 1;

    SSmallInt4 siteA = sSite4;
    siteA.x += ax;
    siteA.y += ay;
    siteA.z += az;
    SSmallInt4 siteApd = sSite4;
    siteApd.x += apdx;
    siteApd.y += apdy;
    siteApd.z += apdz;

    const UINT uiSiteA = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siteA)].m_uiSiteIndex;
    const UINT uiSiteApd = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siteApd)].m_uiSiteIndex;

    const SCHAR byB_p_delta = byB ^ byDelta;
    cuDoubleComplex sum = _zerocd;
    #pragma unroll
    for (BYTE c1 = 0; c1 < 3; ++c1)
    {
        #pragma unroll
        for (BYTE c2 = 0; c2 < 3; ++c2)
        {
            // w1 at source B, sink A (DSource)
            // w2 at source B^d, sink A^d (DdSource)
            sum = cuCadd(sum,
                cuCmul(
                    _cToDouble(w1[byB * 3 + c2][uiSiteA].m_ve[c1]),
                    cuConj(_cToDouble(w2[byB_p_delta * 3 + c2][uiSiteApd].m_ve[c1]))
                )
            );
        }
    }
    res[uiSiteIndex] = sum;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPickEveryTimeSlice(
    const cuDoubleComplex* __restrict__ pAll,
    BYTE uiT,
    INT uiChannel,
    cuDoubleComplex* res
)
{
    intokernalInt4_S(uiT);
    res[uiSiteIndex3D] = pAll[uiChannel * _DC_Volume + uiSiteIndex];
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CMeasureMesonCorrelatorStaggered)

#pragma region Help functions

/**
 * Initialize sign table and delta table matching MesonStructures.py.
 * Sign table: [type][sub][A] = (s_x*a_x + s_y*a_y + s_z*a_z) mod 2
 * Delta table: [type][sub] = spatial displacement index (0-7)
 * Values computed from eta_i/xi_i functions in MesonStructures.py.
 */
void CMeasureMesonCorrelatorStaggered::InitialSignTable(
    BYTE signtable[_kMesonCorrelatorType][_kMaxSubChannels][8],
    BYTE deltatable[_kMesonCorrelatorType][_kMaxSubChannels],
    BYTE nsubchannels[_kMesonCorrelatorType])
{
    // Sub-channel counts: types 0,1,16,17 have 1 sub-channel; rest have 3
    static const BYTE s_nSub[20] = {
        1, 1, 3, 3, 3, 3, 3, 3, 3, 3,
        3, 3, 3, 3, 3, 3, 1, 1, 3, 3
    };
    memcpy(nsubchannels, s_nSub, sizeof(s_nSub));

    // Sign table values from PyQUDA all_signs()
    static const BYTE s_sign[20][3][8] = {
        // Type 0: gamma_S (pseudoscalar), s=(0,0,0), delta=0
        {{0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0}},
        // Type 1: gamma_0 gamma_S, s=(1,1,1), delta=0
        {{0,1,1,0,1,0,0,1}, {0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0}},
        // Type 2: gamma_k gamma_5, 3 sub-channels
        {{0,1,0,1,0,1,0,1}, {0,0,1,1,0,0,1,1}, {0,0,0,0,1,1,1,1}},
        // Type 3: gamma_k gamma_0 gamma_5, 3 sub-channels
        {{0,0,1,1,1,1,0,0}, {0,1,0,1,1,0,1,0}, {0,1,1,0,0,1,1,0}},
        // Type 4: gamma_k, 3 sub-channels
        {{0,0,0,0,0,0,0,0}, {0,1,0,1,0,1,0,1}, {0,1,1,0,0,1,1,0}},
        // Type 5: gamma_0 gamma_k, 3 sub-channels
        {{0,1,1,0,1,0,0,1}, {0,0,1,1,1,1,0,0}, {0,0,0,0,1,1,1,1}},
        // Type 6: gamma_k gamma_5 gamma_j (k!=j), 3 sub-channels
        {{0,1,0,1,0,1,0,1}, {0,1,1,0,0,1,1,0}, {0,1,1,0,1,0,0,1}},
        // Type 7: gamma_0 gamma_k gamma_5 gamma_j, 3 sub-channels
        {{0,0,1,1,1,1,0,0}, {0,0,0,0,1,1,1,1}, {0,0,0,0,0,0,0,0}},
        // Type 8: gamma_k gamma_l (k!=l), 3 sub-channels
        {{0,0,0,0,0,0,0,0}, {0,1,0,1,0,1,0,1}, {0,0,0,0,1,1,1,1}},
        // Type 9: gamma_0 gamma_k gamma_l, 3 sub-channels
        {{0,1,1,0,1,0,0,1}, {0,0,1,1,1,1,0,0}, {0,1,1,0,0,1,1,0}},
        // Type 10: gamma_k gamma_l gamma_5, 3 sub-channels
        {{0,1,0,1,0,1,0,1}, {0,0,1,1,0,0,1,1}, {0,1,1,0,0,1,1,0}},
        // Type 11: gamma_0 gamma_k gamma_l gamma_5, 3 sub-channels
        {{0,0,1,1,1,1,0,0}, {0,1,0,1,1,0,1,0}, {0,0,0,0,1,1,1,1}},
        // Type 12: gamma_k gamma_j (k!=j), 3 sub-channels
        {{0,0,1,1,0,0,1,1}, {0,0,0,0,1,1,1,1}, {0,0,1,1,1,1,0,0}},
        // Type 13: gamma_0 gamma_k gamma_j, 3 sub-channels
        {{0,1,0,1,1,0,1,0}, {0,1,1,0,0,1,1,0}, {0,1,0,1,0,1,0,1}},
        // Type 14: gamma_k gamma_l gamma_5 gamma_j, 3 sub-channels
        {{0,1,1,0,1,0,0,1}, {0,1,1,0,1,0,0,1}, {0,0,0,0,0,0,0,0}},
        // Type 15: gamma_0 gamma_k gamma_l gamma_5 gamma_j, 3 sub-channels
        {{0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0}, {0,1,1,0,1,0,0,1}},
        // Type 16: gamma_1 gamma_2 gamma_3, s=(0,1,0), delta=7
        {{0,0,1,1,0,0,1,1}, {0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0}},
        // Type 17: gamma_0 gamma_1 gamma_2 gamma_3, s=(1,0,1), delta=7
        {{0,1,0,1,1,0,1,0}, {0,0,0,0,0,0,0,0}, {0,0,0,0,0,0,0,0}},
        // Type 18: gamma_k gamma_1 gamma_2 gamma_3 gamma_5, 3 sub-channels
        {{0,1,1,0,0,1,1,0}, {0,0,0,0,0,0,0,0}, {0,0,1,1,1,1,0,0}},
        // Type 19: gamma_0 gamma_k gamma_1 gamma_2 gamma_3 gamma_5, 3 sub-channels
        {{0,0,0,0,1,1,1,1}, {0,1,1,0,1,0,0,1}, {0,1,0,1,0,1,0,1}},
    };
    memcpy(signtable, s_sign, sizeof(s_sign));

    // Delta table values from all_deltas()
    // delta_to_deltaidx: (dz<<2)|(dy<<1)|dx
    static const BYTE s_delta[20][3] = {
        {0, 0, 0},  // Type 0
        {0, 0, 0},  // Type 1
        {0, 0, 0},  // Type 2
        {0, 0, 0},  // Type 3
        {1, 2, 4},  // Type 4: dx, dy, dz
        {1, 2, 4},  // Type 5
        {1, 2, 4},  // Type 6
        {1, 2, 4},  // Type 7
        {2, 4, 1},  // Type 8: dy, dz, dx
        {2, 4, 1},  // Type 9
        {3, 6, 5},  // Type 10: dx+dy, dy+dz, dx+dz
        {3, 6, 5},  // Type 11
        {3, 6, 5},  // Type 12
        {3, 6, 5},  // Type 13
        {3, 6, 5},  // Type 14
        {3, 6, 5},  // Type 15
        {7, 0, 0},  // Type 16: dx+dy+dz
        {7, 0, 0},  // Type 17
        {7, 7, 7},  // Type 18
        {7, 7, 7},  // Type 19
    };
    memcpy(deltatable, s_delta, sizeof(s_delta));
}

void CMeasureMesonCorrelatorStaggered::CalculateSources(INT gn, INT bn, const CFieldGauge* const* gs, const CFieldBoson* const* bs)
{
    for (BYTE shift = 0; shift < 8; ++shift)
    {
        for (BYTE c = 0; c < 3; ++c)
        {
            const INT idx = shift * 3 + c;
            SFermionBosonSource source;
            source.m_bySpinIndex = shift;
            source.m_byColorIndex = c;
            source.m_eSourceType = EFS_StaggeredWall;
            source.m_sSourcePoint = SSmallInt4(0, 0, 0, 0);
            m_pW1[idx]->InitialAsSource(source);
            appParanoiac(_T("generating source(%d)...\n"), idx);
            m_pW1[idx]->InverseDDdagger(gn, bn, 0, gs, bs, NULL);
            m_pW1[idx]->CopyTo(m_pW2[idx]);

            m_pW1[idx]->Ddagger(gn, bn, 0, gs, bs, NULL);
            m_pW2[idx]->D(gn, bn, 0, gs, bs, NULL);
        }
    }
}

/**
 * W2W contraction: for each (A, B), run kernel twice (DSource + DdSource),
 * reduce per time slice, then combine on CPU with delta.
 */
void CMeasureMesonCorrelatorStaggered::CalculateW2W()
{
    // 1. Build device pointer arrays
    deviceSU3Vector* w1[24];
    deviceSU3Vector* w2[24];
    for (BYTE shift = 0; shift < 24; ++shift)
    {
        w1[shift] = m_pW1[shift]->m_pDeviceData;
        w2[shift] = m_pW2[shift]->m_pDeviceData;
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceW1, w1, sizeof(deviceSU3Vector*) * 24, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(m_pDeviceW2, w2, sizeof(deviceSU3Vector*) * 24, cudaMemcpyHostToDevice));

    const INT nt = _HC_Lti;
    const BYTE fieldId = GetFermionFieldId();

    // 2. Allocate temporary per-(A,B) arrays:
    //    flat indexing [(A*8+B)*nt*9 + t*9 + colorPair]
    cuDoubleComplex* tmp1 = (cuDoubleComplex*)malloc(sizeof(cuDoubleComplex) * 64 * nt * 9);
    cuDoubleComplex* tmp2 = (cuDoubleComplex*)malloc(sizeof(cuDoubleComplex) * 64 * nt * 9);

    // 3. For each (A, B): compute tmp1 (DSource) and tmp2 (DdSource)
    for (INT A = 0; A < 8; ++A)
    {
        for (INT B = 0; B < 8; ++B)
        {
            // Kernel call 1: DSource -> tmp1
            {
                preparethread;
                _LAUNCH_KERNEL(_kernelPickPropagatorsW2W, block, threads,
                    m_pDeviceW1,
                    static_cast<BYTE>(A),
                    static_cast<BYTE>(B),
                    fieldId,
                    m_pDeviceKernelBuffer);

                for (INT ch = 0; ch < 9; ++ch)
                {
                    for (INT t = 0; t < nt; ++t)
                    {
                        preparethread_S;
                        _LAUNCH_KERNEL(_kernelPickEveryTimeSlice, block3d, threads3d,
                            m_pDeviceKernelBuffer,
                            static_cast<BYTE>(t),
                            ch,
                            m_pDeviceKernelEveryTimeSlice);
                        tmp1[(((A * 8 + B) * nt + t) * 9) + ch] = CCudaHelper::ReduceComplex(m_pDeviceKernelEveryTimeSlice, _HC_Volume_xyz);
                    }
                }
            }

            // Kernel call 2: DdSource -> tmp2
            {
                preparethread;
                _LAUNCH_KERNEL(_kernelPickPropagatorsW2W, block, threads,
                    m_pDeviceW2,
                    static_cast<BYTE>(A),
                    static_cast<BYTE>(B),
                    fieldId,
                    m_pDeviceKernelBuffer);

                for (INT ch = 0; ch < 9; ++ch)
                {
                    for (INT t = 0; t < nt; ++t)
                    {
                        preparethread_S;
                        _LAUNCH_KERNEL(_kernelPickEveryTimeSlice, block3d, threads3d,
                            m_pDeviceKernelBuffer,
                            static_cast<BYTE>(t),
                            ch,
                            m_pDeviceKernelEveryTimeSlice);
                        tmp2[(((A * 8 + B) * nt + t) * 9) + ch] = CCudaHelper::ReduceComplex(m_pDeviceKernelEveryTimeSlice, _HC_Volume_xyz);
                    }
                }
            }
        }
    }

    // 4. Combine p_w2w array on CPU:
    //    p(t, A, B, d) = sum_{c1,c2} tmp1(A,B,t,c1,c2) * conj(tmp2(A^d,B^d,t,c1,c2))
    for (INT t = 0; t < nt; ++t)
    {
        for (INT A = 0; A < 8; ++A)
        {
            for (INT B = 0; B < 8; ++B)
            {
                for (INT d = 0; d < 8; ++d)
                {
                    const INT Ad = A ^ d;
                    const INT Bd = B ^ d;
                    cuDoubleComplex corr = make_cuDoubleComplex(0.0, 0.0);
                    for (INT ch = 0; ch < 9; ++ch)
                    {
                        const cuDoubleComplex& lhs = tmp1[(((A * 8 + B) * nt + t) * 9) + ch];
                        const cuDoubleComplex& rhs = tmp2[(((Ad * 8 + Bd) * nt + t) * 9) + ch];
                        corr = cuCadd(corr, cuCmul(lhs, cuConj(rhs)));
                    }
                    m_pW2WPArray[((t * 8 + A) * 8 + B) * 8 + d] = corr;
                }
            }
        }
    }

    // 5. Project onto meson channels
    ProjectPArray(m_pW2WPArray, m_pW2WCorrelator);

    free(tmp1);
    free(tmp2);
}

/**
 * P2P contraction: compute p(A,B,delta,t) for all (A,B,delta),
 * then project onto meson channels on CPU.
 */
void CMeasureMesonCorrelatorStaggered::CalculateP2P()
{
    // 1. Build device pointer arrays (may already be set by CalculateW2W)
    deviceSU3Vector* w1[24];
    deviceSU3Vector* w2[24];
    for (BYTE shift = 0; shift < 24; ++shift)
    {
        w1[shift] = m_pW1[shift]->m_pDeviceData;
        w2[shift] = m_pW2[shift]->m_pDeviceData;
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceW1, w1, sizeof(deviceSU3Vector*) * 24, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(m_pDeviceW2, w2, sizeof(deviceSU3Vector*) * 24, cudaMemcpyHostToDevice));

    const INT nt = _HC_Lti;
    const BYTE fieldId = GetFermionFieldId();

    // 2. For each (A, B, delta): launch P2P kernel, reduce per time slice
    for (INT A = 0; A < 8; ++A)
    {
        for (INT B = 0; B < 8; ++B)
        {
            for (INT delta = 0; delta < 8; ++delta)
            {
                preparethread;
                _LAUNCH_KERNEL(_kernelPickPropagatorsP2P, block, threads,
                    m_pDeviceW1,
                    m_pDeviceW2,
                    static_cast<BYTE>(A),
                    static_cast<BYTE>(B),
                    static_cast<BYTE>(delta),
                    fieldId,
                    m_pDeviceKernelBuffer);

                for (INT t = 0; t < nt; ++t)
                {
                    preparethread_S;
                    _LAUNCH_KERNEL(_kernelPickEveryTimeSlice, block3d, threads3d,
                        m_pDeviceKernelBuffer,
                        static_cast<BYTE>(t),
                        0,
                        m_pDeviceKernelEveryTimeSlice);
                    m_pP2PPArray[((t * 8 + A) * 8 + B) * 8 + delta] = CCudaHelper::ReduceComplex(m_pDeviceKernelEveryTimeSlice, _HC_Volume_xyz);
                }
            }
        }
    }

    // 3. Project onto meson channels
    ProjectPArray(m_pP2PPArray, m_pP2PCorrelator);
}

/**
 * projection: C_d(t) = sum_{A,B} sign(A)sign(B) p(t,A,B,d)
 */
void CMeasureMesonCorrelatorStaggered::ProjectPArray(const cuDoubleComplex* pArray, cuDoubleComplex* correlator)
{
    const INT nt = _HC_Lti;
    for (INT ty = 0; ty < 20; ++ty)
    {
        for (INT sub = 0; sub < m_nSubChannels[ty]; ++sub)
        {
            const BYTE d = m_pDeltaTable[ty][sub];
            for (INT t = 0; t < nt; ++t)
            {
                cuDoubleComplex corr = make_cuDoubleComplex(0.0, 0.0);
                for (INT A = 0; A < 8; ++A)
                {
                    const BYTE sa = m_pSignTable[ty][sub][A];
                    for (INT B = 0; B < 8; ++B)
                    {
                        const BYTE sb = m_pSignTable[ty][sub][B];
                        const INT signTotal = sa + sb;
                        const cuDoubleComplex& pval = pArray[((t * 8 + A) * 8 + B) * 8 + d];
                        const cuDoubleComplex pd = make_cuDoubleComplex(pval.x, pval.y);
                        if (signTotal & 1)
                        {
                            corr = cuCsub(corr, pd);
                        }
                        else
                        {
                            corr = cuCadd(corr, pd);
                        }
                    }
                }
                correlator[(ty * 3 + sub) * nt + t] = corr;
            }
        }
    }
}

void CMeasureMesonCorrelatorStaggered::InitialBuffers()
{
    const INT nt = _HC_Lti;

    // Device pointer arrays
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceW1, sizeof(deviceSU3Vector*) * 24));
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceW2, sizeof(deviceSU3Vector*) * 24));

    // Shared GPU kernel buffers
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceKernelBuffer, sizeof(cuDoubleComplex) * _HC_Volume * 9));
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceKernelEveryTimeSlice,
        sizeof(cuDoubleComplex) * _HC_Volume_xyz));

    // Host per-config buffers
    m_pW2WPArray = (cuDoubleComplex*)(malloc(sizeof(cuDoubleComplex) * nt * 8 * 8 * 8));
    m_pW2WCorrelator = (cuDoubleComplex*)(malloc(sizeof(cuDoubleComplex) * _kMesonCorrelatorType * _kMaxSubChannels * nt));
    m_pP2PPArray = (cuDoubleComplex*)(malloc(sizeof(cuDoubleComplex) * nt * 8 * 8 * 8));
    m_pP2PCorrelator = (cuDoubleComplex*)(malloc(sizeof(cuDoubleComplex) * _kMesonCorrelatorType * _kMaxSubChannels * nt));
}

#pragma endregion

CMeasureMesonCorrelatorStaggered::~CMeasureMesonCorrelatorStaggered()
{
    checkCudaErrors(__cudaFree(m_pDeviceW1));
    checkCudaErrors(__cudaFree(m_pDeviceW2));
    checkCudaErrors(__cudaFree(m_pDeviceKernelBuffer));
    checkCudaErrors(__cudaFree(m_pDeviceKernelEveryTimeSlice));

    appSafeFree(m_pW2WPArray);
    appSafeFree(m_pW2WCorrelator);
    appSafeFree(m_pP2PPArray);
    appSafeFree(m_pP2PCorrelator);
}

void CMeasureMesonCorrelatorStaggered::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);
    INT iValue = 1;
    param.FetchValueINT(_T("GaugeFixing"), iValue);
    m_bGaugeFixing = iValue != 0;

    InitialBuffers();
    InitialSignTable(m_pSignTable, m_pDeltaTable, m_nSubChannels);
}

void CMeasureMesonCorrelatorStaggered::OnConfigurationAccepted(INT gn, INT bn, INT tensor2Num, const CFieldGauge* const* gs, const CFieldBoson* const* bs, const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* stp)
{
    for (BYTE i = 0; i < 24; ++i)
    {
        m_pW1[i] = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(GetFermionFieldId(), _T(__FILE__), __LINE__));
        m_pW2[i] = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(GetFermionFieldId(), _T(__FILE__), __LINE__));
    }

    if (m_bGaugeFixing && NULL != appGetLattice()->m_pGaugeFixing)
    {
        TArray<CFieldGauge*> fixedgauges;
        for (INT i = 0; i < gn; ++i)
        {
            CFieldGauge* fixedgauge = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(gs[i], _T("CMeasureMesonCorrelatorStaggered"), __LINE__));
            appGetLattice()->m_pGaugeFixing->GaugeFixing(fixedgauge);
            fixedgauges.AddItem(fixedgauge);
        }

        CalculateSources(gn, bn, fixedgauges.GetData(), bs);

        for (INT i = 0; i < gn; ++i)
        {
            fixedgauges[i]->Return();
        }
    }
    else
    {
        CalculateSources(gn, bn, gs, bs);
    }

    // Compute W2W and P2P contractions
    CalculateW2W();
    CalculateP2P();

    // ========== extract results ===========
    const INT nt = _HC_Lti;

    if (m_bShowResult)
    {
        appPushLogDate(FALSE);
        appGeneral(_T("==================== W2W correlators ===============\n"));
    }

    TArray<TArray<TArray<cuDoubleComplex>>> thisConfW2W;
    TArray<TArray<TArray<cuDoubleComplex>>> thisConfP2P;
    TArray<TArray<cuDoubleComplex>> thisConfW2WCombined;
    TArray<TArray<cuDoubleComplex>> thisConfP2PCombined;

    for (INT i = 0; i < _kMesonCorrelatorType; ++i)
    {
        TArray<TArray<cuDoubleComplex>> thisTypeW2W;
        TArray<TArray<cuDoubleComplex>> thisTypeP2P;
        TArray<cuDoubleComplex> thisTypeW2WCombined;
        TArray<cuDoubleComplex> thisTypeP2PCombined;

        for (INT sub = 0; sub < m_nSubChannels[i]; ++sub)
        {
            const INT w2wOffset = (i * _kMaxSubChannels + sub) * nt;
            const INT p2pOffset = (i * _kMaxSubChannels + sub) * nt;

            if (m_bShowResult)
            {
                appGeneral(_T("Type%d_sub%d(W2W):"), i, sub);
            }
            TArray<cuDoubleComplex> thisSubW2W;
            TArray<cuDoubleComplex> thisSubP2P;
            for (INT t = 0; t < nt; ++t)
            {
                const cuDoubleComplex& w2wRes = m_pW2WCorrelator[w2wOffset + t];
                const cuDoubleComplex& p2pRes = m_pP2PCorrelator[p2pOffset + t];

                if (m_bShowResult)
                {
                    LogGeneralComplex(w2wRes);
                }
                thisSubW2W.AddItem(w2wRes);
                thisSubP2P.AddItem(p2pRes);
            }
            thisTypeW2W.AddItem(thisSubW2W);
            thisTypeP2P.AddItem(thisSubP2P);
            if (m_bShowResult)
            {
                appGeneral(_T("\n"));
            }
        }

        // Combined = sum over sub-channels
        for (INT t = 0; t < nt; ++t)
        {
            cuDoubleComplex combined_w2w = make_cuDoubleComplex(0.0, 0.0);
            cuDoubleComplex combined_p2p = make_cuDoubleComplex(0.0, 0.0);
            for (INT sub = 0; sub < m_nSubChannels[i]; ++sub)
            {
                combined_w2w = cuCadd(combined_w2w, thisTypeW2W[sub][t]);
                combined_p2p = cuCadd(combined_p2p, thisTypeP2P[sub][t]);
            }
            thisTypeW2WCombined.AddItem(combined_w2w);
            thisTypeP2PCombined.AddItem(combined_p2p);
        }

        thisConfW2W.AddItem(thisTypeW2W);
        thisConfP2P.AddItem(thisTypeP2P);
        thisConfP2PCombined.AddItem(thisTypeP2PCombined);
        thisConfW2WCombined.AddItem(thisTypeW2WCombined);

        if (m_bShowResult)
        {
            appGeneral(_T("Type%d(combined W2W):"), i);
            for (INT t = 0; t < nt; ++t)
            {
                LogGeneralComplex(thisTypeW2WCombined[t]);
            }
            appGeneral(_T("\nType%d(combined P2P):"), i);
            for (INT t = 0; t < nt; ++t)
            {
                LogGeneralComplex(thisTypeP2PCombined[t]);
            }
            appGeneral(_T("\n"));
        }
    }

    m_lstW2WCorrelator.AddItem(thisConfW2W);
    m_lstW2WCombinedCorrelator.AddItem(thisConfW2WCombined);
    m_lstP2PCorrelator.AddItem(thisConfP2P);
    m_lstP2PCombinedCorrelator.AddItem(thisConfP2PCombined);

    if (NULL != m_pOwner)
    {
        m_pOwner->AddOneConfigurationResult(this, _T("W2WCorrelator"), thisConfW2W);
        m_pOwner->AddOneConfigurationResult(this, _T("W2WCombinedCorrelator"), thisConfW2WCombined);
        m_pOwner->AddOneConfigurationResult(this, _T("P2PCorrelator"), thisConfP2P);
        m_pOwner->AddOneConfigurationResult(this, _T("P2PCombinedCorrelator"), thisConfP2PCombined);
    }

    if (m_bShowResult)
    {
        appPopLogDate();
    }

    for (BYTE i = 0; i < 24; ++i)
    {
        m_pW1[i]->Return();
        m_pW2[i]->Return();
    }
    ++m_uiConfigurationCount;
}

void CMeasureMesonCorrelatorStaggered::Report()
{
    appPushLogDate(FALSE);
    appGeneral(_T(" =====================================================\n"));
    appGeneral(_T(" ========= Staggered Meson (W2W + P2P) ===============\n"));
    appGeneral(_T(" =====================================================\n\n"));

    const INT nConf = m_lstW2WCombinedCorrelator.Num();
    const INT nt = _HC_Lti;

    // ===== W2W combined correlators =====
    m_lstAverageResults.RemoveAll();
    for (INT ty = 0; ty < _kMesonCorrelatorType; ++ty)
    {
        appGeneral(_T("(* === W2W Combined Type:%d === *)\ntabres%d={\n"), ty, ty);
        TArray<DOUBLE> thisType;
        for (INT conf = 0; conf < nConf; ++conf)
        {
            appGeneral(_T("{"));
            for (INT t = 0; t < nt; ++t)
            {
                appGeneral(_T("%2.12f%s"), m_lstW2WCombinedCorrelator[conf][ty][t].x,
                    (t != (nt - 1)) ? _T(",") : _T(""));
                if (0 == conf)
                {
                    thisType.AddItem(m_lstW2WCombinedCorrelator[conf][ty][t].x);
                }
                else
                {
                    thisType[t] = thisType[t] + m_lstW2WCombinedCorrelator[conf][ty][t].x;
                }
            }
            appGeneral(_T("}%s"), (conf == nConf - 1) ? _T("\n};\n") : _T(",\n"));
        }
        for (INT t = 0; t < nt; ++t)
        {
            thisType[t] = thisType[t] / nConf;
        }
        m_lstAverageResults.AddItem(thisType);
    }

    // ===== Per-sub-channel W2W correlators =====
    appGeneral(_T("\n(* === W2W Per-Sub-Channel Correlators === *)\n"));
    for (INT ty = 0; ty < _kMesonCorrelatorType; ++ty)
    {
        for (INT sub = 0; sub < m_nSubChannels[ty]; ++sub)
        {
            appGeneral(_T("(* W2W Type:%d Sub:%d *)\n"), ty, sub);
            for (INT conf = 0; conf < nConf; ++conf)
            {
                appGeneral(_T("{"));
                for (INT t = 0; t < nt; ++t)
                {
                    appGeneral(_T("%2.12f%s"),
                        m_lstW2WCorrelator[conf][ty][sub][t].x,
                        (t != (nt - 1)) ? _T(",") : _T(""));
                }
                appGeneral(_T("}%s"), (conf == nConf - 1) ? _T("\n") : _T(",\n"));
            }
        }
    }

    // ===== Per-sub-channel P2P correlators =====
    appGeneral(_T("\n(* === P2P Per-Sub-Channel Correlators === *)\n"));
    for (INT ty = 0; ty < _kMesonCorrelatorType; ++ty)
    {
        for (INT sub = 0; sub < m_nSubChannels[ty]; ++sub)
        {
            appGeneral(_T("(* P2P Type:%d Sub:%d *)\n"), ty, sub);
            for (INT conf = 0; conf < nConf; ++conf)
            {
                appGeneral(_T("{"));
                for (INT t = 0; t < nt; ++t)
                {
                    appGeneral(_T("%2.12f%s"),
                        m_lstP2PCorrelator[conf][ty][sub][t].x,
                        (t != (nt - 1)) ? _T(",") : _T(""));
                }
                appGeneral(_T("}%s"), (conf == nConf - 1) ? _T("\n") : _T(",\n"));
            }
        }
    }

    // ===== P2P Combined correlators =====
    appGeneral(_T("\n(* === P2P Combined correlators === *)\n"));
    for (INT ty = 0; ty < _kMesonCorrelatorType; ++ty)
    {
        appGeneral(_T("(* P2P Combined Type:%d *)\n"), ty);
        for (INT conf = 0; conf < nConf; ++conf)
        {
            appGeneral(_T("{"));
            for (INT t = 0; t < nt; ++t)
            {
                appGeneral(_T("%2.12f%s"),
                    m_lstP2PCombinedCorrelator[conf][ty][t].x,
                    (t != (nt - 1)) ? _T(",") : _T(""));
            }
            appGeneral(_T("}%s"), (conf == nConf - 1) ? _T("\n") : _T(",\n"));
        }
    }

    // ===== Average combined W2W results =====
    appGeneral(_T("\n(* === All W2W Combined averages === *)\navr_w2w={\n"));
    for (INT ty = 0; ty < _kMesonCorrelatorType; ++ty)
    {
        appGeneral(_T("{"));
        for (INT t = 0; t < nt; ++t)
        {
            appGeneral(_T("%2.12f%s"), m_lstAverageResults[ty][t],
                (t != (nt - 1)) ? _T(",") : _T(""));
        }
        appGeneral(_T("}%s"), ty == 19 ? _T("\n};\n") : _T(",\n"));
    }

    appPopLogDate();
}

void CMeasureMesonCorrelatorStaggered::Reset()
{
    CMeasure::Reset();
    m_lstW2WCorrelator.RemoveAll();
    m_lstW2WCombinedCorrelator.RemoveAll();
    m_lstP2PCorrelator.RemoveAll();
    m_lstP2PCombinedCorrelator.RemoveAll();
    m_lstAverageResults.RemoveAll();
}


__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
