//=============================================================================
// FILENAME : CGaugeFixingMCGIndirect.cu
//
// DESCRIPTION:
// Indirect Maximal Center Gauge (MCG) fixing for SU(3).
//
// Two-step algorithm (Brower et al., hep-lat/9708008; Langfeld et al., hep-lat/0003021):
//   Step 1: Fix to Maximal Abelian Gauge (MAG).
//           This makes the gauge field as diagonal as possible.
//   Step 2: Fix to Maximal Center Gauge (MCG) on the MAG-fixed field.
//           This further aligns the field toward Z_3 center elements.
//
// Each step uses the Cabibbo-Marinari-Okawa method with red-black
// checkerboard sweeps and overrelaxation.
//
// Convergence criterion:
//   Step 1: theta_MAG = sum off-diagonal / (V * D) < accuracy
//   Step 2: theta_MCG = sum (9 - |Tr U|^2) / (9 * V * D) < accuracy
//
// References:
//  hep-lat/9708008 -- Brower et al., "Center vortices in SU(2) lattice gauge fields"
//  hep-lat/0003021 -- Langfeld et al., "SU(N) vortices and Wilson loops"
//  hep-lat/9906010 -- Montero, "Study of SU(3) vortex-like configurations"
//  hep-lat/0110165 -- Tucker & Stack, "The Maximal Abelian Gauge in SU(3)"
//
// ============================================================================
// FORMULA DERIVATION -- Standard IMCG local objective function
// ============================================================================
//
// The MCG functional (global) is:
//   R = sum_{x,mu} |Tr U_mu(x)|^2                         ... (A)
//
// Ref: Montero, hep-lat/9906010, Eq. (2):
//   "The maximal center gauge in SU(N) lattice gauge theory is defined as
//    the gauge which brings link variables U as close as possible to elements
//    of its center Z_N. This can be achieved by maximizing ... |Tr U_mu(n)|^2"
//
// For the local update at site x, Montero (hep-lat/9906010, Eq. (3)) gives
// the direct-MCG local objective for a FULL SU(N) gauge transformation G(x):
//   R_x = sum_mu |Tr{ G(x) U_mu(x) }|^2
//       + sum_mu |Tr{ U_mu(x-mu) G^dagger(x) }|^2          ... (B)
//
// In our code, we restrict to DIAGONAL gauge transformations:
//   g(x) = diag(z1(x), z2(x), z3(x)),  |zi| = 1,  z1*z2*z3 = 1
//
// We use red-black IN-PLACE updates. After each half-sweep, the field pU
// already contains the cumulative gauge transformations from all previous
// steps. When computing a NEW increment g(x) at site x, the intermediate
// state after applying g(x) (but before neighbors are updated in this sweep)
// contributes:
//
//   Forward link from x to x+mu:
//     g(x) * pU_mu(x)          [neighbors fixed]
//
//   Backward link ending at x (from x-mu):
//     pU_mu(x-mu) * g^dagger(x)   [neighbors fixed]
//
// The local objective for diagonal g(x) on the CURRENT field pU is:
//
//   R_x = sum_mu | sum_i zi * pU_mu^{ii}(x) |^2
//       + sum_mu | sum_i pU_mu^{ii}(x-mu) * zi* |^2       ... (C)
//
// where pU_mu^{ii}(x) are the diagonal elements of the current field.
//
// WHY this is correct in the in-place update framework:
//   pU_mu(x) already contains G_cum(x) U_mu^{orig}(x) G_cum^dagger(x+mu).
//   The diagonal elements satisfy:
//     pU_mu^{ii}(x) = G_cum^{ii}(x) * U_mu^{ii,orig}(x) * G_cum^{ii,*}(x+mu)
//   because diagonal gauge transformations commute with the diagonal part:
//     [G U G^dagger]_{ii} = G_{ii} U_{ii} G_{ii}^* = U_{ii}   (for |G_{ii}|=1)
//
//   After applying the new increment g(x), the full transformed link would be:
//     g(x) * pU_mu(x) * g^dagger(x+mu)   [for forward]
//   but g(x+mu) is fixed in this half-sweep, so the trace is:
//     Tr(g(x) * pU_mu(x)) = sum_i g_i(x) * pU_mu^{ii}(x)
//   which is exactly the first term of (C).
//
//   Similarly for the backward link:
//     Tr(pU_mu(x-mu) * g^dagger(x)) = sum_i pU_mu^{ii}(x-mu) * g_i^*(x)
//   which is the second term of (C).
//
// Ref: Tucker & Stack, hep-lat/0110165, Eq. (5):
//   "This gauge is determined by maximizing the functional
//    G_{SU(3)}^{icg} = (1/9N_link) sum_{x,mu} |Tr[U_mu(x)]|^2"
//   Note: Tucker & Stack first PROJECT to U(1)xU(1) before this step.
//   Our variant skips projection and works with the full field restricted
//   to diagonal gauge transformations.
//
// The kernel _kernelCalculateStandardIMCGSU3 implements exactly formula (C).
// The kernel _kernelCalculateStandardIMCGSU2 implements the SU(2) analogue.
//
// REVISION:
//  [05/14/2026 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CGaugeFixingMCGIndirect.h"
#include "CGaugeFixingMAG.h"
#include "CGaugeFixingMCGDirect.h"

__BEGIN_NAMESPACE

#define _CLG_STANDARD_IMCG_GRID 24

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformResidualOdd(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pGx,
    deviceGauge* pGauge)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (sSite4.IsOdd())
    {
        if (!site.IsDirichlet())
        {
            for (BYTE dir = 0; dir < uiDir; ++dir)
            {
                UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
                pGauge[uiLinkDir] = pGx[uiSiteIndex].MulC(pGauge[uiLinkDir]);
            }
        }
    }
    else
    {
        for (BYTE dir = 0; dir < uiDir; ++dir)
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
            {
                UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
                const SSmallInt4 p_p_mu_site = _deviceSmallInt4OffsetC(sSite4, dir + 1);
                const SIndex& site_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_p_mu_site)];
                if (!site_p_mu.IsDirichlet())
                {
                    pGauge[uiLinkDir].MulDagger(pGx[site_p_mu.m_uiSiteIndex]);
                }
            }
        }
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformResidualEven(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pGx,
    deviceGauge* pGauge)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (!sSite4.IsOdd())
    {
        if (!site.IsDirichlet())
        {
            for (BYTE dir = 0; dir < uiDir; ++dir)
            {
                UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
                pGauge[uiLinkDir] = pGx[uiSiteIndex].MulC(pGauge[uiLinkDir]);
            }
        }
    }
    else
    {
        for (BYTE dir = 0; dir < uiDir; ++dir)
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
            {
                UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
                const SSmallInt4 p_p_mu_site = _deviceSmallInt4OffsetC(sSite4, dir + 1);
                const SIndex& site_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_p_mu_site)];
                if (!site_p_mu.IsDirichlet())
                {
                    pGauge[uiLinkDir].MulDagger(pGx[site_p_mu.m_uiSiteIndex]);
                }
            }
        }
    }
}

static __device__ __inline__ CLGComplex _devicePhaseFromGrid(INT i, INT grid)
{
    const Real theta = F(2.0) * F(3.14159265358979323846) * static_cast<Real>(i) / static_cast<Real>(grid);
    return _make_cuComplex(_cos(theta), _sin(theta));
}

static __device__ __inline__ DOUBLE _deviceAbsSq3(
    const CLGComplex& z1, const CLGComplex& d1,
    const CLGComplex& z2, const CLGComplex& d2,
    const CLGComplex& z3, const CLGComplex& d3)
{
    CLGComplex tr = _cuCaddf(_cuCaddf(_cuCmulf(z1, d1), _cuCmulf(z2, d2)), _cuCmulf(z3, d3));
    return static_cast<DOUBLE>(__cuCabsSqf(tr));
}

// ============================================================================
// _kernelCalculateStandardIMCGSU3
//
// Computes the optimal diagonal gauge transformation g(x) = diag(z1, z2, z3)
// by grid search, maximizing the local objective function (formula C above):
//
//   R_x = sum_mu | sum_i zi * pU_mu^{ii}(x) |^2
//       + sum_mu | sum_i pU_mu^{ii}(x-mu) * zi* |^2
//
// The search is over a uniform grid of phase angles:
//   theta = 2*pi*k / _CLG_STANDARD_IMCG_GRID,  k = 0, 1, ..., grid-1
//
// For SU(3), z1 and z2 are searched independently (grid^2 combinations),
// and z3 = conj(z1 * z2) enforces det(g) = 1.
//
// pU is the CURRENT gauge field (already transformed by all previous sweeps).
// pG stores the NEW increment g(x) for sites updated in this half-sweep.
//
// The red-black decomposition ensures that no thread reads pG from a site
// that is being written by another thread in the same kernel launch.
// ============================================================================

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateStandardIMCGSU3(
    BYTE byFieldId,
    UBOOL bOdd,
    const deviceSU3* __restrict__ pU,
    deviceSU3* pG,
    INT iGrid,
    Real fOmega)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (site.IsDirichlet() || (sSite4.IsOdd() != bOdd))
    {
        return;
    }

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    DOUBLE fBest = -1.0;
    INT iBestA = 0;
    INT iBestB = 0;

    for (INT ia = 0; ia < iGrid; ++ia)
    {
        const CLGComplex z1 = _devicePhaseFromGrid(ia, iGrid);
        const CLGComplex z1dag = _cuConjf(z1);
        for (INT ib = 0; ib < iGrid; ++ib)
        {
            const CLGComplex z2 = _devicePhaseFromGrid(ib, iGrid);
            const CLGComplex z2dag = _cuConjf(z2);
            const CLGComplex z3 = _cuConjf(_cuCmulf(z1, z2));
            const CLGComplex z3dag = _cuConjf(z3);

            DOUBLE fObjective = 0.0;
            for (BYTE dir = 0; dir < uiDir; ++dir)
            {
                if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
                {
                    const deviceSU3& Ufwd = pU[_deviceGetLinkIndex(uiSiteIndex, dir)];
                    fObjective += _deviceAbsSq3(z1, Ufwd.m_me[0], z2, Ufwd.m_me[4], z3, Ufwd.m_me[8]);
                }

                const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
                const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];
                const deviceSU3 Ubwd = _deviceGetGaugeBCSU3DirOneSIndex(pU, site_m_mu);
                fObjective += _deviceAbsSq3(z1dag, Ubwd.m_me[0], z2dag, Ubwd.m_me[4], z3dag, Ubwd.m_me[8]);
            }

            if (fObjective > fBest)
            {
                fBest = fObjective;
                iBestA = ia;
                iBestB = ib;
            }
        }
    }

    // Apply overrelaxation: scale phase angles by fOmega
    const Real fTwoPiOverGrid = F(2.0) * PI / static_cast<Real>(iGrid);
    const Real theta1 = fTwoPiOverGrid * static_cast<Real>(iBestA);
    const Real theta2 = fTwoPiOverGrid * static_cast<Real>(iBestB);
    const Real theta1_o = fOmega * theta1;
    const Real theta2_o = fOmega * theta2;

    const CLGComplex z1 = _make_cuComplex(_cos(theta1_o), _sin(theta1_o));
    const CLGComplex z2 = _make_cuComplex(_cos(theta2_o), _sin(theta2_o));
    const CLGComplex z3 = _cuConjf(_cuCmulf(z1, z2));

    deviceSU3 g = deviceSU3::makeSU3Zero();
    g.m_me[0] = z1;
    g.m_me[4] = z2;
    g.m_me[8] = z3;
    pG[uiSiteIndex] = g;
}

static __device__ __inline__ DOUBLE _deviceAbsSq2(
    const CLGComplex& z, const CLGComplex& a,
    const CLGComplex& zdag, const CLGComplex& d)
{
    CLGComplex tr = _cuCaddf(_cuCmulf(z, a), _cuCmulf(zdag, d));
    return static_cast<DOUBLE>(__cuCabsSqf(tr));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateStandardIMCGSU2(
    BYTE byFieldId,
    UBOOL bOdd,
    const deviceSU2* __restrict__ pU,
    deviceSU2* pG,
    INT iGrid,
    Real fOmega)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (site.IsDirichlet() || (sSite4.IsOdd() != bOdd))
    {
        return;
    }

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    DOUBLE fBest = -1.0;
    INT iBest = 0;

    for (INT ia = 0; ia < iGrid; ++ia)
    {
        const CLGComplex z = _devicePhaseFromGrid(ia, iGrid);
        const CLGComplex zdag = _cuConjf(z);
        DOUBLE fObjective = 0.0;

        for (BYTE dir = 0; dir < uiDir; ++dir)
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
            {
                const deviceSU2& Ufwd = pU[_deviceGetLinkIndex(uiSiteIndex, dir)];
                fObjective += _deviceAbsSq2(z, Ufwd.m_me[0], zdag, Ufwd.m_me[3]);
            }

            const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
            const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];
            const deviceSU2 Ubwd = _deviceGetGaugeBCDirOneSIndexT<deviceSU2>(pU, site_m_mu);
            fObjective += _deviceAbsSq2(zdag, Ubwd.m_me[0], z, Ubwd.m_me[3]);
        }

        if (fObjective > fBest)
        {
            fBest = fObjective;
            iBest = ia;
        }
    }

    // Apply overrelaxation: scale phase angle by fOmega
    const Real fTwoPiOverGrid = F(2.0) * PI / static_cast<Real>(iGrid);
    const Real theta = fTwoPiOverGrid * static_cast<Real>(iBest);
    const Real theta_o = fOmega * theta;

    const CLGComplex z = _make_cuComplex(_cos(theta_o), _sin(theta_o));
    deviceSU2 g = deviceSU2::makeSU2Zero();
    g.m_me[0] = z;
    g.m_me[3] = _cuConjf(z);
    pG[uiSiteIndex] = g;
}

// ============================================================================
// _StandardIMCGGaugeFixingSU3 / _StandardIMCGGaugeFixingSU2
//
// Host-side driver for the standard IMCG (residual-Abelian) gauge fixing.
// Implements a red-black checkerboard sweep:
//
//   For each iteration:
//     1. Compute g for odd sites  (_kernelCalculateStandardIMCGSU3/SU2)
//     2. Apply g to odd sites     (_kernelGaugeTransformResidualOdd)
//     3. Compute g for even sites (_kernelCalculateStandardIMCGSU3/SU2)
//     4. Apply g to even sites    (_kernelGaugeTransformResidualEven)
//
// The convergence check uses CGaugeFixingMCGDirect::CheckRes(), which
// measures theta = sum(9 - |Tr U|^2) / (9 * V * D).
//
// Note: pG is allocated but NOT explicitly initialized. However, in the
// red-black pattern, every read from pG accesses a site that was written
// in the immediately preceding kernel launch, so uninitialized memory is
// never read.
// ============================================================================

static void _StandardIMCGGaugeFixingSU3(CGaugeFixingMCGIndirect* pFixer, deviceSU3* pData, BYTE byFieldId)
{
    deviceSU3* pG = NULL;
    checkCudaErrors(__cudaMalloc((void**)&pG, _HC_Volume * sizeof(deviceSU3)));

    preparethread;
    pFixer->m_iIterate = 0;
    CGaugeFixingMCGDirect checker;

    while (pFixer->m_iIterate < pFixer->m_iMaxIterate)
    {
        if (0 == pFixer->m_iIterate % pFixer->m_iCheckErrorStep)
        {
            const Real fTheta = static_cast<Real>(checker.CheckResLocalSU3(pData, byFieldId));
            appGeneral(_T("MCG Indirect residual-Abelian (SU3) Iterate : %d, theta = %2.12f\n"), pFixer->m_iIterate, fTheta);
            if (fTheta < pFixer->m_fAccuracy)
            {
                cudaSafeFree(pG);
                return;
            }
        }

        _LAUNCH_KERNEL(_kernelCalculateStandardIMCGSU3, block, threads,
            byFieldId, TRUE, pData, pG,
            pFixer->m_iIMCGGrid, pFixer->m_fOmega);
        _LAUNCH_KERNEL(_kernelGaugeTransformResidualOdd<deviceSU3>, block, threads,
            byFieldId, pG, pData);
        _LAUNCH_KERNEL(_kernelCalculateStandardIMCGSU3, block, threads,
            byFieldId, FALSE, pData, pG,
            pFixer->m_iIMCGGrid, pFixer->m_fOmega);
        _LAUNCH_KERNEL(_kernelGaugeTransformResidualEven<deviceSU3>, block, threads,
            byFieldId, pG, pData);

        ++pFixer->m_iIterate;
    }

    appGeneral(_T("MCG Indirect residual-Abelian (SU3) gauge fixing failed.\n"));
    cudaSafeFree(pG);
}

static void _StandardIMCGGaugeFixingSU2(CGaugeFixingMCGIndirect* pFixer, deviceSU2* pData, BYTE byFieldId)
{
    deviceSU2* pG = NULL;
    checkCudaErrors(__cudaMalloc((void**)&pG, _HC_Volume * sizeof(deviceSU2)));

    preparethread;
    pFixer->m_iIterate = 0;
    CGaugeFixingMCGDirect checker;

    while (pFixer->m_iIterate < pFixer->m_iMaxIterate)
    {
        if (0 == pFixer->m_iIterate % pFixer->m_iCheckErrorStep)
        {
            const Real fTheta = static_cast<Real>(checker.CheckResLocalSU2(pData, byFieldId));
            appGeneral(_T("MCG Indirect residual-Abelian (SU2) Iterate : %d, theta = %2.12f\n"), pFixer->m_iIterate, fTheta);
            if (fTheta < pFixer->m_fAccuracy)
            {
                cudaSafeFree(pG);
                return;
            }
        }

        _LAUNCH_KERNEL(_kernelCalculateStandardIMCGSU2, block, threads,
            byFieldId, TRUE, pData, pG,
            pFixer->m_iIMCGGrid, pFixer->m_fOmega);
        _LAUNCH_KERNEL(_kernelGaugeTransformResidualOdd<deviceSU2>, block, threads,
            byFieldId, pG, pData);
        _LAUNCH_KERNEL(_kernelCalculateStandardIMCGSU2, block, threads,
            byFieldId, FALSE, pData, pG,
            pFixer->m_iIMCGGrid, pFixer->m_fOmega);
        _LAUNCH_KERNEL(_kernelGaugeTransformResidualEven<deviceSU2>, block, threads,
            byFieldId, pG, pData);

        ++pFixer->m_iIterate;
    }

    appGeneral(_T("MCG Indirect residual-Abelian (SU2) gauge fixing failed.\n"));
    cudaSafeFree(pG);
}

__CLGIMPLEMENT_CLASS(CGaugeFixingMCGIndirect)

void CGaugeFixingMCGIndirect::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;

    if (!params.FetchValueReal(_T("Omega"), m_fOmega))
    {
        appGeneral(_T("CGaugeFixingMCGIndirect: Omega not set, set to 1.5 by default."));
    }

    if (!params.FetchValueReal(_T("OmegaStage1"), m_fOmegaStage1))
    {
        // Backward compatibility: accept OmegaMAG as alias
        if (!params.FetchValueReal(_T("OmegaMAG"), m_fOmegaStage1))
        {
            // Default: use same omega for both stages
            m_fOmegaStage1 = m_fOmega;
        }
    }

    if (!params.FetchValueReal(_T("Accuracy"), m_fAccuracy))
    {
        appGeneral(_T("CGaugeFixingMCGIndirect: Accuracy not set, set to 0.00000000001 by default."));
        if (m_fAccuracy < _CLG_FLT_EPSILON * F(2.0))
        {
            m_fAccuracy = _CLG_FLT_EPSILON * F(2.0);
            appGeneral(_T("Solver accuracy too small, set to be %2.18f\n"), m_fAccuracy);
        }
    }

    INT iValue = static_cast<INT>(m_iMaxIterate);
    if (!params.FetchValueINT(_T("MaxIterate"), iValue))
    {
        appGeneral(_T("CGaugeFixingMCGIndirect: MaxIterate not set, set to 100000 by default."));
    }
    m_iMaxIterate = static_cast<UINT>(iValue);

    iValue = static_cast<INT>(m_iStage1MaxIterate);
    if (!params.FetchValueINT(_T("Stage1MaxIterate"), iValue))
    {
        // Backward compatibility: accept MAGMaxIterate as alias
        if (!params.FetchValueINT(_T("MAGMaxIterate"), iValue))
        {
            appGeneral(_T("CGaugeFixingMCGIndirect: Stage1MaxIterate not set, set to 100000 by default."));
        }
    }
    m_iStage1MaxIterate = static_cast<UINT>(iValue);

    iValue = static_cast<INT>(m_iCheckErrorStep);
    if (!params.FetchValueINT(_T("CheckErrorStep"), iValue))
    {
        appGeneral(_T("CGaugeFixingMCGIndirect: CheckErrorStep not set, set to 1000 by default."));
    }
    m_iCheckErrorStep = static_cast<UINT>(iValue);

    iValue = m_bUseStandardIMCG ? 1 : 0;
    if (!params.FetchValueINT(_T("UseStandardIMCG"), iValue))
    {
        appGeneral(_T("CGaugeFixingMCGIndirect: UseStandardIMCG not set, using MAG-preconditioned Direct MCG by default."));
    }
    m_bUseStandardIMCG = (0 != iValue);

    iValue = m_iIMCGGrid;
    if (!params.FetchValueINT(_T("IMCGGrid"), iValue))
    {
        appGeneral(_T("CGaugeFixingMCGIndirect: IMCGGrid not set, using default %d.\n"), m_iIMCGGrid);
    }
    m_iIMCGGrid = iValue;
    if (m_iIMCGGrid < 4)
    {
        m_iIMCGGrid = 4;
        appGeneral(_T("CGaugeFixingMCGIndirect: IMCGGrid too small, set to minimum 4.\n"));
    }
}

void CGaugeFixingMCGIndirect::GaugeFixing(CFieldGauge* pResGauge)
{
    if (NULL == pResGauge)
    {
        appCrucial(_T("CGaugeFixingMCGIndirect: pResGauge is NULL!\n"));
        return;
    }

#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        //P4-2.5: gather ONCE to rank 0, then run BOTH stages under the temporary
        //GLOBAL lattice context using the iteration helpers (GaugeFixingLoop)
        //directly. The public GaugeFixing() of the nested fixers must NOT be
        //called here - it would gather again and deadlock on the scatter.
        const EFieldType eType = pResGauge->GetFieldType();
        CFieldGaugeSU2* pGaugeSU2 = (EFT_GaugeSU2 == eType) ? dynamic_cast<CFieldGaugeSU2*>(pResGauge) : NULL;
        CFieldGaugeSU3* pGaugeSU3 = (EFT_GaugeSU3 == eType) ? dynamic_cast<CFieldGaugeSU3*>(pResGauge) : NULL;
        if (NULL == pGaugeSU2 && NULL == pGaugeSU3)
        {
            appCrucial(_T("CGaugeFixingMCGIndirect: unsupported field type %d!\n"), pResGauge->GetFieldType());
            return;
        }

        const UINT uiElemSize = (NULL != pGaugeSU2) ? static_cast<UINT>(sizeof(deviceSU2)) : static_cast<UINT>(sizeof(deviceSU3));
        const UINT uiBytesPerSite = uiElemSize * static_cast<UINT>(_HC_Dir);
        const UINT uiLocalBytes = uiElemSize * static_cast<UINT>(_HC_LinkCount);
        const BYTE* pLocalData = (NULL != pGaugeSU2) ? (const BYTE*)pGaugeSU2->m_pDeviceData : (const BYTE*)pGaugeSU3->m_pDeviceData;

        BYTE* pHostLocal = (BYTE*)malloc(uiLocalBytes);
        checkCudaErrors(cudaMemcpy(pHostLocal, pLocalData, uiLocalBytes, cudaMemcpyDeviceToHost));

        UINT uiGlobalBytes = 0;
        BYTE* pGlobal = appGetComm()->GatherFieldToRoot(pHostLocal, uiBytesPerSite, uiGlobalBytes);
        free(pHostLocal);

        if (appGetComm()->IsRoot())
        {
            BYTE* pDevGlobal = NULL;
            checkCudaErrors(__cudaMalloc((void**)&pDevGlobal, uiGlobalBytes));
            checkCudaErrors(cudaMemcpy(pDevGlobal, pGlobal, uiGlobalBytes, cudaMemcpyHostToDevice));
            free(pGlobal);
            pGlobal = NULL;

            MGEnterGlobalFixerContext();

            //Stage 1: MAG on the global buffer.
            CGaugeFixingMAG magFixer;
            magFixer.m_fOmega = m_fOmegaStage1;
            magFixer.m_fAccuracy = m_fAccuracy;
            magFixer.m_iMaxIterate = m_iStage1MaxIterate;
            magFixer.m_iCheckErrorStep = m_iCheckErrorStep;
            magFixer.m_pOwner = m_pOwner;
            if (NULL != pGaugeSU2)
            {
                magFixer.GaugeFixingLoopSU2(reinterpret_cast<deviceSU2*>(pDevGlobal), pGaugeSU2->m_byFieldId);
            }
            else
            {
                magFixer.GaugeFixingLoopSU3(reinterpret_cast<deviceSU3*>(pDevGlobal), pGaugeSU3->m_byFieldId);
            }

            //Stage 2: standard IMCG or Direct MCG on the global buffer.
            UINT uiStage2Iterate = 0;
            if (m_bUseStandardIMCG)
            {
                if (NULL != pGaugeSU2)
                {
                    _StandardIMCGGaugeFixingSU2(this, reinterpret_cast<deviceSU2*>(pDevGlobal), pGaugeSU2->m_byFieldId);
                }
                else
                {
                    _StandardIMCGGaugeFixingSU3(this, reinterpret_cast<deviceSU3*>(pDevGlobal), pGaugeSU3->m_byFieldId);
                }
                uiStage2Iterate = m_iIterate;
            }
            else
            {
                CGaugeFixingMCGDirect mcgFixer;
                mcgFixer.m_fOmega = m_fOmega;
                mcgFixer.m_fAccuracy = m_fAccuracy;
                mcgFixer.m_iMaxIterate = m_iMaxIterate;
                mcgFixer.m_iCheckErrorStep = m_iCheckErrorStep;
                mcgFixer.m_pOwner = m_pOwner;
                if (NULL != pGaugeSU2)
                {
                    mcgFixer.GaugeFixingLoopSU2(reinterpret_cast<deviceSU2*>(pDevGlobal), pGaugeSU2->m_byFieldId);
                }
                else
                {
                    mcgFixer.GaugeFixingLoopSU3(reinterpret_cast<deviceSU3*>(pDevGlobal), pGaugeSU3->m_byFieldId);
                }
                uiStage2Iterate = mcgFixer.m_iIterate;
            }
            m_iIterate = magFixer.m_iIterate + uiStage2Iterate;

            MGExitGlobalFixerContext();

            pGlobal = (BYTE*)malloc(uiGlobalBytes);
            checkCudaErrors(cudaMemcpy(pGlobal, pDevGlobal, uiGlobalBytes, cudaMemcpyDeviceToHost));
            checkCudaErrors(__cudaFree(pDevGlobal));
        }

        BYTE* pLocalOut = (BYTE*)malloc(uiLocalBytes);
        appGetComm()->ScatterFieldFromRoot(pGlobal, uiBytesPerSite, pLocalOut);
        checkCudaErrors(cudaMemcpy((void*)(NULL != pGaugeSU2 ? (void*)pGaugeSU2->m_pDeviceData : (void*)pGaugeSU3->m_pDeviceData), pLocalOut, uiLocalBytes, cudaMemcpyHostToDevice));
        free(pLocalOut);
        if (NULL != pGlobal)
        {
            free(pGlobal);
        }

        //Every local link was rewritten by the scatter memcpy, which bypasses
        //guarded launches; bump the owner handle so the next reader's guard
        //refills the halo (no-op for an unbound handle).
        pResGauge->NotifyWritten();
        return;
    }
#endif

    switch (pResGauge->GetFieldType())
    {
    case EFT_GaugeSU2:
    case EFT_GaugeSU3:
        break;
    default:
        appCrucial(_T("CGaugeFixingMCGIndirect: unsupported field type %d!\n"), pResGauge->GetFieldType());
        return;
    }

    // ================================================================
    // Stage 1: Maximal Abelian Gauge fixing
    // ================================================================
    appGeneral(_T("MCG Indirect: Starting Stage 1 (MAG)...\n"));

    CGaugeFixingMAG magFixer;
    magFixer.m_fOmega = m_fOmegaStage1;
    magFixer.m_fAccuracy = m_fAccuracy;
    magFixer.m_iMaxIterate = m_iStage1MaxIterate;
    magFixer.m_iCheckErrorStep = m_iCheckErrorStep;
    magFixer.m_pOwner = m_pOwner;

    appGeneral(_T("  MAG Omega=%f, Accuracy=%f, MaxIterate=%d\n"),
        magFixer.m_fOmega, magFixer.m_fAccuracy, magFixer.m_iMaxIterate);

    magFixer.GaugeFixing(pResGauge);

    appGeneral(_T("MCG Indirect: Stage 1 (MAG) completed after %d iterations.\n"), magFixer.m_iIterate);

    // ================================================================
    // Stage 2: MCG fixing on the MAG-fixed field
    // ================================================================
    appGeneral(m_bUseStandardIMCG
        ? _T("MCG Indirect: Starting Stage 2 (residual-Abelian standard IMCG)...\n")
        : _T("MCG Indirect: Starting Stage 2 (MAG-preconditioned Direct MCG)...\n"));

    UINT uiStage2Iterate = 0;
    if (m_bUseStandardIMCG)
    {
        if (EFT_GaugeSU2 == pResGauge->GetFieldType())
        {
            _StandardIMCGGaugeFixingSU2(this, dynamic_cast<CFieldGaugeSU2*>(pResGauge)->m_pDeviceData, pResGauge->m_byFieldId);
        }
        else
        {
            _StandardIMCGGaugeFixingSU3(this, dynamic_cast<CFieldGaugeSU3*>(pResGauge)->m_pDeviceData, pResGauge->m_byFieldId);
        }
        uiStage2Iterate = m_iIterate;
    }
    else
    {
        CGaugeFixingMCGDirect mcgFixer;
        mcgFixer.m_fOmega = m_fOmega;
        mcgFixer.m_fAccuracy = m_fAccuracy;
        mcgFixer.m_iMaxIterate = m_iMaxIterate;
        mcgFixer.m_iCheckErrorStep = m_iCheckErrorStep;
        mcgFixer.m_pOwner = m_pOwner;

        appGeneral(_T("  MCG Omega=%f, Accuracy=%f, MaxIterate=%d\n"),
            mcgFixer.m_fOmega, mcgFixer.m_fAccuracy, mcgFixer.m_iMaxIterate);

        mcgFixer.GaugeFixing(pResGauge);
        uiStage2Iterate = mcgFixer.m_iIterate;
    }

    appGeneral(_T("MCG Indirect: Stage 2 (MCG) completed after %d iterations.\n"), uiStage2Iterate);

    // Track total iterations for reporting
    m_iIterate = magFixer.m_iIterate + uiStage2Iterate;
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CGaugeFixingMCGIndirect::CheckRes(const CFieldGauge* pGauge)
#else
Real CGaugeFixingMCGIndirect::CheckRes(const CFieldGauge* pGauge)
#endif
{
    if (NULL == pGauge)
    {
        appCrucial(_T("CGaugeFixingMCGIndirect::CheckRes: pGauge is NULL!\n"));
#if !_CLG_DOUBLEFLOAT
        return 0.0;
#else
        return F(0.0);
#endif
    }

    switch (pGauge->GetFieldType())
    {
    case EFT_GaugeSU2:
    case EFT_GaugeSU3:
        break;
    default:
        appCrucial(_T("CGaugeFixingMCGIndirect::CheckRes: unsupported field type %d!\n"), pGauge->GetFieldType());
#if !_CLG_DOUBLEFLOAT
        return 0.0;
#else
        return F(0.0);
#endif
    }

    // Return the MCG objective (Stage 2 residual).
    // After both MAG and MCG, the field should be well-centered.
    CGaugeFixingMCGDirect mcgFixer;
    return mcgFixer.CheckRes(pGauge);
}

CCString CGaugeFixingMCGIndirect::GetInfos(const CCString& tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeFixingMCGIndirect\n");
    sRet = sRet + tab + _T("Omega (Stage 2 MCG) : ") + appToString(m_fOmega) + _T("\n");
    sRet = sRet + tab + _T("OmegaStage1 (Stage 1 MAG) : ") + appToString(m_fOmegaStage1) + _T("\n");
    sRet = sRet + tab + _T("Stage1 MaxIterate : ") + appToString(m_iStage1MaxIterate) + _T("\n");
    sRet = sRet + tab + _T("UseStandardIMCG : ") + appToString(m_bUseStandardIMCG) + _T("\n");
    sRet = sRet + tab + _T("Total Iterations : ") + appToString(m_iIterate) + _T("\n");
    return sRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
