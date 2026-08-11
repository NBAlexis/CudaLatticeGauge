//=============================================================================
// FILENAME : CGaugeFixingMAG.cu
//
// DESCRIPTION:
// Maximal Abelian Gauge (MAG) fixing for SU(2) and SU(3).
//
// SU(2) MAG algorithm (Cea & Cosmai, hep-lat/9504008, Appendix A):
// The MAG maximizes the functional
//   G_MAG = sum_{x,mu} Tr(U_mu^dag sigma_3 U_mu sigma_3)
// which drives each link toward diagonal form.
//
// Local update at site x constructs the traceless Hermitian matrix:
//   X(x) = sum_mu [ U_mu(x) sigma_3 U_mu^dag(x)
//                 + U_mu^dag(x-mu) sigma_3 U_mu(x-mu) ]
//
// Then V(x) = X(x) sigma_3 / k(x), where k(x) = sqrt(det(X sigma_3)) and V in SU(2).
// Writing V = v0 I + i(v1 sigma_1 + v2 sigma_2), the optimal gauge transformation
// with overrelaxation parameter omega is:
//   g^omega = cos(omega * alpha) I
//           - i(v1 sigma_1 + v2 sigma_2) / sqrt(1-v0^2) * sin(omega * alpha)
// where cos(2*alpha) = v0, i.e. alpha = acos(v0) / 2.
//
// When v0 is numerically >= 1 (already aligned), g = I.
//
// Convergence criterion:
//   theta = sum_{x,mu} |U_{01}(x,mu)|^2 / (V * D)
//   Stop when theta < accuracy.
//
// SU(3) MAG algorithm (Tucker & Stack, hep-lat/0110165):
// Uses Cabibbo-Marinari on SU(2) subgroups.
//
// Maximal Abelian Projection:
//   SU(2): keep diagonal phase, zero off-diagonal.
//   SU(3): keep diagonal elements, zero off-diagonal, adjust phase for det=1.
//
// References:
//  hep-lat/9504008 -- Cea & Cosmai, "Maximal Abelian Gauge in SU(2)"
//  hep-lat/0110165 -- Tucker & Stack, "The Maximal Abelian Gauge in SU(3)"
//
// REVISION:
//  [05/14/2026 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CGaugeFixingMAG.h"

__BEGIN_NAMESPACE

#pragma region SU2 kernels

// ============================================================================
// Compute SU(2) MAG gauge transformation for odd sites.
//
// Algorithm (Cea & Cosmai, hep-lat/9504008, Appendix A):
//   X(x) = sum_mu [ U_mu(x) sigma_3 U_mu^dag(x)
//                 + U_mu^dag(x-mu) sigma_3 U_mu(x-mu) ]
//
// For SU(2) link U = [a, b; -b*, a*]:
//   U sigma_3 U^dag = [|a|^2-|b|^2, -2ab; -2a*b*, |b|^2-|a|^2]
//   U^dag sigma_3 U  = [|a|^2-|b|^2,  2a*b;  2ab*, |b|^2-|a|^2]
//
// X is Hermitian and traceless: X = [X3, X01; X01*, -X3]
// V = X sigma_3 / k, k = sqrt(X3^2 + |X01|^2)
// V = [v0, v2+iv1; -v2+iv1, v0] with v0 = X3/k, v1 = Im(V01), v2 = Re(V01)
//
// Overrelaxation:
//   g^omega = cos(omega*acos(v0)) I
//           + i(v1 sigma_1 + v2 sigma_2)/sqrt(1-v0^2) * sin(omega*acos(v0))
// When v0 >= 1 (numerically), g = I.
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGMAGSU2Odd(
    BYTE byFieldId,
    const deviceSU2* __restrict__ pU,
    Real fOmega,
    deviceSU2* pG)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (site.IsDirichlet() || !sSite4.IsOdd())
    {
        return;
    }

    // X = X3 sigma_3 + X1 sigma_1 + X2 sigma_2 = [X3, X01; X01*, -X3]
    // We accumulate X3 (real), X01 (complex)
    Real fX3 = F(0.0);
    CLGComplex cX01 = _make_cuComplex(F(0.0), F(0.0));

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        // Forward link: U_mu(x) sigma_3 U_mu^dag(x)
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const deviceSU2& Ufwd = pU[uiLinkIndex];
            // U = [a, b; -b*, a*]
            // U sigma_3 U^dag = [|a|^2-|b|^2, -2ab; -2a*b*, |b|^2-|a|^2]
            CLGComplex a = Ufwd.m_me[0];
            CLGComplex b = Ufwd.m_me[1];
            Real absa2 = __cuCabsSqf(a);
            Real absb2 = __cuCabsSqf(b);
            fX3 += (absa2 - absb2);
            CLGComplex minus2ab = _make_cuComplex(
                -F(2.0) * (a.x * b.x - a.y * b.y),
                -F(2.0) * (a.x * b.y + a.y * b.x));
            cX01 = _cuCaddf(cX01, minus2ab);
        }
        else
        {
            // Boundary: U = I, U sigma_3 U^dag = sigma_3 = [1, 0; 0, -1]
            fX3 += F(1.0);
        }

        // Backward link: U_mu^dag(x-mu) sigma_3 U_mu(x-mu)
        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];
        const deviceSU2 Ubwd = _deviceGetGaugeBCDirOneSIndexT<deviceSU2>(pU, site_m_mu);
        // U^dag sigma_3 U = [|a|^2-|b|^2, 2a*b; 2ab*, |b|^2-|a|^2]
        CLGComplex a = Ubwd.m_me[0];
        CLGComplex b = Ubwd.m_me[1];
        Real absa2 = __cuCabsSqf(a);
        Real absb2 = __cuCabsSqf(b);
        fX3 += (absa2 - absb2);
        CLGComplex twoaStarB = _make_cuComplex(
            F(2.0) * (a.x * b.x + a.y * b.y),
            F(2.0) * (a.x * b.y - a.y * b.x));
        cX01 = _cuCaddf(cX01, twoaStarB);
    }

    // X = [fX3, cX01; cX01*, -fX3]
    // X sigma_3 = [fX3, -cX01; cX01*, fX3]
    // det(X sigma_3) = fX3^2 + |cX01|^2
    Real kSq = fX3 * fX3 + __cuCabsSqf(cX01);
    Real k = _sqrt(kSq);

    deviceSU2 su2Id = deviceSU2::makeSU2Id();
    if (k < _CLG_FLT_EPSILON)
    {
        // X is zero, no update needed
        pG[uiSiteIndex] = su2Id;
        return;
    }

    // V = X sigma_3 / k = [v0, v2+iv1; -v2+iv1, v0]
    Real v0 = fX3 / k;
    CLGComplex V01 = _make_cuComplex(-cX01.x / k, -cX01.y / k);
    Real v1 = V01.y;  // Im(V01)
    Real v2 = V01.x;  // Re(V01)

    // Clamp v0 to [-1, 1] to avoid numerical issues
    if (v0 > F(1.0)) v0 = F(1.0);
    if (v0 < F(-1.0)) v0 = F(-1.0);

    // Overrelaxation: g^omega
    // For SU(2), g = cos(alpha) I - i sin(alpha) (n·sigma) rotates vectors by 2*alpha.
    // We need cos(2*alpha) = v0, so alpha = acos(v0) / 2.
    // g^omega = cos(omega*alpha) I
    //         - i(v1 sigma_1 + v2 sigma_2)/sqrt(1-v0^2) * sin(omega*alpha)
    // When v0 is numerically >= 1 (already aligned), g = I.
    Real norm = _sqrt(F(1.0) - v0 * v0);

    if (norm < _CLG_FLT_EPSILON)
    {
        // Already aligned
        pG[uiSiteIndex] = su2Id;
        return;
    }

    Real alpha = acos(v0) / F(2.0);
    Real omegaAlpha = fOmega * alpha;
    Real cosOmegaAlpha = _cos(omegaAlpha);
    Real sinOmegaAlpha = _sin(omegaAlpha);

    // -i(v1 sigma_1 + v2 sigma_2) = [0, -v2-iv1; v2-iv1, 0]
    // Normalized: [0, (-v2-iv1)/norm; (v2-iv1)/norm, 0]
    Real g01x = (-v2 / norm) * sinOmegaAlpha;
    Real g01y = (-v1 / norm) * sinOmegaAlpha;
    Real g10x = (v2 / norm) * sinOmegaAlpha;
    Real g10y = (-v1 / norm) * sinOmegaAlpha;

    deviceSU2 g;
    g.m_me[0] = _make_cuComplex(cosOmegaAlpha, F(0.0));
    g.m_me[1] = _make_cuComplex(g01x, g01y);
    g.m_me[2] = _make_cuComplex(g10x, g10y);
    g.m_me[3] = _make_cuComplex(cosOmegaAlpha, F(0.0));

    pG[uiSiteIndex] = g;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGMAGSU2Even(
    BYTE byFieldId,
    const deviceSU2* __restrict__ pU,
    Real fOmega,
    deviceSU2* pG)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (site.IsDirichlet() || sSite4.IsOdd())
    {
        return;
    }

    Real fX3 = F(0.0);
    CLGComplex cX01 = _make_cuComplex(F(0.0), F(0.0));

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const deviceSU2& Ufwd = pU[uiLinkIndex];
            CLGComplex a = Ufwd.m_me[0];
            CLGComplex b = Ufwd.m_me[1];
            Real absa2 = __cuCabsSqf(a);
            Real absb2 = __cuCabsSqf(b);
            fX3 += (absa2 - absb2);
            CLGComplex minus2ab = _make_cuComplex(
                -F(2.0) * (a.x * b.x - a.y * b.y),
                -F(2.0) * (a.x * b.y + a.y * b.x));
            cX01 = _cuCaddf(cX01, minus2ab);
        }
        else
        {
            fX3 += F(1.0);
        }

        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];
        const deviceSU2 Ubwd = _deviceGetGaugeBCDirOneSIndexT<deviceSU2>(pU, site_m_mu);
        CLGComplex a = Ubwd.m_me[0];
        CLGComplex b = Ubwd.m_me[1];
        Real absa2 = __cuCabsSqf(a);
        Real absb2 = __cuCabsSqf(b);
        fX3 += (absa2 - absb2);
        CLGComplex twoaStarB = _make_cuComplex(
            F(2.0) * (a.x * b.x + a.y * b.y),
            F(2.0) * (a.x * b.y - a.y * b.x));
        cX01 = _cuCaddf(cX01, twoaStarB);
    }

    Real kSq = fX3 * fX3 + __cuCabsSqf(cX01);
    Real k = _sqrt(kSq);

    deviceSU2 su2Id = deviceSU2::makeSU2Id();
    if (k < _CLG_FLT_EPSILON)
    {
        pG[uiSiteIndex] = su2Id;
        return;
    }

    Real v0 = fX3 / k;
    CLGComplex V01 = _make_cuComplex(-cX01.x / k, -cX01.y / k);
    Real v1 = V01.y;
    Real v2 = V01.x;

    if (v0 > F(1.0)) v0 = F(1.0);
    if (v0 < F(-1.0)) v0 = F(-1.0);

    // Overrelaxation: g^omega
    // For SU(2), g = cos(alpha) I - i sin(alpha) (n·sigma) rotates vectors by 2*alpha.
    // We need cos(2*alpha) = v0, so alpha = acos(v0) / 2.
    Real norm = _sqrt(F(1.0) - v0 * v0);

    if (norm < _CLG_FLT_EPSILON)
    {
        pG[uiSiteIndex] = su2Id;
        return;
    }

    Real alpha = acos(v0) / F(2.0);
    Real omegaAlpha = fOmega * alpha;
    Real cosOmegaAlpha = _cos(omegaAlpha);
    Real sinOmegaAlpha = _sin(omegaAlpha);

    // -i(v1 sigma_1 + v2 sigma_2) = [0, -v2-iv1; v2-iv1, 0]
    Real g01x = (-v2 / norm) * sinOmegaAlpha;
    Real g01y = (-v1 / norm) * sinOmegaAlpha;
    Real g10x = (v2 / norm) * sinOmegaAlpha;
    Real g10y = (-v1 / norm) * sinOmegaAlpha;

    deviceSU2 g;
    g.m_me[0] = _make_cuComplex(cosOmegaAlpha, F(0.0));
    g.m_me[1] = _make_cuComplex(g01x, g01y);
    g.m_me[2] = _make_cuComplex(g10x, g10y);
    g.m_me[3] = _make_cuComplex(cosOmegaAlpha, F(0.0));

    pG[uiSiteIndex] = g;
}

// ============================================================================
// Gauge transform kernels for SU(2) (red-black sweep).
// Same pattern as SU(3).
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformMAGSU2Odd(
    BYTE byFieldId,
    const deviceSU2* __restrict__ pGx,
    deviceSU2* pGauge)
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

__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformMAGSU2Even(
    BYTE byFieldId,
    const deviceSU2* __restrict__ pGx,
    deviceSU2* pGauge)
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

// ============================================================================
// MAG convergence check for SU(2): sum of squared off-diagonal elements.
//
// For SU(2) U = [a, b; -b*, a*], the off-diagonal element is b.
// MAG is maximized when |b|^2 = 0 (fully diagonal).
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateMAGObjectiveSU2(
    BYTE byFieldId,
    DOUBLE* pDeviceRes,
    const deviceSU2* __restrict__ pU)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (site.IsDirichlet())
    {
        pDeviceRes[uiSiteIndex] = 0.0;
        return;
    }

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    DOUBLE objective = 0.0;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const deviceSU2& U = pU[uiLinkIndex];
            // Off-diagonal element is U_{01} = b
            objective += static_cast<DOUBLE>(__cuCabsSqf(U.m_me[1]));
        }
    }

    pDeviceRes[uiSiteIndex] = objective;
}

// ============================================================================
// Maximal Abelian Projection for SU(2): project each link to U(1) subgroup.
//
// For SU(2) U = [a, b; -b*, a*], keep the diagonal phase:
//   U_proj = [a/|a|, 0; 0, a*/|a|]
// This preserves det = 1 since |a/|a|| = 1.
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelMaximalAbelianProjectionSU2(
    BYTE byFieldId,
    deviceSU2* pGauge)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (site.IsDirichlet())
    {
        return;
    }

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            deviceSU2& U = pGauge[uiLinkIndex];
            CLGComplex a = U.m_me[0];
            Real absa = _sqrt(__cuCabsSqf(a));

            if (absa > _CLG_FLT_EPSILON)
            {
                CLGComplex phase = _make_cuComplex(a.x / absa, a.y / absa);
                CLGComplex phaseConj = _make_cuComplex(a.x / absa, -a.y / absa);
                U = deviceSU2::makeSU2Zero();
                U.m_me[0] = phase;
                U.m_me[3] = phaseConj;
            }
            // If absa == 0, the matrix is anti-diagonal; leave as is or set to I
        }
    }
}

// ============================================================================
// Device helper: solve the SU(2) MAG local problem and return the gauge
// transformation g in SU(2).
//
// X = X3 sigma_3 + X1 sigma_1 + X2 sigma_2 = [X3, X01; X01*, -X3]
// Same logic as _kernelCalculateGMAGSU2Odd/Even.
// ============================================================================
__device__ __inline__ deviceSU2 _deviceSolveSU2MAG(
    Real fX3,
    CLGComplex cX01,
    Real fOmega)
{
    deviceSU2 su2Id = deviceSU2::makeSU2Id();

    // X sigma_3 = [fX3, -cX01; cX01*, fX3]
    Real kSq = fX3 * fX3 + __cuCabsSqf(cX01);
    Real k = _sqrt(kSq);

    if (k < _CLG_FLT_EPSILON)
    {
        return su2Id;
    }

    // V = X sigma_3 / k = [v0, v2+iv1; -v2+iv1, v0]
    Real v0 = fX3 / k;
    CLGComplex V01 = _make_cuComplex(-cX01.x / k, -cX01.y / k);
    Real v1 = V01.y;
    Real v2 = V01.x;

    if (v0 > F(1.0)) v0 = F(1.0);
    if (v0 < F(-1.0)) v0 = F(-1.0);

    Real norm = _sqrt(F(1.0) - v0 * v0);

    if (norm < _CLG_FLT_EPSILON)
    {
        return su2Id;
    }

    // For SU(2), g = cos(alpha) I - i sin(alpha) (n·sigma) rotates by 2*alpha.
    // Need cos(2*alpha) = v0, so alpha = acos(v0) / 2.
    Real alpha = acos(v0) / F(2.0);
    Real omegaAlpha = fOmega * alpha;
    Real cosOmegaAlpha = _cos(omegaAlpha);
    Real sinOmegaAlpha = _sin(omegaAlpha);

    Real g01x = (-v2 / norm) * sinOmegaAlpha;
    Real g01y = (-v1 / norm) * sinOmegaAlpha;
    Real g10x = (v2 / norm) * sinOmegaAlpha;
    Real g10y = (-v1 / norm) * sinOmegaAlpha;

    deviceSU2 g;
    g.m_me[0] = _make_cuComplex(cosOmegaAlpha, F(0.0));
    g.m_me[1] = _make_cuComplex(g01x, g01y);
    g.m_me[2] = _make_cuComplex(g10x, g10y);
    g.m_me[3] = _make_cuComplex(cosOmegaAlpha, F(0.0));

    return g;
}

#pragma endregion

#pragma region SU3 kernels

// ============================================================================
// SU(3) MAG via Cabibbo-Marinari on SU(2) subgroups.
//
// For each of the three subgroups R12=(0,1), R23=(1,2), R13=(0,2):
//   1. Extract the 2x2 submatrix W from each link.
//   2. Build the SU(2) MAG environment:
//        X = sum_mu [ W_mu(x) sigma_3 W_mu(x)^dag + W_mu(x-mu)^dag sigma_3 W_mu(x-mu) ]
//   3. Solve for the optimal g in SU(2) (same as SU2 MAG local solve).
//   4. Embed g back into SU(3) as a rotation in the subgroup.
//
// Red-black checkerboard sweep, sequential over subgroups.
// ============================================================================

// ----------------------------------------------------------------------------
// Extract the 2x2 submatrix for a given subgroup.
//   R12: indices (0,1)
//   R23: indices (1,2)
//   R13: indices (0,2)
// ----------------------------------------------------------------------------
__device__ __inline__ deviceSU2 _deviceExtractSU2Subgroup(
    const deviceSU3& U,
    BYTE bySubgroup)
{
    deviceSU2 ret;
    switch (bySubgroup)
    {
    case 0: // R12: (0,1)
        ret.m_me[0] = U.m_me[0];  // U00
        ret.m_me[1] = U.m_me[1];  // U01
        ret.m_me[2] = U.m_me[3];  // U10
        ret.m_me[3] = U.m_me[4];  // U11
        break;
    case 1: // R23: (1,2)
        ret.m_me[0] = U.m_me[4];  // U11
        ret.m_me[1] = U.m_me[5];  // U12
        ret.m_me[2] = U.m_me[7];  // U21
        ret.m_me[3] = U.m_me[8];  // U22
        break;
    case 2: // R13: (0,2)
        ret.m_me[0] = U.m_me[0];  // U00
        ret.m_me[1] = U.m_me[2];  // U02
        ret.m_me[2] = U.m_me[6];  // U20
        ret.m_me[3] = U.m_me[8];  // U22
        break;
    }
    return ret;
}

// ----------------------------------------------------------------------------
// Embed an SU(2) matrix into SU(3) for a given subgroup.
//   R12: diag(g, 1)
//   R23: diag(1, g)
//   R13: [g00, 0, g01; 0, 1, 0; g10, 0, g11]
// ----------------------------------------------------------------------------
__device__ __inline__ deviceSU3 _deviceEmbedSU2ToSU3(
    const deviceSU2& g,
    BYTE bySubgroup)
{
    deviceSU3 ret = deviceSU3::makeSU3Id();
    switch (bySubgroup)
    {
    case 0: // R12: (0,1)
        ret.m_me[0] = g.m_me[0];
        ret.m_me[1] = g.m_me[1];
        ret.m_me[3] = g.m_me[2];
        ret.m_me[4] = g.m_me[3];
        break;
    case 1: // R23: (1,2)
        ret.m_me[4] = g.m_me[0];
        ret.m_me[5] = g.m_me[1];
        ret.m_me[7] = g.m_me[2];
        ret.m_me[8] = g.m_me[3];
        break;
    case 2: // R13: (0,2)
        ret.m_me[0] = g.m_me[0];
        ret.m_me[2] = g.m_me[1];
        ret.m_me[6] = g.m_me[2];
        ret.m_me[8] = g.m_me[3];
        break;
    }
    return ret;
}

// ----------------------------------------------------------------------------
// Compute SU(3) MAG gauge transform for one subgroup on odd sites.
// ----------------------------------------------------------------------------
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGMAGSU3SubgroupOdd(
    BYTE byFieldId,
    BYTE bySubgroup,
    const deviceSU3* __restrict__ pU,
    Real fOmega,
    deviceSU3* pG)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (site.IsDirichlet() || !sSite4.IsOdd())
    {
        return;
    }

    // X = X3 sigma_3 + X1 sigma_1 + X2 sigma_2 = [X3, X01; X01*, -X3]
    Real fX3 = F(0.0);
    CLGComplex cX01 = _make_cuComplex(F(0.0), F(0.0));

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        // Forward link: W sigma_3 W^dag
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            deviceSU2 W = _deviceExtractSU2Subgroup(pU[uiLinkIndex], bySubgroup);
            // W is a general SU(3) principal 2x2 block, not necessarily an SU(2)
            // matrix.  For W = [a,b;c,d], keep the traceless part of
            // W sigma_3 W^dag:
            //   X3  = 1/2[(|a|^2-|b|^2) - (|c|^2-|d|^2)]
            //   X01 = a c* - b d*
            CLGComplex a = W.m_me[0];
            CLGComplex b = W.m_me[1];
            CLGComplex c = W.m_me[2];
            CLGComplex d = W.m_me[3];
            Real absa2 = __cuCabsSqf(a);
            Real absb2 = __cuCabsSqf(b);
            Real absc2 = __cuCabsSqf(c);
            Real absd2 = __cuCabsSqf(d);
            fX3 += F(0.5) * ((absa2 - absb2) - (absc2 - absd2));
            CLGComplex x01 = _cuCsubf(_cuCmulf(a, _cuConjf(c)), _cuCmulf(b, _cuConjf(d)));
            cX01 = _cuCaddf(cX01, x01);
        }
        else
        {
            // Boundary: W = I, contribution is sigma_3
            fX3 += F(1.0);
        }

        // Backward link: W^dag sigma_3 W
        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];
        const deviceSU3 Ubwd3 = _deviceGetGaugeBCSU3DirOneSIndex(pU, site_m_mu);
        deviceSU2 Wbwd = _deviceExtractSU2Subgroup(Ubwd3, bySubgroup);
        // W^dag sigma_3 W for a general W = [a,b;c,d]:
        //   X3  = 1/2[(|a|^2-|c|^2) - (|b|^2-|d|^2)]
        //   X01 = a* b - c* d
        CLGComplex a = Wbwd.m_me[0];
        CLGComplex b = Wbwd.m_me[1];
        CLGComplex c = Wbwd.m_me[2];
        CLGComplex d = Wbwd.m_me[3];
        Real absa2 = __cuCabsSqf(a);
        Real absb2 = __cuCabsSqf(b);
        Real absc2 = __cuCabsSqf(c);
        Real absd2 = __cuCabsSqf(d);
        fX3 += F(0.5) * ((absa2 - absc2) - (absb2 - absd2));
        CLGComplex x01 = _cuCsubf(_cuCmulf(_cuConjf(a), b), _cuCmulf(_cuConjf(c), d));
        cX01 = _cuCaddf(cX01, x01);
    }

    deviceSU2 gSU2 = _deviceSolveSU2MAG(fX3, cX01, fOmega);
    pG[uiSiteIndex] = _deviceEmbedSU2ToSU3(gSU2, bySubgroup);
}

// ----------------------------------------------------------------------------
// Compute SU(3) MAG gauge transform for one subgroup on even sites.
// ----------------------------------------------------------------------------
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGMAGSU3SubgroupEven(
    BYTE byFieldId,
    BYTE bySubgroup,
    const deviceSU3* __restrict__ pU,
    Real fOmega,
    deviceSU3* pG)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (site.IsDirichlet() || sSite4.IsOdd())
    {
        return;
    }

    Real fX3 = F(0.0);
    CLGComplex cX01 = _make_cuComplex(F(0.0), F(0.0));

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            deviceSU2 W = _deviceExtractSU2Subgroup(pU[uiLinkIndex], bySubgroup);
            CLGComplex a = W.m_me[0];
            CLGComplex b = W.m_me[1];
            CLGComplex c = W.m_me[2];
            CLGComplex d = W.m_me[3];
            Real absa2 = __cuCabsSqf(a);
            Real absb2 = __cuCabsSqf(b);
            Real absc2 = __cuCabsSqf(c);
            Real absd2 = __cuCabsSqf(d);
            fX3 += F(0.5) * ((absa2 - absb2) - (absc2 - absd2));
            CLGComplex x01 = _cuCsubf(_cuCmulf(a, _cuConjf(c)), _cuCmulf(b, _cuConjf(d)));
            cX01 = _cuCaddf(cX01, x01);
        }
        else
        {
            fX3 += F(1.0);
        }

        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];
        const deviceSU3 Ubwd3 = _deviceGetGaugeBCSU3DirOneSIndex(pU, site_m_mu);
        deviceSU2 Wbwd = _deviceExtractSU2Subgroup(Ubwd3, bySubgroup);
        CLGComplex a = Wbwd.m_me[0];
        CLGComplex b = Wbwd.m_me[1];
        CLGComplex c = Wbwd.m_me[2];
        CLGComplex d = Wbwd.m_me[3];
        Real absa2 = __cuCabsSqf(a);
        Real absb2 = __cuCabsSqf(b);
        Real absc2 = __cuCabsSqf(c);
        Real absd2 = __cuCabsSqf(d);
        fX3 += F(0.5) * ((absa2 - absc2) - (absb2 - absd2));
        CLGComplex x01 = _cuCsubf(_cuCmulf(_cuConjf(a), b), _cuCmulf(_cuConjf(c), d));
        cX01 = _cuCaddf(cX01, x01);
    }

    deviceSU2 gSU2 = _deviceSolveSU2MAG(fX3, cX01, fOmega);
    pG[uiSiteIndex] = _deviceEmbedSU2ToSU3(gSU2, bySubgroup);
}

// ============================================================================
// Gauge transform kernels (red-black sweep).
// Same as MCG.
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformMAGOdd(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
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

__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformMAGEven(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
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

// ============================================================================
// MAG convergence check: compute sum of squared off-diagonal elements.
//
// MAG functional is maximized when off-diagonal elements are zero.
// Convergence measure:
//   theta = sum_{x,\mu} sum_{i \neq j} |U_{ij}|^2 / (V * D)
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateMAGObjective(
    BYTE byFieldId,
    DOUBLE* pDeviceRes,
    const deviceSU3* __restrict__ pU)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (site.IsDirichlet())
    {
        pDeviceRes[uiSiteIndex] = 0.0;
        return;
    }

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    DOUBLE objective = 0.0;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const deviceSU3& U = pU[uiLinkIndex];
            // Off-diagonal elements: 01, 02, 10, 12, 20, 21
            objective += static_cast<DOUBLE>(__cuCabsSqf(U.m_me[1]));
            objective += static_cast<DOUBLE>(__cuCabsSqf(U.m_me[2]));
            objective += static_cast<DOUBLE>(__cuCabsSqf(U.m_me[3]));
            objective += static_cast<DOUBLE>(__cuCabsSqf(U.m_me[5]));
            objective += static_cast<DOUBLE>(__cuCabsSqf(U.m_me[6]));
            objective += static_cast<DOUBLE>(__cuCabsSqf(U.m_me[7]));
        }
    }

    pDeviceRes[uiSiteIndex] = objective;
}

// ============================================================================
// Maximal Abelian Projection: project each SU(3) link to U(1)xU(1) subgroup.
//
// Algorithm (Tucker & Stack, hep-lat/0110165):
//   For each link U, keep diagonal elements, zero off-diagonal elements,
//   normalize each diagonal to unit magnitude, then adjust overall phase
//   to enforce det = 1.
//
// Practical steps:
//   1. d_i = U_{ii}
//   2. d_i' = d_i / |d_i|  (normalize to unit magnitude)
//   3. phi = arg(d_1') + arg(d_2') + arg(d_3')
//   4. U_abelian = diag(d_1', d_2', d_3') * exp(-i*phi/3)
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelMaximalAbelianProjectionSU3(
    BYTE byFieldId,
    deviceSU3* pGauge)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (site.IsDirichlet())
    {
        return;
    }

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            deviceSU3& U = pGauge[uiLinkIndex];

            // Extract diagonal elements and normalize to unit magnitude
            CLGComplex d1 = U.m_me[0];
            CLGComplex d2 = U.m_me[4];
            CLGComplex d3 = U.m_me[8];

            Real abs1 = _sqrt(__cuCabsSqf(d1));
            Real abs2 = _sqrt(__cuCabsSqf(d2));
            Real abs3 = _sqrt(__cuCabsSqf(d3));

            if (abs1 > _CLG_FLT_EPSILON) d1 = _make_cuComplex(d1.x / abs1, d1.y / abs1);
            else d1 = _make_cuComplex(F(1.0), F(0.0));

            if (abs2 > _CLG_FLT_EPSILON) d2 = _make_cuComplex(d2.x / abs2, d2.y / abs2);
            else d2 = _make_cuComplex(F(1.0), F(0.0));

            if (abs3 > _CLG_FLT_EPSILON) d3 = _make_cuComplex(d3.x / abs3, d3.y / abs3);
            else d3 = _make_cuComplex(F(1.0), F(0.0));

            // Adjust overall phase so that det = 1
            CLGComplex det = _cuCmulf(_cuCmulf(d1, d2), d3);
            Real phi = atan2f(det.y, det.x);

            Real phaseReal = cosf(-phi / F(3.0));
            Real phaseImag = sinf(-phi / F(3.0));
            CLGComplex phase = _make_cuComplex(phaseReal, phaseImag);

            // Project to diagonal with normalized and phase-adjusted elements
            U = deviceSU3::makeSU3Zero();
            U.m_me[0] = _cuCmulf(phase, d1);
            U.m_me[4] = _cuCmulf(phase, d2);
            U.m_me[8] = _cuCmulf(phase, d3);
        }
    }
}

#pragma endregion

// ============================================================================
// Host helper: runs the full MAG iteration for SU(2).
// ============================================================================
void CGaugeFixingMAG::GaugeFixingLoopSU2(
    deviceSU2* pDeviceBufferPointer,
    BYTE byFieldId)
{
    // Allocate site buffer G(x)
    deviceSU2* pG = NULL;
    checkCudaErrors(__cudaMalloc((void**)&pG, _HC_Volume * sizeof(deviceSU2)));

    preparethread;
    m_iIterate = 0;
#if !_CLG_DOUBLEFLOAT
    DOUBLE fTheta = 0.0;
#else
    Real fTheta = F(0.0);
#endif

    while (m_iIterate < m_iMaxIterate)
    {
        // Check convergence
        if (0 == m_iIterate % m_iCheckErrorStep)
        {
            _LAUNCH_KERNEL(_kernelCalculateMAGObjectiveSU2, block, threads,
                byFieldId,
                _D_RealThreadBuffer,
                pDeviceBufferPointer);
            fTheta = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer)
#if !_CLG_DOUBLEFLOAT
                   / (static_cast<DOUBLE>(_HC_Volume) * static_cast<DOUBLE>(_HC_Dir));
#else
                   / (_HC_Volume * _HC_Dir);
#endif

            appGeneral(_T("MAG (SU2) Iterate : %d, theta = %2.12f\n"), m_iIterate, fTheta);
            if (fTheta < m_fAccuracy)
            {
                cudaSafeFree(pG);
                return;
            }
        }

        // Red-black sweep
        _LAUNCH_KERNEL(_kernelCalculateGMAGSU2Odd, block, threads,
            byFieldId, pDeviceBufferPointer, m_fOmega, pG);
        _LAUNCH_KERNEL(_kernelGaugeTransformMAGSU2Odd, block, threads,
            byFieldId, pG, pDeviceBufferPointer);
        _LAUNCH_KERNEL(_kernelCalculateGMAGSU2Even, block, threads,
            byFieldId, pDeviceBufferPointer, m_fOmega, pG);
        _LAUNCH_KERNEL(_kernelGaugeTransformMAGSU2Even, block, threads,
            byFieldId, pG, pDeviceBufferPointer);

        ++m_iIterate;
    }

    appGeneral(_T("MAG (SU2) gauge fixing failed with last theta = %f\n"), fTheta);
    cudaSafeFree(pG);
}

// ============================================================================
// Host helper: runs the full MAG iteration for SU(3).
// ============================================================================
void CGaugeFixingMAG::GaugeFixingLoopSU3(
    deviceSU3* pDeviceBufferPointer,
    BYTE byFieldId)
{
    // Allocate site buffer G(x)
    deviceSU3* pG = NULL;
    checkCudaErrors(__cudaMalloc((void**)&pG, _HC_Volume * sizeof(deviceSU3)));

    preparethread;
    m_iIterate = 0;
#if !_CLG_DOUBLEFLOAT
    DOUBLE fTheta = 0.0;
#else
    Real fTheta = F(0.0);
#endif

    while (m_iIterate < m_iMaxIterate)
    {
        // Check convergence
        if (0 == m_iIterate % m_iCheckErrorStep)
        {
            _LAUNCH_KERNEL(_kernelCalculateMAGObjective, block, threads,
                byFieldId,
                _D_RealThreadBuffer,
                pDeviceBufferPointer);
            fTheta = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer)
#if !_CLG_DOUBLEFLOAT
                   / (static_cast<DOUBLE>(_HC_Volume) * static_cast<DOUBLE>(_HC_Dir));
#else
                   / (_HC_Volume * _HC_Dir);
#endif

            appGeneral(_T("MAG (SU3) Iterate : %d, theta = %2.12f\n"), m_iIterate, fTheta);
            if (fTheta < m_fAccuracy)
            {
                cudaSafeFree(pG);
                return;
            }
        }

        // Cabibbo-Marinari on SU(2) subgroups: R12, R23, R13
        // Do multiple passes per iteration for better convergence.
        for (INT pass = 0; pass < 2; ++pass)
        {
            for (BYTE bySubgroup = 0; bySubgroup < 3; ++bySubgroup)
            {
                _LAUNCH_KERNEL(_kernelCalculateGMAGSU3SubgroupOdd, block, threads,
                    byFieldId, bySubgroup, pDeviceBufferPointer, m_fOmega, pG);
                _LAUNCH_KERNEL(_kernelGaugeTransformMAGOdd, block, threads,
                    byFieldId, pG, pDeviceBufferPointer);
                _LAUNCH_KERNEL(_kernelCalculateGMAGSU3SubgroupEven, block, threads,
                    byFieldId, bySubgroup, pDeviceBufferPointer, m_fOmega, pG);
                _LAUNCH_KERNEL(_kernelGaugeTransformMAGEven, block, threads,
                    byFieldId, pG, pDeviceBufferPointer);
            }
        }

        ++m_iIterate;
    }

    appGeneral(_T("MAG (SU3) gauge fixing failed with last theta = %f\n"), fTheta);
    cudaSafeFree(pG);
}

__CLGIMPLEMENT_CLASS(CGaugeFixingMAG)

void CGaugeFixingMAG::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;

    if (!params.FetchValueReal(_T("Omega"), m_fOmega))
    {
        appGeneral(_T("CGaugeFixingMAG: Omega not set, set to 1.5 by default."));
    }

    if (!params.FetchValueReal(_T("Accuracy"), m_fAccuracy))
    {
        appGeneral(_T("CGaugeFixingMAG: Accuracy not set, set to 0.00000000001 by default."));
        if (m_fAccuracy < _CLG_FLT_EPSILON * F(2.0))
        {
            m_fAccuracy = _CLG_FLT_EPSILON * F(2.0);
            appGeneral(_T("Solver accuracy too small, set to be %2.18f\n"), m_fAccuracy);
        }
    }

    INT iValue = static_cast<INT>(m_iMaxIterate);
    if (!params.FetchValueINT(_T("MaxIterate"), iValue))
    {
        appGeneral(_T("CGaugeFixingMAG: MaxIterate not set, set to 100000 by default."));
    }
    m_iMaxIterate = static_cast<UINT>(iValue);

    iValue = static_cast<INT>(m_iCheckErrorStep);
    if (!params.FetchValueINT(_T("CheckErrorStep"), iValue))
    {
        appGeneral(_T("CGaugeFixingMAG: CheckErrorStep not set, set to 1000 by default."));
    }
    m_iCheckErrorStep = static_cast<UINT>(iValue);
}

void CGaugeFixingMAG::GaugeFixing(CFieldGauge* pResGauge)
{
    if (NULL == pResGauge)
    {
        appCrucial(_T("CGaugeFixingMAG: pResGauge is NULL!\n"));
        return;
    }

#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        //P4-2.5: gather -> rank0 global-context fix -> scatter flow, same as the
        //other fixers. The iteration helpers take a raw buffer so the gathered
        //global buffer is passed directly.
        const EFieldType eType = pResGauge->GetFieldType();
        CFieldGaugeSU2* pGaugeSU2 = NULL;
        CFieldGaugeSU3* pGaugeSU3 = NULL;
        UINT uiElemSize = 0;
        if (EFT_GaugeSU2 == eType)
        {
            pGaugeSU2 = dynamic_cast<CFieldGaugeSU2*>(pResGauge);
            uiElemSize = static_cast<UINT>(sizeof(deviceSU2));
        }
        else if (EFT_GaugeSU3 == eType)
        {
            pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pResGauge);
            uiElemSize = static_cast<UINT>(sizeof(deviceSU3));
        }
        else
        {
            appCrucial(_T("CGaugeFixingMAG: unsupported field type %d!\n"), pResGauge->GetFieldType());
            return;
        }

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
            if (NULL != pGaugeSU2)
            {
                GaugeFixingLoopSU2(reinterpret_cast<deviceSU2*>(pDevGlobal), pGaugeSU2->m_byFieldId);
            }
            else
            {
                GaugeFixingLoopSU3(reinterpret_cast<deviceSU3*>(pDevGlobal), pGaugeSU3->m_byFieldId);
            }
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
        GaugeFixingLoopSU2(
            dynamic_cast<CFieldGaugeSU2*>(pResGauge)->m_pDeviceData,
            pResGauge->m_byFieldId);
        break;
    case EFT_GaugeSU3:
        GaugeFixingLoopSU3(
            dynamic_cast<CFieldGaugeSU3*>(pResGauge)->m_pDeviceData,
            pResGauge->m_byFieldId);
        break;
    default:
        appCrucial(_T("CGaugeFixingMAG: unsupported field type %d!\n"), pResGauge->GetFieldType());
        break;
    }
}

void CGaugeFixingMAG::MaximalAbelianProjection(CFieldGaugeSU2* pGauge)
{
    if (NULL == pGauge)
    {
        appCrucial(_T("CGaugeFixingMAG::MaximalAbelianProjection: pGauge is NULL!\n"));
        return;
    }

    preparethread;
    _LAUNCH_KERNEL(_kernelMaximalAbelianProjectionSU2, block, threads,
        pGauge->m_byFieldId, pGauge->m_pDeviceData);
}

void CGaugeFixingMAG::MaximalAbelianProjection(CFieldGaugeSU3* pGauge)
{
    if (NULL == pGauge)
    {
        appCrucial(_T("CGaugeFixingMAG::MaximalAbelianProjection: pGauge is NULL!\n"));
        return;
    }

    preparethread;
    _LAUNCH_KERNEL(_kernelMaximalAbelianProjectionSU3, block, threads,
        pGauge->m_byFieldId, pGauge->m_pDeviceData);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CGaugeFixingMAG::CheckRes(const CFieldGauge* pGauge)
#else
Real CGaugeFixingMAG::CheckRes(const CFieldGauge* pGauge)
#endif
{
    if (NULL == pGauge)
    {
        appCrucial(_T("CGaugeFixingMAG::CheckRes: pGauge is NULL!\n"));
#if !_CLG_DOUBLEFLOAT
        return 0.0;
#else
        return F(0.0);
#endif
    }

#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        //P4-2.5: the deviation must be measured on the WHOLE lattice; same flow
        //as the other fixers. Non-root ranks report their local deviation; the
        //driver only compares rank 0.
        const EFieldType eType = pGauge->GetFieldType();
        const UINT uiElemSize = (EFT_GaugeSU2 == eType) ? static_cast<UINT>(sizeof(deviceSU2)) : static_cast<UINT>(sizeof(deviceSU3));
        const UINT uiBytesPerSite = uiElemSize * static_cast<UINT>(_HC_Dir);
        const UINT uiLocalBytes = uiElemSize * static_cast<UINT>(_HC_LinkCount);
        const BYTE* pLocalData = (const BYTE*)pGauge->GetData();

        BYTE* pHostLocal = (BYTE*)malloc(uiLocalBytes);
        checkCudaErrors(cudaMemcpy(pHostLocal, pLocalData, uiLocalBytes, cudaMemcpyDeviceToHost));

        UINT uiGlobalBytes = 0;
        BYTE* pGlobal = appGetComm()->GatherFieldToRoot(pHostLocal, uiBytesPerSite, uiGlobalBytes);
        free(pHostLocal);
        if (NULL == pGlobal)
        {
            return (EFT_GaugeSU2 == eType)
                ? CheckResLocalSU2(reinterpret_cast<const deviceSU2*>(pLocalData), pGauge->m_byFieldId)
                : CheckResLocalSU3(reinterpret_cast<const deviceSU3*>(pLocalData), pGauge->m_byFieldId);
        }

        BYTE* pDevGlobal = NULL;
        checkCudaErrors(__cudaMalloc((void**)&pDevGlobal, uiGlobalBytes));
        checkCudaErrors(cudaMemcpy(pDevGlobal, pGlobal, uiGlobalBytes, cudaMemcpyHostToDevice));
        free(pGlobal);

        MGEnterGlobalFixerContext();
        const DOUBLE fRet = (EFT_GaugeSU2 == eType)
            ? CheckResLocalSU2(reinterpret_cast<deviceSU2*>(pDevGlobal), pGauge->m_byFieldId)
            : CheckResLocalSU3(reinterpret_cast<deviceSU3*>(pDevGlobal), pGauge->m_byFieldId);
        MGExitGlobalFixerContext();

        checkCudaErrors(__cudaFree(pDevGlobal));
        return fRet;
    }
#endif

    return (EFT_GaugeSU2 == pGauge->GetFieldType())
        ? CheckResLocalSU2(reinterpret_cast<const deviceSU2*>(pGauge->GetData()), pGauge->m_byFieldId)
        : CheckResLocalSU3(reinterpret_cast<const deviceSU3*>(pGauge->GetData()), pGauge->m_byFieldId);
}

DOUBLE CGaugeFixingMAG::CheckResLocalSU2(const deviceSU2* pGaugeData, BYTE byFieldId)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelCalculateMAGObjectiveSU2, block, threads,
        byFieldId,
        _D_RealThreadBuffer,
        pGaugeData);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer)
#if !_CLG_DOUBLEFLOAT
         / (static_cast<DOUBLE>(_HC_Volume) * static_cast<DOUBLE>(_HC_Dir));
#else
         / (_HC_Volume * _HC_Dir);
#endif
}

DOUBLE CGaugeFixingMAG::CheckResLocalSU3(const deviceSU3* pGaugeData, BYTE byFieldId)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelCalculateMAGObjective, block, threads,
        byFieldId,
        _D_RealThreadBuffer,
        pGaugeData);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer)
#if !_CLG_DOUBLEFLOAT
         / (static_cast<DOUBLE>(_HC_Volume) * static_cast<DOUBLE>(_HC_Dir));
#else
         / (_HC_Volume * _HC_Dir);
#endif
}

CCString CGaugeFixingMAG::GetInfos(const CCString& tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeFixingMAG\n");
    sRet = sRet + tab + _T("Omega : ") + appToString(m_fOmega) + _T("\n");
    return sRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
