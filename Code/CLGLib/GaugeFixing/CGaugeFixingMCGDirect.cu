//=============================================================================
// FILENAME : CGaugeFixingMCGDirect.cu
//
// DESCRIPTION:
// Direct Maximal Center Gauge (MCG) fixing for SU(3).
//
// Standard algorithm (Montero, hep-lat/9906010):
// Directly maximizes the MCG functional
//   R = (1/(N_site * N_dim * N^2)) * sum_{x,\mu} |Tr U_\mu(x)|^2
// by local gauge updates using the Cabibbo-Marinari-Okawa method.
//
// The local quantity to maximize at each site x is:
//   R_x = sum_\mu |Tr{G(x) U_\mu(x)}|^2 + |Tr{U_\mu(x-\mu) G^\dagger(x)}|^2
//
// Local staple matrix:
//   A(x) = sum_\mu [ Tr[U_\mu(x)]* * U_\mu(x)
//                   + Tr[U_\mu(x-\mu)]* * U_\mu(x-\mu)^\dagger ]
//
// The optimal G(x) maximizes Re Tr[G(x) A(x)].
// Since CabbiboMarinariProj(M) finds U maximizing Re Tr[U M^\dagger],
// we project A(x)^\dagger instead of A(x).
//
// Overrelaxation: G_\omega = (1-\omega) I + \omega G, then re-project.
// Red-black checkerboard sweep for parallelization.
//
// Convergence criterion:
//   \theta = sum_{x,\mu} (N^2 - |Tr U_\mu(x)|^2) / (N^2 * V * D)
//   Stop when \theta < accuracy.
//
// Center Projection (Tucker & Stack, hep-lat/0110165):
//   After MCG gauge fixing, links can be projected to Z_3 center elements
//   by finding the nearest center phase for each link.
//   Z_3 = { exp(2*pi*i*m/3) * I | m = 0, 1, 2 }
//
// References:
//  hep-lat/9906010 -- Montero, "Study of SU(3) vortex-like configurations
//                     with a new maximal center gauge fixing method"
//  hep-lat/9708008 -- Brower et al., "Center vortices in SU(2) lattice gauge fields"
//  hep-lat/0003021 -- Langfeld et al., "SU(N) vortices and Wilson loops"
//  hep-lat/0110165 -- Tucker & Stack, "The Maximal Abelian Gauge in SU(3)"
//
// REVISION:
//  [05/14/2026 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CGaugeFixingMCGDirect.h"

__BEGIN_NAMESPACE

#pragma region kernels

// ============================================================================
// Compute trace-weighted staple for odd sites and project to SU(3).
//
// A(x) = sum_\mu [ Tr[U_\mu(x)]* * U_\mu(x)
//                 + Tr[U_\mu(x-\mu)]* * U_\mu(x-\mu)^\dagger ]
//
// The trace weight Tr[U]* comes from the derivative of |Tr U|^2.
//
// Overrelaxation: G = (1 - \omega) * I + \omega * A(x)
// Then project G^\dagger to SU(3) using Cabibbo-Marinari.
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGMCGOdd(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pU,
    Real fOmega,
    deviceSU3* pG)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    const deviceSU3 sunId = deviceSU3::makeSU3Id();
    if (site.IsDirichlet() || !sSite4.IsOdd())
    {
        return;
    }

    pG[uiSiteIndex] = deviceSU3::makeSU3Zero();

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        // Forward link: Tr[U_\mu(x)]* * U_\mu(x)
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const deviceSU3& Ufwd = pU[uiLinkIndex];
            CLGComplex trU = Ufwd.Tr();
            trU.y = -trU.y;  // conjugate: Tr[U]*
            pG[uiSiteIndex].AddDagger(Ufwd.MulCompC(trU));
        }
        else
        {
            // Boundary: treat as identity. Tr[I] = 3, so contribution is 3*I.
            pG[uiSiteIndex].Add(sunId.MulRealC(F(3.0)));
        }

        // Backward link: Tr[U_\mu(x-\mu)]* * U_\mu(x-\mu)^\dagger
        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];
        const deviceSU3 Ubwd = _deviceGetGaugeBCSU3DirOneSIndex(pU, site_m_mu);
        CLGComplex trU = Ubwd.Tr();
        trU.y = -trU.y;  // conjugate
        pG[uiSiteIndex].Add(Ubwd.MulCompC(trU));
    }

    // Overrelaxation: G = (1 - omega) * I + omega * G
    pG[uiSiteIndex].MulReal(fOmega);
    pG[uiSiteIndex].Add(sunId.MulRealC(F(1.0) - fOmega));

    // Cabibbo-Marinari projection on M finds U maximizing Re Tr[U M^\dagger].
    // We need g maximizing Re Tr[g A], so project A^\dagger instead of A.
    pG[uiSiteIndex].Dagger();
    pG[uiSiteIndex].CabbiboMarinariProj();
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGMCGEven(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pU,
    Real fOmega,
    deviceSU3* pG)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    const deviceSU3 sunId = deviceSU3::makeSU3Id();
    if (site.IsDirichlet() || sSite4.IsOdd())
    {
        return;
    }

    pG[uiSiteIndex] = deviceSU3::makeSU3Zero();

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const deviceSU3& Ufwd = pU[uiLinkIndex];
            CLGComplex trU = Ufwd.Tr();
            trU.y = -trU.y;
            pG[uiSiteIndex].AddDagger(Ufwd.MulCompC(trU));
        }
        else
        {
            pG[uiSiteIndex].Add(sunId.MulRealC(F(3.0)));
        }

        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];
        const deviceSU3 Ubwd = _deviceGetGaugeBCSU3DirOneSIndex(pU, site_m_mu);
        CLGComplex trU = Ubwd.Tr();
        trU.y = -trU.y;
        pG[uiSiteIndex].Add(Ubwd.MulCompC(trU));
    }

    pG[uiSiteIndex].MulReal(fOmega);
    pG[uiSiteIndex].Add(sunId.MulRealC(F(1.0) - fOmega));

    // Cabibbo-Marinari projection on M finds U maximizing Re Tr[U M^\dagger].
    // We need g maximizing Re Tr[g A], so project A^\dagger instead of A.
    pG[uiSiteIndex].Dagger();
    pG[uiSiteIndex].CabbiboMarinariProj();
}

// ============================================================================
// Gauge transform kernels (red-black sweep).
// Odd sites:  U'_\mu(x) = G(x) * U_\mu(x)
// Even sites: U'_\mu(x) = U_\mu(x) * G(x+\mu)^\dagger
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformMCGOdd(
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
_kernelGaugeTransformMCGEven(
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
// MCG convergence check: compute the MCG objective per site.
//
// MCG functional: F = sum_{x,\mu} |Tr U_\mu(x)|^2
//
// Convergence measure:
//   \theta = sum_{x,\mu} (N^2 - |Tr U_\mu(x)|^2) / (N^2 * V * D)
//
// For SU(3): N^2 = 9, D = 4, so denominator = 36 * V.
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateMCGObjective(
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
            Real trSq = __cuCabsSqf(U.Tr());
            objective += static_cast<DOUBLE>(F(9.0) - trSq);
        }
    }

    pDeviceRes[uiSiteIndex] = objective;
}

// ============================================================================
// Center projection: project each SU(3) link to the nearest Z_3 element.
//
// Z_3 center of SU(3): { exp(2*pi*i*m/3) * I | m = 0, 1, 2 }
//
// For each link U, compute theta = arg(Tr U) and find nearest center phase:
//   m = 0  if |theta| <= pi/3        -> nearest to I
//   m = 1  if  pi/3 < theta <= pi    -> nearest to exp(2*pi*i/3) I
//   m = 2  if -pi < theta < -pi/3    -> nearest to exp(-2*pi*i/3) I
//
// Reference: Tucker & Stack, hep-lat/0110165.
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelCenterProjectionSU3(
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

    const Real fPiOver3 = F(3.14159265358979323846) / F(3.0);
    const Real fTwoPiOver3 = F(2.0) * F(3.14159265358979323846) / F(3.0);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            deviceSU3& U = pGauge[uiLinkIndex];
            CLGComplex tr = U.Tr();
            Real theta = atan2f(tr.y, tr.x);  // arg(Tr U) in (-pi, pi]

            // Find nearest Z_3 element
            CLGComplex phase;
            if (theta > fPiOver3)
            {
                // Nearest to exp(2*pi*i/3)
                phase = _make_cuComplex(cosf(fTwoPiOver3), sinf(fTwoPiOver3));
            }
            else if (theta < -fPiOver3)
            {
                // Nearest to exp(-2*pi*i/3)
                phase = _make_cuComplex(cosf(-fTwoPiOver3), sinf(-fTwoPiOver3));
            }
            else
            {
                // Nearest to I
                phase = _make_cuComplex(F(1.0), F(0.0));
            }

            // Project to phase * I
            U = deviceSU3::makeSU3Zero();
            U.m_me[0] = phase;
            U.m_me[4] = phase;
            U.m_me[8] = phase;
        }
    }
}


/**
* Remove the Z_3 center from an SU(3) gauge field, in place.
* For each link U find the nearest center element Z = exp(2*pi*i*m/3)*I
* (same m selection as _kernelCenterProjectionSU3) and set
*   U <- Z^dag * U.
* The resulting links stay in SU(3); the center phase is projected out while
* the coset part SU(3)/Z_3 is kept.
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCenterRemoveSU3(
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

    const Real fPiOver3 = F(3.14159265358979323846) / F(3.0);
    const Real fTwoPiOver3 = F(2.0) * F(3.14159265358979323846) / F(3.0);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            deviceSU3& U = pGauge[uiLinkIndex];
            CLGComplex tr = U.Tr();
            Real theta = atan2f(tr.y, tr.x);

            // Nearest Z_3 element Z = exp(2*pi*i*m/3)*I; remove it: U <- Z^dag * U
            CLGComplex conjPhase;
            if (theta > fPiOver3)
            {
                // m = 1, Z = exp(2*pi*i/3), Z^dag = exp(-2*pi*i/3)
                conjPhase = _make_cuComplex(cosf(-fTwoPiOver3), sinf(-fTwoPiOver3));
            }
            else if (theta < -fPiOver3)
            {
                // m = 2, Z = exp(-2*pi*i/3), Z^dag = exp(2*pi*i/3)
                conjPhase = _make_cuComplex(cosf(fTwoPiOver3), sinf(fTwoPiOver3));
            }
            else
            {
                // m = 0, Z = I, Z^dag = I
                conjPhase = _make_cuComplex(F(1.0), F(0.0));
            }

            U.MulComp(conjPhase);
        }
    }
}

#pragma region SU2 kernels

// ============================================================================
// Compute trace-weighted staple for SU(2) odd sites and project to SU(2).
//
// Algorithm (Del Debbio et al., hep-lat/9610005):
// For SU(2) the center is Z_2 = {+I, -I}. The MCG maximizes
//   R = sum_{x,mu} |Tr U_mu(x)|^2.
// Since Tr(U) is real for SU(2), the local staple is:
//   A(x) = sum_mu [ Tr[U_mu(x)] * U_mu(x)
//                 + Tr[U_mu(x-mu)] * U_mu(x-mu)^dag ]
//
// For SU(2) the staple has the form A = [[a, b], [-b*, a*]].
// The optimal g(x) in SU(2) maximizing Re Tr[g A] is:
//   g = A^dag / |A|,  where |A| = sqrt(|a|^2 + |b|^2).
//
// Overrelaxation: G = (1 - omega) * I + omega * A, then project.
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGMCGSU2Odd(
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

    deviceSU2 A = deviceSU2::makeSU2Zero();

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        // Forward link: Tr[U_mu(x)] * U_mu(x)
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const deviceSU2& Ufwd = pU[uiLinkIndex];
            Real trU = Ufwd.ReTr();
            A.Add(Ufwd.MulRealC(trU));
        }
        else
        {
            // Boundary: U = I, Tr[I] = 2
            A.AddReal(F(2.0));
        }

        // Backward link: Tr[U_mu(x-mu)] * U_mu(x-mu)^dag
        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];
        const deviceSU2 Ubwd = _deviceGetGaugeBCDirOneSIndexT<deviceSU2>(pU, site_m_mu);
        Real trU = Ubwd.ReTr();
        A.AddDagger(Ubwd.MulRealC(trU));
    }

    // Overrelaxation
    A.MulReal(fOmega);
    A.AddReal(F(1.0) - fOmega);

    // Project to SU(2): g = A^dag / |A|.
    A.Dagger();
    Real norm = _sqrt(__cuCabsSqf(A.m_me[0]) + __cuCabsSqf(A.m_me[1]));

    deviceSU2 g = deviceSU2::makeSU2Id();
    if (norm > _CLG_FLT_EPSILON)
    {
        Real invNorm = __rcp(norm);
        g.m_me[0] = cuCmulf_cr(A.m_me[0], invNorm);
        g.m_me[1] = cuCmulf_cr(A.m_me[1], invNorm);
        // Enforce exact SU(2) structure: second row from first row
        g.m_me[2] = _make_cuComplex(-g.m_me[1].x, g.m_me[1].y);  // -g01*
        g.m_me[3] = _cuConjf(g.m_me[0]);  // g00*
    }

    pG[uiSiteIndex] = g;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateGMCGSU2Even(
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

    deviceSU2 A = deviceSU2::makeSU2Zero();

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            const deviceSU2& Ufwd = pU[uiLinkIndex];
            Real trU = Ufwd.ReTr();
            A.Add(Ufwd.MulRealC(trU));
        }
        else
        {
            A.AddReal(F(2.0));
        }

        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& site_m_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];
        const deviceSU2 Ubwd = _deviceGetGaugeBCDirOneSIndexT<deviceSU2>(pU, site_m_mu);
        Real trU = Ubwd.ReTr();
        A.AddDagger(Ubwd.MulRealC(trU));
    }

    A.MulReal(fOmega);
    A.AddReal(F(1.0) - fOmega);

    // Project to SU(2): g = A^dag / |A|.
    A.Dagger();
    Real norm = _sqrt(__cuCabsSqf(A.m_me[0]) + __cuCabsSqf(A.m_me[1]));

    deviceSU2 g = deviceSU2::makeSU2Id();
    if (norm > _CLG_FLT_EPSILON)
    {
        Real invNorm = __rcp(norm);
        g.m_me[0] = cuCmulf_cr(A.m_me[0], invNorm);
        g.m_me[1] = cuCmulf_cr(A.m_me[1], invNorm);
        // Enforce exact SU(2) structure: second row from first row
        g.m_me[2] = _make_cuComplex(-g.m_me[1].x, g.m_me[1].y);  // -g01*
        g.m_me[3] = _cuConjf(g.m_me[0]);  // g00*
    }

    pG[uiSiteIndex] = g;
}

// ============================================================================
// Gauge transform kernels for SU(2) (red-black sweep).
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformMCGSU2Odd(
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
_kernelGaugeTransformMCGSU2Even(
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
// MCG convergence check for SU(2).
//
// MCG functional: F = sum_{x,mu} |Tr U_mu(x)|^2
//
// Convergence measure:
//   theta = sum_{x,mu} (N^2 - |Tr U|^2) / (N^2 * V * D)
// For SU(2): N^2 = 4, D = 4, so denominator = 16 * V.
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateMCGObjectiveSU2(
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
            Real tr = U.ReTr();
            Real trSq = tr * tr;
            objective += static_cast<DOUBLE>(F(4.0) - trSq);
        }
    }

    pDeviceRes[uiSiteIndex] = objective;
}

// ============================================================================
// Center projection for SU(2): project each link to nearest Z_2 element.
//
// Z_2 center of SU(2): { +I, -I }
//
// Algorithm (Del Debbio et al., hep-lat/9610005):
//   For SU(2), Tr(U) is real. Project to sign(Tr U) * I.
// ============================================================================
__global__ void _CLG_LAUNCH_BOUND
_kernelCenterProjectionSU2(
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
            Real tr = U.ReTr();

            // Project to sign(Tr U) * I
            CLGComplex phase;
            if (tr > F(0.0))
            {
                phase = _make_cuComplex(F(1.0), F(0.0));
            }
            else
            {
                phase = _make_cuComplex(F(-1.0), F(0.0));
            }

            U = deviceSU2::makeSU2Zero();
            U.m_me[0] = phase;
            U.m_me[3] = phase;
        }
    }
}


/**
* Remove the Z_2 center from an SU(2) gauge field, in place.
* For each link U find the nearest Z_2 element Z = sign(Tr U)*I
* (Z^dag = Z for Z_2) and set
*   U <- Z * U.
* The resulting links stay in SU(2); the center sign is projected out while
* the coset part SU(2)/Z_2 is kept.
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCenterRemoveSU2(
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
            Real tr = U.ReTr();

            // Nearest Z_2 element Z = sign(Tr U)*I; Z^dag = Z, so U <- Z*U
            U.MulReal((tr > F(0.0)) ? F(1.0) : F(-1.0));
        }
    }
}

#pragma endregion

// ============================================================================
// Host helper: runs the full Direct MCG iteration for SU(2).
// ============================================================================
void CGaugeFixingMCGDirect::GaugeFixingLoopSU2(
        deviceSU2* pDeviceBufferPointer,
    BYTE byFieldId)
{
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
        if (0 == m_iIterate % m_iCheckErrorStep)
        {
            _LAUNCH_KERNEL(_kernelCalculateMCGObjectiveSU2, block, threads,
                byFieldId,
                _D_RealThreadBuffer,
                pDeviceBufferPointer);
            fTheta = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer)
#if !_CLG_DOUBLEFLOAT
                   / (static_cast<DOUBLE>(F(4.0)) * static_cast<DOUBLE>(_HC_Volume) * static_cast<DOUBLE>(_HC_Dir));
#else
                   / (F(4.0) * _HC_Volume * _HC_Dir);
#endif

            appGeneral(_T("MCG Direct (SU2) Iterate : %d, theta = %2.12f\n"), m_iIterate, fTheta);
            if (fTheta < m_fAccuracy)
            {
                cudaSafeFree(pG);
                return;
            }
        }

        _LAUNCH_KERNEL(_kernelCalculateGMCGSU2Odd, block, threads,
            byFieldId, pDeviceBufferPointer, m_fOmega, pG);
        _LAUNCH_KERNEL(_kernelGaugeTransformMCGSU2Odd, block, threads,
            byFieldId, pG, pDeviceBufferPointer);
        _LAUNCH_KERNEL(_kernelCalculateGMCGSU2Even, block, threads,
            byFieldId, pDeviceBufferPointer, m_fOmega, pG);
        _LAUNCH_KERNEL(_kernelGaugeTransformMCGSU2Even, block, threads,
            byFieldId, pG, pDeviceBufferPointer);

        ++m_iIterate;
    }

    appGeneral(_T("MCG Direct (SU2) gauge fixing failed with last theta = %f\n"), fTheta);
    cudaSafeFree(pG);
}

// ============================================================================
// Host helper: runs the full Direct MCG iteration for SU(3).
// ============================================================================
void CGaugeFixingMCGDirect::GaugeFixingLoopSU3(
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
            _LAUNCH_KERNEL(_kernelCalculateMCGObjective, block, threads,
                byFieldId,
                _D_RealThreadBuffer,
                pDeviceBufferPointer);
            fTheta = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer)
#if !_CLG_DOUBLEFLOAT
                   / (static_cast<DOUBLE>(F(9.0)) * static_cast<DOUBLE>(_HC_Volume) * static_cast<DOUBLE>(_HC_Dir));
#else
                   / (F(9.0) * _HC_Volume * _HC_Dir);
#endif

            appGeneral(_T("MCG Direct Iterate : %d, theta = %2.12f\n"), m_iIterate, fTheta);
            if (fTheta < m_fAccuracy)
            {
                cudaSafeFree(pG);
                return;
            }
        }

        // Red-black sweep
        _LAUNCH_KERNEL(_kernelCalculateGMCGOdd, block, threads,
            byFieldId, pDeviceBufferPointer, m_fOmega, pG);
        _LAUNCH_KERNEL(_kernelGaugeTransformMCGOdd, block, threads,
            byFieldId, pG, pDeviceBufferPointer);
        _LAUNCH_KERNEL(_kernelCalculateGMCGEven, block, threads,
            byFieldId, pDeviceBufferPointer, m_fOmega, pG);
        _LAUNCH_KERNEL(_kernelGaugeTransformMCGEven, block, threads,
            byFieldId, pG, pDeviceBufferPointer);

        ++m_iIterate;
    }

    appGeneral(_T("MCG Direct gauge fixing failed with last theta = %f\n"), fTheta);
    cudaSafeFree(pG);
}

__CLGIMPLEMENT_CLASS(CGaugeFixingMCGDirect)

void CGaugeFixingMCGDirect::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;

    if (!params.FetchValueReal(_T("Omega"), m_fOmega))
    {
        appGeneral(_T("CGaugeFixingMCGDirect: Omega not set, set to 1.5 by default."));
    }

    if (!params.FetchValueReal(_T("Accuracy"), m_fAccuracy))
    {
        appGeneral(_T("CGaugeFixingMCGDirect: Accuracy not set, set to 0.00000000001 by default."));
        if (m_fAccuracy < _CLG_FLT_EPSILON * F(2.0))
        {
            m_fAccuracy = _CLG_FLT_EPSILON * F(2.0);
            appGeneral(_T("Solver accuracy too small, set to be %2.18f\n"), m_fAccuracy);
        }
    }

    INT iValue = static_cast<INT>(m_iMaxIterate);
    if (!params.FetchValueINT(_T("MaxIterate"), iValue))
    {
        appGeneral(_T("CGaugeFixingMCGDirect: MaxIterate not set, set to 100000 by default."));
    }
    m_iMaxIterate = static_cast<UINT>(iValue);

    iValue = static_cast<INT>(m_iCheckErrorStep);
    if (!params.FetchValueINT(_T("CheckErrorStep"), iValue))
    {
        appGeneral(_T("CGaugeFixingMCGDirect: CheckErrorStep not set, set to 1000 by default."));
    }
    m_iCheckErrorStep = static_cast<UINT>(iValue);
}

void CGaugeFixingMCGDirect::GaugeFixing(CFieldGauge* pResGauge)
{
    if (NULL == pResGauge)
    {
        appCrucial(_T("CGaugeFixingMCGDirect: pResGauge is NULL!\n"));
        return;
    }

#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        //P4-2.5: gather -> rank0 global-context fix -> scatter flow, same as the
        //other fixers. The static iteration helpers take a raw buffer so the
        //gathered global buffer is passed directly.
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
            appCrucial(_T("CGaugeFixingMCGDirect: unsupported field type %d!\n"), pResGauge->GetFieldType());
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
        appCrucial(_T("CGaugeFixingMCGDirect: unsupported field type %d!\n"), pResGauge->GetFieldType());
        break;
    }
}

void CGaugeFixingMCGDirect::CenterProjection(CFieldGaugeSU2* pGauge)
{
    if (NULL == pGauge)
    {
        appCrucial(_T("CGaugeFixingMCGDirect::CenterProjection: pGauge (SU2) is NULL!\n"));
        return;
    }

    preparethread;
    _LAUNCH_KERNEL(_kernelCenterProjectionSU2, block, threads,
        pGauge->m_byFieldId, pGauge->m_pDeviceData);
}

void CGaugeFixingMCGDirect::CenterProjection(CFieldGaugeSU3* pGauge)
{
    if (NULL == pGauge)
    {
        appCrucial(_T("CGaugeFixingMCGDirect::CenterProjection: pGauge (SU3) is NULL!\n"));
        return;
    }

    preparethread;
    _LAUNCH_KERNEL(_kernelCenterProjectionSU3, block, threads,
        pGauge->m_byFieldId, pGauge->m_pDeviceData);
}

void CGaugeFixingMCGDirect::CenterRemove(CFieldGaugeSU2* pGauge)
{
    if (NULL == pGauge)
    {
        appCrucial(_T("CGaugeFixingMCGDirect::CenterRemove: pGauge (SU2) is NULL!\n"));
        return;
    }

    preparethread;
    _LAUNCH_KERNEL(_kernelCenterRemoveSU2, block, threads,
        pGauge->m_byFieldId, pGauge->m_pDeviceData);
}

void CGaugeFixingMCGDirect::CenterRemove(CFieldGaugeSU3* pGauge)
{
    if (NULL == pGauge)
    {
        appCrucial(_T("CGaugeFixingMCGDirect::CenterRemove: pGauge (SU3) is NULL!\n"));
        return;
    }

    preparethread;
    _LAUNCH_KERNEL(_kernelCenterRemoveSU3, block, threads,
        pGauge->m_byFieldId, pGauge->m_pDeviceData);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CGaugeFixingMCGDirect::CheckRes(const CFieldGauge* pGauge)
#else
Real CGaugeFixingMCGDirect::CheckRes(const CFieldGauge* pGauge)
#endif
{
    if (NULL == pGauge)
    {
        appCrucial(_T("CGaugeFixingMCGDirect::CheckRes: pGauge is NULL!\n"));
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

DOUBLE CGaugeFixingMCGDirect::CheckResLocalSU2(const deviceSU2* pGaugeData, BYTE byFieldId)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelCalculateMCGObjectiveSU2, block, threads,
        byFieldId,
        _D_RealThreadBuffer,
        pGaugeData);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer)
#if !_CLG_DOUBLEFLOAT
         / (static_cast<DOUBLE>(F(4.0)) * static_cast<DOUBLE>(_HC_Volume) * static_cast<DOUBLE>(_HC_Dir));
#else
         / (F(4.0) * _HC_Volume * _HC_Dir);
#endif
}

DOUBLE CGaugeFixingMCGDirect::CheckResLocalSU3(const deviceSU3* pGaugeData, BYTE byFieldId)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelCalculateMCGObjective, block, threads,
        byFieldId,
        _D_RealThreadBuffer,
        pGaugeData);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer)
#if !_CLG_DOUBLEFLOAT
         / (static_cast<DOUBLE>(F(9.0)) * static_cast<DOUBLE>(_HC_Volume) * static_cast<DOUBLE>(_HC_Dir));
#else
         / (F(9.0) * _HC_Volume * _HC_Dir);
#endif
}

CCString CGaugeFixingMCGDirect::GetInfos(const CCString& tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeFixingMCGDirect\n");
    sRet = sRet + tab + _T("Omega : ") + appToString(m_fOmega) + _T("\n");
    return sRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
