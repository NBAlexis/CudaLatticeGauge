//=============================================================================
// FILENAME : CActionGaugePlaquettePSU3WithBoundary.cu
//
// DESCRIPTION:
// Implementation of the fundamental-lift Wilson/Villain PSU(3) action with a
// Z3 2-form field B.
//
//   S(U,B) = -(Beta / 3) * sum_p Re[B_p Tr(U_p)]
//
// Energy / force / Z3 heatbath sweeps:
//   - energy: full plaquette sum with the canonical B_p
//   - force : Y_l = -(Beta/3)/2 * sum_{p contains l} conj(B_oriented(p,l)) * Sigma_p,
//             the integrator applies U*Y^dag + TA afterwards
//   - AllowMonopole=1: independent per-plaquette heatbath (unconstrained)
//   - AllowMonopole=0: link-star coboundary heatbath (4 dirs x 2 parity) plus
//                      global coclosed-sheet heatbath for the H2 twist sectors
//
// REVISION:
//  [08/09/26]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionGaugePlaquettePSU3WithBoundary.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquettePSU3WithBoundary)

#pragma region Z3 device helpers

/**
* Map a deviceZN<3> to its exponent index {0,1,2} by comparing
* Re(z * conj(zeta^k)) and taking the largest (robust to float noise near the
* roots, no atan2 needed).
*/
__device__ __inline__ UINT DeviceZ3ToIndex(const deviceZN<3>& z)
{
    Real fBest = -F(2.0);
    UINT kBest = 0;
    for (UINT k = 0; k < 3; ++k)
    {
        const Real fTheta = F(2.0) * PI * static_cast<Real>(k) / F(3.0);
        const CLGComplex zetaConj = _make_cuComplex(_cos(fTheta), -_sin(fTheta));
        const CLGComplex prod = _cuCmulf(z.m_me, zetaConj);
        if (prod.x > fBest)
        {
            fBest = prod.x;
            kBest = k;
        }
    }
    return kBest;
}

/**
* b + delta (mod 3), result in {0,1,2}.
*/
__device__ __inline__ UINT DeviceZ3Add(UINT a, INT iDelta)
{
    INT r = static_cast<INT>(a) + iDelta;
    r %= 3;
    if (r < 0)
    {
        r += 3;
    }
    return static_cast<UINT>(r);
}

/**
* canonical plaquette index for mu < nu (xy,xz,xt,yz,yt,zt).
*/
__device__ __inline__ BYTE DevicePlaqIndex(BYTE mu, BYTE nu)
{
    return static_cast<BYTE>(mu * (2 * _DC_Dir - mu - 1) / 2 + (nu - mu - 1));
}

/**
* Move a site index by +1/-1 in a given direction with periodic wrap.
*/
__device__ __inline__ UINT DeviceSiteMove(UINT uiSiteIndex, BYTE byDir, INT iDelta)
{
    const SSmallInt4 s = __deviceSiteIndexToInt4(uiSiteIndex);
    const INT iLen = static_cast<INT>(_constIntegers[ECI_Lx + byDir]);
    SSmallInt4 t = s;
    INT iNew = 0;
    if (0 == byDir) iNew = static_cast<INT>(s.x) + iDelta;
    else if (1 == byDir) iNew = static_cast<INT>(s.y) + iDelta;
    else if (2 == byDir) iNew = static_cast<INT>(s.z) + iDelta;
    else iNew = static_cast<INT>(s.w) + iDelta;
    iNew = (iNew + iLen) % iLen;
    if (iNew < 0) iNew += iLen;
    if (0 == byDir) t.x = static_cast<SCHAR>(iNew);
    else if (1 == byDir) t.y = static_cast<SCHAR>(iNew);
    else if (2 == byDir) t.z = static_cast<SCHAR>(iNew);
    else t.w = static_cast<SCHAR>(iNew);
    return _deviceGetSiteIndex(t);
}

/**
* Build the canonical plaquette U_p = U_mu(x) U_nu(x+mu) U_mu(x+nu)^dag U_nu(x)^dag
* for plaqIndex 0..5 at baseSite, following the plaquette cache contract.
*/
template<typename deviceGauge>
__device__ __inline__ deviceGauge DeviceBuildCanonicalPlaquette(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pPlaqCache,
    UINT uiSiteIndex,
    BYTE byPlaqIndex)
{
#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT uiPlaqLength = appGetLattice()->m_pIndexCache->m_uiPlaqutteLength;
    const UINT uiPlaqCountPerSite = appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite;
#else
    const UINT uiPlaqLength = 4;
    const UINT uiPlaqCountPerSite = 6;
#endif
    const UINT uiPlaqCountAllSite = uiPlaqLength * uiPlaqCountPerSite;
    const UINT uiBase = uiSiteIndex * uiPlaqCountAllSite + byPlaqIndex * uiPlaqLength;

    SIndex first = pPlaqCache[uiBase + 0];
    deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
    if (first.NeedToDagger())
    {
        toAdd.Dagger();
    }

    for (BYTE j = 1; j < uiPlaqLength; ++j)
    {
        first = pPlaqCache[uiBase + j];
        deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
        if (first.NeedToDagger())
        {
            toAdd.MulDagger(toMul);
        }
        else
        {
            toAdd.Mul(toMul);
        }
    }
    return toAdd;
}

/**
* Numerical-stable categorical sampling over three log weights.
* Consumes exactly one uniform random from the stream fatIndex.
*/
__device__ __inline__ UINT DeviceSampleCategorical3(const DOUBLE logWeight[3], UINT uiSeed)
{
    DOUBLE dMax = logWeight[0];
    if (logWeight[1] > dMax) dMax = logWeight[1];
    if (logWeight[2] > dMax) dMax = logWeight[2];

    DOUBLE w[3];
    DOUBLE dSum = 0.0;
    for (UINT i = 0; i < 3; ++i)
    {
        w[i] = _exp(logWeight[i] - dMax);
        dSum += w[i];
    }

    const DOUBLE r = static_cast<DOUBLE>(_deviceRandomF(uiSeed)) * dSum;
    if (r < w[0]) return 0;
    if (r < w[0] + w[1]) return 1;
    return 2;
}

struct SZ3OrientedPlaquette
{
    UINT uiSite;      // canonical base site of the plaquette
    BYTE byPlaqIndex; // 0..5
    SCHAR bySign;     // +1: B_oriented = B_p; -1: B_oriented = conj(B_p)
};

/**
* The canonical plaquette associated with one staple (forward/backward) of the
* link (uiSiteIndex, byAlpha) in the other direction byNu.
*
*   forward  | alpha < nu | base x    | (alpha,nu) | +1
*   forward  | alpha > nu | base x    | (nu,alpha) | -1
*   backward | alpha < nu | base x-nu | (alpha,nu) | -1
*   backward | alpha > nu | base x-nu | (nu,alpha) | +1
*
* This is also the d-lambda incidence table used by the link-star move.
*/
__device__ __inline__ SZ3OrientedPlaquette DeviceOrientedPlaquetteOfLink(
    UINT uiSiteIndex, BYTE byAlpha, BYTE byNu, INT iForwardBackward)
{
    SZ3OrientedPlaquette ret;
    const BYTE mu = byAlpha < byNu ? byAlpha : byNu;
    const BYTE nu = byAlpha < byNu ? byNu : byAlpha;
    ret.byPlaqIndex = DevicePlaqIndex(mu, nu);

    if (iForwardBackward > 0)
    {
        // forward: base = x
        ret.uiSite = uiSiteIndex;
        ret.bySign = (byAlpha < byNu) ? 1 : -1;
    }
    else
    {
        // backward: base = x - nu
        ret.uiSite = DeviceSiteMove(uiSiteIndex, byNu, -1);
        ret.bySign = (byAlpha < byNu) ? -1 : 1;
    }
    return ret;
}

/**
* Read the Z3 element B[plaqIndex, site].
*/
__device__ __inline__ deviceZN<3> DeviceReadB(
    const deviceZN<3>* __restrict__ pBoundary, UINT uiSiteCount, UINT uiSiteIndex, BYTE byPlaqIndex)
{
    return pBoundary[byPlaqIndex * uiSiteCount + uiSiteIndex];
}

#pragma endregion

#pragma region kernels

/**
* Energy: S_dyn = -(Beta/3) * sum_p Re[B_p Tr(U_p)]
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelEnergyPSU3WithBoundary(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const deviceZN<3>* __restrict__ pBoundary,
    const SIndex* __restrict__ pPlaqCache,
    Real betaOverN,
    UINT uiSiteCount,
    DOUBLE* results)
{
    intokernal;

    DOUBLE resThisThread = 0.0;

    for (BYTE i = 0; i < 6; ++i)
    {
        const deviceGauge up = DeviceBuildCanonicalPlaquette(byFieldId, pDeviceData, pPlaqCache, uiSiteIndex, i);
        const deviceZN<3> b = pBoundary[i * uiSiteCount + uiSiteIndex];
        const CLGComplex ztr = _cuCmulf(b.m_me, up.Tr());
        resThisThread -= static_cast<DOUBLE>(betaOverN) * ztr.x;
    }

    results[uiSiteIndex] = resThisThread;
}

/**
* Force: for each link, accumulate
*   Y_l = -(Beta/3)/2 * sum_{p contains l} conj(B_oriented(p,l)) * Sigma_p
* into the shared force field. The integrator then applies U*Y^dag + TA.
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelForcePSU3WithBoundary(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const deviceZN<3>* __restrict__ pBoundary,
    const SIndex* __restrict__ pStapleCache,
    Real betaOverN,
    UINT uiSiteCount,
    deviceGauge* pForceData)
{
    intokernalDir_NoDir;

    const BYTE byAlpha = static_cast<BYTE>(uiLinkIndex % _DC_Dir);

    deviceGauge res = _makeZero<deviceGauge>();
    const deviceGauge u = pDeviceData[uiLinkIndex];

    // staple cache: for i != alpha (ascending), [forward(3 links), backward(3 links)]
    UINT uiStapleIndex = 0;
    for (BYTE i = 0; i < _DC_Dir; ++i)
    {
        if (i == byAlpha)
        {
            continue;
        }
        const BYTE byNu = i;
        for (INT fb = 0; fb <= 1; ++fb) // 0: forward, 1: backward
        {
            const UINT uiStapleBase = uiLinkIndex * (_DC_Dir - 1) * 2 * 3 + uiStapleIndex * 3;
            SIndex first = pStapleCache[uiStapleBase + 0];
            deviceGauge staple(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
            if (first.NeedToDagger())
            {
                staple.Dagger();
            }
            for (BYTE j = 1; j < 3; ++j)
            {
                first = pStapleCache[uiStapleBase + j];
                deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
                if (first.NeedToDagger())
                {
                    staple.MulDagger(toMul);
                }
                else
                {
                    staple.Mul(toMul);
                }
            }

            // B_oriented: +1 -> B_p, -1 -> conj(B_p)
            const SZ3OrientedPlaquette op = DeviceOrientedPlaquetteOfLink(uiSiteIndex, byAlpha, byNu, fb == 0 ? 1 : -1);
            deviceZN<3> bOriented = pBoundary[op.byPlaqIndex * uiSiteCount + op.uiSite];
            if (op.bySign < 0)
            {
                bOriented.Dagger();
            }
            // Y += conj(B_oriented) * staple
            bOriented.Dagger();
            const deviceGauge term = staple.MulCompC(bOriented.m_me);
            res.Add(term);

            ++uiStapleIndex;
        }
    }

    res.MulReal(-betaOverN * F(0.5));

    // force is additive
    _add(pForceData[uiLinkIndex], res);
}

/**
* Unconstrained per-plaquette heatbath (AllowMonopole=1).
* One thread owns one site and updates its six plaquettes sequentially from
* the same random stream (consistent with the tensor2 initializer).
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelSweepZ3PlaquetteUnconstrained(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    deviceZN<3>* __restrict__ pBoundary,
    const SIndex* __restrict__ pPlaqCache,
    Real betaOverN,
    UINT uiSiteCount)
{
    intokernal;

    const UINT uiSeedIndex = _deviceGetLinkIndex(uiSiteIndex, 0);

    for (BYTE plaq = 0; plaq < 6; ++plaq)
    {
        const deviceGauge up = DeviceBuildCanonicalPlaquette(byFieldId, pDeviceData, pPlaqCache, uiSiteIndex, plaq);
        DOUBLE logWeight[3];
        for (UINT k = 0; k < 3; ++k)
        {
            const deviceZN<3> z = deviceZN<3>::makeAsK(k);
            const CLGComplex ztr = _cuCmulf(z.m_me, up.Tr());
            logWeight[k] = static_cast<DOUBLE>(betaOverN) * ztr.x;
        }
        const UINT kSel = DeviceSampleCategorical3(logWeight, uiSeedIndex);
        pBoundary[plaq * uiSiteCount + uiSiteIndex] = deviceZN<3>::makeAsK(kSel);
    }
}

/**
* Link-star coboundary heatbath (AllowMonopole=0).
* One thread owns one link (fixed direction byAlpha, fixed parity). The six
* incident plaquettes are updated together by a three-state orbit:
*   B_p(q) = B_p(old) * zeta^(epsilon(p,l) * q),  q = 0,1,2
* which is an exact coboundary update: dB and the global twist are unchanged.
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelSweepZ3LinkStar(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    deviceZN<3>* __restrict__ pBoundary,
    const SIndex* __restrict__ pPlaqCache,
    BYTE byAlpha,
    UBOOL bEven,
    Real betaOverN,
    UINT uiSiteCount)
{
    intokernal;

    const SSmallInt4 sSite = __deviceSiteIndexToInt4(uiSiteIndex);
    const UINT uiParity = static_cast<UINT>((sSite.x + sSite.y + sSite.z + sSite.w) & 1);
    if (uiParity != (bEven ? 0u : 1u))
    {
        return;
    }

    DOUBLE logWeight[3] = { 0.0, 0.0, 0.0 };

    for (BYTE i = 0; i < _DC_Dir; ++i)
    {
        if (i == byAlpha)
        {
            continue;
        }
        for (INT fb = 0; fb <= 1; ++fb)
        {
            const SZ3OrientedPlaquette op = DeviceOrientedPlaquetteOfLink(uiSiteIndex, byAlpha, i, fb == 0 ? 1 : -1);
            const deviceZN<3> bOld = pBoundary[op.byPlaqIndex * uiSiteCount + op.uiSite];
            const deviceGauge up = DeviceBuildCanonicalPlaquette(byFieldId, pDeviceData, pPlaqCache, op.uiSite, op.byPlaqIndex);
            const UINT kOld = DeviceZ3ToIndex(bOld);

            for (UINT q = 0; q < 3; ++q)
            {
                const UINT kNew = DeviceZ3Add(kOld, op.bySign * static_cast<INT>(q));
                const deviceZN<3> z = deviceZN<3>::makeAsK(kNew);
                const CLGComplex ztr = _cuCmulf(z.m_me, up.Tr());
                logWeight[q] += static_cast<DOUBLE>(betaOverN) * ztr.x;
            }
        }
    }

    const UINT qSel = DeviceSampleCategorical3(logWeight, _deviceGetLinkIndex(uiSiteIndex, byAlpha));

    for (BYTE i = 0; i < _DC_Dir; ++i)
    {
        if (i == byAlpha)
        {
            continue;
        }
        for (INT fb = 0; fb <= 1; ++fb)
        {
            const SZ3OrientedPlaquette op = DeviceOrientedPlaquetteOfLink(uiSiteIndex, byAlpha, i, fb == 0 ? 1 : -1);
            const UINT kOld = DeviceZ3ToIndex(pBoundary[op.byPlaqIndex * uiSiteCount + op.uiSite]);
            const UINT kNew = DeviceZ3Add(kOld, op.bySign * static_cast<INT>(qSel));
            pBoundary[op.byPlaqIndex * uiSiteCount + op.uiSite] = deviceZN<3>::makeAsK(kNew);
        }
    }
}

/**
* Global coclosed-sheet delta evaluation.
* Only the plaquettes of the fixed representative sheet Sigma_mu_nu contribute:
*   logWeight[q] = +(Beta/3) Re[(B_p * zeta^q) Tr(U_p)]  over p in sheet
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelEvaluateSheetDelta(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const deviceZN<3>* __restrict__ pBoundary,
    const SIndex* __restrict__ pPlaqCache,
    Real betaOverN,
    BYTE byPlaqIndex,
    BYTE byMu,
    BYTE byNu,
    UINT uiFixedMu,
    UINT uiFixedNu,
    UINT q,
    UINT uiSiteCount,
    DOUBLE* results)
{
    intokernal;

    const SSmallInt4 sSite = __deviceSiteIndexToInt4(uiSiteIndex);
    UINT uiCoord[4];
    uiCoord[0] = static_cast<UINT>(sSite.x);
    uiCoord[1] = static_cast<UINT>(sSite.y);
    uiCoord[2] = static_cast<UINT>(sSite.z);
    uiCoord[3] = static_cast<UINT>(sSite.w);

    if (uiCoord[byMu] != uiFixedMu || uiCoord[byNu] != uiFixedNu)
    {
        results[uiSiteIndex] = 0.0;
        return;
    }

    const deviceZN<3> b = pBoundary[byPlaqIndex * uiSiteCount + uiSiteIndex];
    const deviceGauge up = DeviceBuildCanonicalPlaquette(byFieldId, pDeviceData, pPlaqCache, uiSiteIndex, byPlaqIndex);

    const deviceZN<3> zq = deviceZN<3>::makeAsK(q);
    const CLGComplex zb = _cuCmulf(zq.m_me, b.m_me);
    const CLGComplex ztr = _cuCmulf(zb, up.Tr());
    results[uiSiteIndex] = static_cast<DOUBLE>(betaOverN) * ztr.x;
}

/**
* Sample the sheet q from the three accumulated log weights with a single
* device thread (consumes exactly one uniform random from the stream).
*/
__global__ void _kernelSampleSheetQ(
    DOUBLE dW0, DOUBLE dW1, DOUBLE dW2,
    UINT uiSeed,
    UINT* pOutQ)
{
    if (0 != (threadIdx.x + blockIdx.x * blockDim.x))
    {
        return;
    }
    DOUBLE logWeight[3] = { dW0, dW1, dW2 };
    *pOutQ = DeviceSampleCategorical3(logWeight, uiSeed);
}

/**
* Apply the chosen twist q to every plaquette of the fixed sheet.
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelApplySheet(
    deviceZN<3>* __restrict__ pBoundary,
    BYTE byPlaqIndex,
    BYTE byMu,
    BYTE byNu,
    UINT uiFixedMu,
    UINT uiFixedNu,
    UINT q,
    UINT uiSiteCount)
{
    intokernal;

    const SSmallInt4 sSite = __deviceSiteIndexToInt4(uiSiteIndex);
    UINT uiCoord[4];
    uiCoord[0] = static_cast<UINT>(sSite.x);
    uiCoord[1] = static_cast<UINT>(sSite.y);
    uiCoord[2] = static_cast<UINT>(sSite.z);
    uiCoord[3] = static_cast<UINT>(sSite.w);

    if (uiCoord[byMu] != uiFixedMu || uiCoord[byNu] != uiFixedNu)
    {
        return;
    }

    const UINT k = DeviceZ3Add(DeviceZ3ToIndex(pBoundary[byPlaqIndex * uiSiteCount + uiSiteIndex]), static_cast<INT>(q));
    pBoundary[byPlaqIndex * uiSiteCount + uiSiteIndex] = deviceZN<3>::makeAsK(k);
}

/**
* Monopole diagnostic: for every cube (x, mu<nu<rho), compute the Z3 charge
*   m = b_nu_rho(x+mu) - b_nu_rho(x) - b_mu_rho(x+nu) + b_mu_rho(x)
*     + b_mu_nu(x+rho) - b_mu_nu(x)   (mod 3)
* and count the non-zero cubes of each site (4 orientations: xyz, xyt, xzt, yzt).
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCountMonopole(
    const deviceZN<3>* __restrict__ pBoundary,
    UINT uiSiteCount,
    UINT* pCount)
{
    intokernal;

    UINT uiLocal = 0;
    const UINT aiOrient[4][3] = {
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 2, 3 },
        { 1, 2, 3 },
    };

    for (UINT o = 0; o < 4; ++o)
    {
        const BYTE mu = static_cast<BYTE>(aiOrient[o][0]);
        const BYTE nu = static_cast<BYTE>(aiOrient[o][1]);
        const BYTE rho = static_cast<BYTE>(aiOrient[o][2]);

        const UINT bNr_xm = DeviceZ3ToIndex(pBoundary[DevicePlaqIndex(nu, rho) * uiSiteCount + DeviceSiteMove(uiSiteIndex, mu, 1)]);
        const UINT bNr_x  = DeviceZ3ToIndex(pBoundary[DevicePlaqIndex(nu, rho) * uiSiteCount + uiSiteIndex]);
        const UINT bMr_xn = DeviceZ3ToIndex(pBoundary[DevicePlaqIndex(mu, rho) * uiSiteCount + DeviceSiteMove(uiSiteIndex, nu, 1)]);
        const UINT bMr_x  = DeviceZ3ToIndex(pBoundary[DevicePlaqIndex(mu, rho) * uiSiteCount + uiSiteIndex]);
        const UINT bMn_xr = DeviceZ3ToIndex(pBoundary[DevicePlaqIndex(mu, nu) * uiSiteCount + DeviceSiteMove(uiSiteIndex, rho, 1)]);
        const UINT bMn_x  = DeviceZ3ToIndex(pBoundary[DevicePlaqIndex(mu, nu) * uiSiteCount + uiSiteIndex]);

        INT m = static_cast<INT>(bNr_xm) - static_cast<INT>(bNr_x)
              - static_cast<INT>(bMr_xn) + static_cast<INT>(bMr_x)
              + static_cast<INT>(bMn_xr) - static_cast<INT>(bMn_x);
        m %= 3;
        if (m != 0)
        {
            ++uiLocal;
        }
    }

    pCount[uiSiteIndex] = uiLocal;
}

#pragma endregion


//=============================================================================
// class implementation
//=============================================================================

CActionGaugePlaquettePSU3WithBoundary::CActionGaugePlaquettePSU3WithBoundary()
    : CAction()
    , m_byBoundaryTensor2FieldId(0)
    , m_bAllowMonopole(FALSE)
    , m_uiZ3SweepsPerTrajectory(1)
    , m_uiGlobalTwistSweepsPerTrajectory(1)
    , m_bCheckMonopole(FALSE)
    , m_pBoundaryField(NULL)
    , m_fPostSweepEnergy(0.0)
    , m_bPostSweepEnergyValid(FALSE)
{
    for (UINT i = 0; i < 4; ++i)
    {
        m_aiTwistSheetPosition[i] = 0;
    }
}

INT CActionGaugePlaquettePSU3WithBoundary::FindGaugeIndexById(INT num, const CFieldGauge* const* gaugeFields, BYTE byFieldId)
{
    for (INT i = 0; i < num; ++i)
    {
        if (NULL != gaugeFields[i] && gaugeFields[i]->m_byFieldId == byFieldId)
        {
            return i;
        }
    }
    return -1;
}

INT CActionGaugePlaquettePSU3WithBoundary::FindTensor2IndexById(INT num, const CFieldTensor2* const* tensor2Fields, BYTE byFieldId)
{
    for (INT i = 0; i < num; ++i)
    {
        if (NULL != tensor2Fields[i] && tensor2Fields[i]->m_byFieldId == byFieldId)
        {
            return i;
        }
    }
    return -1;
}

void CActionGaugePlaquettePSU3WithBoundary::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    //---------------------------------------------------------------------
    // parameters
    //---------------------------------------------------------------------
    INT iValue = 0;
    if (!param.FetchValueINT(_T("BoundaryFieldId"), iValue) || iValue <= 0)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary: BoundaryFieldId must be a positive field id.\n"));
        _FAIL_EXIT;
        return;
    }
    m_byBoundaryTensor2FieldId = static_cast<BYTE>(iValue);

    iValue = 0;
    param.FetchValueINT(_T("AllowMonopole"), iValue);
    m_bAllowMonopole = (0 != iValue);

    UINT uiValue = 1;
    if (param.FetchValueINT(_T("Z3Sweeps"), iValue))
    {
        m_uiZ3SweepsPerTrajectory = static_cast<UINT>(iValue < 0 ? 0 : iValue);
    }
    if (param.FetchValueINT(_T("GlobalTwistSweeps"), iValue))
    {
        m_uiGlobalTwistSweepsPerTrajectory = static_cast<UINT>(iValue < 0 ? 0 : iValue);
    }
    iValue = 0;
    param.FetchValueINT(_T("CheckMonopole"), iValue);
    m_bCheckMonopole = (0 != iValue);

    iValue = 0;
    param.FetchValueINT(_T("DiscreteTheta"), iValue);
    if (0 != iValue)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary: only DiscreteTheta=0 is supported in this version.\n"));
        _FAIL_EXIT;
        return;
    }

    TArray<INT> aiTwist;
    if (param.FetchValueArrayINT(_T("TwistSheetPosition"), aiTwist) && aiTwist.Num() == 4)
    {
        for (UINT i = 0; i < 4; ++i)
        {
            m_aiTwistSheetPosition[i] = static_cast<UINT>(aiTwist[i] < 0 ? 0 : aiTwist[i]);
        }
    }

    //---------------------------------------------------------------------
    // gauge field checks: exactly one SU(3) gauge field, not the boundary id
    //---------------------------------------------------------------------
    if (1 != m_byGaugeFieldIds.Num())
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary: exactly one GaugeFields entry is required.\n"));
        _FAIL_EXIT;
        return;
    }
    if (m_byGaugeFieldIds[0] == m_byBoundaryTensor2FieldId)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary: gauge field id and boundary field id must differ.\n"));
        _FAIL_EXIT;
        return;
    }
    CField* pGaugeField = pOwner->GetFieldById(m_byGaugeFieldIds[0]);
    CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pGaugeField);
    if (NULL == pGaugeSU3)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary: GaugeFields[0] must be an EFT_GaugeSU3 field.\n"));
        _FAIL_EXIT;
        return;
    }

    //---------------------------------------------------------------------
    // boundary field checks
    //---------------------------------------------------------------------
    CField* pBoundaryField = pOwner->GetFieldById(m_byBoundaryTensor2FieldId);
    CFieldTensor2Z3* pBoundary = dynamic_cast<CFieldTensor2Z3*>(pBoundaryField);
    if (NULL == pBoundary)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary: BoundaryFieldId must point to a CFieldTensor2Z3 field.\n"));
        _FAIL_EXIT;
        return;
    }
    if (!pBoundary->IsDynamic())
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary: the boundary tensor2 field must be Dynamic: 1 (the integrator working copy is required).\n"));
        _FAIL_EXIT;
        return;
    }

    //---------------------------------------------------------------------
    // lattice checks
    //---------------------------------------------------------------------
    if (4 != _HC_Dim)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary: only 4D lattices are supported.\n"));
        _FAIL_EXIT;
        return;
    }

    if (m_bAllowMonopole)
    {
        if (0 != m_uiGlobalTwistSweepsPerTrajectory)
        {
            appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary: AllowMonopole=1 requires GlobalTwistSweeps=0 (global twist is not a protected sector in the unrestricted ensemble).\n"));
            _FAIL_EXIT;
            return;
        }
    }
    else
    {
        // flat mode requires even extents for the checkerboard link-star schedule
        if ((_HC_Lx % 2) || (_HC_Ly % 2) || (_HC_Lz % 2) || (_HC_Lt % 2))
        {
            appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary: AllowMonopole=0 requires all lattice extents to be even (checkerboard link-star schedule).\n"));
            _FAIL_EXIT;
            return;
        }
        if (0 == m_uiGlobalTwistSweepsPerTrajectory)
        {
            appGeneral(_T("CActionGaugePlaquettePSU3WithBoundary: GlobalTwistSweeps=0, the global twist sector is FIXED (not summing all PSU(3) bundles).\n"));
        }
        if (0 == m_uiZ3SweepsPerTrajectory)
        {
            appGeneral(_T("CActionGaugePlaquettePSU3WithBoundary: Z3Sweeps=0, the local B degrees of freedom are FROZEN (not a full joint PSU(3) Monte Carlo).\n"));
        }
    }

    //---------------------------------------------------------------------
    // B validity: flatness when AllowMonopole=0
    //---------------------------------------------------------------------
    if (!m_bAllowMonopole)
    {
        const UINT uiMonopoles = CountMonopoles(pBoundary);
        if (0 != uiMonopoles)
        {
            appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary: initial B is not flat (dB != 0, %u monopole cubes), rejected.\n"), uiMonopoles);
            _FAIL_EXIT;
            return;
        }
    }

    m_pBoundaryField = pBoundary;

    appGeneral(_T("CActionGaugePlaquettePSU3WithBoundary: Beta=%f (fundamental-lift coupling), AllowMonopole=%d, Z3Sweeps=%u, GlobalTwistSweeps=%u, BoundaryFieldId=%d\n"),
        m_fBetaOverN * 3.0, m_bAllowMonopole ? 1 : 0, m_uiZ3SweepsPerTrajectory, m_uiGlobalTwistSweepsPerTrajectory, m_byBoundaryTensor2FieldId);
}

DOUBLE CActionGaugePlaquettePSU3WithBoundary::CalculateEnergyNow(const CFieldGaugeSU3* pGauge, const CFieldTensor2Z3* pBoundary)
{
    preparethread;
    appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelEnergyPSU3WithBoundary<deviceSU3>, block, threads,
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        pBoundary->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId],
        m_fBetaOverNR,
        _HC_Volume,
        _D_RealThreadBuffer);
#else
    _LAUNCH_KERNEL(_kernelEnergyPSU3WithBoundary<deviceSU3>, block, threads,
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        pBoundary->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId],
        m_fBetaOverNR,
        _HC_Volume,
        _D_RealThreadBuffer);
#endif

    m_fNewEnergy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
    return m_fNewEnergy;
}

DOUBLE CActionGaugePlaquettePSU3WithBoundary::Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* stapleFields)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

    const INT gi = FindGaugeIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
    const INT ti = FindTensor2IndexById(tensor2Num, tensor2Fields, m_byBoundaryTensor2FieldId);
    if (gi < 0 || ti < 0)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary::Energy: cannot find gauge/boundary field (gi=%d, ti=%d).\n"), gi, ti);
        _FAIL_EXIT;
        return 0.0;
    }

    const CFieldGaugeSU3* pGauge = dynamic_cast<const CFieldGaugeSU3*>(gaugeFields[gi]);
    const CFieldTensor2Z3* pBoundary = dynamic_cast<const CFieldTensor2Z3*>(tensor2Fields[ti]);
    if (NULL == pGauge || NULL == pBoundary)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary::Energy: wrong field types (gauge or boundary).\n"));
        _FAIL_EXIT;
        return 0.0;
    }

    return CalculateEnergyNow(pGauge, pBoundary);
}

void CActionGaugePlaquettePSU3WithBoundary::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 != uiUpdateIterate)
    {
        return;
    }

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    if (NULL == pGaugeSU3 || NULL == m_pBoundaryField)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary::PrepareForHMC: gauge/boundary field missing.\n"));
        _FAIL_EXIT;
        return;
    }

    m_fLastEnergy = CalculateEnergyNow(pGaugeSU3, m_pBoundaryField);
}

UBOOL CActionGaugePlaquettePSU3WithBoundary::CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const
{
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    if (NULL == pGaugeSU3 || NULL == pForceSU3)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary only works with SU3 gauge fields.\n"));
        return FALSE;
    }
    if (NULL == m_pBoundaryField)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary::CalculateForce: boundary field not initialized.\n"));
        return FALSE;
    }

    preparethreadDir;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelForcePSU3WithBoundary<deviceSU3>, block, threads,
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        m_pBoundaryField->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[pGaugeSU3->m_byFieldId],
        m_fBetaOverNR,
        _HC_Volume,
        pForceSU3->m_pDeviceData);
#else
    _LAUNCH_KERNEL(_kernelForcePSU3WithBoundary<deviceSU3>, block, threads,
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        m_pBoundaryField->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[pGaugeSU3->m_byFieldId],
        m_fBetaOverNR,
        _HC_Volume,
        pForceSU3->m_pDeviceData);
#endif
    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

void CActionGaugePlaquettePSU3WithBoundary::SweepZ3(CFieldGaugeSU3* pGauge, CFieldTensor2Z3* pBoundary)
{
    for (UINT i = 0; i < m_uiZ3SweepsPerTrajectory; ++i)
    {
        if (m_bAllowMonopole)
        {
            SweepZ3Unconstrained(pGauge, pBoundary);
        }
        else
        {
            for (BYTE alpha = 0; alpha < _DC_Dir; ++alpha)
            {
                SweepZ3LinkStar(pGauge, pBoundary, alpha, TRUE);
                SweepZ3LinkStar(pGauge, pBoundary, alpha, FALSE);
            }
        }
    }

    if (!m_bAllowMonopole)
    {
        for (UINT i = 0; i < m_uiGlobalTwistSweepsPerTrajectory; ++i)
        {
            for (BYTE pi = 0; pi < 6; ++pi)
            {
                SweepZ3GlobalTwist(pGauge, pBoundary, pi);
            }
        }
    }

    if (m_bCheckMonopole)
    {
        const UINT uiMonopoles = CountMonopoles(pBoundary);
        appGeneral(_T("CActionGaugePlaquettePSU3WithBoundary: post-sweep monopole count = %u (AllowMonopole=%d)\n"), uiMonopoles, m_bAllowMonopole ? 1 : 0);
        if (!m_bAllowMonopole && 0 != uiMonopoles)
        {
            appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary: monopole appeared in flat mode (dB != 0)!\n"));
            _FAIL_EXIT;
            return;
        }
    }
}

void CActionGaugePlaquettePSU3WithBoundary::SweepZ3Unconstrained(CFieldGaugeSU3* pGauge, CFieldTensor2Z3* pBoundary)
{
    preparethread;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelSweepZ3PlaquetteUnconstrained<deviceSU3>, block, threads,
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        pBoundary->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId],
        m_fBetaOverNR,
        _HC_Volume);
#else
    _LAUNCH_KERNEL(_kernelSweepZ3PlaquetteUnconstrained<deviceSU3>, block, threads,
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        pBoundary->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId],
        m_fBetaOverNR,
        _HC_Volume);
#endif
    checkCudaErrors(cudaDeviceSynchronize());
}

void CActionGaugePlaquettePSU3WithBoundary::SweepZ3LinkStar(CFieldGaugeSU3* pGauge, CFieldTensor2Z3* pBoundary, BYTE byAlpha, UBOOL bEven)
{
    preparethread;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelSweepZ3LinkStar<deviceSU3>, block, threads,
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        pBoundary->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId],
        byAlpha, bEven,
        m_fBetaOverNR,
        _HC_Volume);
#else
    _LAUNCH_KERNEL(_kernelSweepZ3LinkStar<deviceSU3>, block, threads,
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        pBoundary->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId],
        byAlpha, bEven,
        m_fBetaOverNR,
        _HC_Volume);
#endif
    checkCudaErrors(cudaDeviceSynchronize());
}

void CActionGaugePlaquettePSU3WithBoundary::SweepZ3GlobalTwist(CFieldGaugeSU3* pGauge, CFieldTensor2Z3* pBoundary, BYTE byPlaqIndex)
{
    // orientation pair of plaqIndex 0..5: (0,1),(0,2),(0,3),(1,2),(1,3),(2,3)
    const BYTE byMu = byPlaqIndex <= 2 ? 0 : (byPlaqIndex <= 4 ? 1 : 2);
    const BYTE byNu = byPlaqIndex <= 2 ? static_cast<BYTE>(byPlaqIndex + 1) : (byPlaqIndex == 3 ? 2 : (byPlaqIndex == 4 ? 3 : 3));
    const UINT uiFixedMu = m_aiTwistSheetPosition[byMu] % (_HC_Lx + byMu);
    const UINT uiFixedNu = m_aiTwistSheetPosition[byNu] % (_HC_Lx + byNu);

    DOUBLE dW[3] = { 0.0, 0.0, 0.0 };

    preparethread;
    for (UINT q = 0; q < 3; ++q)
    {
        appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);
#if !_CLG_ASSUME_SQUARE_LATTICE
        _LAUNCH_KERNEL(_kernelEvaluateSheetDelta<deviceSU3>, block, threads,
            pGauge->m_byFieldId,
            pGauge->m_pDeviceData,
            pBoundary->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId],
            m_fBetaOverNR,
            byPlaqIndex, byMu, byNu, uiFixedMu, uiFixedNu, q,
            _HC_Volume,
            _D_RealThreadBuffer);
#else
        _LAUNCH_KERNEL(_kernelEvaluateSheetDelta<deviceSU3>, block, threads,
            pGauge->m_byFieldId,
            pGauge->m_pDeviceData,
            pBoundary->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId],
            m_fBetaOverNR,
            byPlaqIndex, byMu, byNu, uiFixedMu, uiFixedNu, q,
            _HC_Volume,
            _D_RealThreadBuffer);
#endif
        dW[q] = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
    }

    // sample q on a single device thread, consuming one uniform random
    UINT* pDeviceQ = NULL;
    checkCudaErrors(cudaMalloc((void**)&pDeviceQ, sizeof(UINT)));
    // stream index 0 (link 0) is always a legal random table index
    const UINT uiSeed = 0;
    _kernelSampleSheetQ<<<1, 1>>>(dW[0], dW[1], dW[2], uiSeed, pDeviceQ);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    UINT qSel = 0;
    checkCudaErrors(cudaMemcpy(&qSel, pDeviceQ, sizeof(UINT), cudaMemcpyDeviceToHost));

    // apply
    _LAUNCH_KERNEL(_kernelApplySheet, block, threads,
        pBoundary->m_pDeviceData,
        byPlaqIndex, byMu, byNu, uiFixedMu, uiFixedNu, qSel,
        _HC_Volume);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    checkCudaErrors(cudaFree(pDeviceQ));
}

UINT CActionGaugePlaquettePSU3WithBoundary::CountMonopoles(const CFieldTensor2Z3* pBoundary) const
{
    preparethread;
    UINT* pDeviceCount = NULL;
    checkCudaErrors(cudaMalloc((void**)&pDeviceCount, sizeof(UINT) * _HC_Volume));

#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelCountMonopole, block, threads,
        pBoundary->m_pDeviceData, _HC_Volume, pDeviceCount);
#else
    _LAUNCH_KERNEL(_kernelCountMonopole, block, threads,
        pBoundary->m_pDeviceData, _HC_Volume, pDeviceCount);
#endif
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    UINT* pHostCount = new UINT[_HC_Volume];
    checkCudaErrors(cudaMemcpy(pHostCount, pDeviceCount, sizeof(UINT) * _HC_Volume, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(pDeviceCount));

    UINT uiTotal = 0;
    for (UINT i = 0; i < _HC_Volume; ++i)
    {
        uiTotal += pHostCount[i];
    }
    delete[] pHostCount;
    return uiTotal;
}

void CActionGaugePlaquettePSU3WithBoundary::OnFinishTrajectory(UBOOL bWillBeAccept, INT gaugeNum, INT bosonNum, INT tensor2Num, CFieldGauge* const* gaugeFields, CFieldBoson* const* bosonFields, CFieldTensor2* const* tensor2Fields)
{
    CFieldGaugeSU3* pGauge = NULL;
    CFieldTensor2Z3* pBoundary = NULL;

    if (bWillBeAccept)
    {
        const INT gi = FindGaugeIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
        const INT ti = FindTensor2IndexById(tensor2Num, tensor2Fields, m_byBoundaryTensor2FieldId);
        pGauge = (gi >= 0) ? dynamic_cast<CFieldGaugeSU3*>(gaugeFields[gi]) : NULL;
        pBoundary = (ti >= 0) ? dynamic_cast<CFieldTensor2Z3*>(tensor2Fields[ti]) : NULL;
    }
    else
    {
        // rejected: sweep the lattice current state (U_old, B_old)
        pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(m_byGaugeFieldIds[0]));
        pBoundary = const_cast<CFieldTensor2Z3*>(m_pBoundaryField);
    }

    if (NULL == pGauge || NULL == pBoundary)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary::OnFinishTrajectory: cannot resolve the (U,B) state to sweep.\n"));
        _FAIL_EXIT;
        return;
    }

    SweepZ3(pGauge, pBoundary);
    checkCudaErrors(cudaDeviceSynchronize());

    m_fPostSweepEnergy = CalculateEnergyNow(pGauge, pBoundary);
    m_bPostSweepEnergyValid = TRUE;
}

void CActionGaugePlaquettePSU3WithBoundary::OnFinishTrajectory(UBOOL bAccepted)
{
    UN_USE(bAccepted);

    if (!m_bPostSweepEnergyValid)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3WithBoundary::OnFinishTrajectory: no valid post-sweep energy recorded.\n"));
        _FAIL_EXIT;
        return;
    }

    m_fLastEnergy = m_fPostSweepEnergy;
    m_bPostSweepEnergyValid = FALSE;
}

CCString CActionGaugePlaquettePSU3WithBoundary::GetInfos(const CCString& tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("PSU3WithBoundary Fundamental-Lift Plaquette Action (Z3 2-form B)\n");
    sRet = sRet + tab + _T("BoundaryFieldId : ") + appToString(m_byBoundaryTensor2FieldId) + _T("\n");
    sRet = sRet + tab + _T("AllowMonopole : ") + appToString(m_bAllowMonopole ? 1 : 0) + _T("\n");
    sRet = sRet + tab + _T("Z3Sweeps : ") + appToString(m_uiZ3SweepsPerTrajectory) + _T("\n");
    sRet = sRet + tab + _T("GlobalTwistSweeps : ") + appToString(m_uiGlobalTwistSweepsPerTrajectory) + _T("\n");
    return sRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
