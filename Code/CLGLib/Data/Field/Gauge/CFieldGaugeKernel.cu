//=============================================================================
// FILENAME : CFieldGaugeKernel.cu
// 
// DESCRIPTION:
//
// REVISION:
//  [mm/dd/yy]
//  [07/24/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldGaugeKernel.h"
#include "CFieldGaugeLink.h"
#include "Tools/Math/DeviceTemplates/DeviceInlineStaggeredRotation.h"

__BEGIN_NAMESPACE

#pragma region gauge field

#pragma region Kernels

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteGaugeCacheIndex(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceGauge* pStapleData, //can be NULL
    deviceGauge* pForceData,
    Real betaOverN)
{
    intokernalDir_NoDir;

    //Real test_force = F(0.0);
    betaOverN = betaOverN * F(-0.5);
#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1;
#endif

    deviceGauge res = _makeZero<deviceGauge>();

    //there are 6 staples, each is sum of two plaquttes
    for (INT i = 0; i < plaqCountPerLink; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLengthm1 + uiLinkIndex * plaqCountAllLink];
        deviceGauge toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }

        for (INT j = 1; j < plaqLengthm1; ++j)
        {
            SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + uiLinkIndex * plaqCountAllLink];
            const deviceGauge& toMul = pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)];

            if (nextlink.NeedToDagger())
            {
                _muldag(toAdd, toMul);
            }
            else
            {
                _mul(toAdd, toMul);
            }
        }
        _add(res, toAdd);
    }
    if (NULL != pStapleData)
    {
        pStapleData[uiLinkIndex] = res;
    }

    //staple calculated
    //deviceGauge force = pDeviceData[uiLinkIndex];
    //_muldag(force, res);
    //_ta(force);
    _mul(res, betaOverN);

    //force is additive
    _add(pForceData[uiLinkIndex], res);
}


template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeForceCacheIndexAnisotropy(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceGauge* pForceData,
    DOUBLE betaOverN,
    DOUBLE xi)
{
    intokernalDir;

    betaOverN = betaOverN * -0.5;
#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1;
#endif
    const DOUBLE oneoverxi = 1.0 / xi;

    deviceGauge res = _makeZero<deviceGauge>();

    //there are 6 staples, each is sum of two plaquttes
    for (INT i = 0; i < plaqCountPerLink; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLengthm1 + uiLinkIndex * plaqCountAllLink];
        deviceGauge toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }

        for (INT j = 1; j < plaqLengthm1; ++j)
        {
            SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + uiLinkIndex * plaqCountAllLink];
            const deviceGauge& toMul = pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)];

            if (nextlink.NeedToDagger())
            {
                _muldag(toAdd, toMul);
            }
            else
            {
                _mul(toAdd, toMul);
            }
        }
        
        if (3 == dir || i > 3)
        {
            _mul(toAdd, static_cast<Real>(betaOverN * oneoverxi));
        }
        else
        {
            _mul(toAdd, static_cast<Real>(betaOverN * xi));
        }
        _add(res, toAdd);
    }

    //force is additive
    _add(pForceData[uiLinkIndex], res);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateOnlyStapleGauge(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceGauge* pStapleData)
{
    intokernalDir_NoDir;

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1;
#endif

    deviceGauge res = _makeZero<deviceGauge>();

    //there are 6 staples, each is sum of two plaquttes
    for (INT i = 0; i < plaqCountPerLink; ++i)
    {
        const SIndex& first = pCachedIndex[i * plaqLengthm1 + uiLinkIndex * plaqCountAllLink];
        deviceGauge toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }

        for (INT j = 1; j < plaqLengthm1; ++j)
        {
            const SIndex& nextlink = pCachedIndex[i * plaqLengthm1 + j + uiLinkIndex * plaqCountAllLink];
            const deviceGauge& toMul = pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)];

            if (nextlink.NeedToDagger())
            {
                _muldag(toAdd, toMul);
            }
            else
            {
                _mul(toAdd, toMul);
            }
        }
        _add(res, toAdd);
    }
    pStapleData[uiLinkIndex] = res;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAllStaplesGauge(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceGauge* const * __restrict__ ppStapleData)
{
    intokernalDir_NoDir;

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1;
#endif

    //there are 6 staples, each is sum of two plaquttes, plaqCount should be 12
    for (INT i = 0; i < plaqCountPerLink; ++i)
    {
        const SIndex& first = pCachedIndex[i * plaqLengthm1 + uiLinkIndex * plaqCountAllLink];
        deviceGauge toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }

        for (INT j = 1; j < plaqLengthm1; ++j)
        {
            const SIndex& nextlink = pCachedIndex[i * plaqLengthm1 + j + uiLinkIndex * plaqCountAllLink];
            const deviceGauge& toMul = pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)];

            if (nextlink.NeedToDagger())
            {
                _muldag(toAdd, toMul);
            }
            else
            {
                _mul(toAdd, toMul);
            }
        }
        ppStapleData[i][uiLinkIndex] = toAdd;
    }
}

/**
* 
* index is only cached for +mu, +nu, not mu, -nu, so we do not use cached index
* need to cache:
* (x,y), (x,-y), (-x, y), (-x, -y)
* (x,z), (x,-z), (-x, z), (-x, -z)
* (x,t), (x,-t), (-x, t), (-x, -t)
* 
* (y,z), (y,-z), (-y, z), (-y, -z)
* (y,t), (y,-t), (-y, t), (-y, -t)
* 
* (z,t), (z,-t), (-z, t), (-z, -t)
* 
* link index mu to mimic the signs, 0:00, 1:01, 2:10, 3:11
* so:
* mu,    0,        1,        2         3         4         5
* 0   ( x, y)   ( x, z)   ( x, t)   ( y, z)   ( y, t)   ( z, t)
* 1   (-x, y)   (-x, z)   (-x, t)   (-y, z)   (-y, t)   (-z, t)
* 2   ( x,-y)   ( x,-z)   ( x,-t)   ( y,-z)   ( y,-t)   ( z,-t)
* 3   (-x,-y)   (-x,-z)   (-x,-t)   (-y,-z)   (-y,-t)   (-z,-t)
* 
* for example (-t, y) is (y, -t)^dagger
* 
* 
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAllPlaquttesGauge(
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* const* __restrict__ ppPlaqutteData,
    BYTE byFieldId)
{
    intokernalE(24);
    const BYTE sign = static_cast<BYTE>(elementIdx & 3U);
    const BYTE munu = elementIdx >> 2U;
    const SCHAR mu = _plaq_idx[munu][0];
    const SCHAR nu = _plaq_idx[munu][1];

    const SCHAR sign1 = 1 - static_cast<SCHAR>((sign & 1) << 1);
    const SCHAR sign2 = 1 - static_cast<SCHAR>(sign & 2);

    SCHAR path[4] = { 
        static_cast<SCHAR>(mu * sign1),
        static_cast<SCHAR>(nu * sign2),
        static_cast<SCHAR>(-mu * sign1),
        static_cast<SCHAR>(-nu * sign2)
    };
    ppPlaqutteData[munu][_deviceGetLinkIndex(uiSiteIndex, sign)] = _deviceLinkT(pDeviceData, __deviceSiteIndexToInt4(uiSiteIndex), 4, byFieldId, path);
}

/**
* Fmunu = (Pmunu - Pmunu^+) / 8i, where sum is mu<nu, Ocl = csw * (-1/4) * sigmamunu Fmunu (arXiv:1311.6312)
* 
* Note that, if we use sum_{mu>nu}, we need also make Fmunu traceless
*
* Pmunu  = + (mu, nu) + (nu, -mu) - (-nu, -mu) - (mu, -nu) 10.1016/0550-3213(85)90002-1
* Pmunu+ = + (nu, mu) + (-mu, nu) - (-mu, -nu) - (-nu, mu)
* Pmunu - Pmunu+  = Qmunu - Qmunu+
* Qmunu = + (mu, nu) - (-mu, nu) + (-mu, -nu) - (mu, -nu)
* 
* Another notation:
* Pmunu - Pnumu, with 
* Pmunu =  (mu, nu) + (nu, -mu) + (-mu, -nu) + (-nu, mu)
*-Pnumu = -(nu, mu) - (mu, -nu) - (-nu, -mu) - (-mu, nu)
* 
* So, it is Qmunu  = (mu, nu) - (mu, -nu) - (-mu, nu) + (-mu, -nu)
*          -Qmunu+ =-(nu, mu) + (-nu, mu) + (nu, -mu) - (-nu, -mu)
* 
*
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateFmunuGauge(
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pFmunu,
    BYTE byFieldId)
{
    intokernalE(24);
    const BYTE sign = static_cast<BYTE>(elementIdx & 3U);
    const BYTE munu = elementIdx >> 2U;
    const SCHAR mu = _plaq_idx[munu][0];
    const SCHAR nu = _plaq_idx[munu][1];

    const SCHAR sign1 = 1 - static_cast<SCHAR>((sign & 1) << 1);
    const SCHAR sign2 = 1 - static_cast<SCHAR>(sign & 2);

    //sign = 0:  mu  nu    sign1 = 1  sign2 = 1
    //       1: -mu  nu    sign1 =-1  sign2 = 1
    //       2:  mu -nu    sign1 = 1  sign2 =-1
    //       3: -mu -nu    sign1 =-1  sign2 =-1

    SCHAR path[4] = {
        static_cast<SCHAR>(mu * sign1),
        static_cast<SCHAR>(nu * sign2),
        static_cast<SCHAR>(-mu * sign1),
        static_cast<SCHAR>(-nu * sign2)
    };

    deviceGauge toAdd = _deviceLinkT(pDeviceData, __deviceSiteIndexToInt4(uiSiteIndex), 4, byFieldId, path);

    const UINT uiIdx = munu * _DC_Volume + uiSiteIndex;
    if (0 == sign)
    {
        pFmunu[uiIdx] = toAdd;
    }
    __syncthreads();

    if (1 == sign)
    {
        _sub(pFmunu[uiIdx], toAdd);
    }
    __syncthreads();

    if (2 == sign)
    {
        _sub(pFmunu[uiIdx], toAdd);
    }
    __syncthreads();

    if (3 == sign)
    {
        _add(pFmunu[uiIdx], toAdd);

        //_ta(pFmunu[uiIdx]);
        _sub(pFmunu[uiIdx], _daggerC(pFmunu[uiIdx]));
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyGaugeCacheIndex(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerSite,
#endif
    DOUBLE betaOverN,
    DOUBLE* results
)
{
    intokernalE(plaqCountPerSite);

    const UINT indexSkip = plaqCountPerSite * plaqLength * uiSiteIndex;

    SIndex first = pCachedIndex[elementIdx * plaqLength + indexSkip];
    deviceGauge toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

    if (first.NeedToDagger())
    {
        _dagger(toAdd);
    }

    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedIndex[elementIdx * plaqLength + j + indexSkip];
        const deviceGauge& toMul = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];
        if (first.NeedToDagger())
        {
            _muldag(toAdd, toMul);
        }
        else
        {
            _mul(toAdd, toMul);
        }
    }

    const DOUBLE res = (_dim<deviceGauge>() - _retr(toAdd));
    for (BYTE i = 0; i < plaqCountPerSite; ++i)
    {
        if (i == elementIdx)
        {
            if (0 == i)
            {
                results[uiSiteIndex] = res;
            }
            else
            {
                results[uiSiteIndex] += res;
            }
        }
        __syncthreads();
    }

    if (0 == elementIdx)
    {
        results[uiSiteIndex] = results[uiSiteIndex] * betaOverN;
    }
}

/**
* the plaquttes from 0 to 5 are:
* xy
* xz
* xt
* yz
* yt
* zt
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyGaugeCacheIndexAnisotropy(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerSite,
#endif
    DOUBLE betaOverN,
    DOUBLE xi,
    DOUBLE* results
)
{
    intokernalE(plaqCountPerSite);

    const UINT indexSkip = plaqCountPerSite * plaqLength * uiSiteIndex;

    const DOUBLE oneoverxi = 1.0 / xi;

    SIndex first = pCachedIndex[elementIdx * plaqLength + indexSkip];
    deviceGauge toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

    if (first.NeedToDagger())
    {
        _dagger(toAdd);
    }

    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedIndex[elementIdx * plaqLength + j + indexSkip];
        const deviceGauge& toMul = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];
        if (first.NeedToDagger())
        {
            _muldag(toAdd, toMul);
        }
        else
        {
            _mul(toAdd, toMul);
        }
    }

    const DOUBLE res = (_dim<deviceGauge>() - _retr(toAdd));
    for (BYTE i = 0; i < plaqCountPerSite; ++i)
    {
        if (i == elementIdx)
        {
            //0, 1, 3 times xi
            //2, 4, 5 times 1/xi
            if (0 == i)
            {
                results[uiSiteIndex] = res * xi;
            }
            else
            {
                
                if (i & (i + 1))
                {
                    //2, 4, 5
                    results[uiSiteIndex] += res * oneoverxi;
                }
                else
                {
                    //0, 1, 3
                    results[uiSiteIndex] += res * xi;
                }
            }
        }
        __syncthreads();
    }

    if (0 == elementIdx)
    {
        results[uiSiteIndex] = results[uiSiteIndex] * betaOverN;
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyGauge_UseClover(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE stapleConstant,
    DOUBLE fBetaOverN,
    DOUBLE* results
)
{
    intokernalDirInt4;

    DOUBLE fRes = 0.0;
    for (BYTE byDir2 = dir + 1; byDir2 < _DC_Dir; ++byDir2)
    {
        fRes += _deviceCloverRetrT(pDeviceData, sSite4, __bi(sSite4), dir, byDir2, byFieldId);
    }
    for (BYTE i = 0; i < _DC_Dir; ++i)
    {
        if (i == dir)
        {
            if (0 == i)
            {
                results[uiSiteIndex] = fRes;
            }
            else
            {
                results[uiSiteIndex] += fRes;
            }
        }
        __syncthreads();
    }

    if (0 == dir)
    {
        results[uiSiteIndex] = (stapleConstant - 0.25 * results[uiSiteIndex]) * fBetaOverN;
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyGauge_UseCloverAnisotropy(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE fBetaOverN,
    DOUBLE xi,
    DOUBLE* results
)
{
    intokernalDirInt4;

    const DOUBLE oneoverxi = 1.0 / xi;

    DOUBLE fRes = 0.0;
    for (BYTE byDir2 = dir + 1; byDir2 < _DC_Dir; ++byDir2)
    {
        if (3 == dir || 3 == byDir2)
        {
            fRes += (_dim<deviceGauge>() - 0.25 * _deviceCloverRetrT(pDeviceData, sSite4, __bi(sSite4), dir, byDir2, byFieldId)) * oneoverxi;
        }
        else
        {
            fRes += (_dim<deviceGauge>() - 0.25 * _deviceCloverRetrT(pDeviceData, sSite4, __bi(sSite4), dir, byDir2, byFieldId)) * xi;
        }
    }
    for (BYTE i = 0; i < _DC_Dir; ++i)
    {
        if (i == dir)
        {
            if (0 == i)
            {
                results[uiSiteIndex] = fRes;
            }
            else
            {
                results[uiSiteIndex] += fRes;
            }
        }
        __syncthreads();
    }

    if (0 == dir)
    {
        results[uiSiteIndex] = results[uiSiteIndex] * fBetaOverN;
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyUsingStapleGauge(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pStapleData,
    DOUBLE stapleConstant,
    DOUBLE betaOverN,
    DOUBLE* results
)
{
    intokernalDir;
    const DOUBLE res = (stapleConstant - _retr(_muldagC(pDeviceData[uiLinkIndex], pStapleData[uiLinkIndex])));
    for (BYTE i = 0; i < _DC_Dir; ++i)
    {
        if (i == dir)
        {
            if (0 == i)
            {
                results[uiSiteIndex] = res;
            }
            else
            {
                results[uiSiteIndex] += res;
            }
        }
        __syncthreads();
    }

    if (0 == dir)
    {
        results[uiSiteIndex] = results[uiSiteIndex] * betaOverN * 0.25;
    }
    
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteCacheIndexT_D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceGauge* pStapleData, //can be NULL
    deviceGauge* pForceData,
    Real betaOverN)
{
    intokernalDir_NoDir;

    //Real test_force = F(0.0);
    betaOverN = betaOverN * F(-0.5);
#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1;
#endif

    deviceGauge res = _makeZero<deviceGauge>();

    //there are 6 staples, each is sum of two plaquttes
    for (BYTE i = 0; i < plaqCountPerLink; ++i)
    {
        BYTE diricCount = 0;
        const SIndex& first = pCachedIndex[i * plaqLengthm1 + uiLinkIndex * plaqCountAllLink];
        if (first.IsDirichlet())
        {
            ++diricCount;
        }
        //deviceGauge toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
        deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }

        for (BYTE j = 1; j < plaqLengthm1; ++j)
        {
            const SIndex& nextlink = pCachedIndex[i * plaqLengthm1 + j + uiLinkIndex * plaqCountAllLink];
            if (nextlink.IsDirichlet())
            {
                ++diricCount;
            }
            //deviceGauge toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);
            deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink));

            if (nextlink.NeedToDagger())
            {
                _muldag(toAdd, toMul);
            }
            else
            {
                _mul(toAdd, toMul);
            }
        }
        if (diricCount < plaqLength - 1)
        {
            // If more than 3(including 3) of the edges are Dirichlet, 
            // the plaqutte dose NOT exist.
            _add(res, toAdd);
        }
    }
    if (NULL != pStapleData)
    {
        pStapleData[uiLinkIndex] = res;
    }

    //staple calculated
    //deviceGauge force(pDeviceData[uiLinkIndex]);
    //_muldag(force, res);
    //test_force += F(-2.0) * betaOverN * __SU3Generators[8].MulC(force).ImTr();
    //_ta(force);
    _mul(res, betaOverN);

    //force is additive
    _add(pForceData[uiLinkIndex], res);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelForceAtSiteCacheIndexAnisotropyT_D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceGauge* pForceData,
    DOUBLE betaOverN,
    DOUBLE xi)
{
    intokernalDir;

    //Real test_force = F(0.0);
    betaOverN = betaOverN * -0.5;
    const DOUBLE oneoverxi = 1.0 / xi;
#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1;
#endif

    deviceGauge res = _makeZero<deviceGauge>();

    //there are 6 staples, each is sum of two plaquttes
    for (BYTE i = 0; i < plaqCountPerLink; ++i)
    {
        BYTE diricCount = 0;
        const SIndex& first = pCachedIndex[i * plaqLengthm1 + uiLinkIndex * plaqCountAllLink];
        if (first.IsDirichlet())
        {
            ++diricCount;
        }
        //deviceGauge toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
        deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }

        for (BYTE j = 1; j < plaqLengthm1; ++j)
        {
            const SIndex& nextlink = pCachedIndex[i * plaqLengthm1 + j + uiLinkIndex * plaqCountAllLink];
            if (nextlink.IsDirichlet())
            {
                ++diricCount;
            }
            //deviceGauge toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);
            deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink));

            if (nextlink.NeedToDagger())
            {
                _muldag(toAdd, toMul);
            }
            else
            {
                _mul(toAdd, toMul);
            }
        }
        if (diricCount < plaqLength - 1)
        {
            // If more than 3(including 3) of the edges are Dirichlet, 
            // the plaqutte dose NOT exist.
            if (3 == dir || i > 3)
            {
                _mul(toAdd, static_cast<Real>(betaOverN * oneoverxi));
            }
            else
            {
                _mul(toAdd, static_cast<Real>(betaOverN * xi));
            }
            _add(res, toAdd);
        }
    }

    //force is additive
    _add(pForceData[uiLinkIndex], res);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteCacheIndexCloverT_D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceGauge* pStapleData, //can be NULL
    deviceGauge* pForceData,
    Real betaOverN)
{
    intokernalDir_NoDir;

    //Real test_force = F(0.0);
    betaOverN = betaOverN * F(-0.5);
#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1;
#endif

    deviceGauge res = _makeZero<deviceGauge>();
    //there are 6 staples, each is sum of two plaquttes
    for (BYTE i = 0; i < plaqCountPerLink; ++i)
    {
        BYTE diricCount = 0;
        const SIndex& first = pCachedIndex[i * plaqLengthm1 + uiLinkIndex * plaqCountAllLink];
        if (first.IsOutside())
        {
            ++diricCount;
        }
        
        //deviceGauge toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
        deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }

        for (BYTE j = 1; j < plaqLengthm1; ++j)
        {
            const SIndex& nextlink = pCachedIndex[i * plaqLengthm1 + j + uiLinkIndex * plaqCountAllLink];
            if (nextlink.IsOutside())
            {
                ++diricCount;
            }
            //deviceGauge toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);
            deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink));

            if (nextlink.NeedToDagger())
            {
                _muldag(toAdd, toMul);
            }
            else
            {
                _mul(toAdd, toMul);
            }
        }
        if (diricCount < plaqLengthm1)
        {
            // If more than 3(including 3) of the edges are Dirichlet, 
            // the plaqutte dose NOT exist.
            if (1 == diricCount)
            {
                _mul(toAdd, F(0.5));
            }
            else if (2 == diricCount)
            {
                _mul(toAdd, F(0.25));
            }
            _add(res, toAdd);
        }
    }
    if (NULL != pStapleData)
    {
        pStapleData[uiLinkIndex] = res;
    }

    //staple calculated
    //deviceGauge force(pDeviceData[uiLinkIndex]);
    //_muldag(force, res);
    //test_force += F(-2.0) * betaOverN * __SU3Generators[8].MulC(force).ImTr();
    //_ta(force);
    _mul(res, betaOverN);

    //force is additive
    _add(pForceData[uiLinkIndex], res);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelForceAtSiteCacheIndexCloverAnisotropyT_D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceGauge* pForceData,
    DOUBLE betaOverN,
    DOUBLE xi)
{
    intokernalDir;

    //Real test_force = F(0.0);
    betaOverN = betaOverN * -0.5;
    const DOUBLE oneoverxi = 1.0 / xi;
#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1;
#endif

    deviceGauge res = _makeZero<deviceGauge>();
    //there are 6 staples, each is sum of two plaquttes
    for (BYTE i = 0; i < plaqCountPerLink; ++i)
    {
        BYTE diricCount = 0;
        const SIndex& first = pCachedIndex[i * plaqLengthm1 + uiLinkIndex * plaqCountAllLink];
        if (first.IsOutside())
        {
            ++diricCount;
        }

        //deviceGauge toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
        deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }

        for (BYTE j = 1; j < plaqLengthm1; ++j)
        {
            const SIndex& nextlink = pCachedIndex[i * plaqLengthm1 + j + uiLinkIndex * plaqCountAllLink];
            if (nextlink.IsOutside())
            {
                ++diricCount;
            }
            //deviceGauge toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);
            deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink));

            if (nextlink.NeedToDagger())
            {
                _muldag(toAdd, toMul);
            }
            else
            {
                _mul(toAdd, toMul);
            }
        }
        if (diricCount < plaqLengthm1)
        {
            // If more than 3(including 3) of the edges are Dirichlet, 
            // the plaqutte dose NOT exist.
            if (3 == dir || i > 3)
            {
                _mul(toAdd, static_cast<Real>(betaOverN * oneoverxi));
            }
            else
            {
                _mul(toAdd, static_cast<Real>(betaOverN * xi));
            }

            if (1 == diricCount)
            {
                _mul(toAdd, F(0.5));
            }
            else if (2 == diricCount)
            {
                _mul(toAdd, F(0.25));
            }
            _add(res, toAdd);
        }
    }

    //staple calculated
    //deviceGauge force(pDeviceData[uiLinkIndex]);
    //_muldag(force, res);
    //test_force += F(-2.0) * betaOverN * __SU3Generators[8].MulC(force).ImTr();
    //_ta(force);
    //_mul(res, betaOverN);

    //force is additive
    _add(pForceData[uiLinkIndex], res);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyCacheIndexT_D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerSite,
#endif
    DOUBLE betaOverN,
    DOUBLE* results
)
{
    intokernalE(plaqCountPerSite);

#if !_CLG_ASSUME_SQUARE_LATTICE
    UINT plaqCountAllSite = plaqCountPerSite * plaqLength;
#endif

    SIndex first = pCachedIndex[elementIdx * plaqLength + uiSiteIndex * plaqCountAllSite];
    deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

    if (first.NeedToDagger())
    {
        _dagger(toAdd);
    }

    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedIndex[elementIdx * plaqLength + j + uiSiteIndex * plaqCountAllSite];
        deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
        if (first.NeedToDagger())
        {
            _muldag(toAdd, toMul);
        }
        else
        {
            _mul(toAdd, toMul);
        }
    }
    const DOUBLE res = (_dim<deviceGauge>() - _retr(toAdd));
    for (BYTE i = 0; i < plaqCountPerSite; ++i)
    {
        if (i == elementIdx)
        {
            if (0 == i)
            {
                results[uiSiteIndex] = res;
            }
            else
            {
                results[uiSiteIndex] += res;
            }
        }
        __syncthreads();
    }

    if (0 == elementIdx)
    {
        results[uiSiteIndex] = results[uiSiteIndex] * betaOverN;
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyCacheIndexAnisotropyT_D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerSite,
#endif
    DOUBLE betaOverN,
    DOUBLE xi,
    DOUBLE* results
)
{
    intokernalE(plaqCountPerSite);

#if !_CLG_ASSUME_SQUARE_LATTICE
    UINT plaqCountAllSite = plaqCountPerSite * plaqLength;
#endif

    const DOUBLE oneoverxi = 1.0 / xi;

    SIndex first = pCachedIndex[elementIdx * plaqLength + uiSiteIndex * plaqCountAllSite];
    deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

    if (first.NeedToDagger())
    {
        _dagger(toAdd);
    }

    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedIndex[elementIdx * plaqLength + j + uiSiteIndex * plaqCountAllSite];
        deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
        if (first.NeedToDagger())
        {
            _muldag(toAdd, toMul);
        }
        else
        {
            _mul(toAdd, toMul);
        }
    }
    const DOUBLE res = (_dim<deviceGauge>() - _retr(toAdd));
    for (BYTE i = 0; i < plaqCountPerSite; ++i)
    {
        if (i == elementIdx)
        {
            if (0 == i)
            {
                results[uiSiteIndex] = res * xi;
            }
            else
            {
                if (i & (i + 1))
                {
                    //2, 4, 5
                    results[uiSiteIndex] += res * oneoverxi;
                }
                else
                {
                    //0, 1, 3
                    results[uiSiteIndex] += res * xi;
                }
            }
        }
        __syncthreads();
    }

    if (0 == elementIdx)
    {
        results[uiSiteIndex] = results[uiSiteIndex] * betaOverN;
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateOnlyStapleT_D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceGauge* pStapleData)
{
    intokernalDir_NoDir;

    //Real test_force = F(0.0);
#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1;
#endif

    deviceGauge res = _makeZero<deviceGauge>();

    //there are 6 staples, each is sum of two plaquttes
    for (INT i = 0; i < plaqCountPerLink; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLengthm1 + uiLinkIndex * plaqCountAllLink];
        deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }

        for (int j = 1; j < plaqLengthm1; ++j)
        {
            SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + uiLinkIndex * plaqCountAllLink];
            deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink));

            if (nextlink.NeedToDagger())
            {
                _muldag(toAdd, toMul);
            }
            else
            {
                _mul(toAdd, toMul);
            }
        }
        _add(res, toAdd);
    }
    pStapleData[uiLinkIndex] = res;
}

#pragma endregion

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceAndStaple(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, deviceGauge* pStaple, Real betaOverN)
{
    preparethreadDir;
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId]);
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelStapleAtSiteGaugeCacheIndex<deviceGauge>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStaple,
        pForce,
        betaOverN);
#else
    _LAUNCH_KERNEL(_kernelStapleAtSiteGaugeCacheIndex<deviceGauge>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        pStaple,
        pForce,
        betaOverN);
#endif
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceAndStaple(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, deviceGauge* pStaple, CLGComplex betaOverN)
{
    appCrucial(_T("Complex version CalculateForceAndStaple for this gauge group not implemented!\n"));
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverN, DOUBLE xi)
{
    preparethreadDir;
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId]);
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelGaugeForceCacheIndexAnisotropy<deviceGauge>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pForce,
        betaOverN,
        xi);
#else
    _LAUNCH_KERNEL(_kernelGaugeForceCacheIndexAnisotropy<deviceGauge>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        pForce,
        betaOverN,
        xi);
#endif
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateOnlyStaple(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pStaple)
{
    preparethreadDir;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelCalculateOnlyStapleGauge<deviceGauge>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStaple);
#else
    _LAUNCH_KERNEL(_kernelCalculateOnlyStapleGauge<deviceGauge>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        pStaple);
#endif
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateAllStaples(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* const * ppStaples)
{
    preparethreadDir;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelCalculateAllStaplesGauge<deviceGauge>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        ppStaples);
#else
    _LAUNCH_KERNEL(_kernelCalculateAllStaplesGauge<deviceGauge>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        ppStaples);
#endif
    _CHECKCUDA;
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateAllPlaquttes(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* const* ppPlaquettes)
{
    preparethreadE(24);
    _LAUNCH_KERNEL(_kernelCalculateAllPlaquttesGauge<deviceGauge>, block, threads, deviceData, ppPlaquettes, byFieldId);
    _CHECKCUDA;
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateFmunu(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pFmunu)
{
    preparethreadE(24);
    _LAUNCH_KERNEL(_kernelCalculateFmunuGauge<deviceGauge>, block, threads, deviceData, pFmunu, byFieldId);
    _CHECKCUDA;
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN)
{
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId]);

#if !_CLG_ASSUME_SQUARE_LATTICE
    preparethreadE(appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyGaugeCacheIndex<deviceGauge>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer
    );
#else
    preparethreadE(plaqCountPerSite);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyGaugeCacheIndex<deviceGauge>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId],
        betaOverN,
        _D_RealThreadBuffer
    );
#endif


    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergyAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN, DOUBLE xi)
{
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId]);

#if !_CLG_ASSUME_SQUARE_LATTICE
    preparethreadE(appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyGaugeCacheIndexAnisotropy<deviceGauge>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        xi,
        _D_RealThreadBuffer
    );
#else
    preparethreadE(plaqCountPerSite);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyGaugeCacheIndexAnisotropy<deviceGauge>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId],
        betaOverN,
        xi,
        _D_RealThreadBuffer
    );
#endif


    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergyUseClover(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN)
{
    //appGeneral(_T("const %f\n"), 3.0 * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite);
#if !_CLG_ASSUME_SQUARE_LATTICE
    const DOUBLE stpc = matrixN * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite;
#else
    const DOUBLE stpc = matrixN * plaqCountPerSite;
#endif
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyGauge_UseClover<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        stpc,
        betaOverN,
        _D_RealThreadBuffer);
    //cudaDeviceSynchronize();
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergyUseCloverAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN, DOUBLE xi)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyGauge_UseCloverAnisotropy<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        betaOverN,
        xi,
        _D_RealThreadBuffer);
    //cudaDeviceSynchronize();
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergyUsingStaple(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN, const deviceGauge* pStaple)
{
#if !_CLG_ASSUME_SQUARE_LATTICE
    const DOUBLE stpc = matrixN * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink;
#else
    const DOUBLE stpc = matrixN * plaqCountPerLink;
#endif
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyUsingStapleGauge<deviceGauge>, block, threads,
        deviceData,
        pStaple,
        stpc,
        betaOverN,
        _D_RealThreadBuffer);
    //cudaDeviceSynchronize();
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceAndStaple_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, deviceGauge* pStaple, Real betaOverN)
{
    preparethreadDir;
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId]);
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelStapleAtSiteCacheIndexT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStaple,
        pForce,
        betaOverN);
#else
    _LAUNCH_KERNEL(_kernelStapleAtSiteCacheIndexT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        pStaple,
        pForce,
        betaOverN);
#endif
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceAnisotropy_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverN, DOUBLE xi)
{
    preparethreadDir;
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId]);
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelForceAtSiteCacheIndexAnisotropyT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pForce,
        betaOverN,
        xi);
#else
    _LAUNCH_KERNEL(_kernelForceAtSiteCacheIndexAnisotropyT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        pForce,
        betaOverN,
        xi);
#endif
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceAndStapleClover_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, deviceGauge* pStaple, Real betaOverN)
{
    preparethreadDir;
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId]);
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelStapleAtSiteCacheIndexCloverT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStaple,
        pForce,
        betaOverN);
#else
    _LAUNCH_KERNEL(_kernelStapleAtSiteCacheIndexCloverT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        pStaple,
        pForce,
        betaOverN);
#endif
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceCloverAnisotropy_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverN, DOUBLE xi)
{
    preparethreadDir;
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId]);
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelForceAtSiteCacheIndexCloverAnisotropyT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pForce,
        betaOverN,
        xi);
#else
    _LAUNCH_KERNEL(_kernelForceAtSiteCacheIndexCloverAnisotropyT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        pForce,
        betaOverN,
        xi);
#endif
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateOnlyStaple_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pStaple)
{
    preparethreadDir;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelCalculateOnlyStapleT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStaple);
#else
    _LAUNCH_KERNEL(_kernelCalculateOnlyStapleT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        pStaple);
#endif
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergy_D(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN)
{
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId]);

#if !_CLG_ASSUME_SQUARE_LATTICE
    preparethreadE(appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyCacheIndexT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer
    );
#else
    preparethreadE(plaqCountPerSite);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyCacheIndexT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId],
        betaOverN,
        _D_RealThreadBuffer
    );
#endif

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergyAnisotropy_D(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN, DOUBLE xi)
{
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId]);

#if !_CLG_ASSUME_SQUARE_LATTICE
    preparethreadE(appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyCacheIndexAnisotropyT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        xi,
        _D_RealThreadBuffer
    );
#else
    preparethreadE(plaqCountPerSite);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyCacheIndexAnisotropyT_D<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId],
        betaOverN,
        xi,
        _D_RealThreadBuffer
    );
#endif

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

#pragma endregion


#pragma region rotation - related

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_Cache_XYTau_TermT(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pToBeCached,
    BYTE byGaugeFieldId)
{
    intokernalE(8);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    const UBOOL bPlusX = (0 != (elementIdx & 1));
    const UBOOL bPlusY = (0 != (elementIdx & 2));
    const UBOOL bPlusT = (0 != (elementIdx & 4));

    //same check as _kernelDFermionKS_PR_XYTau_TermT, do not access gauge links out of the Dirichlet region
    SSmallInt4 sOffset = sSite4;
    sOffset.x = sOffset.x + (bPlusX ? 1 : -1);
    sOffset.y = sOffset.y + (bPlusY ? 1 : -1);
    sOffset.w = sOffset.w + (bPlusT ? 1 : -1);
    if (__idx->m_pDeviceIndexPositionToSIndex[byGaugeFieldId][__bi(sOffset)].IsDirichlet())
    {
        //gauge field is identity across the boundary, and this entry is always skipped by the Dirichlet check of the D operator kernels
        pToBeCached[uiSiteIndex * 8 + elementIdx] = _makeId<deviceGauge>();
        return;
    }

    pToBeCached[uiSiteIndex * 8 + elementIdx] = _deviceVXYTOptimizedT(pGauge, sSite4, byGaugeFieldId, bPlusX, bPlusY, bPlusT);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_Cache_XYTermT(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pToBeCached,
    BYTE byGaugeFieldId)
{
    intokernalE(8);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    const UBOOL bPlusMu = elementIdx & 2;
    const UBOOL bPlusTau = elementIdx & 4;
    const UINT bXorY = elementIdx & 1;

    //same check as _kernelDFermionKS_PR_XYTermT, do not access gauge links out of the Dirichlet region
    const UINT bYorX = 1 - bXorY;
    SSmallInt4 sTargetSite = sSite4;
    sTargetSite.m_byData4[bYorX] = sTargetSite.m_byData4[bYorX] + (bPlusMu ? 2 : -2);
    sTargetSite.w = sTargetSite.w + (bPlusTau ? 1 : -1);
    if (__idx->m_pDeviceIndexPositionToSIndex[byGaugeFieldId][__bi(sTargetSite)].IsDirichlet())
    {
        //gauge field is identity across the boundary, and this entry is always skipped by the Dirichlet check of the D operator kernels
        pToBeCached[uiSiteIndex * 8 + elementIdx] = _makeId<deviceGauge>();
        return;
    }

    pToBeCached[uiSiteIndex * 8 + elementIdx] = _deviceVXXTauOptimizedT(pGauge, sSite4, byGaugeFieldId, bXorY, bPlusMu, bPlusTau);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_Cache_XYTau_TermEMT(
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pPhase,
    Real fCharge,
    deviceGauge* pToBeCached,
    BYTE byGaugeFieldId)
{
    intokernalE(8);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    const UBOOL bPlusX = (0 != (elementIdx & 1));
    const UBOOL bPlusY = (0 != (elementIdx & 2));
    const UBOOL bPlusT = (0 != (elementIdx & 4));

    //same check as _kernelDFermionKS_Cache_XYTau_TermT, do not access gauge links out of the Dirichlet region
    SSmallInt4 sOffset = sSite4;
    sOffset.x = sOffset.x + (bPlusX ? 1 : -1);
    sOffset.y = sOffset.y + (bPlusY ? 1 : -1);
    sOffset.w = sOffset.w + (bPlusT ? 1 : -1);
    if (__idx->m_pDeviceIndexPositionToSIndex[byGaugeFieldId][__bi(sOffset)].IsDirichlet())
    {
        //gauge field is identity across the boundary
        pToBeCached[uiSiteIndex * 8 + elementIdx] = _makeId<deviceGauge>();
        return;
    }

    pToBeCached[uiSiteIndex * 8 + elementIdx] = _deviceVXYTOptimizedEMT(pGauge, pPhase, sSite4, fCharge, byGaugeFieldId, bPlusX, bPlusY, bPlusT);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_Cache_XYTermEMT(
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pPhase,
    Real fCharge,
    deviceGauge* pToBeCached,
    BYTE byGaugeFieldId)
{
    intokernalE(8);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    const UBOOL bPlusMu = elementIdx & 2;
    const UBOOL bPlusTau = elementIdx & 4;
    const UINT bXorY = elementIdx & 1;

    //same check as _kernelDFermionKS_Cache_XYTermT, do not access gauge links out of the Dirichlet region
    const UINT bYorX = 1 - bXorY;
    SSmallInt4 sTargetSite = sSite4;
    sTargetSite.m_byData4[bYorX] = sTargetSite.m_byData4[bYorX] + (bPlusMu ? 2 : -2);
    sTargetSite.w = sTargetSite.w + (bPlusTau ? 1 : -1);
    if (__idx->m_pDeviceIndexPositionToSIndex[byGaugeFieldId][__bi(sTargetSite)].IsDirichlet())
    {
        //gauge field is identity across the boundary
        pToBeCached[uiSiteIndex * 8 + elementIdx] = _makeId<deviceGauge>();
        return;
    }

    pToBeCached[uiSiteIndex * 8 + elementIdx] = _deviceVXXTauOptimizedEMT(pGauge, pPhase, sSite4, fCharge, byGaugeFieldId, bXorY, bPlusMu, bPlusTau);
}


template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CacheKSRotationGaugeBuffer(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* res, const Real* pPhase, Real fCharge, UBOOL bHasCharge)
{
    preparethreadE(8);
    if (bHasCharge)
    {
        _LAUNCH_KERNEL(_kernelDFermionKS_Cache_XYTermEMT<deviceGauge>, block, threads,
            deviceData,
            pPhase,
            fCharge,
            res,
            byFieldId
        );

        _LAUNCH_KERNEL(_kernelDFermionKS_Cache_XYTau_TermEMT<deviceGauge>, block, threads,
            deviceData,
            pPhase,
            fCharge,
            res + _HC_Volume * 8,
            byFieldId
        );
    }
    else
    {
        _LAUNCH_KERNEL(_kernelDFermionKS_Cache_XYTermT<deviceGauge>, block, threads,
            deviceData,
            res,
            byFieldId
        );

        _LAUNCH_KERNEL(_kernelDFermionKS_Cache_XYTau_TermT<deviceGauge>, block, threads,
            deviceData,
            res + _HC_Volume * 8,
            byFieldId
        );
    }

}

#pragma endregion


#pragma region tree level - related

#pragma region kernel

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_Rectangular_NormalE(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const DOUBLE betaOverNTimesCRect,
    DOUBLE* results
)
{
    intokernalE(6);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    const SCHAR fwdir = _plaq_idx[elementIdx][0];
    const SCHAR fwnu = _plaq_idx[elementIdx][1];

    // | <-- <-- ^
    // v         |
    // o --> --> |
    SCHAR path[6] = { fwdir, fwdir, fwnu, static_cast<SCHAR>(-fwdir), static_cast<SCHAR>(-fwdir), static_cast<SCHAR>(-fwnu) };
    DOUBLE res = _retr(_deviceLinkT(pDeviceData, sSite4, 6, byFieldId, path));

    // |<--  ^
    // v     |
    // |     ^
    // v     |
    // o --> |
    path[0] = fwdir;
    path[1] = fwnu;
    path[2] = fwnu;
    path[3] = -fwdir;
    path[4] = -fwnu;
    path[5] = -fwnu;
    res += _retr(_deviceLinkT(pDeviceData, sSite4, 6, byFieldId, path));

    res = (6.0 - res) * betaOverNTimesCRect;
    if (0 == elementIdx)
    {
        results[uiSiteIndex] = res;
    }
    __syncthreads();

    #pragma unroll
    for (BYTE i = 1U; i < 6U; ++i)
    {
        if (i == elementIdx)
        {
            results[uiSiteIndex] += res;
        }
        __syncthreads();
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_Rectangular_Clover(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const DOUBLE betaOverNTimesCRect,
    DOUBLE* results
)
{
    intokernalDirInt4;

    DOUBLE res = 0.0;
    for (BYTE byNu = dir + 1; byNu < _DC_Dir; ++byNu)
    {
        res += 12.0 - _deviceCloverRectangleRetrT(byFieldId, pDeviceData, sSite4, dir, byNu);
    }
    //in total, 4*6=24 rectangular are calculated, but we need 12 rectangulars
    res = res * betaOverNTimesCRect * 0.5;
    if (0 == dir)
    {
        results[uiSiteIndex] = res;
    }
    __syncthreads();
    if (1 == dir)
    {
        results[uiSiteIndex] += res;
    }
    __syncthreads();
    if (2 == dir)
    {
        results[uiSiteIndex] += res;
    }
    __syncthreads();
    if (3 == dir)
    {
        results[uiSiteIndex] += res;
    }
}


// { fwnu, fwdir, fwdir, bcnu, bcdir };
//rectangluar force idx
// ^ nu
// |
// --->mu
__device__ __constant__ constexpr BYTE _rectangluar_force_idx[6][5] = {

    // ^ -->--> |
    // |        |
    // o    <-- v
    {2, 0, 0, 3, 1},

    // o    <-- ^
    // |        |
    // v -->--> |
    {3, 0, 0, 2, 1},

    // ^ -->--> |
    // |        |
    // <--o     v
    {1, 2, 0, 0, 3},

    // <--o     ^
    // |        |
    // v -->--> |
    {1, 3, 0, 0, 2},

    // ^ --> |
    // |     |
    // ^     v
    // |     |
    // o     v
    {2, 2, 0, 3, 3},

    // o     ^
    // |     |
    // v     ^
    // |     |
    // v --> |
    {3, 3, 0, 2, 2},
};

__device__ __constant__ constexpr SCHAR _rectangluar_force_corner_shift[6][2] = {

    // ^ -->--> |
    // |        |
    // o    <-- v
    {0, 0},

    // o    <-- ^
    // |        |
    // v -->--> |
    {0, -1},

    // ^ -->--> |
    // |        |
    // <--o     v
    {-1, 0},

    // <--o     ^
    // |        |
    // v -->--> |
    {-1, -1},

    // ^ --> |
    // |     |
    // ^     v
    // |     |
    // o     v
    {0, 0},

    // o     ^
    // |     |
    // v     ^
    // |     |
    // v --> |
    {0, -2},
};

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_Rectangular_ForceE(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE fCoeff,
    deviceGauge* pForce)
{
    intokernalEDir(18U);
    //const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);
    //const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
    //{
    //    return;
    //}

    /*old
    const BYTE dirid = elementIdx % 12U; //0-11
    const BYTE termid = elementIdx / 12U; //0-5
    //const BYTE nu = dirid >> 2U; //0-2
    const BYTE mu = dirid & 3U;  //0-3
    */
    const BYTE dirid = elementIdx % 3U; //0-2
    const BYTE termid = elementIdx / 3U; //0-5

    //mu = 0,1,2,3
    //nu != mu, so, 
    //0,1,2 -> 
    //0: 1,2,3
    //1: 2,3,0
    //2: 3,0,1
    //3: 0,1,2
    //which is 
    // (mu + nu + 1U) & 3U
    const BYTE nu = (dir + dirid + 1U) & 3U;
    const SCHAR dirs[4] = {
        static_cast<SCHAR>(dir + 1),
        static_cast<SCHAR>(-static_cast<SCHAR>(dir) - 1),
        static_cast<SCHAR>(nu + 1),
        static_cast<SCHAR>(-static_cast<SCHAR>(nu) - 1) };
    const SCHAR path[5] = {
        dirs[_rectangluar_force_idx[termid][0]],
        dirs[_rectangluar_force_idx[termid][1]],
        dirs[_rectangluar_force_idx[termid][2]],
        dirs[_rectangluar_force_idx[termid][3]],
        dirs[_rectangluar_force_idx[termid][4]]
    };

    //deviceGauge res = _deviceLinkT(pDeviceData, sSite4, 5, byFieldId, path);
    deviceGauge res = _deviceLinkT(pDeviceData, __deviceSiteIndexToInt4(uiSiteIndex), 5, byFieldId, path);
    _mul(res, fCoeff);
    for (BYTE i = 0U; i < 18U; ++i)
    {
        if (i == elementIdx)
        {
            _add(pForce[uiLinkIndex], res);
        }
        __syncthreads();
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_Rectangular_ForceE_D(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE fCoeff,
    deviceGauge* pForce,
    SSmallInt4 sPeriod)
{
    intokernalEDir(18U);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);
    //const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
    //{
    //    return;
    //}


    /*old
    const BYTE dirid = elementIdx % 12U; //0-11
    const BYTE termid = elementIdx / 12U; //0-5
    //const BYTE nu = dirid >> 2U; //0-2
    const BYTE mu = dirid & 3U;  //0-3
    */
    const BYTE dirid = elementIdx % 3U; //0-2
    const BYTE termid = elementIdx / 3U; //0-5

    //mu = 0,1,2,3
    //nu != mu, so, 
    //0,1,2 -> 
    //0: 1,2,3
    //1: 2,3,0
    //2: 3,0,1
    //3: 0,1,2
    //which is 
    // (mu + nu + 1U) & 3U
    const BYTE nu = (dir + dirid + 1U) & 3U;
    SSmallInt4 sCorner = sSite4;
    sCorner.m_byData4[dir] += _rectangluar_force_corner_shift[termid][0];
    sCorner.m_byData4[nu] += _rectangluar_force_corner_shift[termid][1];
    if (sCorner.m_byData4[dir] < 0 && 0 == sPeriod.m_byData4[dir])
    {
        return;
    }
    if (sCorner.m_byData4[nu] < 0 && 0 == sPeriod.m_byData4[nu])
    {
        return;
    }

    const SCHAR dirs[4] = {
        static_cast<SCHAR>(dir + 1),
        static_cast<SCHAR>(-static_cast<SCHAR>(dir) - 1),
        static_cast<SCHAR>(nu + 1),
        static_cast<SCHAR>(-static_cast<SCHAR>(nu) - 1) };

    const SCHAR path[5] = {
        dirs[_rectangluar_force_idx[termid][0]],
        dirs[_rectangluar_force_idx[termid][1]],
        dirs[_rectangluar_force_idx[termid][2]],
        dirs[_rectangluar_force_idx[termid][3]],
        dirs[_rectangluar_force_idx[termid][4]]
    };

    //deviceGauge res = _deviceLinkT(pDeviceData, sSite4, 5, byFieldId, path);
    deviceGauge res = _deviceLinkT(pDeviceData, sSite4, 5, byFieldId, path);
    _mul(res, fCoeff);
    for (BYTE i = 0U; i < 18U; ++i)
    {
        if (i == elementIdx)
        {
            _add(pForce[uiLinkIndex], res);
        }
        __syncthreads();
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_Rectangular_ForceE_CloverD(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE fCoeff,
    deviceGauge* pForce,
    SSmallInt4 sPeriod)
{
    intokernalEDir(18U);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);
    //const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
    //{
    //    return;
    //}


    /*old
    const BYTE dirid = elementIdx % 12U; //0-11
    const BYTE termid = elementIdx / 12U; //0-5
    //const BYTE nu = dirid >> 2U; //0-2
    const BYTE mu = dirid & 3U;  //0-3
    */
    const BYTE dirid = elementIdx % 3U; //0-2
    const BYTE termid = elementIdx / 3U; //0-5

    //mu = 0,1,2,3
    //nu != mu, so, 
    //0,1,2 -> 
    //0: 1,2,3
    //1: 2,3,0
    //2: 3,0,1
    //3: 0,1,2
    //which is 
    // (mu + nu + 1U) & 3U
    const BYTE nu = (dir + dirid + 1U) & 3U;
    if (4U == termid && sSite4.m_byData4[nu] == (_constIntegers[ECI_Lx + nu] - 1) && 0 == sPeriod.m_byData4[nu])
    {
        return;
    }

    if ((0U == termid || 1U == termid) && sSite4.m_byData4[dir] == (_constIntegers[ECI_Lx + dir] - 1) && 0 == sPeriod.m_byData4[dir])
    {
        return;
    }

    const SCHAR dirs[4] = {
        static_cast<SCHAR>(dir + 1),
        static_cast<SCHAR>(-static_cast<SCHAR>(dir) - 1),
        static_cast<SCHAR>(nu + 1),
        static_cast<SCHAR>(-static_cast<SCHAR>(nu) - 1) };

    const SCHAR path[5] = {
        dirs[_rectangluar_force_idx[termid][0]],
        dirs[_rectangluar_force_idx[termid][1]],
        dirs[_rectangluar_force_idx[termid][2]],
        dirs[_rectangluar_force_idx[termid][3]],
        dirs[_rectangluar_force_idx[termid][4]]
    };

    //deviceGauge res = _deviceLinkT(pDeviceData, sSite4, 5, byFieldId, path);
    deviceGauge res = _deviceLinkT(pDeviceData, sSite4, 5, byFieldId, path);
    _mul(res, fCoeff);
    if ((4U == termid || 5U == termid) && sSite4.m_byData4[dir] == (_constIntegers[ECI_Lx + dir] - 1) && 0 == sPeriod.m_byData4[dir])
    {
        _mul(res, F(0.5));
    }
    if ((0U == termid || 2U == termid) && sSite4.m_byData4[nu] == (_constIntegers[ECI_Lx + nu] - 1) && 0 == sPeriod.m_byData4[nu])
    {
        _mul(res, F(0.5));
    }
    for (BYTE i = 0U; i < 18U; ++i)
    {
        if (i == elementIdx)
        {
            _add(pForce[uiLinkIndex], res);
        }
        __syncthreads();
    }
}

#pragma endregion

#pragma region rectangular anisotropy

//the anisotropic counterparts of the rectangular kernels
//the plane weight convention is the same as the plaquette anisotropy kernels:
//spatial planes are multiplied by xi, planes containing the temporal direction by 1/xi

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_Rectangular_NormalEAnisotropy(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const DOUBLE betaOverNTimesCRect,
    const DOUBLE xi,
    DOUBLE* results
)
{
    intokernalE(6);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    const SCHAR fwdir = _plaq_idx[elementIdx][0];
    const SCHAR fwnu = _plaq_idx[elementIdx][1];

    // | <-- <-- ^
    // v         |
    // o --> --> |
    SCHAR path[6] = { fwdir, fwdir, fwnu, static_cast<SCHAR>(-fwdir), static_cast<SCHAR>(-fwdir), static_cast<SCHAR>(-fwnu) };
    DOUBLE res = _retr(_deviceLinkT(pDeviceData, sSite4, 6, byFieldId, path));

    // |<--  ^
    // v     |
    // |     ^
    // v     |
    // o --> |
    path[0] = fwdir;
    path[1] = fwnu;
    path[2] = fwnu;
    path[3] = -fwdir;
    path[4] = -fwnu;
    path[5] = -fwnu;
    res += _retr(_deviceLinkT(pDeviceData, sSite4, 6, byFieldId, path));

    //both loops are in the same plane, 4 is the temporal direction
    const DOUBLE fWeight = (4 == fwdir || 4 == fwnu) ? (1.0 / xi) : xi;
    res = (6.0 - res) * betaOverNTimesCRect * fWeight;
    if (0 == elementIdx)
    {
        results[uiSiteIndex] = res;
    }
    __syncthreads();

    #pragma unroll
    for (BYTE i = 1U; i < 6U; ++i)
    {
        if (i == elementIdx)
        {
            results[uiSiteIndex] += res;
        }
        __syncthreads();
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_Rectangular_CloverAnisotropy(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const DOUBLE betaOverNTimesCRect,
    const DOUBLE xi,
    DOUBLE* results
)
{
    intokernalDirInt4;

    DOUBLE res = 0.0;
    for (BYTE byNu = dir + 1; byNu < _DC_Dir; ++byNu)
    {
        //dir and byNu are 0 based, 3 is the temporal direction
        const DOUBLE fWeight = (3 == dir || 3 == byNu) ? (1.0 / xi) : xi;
        res += (12.0 - _deviceCloverRectangleRetrT(byFieldId, pDeviceData, sSite4, dir, byNu)) * fWeight;
    }
    //in total, 4*6=24 rectangular are calculated, but we need 12 rectangulars
    res = res * betaOverNTimesCRect * 0.5;
    if (0 == dir)
    {
        results[uiSiteIndex] = res;
    }
    __syncthreads();
    if (1 == dir)
    {
        results[uiSiteIndex] += res;
    }
    __syncthreads();
    if (2 == dir)
    {
        results[uiSiteIndex] += res;
    }
    __syncthreads();
    if (3 == dir)
    {
        results[uiSiteIndex] += res;
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_Rectangular_ForceEAnisotropy(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE fCoeff,
    const DOUBLE xi,
    deviceGauge* pForce)
{
    intokernalEDir(18U);

    const BYTE dirid = elementIdx % 3U; //0-2
    const BYTE termid = elementIdx / 3U; //0-5

    //nu != dir, (dir + dirid + 1U) & 3U
    const BYTE nu = (dir + dirid + 1U) & 3U;
    const SCHAR dirs[4] = {
        static_cast<SCHAR>(dir + 1),
        static_cast<SCHAR>(-static_cast<SCHAR>(dir) - 1),
        static_cast<SCHAR>(nu + 1),
        static_cast<SCHAR>(-static_cast<SCHAR>(nu) - 1) };
    const SCHAR path[5] = {
        dirs[_rectangluar_force_idx[termid][0]],
        dirs[_rectangluar_force_idx[termid][1]],
        dirs[_rectangluar_force_idx[termid][2]],
        dirs[_rectangluar_force_idx[termid][3]],
        dirs[_rectangluar_force_idx[termid][4]]
    };

    deviceGauge res = _deviceLinkT(pDeviceData, __deviceSiteIndexToInt4(uiSiteIndex), 5, byFieldId, path);
    //dir and nu are 0 based, 3 is the temporal direction
    _mul(res, fCoeff * ((3 == dir || 3 == nu) ? (1.0 / xi) : xi));
    for (BYTE i = 0U; i < 18U; ++i)
    {
        if (i == elementIdx)
        {
            _add(pForce[uiLinkIndex], res);
        }
        __syncthreads();
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_Rectangular_ForceEAnisotropy_D(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE fCoeff,
    const DOUBLE xi,
    deviceGauge* pForce,
    SSmallInt4 sPeriod)
{
    intokernalEDir(18U);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    const BYTE dirid = elementIdx % 3U; //0-2
    const BYTE termid = elementIdx / 3U; //0-5

    //nu != dir, (dir + dirid + 1U) & 3U
    const BYTE nu = (dir + dirid + 1U) & 3U;
    SSmallInt4 sCorner = sSite4;
    sCorner.m_byData4[dir] += _rectangluar_force_corner_shift[termid][0];
    sCorner.m_byData4[nu] += _rectangluar_force_corner_shift[termid][1];
    if (sCorner.m_byData4[dir] < 0 && 0 == sPeriod.m_byData4[dir])
    {
        return;
    }
    if (sCorner.m_byData4[nu] < 0 && 0 == sPeriod.m_byData4[nu])
    {
        return;
    }

    const SCHAR dirs[4] = {
        static_cast<SCHAR>(dir + 1),
        static_cast<SCHAR>(-static_cast<SCHAR>(dir) - 1),
        static_cast<SCHAR>(nu + 1),
        static_cast<SCHAR>(-static_cast<SCHAR>(nu) - 1) };

    const SCHAR path[5] = {
        dirs[_rectangluar_force_idx[termid][0]],
        dirs[_rectangluar_force_idx[termid][1]],
        dirs[_rectangluar_force_idx[termid][2]],
        dirs[_rectangluar_force_idx[termid][3]],
        dirs[_rectangluar_force_idx[termid][4]]
    };

    deviceGauge res = _deviceLinkT(pDeviceData, sSite4, 5, byFieldId, path);
    //dir and nu are 0 based, 3 is the temporal direction
    _mul(res, fCoeff * ((3 == dir || 3 == nu) ? (1.0 / xi) : xi));
    for (BYTE i = 0U; i < 18U; ++i)
    {
        if (i == elementIdx)
        {
            _add(pForce[uiLinkIndex], res);
        }
        __syncthreads();
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_Rectangular_ForceE_CloverAnisotropy_D(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE fCoeff,
    const DOUBLE xi,
    deviceGauge* pForce,
    SSmallInt4 sPeriod)
{
    intokernalEDir(18U);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    const BYTE dirid = elementIdx % 3U; //0-2
    const BYTE termid = elementIdx / 3U; //0-5

    //nu != dir, (dir + dirid + 1U) & 3U
    const BYTE nu = (dir + dirid + 1U) & 3U;
    if (4U == termid && sSite4.m_byData4[nu] == (_constIntegers[ECI_Lx + nu] - 1) && 0 == sPeriod.m_byData4[nu])
    {
        return;
    }

    if ((0U == termid || 1U == termid) && sSite4.m_byData4[dir] == (_constIntegers[ECI_Lx + dir] - 1) && 0 == sPeriod.m_byData4[dir])
    {
        return;
    }

    const SCHAR dirs[4] = {
        static_cast<SCHAR>(dir + 1),
        static_cast<SCHAR>(-static_cast<SCHAR>(dir) - 1),
        static_cast<SCHAR>(nu + 1),
        static_cast<SCHAR>(-static_cast<SCHAR>(nu) - 1) };

    const SCHAR path[5] = {
        dirs[_rectangluar_force_idx[termid][0]],
        dirs[_rectangluar_force_idx[termid][1]],
        dirs[_rectangluar_force_idx[termid][2]],
        dirs[_rectangluar_force_idx[termid][3]],
        dirs[_rectangluar_force_idx[termid][4]]
    };

    deviceGauge res = _deviceLinkT(pDeviceData, sSite4, 5, byFieldId, path);
    //dir and nu are 0 based, 3 is the temporal direction
    _mul(res, fCoeff * ((3 == dir || 3 == nu) ? (1.0 / xi) : xi));
    if ((4U == termid || 5U == termid) && sSite4.m_byData4[dir] == (_constIntegers[ECI_Lx + dir] - 1) && 0 == sPeriod.m_byData4[dir])
    {
        _mul(res, F(0.5));
    }
    if ((0U == termid || 2U == termid) && sSite4.m_byData4[nu] == (_constIntegers[ECI_Lx + nu] - 1) && 0 == sPeriod.m_byData4[nu])
    {
        _mul(res, F(0.5));
    }
    for (BYTE i = 0U; i < 18U; ++i)
    {
        if (i == elementIdx)
        {
            _add(pForce[uiLinkIndex], res);
        }
        __syncthreads();
    }
}

#pragma endregion

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculateRectangularEnergy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverNtimeRect)
{
    preparethreadE(6);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_Rectangular_NormalE<deviceGauge>, block, threads, byFieldId, deviceData, betaOverNtimeRect, _D_RealThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculateRectangularEnergyUseClover(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverNtimeRect)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_Rectangular_Clover<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        betaOverNtimeRect,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceRectangular(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeRect)
{
    preparethreadEDir(18);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_Rectangular_ForceE<deviceGauge>, block, threads, byFieldId, deviceData, betaOverNtimeRect, pForce);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceRectangular_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeRect, SSmallInt4 sPeriod)
{
    preparethreadEDir(18);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_Rectangular_ForceE_D<deviceGauge>, block, threads, byFieldId, deviceData, betaOverNtimeRect, pForce, sPeriod);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceRectangularClover_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeRect, SSmallInt4 sPeriod)
{
    preparethreadEDir(18);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_Rectangular_ForceE_CloverD<deviceGauge>, block, threads, byFieldId, deviceData, betaOverNtimeRect, pForce, sPeriod);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculateRectangularEnergyAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverNtimeRect, DOUBLE xi)
{
    preparethreadE(6);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_Rectangular_NormalEAnisotropy<deviceGauge>, block, threads, byFieldId, deviceData, betaOverNtimeRect, xi, _D_RealThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculateRectangularEnergyUseCloverAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverNtimeRect, DOUBLE xi)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_Rectangular_CloverAnisotropy<deviceGauge>, block, threads,
        byFieldId,
        deviceData,
        betaOverNtimeRect,
        xi,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceRectangularAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeRect, DOUBLE xi)
{
    preparethreadEDir(18);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_Rectangular_ForceEAnisotropy<deviceGauge>, block, threads, byFieldId, deviceData, betaOverNtimeRect, xi, pForce);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceRectangularAnisotropy_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeRect, DOUBLE xi, SSmallInt4 sPeriod)
{
    preparethreadEDir(18);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_Rectangular_ForceEAnisotropy_D<deviceGauge>, block, threads, byFieldId, deviceData, betaOverNtimeRect, xi, pForce, sPeriod);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceRectangularCloverAnisotropy_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeRect, DOUBLE xi, SSmallInt4 sPeriod)
{
    preparethreadEDir(18);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_Rectangular_ForceE_CloverAnisotropy_D<deviceGauge>, block, threads, byFieldId, deviceData, betaOverNtimeRect, xi, pForce, sPeriod);
}

/**
*  |
* ||
*  |
* is like
*
* ____
* |  |
* |__|
*
* so just calculate all "staples" with one of the links replaced by f0
* o is start, x is end
*
*   /|
*  / |
* |  |
* |  |
* x  |
*    |
*    o
*
*
*   /|
*  / |
* x  |
*    |
*    |
* o  |
*  \ |
*   \|
*
*
*    x
*    |
* o  |
* |  |
* |  |
*  \ |
*   \|
*
*
* We loop for every "long" link, and add force for its three contributions
*
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelNaikForce_Old(
    const deviceGauge* __restrict__ f0,
    const deviceGauge* __restrict__ gauge,
    deviceGauge* force1,
    deviceGauge* force2,
    deviceGauge* force3,
    BYTE byGaugeFieldId)
{
    intokernalInt4;

    for (BYTE dir = 0; dir < _DC_Dir; ++dir)
    {
        const BYTE fwddir = dir + 1;
        const SSmallInt4 n2 = _deviceSmallInt4OffsetC(sSite4, fwddir);
        const SSmallInt4 n3 = _deviceSmallInt4OffsetC(n2, fwddir);

        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        //const SIndex& n1__mu = __idx->m_pDeviceIndexLinkToSIndex[byGaugeFieldId][__bi4(sSite4) + dir];
        const SIndex n2__mu = __idx->m_pDeviceIndexLinkToSIndex[byGaugeFieldId][__bi4(n2) + dir];
        const SIndex n3__mu = __idx->m_pDeviceIndexLinkToSIndex[byGaugeFieldId][__bi4(n3) + dir];
        const UINT linkIndex2 = _deviceGetLinkIndex(n2__mu.m_uiSiteIndex, dir);
        const UINT linkIndex3 = _deviceGetLinkIndex(n3__mu.m_uiSiteIndex, dir);

        deviceGauge f1 = _muldagC(f0[uiLinkIndex], gauge[linkIndex3]);
        _muldag(f1, gauge[linkIndex2]);

        deviceGauge f2 = _dagmulC(gauge[uiLinkIndex], f0[uiLinkIndex]);
        _muldag(f2, gauge[linkIndex3]);

        deviceGauge f3 = _mulC(gauge[uiLinkIndex], gauge[linkIndex2]);
        _dagmul(f3, f0[uiLinkIndex]);

        force1[uiLinkIndex] = f1;
        force2[linkIndex2] = f2;
        force3[linkIndex3] = f3;
    }

}

/**
* very IMPORTANT NOTE here:
* we usually calculate force on a SINGLE link at a time
* however, it seems cheaper to calculate force for different links here
* so we modify the values of three links here
*
* This will cause problem when crossing different blocks
* Different blocks may modify a same link at the same time!
* this cannot be solved even with __syncthreads()
* we have to use three different buffers!
*
* Note2:
* If do not use even-odd, it is just -(f)^+
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelNaikForce(
    const deviceGauge* __restrict__ f0,
    const deviceGauge* __restrict__ gauge,
    deviceGauge* force1,
    deviceGauge* force2,
    deviceGauge* force3,
    BYTE byGaugeFieldId)
{
    intokernalDirInt4;

    const BYTE fwddir = dir + 1;
    const SSmallInt4 p1 = _deviceSmallInt4OffsetC(sSite4, static_cast<SCHAR>(fwddir));
    const SSmallInt4 p2 = _deviceSmallInt4OffsetC(p1, static_cast<SCHAR>(fwddir));
    const SSmallInt4 m1 = _deviceSmallInt4OffsetC(sSite4, static_cast<SCHAR>(-static_cast<INT>(fwddir)));
    const SSmallInt4 m2 = _deviceSmallInt4OffsetC(m1, static_cast<SCHAR>(-static_cast<INT>(fwddir)));

    const SIndex p1__mu = __idx->m_pDeviceIndexLinkToSIndex[byGaugeFieldId][__bi4(p1) + dir];
    const SIndex p2__mu = __idx->m_pDeviceIndexLinkToSIndex[byGaugeFieldId][__bi4(p2) + dir];
    const SIndex m1__mu = __idx->m_pDeviceIndexLinkToSIndex[byGaugeFieldId][__bi4(m1) + dir];
    const SIndex m2__mu = __idx->m_pDeviceIndexLinkToSIndex[byGaugeFieldId][__bi4(m2) + dir];
    const UINT linkIndexP1 = _deviceGetLinkIndex(p1__mu.m_uiSiteIndex, dir);
    const UINT linkIndexP2 = _deviceGetLinkIndex(p2__mu.m_uiSiteIndex, dir);
    const UINT linkIndexM1 = _deviceGetLinkIndex(m1__mu.m_uiSiteIndex, dir);
    const UINT linkIndexM2 = _deviceGetLinkIndex(m2__mu.m_uiSiteIndex, dir);

    //Multi-GPU: this kernel is written in GATHER form (every local link computes
    //its own f1/f2/f3). The original SCATTER form wrote f2 at site+1 and f3 at
    //site+2, which under domain decomposition lands in the halo tail for
    //cross-rank neighbours and leaves the owner rank's boundary links unwritten.
    //Each link's expression is the SAME sequence of matrix ops as the scatter
    //version (evaluated at the same source site), so single-GPU results are
    //bit-identical. Gather reads reach up to 2 hops, so the caller must refill
    //the f0 and gauge halos (width >= 2) before launch.
    //f1[s] = f0[s] g[s+2]^dag g[s+1]^dag (was already local)
    deviceGauge f1 = _muldagC(f0[uiLinkIndex], gauge[linkIndexP2]);
    _muldag(f1, gauge[linkIndexP1]);

    //f2[L] = g[L-1]^dag f0[L-1] g[L+1]^dag (scatter wrote this from s=L-1 into L)
    deviceGauge f2 = _dagmulC(gauge[linkIndexM1], f0[linkIndexM1]);
    _muldag(f2, gauge[linkIndexP1]);

    //f3[L] = (g[L-2] g[L-1])^dag f0[L-2] (scatter wrote this from s=L-2 into L)
    deviceGauge f3 = _mulC(gauge[linkIndexM2], gauge[linkIndexM1]);
    _dagmul(f3, f0[linkIndexM2]);

    force1[uiLinkIndex] = f1;

    force2[uiLinkIndex] = f2;

    force3[uiLinkIndex] = f3;
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::NaikForce(const CFieldGauge* naikf0, CFieldGauge* naikforce)
{
    _RECORD(CFieldGaugeKernel::NaikForce);
    appParanoiac(_T("CFieldFermionHISQSU3 NaikForce\n"));
    CFieldGauge* pNaikForce2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(naikforce->m_byFieldId, _T(__FILE__), __LINE__));
    CFieldGauge* pNaikForce3 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(naikforce->m_byFieldId, _T(__FILE__), __LINE__));
    const deviceGauge* pEffectiveLevel1 = (const deviceGauge*)(appGetGaugeSmearing(naikf0->m_byFieldId)->GetEffectiveGaugeLevel1()->GetData());
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelNaikForce<deviceGauge>, block, threads,
        (const deviceGauge*)(naikf0->GetData()),
        pEffectiveLevel1,
        (deviceGauge*)(naikforce->GetData()),
        (deviceGauge*)(pNaikForce2->GetData()),
        (deviceGauge*)(pNaikForce3->GetData()),
        naikf0->m_byFieldId
    );
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    naikforce->AxpyPlus(pNaikForce2);
    naikforce->AxpyPlus(pNaikForce3);
    pNaikForce2->Return();
    pNaikForce3->Return();
}

#pragma endregion

#pragma region one-loop gauge

#pragma region kernel

/**
* copy from pyQuda
* 
            [X, Y, Z, -X, -Y, -Z],
            [X, Y, -Z, -X, -Y, Z],
            [X, -Y, Z, -X, Y, -Z],
            [X, -Y, -Z, -X, Y, Z],
            [X, Y, T, -X, -Y, -T],
            [X, Y, -T, -X, -Y, T],
            [X, -Y, T, -X, Y, -T],
            [X, -Y, -T, -X, Y, T],
            [X, Z, T, -X, -Z, -T],
            [X, Z, -T, -X, -Z, T],
            [X, -Z, T, -X, Z, -T],
            [X, -Z, -T, -X, Z, T],
            [Y, Z, T, -Y, -Z, -T],
            [Y, Z, -T, -Y, -Z, T],
            [Y, -Z, T, -Y, Z, -T],
            [Y, -Z, -T, -Y, Z, T],
*
* so, it is xyz, xyt, xzt, yzt
* each has four directions (E=16)
* 0 1 2
* 0 1 3
* 0 2 3
* 1 2 3
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_OneLoop_NormalE(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const DOUBLE betaOverNTimesCRect,
    DOUBLE* results
)
{
    intokernalE(16U);
    const BYTE dirs  = elementIdx >> 2U;
    const SCHAR signs = static_cast<SCHAR>(elementIdx & 3U);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    const SCHAR dir1 = static_cast<SCHAR>(1U + (dirs / 3U));
    SCHAR dir2 = static_cast<SCHAR>(2U + (dirs >> 1U));
    SCHAR dir3 = static_cast<SCHAR>(3U + ((dirs + 2U) / 3U));

    //printf("%d %d %d %d, signs=%d sl=%d, sr=%d\n", dirs, dir1, dir2, dir3, signs, 1 - (signs&2), 1 - ((signs & 1) << 1));

    dir2 = dir2 * (1 - (signs&2));
    dir3 = dir3 * (1 - ((signs&1) << 1));

    
    SCHAR path[6] = { dir1, dir2, dir3, static_cast<SCHAR>(-dir1), static_cast<SCHAR>(-dir2), static_cast<SCHAR>(-dir3) };
    DOUBLE res = _retr(_deviceLinkT(pDeviceData, sSite4, 6, byFieldId, path));

    res = (3.0 - res) * betaOverNTimesCRect;
    if (0 == elementIdx)
    {
        results[uiSiteIndex] = res;
    }
    __syncthreads();

    #pragma unroll
    for (BYTE i = 1U; i < 16U; ++i)
    {
        if (i == elementIdx)
        {
            results[uiSiteIndex] += res;
        }
        __syncthreads();
    }
}

/**
* copy from pyQuda
*
            [X, Y, Z, -X, -Y, -Z],
            [X, Y, -Z, -X, -Y, Z],
            [X, -Y, Z, -X, Y, -Z],
            [X, -Y, -Z, -X, Y, Z],
            [X, Y, T, -X, -Y, -T],
            [X, Y, -T, -X, -Y, T],
            [X, -Y, T, -X, Y, -T],
            [X, -Y, -T, -X, Y, T],
            [X, Z, T, -X, -Z, -T],
            [X, Z, -T, -X, -Z, T],
            [X, -Z, T, -X, Z, -T],
            [X, -Z, -T, -X, Z, T],
            [Y, Z, T, -Y, -Z, -T],
            [Y, Z, -T, -Y, -Z, T],
            [Y, -Z, T, -Y, Z, -T],
            [Y, -Z, -T, -Y, Z, T],
*
* so, taking X for example, the left-right directions are:
* +-Z X +-Y
* +-T X +-Y
* +-T X +-Z
* 
* taking Z for example,
* 
* +-Y Z +-X
* +-X Z +-T
* +-Y Z +-T
* 
* so, if arrange as [a,b,c,d,e]
* it was always bcd, acd, bce
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_OneLoop_ForceE(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE fCoeff,
    deviceGauge* pForce)
{
    intokernalEDir(12U);
    const BYTE dirs = elementIdx >> 2U;
    const SCHAR signs = static_cast<SCHAR>(elementIdx & 3U);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    constexpr SCHAR dirtableperodic[8] = {4, 3, 2, 1, 4, 3, 2, 1};

    //x->3, y->2, z->5, t->4
    const BYTE middle = (3U - dir) + ((dir >> 1U) << 2U);
    const BYTE left = middle - 1U - (dirs >> 1U);
    const BYTE right = middle + 1U + (dirs & 1U);
    
    const SCHAR middir = dirtableperodic[middle];
    const SCHAR leftdir = dirtableperodic[left] * (1 - (signs&2));
    const SCHAR rightdir = dirtableperodic[right] * (1 - ((signs&1) << 1));

    SCHAR path[5] = { rightdir, leftdir, middir, static_cast<SCHAR>(-rightdir), static_cast<SCHAR>(-leftdir)};

    deviceGauge res = _deviceLinkT(pDeviceData, __deviceSiteIndexToInt4(uiSiteIndex), 5, byFieldId, path);

    path[0] = leftdir;
    path[1] = rightdir;
    path[2] = middir;
    path[3] = static_cast<SCHAR>(-leftdir);
    path[4] = static_cast<SCHAR>(-rightdir);
    _add(res, _deviceLinkT(pDeviceData, __deviceSiteIndexToInt4(uiSiteIndex), 5, byFieldId, path));

    _mul(res, fCoeff);
    for (BYTE i = 0U; i < 12U; ++i)
    {
        if (i == elementIdx)
        {
            _add(pForce[uiLinkIndex], res);
        }
        __syncthreads();
    }
}

#pragma endregion

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculateTwistedLoopEnergy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverNtimeCT)
{
    preparethreadE(16);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_OneLoop_NormalE<deviceGauge>, block, threads, byFieldId, deviceData, betaOverNtimeCT, _D_RealThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceTwistedLoop(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeCT)
{
    preparethreadEDir(12);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_OneLoop_ForceE<deviceGauge>, block, threads, byFieldId, deviceData, betaOverNtimeCT, pForce);
}

//the anisotropic counterparts of the twisted loop kernels
//every loop spans a triple of directions, triples containing the temporal
//direction (4) are multiplied by 1/xi, spatial triples by xi

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_OneLoop_NormalEAnisotropy(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const DOUBLE betaOverNTimesCRect,
    const DOUBLE xi,
    DOUBLE* results
)
{
    intokernalE(16U);
    const BYTE dirs  = elementIdx >> 2U;
    const SCHAR signs = static_cast<SCHAR>(elementIdx & 3U);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    const SCHAR dir1 = static_cast<SCHAR>(1U + (dirs / 3U));
    SCHAR dir2 = static_cast<SCHAR>(2U + (dirs >> 1U));
    SCHAR dir3 = static_cast<SCHAR>(3U + ((dirs + 2U) / 3U));

    dir2 = dir2 * (1 - (signs&2));
    dir3 = dir3 * (1 - ((signs & 1) << 1));

    SCHAR path[6] = { dir1, dir2, dir3, static_cast<SCHAR>(-dir1), static_cast<SCHAR>(-dir2), static_cast<SCHAR>(-dir3) };
    DOUBLE res = _retr(_deviceLinkT(pDeviceData, sSite4, 6, byFieldId, path));

    //dirs == 0 is the purely spatial triple (1,2,3), the other triples contain the temporal direction.
    //In the weak field expansion a twisted loop [mu,nu,rho,-mu,-nu,-rho] contains all three
    //planes F_mn, F_mr, F_nr, so the purely spatial triple carries (2 xi - 1/xi) while
    //time containing triples carry 1/xi (continuum-matched, see Docs anisotropy notes).
    const DOUBLE fWeight = (0 == dirs) ? (2.0 * xi - 1.0 / xi) : (1.0 / xi);
    res = (3.0 - res) * betaOverNTimesCRect * fWeight;
    if (0 == elementIdx)
    {
        results[uiSiteIndex] = res;
    }
    __syncthreads();

    #pragma unroll
    for (BYTE i = 1U; i < 16U; ++i)
    {
        if (i == elementIdx)
        {
            results[uiSiteIndex] += res;
        }
        __syncthreads();
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergy_OneLoop_ForceEAnisotropy(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE fCoeff,
    const DOUBLE xi,
    deviceGauge* pForce)
{
    intokernalEDir(12U);
    const BYTE dirs = elementIdx >> 2U;
    const SCHAR signs = static_cast<SCHAR>(elementIdx & 3U);
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    constexpr SCHAR dirtableperodic[8] = {4, 3, 2, 1, 4, 3, 2, 1};

    //x->3, y->2, z->5, t->4
    const BYTE middle = (3U - dir) + ((dir >> 1U) << 2U);
    const BYTE left = middle - 1U - (dirs >> 1U);
    const BYTE right = middle + 1U + (dirs & 1U);

    const SCHAR middir = dirtableperodic[middle];
    const SCHAR leftdir = dirtableperodic[left] * (1 - (signs&2));
    const SCHAR rightdir = dirtableperodic[right] * (1 - ((signs&1) << 1));

    SCHAR path[5] = { rightdir, leftdir, middir, static_cast<SCHAR>(-rightdir), static_cast<SCHAR>(-leftdir)};

    deviceGauge res = _deviceLinkT(pDeviceData, __deviceSiteIndexToInt4(uiSiteIndex), 5, byFieldId, path);

    path[0] = leftdir;
    path[1] = rightdir;
    path[2] = middir;
    path[3] = static_cast<SCHAR>(-leftdir);
    path[4] = static_cast<SCHAR>(-rightdir);
    _add(res, _deviceLinkT(pDeviceData, __deviceSiteIndexToInt4(uiSiteIndex), 5, byFieldId, path));

    //4 is the temporal direction
    //the weight must be the same function of the loop triple as in the energy kernel:
    //purely spatial triples carry (2 xi - 1/xi), time containing triples carry 1/xi
    const UBOOL bTemporal = (4 == middir) || (4 == abs(leftdir)) || (4 == abs(rightdir));
    _mul(res, fCoeff * (bTemporal ? (1.0 / xi) : (2.0 * xi - 1.0 / xi)));
    for (BYTE i = 0U; i < 12U; ++i)
    {
        if (i == elementIdx)
        {
            _add(pForce[uiLinkIndex], res);
        }
        __syncthreads();
    }
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculateTwistedLoopEnergyAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverNtimeCT, DOUBLE xi)
{
    preparethreadE(16);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_OneLoop_NormalEAnisotropy<deviceGauge>, block, threads, byFieldId, deviceData, betaOverNtimeCT, xi, _D_RealThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceTwistedLoopAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeCT, DOUBLE xi)
{
    preparethreadEDir(12);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergy_OneLoop_ForceEAnisotropy<deviceGauge>, block, threads, byFieldId, deviceData, betaOverNtimeCT, xi, pForce);
}

#pragma endregion

#pragma region ZN partial specialization

template<INT N>
void CFieldGaugeKernel<deviceZN<N>, 1>::CalculateOnlyStaple(const deviceZN<N>* deviceData, BYTE byFieldId, deviceZN<N>* pStaple)
{
    preparethreadDir;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelCalculateOnlyStapleGauge<deviceZN<N>>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStaple);
#else
    _LAUNCH_KERNEL(_kernelCalculateOnlyStapleGauge<deviceZN<N>>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        pStaple);
#endif
}

// Note: parameter name is betaOverN for interface consistency, but for discrete
// groups (Z_N, D_N) the caller passes the full beta (not beta/N).
template<INT N>
DOUBLE CFieldGaugeKernel<deviceZN<N>, 1>::CalculatePlaqutteEnergy(const deviceZN<N>* deviceData, BYTE byFieldId, DOUBLE betaOverN)
{
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId]);

#if !_CLG_ASSUME_SQUARE_LATTICE
    preparethreadE(appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyGaugeCacheIndex<deviceZN<N>>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer
    );
#else
    preparethreadE(plaqCountPerSite);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyGaugeCacheIndex<deviceZN<N>>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId],
        betaOverN,
        _D_RealThreadBuffer
    );
#endif

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}


#pragma endregion

template class CFieldGaugeKernel<CLGComplex, 1>;
template class CFieldGaugeKernel<deviceSU2, 2>;
template class CFieldGaugeKernel<deviceSU3, 3>;

#if _CLG_SU4_GAUGE
template class CFieldGaugeKernel<deviceSU4, 4>;
#endif
#if _CLG_SU5_GAUGE
template class CFieldGaugeKernel<deviceSU5, 5>;
#endif
#if _CLG_SU6_GAUGE
template class CFieldGaugeKernel<deviceSU6, 6>;
#endif
#if _CLG_SU7_GAUGE
template class CFieldGaugeKernel<deviceSU7, 7>;
#endif
#if _CLG_SU8_GAUGE
template class CFieldGaugeKernel<deviceSU8, 8>;
#endif

#if _CLG_Z2_GAUGE
template class CFieldGaugeKernel<deviceZN<2>, 1>;
#endif
#if _CLG_Z3_GAUGE
template class CFieldGaugeKernel<deviceZN<3>, 1>;
#endif
#if _CLG_Z4_GAUGE
template class CFieldGaugeKernel<deviceZN<4>, 1>;
#endif
#if _CLG_Z5_GAUGE
template class CFieldGaugeKernel<deviceZN<5>, 1>;
#endif
#if _CLG_Z6_GAUGE
template class CFieldGaugeKernel<deviceZN<6>, 1>;
#endif

#pragma region DN partial specialization

template<INT N>
void CFieldGaugeKernel<deviceDN<N>, 2>::CalculateOnlyStaple(const deviceDN<N>* deviceData, BYTE byFieldId, deviceDN<N>* pStaple)
{
    preparethreadDir;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelCalculateOnlyStapleGauge<deviceDN<N>>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStaple);
#else
    _LAUNCH_KERNEL(_kernelCalculateOnlyStapleGauge<deviceDN<N>>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byFieldId],
        pStaple);
#endif
}

// Note: parameter name is betaOverN for interface consistency, but for discrete
// groups (Z_N, D_N) the caller passes the full beta (not beta/N).
template<INT N>
DOUBLE CFieldGaugeKernel<deviceDN<N>, 2>::CalculatePlaqutteEnergy(const deviceDN<N>* deviceData, BYTE byFieldId, DOUBLE betaOverN)
{
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId]);

#if !_CLG_ASSUME_SQUARE_LATTICE
    preparethreadE(appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyGaugeCacheIndex<deviceDN<N>>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer
    );
#else
    preparethreadE(plaqCountPerSite);
    _LAUNCH_KERNEL(_kernelPlaqutteEnergyGaugeCacheIndex<deviceDN<N>>, block, threads,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[byFieldId],
        betaOverN,
        _D_RealThreadBuffer
    );
#endif

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

#pragma endregion

#if _CLG_D3_GAUGE
template class CFieldGaugeKernel<deviceDN<3>, 2>;
#endif
#if _CLG_D4_GAUGE
template class CFieldGaugeKernel<deviceDN<4>, 2>;
#endif
#if _CLG_D8_GAUGE
template class CFieldGaugeKernel<deviceDN<8>, 2>;
#endif

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================