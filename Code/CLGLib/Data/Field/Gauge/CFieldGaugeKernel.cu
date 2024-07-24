//=============================================================================
// FILENAME : CFieldGaugeKernel.cu
// 
// DESCRIPTION:
//
// REVISION:
//  [07/24/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldGaugeKernel.h"
#include "CFieldGaugeLink.h"

__BEGIN_NAMESPACE

#pragma region gauge field

#pragma region Kernels

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteGaugeCacheIndex(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceGauge* pStapleData, //can be NULL
    deviceGauge* pForceData,
    Real betaOverN)
{
    intokernaldir;

    //Real test_force = F(0.0);
    betaOverN = betaOverN * F(-0.5);
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceGauge res = _makeZero<deviceGauge>();

        //there are 6 staples, each is sum of two plaquttes
        for (INT i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceGauge toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

            if (first.NeedToDagger())
            {
                _dagger(toAdd);
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
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
            pStapleData[linkIndex] = res;
        }

        //staple calculated
        deviceGauge force = pDeviceData[linkIndex];
        _muldag(force, res);
        _ta(force);
        _mul(force, betaOverN);

        //force is additive
        _add(pForceData[linkIndex], force);
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateOnlyStapleGauge(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceGauge* pStapleData)
{
    intokernaldir;

    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceGauge res = _makeZero<deviceGauge>();

        //there are 6 staples, each is sum of two plaquttes
        for (INT i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceGauge toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

            if (first.NeedToDagger())
            {
                _dagger(toAdd);
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
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
        pStapleData[linkIndex] = res;
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyGaugeCacheIndex(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    DOUBLE betaOverN,
    DOUBLE* results
)
{
    intokernal;

    DOUBLE resThisThread = 0.0;
    const UINT indexSkip = plaqCount * plaqLength * uiSiteIndex;
    for (BYTE i = 0; i < plaqCount; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLength + indexSkip];
        deviceGauge toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }

        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedIndex[i * plaqLength + j + indexSkip];
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

        resThisThread += (_dim<deviceGauge>() - _retr(toAdd));
    }

    results[uiSiteIndex] = resThisThread * betaOverN;
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
    intokernalInt4;

    DOUBLE fRes = 0.0;
    for (BYTE byDir1 = 0; byDir1 < _DC_Dir; ++byDir1)
    {
        for (BYTE byDir2 = byDir1 + 1; byDir2 < _DC_Dir; ++byDir2)
        {
            fRes += _deviceCloverRetrT(pDeviceData, sSite4, __bi(sSite4), byDir1, byDir2, byFieldId);
        }
    }
    fRes = stapleConstant - 0.25 * fRes;
    results[uiSiteIndex] = fRes * fBetaOverN;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyUsingStableGauge(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pStableData,
    DOUBLE stapleConstant,
    DOUBLE betaOverN,
    DOUBLE* results
)
{
    intokernaldir;

    DOUBLE resThisThread = 0.0;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        resThisThread += (stapleConstant - _retr(_muldagC(pDeviceData[linkIndex], pStableData[linkIndex])));
    }

    results[uiSiteIndex] = resThisThread * betaOverN * 0.25;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteCacheIndexT_D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceGauge* pStapleData, //can be NULL
    deviceGauge* pForceData,
    Real betaOverN)
{
    intokernaldir;

    //Real test_force = F(0.0);
    betaOverN = betaOverN * F(-0.5);
    const UINT plaqLengthm1 = plaqLength - 1;
    const UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceGauge res = _makeZero<deviceGauge>();

        //there are 6 staples, each is sum of two plaquttes
        for (BYTE i = 0; i < plaqCount; ++i)
        {
            BYTE diricCount = 0;
            const SIndex& first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
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
                const SIndex& nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
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
            pStapleData[linkIndex] = res;
        }

        //staple calculated
        deviceGauge force(pDeviceData[linkIndex]);
        _muldag(force, res);
        //test_force += F(-2.0) * betaOverN * __SU3Generators[8].MulC(force).ImTr();
        _ta(force);
        _mul(force, betaOverN);

        //force is additive
        _add(pForceData[linkIndex], force);
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyCacheIndexT_D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    DOUBLE betaOverN,
    DOUBLE* results
)
{
    intokernal;

    DOUBLE resThisThread = 0.0;
    UINT plaqCountAll = plaqCount * plaqLength;
    for (BYTE i = 0; i < plaqCount; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLength + uiSiteIndex * plaqCountAll];
        //deviceGauge toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
        deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }

        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedIndex[i * plaqLength + j + uiSiteIndex * plaqCountAll];
            //deviceGauge toMul(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
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

        resThisThread += (static_cast<DOUBLE>(_dim<deviceGauge>()) - _retr(toAdd));
    }

    results[uiSiteIndex] = resThisThread * betaOverN;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateOnlyStapleT_D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceGauge* pStapleData)
{
    intokernaldir;

    //Real test_force = F(0.0);
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceGauge res = _makeZero<deviceGauge>();

        //there are 6 staples, each is sum of two plaquttes
        for (int i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

            if (first.NeedToDagger())
            {
                _dagger(toAdd);
            }

            for (int j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
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
        pStapleData[linkIndex] = res;
    }
}

#pragma endregion

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceAndStaple(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, deviceGauge* pStaple, Real betaOverN)
{
    preparethread;
    assert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache);
    _kernelStapleAtSiteGaugeCacheIndex << <block, threads >> > (
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStaple,
        pForce,
        betaOverN);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateOnlyStaple(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pStaple)
{
    preparethread;
    _kernelCalculateOnlyStapleGauge << <block, threads >> > (
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStaple);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN)
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);

    preparethread;
    _kernelPlaqutteEnergyGaugeCacheIndex << <block, threads >> > (
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer
        );

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergyUseClover(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN)
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);
    //appGeneral(_T("const %f\n"), 3.0 * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite);
    preparethread;
    _kernelPlaqutteEnergyGauge_UseClover << <block, threads >> > (
        byFieldId,
        deviceData,
        matrixN * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergyUsingStable(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN, const deviceGauge* pStaple)
{
    preparethread;
    _kernelPlaqutteEnergyUsingStableGauge << <block, threads >> > (
        deviceData,
        pStaple,
        matrixN * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        betaOverN,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateForceAndStaple_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, deviceGauge* pStaple, Real betaOverN)
{
    preparethread;
    assert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache);
    _kernelStapleAtSiteCacheIndexT_D << <block, threads >> > (
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStaple,
        pForce,
        betaOverN);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeKernel<deviceGauge, matrixN>::CalculateOnlyStaple_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pStaple)
{
    preparethread;
    _kernelCalculateOnlyStapleT_D << <block, threads >> > (
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStaple);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeKernel<deviceGauge, matrixN>::CalculatePlaqutteEnergy_D(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN)
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);

    preparethread;
    _kernelPlaqutteEnergyCacheIndexT_D << <block, threads >> > (
        byFieldId,
        deviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer
        );

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

#pragma endregion

template class CFieldGaugeKernel<CLGComplex, 1>;
template class CFieldGaugeKernel<deviceSU2, 2>;
template class CFieldGaugeKernel<deviceSU4, 4>;
template class CFieldGaugeKernel<deviceSU5, 5>;
template class CFieldGaugeKernel<deviceSU6, 6>;
template class CFieldGaugeKernel<deviceSU7, 7>;
template class CFieldGaugeKernel<deviceSU8, 8>;

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================