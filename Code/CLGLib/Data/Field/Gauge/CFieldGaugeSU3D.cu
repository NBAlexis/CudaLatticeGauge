//=============================================================================
// FILENAME : CFieldGaugeSU3D.cu
// 
// DESCRIPTION:
//
// REVISION:
//  [05/17/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldGaugeSU3D)

#pragma region Kernels

/**
* Initial SU3 Field with a value
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialSU3Generator_D(deviceSU3 *pDevicePtr)
{
    deviceSU3 zero = deviceSU3::makeSU3Zero();

    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, idir))
        {
            pDevicePtr[uiLinkIndex] = zero;
        }
        else
        {
            pDevicePtr[uiLinkIndex] = deviceSU3::makeSU3RandomGenerator(_deviceGetFatIndex(uiSiteIndex, idir + 1));
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteSU3CacheIndex_D(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSU3 *pStapleData, //can be NULL
    deviceSU3 *pForceData,
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
        deviceSU3 res = deviceSU3::makeSU3Zero();

        //there are 6 staples, each is sum of two plaquttes
        for (BYTE i = 0; i < plaqCount; ++i)
        {
            BYTE diricCount = 0;
            const SIndex& first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            if (first.IsDirichlet())
            {
                ++diricCount;
            }
            //deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
            deviceSU3 toAdd(_deviceGetGaugeBCSU3(pDeviceData, first));

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (BYTE j = 1; j < plaqLengthm1; ++j)
            {
                const SIndex& nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                if (nextlink.IsDirichlet())
                {
                    ++diricCount;
                }
                //deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);
                deviceSU3 toMul(_deviceGetGaugeBCSU3(pDeviceData, nextlink));

                if (nextlink.NeedToDagger())
                {
                    toAdd.MulDagger(toMul);
                }
                else
                {
                    toAdd.Mul(toMul);
                }
            }
            if (diricCount < plaqLength - 1)
            {
                // If more than 3(including 3) of the edges are Dirichlet, 
                // the plaqutte dose NOT exist.
                res.Add(toAdd);
            }
        }
        if (NULL != pStapleData)
        {
            pStapleData[linkIndex] = res;
        }

        //staple calculated
        deviceSU3 force(pDeviceData[linkIndex]);
        force.MulDagger(res);
        //test_force += F(-2.0) * betaOverN * __SU3Generators[8].MulC(force).ImTr();
        force.Ta();
        force.MulReal(betaOverN);

        //force is additive
        pForceData[linkIndex].Add(force);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergySU3CacheIndex_D(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN,
    DOUBLE* results
#else
    Real betaOverN,
    Real* results
#endif
)
{
    intokernal;

    Real resThisThread = F(0.0);
    UINT plaqCountAll = plaqCount * plaqLength;
    for (BYTE i = 0; i < plaqCount; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLength + uiSiteIndex * plaqCountAll];
        //deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
        deviceSU3 toAdd(_deviceGetGaugeBCSU3(pDeviceData, first));

        if (first.NeedToDagger())
        {
            toAdd.Dagger();
        }

        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedIndex[i * plaqLength + j + uiSiteIndex * plaqCountAll];
            //deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
            deviceSU3 toMul(_deviceGetGaugeBCSU3(pDeviceData, first));
            if (first.NeedToDagger())
            {
                toAdd.MulDagger(toMul);
            }
            else
            {
                toAdd.Mul(toMul);
            }
        }

#if _CLG_DEBUG
        Real reTr = toAdd.ReTr();
        assert(reTr > -F(1.50001));
        assert(reTr < F(3.00001));
        //printf("  ---- energy: thread=%d, res=%f\n", uiSiteIndex, reTr);
#endif
        resThisThread += (F(3.0) - toAdd.ReTr());
    }

    results[uiSiteIndex] = resThisThread * betaOverN;

    //printf("  ---- energy: thread=%d, res=%f\n", uiSiteIndex, results[uiSiteIndex]);
}


__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateOnlyStaple_D(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSU3 *pStapleData)
{
    intokernaldir;

    //Real test_force = F(0.0);
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU3 res = deviceSU3::makeSU3Zero();

        //there are 6 staples, each is sum of two plaquttes
        for (int i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            //deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
            deviceSU3 toAdd(_deviceGetGaugeBCSU3(pDeviceData, first));

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (int j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                //deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);
                deviceSU3 toMul(_deviceGetGaugeBCSU3(pDeviceData, nextlink));

                if (nextlink.NeedToDagger())
                {
                    toAdd.MulDagger(toMul);
                }
                else
                {
                    toAdd.Mul(toMul);
                }
            }
            res.Add(toAdd);
        }
        pStapleData[linkIndex] = res;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelExpMultSU3RealQ_D(
    const deviceSU3 * __restrict__ pMyDeviceData,
    Real a,
    deviceSU3 *pU)
{
    intokernalInt4;
    const UINT uiDir = _DC_Dir;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE idir = 0; idir < uiDir; ++idir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, idir))
        {
            UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
            deviceSU3 expP = pMyDeviceData[linkIndex].QuickExp(a);
            expP.Mul(pU[linkIndex]);
            pU[linkIndex] = expP;
        }
        else
        {
            UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
            pU[linkIndex] = deviceSU3::makeSU3Id();
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelExpMultSU3Real_D(
    const deviceSU3 * __restrict__ pMyDeviceData,
    Real a,
    deviceSU3 *pU,
    BYTE prec)
{
    intokernalInt4;
    const UINT uiDir = _DC_Dir;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE idir = 0; idir < uiDir; ++idir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, idir))
        {
            UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
            deviceSU3 expP = pMyDeviceData[linkIndex].ExpReal(a, prec);
            expP.Mul(pU[linkIndex]);
            pU[linkIndex] = expP;
        }
        else
        {
            UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
            pU[linkIndex] = deviceSU3::makeSU3Id();
        }
    }
}

/**
* Trace (P^2)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateKinematicEnergySU3_D(const deviceSU3 * __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE* results
#else
    Real* results
#endif
)
{
    intokernalInt4;
    const UINT uiDir = _DC_Dir;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    Real resThisThread = F(0.0);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, idir))
        {
            UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
            resThisThread += pDeviceData[linkIndex].DaggerMulC(pDeviceData[linkIndex]).ReTr();
        }
    }
    results[uiSiteIndex] = resThisThread;
}


__global__ void _CLG_LAUNCH_BOUND
_kernelFixBoundarySU3_D(deviceSU3 * pDeviceData)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const UINT uiDir = _DC_Dir;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, idir))
        {
            SIndex idx = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
            pDeviceData[_deviceGetLinkIndex(uiSiteIndex, idir)] =
                ((CFieldBoundaryGaugeSU3*)__boundaryFieldPointers[1])->m_pDeviceData
                [
                    __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idir
                ];
        }

    }
}


/**
 * iA = U.TA() 
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToIA_D(
    deviceSU3* pDeviceData)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex].Ta();
        }
    }
}

#if _CLG_DEBUG
__global__ void _CLG_LAUNCH_BOUND_HALF
#else
__global__ void _CLG_LAUNCH_BOUND
#endif
_kernelTransformToIALog_D(
    deviceSU3* pDeviceData)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].Log();
        }
    }
}

/**
 * U = exp(A)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToU_D(
    deviceSU3* pDeviceData)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            //pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].ExpReal(F(1.0), 8);
            pDeviceData[uiLinkIndex] = (0 == _DC_ExpPrecision)
                ? pDeviceData[uiLinkIndex].QuickExp(F(1.0))
                : pDeviceData[uiLinkIndex].ExpReal(F(1.0), _DC_ExpPrecision);
        }
    }
}

#if _CLG_DEBUG
__global__ void _CLG_LAUNCH_BOUND_HALF
#else
__global__ void _CLG_LAUNCH_BOUND
#endif
_kernelTransformToU_DLog(
    deviceSU3* pDeviceData)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            //pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].ExpReal(F(1.0), 8);
            pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].StrictExp();
        }
    }
}

#pragma endregion


void CFieldGaugeSU3D::MakeRandomGenerator()
{
    preparethread;
    _kernelInitialSU3Generator_D << <block, threads >> > (m_pDeviceData);
}

/**
* (1) calculate staples
* (2) calculate force(additive)
* (3) calculate energy
*/
void CFieldGaugeSU3D::CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, Real betaOverN) const
{
    if (NULL == pForce || EFT_GaugeSU3 != pForce->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: force field is not SU3");
        return;
    }
    if (NULL != pStable && EFT_GaugeSU3 != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stape field is not SU3");
        return;
    }

    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    CFieldGaugeSU3* pStableSU3 = NULL == pStable ? NULL : dynamic_cast<CFieldGaugeSU3*>(pStable);

    preparethread;

    assert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache);

    _kernelStapleAtSiteSU3CacheIndex_D << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData,
        betaOverN);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeSU3D::CalculatePlaqutteEnergy(DOUBLE betaOverN) const
#else
Real CFieldGaugeSU3D::CalculatePlaqutteEnergy(Real betaOverN) const
#endif
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);

    preparethread;
    _kernelPlaqutteEnergySU3CacheIndex_D << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeSU3D::CalculateKinematicEnergy() const
#else
Real CFieldGaugeSU3D::CalculateKinematicEnergy() const
#endif
{
    preparethread;
    _kernelCalculateKinematicEnergySU3_D << <block, threads >> > (m_pDeviceData, _D_RealThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

void CFieldGaugeSU3D::CalculateOnlyStaple(CFieldGauge* pStaple) const
{
    if (NULL == pStaple || EFT_GaugeSU3 != pStaple->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stable field is not SU3");
        return;
    }
    CFieldGaugeSU3* pStapleSU3 = dynamic_cast<CFieldGaugeSU3*>(pStaple);

    preparethread;
    _kernelCalculateOnlyStaple_D << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStapleSU3->m_pDeviceData);
}

void CFieldGaugeSU3D::ExpMult(Real a, CField* U) const
{
    if (NULL == U || EFT_GaugeSU3 != U->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: U field is not SU3");
        return;
    }

    CFieldGaugeSU3* pUField = dynamic_cast<CFieldGaugeSU3*>(U);

    preparethread;
    if (0 == _HC_ExpPrecision)
    {
        _kernelExpMultSU3RealQ_D << < block, threads >> > (m_pDeviceData, a, pUField->m_pDeviceData);
    }
    else
    {
        _kernelExpMultSU3Real_D << < block, threads >> > (m_pDeviceData, a, pUField->m_pDeviceData, static_cast<BYTE>(_HC_ExpPrecision));
    }
    
}

void CFieldGaugeSU3D::FixBoundary()
{
    appDetailed(_T("CFieldGaugeSU3::FixBoundary()\n"));

    preparethread;
    _kernelFixBoundarySU3_D << <block, threads >> > (m_pDeviceData);

    //DebugPrintMe();
    //INT i = 0;
}

void CFieldGaugeSU3D::TransformToIA()
{
    preparethread;
    if (0 == _HC_ALog)
    {
        _kernelTransformToIA_D << <block, threads >> > (m_pDeviceData);
    }
    else
    {
        _kernelTransformToIALog_D << <block, threads >> > (m_pDeviceData);
    }
}

void CFieldGaugeSU3D::TransformToU()
{
    preparethread;
    if (0 == _HC_ALog)
    {
        _kernelTransformToU_D << <block, threads >> > (m_pDeviceData);
    }
    else
    {
        _kernelTransformToU_DLog << <block, threads >> > (m_pDeviceData);
    }
    
}

void CFieldGaugeSU3D::CopyTo(CField* pTarget) const
{
    CFieldGaugeSU3::CopyTo(pTarget);
}

CCString CFieldGaugeSU3D::GetInfos(const CCString &tab) const
{
    CCString sRet = tab + _T("Name : CFieldGaugeSU3 Dirichlet\n");
    sRet = sRet + CFieldGauge::GetInfos(tab);
    SSmallInt4 boundary = appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(m_byFieldId);
    sRet = sRet + tab + appToString(boundary) + _T("\n");

    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================