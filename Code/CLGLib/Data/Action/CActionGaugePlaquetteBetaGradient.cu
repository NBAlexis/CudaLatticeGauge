//=============================================================================
// FILENAME : CActionGaugePlaquetteGradient.cpp
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [08/15/2022 nbale]
//=============================================================================
#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteGradient)

#pragma region kernel

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergySU3_UseCloverGradient(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    const DOUBLE* __restrict__ fBetaOverN,
    DOUBLE* results
#else
    const Real* __restrict__ fBetaOverN,
    Real* results
#endif
)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        //printf("x=%d y=%d z=%d w=%d\n", sSite4.x, sSite4.y, sSite4.z, sSite4.w);
        return;
    }

#if !_CLG_DOUBLEFLOAT
    DOUBLE fRes = 0.0;
#else
    Real fRes = F(0.0);
#endif
    for (BYTE byDir1 = 0; byDir1 < _DC_Dir; ++byDir1)
    {
        for (BYTE byDir2 = byDir1 + 1; byDir2 < _DC_Dir; ++byDir2)
        {
            fRes += _deviceCloverRetr(pDeviceData, sSite4, __bi(sSite4), byDir1, byDir2, byFieldId);
        }
    }
#if !_CLG_DOUBLEFLOAT
    fRes = 18.0 - 0.25 * fRes;
#else
    fRes = F(18.0) - F(0.25) * fRes;
#endif
    results[uiSiteIndex] = fRes * fBetaOverN[sSite4.z];
}

__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteSU3CacheIndexGradient(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSU3* pStapleData, //can be NULL
    deviceSU3* pForceData,
#if !_CLG_DOUBLEFLOAT
    const DOUBLE* __restrict__ betaOverN
#else
    const Real* __restrict__ betaOverN
#endif
)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    //UBOOL bBoundStartFromDirichlet = FALSE;
    //if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].IsDirichlet())
    //{
    //    bBoundStartFromDirichlet = TRUE;
    //}

    //Real test_force = F(0.0);
    //betaOverN = betaOverN * F(-0.5);
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, static_cast<BYTE>(idir)))
        {
            continue;
        }

        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU3 res = deviceSU3::makeSU3Zero();

        //there are 6 staples, each is sum of two plaquttes
        SBYTE z1 = 0;
        SBYTE z2 = 0;
        BYTE zNotFound = 1;
        for (INT i = 0; i < plaqCount; ++i)
        {
            BYTE diricCount = 0;
            const SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            if (first.IsDirichlet())
            {
                ++diricCount;
            }
            deviceSU3 toAdd(_deviceGetGaugeBCSU3(pDeviceData, first));

            const SSmallInt4 firstn = __deviceSiteIndexToInt4(first.m_uiSiteIndex);
            z1 = firstn.z;
            if (2 == first.m_byDir)
            {
                z2 = z1 + 1;
                zNotFound = 0;
            }

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                const SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
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

                if (zNotFound && nextlink.m_byDir == 2)
                {
                    zNotFound = 0;
                    const SSmallInt4 nextlinkn = __deviceSiteIndexToInt4(nextlink.m_uiSiteIndex);
                    z1 = nextlinkn.z;
                    z2 = z1 + 1;
                }
            }

            if (diricCount < plaqLengthm1)
            {
                // If more than 3(including 3) of the edges are Dirichlet, 
                // the plaqutte dose NOT exist.
                if (zNotFound)
                {
                    //zNotFound means this is a plaqutte not have z bounds
                    toAdd.MulReal(betaOverN[z1] * F(-0.5));
                }
                else
                {
                    if (0 == z1)
                    {
                        if (diricCount > 0)
                        {
                            z1 = z2;
                        }
                    }

                    if (z2 >= _DC_Lz)
                    {
                        if (diricCount > 0)
                        {
                            z2 = _DC_Lz - 1;
                        }
                        else
                        {
                            z2 = 0;
                        }
                    }
                    toAdd.MulReal(F(-0.25) * (betaOverN[z1] + betaOverN[z2]));
                }
                if (diricCount > 0)
                {
                    toAdd.MulReal(F(0.5));
                }
                res.Add(toAdd);
            }
            //else
            //{
            //    printf("do we have this?\n");
            //}
            //printf("diricCount=%d\n", diricCount);
        }
        if (NULL != pStapleData)
        {
            pStapleData[linkIndex] = res;
        }

        //staple calculated
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, static_cast<BYTE>(idir)))
        {
            deviceSU3 force(pDeviceData[linkIndex]);
            force.MulDagger(res);
            //test_force += F(-2.0) * betaOverN * __SU3Generators[8].MulC(force).ImTr();
            force.Ta();

            //this is the average over 4 cornels, so this is different for different dirs
            //force.MulReal(betaOverN);

            //force is additive
            pForceData[linkIndex].Add(force);
        }
    }
}

#pragma endregion

CActionGaugePlaquetteGradient::CActionGaugePlaquetteGradient()
    : CAction()
    , m_pDeviceBetaArray(NULL)
    , m_uiPlaqutteCount(0)
{
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CActionGaugePlaquetteGradient::CalculatePlaqutteEnergyUseClover(const CFieldGaugeSU3* pGauge) const
#else
Real CActionGaugePlaquetteGradient::CalculatePlaqutteEnergyUseClover(const CFieldGaugeSU3* pGauge) const
#endif
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);
    //pGauge->FixBoundary();
    //pGauge->DebugPrintMe();

    preparethread;
    _kernelPlaqutteEnergySU3_UseCloverGradient << <block, threads >> > (
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        m_pDeviceBetaArray,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

void CActionGaugePlaquetteGradient::CalculateForceAndStaple(const CFieldGaugeSU3* pGauge, CFieldGaugeSU3* pForce) const
{
    preparethread;

    assert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache);

    _kernelStapleAtSiteSU3CacheIndexGradient << <block, threads >> > (
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        NULL,
        pForce->m_pDeviceData,
        m_pDeviceBetaArray);
}

void CActionGaugePlaquetteGradient::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
        if (NULL == pGaugeSU3)
        {
            appCrucial(_T("CActionGaugePlaquetteGradient must be CFieldGaugeSU3!\n"));
            return;
        }

        assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);
        m_fLastEnergy = CalculatePlaqutteEnergyUseClover(pGaugeSU3);
    }
}

void CActionGaugePlaquetteGradient::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

#if !_CLG_DOUBLEFLOAT
    param.FetchValueArrayDOUBLE(_T("Beta"), m_fBetaArray);
    for (INT i = 0; i < _HC_Lzi; ++i)
    {
        if (i < m_fBetaArray.Num())
        {
            m_fBetaArray[i] = m_fBetaArray[i] / static_cast<DOUBLE>(_HC_SUN);
        }
        else
        {
            if (0 == i)
            {
                m_fBetaArray.AddItem(5.0 / static_cast<DOUBLE>(_HC_SUN));
            }
            else
            {
                m_fBetaArray.AddItem(m_fBetaArray.Num() - 1);
            }
        }
    }
#else
    param.FetchValueArrayReal(_T("Beta"), m_fBetaArray);
    for (INT i = 0; i < _HC_Lzi; ++i)
    {
        if (i < m_fBetaArray.Num())
        {
            m_fBetaArray[i] = m_fBetaArray[i] / static_cast<Real>(_HC_SUN);
        }
        else
        {
            if (0 == i)
            {
                m_fBetaArray.AddItem(F(5.0) / static_cast<Real>(_HC_SUN));
            }
            else
            {
                m_fBetaArray.AddItem(m_fBetaArray.Num() - 1);
            }
        }
    }
#endif
    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

#if !_CLG_DOUBLEFLOAT
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceBetaArray, sizeof(DOUBLE) * _HC_Lz));
    checkCudaErrors(cudaMemcpy(m_pDeviceBetaArray, m_fBetaArray.GetData(), sizeof(DOUBLE) * _HC_Lz, cudaMemcpyHostToDevice));
#else
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceBetaArray, sizeof(Real) * _HC_Lz));
    checkCudaErrors(cudaMemcpy(m_pDeviceBetaArray, m_fBetaArray.GetData(), sizeof(Real) * _HC_Lz, cudaMemcpyHostToDevice));
#endif


}

void CActionGaugePlaquetteGradient::SetBeta(const TArray<DOUBLE>& fBeta)
{
    m_fBetaArray = fBeta;
#if !_CLG_DOUBLEFLOAT
    for (INT i = 0; i < _HC_Lzi; ++i)
    {
        if (i < m_fBetaArray.Num())
        {
            m_fBetaArray[i] = m_fBetaArray[i] / static_cast<DOUBLE>(_HC_SUN);
        }
        else
        {
            if (0 == i)
            {
                m_fBetaArray.AddItem(5.0 / static_cast<DOUBLE>(_HC_SUN));
            }
            else
            {
                m_fBetaArray.AddItem(m_fBetaArray.Num() - 1);
            }
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceBetaArray, m_fBetaArray.GetData(), sizeof(DOUBLE) * _HC_Lz, cudaMemcpyHostToDevice));
#else
    for (INT i = 0; i < _HC_Lzi; ++i)
    {
        if (i < m_fBetaArray.Num())
        {
            m_fBetaArray[i] = m_fBetaArray[i] / static_cast<Real>(_HC_SUN);
        }
        else
        {
            if (0 == i)
            {
                m_fBetaArray.AddItem(F(5.0) / static_cast<Real>(_HC_SUN));
            }
            else
            {
                m_fBetaArray.AddItem(m_fBetaArray.Num() - 1);
            }
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceBetaArray, m_fBetaArray.GetData(), sizeof(Real)* _HC_Lz, cudaMemcpyHostToDevice));
#endif

}

UBOOL CActionGaugePlaquetteGradient::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    if (NULL == pGaugeSU3 || NULL == pForceSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteGradient must be CFieldGaugeSU3!\n"));
        return FALSE;
    }

    CalculateForceAndStaple(pGaugeSU3, pForceSU3);
    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
DOUBLE CActionGaugePlaquetteGradient::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    if (NULL == pGaugeSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteGradient must be CFieldGaugeSU3!\n"));
        return m_fNewEnergy;
    }
    m_fNewEnergy = CalculatePlaqutteEnergyUseClover(pGaugeSU3);
    return m_fNewEnergy;
}

CCString CActionGaugePlaquetteGradient::GetInfos(const CCString &tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("Beta : ");
    for (INT i = 0; i < _HC_Lzi; ++i)
    {
        sRet = sRet + appToString(m_fBetaArray[i]) + _T(", ");
    }
    sRet = sRet + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================