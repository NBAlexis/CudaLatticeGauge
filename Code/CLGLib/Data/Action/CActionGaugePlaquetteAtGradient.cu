//=============================================================================
// FILENAME : CActionGaugePlaquetteAtGradient.cu
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [mm/dd/yy]
//  [08/15/2022 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionGaugePlaquetteAtGradient.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteAtGradient)

#pragma region kernel

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergySU3_UseCloverAtGradient(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const DOUBLE* __restrict__ fXiList,
    const DOUBLE fBetaOverN,
    DOUBLE* results
)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        return;
    }

    DOUBLE fRes = 0.0;
    const DOUBLE fXi = fXiList[sSite4.z];
    for (BYTE byDir1 = 0; byDir1 < _DC_Dir; ++byDir1)
    {
        for (BYTE byDir2 = byDir1 + 1; byDir2 < _DC_Dir; ++byDir2)
        {
            if (3 == byDir2)
            {
                fRes += (3.0 - 0.25 * _deviceCloverRetrT(pDeviceData, sSite4, __bi(sSite4), byDir1, byDir2, byFieldId)) / fXi;
            }
            else
            {
                fRes += (3.0 - 0.25 * _deviceCloverRetrT(pDeviceData, sSite4, __bi(sSite4), byDir1, byDir2, byFieldId)) * fXi;
            }
        }
    }
    results[uiSiteIndex] = fRes * fBetaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteSU3CacheIndexAtGradient(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceSU3* pStapleData, //can be NULL
    deviceSU3* pForceData,
    DOUBLE fBetaOverN,
    const DOUBLE* __restrict__ fXiList
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
#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1;
#endif

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, static_cast<BYTE>(idir)))
        {
            continue;
        }

        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU3 res = deviceSU3::makeSU3Zero();

        //there are 6 staples, each is sum of two plaquttes
        SCHAR z1 = 0;
        SCHAR z2 = 0;
        BYTE zNotFound = 1;
        BYTE hasT = 0;
        for (INT i = 0; i < plaqCountPerLink; ++i)
        {
            BYTE diricCount = 0;
            const SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAllLink];
            if (first.IsDirichlet())
            {
                ++diricCount;
            }
            deviceSU3 toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

            const SSmallInt4 firstn = __deviceSiteIndexToInt4(first.m_uiSiteIndex);
            z1 = firstn.z;
            if (2 == first.m_byDir)
            {
                z2 = z1 + 1;
                zNotFound = 0;
            }

            if (3 == first.m_byDir)
            {
                hasT = 1;
            }

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                const SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAllLink];
                if (nextlink.IsDirichlet())
                {
                    ++diricCount;
                }
                //deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);
                deviceSU3 toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink));

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

                if (3 == nextlink.m_byDir)
                {
                    hasT = 1;
                }
            }

            if (diricCount < plaqLengthm1)
            {
                // If more than 3(including 3) of the edges are Dirichlet, 
                // the plaqutte dose NOT exist.
                if (zNotFound)
                {
                    //zNotFound means this is a plaqutte not have z bounds
                    if (hasT)
                    {
                        toAdd.MulReal(fBetaOverN * F(-0.5) / fXiList[z1]);
                    }
                    else
                    {
                        toAdd.MulReal(fBetaOverN * F(-0.5) * fXiList[z1]);
                    }
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
                    if (hasT)
                    {
                        toAdd.MulReal(fBetaOverN * F(-0.25) * (__rcp(fXiList[z1]) + __rcp(fXiList[z2])));
                    }
                    else
                    {
                        toAdd.MulReal(fBetaOverN * F(-0.25) * (fXiList[z1] + fXiList[z2]));
                    }
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
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, static_cast<BYTE>(idir)))
        {
            //deviceSU3 force(pDeviceData[linkIndex]);
            //force.MulDagger(res);
            //test_force += F(-2.0) * betaOverN * __SU3Generators[8].MulC(force).ImTr();
            //force.Ta();

            //this is the average over 4 cornels, so this is different for different dirs
            //force.MulReal(betaOverN);

            //force is additive
            pForceData[linkIndex].Add(res);
        }
    }
}

#pragma endregion

CActionGaugePlaquetteAtGradient::CActionGaugePlaquetteAtGradient()
    : CAction()
    , m_pDeviceXiArray(NULL)
    , m_uiPlaqutteCount(0)
{
}

DOUBLE CActionGaugePlaquetteAtGradient::CalculatePlaqutteEnergyUseClover(const CFieldGaugeSU3* pGauge) const
{
    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId]);
    //pGauge->FixBoundary();
    //pGauge->DebugPrintMe();

    preparethread;
    _LAUNCH_KERNEL(_kernelPlaqutteEnergySU3_UseCloverAtGradient, block, threads, 
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        m_pDeviceXiArray,
        m_fBetaOverN,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

void CActionGaugePlaquetteAtGradient::CalculateForceAndStaple(const CFieldGaugeSU3* pGauge, CFieldGaugeSU3* pForce) const
{
    preparethread;

    appAssert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache[pGauge->m_byFieldId]);
    deviceSU3* emptystaple = NULL;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelStapleAtSiteSU3CacheIndexAtGradient, block, threads,
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[pGauge->m_byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        emptystaple,
        pForce->m_pDeviceData,
        m_fBetaOverN,
        m_pDeviceXiArray);
#else
    _LAUNCH_KERNEL(_kernelStapleAtSiteSU3CacheIndexAtGradient, block, threads, 
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[pGauge->m_byFieldId],
        emptystaple,
        pForce->m_pDeviceData,
        m_fBetaOverN,
        m_pDeviceXiArray);
#endif
}

void CActionGaugePlaquetteAtGradient::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
        if (NULL == pGaugeSU3)
        {
            appCrucial(_T("CActionGaugePlaquetteGradient must be CFieldGaugeSU3!\n"));
            return;
        }

        appAssert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGauge->m_byFieldId]);
        m_fLastEnergy = CalculatePlaqutteEnergyUseClover(pGaugeSU3);
        //I9: the cached "before" energy must be the GLOBAL action (the clover
        //energy above sums only this rank's local sub-lattice), or HMC sees a
        //constant ΔH offset. No-op on a lone rank.
        appGlobalSum(m_fLastEnergy);
    }
}

void CActionGaugePlaquetteAtGradient::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    param.FetchValueArrayDOUBLE(_T("Xi"), m_fXiArray);
    for (INT i = 0; i < _HC_Lzi; ++i)
    {
        if (i >= m_fXiArray.Num())
        {
            if (0 == i)
            {
                m_fXiArray.AddItem(1.0);
            }
            else
            {
                m_fXiArray.AddItem(m_fXiArray[m_fXiArray.Num() - 1]);
            }
        }
    }

    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceXiArray, sizeof(DOUBLE) * _HC_Lz));
    checkCudaErrors(cudaMemcpy(m_pDeviceXiArray, m_fXiArray.GetData(), sizeof(DOUBLE) * _HC_Lz, cudaMemcpyHostToDevice));
}

void CActionGaugePlaquetteAtGradient::SetXiList(const TArray<DOUBLE>& fXi)
{
    m_fXiArray = fXi;
    for (INT i = 0; i < _HC_Lzi; ++i)
    {
        if (i >= m_fXiArray.Num())
        {
            if (0 == i)
            {
                m_fXiArray.AddItem(1.0);
            }
            else
            {
                m_fXiArray.AddItem(m_fXiArray[m_fXiArray.Num() - 1]);
            }
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceXiArray, m_fXiArray.GetData(), sizeof(DOUBLE) * _HC_Lz, cudaMemcpyHostToDevice));
}

UBOOL CActionGaugePlaquetteAtGradient::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
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
DOUBLE CActionGaugePlaquetteAtGradient::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    if (NULL == pGaugeSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteAtGradient must be CFieldGaugeSU3!\n"));
        return m_fNewEnergy;
    }
    m_fNewEnergy = CalculatePlaqutteEnergyUseClover(pGaugeSU3);
    //I9: local partial sum -> global action (see PrepareForHMCSingleField).
    appGlobalSum(m_fNewEnergy);
    return m_fNewEnergy;
}

CCString CActionGaugePlaquetteAtGradient::GetInfos(const CCString &tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("Beta : ") + appToString(m_fXiArray) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================