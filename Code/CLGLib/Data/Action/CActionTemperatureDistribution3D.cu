//=============================================================================
// FILENAME : CActionTemperatureDistribution.cu
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [mm/dd/yy]
//  [08/15/2022 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionTemperatureDistribution.h"
#include "Data/Field/Boson/CFieldBosonReal.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionTemperatureDistribution)

#pragma region kernel

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergySU3_UseClover_TemperatureDistribution(
    BYTE byGaugeFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const Real* __restrict__ pBoson,
    const DOUBLE fBetaOverN,
    DOUBLE* results
)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //if (__idx->m_pDeviceIndexPositionToSIndex[byGaugeFieldId][uiBigIdx].IsDirichlet())
    //{
    //    results[uiSiteIndex] = F(0.0);
    //    return;
    //}

    DOUBLE fRes = 0.0;
    const Real fXi = pBoson[uiSiteIndex];
    for (BYTE byDir1 = 0; byDir1 < _DC_Dir; ++byDir1)
    {
        for (BYTE byDir2 = byDir1 + 1; byDir2 < _DC_Dir; ++byDir2)
        {
            if (3 == byDir2)
            {
                fRes += (3.0 - 0.25 * _deviceCloverRetrT(pDeviceData, sSite4, uiBigIdx, byDir1, byDir2, byGaugeFieldId)) / fXi;
            }
            else
            {
                fRes += (3.0 - 0.25 * _deviceCloverRetrT(pDeviceData, sSite4, uiBigIdx, byDir1, byDir2, byGaugeFieldId)) * fXi;
            }
        }
    }
    //printf("%f", fXi);
    results[uiSiteIndex] = fRes * fBetaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergySU3_BosonForceUseClover_TemperatureDistribution(
    BYTE byGaugeFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const Real* __restrict__ pBoson,
    const DOUBLE fBetaOverN,
    Real* bosonForce
)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //if (__idx->m_pDeviceIndexPositionToSIndex[byGaugeFieldId][uiBigIdx].IsDirichlet())
    //{
    //    results[uiSiteIndex] = F(0.0);
    //    return;
    //}

    Real fRes = F(0.0);
    Real fXi = pBoson[uiSiteIndex];
    fXi = fXi * fXi;
    for (BYTE byDir1 = 0; byDir1 < _DC_Dir; ++byDir1)
    {
        for (BYTE byDir2 = byDir1 + 1; byDir2 < _DC_Dir; ++byDir2)
        {
            if (3 == byDir2)
            {
                fRes -= (F(3.0) - F(0.25) * _deviceCloverRetrT(pDeviceData, sSite4, uiBigIdx, byDir1, byDir2, byGaugeFieldId)) / fXi;
            }
            else
            {
                fRes += (F(3.0) - F(0.25) * _deviceCloverRetrT(pDeviceData, sSite4, uiBigIdx, byDir1, byDir2, byGaugeFieldId));
            }
        }
    }
    //printf("%f", fXi);
    bosonForce[uiSiteIndex] -= fRes * fBetaOverN * F(0.5);
}

/**
* we run over all sites and all links of a site
* therefore the force of boson is calculated at the same time
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteSU3CacheIndex_TemperatureDistribution(
    BYTE byGaugeFieldId,
    BYTE byBosonFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const Real* __restrict__ pBoson,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceSU3* pStapleData, //can be NULL
    deviceSU3* pForceData,
    //Real* pForceBoson,
    DOUBLE fBetaOverN
)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
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

    //count each site index twice
    UINT eightSites[8];
    UINT maxSiteCount = 0;
    eightSites[0] = uiSiteIndex;

    for (SCHAR idir = 0; idir < uiDir; ++idir)
    {
        //if (__idx->_deviceIsBondOnSurface(uiBigIdx, byGaugeFieldId, static_cast<BYTE>(idir)))
        //{
        //    continue;
        //}

        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU3 res = deviceSU3::makeSU3Zero();

        //there are 6 staples, each is sum of two plaquttes
        SSmallInt4 sSite4Move = sSite4;
        //if (0 == uiSiteIndex)
        //{
        //    printf("first idx = %d, %d, %d, %d\n", sSite4Move.x, sSite4Move.y, sSite4Move.z, sSite4Move.w);
        //}
        _deviceSmallInt4Offset(sSite4Move, idir + 1);
        //if (0 == uiSiteIndex)
        //{
        //    printf("second idx = %d, %d, %d, %d\n", sSite4Move.x, sSite4Move.y, sSite4Move.z, sSite4Move.w);
        //}
        eightSites[1] = __idx->m_pDeviceIndexPositionToSIndex[byBosonFieldId][__bi(sSite4Move)].m_uiSiteIndex;

        for (INT i = 0; i < plaqCountPerLink; ++i)
        {
            BYTE hasT = 0;
            //BYTE diricCount = 0;
            const SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAllLink];
            //if (first.IsDirichlet())
            //{
            //    ++diricCount;
            //}
            deviceSU3 toAdd(_deviceGetGaugeBCT(byGaugeFieldId, pDeviceData, first));

            eightSites[2] = first.m_uiSiteIndex;
            sSite4Move = __deviceSiteIndexToInt4(first.m_uiSiteIndex);
            _deviceSmallInt4Offset(sSite4Move, first.m_byDir + 1);
            eightSites[3] = __idx->m_pDeviceIndexPositionToSIndex[byBosonFieldId][__bi(sSite4Move)].m_uiSiteIndex;
            maxSiteCount = 4;

            if (3 == first.m_byDir)
            {
                hasT = 1;
            }

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            //Here we assumes 3 == plaqLengthm1
            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                const SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAllLink];
                //if (nextlink.IsDirichlet())
                //{
                //    ++diricCount;
                //}
                //deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);
                deviceSU3 toMul(_deviceGetGaugeBCT(byGaugeFieldId, pDeviceData, nextlink));

                if (nextlink.NeedToDagger())
                {
                    toAdd.MulDagger(toMul);
                }
                else
                {
                    toAdd.Mul(toMul);
                }

                if (j < 3) //0,1,2
                {
                    eightSites[2 * j + 2] = nextlink.m_uiSiteIndex;
                    sSite4Move = __deviceSiteIndexToInt4(nextlink.m_uiSiteIndex);
                    _deviceSmallInt4Offset(sSite4Move, nextlink.m_byDir + 1);
                    eightSites[2 * j + 3] = __idx->m_pDeviceIndexPositionToSIndex[byBosonFieldId][__bi(sSite4Move)].m_uiSiteIndex;

                    maxSiteCount += 2;
                }
                else
                {
                    printf("not support plaq with more than 4 links\n");
                }


                if (3 == nextlink.m_byDir)
                {
                    hasT = 1;
                }
            }

            //deviceSU3 forceBoson(pDeviceData[linkIndex]);
            //forceBoson.MulDagger(toAdd);
            //const Real fForceBoson = forceBoson.ReTr(); //we don't care the direction of a closed loop
            Real fFactor = F(0.0);
            const Real fRepeatFactor = __rcp(maxSiteCount);
            for (UINT j = 0; j < maxSiteCount; ++j)
            {
                //if we run over all staple, every plaqutee is considered four times.
                //N * 4 * 6 plaquttes are calculatted
                //N * 4 * 3 / 2 plaquttes in all for this lattice

                //for each site, 4 * 6 plaquttes should contribute
                //if we run over all staple, the contribution is added for N * 4 * 6 * 4 times

                //if (NULL != pForceBoson) //Whether boson is dynamic is decided by whether this is NULL
                //{
                //    if (hasT)
                //    {
                //        pForceBoson[eightSites[j]] += F(0.125) * fBetaOverN * fForceBoson * fRepeatFactor / (pBoson[eightSites[j]] * pBoson[eightSites[j]]);
                //    }
                //    else
                //    {
                //        pForceBoson[eightSites[j]] -= F(0.125) * fBetaOverN * fForceBoson * fRepeatFactor;
                //    }
                //}

                if (hasT)
                {
                    fFactor = fFactor + __rcp(pBoson[eightSites[j]]);
                }
                else
                {
                    fFactor = fFactor + pBoson[eightSites[j]];
                }

                //if (0 == uiSiteIndex)
                //{
                //    SSmallInt4 siteshow = __deviceSiteIndexToInt4(eightSites[j]);
                //    printf("j=%d, eightsitej=(%d %d %d %d)\n", j, siteshow.x, siteshow.y, siteshow.z, siteshow.w);
                //}
            }
            toAdd.MulReal(fBetaOverN * F(-0.5) * fRepeatFactor * fFactor);
            res.Add(toAdd);
        }
        if (NULL != pStapleData)
        {
            pStapleData[linkIndex] = res;
        }

        //staple calculated
        //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byGaugeFieldId, static_cast<BYTE>(idir)))
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

CActionTemperatureDistribution::CActionTemperatureDistribution()
    : CAction()
    , m_fBosonKinetic(1.0)
    , m_pStaticBosonField(NULL)
{
}

void CActionTemperatureDistribution::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    INT iVal = 0;
    param.FetchValueINT(_T("StaticBosonFieldId"), iVal);
    if (iVal > 0)
    {
        BYTE byFieldId = static_cast<BYTE>(iVal);
        //CField* pField = appGetLattice()->GetFieldById(byFieldId);
        m_pStaticBosonField = dynamic_cast<CFieldBosonReal*>(appGetLattice()->GetFieldById(byFieldId));
    }

    DOUBLE bosonKinetic = 1.0;
    if (param.FetchValueDOUBLE(_T("BosonKinetic"), bosonKinetic))
    {
        m_fBosonKinetic = bosonKinetic;
    }
}

DOUBLE CActionTemperatureDistribution::Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldGauge* const* stapleFields)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

    if (NULL != m_pStaticBosonField)
    {
        if (m_byGaugeFieldIds.Num() < 1)
        {
            appCrucial(_T("CActionTemperatureDistribution needs a gauge field!\n"));
            return 0.0;
        }

        INT iGaugeidx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
        if (iGaugeidx < 0)
        {
            appCrucial(_T("CActionTemperatureDistribution needs a gauge field and a boson field!\n"));
            return 0.0;
        }
        const CFieldGaugeSU3* pGauge = dynamic_cast<const CFieldGaugeSU3*>(gaugeFields[iGaugeidx]);
        if (NULL == pGauge)
        {
            appCrucial(_T("CActionTemperatureDistribution needs a SU3 gauge field and a Real boson field!\n"));
            return 0.0;
        }

        //If it was static, it did not fix boundary
        m_pStaticBosonField->FixBoundary(EFB_Field);

        preparethread;
        _LAUNCH_KERNEL(_kernelPlaqutteEnergySU3_UseClover_TemperatureDistribution, block, threads, 
            pGauge->m_byFieldId,
            pGauge->m_pDeviceData,
            m_pStaticBosonField->m_pDeviceData,
            m_fBetaOverN,
            _D_RealThreadBuffer);

        m_fNewEnergy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        return m_fNewEnergy;
    }

    if (m_byGaugeFieldIds.Num() < 1 || m_byBosonFieldIds.Num() < 1)
    {
        appCrucial(_T("CActionTemperatureDistribution needs a gauge field and a boson field!\n"));
        return 0.0;
    }

    INT iGaugeidx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
    INT iBosonidx = CLatticeData::GetBosonFieldIndexById(bosonNum, bosonFields, m_byBosonFieldIds[0]);
    if (iGaugeidx < 0 || iBosonidx < 0)
    {
        appCrucial(_T("CActionTemperatureDistribution needs a gauge field and a boson field!\n"));
        return 0.0;
    }
    const CFieldGaugeSU3* pGauge = dynamic_cast<const CFieldGaugeSU3*>(gaugeFields[iGaugeidx]);
    const CFieldBosonReal* pBoson = dynamic_cast<const CFieldBosonReal*>(bosonFields[iBosonidx]);
    if (NULL == pGauge || NULL == pBoson)
    {
        appCrucial(_T("CActionTemperatureDistribution needs a SU3 gauge field and a Real boson field!\n"));
        return 0.0;
    }

    //pBoson->DebugPrintMe();

    preparethread;
    _LAUNCH_KERNEL(_kernelPlaqutteEnergySU3_UseClover_TemperatureDistribution, block, threads, 
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        pBoson->m_pDeviceData,
        m_fBetaOverN,
        _D_RealThreadBuffer);

    m_fNewEnergy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    //add kinetic term
    if (abs(m_fBosonKinetic) > _CLG_FLT_MIN_)
    {
        /**
    INT ibosonidx = CLatticeData::GetBosonFieldIndexById(bosonNum, bosonFields, m_byBosonFieldIds[0]);
    const CFieldBoson* bosonfield = bosonFields[ibosonidx];
    CFieldBoson* dboson = dynamic_cast<CFieldBoson*>(appGetLattice()->GetPooledCopy(bosonfield));
    dboson->D(gaugeNum, bosonNum, gaugeFields, bosonFields);
    checkCudaErrors(cudaDeviceSynchronize());
    cuDoubleComplex partiald = bosonfield->Dot(dboson);
    checkCudaErrors(cudaDeviceSynchronize());
    DOUBLE phi2 = bosonfield->GetLength();
    checkCudaErrors(cudaDeviceSynchronize());
    bosonfield->CopyTo(dboson);
    dboson->Mul(bosonfield);
    checkCudaErrors(cudaDeviceSynchronize());
    DOUBLE phi4 = dboson->GetLength();
    checkCudaErrors(cudaDeviceSynchronize());
    m_fLastEnergy = (8.0 + m_fM) * phi2 + m_fLambda * phi4 - partiald.x;
        */
        CFieldBoson* dboson = dynamic_cast<CFieldBoson*>(appGetLattice()->GetPooledCopy(pBoson, _T(__FILE__), __LINE__));
        dboson->D(0, 0, NULL, NULL);
        cuDoubleComplex partiald = pBoson->Dot(dboson);
        DOUBLE phi2 = pBoson->GetLength();
        //pBoson->DebugPrintMe();
        //appGeneral(_T("================="));
        //dboson->DebugPrintMe();
        m_fNewEnergy += (8.0 * phi2 - partiald.x) * m_fBosonKinetic;
        dboson->Return();
    }

    return m_fNewEnergy;
}

UBOOL CActionTemperatureDistribution::CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
    CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
    CFieldGauge* const* stapleFields, ESolverPhase ePhase) const
{
    if (NULL != m_pStaticBosonField)
    {
        if (m_byGaugeFieldIds.Num() < 1)
        {
            appCrucial(_T("CActionTemperatureDistribution needs a gauge field and a boson field!\n"));
            return FALSE;
        }

        INT iGaugeidx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
        if (iGaugeidx < 0)
        {
            appCrucial(_T("CActionTemperatureDistribution needs a gauge field and a boson field!\n"));
            return FALSE;
        }

        const CFieldGaugeSU3* pGauge = dynamic_cast<const CFieldGaugeSU3*>(gaugeFields[iGaugeidx]);
        if (NULL == pGauge)
        {
            appCrucial(_T("CActionTemperatureDistribution needs a SU3 gauge field and a Real boson field!\n"));
            return FALSE;
        }

        CFieldGaugeSU3* pGaugeForce = dynamic_cast<CFieldGaugeSU3*>(gaugeForces[iGaugeidx]);
        CFieldGaugeSU3* pStaple = (NULL != stapleFields) ? dynamic_cast<CFieldGaugeSU3*>(stapleFields[iGaugeidx]) : NULL;
        if (NULL == pGaugeForce)
        {
            appCrucial(_T("CActionTemperatureDistribution both field should be dynamic\n"));
            return FALSE;
        }

        //If it was static, it did not fix boundary
        m_pStaticBosonField->FixBoundary(EFB_Field);

        preparethread;
#if !_CLG_ASSUME_SQUARE_LATTICE
        _LAUNCH_KERNEL(_kernelStapleAtSiteSU3CacheIndex_TemperatureDistribution, block, threads, 
            pGauge->m_byFieldId,
            m_pStaticBosonField->m_byFieldId,
            pGauge->m_pDeviceData,
            m_pStaticBosonField->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pStappleCache[pGauge->m_byFieldId],
            appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
            appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
            (NULL == pStaple) ? NULL : pStaple->m_pDeviceData,
            pGaugeForce->m_pDeviceData,
            //NULL,
            m_fBetaOverN
            );
#else
        _LAUNCH_KERNEL(_kernelStapleAtSiteSU3CacheIndex_TemperatureDistribution, block, threads,
            pGauge->m_byFieldId,
            m_pStaticBosonField->m_byFieldId,
            pGauge->m_pDeviceData,
            m_pStaticBosonField->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pStappleCache[pGauge->m_byFieldId],
            (NULL == pStaple) ? NULL : pStaple->m_pDeviceData,
            pGaugeForce->m_pDeviceData,
            //NULL,
            m_fBetaOverN
        );
#endif
        return TRUE;
    }

    if (m_byGaugeFieldIds.Num() < 1 || m_byBosonFieldIds.Num() < 1)
    {
        appCrucial(_T("CActionTemperatureDistribution needs a gauge field and a boson field!\n"));
        return FALSE;
    }

    INT iGaugeidx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
    INT iBosonidx = CLatticeData::GetBosonFieldIndexById(bosonNum, bosonFields, m_byBosonFieldIds[0]);
    if (iGaugeidx < 0 || iBosonidx < 0)
    {
        appCrucial(_T("CActionTemperatureDistribution needs a gauge field and a boson field!\n"));
        return FALSE;
    }

    const CFieldGaugeSU3* pGauge = dynamic_cast<const CFieldGaugeSU3*>(gaugeFields[iGaugeidx]);
    const CFieldBosonReal* pBoson = dynamic_cast<const CFieldBosonReal*>(bosonFields[iBosonidx]);
    if (NULL == pGauge || NULL == pBoson)
    {
        appCrucial(_T("CActionTemperatureDistribution needs a SU3 gauge field and a Real boson field!\n"));
        return FALSE;
    }

    CFieldGaugeSU3* pGaugeForce = dynamic_cast<CFieldGaugeSU3*>(gaugeForces[iGaugeidx]);
    CFieldGaugeSU3* pStaple = (NULL != stapleFields) ? dynamic_cast<CFieldGaugeSU3*>(stapleFields[iGaugeidx]) : NULL;
    CFieldBosonReal* pBosonForce = dynamic_cast<CFieldBosonReal*>(bosonForces[iBosonidx]);
    if (NULL == pGaugeForce || NULL == pBosonForce)
    {
        appCrucial(_T("CActionTemperatureDistribution both field should be dynamic\n"));
        return FALSE;
    }

    preparethread;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelStapleAtSiteSU3CacheIndex_TemperatureDistribution, block, threads, 
        pGauge->m_byFieldId,
        pBoson->m_byFieldId,
        pGauge->m_pDeviceData,
        pBoson->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[pGauge->m_byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        (NULL == pStaple) ? NULL : pStaple->m_pDeviceData,
        pGaugeForce->m_pDeviceData,
        //pBosonForce->m_pDeviceData,
        m_fBetaOverN
        );
#else
    _LAUNCH_KERNEL(_kernelStapleAtSiteSU3CacheIndex_TemperatureDistribution, block, threads,
        pGauge->m_byFieldId,
        pBoson->m_byFieldId,
        pGauge->m_pDeviceData,
        pBoson->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[pGauge->m_byFieldId],
        (NULL == pStaple) ? NULL : pStaple->m_pDeviceData,
        pGaugeForce->m_pDeviceData,
        //pBosonForce->m_pDeviceData,
        m_fBetaOverN
    );
#endif
    _LAUNCH_KERNEL(_kernelPlaqutteEnergySU3_BosonForceUseClover_TemperatureDistribution, block, threads, 
        pGauge->m_byFieldId,
        pGauge->m_pDeviceData,
        pBoson->m_pDeviceData,
        m_fBetaOverN,
        pBosonForce->m_pDeviceData);

    //add kinetic term force
    if (abs(m_fBosonKinetic) > _CLG_FLT_MIN_)
    {
        //add to boson force
        CFieldBoson* dboson = dynamic_cast<CFieldBoson*>(appGetLattice()->GetPooledCopy(pBoson, _T(__FILE__), __LINE__));
        dboson->D(0, 0, NULL, NULL);
        /**
        dboson->D(gaugeNum, bosonNum, gaugeFields, bosonFields);
        checkCudaErrors(cudaDeviceSynchronize());

        bosonForces[ibosonidx]->Axpy(_make_cuComplex(F(1.0), F(0.0)), dboson);
        checkCudaErrors(cudaDeviceSynchronize());
        bosonForces[ibosonidx]->Axpy(_make_cuComplex(-F(1.0) * (F(8.0) + m_fM), F(0.0)), bosonfield);
        checkCudaErrors(cudaDeviceSynchronize());
        bosonfield->CopyTo(dboson);
        dboson->Mul(bosonfield);
        checkCudaErrors(cudaDeviceSynchronize());
        dboson->Mul(bosonfield, FALSE);
        checkCudaErrors(cudaDeviceSynchronize());
        bosonForces[ibosonidx]->Axpy(_make_cuComplex(-F(2.0) * m_fLambda, F(0.0)), dboson);
        checkCudaErrors(cudaDeviceSynchronize());
        */


        pBosonForce->Axpy(static_cast<Real>(m_fBosonKinetic), dboson);
        pBosonForce->Axpy(-static_cast<Real>(m_fBosonKinetic * 8.0), pBoson);
        //pBoson->DebugPrintMe();
        //appGeneral(_T("===================="));
        //dboson->DebugPrintMe();
        //appGeneral(_T("===================="));
        //dboson->Axpy(-F(8.0), pBoson);
        //dboson->DebugPrintMe();

        //appGeneral(_T("boson kinetic = %f\n"), m_fBosonKinetic);

        dboson->Return();
    }

    return TRUE;
}

void CActionTemperatureDistribution::PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate)
{
    if (0 == iUpdateIterate)
    {
        m_fLastEnergy = Energy(FALSE, gaugeNum, bosonNum, gaugeFields, bosonFields, NULL);
    }
}


CCString CActionTemperatureDistribution::GetInfos(const CCString& tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    if (NULL == m_pStaticBosonField)
    {
        sRet = sRet + tab + _T("With dynamic boson field\n");
    }
    else
    {
        sRet = sRet + tab + _T("Static Boson FieldId : ") + appToString(m_pStaticBosonField->m_byFieldId) + _T("\n");
    }
    return sRet;
}


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================