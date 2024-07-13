//=============================================================================
// FILENAME : CMeasureMesonCorrelator.cu
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [02/22/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3.h"
#include "CMeasureMesonCorrelator.h"

//c2^*: if c2=i^0 or i^2, c2*=i^0 or i^2, if c2=i^1, c2*=i^3, if c2=i^3, c2*=i^1
//so, (((c2 & 1) << 1) + c2) & 3
//for c2 = 0, 2, (c2 & 1 << 1)=0, so it is (c2 & 3) = 0, 2
//for c2 = 1, 3, (c2 & 1 << 1)=2, so it is (3 or 5)&3 = 3 or 1
//& 3 is not need because i^1 = i^5
#define __z4h(v) ((((v) & 1) << 1) + (v))

__BEGIN_NAMESPACE

#pragma region kernels

__global__ void
_CLG_LAUNCH_BOUND
_kernelCalculateCorrelatorSU3(
    const deviceWilsonVectorSU3** __restrict__ sources,
    EGammaMatrix eGamma,
#if !_CLG_DOUBLEFLOAT
    DOUBLE* correlatorRes
#else
    Real* correlatorRes
#endif
)
{
    intokernal;

    /*
    * No comments.... see the detail.pdf for detail
    */

#if !_CLG_DOUBLEFLOAT
    DOUBLE res = 0.0;
#else
    Real res = F(0.0);
#endif
    //Constract U_{ij}
    //Constract U_{p(i)p(j)}
    CLGComplex uijdagger[3];
    CLGComplex upipj[3];
    gammaMatrix gm = __chiralGamma[eGamma];
    for (BYTE i = 1; i < 4; ++i)
    {
        for (BYTE j = i + 1; j < 4; ++j)
        {
            Real tmpRes = F(0.0);

            //UijT 0, 1, 2 and Upipj 0, 3, 6
            memcpy(uijdagger, sources[j * 3 + 0][uiSiteIndex].m_d[i].m_ve, sizeof(CLGComplex) * 3);
            memcpy(upipj, sources[gm.m_uiIndex[j] * 3 + 0][uiSiteIndex].m_d[gm.m_uiIndex[i]].m_ve, sizeof(CLGComplex) * 3);

            //uijdagger^* x upipj
            tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
            tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
            tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;

            //UijT 3, 4, 5 and Upipj 1, 4, 7
            memcpy(uijdagger, sources[j * 3 + 1][uiSiteIndex].m_d[i].m_ve, sizeof(CLGComplex) * 3);
            memcpy(upipj, sources[gm.m_uiIndex[j] * 3 + 1][uiSiteIndex].m_d[gm.m_uiIndex[i]].m_ve, sizeof(CLGComplex) * 3);
            tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
            tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
            tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
            memcpy(uijdagger, sources[j * 3 + 2][uiSiteIndex].m_d[i].m_ve, sizeof(CLGComplex) * 3);
            memcpy(upipj, sources[gm.m_uiIndex[j] * 3 + 2][uiSiteIndex].m_d[gm.m_uiIndex[i]].m_ve, sizeof(CLGComplex) * 3);
            tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
            tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
            tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;


            //coefficient c1c2*=i^(c1+c2*)
            //it can only be 0,2,4,6,8...
            //for 2,6,10,... it is -1, otherwise, it is 1.
            //for 2,6,10, gm.m_byZ4[i] + gm.m_byZ4[j]) & 0x02 = 2, otherwise is 0
            //res += (F(1.0) - ((gm.m_byZ4[i] + __z4h(gm.m_byZ4[gm.m_uiIndex[j]])) & 0x02)) * tmpRes;
            res += (F(1.0) - ((gm.m_byZ4[i] + __z4h(gm.m_byZ4[j])) & 0x02)) * tmpRes;
        }
    }

    //diagnal
    if (0 == gm.m_byNextSymmetricIndex)
    {
        
        for (BYTE i = 0; i < 4; ++i)
        {
            Real tmpRes = F(0.0);
            //p(i)=i
            memcpy(uijdagger, sources[i * 3 + 0][uiSiteIndex].m_d[i].m_ve, sizeof(CLGComplex) * 3);
            tmpRes += uijdagger[0].x * uijdagger[0].x + uijdagger[0].y * uijdagger[0].y;
            tmpRes += uijdagger[1].x * uijdagger[1].x + uijdagger[1].y * uijdagger[1].y;
            tmpRes += uijdagger[2].x * uijdagger[2].x + uijdagger[2].y * uijdagger[2].y;
            memcpy(uijdagger, sources[i * 3 + 1][uiSiteIndex].m_d[i].m_ve, sizeof(CLGComplex) * 3);
            tmpRes += uijdagger[0].x * uijdagger[0].x + uijdagger[0].y * uijdagger[0].y;
            tmpRes += uijdagger[1].x * uijdagger[1].x + uijdagger[1].y * uijdagger[1].y;
            tmpRes += uijdagger[2].x * uijdagger[2].x + uijdagger[2].y * uijdagger[2].y;
            memcpy(uijdagger, sources[i * 3 + 2][uiSiteIndex].m_d[i].m_ve, sizeof(CLGComplex) * 3);
            tmpRes += uijdagger[0].x * uijdagger[0].x + uijdagger[0].y * uijdagger[0].y;
            tmpRes += uijdagger[1].x * uijdagger[1].x + uijdagger[1].y * uijdagger[1].y;
            tmpRes += uijdagger[2].x * uijdagger[2].x + uijdagger[2].y * uijdagger[2].y;

            //res += (F(1.0) - ((gm.m_byZ4[i] + __z4h(gm.m_byZ4[i])) & 0x02)) * F(0.5) *tmpRes;
            res += F(0.5) * tmpRes;
        }
    }
    else
    {
        Real tmpRes = F(0.0);

        //U_00 and U_p0p0
        memcpy(uijdagger, sources[0 * 3 + 0][uiSiteIndex].m_d[0].m_ve, sizeof(CLGComplex) * 3);
        memcpy(upipj, sources[gm.m_uiIndex[0] * 3 + 0][uiSiteIndex].m_d[gm.m_uiIndex[0]].m_ve, sizeof(CLGComplex) * 3);
        tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
        tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
        tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
        memcpy(uijdagger, sources[0 * 3 + 1][uiSiteIndex].m_d[0].m_ve, sizeof(CLGComplex) * 3);
        memcpy(upipj, sources[gm.m_uiIndex[0] * 3 + 1][uiSiteIndex].m_d[gm.m_uiIndex[0]].m_ve, sizeof(CLGComplex) * 3);
        tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
        tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
        tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
        memcpy(uijdagger, sources[0 * 3 + 2][uiSiteIndex].m_d[0].m_ve, sizeof(CLGComplex) * 3);
        memcpy(upipj, sources[gm.m_uiIndex[0] * 3 + 2][uiSiteIndex].m_d[gm.m_uiIndex[0]].m_ve, sizeof(CLGComplex) * 3);
        tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
        tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
        tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
        //res += (F(1.0) - ((gm.m_byZ4[0] + __z4h(gm.m_byZ4[gm.m_uiIndex[0]])) & 0x02)) * tmpRes;
        res += tmpRes;

        //U_kk and U_pkpk
        BYTE byNextIdx = gm.m_byNextSymmetricIndex;
        tmpRes = F(0.0);
        memcpy(uijdagger, sources[byNextIdx * 3 + 0][uiSiteIndex].m_d[byNextIdx].m_ve, sizeof(CLGComplex) * 3);
        memcpy(upipj, sources[gm.m_uiIndex[byNextIdx] * 3 + 0][uiSiteIndex].m_d[gm.m_uiIndex[byNextIdx]].m_ve, sizeof(CLGComplex) * 3);
        tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
        tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
        tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
        memcpy(uijdagger, sources[byNextIdx * 3 + 1][uiSiteIndex].m_d[byNextIdx].m_ve, sizeof(CLGComplex) * 3);
        memcpy(upipj, sources[gm.m_uiIndex[byNextIdx] * 3 + 1][uiSiteIndex].m_d[gm.m_uiIndex[byNextIdx]].m_ve, sizeof(CLGComplex) * 3);
        tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
        tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
        tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
        memcpy(uijdagger, sources[byNextIdx * 3 + 2][uiSiteIndex].m_d[byNextIdx].m_ve, sizeof(CLGComplex) * 3);
        memcpy(upipj, sources[gm.m_uiIndex[byNextIdx] * 3 + 2][uiSiteIndex].m_d[gm.m_uiIndex[byNextIdx]].m_ve, sizeof(CLGComplex) * 3);
        tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
        tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
        tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
        //res += (F(1.0) - ((gm.m_byZ4[byNextIdx] + __z4h(gm.m_byZ4[gm.m_uiIndex[byNextIdx]])) & 0x02)) * tmpRes;
        res += tmpRes;
    }

    //site index = xyz * lt + t
    //change it to t * v_xyz + xyz
    UINT uiXYZ = uiSiteIndex / _DC_Lt;
    UINT uiT = uiSiteIndex % _DC_Lt;
    correlatorRes[uiT * _DC_Volume_xyz + uiXYZ] = res;
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CMeasureMesonCorrelator)

void CMeasureMesonCorrelator::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    //find the fermion field id and gamma
    TArray<CCString> allGammas;
    param.FetchStringVectorValue(_T("GammaMatrix"), allGammas);
    if (0 == allGammas.Num())
    {
        allGammas.AddItem(_T("UNITY"));
    }
    m_lstGammas.RemoveAll();
    for (INT i = 0; i < allGammas.Num(); ++i)
    {
        EGammaMatrix eGamma = __STRING_TO_ENUM(EGammaMatrix, allGammas[i]);
        if (m_lstGammas.FindItemIndex(eGamma) < 0)
        {
            m_lstGammas.AddItem(eGamma);
        }
    }
    if (0 == m_lstGammas.Num())
    {
        m_lstGammas.AddItem(UNITY);
    }

    m_uiLt = _HC_Lt;
#if !_CLG_DOUBLEFLOAT
    m_f2OverVolumnSqrt = 2.0 / _hostsqrtd(static_cast<DOUBLE>(_HC_Volume));
#else
    m_f2OverVolumnSqrt = F(2.0) / _hostsqrt(static_cast<Real>(_HC_Volume));
#endif
    m_uiConfigurationCount = 0;
}

void CMeasureMesonCorrelator::OnConfigurationAccepted(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const CFieldGauge* const* pStaples)
{
    CalculateCorrelator(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson, pStaples);
}

void CMeasureMesonCorrelator::Report()
{
    if (0 == m_uiConfigurationCount)
    {
        appGeneral(_T("Not measured yet.\n"));
        return;
    }
    if (0 == m_lstGammas.Num())
    {
        appGeneral(_T("No Gamma matrix is measured.\n"));
        return;
    }

    appPushLogDate(FALSE);

    assert(m_lstResults.Num() == m_lstGammas.Num());
    appGeneral(_T("CMeasureMesonCorrelator final report: Number of configurations = %d\n"), m_uiConfigurationCount);
    for (INT i = 0; i < m_lstResults.Num(); ++i)
    {
        assert(m_lstResults[i].Num() == static_cast<INT>(m_uiLt));
        appGeneral(_T("CMeasureMesonCorrelator final report Gamma = %s, C(nt=0 to %d)=\n"),
            __ENUM_TO_STRING(EGammaMatrix, m_lstGammas[i]).c_str(),
            m_lstResults[i].Num() - 1);
        for (UINT j = 0; j < m_uiLt; ++j)
        {
            appGeneral(_T("%8.12f, "), m_lstResults[i][j]);
        }
        appGeneral(_T("\n log10(C(nt))=\n"));
        for (UINT j = 0; j < m_uiLt; ++j)
        {
#if !_CLG_DOUBLEFLOAT
            appGeneral(_T("%8.12f, "), _hostlog10d(appAbs(m_lstResults[i][j])));
#else
            appGeneral(_T("%8.12f, "), _hostlog10(appAbs(m_lstResults[i][j])));
#endif
        }
        appGeneral(_T("\n log10(C(nt)/C(0))=\n"));
        for (UINT j = 1; j < m_uiLt; ++j)
        {
#if !_CLG_DOUBLEFLOAT
            appGeneral(_T("%8.12f, "), _hostlog10d(appAbs(m_lstResults[i][j] / m_lstResults[i][0])));
#else
            appGeneral(_T("%8.12f, "), _hostlog10(appAbs(m_lstResults[i][j] / m_lstResults[i][0])));
#endif
        }
        appGeneral(_T("\n"));
    }

    appPopLogDate();
}

void CMeasureMesonCorrelator::Reset()
{
    CMeasure::Reset();
    m_uiLt = _HC_Lt;
    m_f2OverVolumnSqrt = F(2.0) / _hostsqrt(static_cast<Real>(_HC_Volume));
    m_lstResults.RemoveAll();
    m_lstResultsLastConf.RemoveAll();
}

void CMeasureMesonCorrelator::CalculateCorrelator(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const CFieldGauge* const* pStaple)
{
    //CFieldGauge* pCopyGauge = NULL;
    //const CFieldGauge* pGaugeField = NULL;
    TArray<const CFieldGauge*> gauges;
    TArray<CFieldGauge*> returngauges;
    if (m_bNeedSmearing && NULL != appGetGaugeSmearing())
    {
        for (INT i = 0; i < gaugeNum; ++i)
        {
            CFieldGauge* pCopyGauge = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pAcceptGauge[i]));
            returngauges.AddItem(pCopyGauge);
            CFieldGauge* pCopyStaple = NULL;
            if (NULL != pStaple)
            {
                pCopyStaple = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pStaple[i]));
            }
            else
            {
                pCopyStaple = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(pStaple[i]->m_byFieldId));
                pCopyGauge->CalculateOnlyStaple(pCopyStaple);
            }
            appGetGaugeSmearing()->GaugeSmearing(pCopyGauge, pCopyStaple);
            gauges.AddItem(pCopyGauge);
            pCopyStaple->Return();
        }
    }
    else
    {
        for (INT i = 0; i < gaugeNum; ++i)
        {
            gauges.AddItem(pAcceptGauge[i]);
        }
    }

    CFieldFermionWilsonSquareSU3* pFermionSources[12];
    for (UINT i = 0; i < 12; ++i)
    {
        pFermionSources[i] = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(GetFermionFieldId()));
        if (NULL == pFermionSources[i])
        {
            appCrucial(_T("Meson correlator only implemented with Wilson SU3\n"));
            _FAIL_EXIT;
        }
    }

    deviceWilsonVectorSU3* pDevicePtr[12];
    for (BYTE s = 0; s < 4; ++s)
    {
        for (BYTE c = 0; c < 3; ++c)
        {
            SFermionBosonSource sourceData;
            sourceData.m_eSourceType = EFS_Point;
            sourceData.m_sSourcePoint.x = 0;
            sourceData.m_sSourcePoint.y = 0;
            sourceData.m_sSourcePoint.z = 0;
            sourceData.m_sSourcePoint.w = 0;
            sourceData.m_byColorIndex = c;
            sourceData.m_bySpinIndex = s;
            
            pFermionSources[s * 3 + c]->InitialAsSource(sourceData);

            if (NULL != appGetFermionSolver(GetFermionFieldId()) && !appGetFermionSolver(GetFermionFieldId())->IsAbsoluteAccuracy())
            {
                pFermionSources[s * 3 + c]->m_fLength = pFermionSources[s * 3 + c]->Dot(pFermionSources[s * 3 + c]).x;
            }
            pFermionSources[s * 3 + c]->InverseD(gaugeNum, bosonNum, gauges.GetData(), pAcceptBoson);
            pDevicePtr[s * 3 + c] = pFermionSources[s * 3 + c]->m_pDeviceData;
        }
    }

    deviceWilsonVectorSU3** ppDevicePtr;
    checkCudaErrors(cudaMalloc((void**)&ppDevicePtr, sizeof(deviceWilsonVectorSU3*) * 12));
    checkCudaErrors(cudaMemcpy(ppDevicePtr, pDevicePtr, sizeof(deviceWilsonVectorSU3*) * 12, cudaMemcpyHostToDevice));

    preparethread;
    for (INT i = 0; i < m_lstGammas.Num(); ++i)
    {
        _kernelCalculateCorrelatorSU3 << <block, threads >> > (
            (const deviceWilsonVectorSU3**)ppDevicePtr, 
            m_lstGammas[i],
            _D_RealThreadBuffer);

        //if (0 == i && 0 == m_uiResoultCount)
        //{
        //    Real fSum = F(0.0);
        //    Real* hostRes = (Real*)malloc(sizeof(Real) * _HC_Volume_xyz);
        //    checkCudaErrors(cudaMemcpy(hostRes, _D_RealThreadBuffer, sizeof(Real) * _HC_Volume_xyz, cudaMemcpyDeviceToHost));
        //    for (int j = 0; j < _HC_Volume_xyz; ++j)
        //    {
        //        //appGeneral(_T("%f\n"), hostRes[j]);
        //        fSum += hostRes[j];
        //    }
        //    free(hostRes);
        //    appGeneral(_T("sum of t=0 is %f\n"), fSum * m_f2OverVolumnSqrt);
        //}

        appParanoiac(_T("CMeasureMesonCorrelator for %s: "), __ENUM_TO_STRING(EGammaMatrix, m_lstGammas[i]).c_str());
#if !_CLG_DOUBLEFLOAT
        DOUBLE* sumSpatial = (DOUBLE*)appAlloca(sizeof(DOUBLE) * m_uiLt);
#else
        Real* sumSpatial = (Real*)appAlloca(sizeof(Real) * m_uiLt);
#endif
        for (UINT j = 0; j < m_uiLt; ++j)
        {
            sumSpatial[j] = m_f2OverVolumnSqrt * CCudaHelper::ReduceReal(_D_RealThreadBuffer + j * _HC_Volume_xyz, _HC_Volume_xyz);
#if !_CLG_DOUBLEFLOAT
            appParanoiac(_T("C(nt=%d)=%f, log10(C(nt))=%f, \n"), j, sumSpatial[j], _hostlog10d(appAbs(sumSpatial[j])));
#else
            appParanoiac(_T("C(nt=%d)=%f, log10(C(nt))=%f, \n"), j, sumSpatial[j], _hostlog10(appAbs(sumSpatial[j])));
#endif
        }
        appParanoiac(_T("\n"));
        //anverage
        if (m_uiConfigurationCount == 0)
        {
#if !_CLG_DOUBLEFLOAT
            TArray<DOUBLE> thisGammaResult;
#else
            TArray<Real> thisGammaResult;
#endif
            for (UINT j = 0; j < m_uiLt; ++j)
            {
                thisGammaResult.AddItem(sumSpatial[j]);
            }
            m_lstResults.AddItem(thisGammaResult);
            m_lstResultsLastConf.AddItem(thisGammaResult);
        }
        else
        {
            assert(m_lstResults[i].Num() == static_cast<INT>(m_uiLt));
            for (UINT j = 0; j < m_uiLt; ++j)
            {
                m_lstResultsLastConf[i][j] = sumSpatial[j];
                m_lstResults[i][j] = 
                    (m_lstResults[i][j] * m_uiConfigurationCount + sumSpatial[j]) / (m_uiConfigurationCount + 1);
            }
        }
    }
    ++m_uiConfigurationCount;

    for (UINT i = 0; i < 12; ++i)
    {
        pFermionSources[i]->Return();
    }
    checkCudaErrors(cudaFree(ppDevicePtr));

    for (INT i = 0; i < returngauges.Num(); ++i)
    {
        returngauges[i]->Return();
    }
}


__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================