//=============================================================================
// FILENAME : CMeasurePlaqutteEnergy.cpp
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [02/22/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernels

__global__ void
_CLG_LAUNCH_BOUND
_kernelCalculateCorrelatorSU3(
    const deviceWilsonVectorSU3** __restrict__ sources,
    EGammaMatrix eGamma,
    UBOOL bDirac,
    Real* correlatorRes
)
{
    intokernal;

    /*
    * No comments.... see the detail.pdf for detail
    */

    Real res = F(0.0);
    //Constract U_{ij}
    //Constract U_{p(i)p(j)}
    _Complex uijdagger[3];
    _Complex upipj[3];
    gammaMatrix gm = bDirac ? gammaMatrix(__diracGamma[eGamma]) : gammaMatrix(__chiralGamma[eGamma]);
    for (BYTE i = 1; i < 4; ++i)
    {
        for (BYTE j = i + 1; j < 4; ++j)
        {
            Real tmpRes = F(0.0);

            //UijT 0, 1, 2 and Upipj 0, 3, 6
            memcpy(uijdagger, sources[j * 3 + 0][uiSiteIndex].m_d[i].m_ve, sizeof(_Complex) * 3);
            memcpy(upipj, sources[gm.m_uiIndex[j] * 3 + 0][uiSiteIndex].m_d[gm.m_uiIndex[i]].m_ve, sizeof(_Complex) * 3);

            //uijdagger^* x upipj
            tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
            tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
            tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;

            //UijT 3, 4, 5 and Upipj 1, 4, 7
            memcpy(uijdagger, sources[j * 3 + 1][uiSiteIndex].m_d[i].m_ve, sizeof(_Complex) * 3);
            memcpy(upipj, sources[gm.m_uiIndex[j] * 3 + 1][uiSiteIndex].m_d[gm.m_uiIndex[i]].m_ve, sizeof(_Complex) * 3);
            tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
            tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
            tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
            memcpy(uijdagger, sources[j * 3 + 2][uiSiteIndex].m_d[i].m_ve, sizeof(_Complex) * 3);
            memcpy(upipj, sources[gm.m_uiIndex[j] * 3 + 2][uiSiteIndex].m_d[gm.m_uiIndex[i]].m_ve, sizeof(_Complex) * 3);
            tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
            tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
            tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;

            //coefficient c1c2
            //it can only be 0,2,4,6,8...
            //for 2,6,10,... it is -1, otherwise, it is 1.
            //for 2,6,10, gm.m_byZ4[i] + gm.m_byZ4[j]) & 0x02 = 2, otherwise is 0
            res += (F(1.0) - ((gm.m_byZ4[i] + gm.m_byZ4[j]) & 0x02)) * tmpRes;
        }
    }

    //diagnal
    if (0 == gm.m_byNextSymmetricIndex)
    {
        
        for (BYTE i = 0; i < 4; ++i)
        {
            Real tmpRes = F(0.0);
            //p(i)=i
            memcpy(uijdagger, sources[i * 3 + 0][uiSiteIndex].m_d[i].m_ve, sizeof(_Complex) * 3);
            tmpRes += uijdagger[0].x * uijdagger[0].x + uijdagger[0].y * uijdagger[0].y;
            tmpRes += uijdagger[1].x * uijdagger[1].x + uijdagger[1].y * uijdagger[1].y;
            tmpRes += uijdagger[2].x * uijdagger[2].x + uijdagger[2].y * uijdagger[2].y;
            memcpy(uijdagger, sources[i * 3 + 1][uiSiteIndex].m_d[i].m_ve, sizeof(_Complex) * 3);
            tmpRes += uijdagger[0].x * uijdagger[0].x + uijdagger[0].y * uijdagger[0].y;
            tmpRes += uijdagger[1].x * uijdagger[1].x + uijdagger[1].y * uijdagger[1].y;
            tmpRes += uijdagger[2].x * uijdagger[2].x + uijdagger[2].y * uijdagger[2].y;
            memcpy(uijdagger, sources[i * 3 + 2][uiSiteIndex].m_d[i].m_ve, sizeof(_Complex) * 3);
            tmpRes += uijdagger[0].x * uijdagger[0].x + uijdagger[0].y * uijdagger[0].y;
            tmpRes += uijdagger[1].x * uijdagger[1].x + uijdagger[1].y * uijdagger[1].y;
            tmpRes += uijdagger[2].x * uijdagger[2].x + uijdagger[2].y * uijdagger[2].y;

            res += (F(1.0) - ((gm.m_byZ4[i] + gm.m_byZ4[i]) & 0x02)) * F(0.5) *tmpRes;
        }
    }
    else
    {
        Real tmpRes = F(0.0);

        //U_00 and U_p0p0
        memcpy(uijdagger, sources[0 * 3 + 0][uiSiteIndex].m_d[0].m_ve, sizeof(_Complex) * 3);
        memcpy(upipj, sources[gm.m_uiIndex[0] * 3 + 0][uiSiteIndex].m_d[gm.m_uiIndex[0]].m_ve, sizeof(_Complex) * 3);
        tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
        tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
        tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
        memcpy(uijdagger, sources[0 * 3 + 1][uiSiteIndex].m_d[0].m_ve, sizeof(_Complex) * 3);
        memcpy(upipj, sources[gm.m_uiIndex[0] * 3 + 1][uiSiteIndex].m_d[gm.m_uiIndex[0]].m_ve, sizeof(_Complex) * 3);
        tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
        tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
        tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
        memcpy(uijdagger, sources[0 * 3 + 2][uiSiteIndex].m_d[0].m_ve, sizeof(_Complex) * 3);
        memcpy(upipj, sources[gm.m_uiIndex[0] * 3 + 2][uiSiteIndex].m_d[gm.m_uiIndex[0]].m_ve, sizeof(_Complex) * 3);
        tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
        tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
        tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
        res += (F(1.0) - ((gm.m_byZ4[0] + gm.m_byZ4[gm.m_uiIndex[0]]) & 0x02)) * tmpRes;

        //U_kk and U_pkpk
        BYTE byNextIdx = gm.m_byNextSymmetricIndex;
        tmpRes = F(0.0);
        memcpy(uijdagger, sources[byNextIdx * 3 + 0][uiSiteIndex].m_d[byNextIdx].m_ve, sizeof(_Complex) * 3);
        memcpy(upipj, sources[gm.m_uiIndex[byNextIdx] * 3 + 0][uiSiteIndex].m_d[gm.m_uiIndex[byNextIdx]].m_ve, sizeof(_Complex) * 3);
        tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
        tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
        tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
        memcpy(uijdagger, sources[byNextIdx * 3 + 1][uiSiteIndex].m_d[byNextIdx].m_ve, sizeof(_Complex) * 3);
        memcpy(upipj, sources[gm.m_uiIndex[byNextIdx] * 3 + 1][uiSiteIndex].m_d[gm.m_uiIndex[byNextIdx]].m_ve, sizeof(_Complex) * 3);
        tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
        tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
        tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
        memcpy(uijdagger, sources[byNextIdx * 3 + 2][uiSiteIndex].m_d[byNextIdx].m_ve, sizeof(_Complex) * 3);
        memcpy(upipj, sources[gm.m_uiIndex[byNextIdx] * 3 + 2][uiSiteIndex].m_d[gm.m_uiIndex[byNextIdx]].m_ve, sizeof(_Complex) * 3);
        tmpRes += uijdagger[0].x * upipj[0].x + uijdagger[0].y * upipj[0].y;
        tmpRes += uijdagger[1].x * upipj[1].x + uijdagger[1].y * upipj[1].y;
        tmpRes += uijdagger[2].x * upipj[2].x + uijdagger[2].y * upipj[2].y;
        res += (F(1.0) - ((gm.m_byZ4[byNextIdx] + gm.m_byZ4[gm.m_uiIndex[byNextIdx]]) & 0x02)) * tmpRes;
    }

    //site index = xyz * lt + t
    //change it to t * v_xyz + xyz
    UINT uiXYZ = uiSiteIndex / _DC_Lt;
    UINT uiT = uiSiteIndex % _DC_Lt;
    correlatorRes[uiT * _DC_Volumn_xyz + uiXYZ] = res;
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CMeasureMesonCorrelator)

void CMeasureMesonCorrelator::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    m_pOwner = pOwner;
    m_pLatticeData = pLatticeData;
    m_byId = byId;

    //find the fermion field id and gamma
    INT iValue = 0;
    param.FetchValueINT(_T("FieldId"), iValue);
    if (iValue < 1 || iValue >= kMaxFieldCount)
    {
        appCrucial(_T("CMeasureMesonCorrelator: must give a fermion field id, but fetch %d"), iValue);
    }
    m_byFermionFieldId = static_cast<BYTE>(iValue);
    CCString sValue = _T("UNITY");
    param.FetchStringValue(_T("GammaMatrix"), sValue);
    m_eGamma = __STRING_TO_ENUM(EGammaMatrix, sValue);

    m_uiLt = _HC_Lt;
    m_f2OverVolumnSqrt = F(2.0) / _hostsqrt(static_cast<Real>(_HC_Volumn));
    m_uiResoultCount = 0;
}

void CMeasureMesonCorrelator::OnConfigurationAccepted()
{
    CalculateCorrelator(appGetLattice()->m_pGaugeField);
}

void CMeasureMesonCorrelator::Average(UINT)
{
    //do nothing
}

void CMeasureMesonCorrelator::Report()
{
    if (0 == m_uiResoultCount)
    {
        appGeneral(_T("Not measured yet.\n"));
        return;
    }
    assert(m_lstResults.Num() == static_cast<INT>(m_uiLt));
    for (UINT i = 0; i < m_uiLt; ++i)
    {
        appGeneral(_T("CMeasureMesonCorrelator: C(nt=%d)=%f, log(C(nt))=%f\n"), i, m_lstResults[i], _hostlog(appAbs(m_lstResults[i])));
    }
}

void CMeasureMesonCorrelator::Reset()
{
    m_uiLt = _HC_Lt;
    m_f2OverVolumnSqrt = F(2.0) / _hostsqrt(static_cast<Real>(_HC_Volumn));
    m_uiResoultCount = 0;
    m_lstResults.RemoveAll();
}

void CMeasureMesonCorrelator::CalculateCorrelator(const CFieldGauge* pGaugeField)
{
    if (NULL == pGaugeField || EFT_GaugeSU3 != pGaugeField->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    CFieldFermionWilsonSquareSU3* pFermionSources[12];
    for (UINT i = 0; i < 12; ++i)
    {
        pFermionSources[i] = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(m_byFermionFieldId));
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
            SFermionSource sourceData;
            sourceData.m_eSourceType = EFS_Point;
            sourceData.m_sSourcePoint.x = 0;
            sourceData.m_sSourcePoint.y = 0;
            sourceData.m_sSourcePoint.z = 0;
            sourceData.m_sSourcePoint.w = 0;
            sourceData.m_byColorIndex = c;
            sourceData.m_bySpinIndex = s;
            
            pFermionSources[s * 3 + c]->InitialAsSource(sourceData);

            if (NULL != appGetFermionSolver() && !appGetFermionSolver()->IsAbsoluteAccuracy())
            {
                pFermionSources[s * 3 + c]->m_fLength = pFermionSources[s * 3 + c]->Dot(pFermionSources[s * 3 + c]).x;
            }
            pFermionSources[s * 3 + c]->InverseD(pGaugeField);
            pDevicePtr[s * 3 + c] = pFermionSources[s * 3 + c]->m_pDeviceData;
        }
    }

    deviceWilsonVectorSU3** ppDevicePtr;
    checkCudaErrors(cudaMalloc((void**)&ppDevicePtr, sizeof(deviceWilsonVectorSU3*) * 12));
    checkCudaErrors(cudaMemcpy(ppDevicePtr, pDevicePtr, sizeof(deviceWilsonVectorSU3*) * 12, cudaMemcpyHostToDevice));

    preparethread;
    _kernelCalculateCorrelatorSU3 << <block, threads >> > ((const deviceWilsonVectorSU3**)ppDevicePtr, m_eGamma, TRUE, _D_RealThreadBuffer);

    Real* sumSpatial = (Real*)appAlloca(sizeof(Real) * m_uiLt);
    for (UINT i = 0; i < m_uiLt; ++i)
    {
        sumSpatial[i] = m_f2OverVolumnSqrt * CCudaHelper::ReduceReal(_D_RealThreadBuffer + i * _HC_Volumn_xyz, _HC_Volumn_xyz);
        appParanoiac(_T("CMeasureMesonCorrelator: C(nt=%d)=%f, log(C(nt))=%f\n"), i, sumSpatial[i], _hostlog(appAbs(sumSpatial[i])));
    }

    //anverage
    if (m_uiResoultCount == 0)
    {
        for (UINT i = 0; i < m_uiLt; ++i)
        {
            m_lstResults.AddItem(sumSpatial[i]);
        }
        ++m_uiResoultCount;
    }
    else
    {
        assert(m_lstResults.Num() == static_cast<INT>(m_uiLt));
        for (UINT i = 0; i < m_uiLt; ++i)
        {
            m_lstResults[i] = (m_lstResults[i] * m_uiResoultCount + sumSpatial[i]) / (m_uiResoultCount + 1);
        }
        ++m_uiResoultCount;
    }

    for (UINT i = 0; i < 12; ++i)
    {
        pFermionSources[i]->Return();
    }
    checkCudaErrors(cudaFree(ppDevicePtr));
}


__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================