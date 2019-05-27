//=============================================================================
// FILENAME : CMeasureTopologicChargeXLine.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/28/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureTopologicChargeXLine)

#pragma region kernles

__global__ void
_CLG_LAUNCH_(12, 1)
_kernel_TopoChargeXL(
    deviceWilsonVectorSU3** pSources,
    SSmallInt4 sSite4,
    BYTE byArrayIdx,
    BYTE byFieldId,
    deviceWilsonVectorSU3* res)
{
    //s * 3 + c
    UINT uiC = threadIdx.x;
    UINT uiS = threadIdx.y;
    UINT uiCS = uiS * 3 + uiC;

    UINT uiSiteIndex = _deviceGetSiteIndex(sSite4);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        res[byArrayIdx] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    deviceWilsonVectorSU3 right_element(pSources[uiCS][uiSiteIndex]);
    right_element = __chiralGamma[GAMMA5].MulWilsonC(right_element);

    //Note that, it is not res[byArrayIdx] = result
    //It is res[c,s] = delta_{cc}delta_ss result[c,s]
    res[byArrayIdx].m_d[uiS].m_ve[uiC] = right_element.m_d[uiS].m_ve[uiC];
}

/**
*
*/
__global__ void
_CLG_LAUNCH_(128, 1)
_kernel_Trace_TopoChargeXL(
    const deviceWilsonVectorSU3* __restrict__ pOperator,
    CLGComplex* pResLine)
{
    UINT uiIdx = threadIdx.x;
    pResLine[uiIdx] = pOperator[uiIdx].Sum();
}

#pragma endregion

CMeasureTopologicChargeXLine::~CMeasureTopologicChargeXLine()
{
    if (NULL != m_pHostDataBuffer)
    {
        free(m_pHostDataBuffer);
    }
    if (NULL != m_pDeviceDataBuffer)
    {
        checkCudaErrors(cudaFree(m_pDeviceDataBuffer));
    }
    if (NULL != m_pOperatorData)
    {
        checkCudaErrors(cudaFree(m_pOperatorData));
    }

}

void CMeasureTopologicChargeXLine::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    m_pHostDataBuffer = (CLGComplex*)malloc(sizeof(CLGComplex) * (_HC_Lx - 1));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDataBuffer, sizeof(CLGComplex) * (_HC_Lx - 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pOperatorData, sizeof(deviceWilsonVectorSU3) * (_HC_Lx - 1)));

    Reset();

    INT iValue = 1;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;
}

void CMeasureTopologicChargeXLine::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    dim3 _blocks(1, 1, 1);
    dim3 _thread1(12, 1, 1);
    for (UINT i = 1; i < _HC_Lx; ++i)
    {
        CFieldFermionWilsonSquareSU3DR* pFermionSources[12];
        for (UINT j = 0; j < 12; ++j)
        {
            pFermionSources[j] = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
            if (NULL == pFermionSources[j])
            {
                appCrucial(_T("Meson correlator only implemented with Wilson SU3 Dirichlet Rotating\n"));
                _FAIL_EXIT;
            }
        }

        deviceWilsonVectorSU3* pDevicePtr[12];
        SSmallInt4 sourceSite;
        sourceSite.x = static_cast<SBYTE>(i);
        sourceSite.y = CCommonData::m_sCenter.y;
        sourceSite.z = CCommonData::m_sCenter.z;
        sourceSite.w = CCommonData::m_sCenter.w;
        for (BYTE s = 0; s < 4; ++s)
        {
            for (BYTE c = 0; c < 3; ++c)
            {
                SFermionSource sourceData;
                sourceData.m_eSourceType = EFS_Point;
                sourceData.m_sSourcePoint = sourceSite;
                sourceData.m_byColorIndex = c;
                sourceData.m_bySpinIndex = s;

                pFermionSources[s * 3 + c]->InitialAsSource(sourceData);

                if (NULL != appGetFermionSolver() && !appGetFermionSolver()->IsAbsoluteAccuracy())
                {
                    pFermionSources[s * 3 + c]->m_fLength = pFermionSources[s * 3 + c]->Dot(pFermionSources[s * 3 + c]).x;
                }
                //================NOTE===================
                //This is different from the other charges
                pFermionSources[s * 3 + c]->D(pGaugeSU3);
                pDevicePtr[s * 3 + c] = pFermionSources[s * 3 + c]->m_pDeviceData;
            }
        }

        deviceWilsonVectorSU3** ppDevicePtr;
        checkCudaErrors(cudaMalloc((void**)&ppDevicePtr, sizeof(deviceWilsonVectorSU3*) * 12));
        checkCudaErrors(cudaMemcpy(ppDevicePtr, pDevicePtr, sizeof(deviceWilsonVectorSU3*) * 12, cudaMemcpyHostToDevice));

        _kernel_TopoChargeXL << <_blocks, _thread1 >> > (
            ppDevicePtr,
            sourceSite,
            static_cast<BYTE>(i - 1),
            m_byFieldId,
            m_pOperatorData
            );

        checkCudaErrors(cudaFree(ppDevicePtr));
        for (UINT j = 0; j < 12; ++j)
        {
            pFermionSources[j]->Return();
        }
    }

    ++m_uiConfigurationCount;
    dim3 _thread2(_HC_Lx - 1, 1, 1);
    _kernel_Trace_TopoChargeXL << <_blocks, _thread2 >> > (m_pOperatorData, m_pDeviceDataBuffer);

    if (m_bShowResult)
    {
        appDetailed(_T("\n\n ==================== Topological Charge (%d con)============================ \n\n"), m_uiConfigurationCount);
    }
    checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBuffer, sizeof(CLGComplex) * (_HC_Lx - 1), cudaMemcpyDeviceToHost));

    for (UINT i = 0; i < _HC_Lx - 1; ++i)
    {
        m_lstAllRes.AddItem(m_pHostDataBuffer[i].x);
        if (m_bShowResult)
        {
            appDetailed(_T("%d=(%1.6f,%1.6f)   "), i, m_pHostDataBuffer[i].x, m_pHostDataBuffer[i].y);
        }
    }

    if (m_bShowResult)
    {
        appDetailed(_T("\n=====================================================\n"), m_uiConfigurationCount);
    }
}

void CMeasureTopologicChargeXLine::Average(UINT )
{
    //nothing to do
}

void CMeasureTopologicChargeXLine::Report()
{
    assert(m_uiConfigurationCount * (_HC_Lx - 1) == static_cast<UINT>(m_lstAllRes.Num()));
    appSetLogDate(FALSE);
    TArray<Real> tmpSum;

    appGeneral(_T("\n\n==========================================================================\n"));
    appGeneral(_T("==================== Topological Charge (%d con)============================\n"), m_uiConfigurationCount);

    appGeneral(_T("{\n"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        appGeneral(_T("{"));
        for (UINT j = 0; j < _HC_Lx - 1; ++j)
        {
            appGeneral(_T("%2.12f, "), m_lstAllRes[i * (_HC_Lx - 1) + j]);
            if (0 == i)
            {
                tmpSum.AddItem(m_lstAllRes[i* (_HC_Lx - 1) + j]);
            }
            else
            {
                tmpSum[j] += m_lstAllRes[i * (_HC_Lx - 1) + j];
            }
        }
        appGeneral(_T("},\n"));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n ----------- average ------------- \n"));
    for (UINT j = 0; j < _HC_Lx - 1; ++j)
    {
        appGeneral(_T("%2.12f, "), tmpSum[j] / m_uiConfigurationCount);
    }

    appGeneral(_T("\n==========================================================================\n"));
    appGeneral(_T("==========================================================================\n\n"));
    appSetLogDate(TRUE);
}

void CMeasureTopologicChargeXLine::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstAllRes.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================