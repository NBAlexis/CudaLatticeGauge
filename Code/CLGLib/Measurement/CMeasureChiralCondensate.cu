//=============================================================================
// FILENAME : CMeasureChiralCondensate.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [06/13/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureChiralCondensate)

#pragma region kernels

__global__ void _CLG_LAUNCH_BOUND
_kernelDotAndGatherXY(
    const deviceWilsonVectorSU3 * __restrict__ pMe,
    const deviceWilsonVectorSU3 * __restrict__ pOther,
    CLGComplex * resultXYPlane,
    CLGComplex * result)
{
    intokernal;

    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
    result[uiSiteIndex] = pMe[uiSiteIndex].ConjugateDotC(pOther[uiSiteIndex]);
    atomicAdd(&resultXYPlane[uiXY].x, result[uiSiteIndex].x);
    atomicAdd(&resultXYPlane[uiXY].y, result[uiSiteIndex].y);
}

#pragma endregion

CMeasureChiralCondensate::~CMeasureChiralCondensate()
{
    if (NULL != m_pDeviceXYBuffer)
    {
        checkCudaErrors(cudaFree(m_pDeviceXYBuffer));
        free(m_pHostXYBuffer);
    }
}

void CMeasureChiralCondensate::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceXYBuffer, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    m_pHostXYBuffer = (CLGComplex*)malloc(sizeof(CLGComplex) * _HC_Lx * _HC_Ly);

    Reset();

    INT iValue = 2;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("FieldCount"), iValue);
    m_uiFieldCount = static_cast<UINT>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;
}

void CMeasureChiralCondensate::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{
    CLGComplex res = _make_cuComplex(F(0.0), F(0.0));

    CFieldFermion* pF1 = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    CFieldFermion* pF2 = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    UINT uiVolume = appGetLattice()->m_pIndexCache->m_uiSiteNumber[m_byFieldId];
    CFieldFermionWilsonSquareSU3 * pF1W = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pF1);
    CFieldFermionWilsonSquareSU3 * pF2W = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pF2);

    _ZeroXYPlaneC(m_pDeviceXYBuffer);
    for (UINT i = 0; i < m_uiFieldCount; ++i)
    {
        pF1->InitialField(EFIT_RandomZ4);
        pF1->FixBoundary();
        //pF1->DebugPrintMe();
        //pGauge->DebugPrintMe();
        pF1->CopyTo(pF2);
        pF1->InverseD(pGauge);       

#pragma region Dot

        preparethread;
        _kernelDotAndGatherXY << <block, threads >> > (
            pF1W->m_pDeviceData,
            pF2W->m_pDeviceData,
            m_pDeviceXYBuffer,
            _D_ComplexThreadBuffer);

        CLGComplex thisSum = appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
#pragma endregion

        res.x = res.x + thisSum.x / uiVolume;
        res.y = res.y + thisSum.y / uiVolume;

    }

    checkCudaErrors(cudaMemcpy(m_pHostXYBuffer, m_pDeviceXYBuffer, sizeof(CLGComplex) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
    if (m_bShowResult)
    {
        appDetailed(_T("\n ------ Densisty -----\n"));
    }
    
    for (UINT i = static_cast<UINT>(CCommonData::m_sCenter.x); i < _HC_Lx; ++i)
    {
        CLGComplex cvalue = m_pHostXYBuffer[i * _HC_Ly + CCommonData::m_sCenter.y];
        cvalue.x = cvalue.x / (m_uiFieldCount * _HC_Lz * _HC_Lt);
        cvalue.y = cvalue.y / (m_uiFieldCount * _HC_Lz * _HC_Lt);
        m_lstCondensateDensity.AddItem(cvalue);
        if (m_bShowResult)
        {
            appDetailed(_T("(%d,%d)=%1.6f %s %1.6f I   "), i, CCommonData::m_sCenter.y,
                cvalue.x,
                cvalue.y < F(0.0) ? _T("") : _T("+"),
                appAbs(cvalue.y));
        }
    }
    if (m_bShowResult)
    {
        appDetailed(_T("\n ------ Densisty -----\n"));
    }

    pF1->Return();
    pF2->Return();
    res.x = res.x / m_uiFieldCount;
    res.y = res.y / m_uiFieldCount;
    appDetailed(_T("\nChiral Condensate = %2.12f + %2.12f\n"), res.x, res.y);
    ++m_uiConfigurationCount;
    m_lstCondensate.AddItem(res);
}

void CMeasureChiralCondensate::Average(UINT )
{
    //nothing to do
}

void CMeasureChiralCondensate::Report()
{
    assert(m_uiConfigurationCount == static_cast<UINT>(m_lstCondensate.Num()));
    assert(static_cast<UINT>(m_uiConfigurationCount * CCommonData::m_sCenter.x)
        == static_cast<UINT>(m_lstCondensateDensity.Num()));

    appSetLogDate(FALSE);
    CLGComplex tmpChargeSum = _make_cuComplex(F(0.0), F(0.0));
    m_lstAverageCondensateDensity.RemoveAll();

    appGeneral(_T("\n==========================================================================\n"));
    appGeneral(_T("==================== Chiral Condensate (%d con)============================\n"), m_uiConfigurationCount);

    if (m_uiConfigurationCount > 1)
    {
        appGeneral(_T("\n ----------- each configuration ------------- \n"));
        appGeneral(_T("{"));

        for (UINT i = 0; i < m_uiConfigurationCount; ++i)
        {
            tmpChargeSum.x += m_lstCondensate[i].x;
            tmpChargeSum.y += m_lstCondensate[i].y;
            LogGeneralComplex(m_lstCondensate[i]);
        }
        appGeneral(_T("}\n"));

        tmpChargeSum.x = tmpChargeSum.x / m_uiConfigurationCount;
        tmpChargeSum.y = tmpChargeSum.y / m_uiConfigurationCount;
        appGeneral(_T("\n ----------- average condensate = %2.12f + %2.12f ------------- \n"),
            tmpChargeSum.x, tmpChargeSum.y);

        m_cAverageCondensate = tmpChargeSum;
    }
    else
    {
        appGeneral(_T("\n ----------- average condensate = %2.12f + %2.12f ------------- \n"),
            m_lstCondensate[0].x,
            m_lstCondensate[0].y);

        m_cAverageCondensate = m_lstCondensate[0];
    }

    appGeneral(_T("\n ----------- condensate density------------- \n"));
    appGeneral(_T("{\n"));
    for (UINT k = 0; k < m_uiConfigurationCount; ++k)
    {
        appGeneral(_T("{"));
        for (UINT i = 0; i < static_cast<UINT>(CCommonData::m_sCenter.x); ++i)
        {
            LogGeneralComplex(m_lstCondensateDensity[k * CCommonData::m_sCenter.x + i]);

            if (0 == k)
            {
                m_lstAverageCondensateDensity.AddItem(m_lstCondensateDensity[k * CCommonData::m_sCenter.x + i]);
            }
            else
            {
                m_lstAverageCondensateDensity[i] = _cuCaddf(m_lstAverageCondensateDensity[i], m_lstCondensateDensity[k * CCommonData::m_sCenter.x + i]);
            }

            if (k == m_uiConfigurationCount - 1)
            {
                m_lstAverageCondensateDensity[i].x = m_lstAverageCondensateDensity[i].x / m_uiConfigurationCount;
                m_lstAverageCondensateDensity[i].y = m_lstAverageCondensateDensity[i].y / m_uiConfigurationCount;
            }
        }
        appGeneral(_T("},\n"));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("==========================================================================\n"));
    appSetLogDate(TRUE);
}

void CMeasureChiralCondensate::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstCondensate.RemoveAll();
    m_lstCondensateDensity.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================