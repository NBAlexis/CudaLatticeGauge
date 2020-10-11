//=============================================================================
// FILENAME : CMeasureTopologicChargeXY.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/28/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureTopologicChargeXY)

#pragma region kernles

__global__ void
_CLG_LAUNCH_BOUND
_kernelTopoChargeClover(
    const deviceSU3* __restrict__ pDeviceBuffer,
    BYTE byFieldId,
#if !_CLG_DOUBLEFLOAT
    DOUBLE* pResBuffer
#else
    Real* pResBuffer
#endif
)
{
    intokernalInt4;
    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);
#if !_CLG_DOUBLEFLOAT
    DOUBLE fRes = 0.0;
#else
    Real fRes = F(0.0);
#endif
    if (!__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet())
    {
        fRes = _deviceTopologicalCharge(pDeviceBuffer, byFieldId, sSite4, uiN);
    }

    pResBuffer[uiSiteIndex] = fRes;
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelTopoChargeSumOverZT(
#if !_CLG_DOUBLEFLOAT
    const DOUBLE* __restrict__ pChargeDensity,
#else
    const Real* __restrict__ pChargeDensity,   
#endif
    Real* densityXYPlane
)
{
    intokernalInt4;
    const UINT uiXYIndex = sSite4.x * _DC_Ly + sSite4.y;
#if !_CLG_DOUBLEFLOAT
    atomicAdd(&densityXYPlane[uiXYIndex], static_cast<Real>(pChargeDensity[uiSiteIndex]));
#else
    atomicAdd(&densityXYPlane[uiXYIndex], pChargeDensity[uiSiteIndex]);
#endif
}

#pragma endregion

CMeasureTopologicChargeXY::~CMeasureTopologicChargeXY()
{
    if (NULL != m_pXYHostDensity)
    {
        free(m_pXYHostDensity);
    }

    if (NULL != m_pXYDeviceDensity)
    {
        checkCudaErrors(cudaFree(m_pXYDeviceDensity));
    }
}

void CMeasureTopologicChargeXY::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    m_pXYHostDensity = (Real*)malloc(sizeof(Real) * _HC_Lx * _HC_Ly);
    checkCudaErrors(cudaMalloc((void**)&m_pXYDeviceDensity, sizeof(Real) * _HC_Lx * _HC_Ly));

    Reset();

    INT iValue = 1;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;
}

void CMeasureTopologicChargeXY::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    preparethread;
    _kernelTopoChargeClover << <block, threads >>>(pGaugeSU3->m_pDeviceData, pGaugeSU3->m_byFieldId, _D_RealThreadBuffer);

    _ZeroXYPlane(m_pXYDeviceDensity);
    _kernelTopoChargeSumOverZT<<<block, threads>>>(_D_RealThreadBuffer, m_pXYDeviceDensity);

    checkCudaErrors(cudaMemcpy(m_pXYHostDensity, m_pXYDeviceDensity, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));

    const Real fCharge = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    ++m_uiConfigurationCount;
    if (m_bShowResult)
    {
        appDetailed(_T("\n\n ==================== Topological Charge (%d con)============================ \n\n"), m_uiConfigurationCount);
    }
    m_lstCharge.AddItem(fCharge);
    if (m_bShowResult)
    {
        appDetailed(_T("Charge is %f\n"), fCharge);
    }
    for (UINT i = 1; i < _HC_Ly; ++i)
    {
        for (UINT j = 1; j < _HC_Lx; ++j)
        {
            m_lstXYDensity.AddItem(m_pXYHostDensity[j * _HC_Ly + i]);
            if (m_bShowResult)
            {
                appDetailed(_T("(%d,%d)=%1.6f   "), j, i, m_pXYHostDensity[j * _HC_Ly + i]);
            }
        }
        appDetailed(_T("\n"));
    }

    if (m_bShowResult)
    {
        appDetailed(_T("\n=====================================================\n"), m_uiConfigurationCount);
    }
}

void CMeasureTopologicChargeXY::Average(UINT )
{
    //nothing to do
}

void CMeasureTopologicChargeXY::Report()
{
    assert(m_uiConfigurationCount == static_cast<UINT>(m_lstCharge.Num()));
    assert(m_uiConfigurationCount * (_HC_Lx - 1) * (_HC_Ly - 1) == static_cast<UINT>(m_lstXYDensity.Num()));

    appSetLogDate(FALSE);
    Real tmpChargeSum = F(0.0);

    appGeneral(_T("\n\n==========================================================================\n"));
    appGeneral(_T("==================== Topological Charge (%d con)============================\n"), m_uiConfigurationCount);

    appGeneral(_T("\n ----------- charge ------------- \n"));

    appGeneral(_T("{"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        tmpChargeSum += m_lstCharge[i];
        appGeneral(_T("%2.12f,  "), m_lstCharge[i]);
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n ----------- average charge = %2.12f ------------- \n"), tmpChargeSum / m_uiConfigurationCount);

    appGeneral(_T("\n ----------- charge density ------------- \n"));

    TArray<Real> tmp;
    appGeneral(_T("{\n"));
    for (UINT k = 0; k < m_uiConfigurationCount; ++k)
    {
        appGeneral(_T("{"));
        for (UINT i = 0; i < _HC_Ly - 1; ++i)
        {
            appGeneral(_T("{"));
            for (UINT j = 0; j < _HC_Lx - 1; ++j)
            {
                appGeneral(_T("%1.6f, "), m_lstXYDensity[k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j]);
                if (0 == k)
                {
                    tmp.AddItem(m_lstXYDensity[k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j]);
                }
                else
                {
                    tmp[i * (_HC_Lx - 1) + j] += m_lstXYDensity[k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j];
                }
            }
            appGeneral(_T("}, "));
        }
        appGeneral(_T("}\n"));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n ----------- charge density average ------------- \n"));

    for (UINT i = 0; i < _HC_Ly - 1; ++i)
    {
        for (UINT j = 0; j < _HC_Lx - 1; ++j)
        {
            appGeneral(_T("(x=%d,y=%d)%2.8f,   "), j + 1, i + 1, tmp[i * (_HC_Lx - 1) + j] / m_uiConfigurationCount);
        }
        appGeneral(_T("\n"));
    }

    appGeneral(_T("\n==========================================================================\n"));
    appGeneral(_T("==========================================================================\n\n"));
    appSetLogDate(TRUE);
}

void CMeasureTopologicChargeXY::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstCharge.RemoveAll();
    m_lstXYDensity.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================