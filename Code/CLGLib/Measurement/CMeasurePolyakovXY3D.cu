//=============================================================================
// FILENAME : CMeasurePolyakovXY3D.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [29/10/2022 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasurePolyakovXY3D)

#pragma region kernles 

__global__ void
_CLG_LAUNCH_BOUND 
_kernelPolyakovLoopOfSite3D_SU3(
    const deviceSU3* __restrict__ pDeviceBuffer,
    CLGComplex* res,
    Real* resAbs)
{
    UINT uiXY = (threadIdx.x + blockIdx.x * blockDim.x);
    deviceSU3 resZ = deviceSU3::makeSU3Zero();

    for (UINT uiZ = 0; uiZ < _DC_Lz; ++uiZ)
    {
        const UINT uiSiteIndex = uiXY * _DC_GridDimZT + uiZ * _DC_Lt;
        const UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, _DC_Dir - 1);
        const SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
        const UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

        if (0 == uiZ)
        {
            if (__idx->_deviceIsBondOnSurface(uiBigIdx, _DC_Dir - 1))
            {
                resZ = deviceSU3::makeSU3Zero();
            }
            else
            {
                resZ = pDeviceBuffer[uiLinkIdx];
            }
        }
        else
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, _DC_Dir - 1))
            {
                resZ.Mul(pDeviceBuffer[uiLinkIdx]);
            }
        }
    }
    res[uiXY] = resZ.Tr();
    resAbs[uiXY] = _cuCabsf(res[uiXY]);
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelPolyakovLoopOfSite3D_U1(
    const CLGComplex* __restrict__ pDeviceBuffer,
    CLGComplex* res,
    Real* resAbs)
{
    UINT uiXY = (threadIdx.x + blockIdx.x * blockDim.x);

    for (UINT uiZ = 0; uiZ < _DC_Lz; ++uiZ)
    {
        const UINT uiSiteIndex = uiXY * _DC_GridDimZT + uiZ * _DC_Lt;
        const UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, _DC_Dir - 1);
        const SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
        const UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

        if (0 == uiZ)
        {
            if (__idx->_deviceIsBondOnSurface(uiBigIdx, _DC_Dir - 1))
            {
                res[uiXY] = _zeroc;
            }
            else
            {
                res[uiXY] = pDeviceBuffer[uiLinkIdx];
            }
        }
        else
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, _DC_Dir - 1))
            {
                res[uiXY] = _cuCmulf(res[uiXY], pDeviceBuffer[uiLinkIdx]);
            }
        }
    }
    resAbs[uiXY] = _cuCabsf(res[uiXY]);
}

#pragma endregion

CMeasurePolyakovXY3D::~CMeasurePolyakovXY3D()
{
    if (NULL != m_pXYHostLoopDensity)
    {
        free(m_pXYHostLoopDensity);
    }

    if (NULL != m_pXYDeviceLoopDensity)
    {
        checkCudaErrors(cudaFree(m_pXYDeviceLoopDensity));
    }

    if (NULL != m_pXYHostLoopDensityAbs)
    {
        free(m_pXYHostLoopDensityAbs);
    }

    if (NULL != m_pXYDeviceLoopDensityAbs)
    {
        checkCudaErrors(cudaFree(m_pXYDeviceLoopDensityAbs));
    }

    if (NULL != m_pDistributionR)
    {
        checkCudaErrors(cudaFree(m_pDistributionR));
    }

    if (NULL != m_pDistributionP)
    {
        checkCudaErrors(cudaFree(m_pDistributionP));
    }

    if (NULL != m_pDistributionPAbs)
    {
        checkCudaErrors(cudaFree(m_pDistributionPAbs));
    }

    if (NULL != m_pHostDistributionR)
    {
        free(m_pHostDistributionR);
    }

    if (NULL != m_pHostDistributionP)
    {
        free(m_pHostDistributionP);
    }

    if (NULL != m_pHostDistributionPAbs)
    {
        free(m_pHostDistributionPAbs);
    }
}

void CMeasurePolyakovXY3D::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    m_pXYHostLoopDensity = (CLGComplex*)malloc(sizeof(CLGComplex) * _HC_Lx * _HC_Ly);
    checkCudaErrors(cudaMalloc((void**)&m_pXYDeviceLoopDensity, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    m_pXYHostLoopDensityAbs = (Real*)malloc(sizeof(Real) * _HC_Lx * _HC_Ly);
    checkCudaErrors(cudaMalloc((void**)&m_pXYDeviceLoopDensityAbs, sizeof(Real) * _HC_Lx * _HC_Ly));
    Reset();

    INT iValue = 1;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;

    iValue = 0;
    param.FetchValueINT(_T("U1"), iValue);
    m_bU1 = iValue != 0;

    m_bMeasureDistribution = TRUE;

    iValue = 0;
    param.FetchValueINT(_T("ShiftCenter"), iValue);
    m_bShiftCenter = iValue != 0;

    //assuming the center is really at center
    SetMaxAndEdge(&m_uiMaxR, &m_uiEdgeR, m_bShiftCenter);
    checkCudaErrors(cudaMalloc((void**)&m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pDistributionP, sizeof(CLGComplex) * (m_uiMaxR + 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pDistributionPAbs, sizeof(Real) * (m_uiMaxR + 1)));

    m_pHostDistributionR = (UINT*)malloc(sizeof(UINT) * (m_uiMaxR + 1));
    m_pHostDistributionP = (CLGComplex*)malloc(sizeof(CLGComplex) * (m_uiMaxR + 1));
    m_pHostDistributionPAbs = (Real*)malloc(sizeof(Real) * (m_uiMaxR + 1));

}

void CMeasurePolyakovXY3D::OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple)
{

    if (!m_bU1)
    {
        if (NULL == pAcceptGauge || EFT_GaugeSU3 != pAcceptGauge->GetFieldType())
        {
            appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
            return;
        }
        const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pAcceptGauge);
        dim3 block1(_HC_DecompX, 1, 1);
        dim3 threads1(_HC_DecompLx, 1, 1);
        _kernelPolyakovLoopOfSite3D_SU3 << <block1, threads1 >> > (pGaugeSU3->m_pDeviceData, m_pXYDeviceLoopDensity, m_pXYDeviceLoopDensityAbs);
    }
    else
    {
        if (NULL == pAcceptGauge || EFT_GaugeU1 != pAcceptGauge->GetFieldType())
        {
            appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge U1!\n"));
            return;
        }
        const CFieldGaugeU1* pGaugeU1 = dynamic_cast<const CFieldGaugeU1*>(pAcceptGauge);
        dim3 block1(_HC_DecompX, 1, 1);
        dim3 threads1(_HC_DecompLx, 1, 1);
        _kernelPolyakovLoopOfSite3D_U1 << <block1, threads1 >> > (pGaugeU1->m_pDeviceData, m_pXYDeviceLoopDensity, m_pXYDeviceLoopDensityAbs);
    }
    
    checkCudaErrors(cudaMemcpy(m_pXYHostLoopDensity, m_pXYDeviceLoopDensity, sizeof(CLGComplex) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pXYHostLoopDensityAbs, m_pXYDeviceLoopDensityAbs, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
    for (UINT i = CCommonData::m_sCenter.x; i < _HC_Lx; ++i)
    {
        m_lstLoopDensity.AddItem(m_pXYHostLoopDensity[i * _HC_Ly + CCommonData::m_sCenter.y]);
    }

    TransformFromXYDataToRDataOnce_C(
        m_bShiftCenter,
        m_pXYDeviceLoopDensity,
        m_pDistributionR,
        m_pDistributionP,
        m_pHostDistributionR,
        m_pHostDistributionP,
        m_uiMaxR,
        m_uiEdgeR,
        TRUE,
        m_byFieldId,
        m_lstP,
        &m_lstLoopInner,
        m_lstLoop,
        m_lstR,
        m_uiConfigurationCount,
        F(1.0) / static_cast<Real>(_HC_Lz)
    );

    TransformFromXYDataToRDataOnce_R(
        m_bShiftCenter,
        m_pXYDeviceLoopDensityAbs,
        m_pDistributionR,
        m_pDistributionPAbs,
        m_pHostDistributionR,
        m_pHostDistributionPAbs,
        m_uiMaxR,
        m_uiEdgeR,
        FALSE,
        m_byFieldId,
        m_lstPAbs,
        &m_lstLoopAbsInner,
        m_lstLoopAbs,
        m_lstR,
        m_uiConfigurationCount,
        F(1.0) / static_cast<Real>(_HC_Lz)
    );

    if (m_bShowResult)
    {
        appDetailed(_T("\n\n ==================== Polyakov Loop (%d con)============================ \n\n"), m_uiConfigurationCount);
    }

    if (m_bShowResult)
    {
        appSetLogDate(FALSE);
        appGeneral(_T("Loop is "));
        LogGeneralComplex(m_lstLoop[m_lstLoop.GetCount() - 1]);
        appGeneral(_T(" Abs is %f\n"), m_lstLoopAbs[m_lstLoopAbs.GetCount() - 1]);
        appSetLogDate(TRUE);
    }

    if (m_bShowResult)
    {
        for (UINT i = 1; i < _HC_Lx; ++i)
        {
            appDetailed(_T("{"));
            for (UINT j = 1; j < _HC_Ly; ++j)
            {
                appDetailed(_T("%1.12f %s %1.12f I%s"),
                    m_pXYHostLoopDensity[i * _HC_Ly + j].x,
                    m_pXYHostLoopDensity[i * _HC_Ly + j].y < F(0.0) ? _T("-") : _T("+"),
                    appAbs(m_pXYHostLoopDensity[i * _HC_Ly + j].y),
                    (j == _HC_Ly - 1) ? _T("},\n") : _T(",   ")
                );
            }
        }
    }

    if (m_bShowResult)
    {
        appGeneral(_T("\n"));
    }

    if (m_bShowResult)
    {
        appDetailed(_T("\n=====================================================\n"), m_uiConfigurationCount);
    }

    ++m_uiConfigurationCount;
}

void CMeasurePolyakovXY3D::Average(UINT )
{
    //nothing to do
}

void CMeasurePolyakovXY3D::Report()
{
    assert(m_uiConfigurationCount == static_cast<UINT>(m_lstLoop.Num()));
    assert(static_cast<UINT>(m_uiConfigurationCount * CCommonData::m_sCenter.x)
        == static_cast<UINT>(m_lstLoopDensity.Num()));

    appSetLogDate(FALSE);
    CLGComplex tmpChargeSum = _make_cuComplex(F(0.0), F(0.0));
    m_lstAverageLoopDensity.RemoveAll();

    appGeneral(_T("\n\n==========================================================================\n"));
    appGeneral(_T("==================== Polyakov Loop (%d con)============================\n"), m_uiConfigurationCount);

    appGeneral(_T("\n ----------- Loop ------------- \n"));

    appGeneral(_T("{"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        tmpChargeSum.x += m_lstLoop[i].x;
        tmpChargeSum.y += m_lstLoop[i].y;
        LogGeneralComplex(m_lstLoop[i]);
    }
    appGeneral(_T("}\n"));

    tmpChargeSum.x = tmpChargeSum.x / m_uiConfigurationCount;
    tmpChargeSum.y = tmpChargeSum.y / m_uiConfigurationCount;
    m_cAverageLoop = tmpChargeSum;
    appGeneral(_T("\n ----------- average Loop |<P>| = %2.12f arg(P) = %2.12f ------------- \n"), _cuCabsf(tmpChargeSum), __cuCargf(tmpChargeSum));

    appGeneral(_T("\n ----------- Loop density ------------- \n"));

    appGeneral(_T("{\n"));
    for (UINT k = 0; k < m_uiConfigurationCount; ++k)
    {
        appGeneral(_T("{"));
        for (UINT i = 0; i < static_cast<UINT>(CCommonData::m_sCenter.x); ++i)
        {
            LogGeneralComplex(m_lstLoopDensity[k * CCommonData::m_sCenter.x + i]);

            if (0 == k)
            {
                m_lstAverageLoopDensity.AddItem(m_lstLoopDensity[k * CCommonData::m_sCenter.x + i]);
            }
            else
            {
                m_lstAverageLoopDensity[i] = _cuCaddf(m_lstAverageLoopDensity[i], m_lstLoopDensity[k * CCommonData::m_sCenter.x + i]);
            }

            if (k == m_uiConfigurationCount - 1)
            {
                m_lstAverageLoopDensity[i].x = m_lstAverageLoopDensity[i].x / m_uiConfigurationCount;
                m_lstAverageLoopDensity[i].y = m_lstAverageLoopDensity[i].y / m_uiConfigurationCount;
            }
        }
        appGeneral(_T("}\n"));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n==========================================================================\n"));
    appGeneral(_T("==========================================================================\n\n"));
    appSetLogDate(TRUE);
}

void CMeasurePolyakovXY3D::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstLoop.RemoveAll();
    m_lstLoopInner.RemoveAll();
    m_lstLoopAbs.RemoveAll();
    m_lstLoopAbsInner.RemoveAll();
    m_lstLoopDensity.RemoveAll();

    m_lstR.RemoveAll();
    m_lstP.RemoveAll();
    m_lstPAbs.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================