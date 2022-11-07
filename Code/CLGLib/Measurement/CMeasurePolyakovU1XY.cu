//=============================================================================
// FILENAME : CMeasurePolyakovU1XY.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [29/10/2022 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasurePolyakovU1XY)

#pragma region kernles 

__global__ void
_CLG_LAUNCH_BOUND 
_kernelPolyakovLoopOfSiteU1(
    const CLGComplex* __restrict__ pDeviceBuffer,
    UINT uiT,
    CLGComplex* res)
{
    UINT uiXYZ = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    const UINT uiSiteIndex = uiXYZ * _DC_Lt + uiT;
    UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, _DC_Dir - 1);
    //(uiSiteIndex + 1) * _DC_Dir - 1;//uiSiteIndex * _DC_Dir + (_DC_Dir - 1);
    //if (0 == uiXYZ)
    //{
    //    printf("t=%d, site=%d, linkidx=%d\n", uiT, uiSiteIndex, uiLinkIdx);
    //}

    const SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

    if (0 == uiT)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, _DC_Dir - 1))
        {
            res[uiXYZ] = _zeroc;
        }
        else
        {
            res[uiXYZ] = pDeviceBuffer[uiLinkIdx];
        }
    }
    else
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, _DC_Dir - 1))
        {
            res[uiXYZ] = _zeroc;
        }
        else
        {
            res[uiXYZ] = _cuCmulf(res[uiXYZ], pDeviceBuffer[uiLinkIdx]);
        }
    }
}

/**
 * Before call me, set block dim thread dim.y = 1
 */
__global__ void
_CLG_LAUNCH_BOUND
_kernelPolyakovLoopOfSiteZU1(
    const CLGComplex* __restrict__ pDeviceBuffer,
    CLGComplex* res)
{
    UINT uiXYT = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z);

    for (UINT z = 0; z < _DC_Lz; ++z)
    {
        const UINT uiSiteIndex = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + z * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z);
        UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, 2);
        const SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
        const UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

        if (0 == z)
        {
            if (__idx->_deviceIsBondOnSurface(uiBigIdx, 2))
            {
                res[uiXYT] = _zeroc;
            }
            else
            {
                res[uiXYT] = pDeviceBuffer[uiLinkIdx];
            }
        }
        else
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 2))
            {
                res[uiXYT] = _cuCmulf(res[uiXYT], pDeviceBuffer[uiLinkIdx]);
            }
        }
    }
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelPolyakovTraceOfSiteXYU1(
    const CLGComplex* __restrict__ resXYZ,
    CLGComplex* resXY,
    CLGComplex* resZ,
    Real* resXYAbs,
    Real* resZAbs)
{
    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
    UINT uiZ = threadIdx.y + blockIdx.y * blockDim.y;
    UINT uiXYZ = uiXY * _DC_Lz + uiZ;
    const CLGComplex trres = resXYZ[uiXYZ];
    atomicAdd(&resXY[uiXY].x, trres.x);
    atomicAdd(&resXY[uiXY].y, trres.y);
    atomicAdd(&resXYAbs[uiXY], _cuCabsf(trres));
    if (NULL != resZ)
    {
        atomicAdd(&resZ[uiZ].x, trres.x);
        atomicAdd(&resZ[uiZ].y, trres.y);
    }

    if (NULL != resZAbs)
    {
        atomicAdd(&resZAbs[uiZ], _cuCabsf(trres));
    }
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelInitialZSliceU1(CLGComplex* resZ, Real* resZAbs)
{
    resZ[threadIdx.x + blockIdx.x * blockDim.x] = _zeroc;
    resZAbs[threadIdx.x + blockIdx.x * blockDim.x] = F(0.0);
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelPolyakovZTraceOfSiteXYU1(
    const CLGComplex* __restrict__ resXYT,
    CLGComplex* resXY)
{
    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
    UINT uiXYT = uiXY * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z);
    CLGComplex trres = resXYT[uiXYT];
    atomicAdd(&resXY[uiXY].x, trres.x);
    atomicAdd(&resXY[uiXY].y, trres.y);
}

#pragma endregion

CLGAPI void _PolyakovAtSiteU1(const CLGComplex* __restrict__ pDeviceBuffer, CLGComplex* pRes)
{
    dim3 block1(_HC_DecompX, _HC_DecompY, 1);
    dim3 threads1(_HC_DecompLx, _HC_DecompLy, 1);
    for (UINT uiT = 0; uiT < _HC_Lt; ++uiT)
    {
        _kernelPolyakovLoopOfSiteU1 << <block1, threads1 >> >(pDeviceBuffer, uiT, pRes);
    }
}

CMeasurePolyakovU1XY::~CMeasurePolyakovU1XY()
{
    if (NULL != m_pXYHostLoopDensity)
    {
        free(m_pXYHostLoopDensity);
    }

    if (NULL != m_pXYHostLoopDensityAbs)
    {
        free(m_pXYHostLoopDensityAbs);
    }

    if (NULL != m_pTmpDeviceSum)
    {
        checkCudaErrors(cudaFree(m_pTmpDeviceSum));
    }

    if (NULL != m_pXYDeviceLoopDensity)
    {
        checkCudaErrors(cudaFree(m_pXYDeviceLoopDensity));
    }

    if (NULL != m_pXYDeviceLoopDensityAbs)
    {
        checkCudaErrors(cudaFree(m_pXYDeviceLoopDensityAbs));
    }

    if (NULL != m_pTmpLoop)
    {
        checkCudaErrors(cudaFree(m_pTmpLoop));
    }

    if (NULL != m_pTmpLoopZ)
    {
        checkCudaErrors(cudaFree(m_pTmpLoopZ));
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

    if (NULL != m_pZDeviceLoopDensity)
    {
        checkCudaErrors(cudaFree(m_pZDeviceLoopDensity));
    }

    if (NULL != m_pZHostLoopDensity)
    {
        free(m_pZHostLoopDensity);
    }

    if (NULL != m_pZDeviceLoopDensityAbs)
    {
        checkCudaErrors(cudaFree(m_pZDeviceLoopDensityAbs));
    }

    if (NULL != m_pZHostLoopDensityAbs)
    {
        free(m_pZHostLoopDensityAbs);
    }
}

void CMeasurePolyakovU1XY::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    m_pXYHostLoopDensity = (CLGComplex*)malloc(sizeof(CLGComplex) * _HC_Lx * _HC_Ly);
    m_pXYHostLoopDensityAbs = (Real*)malloc(sizeof(Real) * _HC_Lx * _HC_Ly);
    checkCudaErrors(cudaMalloc((void**)&m_pTmpDeviceSum, sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)&m_pXYDeviceLoopDensity, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    checkCudaErrors(cudaMalloc((void**)&m_pXYDeviceLoopDensityAbs, sizeof(Real) * _HC_Lx * _HC_Ly));
    checkCudaErrors(cudaMalloc((void**)&m_pTmpLoop, sizeof(CLGComplex) * _HC_Lx * _HC_Ly * _HC_Lz));
    Reset();

    INT iValue = 1;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;

    iValue = 0;
    param.FetchValueINT(_T("MeasureZ"), iValue);
    m_bMeasureLoopZ = iValue != 0;

    iValue = 0;
    param.FetchValueINT(_T("ZSlice"), iValue);
    m_bMeasureZSlice = iValue != 0;
    if (m_bMeasureZSlice)
    {
        checkCudaErrors(cudaMalloc((void**)&m_pZDeviceLoopDensity, sizeof(CLGComplex) * _HC_Lz));
        m_pZHostLoopDensity = (CLGComplex*)malloc(sizeof(CLGComplex) * _HC_Lz);
        checkCudaErrors(cudaMalloc((void**)&m_pZDeviceLoopDensityAbs, sizeof(Real) * _HC_Lz));
        m_pZHostLoopDensityAbs = (Real*)malloc(sizeof(Real) * _HC_Lz);
    }

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

    if (m_bMeasureLoopZ)
    {
        checkCudaErrors(cudaMalloc((void**)&m_pTmpLoopZ, sizeof(CLGComplex) * _HC_Lx * _HC_Ly * _HC_Lt));
    }
}

void CMeasurePolyakovU1XY::OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pAcceptGauge || EFT_GaugeU1 != pAcceptGauge->GetFieldType())
    {
        appCrucial(_T("CMeasurePolyakovU1XY only implemented with gauge U1!\n"));
        return;
    }
    const CFieldGaugeU1* pGaugeU1 = dynamic_cast<const CFieldGaugeU1*>(pAcceptGauge);

    dim3 block1(_HC_DecompX, _HC_DecompY, 1); 
    dim3 threads1(_HC_DecompLx, _HC_DecompLy, 1);
    for (UINT uiT = 0; uiT < _HC_Lt; ++uiT)
    {
        _kernelPolyakovLoopOfSiteU1 << <block1, threads1 >> > (pGaugeU1->m_pDeviceData, uiT, m_pTmpLoop);
    }
    _ZeroXYPlaneC(m_pXYDeviceLoopDensity);
    _ZeroXYPlane(m_pXYDeviceLoopDensityAbs);
    if (m_bMeasureZSlice)
    {
        dim3 blockz(_HC_DecompY, 1, 1);
        dim3 threadz(_HC_DecompLy, 1, 1);
        _kernelInitialZSliceU1 << <blockz , threadz >> > (m_pZDeviceLoopDensity, m_pZDeviceLoopDensityAbs);
    }
    _kernelPolyakovTraceOfSiteXYU1 << <block1, threads1 >> > (m_pTmpLoop, 
        m_pXYDeviceLoopDensity, m_pZDeviceLoopDensity,
        m_pXYDeviceLoopDensityAbs, m_pZDeviceLoopDensityAbs);
    checkCudaErrors(cudaMemcpy(m_pXYHostLoopDensity, m_pXYDeviceLoopDensity, sizeof(CLGComplex) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pXYHostLoopDensityAbs, m_pXYDeviceLoopDensityAbs, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
    if (m_bMeasureZSlice)
    {
        const Real fFactor = F(1.0) / static_cast<Real>(_HC_Lx * _HC_Ly);
        checkCudaErrors(cudaMemcpy(m_pZHostLoopDensity, m_pZDeviceLoopDensity, sizeof(CLGComplex) * _HC_Lz, cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(m_pZHostLoopDensityAbs, m_pZDeviceLoopDensityAbs, sizeof(Real) * _HC_Lz, cudaMemcpyDeviceToHost));
        for (UINT i = 0; i < _HC_Lz; ++i)
        {
            m_lstPZSlice.AddItem(cuCmulf_cr(m_pZHostLoopDensity[i], fFactor));
            m_lstPZSliceAbs.AddItem(m_pZHostLoopDensityAbs[i] / fFactor);
        }
    }
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

    if (m_bMeasureLoopZ)
    {
        dim3 block3(_HC_DecompX, 1, _HC_DecompZ);
        dim3 threads3(_HC_DecompLx, 1, _HC_DecompLz);
        _kernelPolyakovLoopOfSiteZU1 << <block3, threads3 >> > (pGaugeU1->m_pDeviceData, m_pTmpLoopZ);
        _ZeroXYPlaneC(m_pXYDeviceLoopDensity);
        _kernelPolyakovZTraceOfSiteXYU1 << <block3, threads3 >> > (m_pTmpLoopZ, m_pXYDeviceLoopDensity);

        checkCudaErrors(cudaMemcpy(m_pXYHostLoopDensity, m_pXYDeviceLoopDensity, sizeof(CLGComplex) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
        for (UINT i = CCommonData::m_sCenter.x; i < _HC_Lx; ++i)
        {
            m_lstLoopZDensity.AddItem(m_pXYHostLoopDensity[i * _HC_Ly + CCommonData::m_sCenter.y]);
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
            FALSE,
            m_byFieldId,
            m_lstPZ,
            &m_lstLoopZInner,
            m_lstLoopZ,
            m_lstR,
            m_uiConfigurationCount,
            F(1.0) / static_cast<Real>(_HC_Lt)
        );

        if (m_bShowResult)
        {
            appDetailed(_T("\n\n ==================== Polyakov LoopZ (%d con)============================ \n\n"), m_uiConfigurationCount);
        }

        if (m_bShowResult)
        {
            appSetLogDate(FALSE);
            appGeneral(_T("Loop Z is "));
            LogGeneralComplex(m_lstLoopZ[m_lstLoopZ.GetCount() - 1]);
            appGeneral(_T("\n"));
            appSetLogDate(TRUE);
            //appGeneral(_T("Loop is %f + %f I\n"), res[0].x, res[0].y);
        }
    }

    ++m_uiConfigurationCount;
}

void CMeasurePolyakovU1XY::Average(UINT )
{
    //nothing to do
}

void CMeasurePolyakovU1XY::Report()
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

void CMeasurePolyakovU1XY::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstLoop.RemoveAll();
    m_lstLoopInner.RemoveAll();
    m_lstLoopAbs.RemoveAll();
    m_lstLoopAbsInner.RemoveAll();
    m_lstLoopDensity.RemoveAll();
    m_lstLoopZ.RemoveAll();
    m_lstLoopZInner.RemoveAll();
    m_lstLoopZDensity.RemoveAll();

    m_lstR.RemoveAll();
    m_lstP.RemoveAll();
    m_lstPAbs.RemoveAll();
    m_lstPZ.RemoveAll();
    m_lstPZSlice.RemoveAll();
    m_lstPZSliceAbs.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================