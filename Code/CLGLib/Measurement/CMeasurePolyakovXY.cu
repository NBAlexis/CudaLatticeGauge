//=============================================================================
// FILENAME : CMeasurePolyakovXY.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/29/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasurePolyakovXY)

#pragma region kernles 

__global__ void
_CLG_LAUNCH_BOUND 
_kernelPolyakovLoopOfSite(
    const deviceSU3* __restrict__ pDeviceBuffer,
    UINT uiT,
    deviceSU3* res)
{
    UINT uiXYZ = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    UINT uiLinkIdx = (uiXYZ * _DC_Lt + uiT + 1) * _DC_Dir - 1;//(uiXYZ * _DC_Lt + uiT) * _DC_Dir + (_DC_Dir - 1);

    SSmallInt4 site4 = __deviceSiteIndexToInt4(uiXYZ * _DC_Lt + uiT);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

    if (0 == uiT)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, _DC_Dir - 1))
        {
            res[uiXYZ] = deviceSU3::makeSU3Zero();
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
            res[uiXYZ] = deviceSU3::makeSU3Zero();
        }
        else
        {
            res[uiXYZ].Mul(pDeviceBuffer[uiLinkIdx]);
        }
    }
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelPolyakovTraceOfSiteXY(
    const deviceSU3* __restrict__ resXYZ,
    CLGComplex* resXY)
{
    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
    UINT uiXYZ = uiXY * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    CLGComplex trres = resXYZ[uiXYZ].Tr();
    atomicAdd(&resXY[uiXY].x, trres.x);
    atomicAdd(&resXY[uiXY].y, trres.y);
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelPolyakovZeroXYPlane(
    CLGComplex* resXY,
    CLGComplex* total)
{
    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
    resXY[uiXY] = _make_cuComplex(F(0.0), F(0.0));

    if (0 == threadIdx.x && 0 == blockIdx.x && NULL != total)
    {
        total[0] = _make_cuComplex(F(0.0), F(0.0));
    }
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelPolyakovAverageOverZAndSum(
    CLGComplex* resXY,
    CLGComplex* total)
{
    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
    atomicAdd(&total[0].x, resXY[uiXY].x);
    atomicAdd(&total[0].y, resXY[uiXY].y);
    resXY[uiXY] = cuCdivf_cr(resXY[uiXY], _DC_Lz);
}

#pragma endregion

CLGAPI void _ZeroXYPlaneC(CLGComplex* pDeviceRes)
{
    preparethread;
    _kernelPolyakovZeroXYPlane << <block, threads >> > (pDeviceRes, NULL);
}

CMeasurePolyakovXY::~CMeasurePolyakovXY()
{
    if (NULL != m_pXYHostLoopDensity)
    {
        free(m_pXYHostLoopDensity);
    }

    if (NULL != m_pTmpDeviceSum)
    {
        checkCudaErrors(cudaFree(m_pTmpDeviceSum));
    }

    if (NULL != m_pXYDeviceLoopDensity)
    {
        checkCudaErrors(cudaFree(m_pXYDeviceLoopDensity));
    }

    if (NULL != m_pTmpLoop)
    {
        checkCudaErrors(cudaFree(m_pTmpLoop));
    }
}

void CMeasurePolyakovXY::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    m_pXYHostLoopDensity = (CLGComplex*)malloc(sizeof(CLGComplex) * _HC_Lx * _HC_Ly);
    checkCudaErrors(cudaMalloc((void**)&m_pTmpDeviceSum, sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)&m_pXYDeviceLoopDensity, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    checkCudaErrors(cudaMalloc((void**)&m_pTmpLoop, sizeof(deviceSU3) * _HC_Lx * _HC_Ly * _HC_Lz));
    Reset();

    INT iValue = 1;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;
}

void CMeasurePolyakovXY::OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pAcceptGauge || EFT_GaugeSU3 != pAcceptGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pAcceptGauge);

    dim3 block1(_HC_DecompX, _HC_DecompY, 1); 
    dim3 threads1(_HC_DecompLx, _HC_DecompLy, 1);

    dim3 block2(block1.x, 1, 1);
    dim3 threads2(threads1.x, 1, 1);

    _kernelPolyakovZeroXYPlane<<<block2, threads2 >>>(m_pXYDeviceLoopDensity, m_pTmpDeviceSum);

    for (UINT uiT = 0; uiT < _HC_Lt; ++uiT)
    {
        _kernelPolyakovLoopOfSite << <block1, threads1 >> >(pGaugeSU3->m_pDeviceData, uiT, m_pTmpLoop);
    }

    _kernelPolyakovTraceOfSiteXY << <block1, threads1 >> >(m_pTmpLoop, m_pXYDeviceLoopDensity);

    _kernelPolyakovAverageOverZAndSum << <block2, threads2 >> >(m_pXYDeviceLoopDensity, m_pTmpDeviceSum);

    //extract res
    CLGComplex res[1];
    checkCudaErrors(cudaMemcpy(res, m_pTmpDeviceSum, sizeof(CLGComplex), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pXYHostLoopDensity, m_pXYDeviceLoopDensity, sizeof(CLGComplex) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));

    ++m_uiConfigurationCount;
    if (m_bShowResult)
    {
        appDetailed(_T("\n\n ==================== Polyakov Loop (%d con)============================ \n\n"), m_uiConfigurationCount);
    }
    
    UINT uiSiteNumber = appGetLattice()->m_pIndexCache->m_uiSiteXYZ;
    res[0].x = res[0].x / uiSiteNumber;
    res[0].y = res[0].y / uiSiteNumber;
    m_lstLoop.AddItem(res[0]);
    if (m_bShowResult)
    {
        appGeneral(_T("Loop is %f + %f i\n"), res[0].x, res[0].y);
    }


    for (UINT i = CCommonData::m_sCenter.x; i < _HC_Lx; ++i)
    {
        m_lstLoopDensity.AddItem(m_pXYHostLoopDensity[
            i * _HC_Ly + CCommonData::m_sCenter.y]);
        if (m_bShowResult)
        {
            appDetailed(_T("(%d,%d)=%1.6f %s %1.6f I   "), i, CCommonData::m_sCenter.y,
                m_pXYHostLoopDensity[i * _HC_Ly + CCommonData::m_sCenter.y].x,
                m_pXYHostLoopDensity[i * _HC_Ly + CCommonData::m_sCenter.y].y < F(0.0) ? _T("") : _T("+"),
                appAbs(m_pXYHostLoopDensity[i * _HC_Ly + CCommonData::m_sCenter.y].y));
        }
    }
    if (m_bShowResult)
    {
        appDetailed(_T("\n"));
    }

    if (m_bShowResult)
    {
        appDetailed(_T("\n=====================================================\n"), m_uiConfigurationCount);
    }
}

void CMeasurePolyakovXY::Average(UINT )
{
    //nothing to do
}

void CMeasurePolyakovXY::Report()
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
    appGeneral(_T("\n ----------- average Loop |<P>| = %2.12f ------------- \n"), _cuCabsf(tmpChargeSum));

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

void CMeasurePolyakovXY::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstLoop.RemoveAll();
    m_lstLoopDensity.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================