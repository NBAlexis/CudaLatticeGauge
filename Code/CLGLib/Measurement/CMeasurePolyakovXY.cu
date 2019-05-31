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
    if (0 == uiT)
    {
        res[uiXYZ] = pDeviceBuffer[uiLinkIdx];
    }
    else
    {
        res[uiXYZ].Mul(pDeviceBuffer[uiLinkIdx]);
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

    if (0 == threadIdx.x && 0 == blockIdx.x)
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
    res[0].x = res[0].x / _HC_Volume_xyz;
    res[0].y = res[0].y / _HC_Volume_xyz;
    m_lstLoop.AddItem(res[0]);
    if (m_bShowResult)
    {
        appDetailed(_T("Loop is %f + %f i\n"), res[0].x, res[0].y);
    }
    for (UINT i = 1; i < _HC_Ly; ++i)
    {
        for (UINT j = 1; j < _HC_Lx; ++j)
        {
            m_lstLoopDensity.AddItem(m_pXYHostLoopDensity[j * _HC_Ly + i]);
            if (m_bShowResult)
            {
                appDetailed(_T("(%d,%d)=%1.6f %s %1.6f I   "), j, i, 
                    m_pXYHostLoopDensity[j * _HC_Ly + i].x,
                    m_pXYHostLoopDensity[j * _HC_Ly + i].y < F(0.0) ? _T("-") : _T("+"),
                    appAbs(m_pXYHostLoopDensity[j * _HC_Ly + i].y));
            }
        }
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
    assert(m_uiConfigurationCount * (_HC_Lx - 1) * (_HC_Ly - 1) == static_cast<UINT>(m_lstLoopDensity.Num()));

    appSetLogDate(FALSE);
    CLGComplex tmpChargeSum = _make_cuComplex(F(0.0), F(0.0));

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
    appGeneral(_T("\n ----------- average Loop |<P>| = %2.12f ------------- \n"), _cuCabsf(tmpChargeSum));

    appGeneral(_T("\n ----------- Loop density ------------- \n"));

    TArray<CLGComplex> tmp;
    appGeneral(_T("{\n"));
    for (UINT k = 0; k < m_uiConfigurationCount; ++k)
    {
        appGeneral(_T("{"));
        for (UINT i = 0; i < _HC_Ly - 1; ++i)
        {
            appGeneral(_T("{"));
            for (UINT j = 0; j < _HC_Lx - 1; ++j)
            {
                LogGeneralComplex(m_lstLoopDensity[k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j]);
                if (0 == k)
                {
                    tmp.AddItem(m_lstLoopDensity[k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j]);
                }
                else
                {
                    tmp[i * (_HC_Lx - 1) + j].x += m_lstLoopDensity[k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j].x;
                    tmp[i * (_HC_Lx - 1) + j].y += m_lstLoopDensity[k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j].y;
                }
            }
            appGeneral(_T("}, "));
        }
        appGeneral(_T("}\n"));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n ----------- loop density average ------------- \n"));

    for (UINT i = 0; i < _HC_Ly - 1; ++i)
    {
        for (UINT j = 0; j < _HC_Lx - 1; ++j)
        {
            tmp[i * (_HC_Lx - 1) + j].x = tmp[i * (_HC_Lx - 1) + j].x / m_uiConfigurationCount;
            tmp[i * (_HC_Lx - 1) + j].y = tmp[i * (_HC_Lx - 1) + j].y / m_uiConfigurationCount;
            appGeneral(_T("(x=%d,y=%d)="), j + 1, i + 1);
            LogGeneralComplex(tmp[i * (_HC_Lx - 1) + j]);
        }
        appGeneral(_T("\n"));
    }

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