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
#include "CMeasurePolyakovXY.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasurePolyakovXY)

#pragma region kernles 

__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovTraceOfSite_X(
    const cuDoubleComplex* __restrict__ resYZT,
    cuDoubleComplex* resY,
    cuDoubleComplex* resZ,
    cuDoubleComplex* resT,
    DOUBLE* resYAbs,
    DOUBLE* resZAbs,
    DOUBLE* resTAbs)
{
    intokernalInt4_Syzt_Only3D(0);

    if (NULL != resY)
    {
        atomicAdd(&resY[sSite4.y].x, resYZT[uiSiteIndex3DYZT].x);
        atomicAdd(&resY[sSite4.y].y, resYZT[uiSiteIndex3DYZT].y);
    }
    if (NULL != resZ)
    {
        atomicAdd(&resZ[sSite4.z].x, resYZT[uiSiteIndex3DYZT].x);
        atomicAdd(&resZ[sSite4.z].y, resYZT[uiSiteIndex3DYZT].y);
    }
    if (NULL != resT)
    {
        atomicAdd(&resT[sSite4.w].x, resYZT[uiSiteIndex3DYZT].x);
        atomicAdd(&resT[sSite4.w].y, resYZT[uiSiteIndex3DYZT].y);
    }

    if (NULL != resYAbs)
    {
        atomicAdd(&resYAbs[sSite4.y], cuCabs(resYZT[uiSiteIndex3DYZT]));
    }
    if (NULL != resZAbs)
    {
        atomicAdd(&resZAbs[sSite4.z], cuCabs(resYZT[uiSiteIndex3DYZT]));
    }
    if (NULL != resTAbs)
    {
        atomicAdd(&resTAbs[sSite4.w], cuCabs(resYZT[uiSiteIndex3DYZT]));
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovTraceOfSite_Y(
    const cuDoubleComplex* __restrict__ resXZT,
    cuDoubleComplex* resX,
    cuDoubleComplex* resZ,
    cuDoubleComplex* resT,
    DOUBLE* resXAbs,
    DOUBLE* resZAbs,
    DOUBLE* resTAbs)
{
    intokernalInt4_Sxzt_Only3D(0);

    if (NULL != resX)
    {
        atomicAdd(&resX[sSite4.x].x, resXZT[uiSiteIndex3DXZT].x);
        atomicAdd(&resX[sSite4.x].y, resXZT[uiSiteIndex3DXZT].y);
    }
    if (NULL != resZ)
    {
        atomicAdd(&resZ[sSite4.z].x, resXZT[uiSiteIndex3DXZT].x);
        atomicAdd(&resZ[sSite4.z].y, resXZT[uiSiteIndex3DXZT].y);
    }
    if (NULL != resT)
    {
        atomicAdd(&resT[sSite4.w].x, resXZT[uiSiteIndex3DXZT].x);
        atomicAdd(&resT[sSite4.w].y, resXZT[uiSiteIndex3DXZT].y);
    }

    if (NULL != resXAbs)
    {
        atomicAdd(&resXAbs[sSite4.x], cuCabs(resXZT[uiSiteIndex3DXZT]));
    }
    if (NULL != resZAbs)
    {
        atomicAdd(&resZAbs[sSite4.z], cuCabs(resXZT[uiSiteIndex3DXZT]));
    }
    if (NULL != resTAbs)
    {
        atomicAdd(&resTAbs[sSite4.w], cuCabs(resXZT[uiSiteIndex3DXZT]));
    }
}


__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovTraceOfSite_Z(
    const cuDoubleComplex* __restrict__ resXYT,
    cuDoubleComplex* resXY,
    cuDoubleComplex* resX,
    cuDoubleComplex* resY,
    cuDoubleComplex* resT,
    DOUBLE* resXYAbs,
    DOUBLE* resXAbs,
    DOUBLE* resYAbs,
    DOUBLE* resTAbs)
{
    intokernalInt4_Sxyt_Only3D(0);

    const UINT uiXY = sSite4.x * _DC_Ly + sSite4.y;

    if (NULL != resXY)
    {
        atomicAdd(&resXY[uiXY].x, resXYT[uiSiteIndex3DXYT].x);
        atomicAdd(&resXY[uiXY].y, resXYT[uiSiteIndex3DXYT].y);
    }
    if (NULL != resX)
    {
        atomicAdd(&resX[sSite4.x].x, resXYT[uiSiteIndex3DXYT].x);
        atomicAdd(&resX[sSite4.x].y, resXYT[uiSiteIndex3DXYT].y);
    }
    if (NULL != resY)
    {
        atomicAdd(&resY[sSite4.y].x, resXYT[uiSiteIndex3DXYT].x);
        atomicAdd(&resY[sSite4.y].y, resXYT[uiSiteIndex3DXYT].y);
    }
    if (NULL != resT)
    {
        atomicAdd(&resT[sSite4.w].x, resXYT[uiSiteIndex3DXYT].x);
        atomicAdd(&resT[sSite4.w].y, resXYT[uiSiteIndex3DXYT].y);
    }

    if (NULL != resXYAbs)
    {
        atomicAdd(&resXYAbs[uiXY], cuCabs(resXYT[uiSiteIndex3DXYT]));
    }
    if (NULL != resXAbs)
    {
        atomicAdd(&resXAbs[sSite4.x], cuCabs(resXYT[uiSiteIndex3DXYT]));
    }
    if (NULL != resYAbs)
    {
        atomicAdd(&resYAbs[sSite4.y], cuCabs(resXYT[uiSiteIndex3DXYT]));
    }
    if (NULL != resTAbs)
    {
        atomicAdd(&resTAbs[sSite4.w], cuCabs(resXYT[uiSiteIndex3DXYT]));
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovTraceOfSite_T(
    const cuDoubleComplex* __restrict__ resXYZ,
    cuDoubleComplex* resXY,
    cuDoubleComplex* resX,
    cuDoubleComplex* resY,
    cuDoubleComplex* resZ,
    DOUBLE* resXYAbs,
    DOUBLE* resXAbs,
    DOUBLE* resYAbs,
    DOUBLE* resZAbs)
{
    intokernalInt4_S_Only3D(0);

    const UINT uiXY = sSite4.x * _DC_Ly + sSite4.y;

    if (NULL != resXY)
    {
        atomicAdd(&resXY[uiXY].x, resXYZ[uiSiteIndex3D].x);
        atomicAdd(&resXY[uiXY].y, resXYZ[uiSiteIndex3D].y);
    }
    if (NULL != resX)
    {
        atomicAdd(&resX[sSite4.x].x, resXYZ[uiSiteIndex3D].x);
        atomicAdd(&resX[sSite4.x].y, resXYZ[uiSiteIndex3D].y);
    }
    if (NULL != resY)
    {
        atomicAdd(&resY[sSite4.y].x, resXYZ[uiSiteIndex3D].x);
        atomicAdd(&resY[sSite4.y].y, resXYZ[uiSiteIndex3D].y);
    }
    if (NULL != resZ)
    {
        atomicAdd(&resZ[sSite4.z].x, resXYZ[uiSiteIndex3D].x);
        atomicAdd(&resZ[sSite4.z].y, resXYZ[uiSiteIndex3D].y);
    }

    if (NULL != resXYAbs)
    {
        atomicAdd(&resXYAbs[uiXY], cuCabs(resXYZ[uiSiteIndex3D]));
    }
    if (NULL != resXAbs)
    {
        atomicAdd(&resXAbs[sSite4.x], cuCabs(resXYZ[uiSiteIndex3D]));
    }
    if (NULL != resYAbs)
    {
        atomicAdd(&resYAbs[sSite4.y], cuCabs(resXYZ[uiSiteIndex3D]));
    }
    if (NULL != resZAbs)
    {
        atomicAdd(&resZAbs[sSite4.z], cuCabs(resXYZ[uiSiteIndex3D]));
    }
}



//__global__ void
//_CLG_LAUNCH_BOUND 
//_kernelPolyakovLoopOfSite(
//    BYTE byFieldId,
//    const deviceSU3* __restrict__ pDeviceBuffer,
//    UINT uiT,
//    deviceSU3* res)
//{
//    UINT uiXYZ = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
//    const UINT uiSiteIndex = uiXYZ * _DC_Lt + uiT;
//    UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, _DC_Dir - 1);
//    //(uiSiteIndex + 1) * _DC_Dir - 1;//uiSiteIndex * _DC_Dir + (_DC_Dir - 1);
//    //if (0 == uiXYZ)
//    //{
//    //    printf("t=%d, site=%d, linkidx=%d\n", uiT, uiSiteIndex, uiLinkIdx);
//    //}
//
//    const SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
//    const UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);
//
//    if (0 == uiT)
//    {
//        if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, _DC_Dir - 1))
//        {
//            res[uiXYZ] = deviceSU3::makeSU3Zero();
//        }
//        else
//        {
//            res[uiXYZ] = pDeviceBuffer[uiLinkIdx];
//        }
//    }
//    else
//    {
//        if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, _DC_Dir - 1))
//        {
//            res[uiXYZ] = deviceSU3::makeSU3Zero();
//        }
//        else
//        {
//            res[uiXYZ].Mul(pDeviceBuffer[uiLinkIdx]);
//        }
//    }
//}
//
///**
// * Before call me, set block dim thread dim.y = 1
// */
//__global__ void
//_CLG_LAUNCH_BOUND
//_kernelPolyakovLoopOfSiteZ(
//    BYTE byFieldId,
//    const deviceSU3* __restrict__ pDeviceBuffer,
//    deviceSU3* res)
//{
//    UINT uiXYT = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z);
//
//    for (UINT z = 0; z < _DC_Lz; ++z)
//    {
//        const UINT uiSiteIndex = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + z * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z);
//        UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, 2);
//        const SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
//        const UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);
//
//        if (0 == z)
//        {
//            if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 2))
//            {
//                res[uiXYT] = deviceSU3::makeSU3Zero();
//            }
//            else
//            {
//                res[uiXYT] = pDeviceBuffer[uiLinkIdx];
//            }
//        }
//        else
//        {
//            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 2))
//            {
//                res[uiXYT].Mul(pDeviceBuffer[uiLinkIdx]);
//            }
//        }
//    }
//}
//
//
//
//__global__ void
//_CLG_LAUNCH_BOUND
//_kernelInitialZSlice(CLGComplex* resZ, Real* resZAbs)
//{
//    resZ[threadIdx.x + blockIdx.x * blockDim.x] = _zeroc;
//    resZAbs[threadIdx.x + blockIdx.x * blockDim.x] = F(0.0);
//}
//
//__global__ void
//_CLG_LAUNCH_BOUND
//_kernelPolyakovZTraceOfSiteXY(
//    const deviceSU3* __restrict__ resXYT,
//    CLGComplex* resXY)
//{
//    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
//    UINT uiXYT = uiXY * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z);
//    const CLGComplex trres = resXYT[uiXYT].Tr();
//    atomicAdd(&resXY[uiXY].x, trres.x);
//    atomicAdd(&resXY[uiXY].y, trres.y);
//}


#pragma endregion


CMeasurePolyakovXY::~CMeasurePolyakovXY()
{
    appSafeFree(m_pXYHostLoopDensity);
    appSafeFree(m_pXHostLoopDensity);
    appSafeFree(m_pYHostLoopDensity);
    appSafeFree(m_pZHostLoopDensity);
    appSafeFree(m_pTHostLoopDensity);
    appSafeCudaFree(m_pXYDeviceLoopDensity);
    appSafeCudaFree(m_pXDeviceLoopDensity);
    appSafeCudaFree(m_pYDeviceLoopDensity);
    appSafeCudaFree(m_pZDeviceLoopDensity);
    appSafeCudaFree(m_pTDeviceLoopDensity);
    appSafeFree(m_pXYHostLoopDensityAbs);
    appSafeFree(m_pXHostLoopDensityAbs);
    appSafeFree(m_pYHostLoopDensityAbs);
    appSafeFree(m_pZHostLoopDensityAbs);
    appSafeFree(m_pTHostLoopDensityAbs);
    appSafeCudaFree(m_pXYDeviceLoopDensityAbs);
    appSafeCudaFree(m_pXDeviceLoopDensityAbs);
    appSafeCudaFree(m_pYDeviceLoopDensityAbs);
    appSafeCudaFree(m_pZDeviceLoopDensityAbs);
    appSafeCudaFree(m_pTDeviceLoopDensityAbs);

    appSafeCudaFree(m_pTmpLoop);
    appSafeCudaFree(m_pTmpLoopX);
    appSafeCudaFree(m_pTmpLoopY);
    appSafeCudaFree(m_pTmpLoopZ);

    appSafeFree(m_pHostDistributionR);
    appSafeFree(m_pHostDistributionP);
    appSafeFree(m_pHostDistributionPAbs);
    appSafeCudaFree(m_pDistributionR);
    appSafeCudaFree(m_pDistributionP);
    appSafeCudaFree(m_pDistributionPAbs);
}

void CMeasurePolyakovXY::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    checkCudaErrors(__cudaMalloc((void**)&m_pTmpLoop, sizeof(cuDoubleComplex) * _HC_Lx * _HC_Ly * _HC_Lz));

    Reset();

    INT iValue = 0;
    param.FetchValueINT(_T("MeasureX"), iValue);
    m_bMeasureLoopX = iValue != 0;
    iValue = 0;
    param.FetchValueINT(_T("MeasureY"), iValue);
    m_bMeasureLoopY = iValue != 0;
    iValue = 0;
    param.FetchValueINT(_T("MeasureZ"), iValue);
    m_bMeasureLoopZ = iValue != 0;

    if (m_bMeasureLoopX)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pTmpLoopX, sizeof(cuDoubleComplex) * _HC_Ly * _HC_Lz * _HC_Lt));
    }
    if (m_bMeasureLoopY)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pTmpLoopY, sizeof(cuDoubleComplex) * _HC_Lx * _HC_Lz * _HC_Lt));
    }
    if (m_bMeasureLoopZ)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pTmpLoopZ, sizeof(cuDoubleComplex) * _HC_Lx * _HC_Ly * _HC_Lt));
    }

    iValue = 0;
    param.FetchValueINT(_T("XSlice"), iValue);
    m_bMeasureXSlice = iValue != 0;
    iValue = 0;
    param.FetchValueINT(_T("YSlice"), iValue);
    m_bMeasureYSlice = iValue != 0;
    iValue = 0;
    param.FetchValueINT(_T("ZSlice"), iValue);
    m_bMeasureZSlice = iValue != 0;
    iValue = 0;
    param.FetchValueINT(_T("TSlice"), iValue);
    m_bMeasureTSlice = iValue != 0;

    iValue = 1;
    param.FetchValueINT(_T("MeasureDist"), iValue);
    m_bMeasureDistribution = iValue != 0;
    iValue = 0;
    param.FetchValueINT(_T("Absolute"), iValue);
    m_bMeasureAbs = iValue != 0;

    if (m_bMeasureDistribution)
    {
        m_pXYHostLoopDensity = (cuDoubleComplex*)malloc(sizeof(cuDoubleComplex) * _HC_Lx * _HC_Ly);
        checkCudaErrors(__cudaMalloc((void**)&m_pXYDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lx * _HC_Ly));
        if (m_bMeasureAbs)
        {
            m_pXYHostLoopDensityAbs = (DOUBLE*)malloc(sizeof(DOUBLE) * _HC_Lx * _HC_Ly);
            checkCudaErrors(__cudaMalloc((void**)&m_pXYDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lx * _HC_Ly));
        }

        iValue = 0;
        param.FetchValueINT(_T("ShiftCenter"), iValue);
        m_bShiftCenter = iValue != 0;

        //assuming the center is really at center
        SetMaxAndEdge(&m_uiMaxR, &m_uiEdgeR, m_bShiftCenter);
        checkCudaErrors(__cudaMalloc((void**)&m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1)));
        checkCudaErrors(__cudaMalloc((void**)&m_pDistributionP, sizeof(cuDoubleComplex) * (m_uiMaxR + 1)));
        if (m_bMeasureAbs)
        {
            checkCudaErrors(__cudaMalloc((void**)&m_pDistributionPAbs, sizeof(DOUBLE) * (m_uiMaxR + 1)));
        }

        m_pHostDistributionR = (UINT*)malloc(sizeof(UINT) * (m_uiMaxR + 1));
        m_pHostDistributionP = (cuDoubleComplex*)malloc(sizeof(cuDoubleComplex) * (m_uiMaxR + 1));
        if (m_bMeasureAbs)
        {
            m_pHostDistributionPAbs = (DOUBLE*)malloc(sizeof(DOUBLE) * (m_uiMaxR + 1));
        }
    }

    if (m_bMeasureXSlice)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pXDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lx));
        m_pXHostLoopDensity = (cuDoubleComplex*)malloc(sizeof(cuDoubleComplex) * _HC_Lx);
        if (m_bMeasureAbs)
        {
            checkCudaErrors(__cudaMalloc((void**)&m_pXDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lx));
            m_pXHostLoopDensityAbs = (DOUBLE*)malloc(sizeof(DOUBLE) * _HC_Lx);
        }
    }

    if (m_bMeasureYSlice)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pYDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Ly));
        m_pYHostLoopDensity = (cuDoubleComplex*)malloc(sizeof(cuDoubleComplex) * _HC_Ly);
        if (m_bMeasureAbs)
        {
            checkCudaErrors(__cudaMalloc((void**)&m_pYDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Ly));
            m_pYHostLoopDensityAbs = (DOUBLE*)malloc(sizeof(DOUBLE) * _HC_Ly);
        }
    }

    if (m_bMeasureZSlice)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pZDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lz));
        m_pZHostLoopDensity = (cuDoubleComplex*)malloc(sizeof(cuDoubleComplex) * _HC_Lz);
        if (m_bMeasureAbs)
        {
            checkCudaErrors(__cudaMalloc((void**)&m_pZDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lz));
            m_pZHostLoopDensityAbs = (DOUBLE*)malloc(sizeof(DOUBLE) * _HC_Lz);
        }
    }

    if (m_bMeasureTSlice)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pTDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lt));
        m_pTHostLoopDensity = (cuDoubleComplex*)malloc(sizeof(cuDoubleComplex) * _HC_Lt);
        if (m_bMeasureAbs)
        {
            checkCudaErrors(__cudaMalloc((void**)&m_pTDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lt));
            m_pTHostLoopDensityAbs = (DOUBLE*)malloc(sizeof(DOUBLE) * _HC_Lt);
        }
    }
}

/**
* Expose one per-slice profile (the last uiSliceCount values of the list, i.e. the
* values of the current configuration) to the measurement dictionary
*/
static void _AddSliceToDictionary(CMeasurementManager* pOwner, CMeasure* pMeasure, const CCString& sName,
    const TArray<cuDoubleComplex>& lstSlice, const TArray<DOUBLE>& lstSliceAbs, UINT uiSliceCount, UBOOL bMeasureAbs)
{
    TArray<CLGComplex> thisConfiguration;
    for (UINT i = 0; i < uiSliceCount; ++i)
    {
        const cuDoubleComplex& c = lstSlice[lstSlice.Num() - uiSliceCount + i];
        thisConfiguration.AddItem(_make_cuComplex(static_cast<Real>(c.x), static_cast<Real>(c.y)));
    }
    pOwner->AddOneConfigurationResult(pMeasure, sName, thisConfiguration);
    if (bMeasureAbs)
    {
        TArray<DOUBLE> thisConfigurationAbs;
        for (UINT i = 0; i < uiSliceCount; ++i)
        {
            thisConfigurationAbs.AddItem(lstSliceAbs[lstSliceAbs.Num() - uiSliceCount + i]);
        }
        pOwner->AddOneConfigurationResult(pMeasure, sName + _T("Abs"), thisConfigurationAbs);
    }
}

void CMeasurePolyakovXY::OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pAcceptGauge)
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }

#if _CLG_MULTI_GPU
    //I9: the T-direction Polyakov chain below is a per-site PRODUCT along t;
    //with a split t direction each rank holds only a partial chain and the
    //loop average (and the t-slice profile) is silently wrong. Reject
    //explicitly, same rule as the loopX/Y/Z chain guards below.
    if (NULL != appGetComm() && appGetComm()->GpuGrid()[3] > 1)
    {
        appCrucial(_T("CMeasurePolyakovXY: not supported on multi-GPU with a split t direction (cross-rank t chain not implemented). Rejected.\n"));
        return;
    }
#endif

    pAcceptGauge->PolyakovOnSpatialSite(m_pTmpLoop, 3);
    if (m_bMeasureDistribution)
    {
        _ZeroXYPlane(m_pXYDeviceLoopDensity);
        if (m_bMeasureAbs)
        {
            _ZeroXYPlane(m_pXYDeviceLoopDensityAbs);
        }
    }
    if (m_bMeasureXSlice)
    {
        _ZeroSlice(m_pXDeviceLoopDensity, 0);
        if (m_bMeasureAbs)
        {
            _ZeroSlice(m_pXDeviceLoopDensityAbs, 0);
        }
    }
    if (m_bMeasureYSlice)
    {
        _ZeroSlice(m_pYDeviceLoopDensity, 1);
        if (m_bMeasureAbs)
        {
            _ZeroSlice(m_pYDeviceLoopDensityAbs, 1);
        }
    }
    if (m_bMeasureZSlice)
    {
        _ZeroSlice(m_pZDeviceLoopDensity, 2);
        if (m_bMeasureAbs)
        {
            _ZeroSlice(m_pZDeviceLoopDensityAbs, 2);
        }
    }
    preparethread_S;
    _LAUNCH_KERNEL(_kernelPolyakovTraceOfSite_T, block3d, threads3d, m_pTmpLoop, 
        m_bMeasureDistribution ? m_pXYDeviceLoopDensity : NULL,
        m_bMeasureXSlice ? m_pXDeviceLoopDensity : NULL,
        m_bMeasureYSlice ? m_pYDeviceLoopDensity : NULL,
        m_bMeasureZSlice ? m_pZDeviceLoopDensity : NULL,
        (m_bMeasureDistribution && m_bMeasureAbs) ? m_pXYDeviceLoopDensityAbs : NULL,
        (m_bMeasureXSlice && m_bMeasureAbs) ? m_pXDeviceLoopDensityAbs : NULL,
        (m_bMeasureYSlice && m_bMeasureAbs) ? m_pYDeviceLoopDensityAbs : NULL,
        (m_bMeasureZSlice && m_bMeasureAbs) ? m_pZDeviceLoopDensityAbs : NULL);

    if (m_bMeasureDistribution)
    {
        checkCudaErrors(cudaMemcpy(m_pXYHostLoopDensity, m_pXYDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
        if (m_bMeasureAbs)
        {
            checkCudaErrors(cudaMemcpy(m_pXYHostLoopDensityAbs, m_pXYDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
        }
#if _CLG_MULTI_GPU
        //P4-3.3: the XY density holds the per-rank partial sum over the LOCAL
        //z extent (the kernel accumulates only the local z-slice of the loop
        //product). Sum the per-(x,y) values across ranks so the R-distribution
        //transform below sees the GLOBAL density; then write back to the device
        //array (the transform kernels read the device buffer).
        if (NULL != appGetComm())
        {
            //P4-3.3: requires x/y NOT split (same (x,y) index set on every rank);
            //a split x/y would need a global-index gather, out of P4-3.3 scope.
            if (appGetComm()->GpuGrid()[0] > 1 || appGetComm()->GpuGrid()[1] > 1)
            {
                appCrucial(_T("CMeasurePolyakovXY: the XY distribution is not supported on multi-GPU with a split x/y direction. Rejected.\n"));
                return;
            }
            appGlobalSum(m_pXYHostLoopDensity, _HC_Lx * _HC_Ly);
            if (m_bMeasureAbs)
            {
                GlobalSumRealArray(m_pXYHostLoopDensityAbs, _HC_Lx * _HC_Ly);
            }
            checkCudaErrors(cudaMemcpy(m_pXYDeviceLoopDensity, m_pXYHostLoopDensity, sizeof(cuDoubleComplex) * _HC_Lx * _HC_Ly, cudaMemcpyHostToDevice));
            if (m_bMeasureAbs)
            {
                checkCudaErrors(cudaMemcpy(m_pXYDeviceLoopDensityAbs, m_pXYHostLoopDensityAbs, sizeof(DOUBLE) * _HC_Lx * _HC_Ly, cudaMemcpyHostToDevice));
            }
        }
#endif

        TransformFromXYDataToRDataOnce(
            m_bShiftCenter,
            m_pXYDeviceLoopDensity,
            m_pDistributionR,
            m_pDistributionP,
            m_pHostDistributionR,
            m_pHostDistributionP,
            m_uiMaxR,
            m_uiEdgeR,
            TRUE,
            GetGaugeFieldIdSingleField(),
            m_lstP,
            &m_lstLoopInner,
            m_lstLoop,
            m_lstR,
            m_uiConfigurationCount,
            F(1.0) / static_cast<Real>(GlobalL(2))
        );

        if (m_bMeasureAbs)
        {
            TransformFromXYDataToRDataOnce(
                m_bShiftCenter,
                m_pXYDeviceLoopDensityAbs,
                m_pDistributionR,
                m_pDistributionPAbs,
                m_pHostDistributionR,
                m_pHostDistributionPAbs,
                m_uiMaxR,
                m_uiEdgeR,
                FALSE,
                GetGaugeFieldIdSingleField(),
                m_lstPAbs,
                &m_lstLoopAbsInner,
                m_lstLoopAbs,
                m_lstR,
                m_uiConfigurationCount,
                F(1.0) / static_cast<Real>(GlobalL(2))
            );
        }
    }
    else
    {
        //P4-3.3: ReduceComplex returns the LOCAL sum over this rank's spatial
        //sites (the T-direction loop product is complete per site when t is not
        //split; for a split t the loop chain itself is incomplete and this
        //measurement is rejected at the top of this function -- see the MG
        //guard there).
        const cuDoubleComplex cSum = appGetCudaHelper()->ReduceComplex(m_pTmpLoop, _HC_Lx * _HC_Ly * _HC_Lz);
        cuDoubleComplex cLoopAvg = cSum;
        appGlobalSum(cLoopAvg);
        //Normalize by the GLOBAL spatial volume (local volume is the per-rank
        //share; single-GPU unchanged since local == global there).
        const DOUBLE fVol = GlobalL(0) * GlobalL(1) * GlobalL(2);
        cLoopAvg = cuCdivf_cd_host(cLoopAvg, fVol);
        m_lstLoop.AddItem(cLoopAvg);
        if (m_bMeasureAbs)
        {
            //also used by the result print and the dictionary
            m_lstLoopAbs.AddItem(cuCabs(cLoopAvg));
        }
    }

    if (m_bMeasureXSlice)
    {
#if _CLG_MULTI_GPU
        //I9: the x-slice profile lives on the LOCAL x extent; with a split x
        //direction each rank holds a different global-x slice and assembling
        //the global profile is a gather, not a sum. Reject explicitly instead
        //of reporting a silently-wrong partial profile (same rule as ZSlice).
        if (NULL != appGetComm() && appGetComm()->GpuGrid()[0] > 1)
        {
            appCrucial(_T("CMeasurePolyakovXY: m_bMeasureXSlice is not supported on multi-GPU with a split x direction (global-x gather not implemented). Rejected.\n"));
            return;
        }
#endif
        //P4-3.3: per-(y,z)-slice sums are per-rank partials; reduce the array.
        const DOUBLE fFactor = 1.0 / static_cast<DOUBLE>(GlobalL(1) * GlobalL(2));
        checkCudaErrors(cudaMemcpy(m_pXHostLoopDensity, m_pXDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lx, cudaMemcpyDeviceToHost));
        if (m_bMeasureAbs)
        {
            checkCudaErrors(cudaMemcpy(m_pXHostLoopDensityAbs, m_pXDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lx, cudaMemcpyDeviceToHost));
        }
#if _CLG_MULTI_GPU
        if (NULL != appGetComm())
        {
            appGlobalSum(m_pXHostLoopDensity, _HC_Lx);
            if (m_bMeasureAbs)
            {
                GlobalSumRealArray(m_pXHostLoopDensityAbs, _HC_Lx);
            }
        }
#endif
        for (UINT i = 0; i < _HC_Lx; ++i)
        {
            m_lstP_XSlice.AddItem(cuCmulf_cd(m_pXHostLoopDensity[i], fFactor));
            if (m_bMeasureAbs)
            {
                m_lstP_XSliceAbs.AddItem(m_pXHostLoopDensityAbs[i] * fFactor);
            }
        }
    }
    if (m_bMeasureYSlice)
    {
#if _CLG_MULTI_GPU
        //I9: the y-slice profile lives on the LOCAL y extent; with a split y
        //direction each rank holds a different global-y slice and assembling
        //the global profile is a gather, not a sum. Reject explicitly instead
        //of reporting a silently-wrong partial profile (same rule as ZSlice).
        if (NULL != appGetComm() && appGetComm()->GpuGrid()[1] > 1)
        {
            appCrucial(_T("CMeasurePolyakovXY: m_bMeasureYSlice is not supported on multi-GPU with a split y direction (global-y gather not implemented). Rejected.\n"));
            return;
        }
#endif
        //P4-3.3: see X slice; per-(x,z)-slice sums reduced across ranks.
        const DOUBLE fFactor = 1.0 / static_cast<DOUBLE>(GlobalL(0) * GlobalL(2));
        checkCudaErrors(cudaMemcpy(m_pYHostLoopDensity, m_pYDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Ly, cudaMemcpyDeviceToHost));
        if (m_bMeasureAbs)
        {
            checkCudaErrors(cudaMemcpy(m_pYHostLoopDensityAbs, m_pYDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Ly, cudaMemcpyDeviceToHost));
        }
#if _CLG_MULTI_GPU
        if (NULL != appGetComm())
        {
            appGlobalSum(m_pYHostLoopDensity, _HC_Ly);
            if (m_bMeasureAbs)
            {
                GlobalSumRealArray(m_pYHostLoopDensityAbs, _HC_Ly);
            }
        }
#endif
        for (UINT i = 0; i < _HC_Ly; ++i)
        {
            m_lstP_YSlice.AddItem(cuCmulf_cd(m_pYHostLoopDensity[i], fFactor));
            if (m_bMeasureAbs)
            {
                m_lstP_YSliceAbs.AddItem(m_pYHostLoopDensityAbs[i] * fFactor);
            }
        }
    }
    if (m_bMeasureZSlice)
    {
#if _CLG_MULTI_GPU
        //P4-3.3: the z-slice profile lives on the LOCAL z extent; with a split z
        //direction (GpuGridZ > 1) each rank holds a different global-z slice and
        //assembling the global profile is a gather, not a sum. Reject explicitly
        //instead of reporting a silently-wrong partial profile.
        if (NULL != appGetComm() && appGetComm()->GpuGrid()[2] > 1)
        {
            appCrucial(_T("CMeasurePolyakovXY: m_bMeasureZSlice is not supported on multi-GPU with a split z direction (global-z gather not implemented in P4-3.3). Rejected.\n"));
            return;
        }
#endif
        const DOUBLE fFactor = 1.0 / static_cast<DOUBLE>(GlobalL(0) * GlobalL(1));
        checkCudaErrors(cudaMemcpy(m_pZHostLoopDensity, m_pZDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lz, cudaMemcpyDeviceToHost));
        if (m_bMeasureAbs)
        {
            checkCudaErrors(cudaMemcpy(m_pZHostLoopDensityAbs, m_pZDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lz, cudaMemcpyDeviceToHost));
        }
#if _CLG_MULTI_GPU
        if (NULL != appGetComm())
        {
            appGlobalSum(m_pZHostLoopDensity, _HC_Lz);
            if (m_bMeasureAbs)
            {
                GlobalSumRealArray(m_pZHostLoopDensityAbs, _HC_Lz);
            }
        }
#endif
        for (UINT i = 0; i < _HC_Lz; ++i)
        {
            m_lstP_ZSlice.AddItem(cuCmulf_cd(m_pZHostLoopDensity[i], fFactor));
            if (m_bMeasureAbs)
            {
                m_lstP_ZSliceAbs.AddItem(m_pZHostLoopDensityAbs[i] * fFactor);
            }
        }
    }

    if (m_bMeasureLoopX)
    {
#if _CLG_MULTI_GPU
        //I9: the X-direction Polyakov chain is a per-site product along x;
        //with a split x direction each rank holds only a partial chain and the
        //loop average is silently wrong. Reject (same rule as loopZ below).
        if (NULL != appGetComm() && appGetComm()->GpuGrid()[0] > 1)
        {
            appCrucial(_T("CMeasurePolyakovXY: loopX is not supported on multi-GPU with a split x direction (cross-rank x chain not implemented). Rejected.\n"));
            return;
        }
#endif
        pAcceptGauge->PolyakovOnSpatialSite(m_pTmpLoopX, 0);
        if (m_bMeasureYSlice)
        {
            _ZeroSlice(m_pYDeviceLoopDensity, 1);
            if (m_bMeasureAbs)
            {
                _ZeroSlice(m_pYDeviceLoopDensityAbs, 1);
            }
        }
        if (m_bMeasureZSlice)
        {
            _ZeroSlice(m_pZDeviceLoopDensity, 2);
            if (m_bMeasureAbs)
            {
                _ZeroSlice(m_pZDeviceLoopDensityAbs, 2);
            }
        }
        if (m_bMeasureTSlice)
        {
            _ZeroSlice(m_pTDeviceLoopDensity, 3);
            if (m_bMeasureAbs)
            {
                _ZeroSlice(m_pTDeviceLoopDensityAbs, 3);
            }
        }
        preparethread_Syzt;
        _LAUNCH_KERNEL(_kernelPolyakovTraceOfSite_X, block3dyzt, threads3dyzt, m_pTmpLoopX,
            m_bMeasureYSlice ? m_pYDeviceLoopDensity : NULL,
            m_bMeasureZSlice ? m_pZDeviceLoopDensity : NULL,
            m_bMeasureTSlice ? m_pTDeviceLoopDensity : NULL,
            (m_bMeasureYSlice && m_bMeasureAbs) ? m_pYDeviceLoopDensityAbs : NULL,
            (m_bMeasureZSlice && m_bMeasureAbs) ? m_pZDeviceLoopDensityAbs : NULL,
            (m_bMeasureTSlice && m_bMeasureAbs) ? m_pTDeviceLoopDensityAbs : NULL);

        if (m_bMeasureYSlice)
        {
            //P4-3.3: see main-section slices; per-(x,z)-slice partial sums reduced.
            const DOUBLE fFactor = 1.0 / static_cast<DOUBLE>(GlobalL(2) * GlobalL(3));
            checkCudaErrors(cudaMemcpy(m_pYHostLoopDensity, m_pYDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Ly, cudaMemcpyDeviceToHost));
            if (m_bMeasureAbs)
            {
                checkCudaErrors(cudaMemcpy(m_pYHostLoopDensityAbs, m_pYDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Ly, cudaMemcpyDeviceToHost));
            }
#if _CLG_MULTI_GPU
            if (NULL != appGetComm())
            {
                appGlobalSum(m_pYHostLoopDensity, _HC_Ly);
                if (m_bMeasureAbs)
                {
                    GlobalSumRealArray(m_pYHostLoopDensityAbs, _HC_Ly);
                }
            }
#endif
            for (UINT i = 0; i < _HC_Ly; ++i)
            {
                m_lstPX_YSlice.AddItem(cuCmulf_cd(m_pYHostLoopDensity[i], fFactor));
                if (m_bMeasureAbs)
                {
                    m_lstPX_YSliceAbs.AddItem(m_pYHostLoopDensityAbs[i] * fFactor);
                }
            }
        }
        if (m_bMeasureZSlice)
        {
#if _CLG_MULTI_GPU
            //P4-3.3: same global-z gather limitation as the main-section Z slice.
            if (NULL != appGetComm() && appGetComm()->GpuGrid()[2] > 1)
            {
                appCrucial(_T("CMeasurePolyakovXY: loopX Z-slice is not supported on multi-GPU with a split z direction. Rejected.\n"));
                return;
            }
#endif
            const DOUBLE fFactor = 1.0 / static_cast<DOUBLE>(GlobalL(1) * GlobalL(3));
            checkCudaErrors(cudaMemcpy(m_pZHostLoopDensity, m_pZDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lz, cudaMemcpyDeviceToHost));
            if (m_bMeasureAbs)
            {
                checkCudaErrors(cudaMemcpy(m_pZHostLoopDensityAbs, m_pZDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lz, cudaMemcpyDeviceToHost));
            }
#if _CLG_MULTI_GPU
            if (NULL != appGetComm())
            {
                appGlobalSum(m_pZHostLoopDensity, _HC_Lz);
                if (m_bMeasureAbs)
                {
                    GlobalSumRealArray(m_pZHostLoopDensityAbs, _HC_Lz);
                }
            }
#endif
            for (UINT i = 0; i < _HC_Lz; ++i)
            {
                m_lstPX_ZSlice.AddItem(cuCmulf_cd(m_pZHostLoopDensity[i], fFactor));
                if (m_bMeasureAbs)
                {
                    m_lstPX_ZSliceAbs.AddItem(m_pZHostLoopDensityAbs[i] * fFactor);
                }
            }
        }
        if (m_bMeasureTSlice)
        {
            //P4-3.3: t is not split in the P4-3.3 test grids; reduce anyway.
            const DOUBLE fFactor = 1.0 / static_cast<DOUBLE>(GlobalL(1) * GlobalL(2));
            checkCudaErrors(cudaMemcpy(m_pTHostLoopDensity, m_pTDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lt, cudaMemcpyDeviceToHost));
            if (m_bMeasureAbs)
            {
                checkCudaErrors(cudaMemcpy(m_pTHostLoopDensityAbs, m_pTDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lt, cudaMemcpyDeviceToHost));
            }
#if _CLG_MULTI_GPU
            if (NULL != appGetComm())
            {
                appGlobalSum(m_pTHostLoopDensity, _HC_Lt);
                if (m_bMeasureAbs)
                {
                    GlobalSumRealArray(m_pTHostLoopDensityAbs, _HC_Lt);
                }
            }
#endif
            for (UINT i = 0; i < _HC_Lt; ++i)
            {
                m_lstPX_TSlice.AddItem(cuCmulf_cd(m_pTHostLoopDensity[i], fFactor));
                if (m_bMeasureAbs)
                {
                    m_lstPX_TSliceAbs.AddItem(m_pTHostLoopDensityAbs[i] * fFactor);
                }
            }
        }

        //P4-3.3: local spatial sum -> global, normalized by the GLOBAL volume.
        const cuDoubleComplex cSumX = appGetCudaHelper()->ReduceComplex(m_pTmpLoopX, _HC_Ly * _HC_Lz * _HC_Lt);
        cuDoubleComplex cLoopAvgX = cSumX;
        appGlobalSum(cLoopAvgX);
        const DOUBLE fVolX = GlobalL(1) * GlobalL(2) * GlobalL(3);
        m_lstLoopX.AddItem(cuCdivf_cd_host(cLoopAvgX, fVolX));
    }

    if (m_bMeasureLoopY)
    {
#if _CLG_MULTI_GPU
        //I9: the Y-direction Polyakov chain is a per-site product along y;
        //with a split y direction each rank holds only a partial chain and the
        //loop average is silently wrong. Reject (same rule as loopZ below).
        if (NULL != appGetComm() && appGetComm()->GpuGrid()[1] > 1)
        {
            appCrucial(_T("CMeasurePolyakovXY: loopY is not supported on multi-GPU with a split y direction (cross-rank y chain not implemented). Rejected.\n"));
            return;
        }
#endif
        pAcceptGauge->PolyakovOnSpatialSite(m_pTmpLoopY, 1);
        if (m_bMeasureXSlice)
        {
            _ZeroSlice(m_pXDeviceLoopDensity, 0);
            if (m_bMeasureAbs)
            {
                _ZeroSlice(m_pXDeviceLoopDensityAbs, 0);
            }
        }
        if (m_bMeasureZSlice)
        {
            _ZeroSlice(m_pZDeviceLoopDensity, 2);
            if (m_bMeasureAbs)
            {
                _ZeroSlice(m_pZDeviceLoopDensityAbs, 2);
            }
        }
        if (m_bMeasureTSlice)
        {
            _ZeroSlice(m_pTDeviceLoopDensity, 3);
            if (m_bMeasureAbs)
            {
                _ZeroSlice(m_pTDeviceLoopDensityAbs, 3);
            }
        }
        preparethread_Sxzt;
        _LAUNCH_KERNEL(_kernelPolyakovTraceOfSite_Y, block3dxzt, threads3dxzt, m_pTmpLoopY,
            m_bMeasureXSlice ? m_pXDeviceLoopDensity : NULL,
            m_bMeasureZSlice ? m_pZDeviceLoopDensity : NULL,
            m_bMeasureTSlice ? m_pTDeviceLoopDensity : NULL,
            (m_bMeasureXSlice && m_bMeasureAbs) ? m_pXDeviceLoopDensityAbs : NULL,
            (m_bMeasureZSlice && m_bMeasureAbs) ? m_pZDeviceLoopDensityAbs : NULL,
            (m_bMeasureTSlice && m_bMeasureAbs) ? m_pTDeviceLoopDensityAbs : NULL);

        if (m_bMeasureXSlice)
        {
            //P4-3.3: see main-section slices.
            const DOUBLE fFactor = 1.0 / static_cast<DOUBLE>(GlobalL(2) * GlobalL(3));
            checkCudaErrors(cudaMemcpy(m_pXHostLoopDensity, m_pXDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lx, cudaMemcpyDeviceToHost));
            if (m_bMeasureAbs)
            {
                checkCudaErrors(cudaMemcpy(m_pXHostLoopDensityAbs, m_pXDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lx, cudaMemcpyDeviceToHost));
            }
#if _CLG_MULTI_GPU
            if (NULL != appGetComm())
            {
                appGlobalSum(m_pXHostLoopDensity, _HC_Lx);
                if (m_bMeasureAbs)
                {
                    GlobalSumRealArray(m_pXHostLoopDensityAbs, _HC_Lx);
                }
            }
#endif
            for (UINT i = 0; i < _HC_Lx; ++i)
            {
                m_lstPY_XSlice.AddItem(cuCmulf_cd(m_pXHostLoopDensity[i], fFactor));
                if (m_bMeasureAbs)
                {
                    m_lstPY_XSliceAbs.AddItem(m_pXHostLoopDensityAbs[i] * fFactor);
                }
            }
        }
        if (m_bMeasureZSlice)
        {
#if _CLG_MULTI_GPU
            if (NULL != appGetComm() && appGetComm()->GpuGrid()[2] > 1)
            {
                appCrucial(_T("CMeasurePolyakovXY: loopY Z-slice is not supported on multi-GPU with a split z direction. Rejected.\n"));
                return;
            }
#endif
            const DOUBLE fFactor = 1.0 / static_cast<DOUBLE>(GlobalL(0) * GlobalL(3));
            checkCudaErrors(cudaMemcpy(m_pZHostLoopDensity, m_pZDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lz, cudaMemcpyDeviceToHost));
            if (m_bMeasureAbs)
            {
                checkCudaErrors(cudaMemcpy(m_pZHostLoopDensityAbs, m_pZDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lz, cudaMemcpyDeviceToHost));
            }
#if _CLG_MULTI_GPU
            if (NULL != appGetComm())
            {
                appGlobalSum(m_pZHostLoopDensity, _HC_Lz);
                if (m_bMeasureAbs)
                {
                    GlobalSumRealArray(m_pZHostLoopDensityAbs, _HC_Lz);
                }
            }
#endif
            for (UINT i = 0; i < _HC_Lz; ++i)
            {
                m_lstPY_ZSlice.AddItem(cuCmulf_cd(m_pZHostLoopDensity[i], fFactor));
                if (m_bMeasureAbs)
                {
                    m_lstPY_ZSliceAbs.AddItem(m_pZHostLoopDensityAbs[i] * fFactor);
                }
            }
        }
        if (m_bMeasureTSlice)
        {
            //P4-3.3: see main-section slices.
            const DOUBLE fFactor = 1.0 / static_cast<DOUBLE>(GlobalL(0) * GlobalL(2));
            checkCudaErrors(cudaMemcpy(m_pTHostLoopDensity, m_pTDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lt, cudaMemcpyDeviceToHost));
            if (m_bMeasureAbs)
            {
                checkCudaErrors(cudaMemcpy(m_pTHostLoopDensityAbs, m_pTDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lt, cudaMemcpyDeviceToHost));
            }
#if _CLG_MULTI_GPU
            if (NULL != appGetComm())
            {
                appGlobalSum(m_pTHostLoopDensity, _HC_Lt);
                if (m_bMeasureAbs)
                {
                    GlobalSumRealArray(m_pTHostLoopDensityAbs, _HC_Lt);
                }
            }
#endif
            for (UINT i = 0; i < _HC_Lt; ++i)
            {
                m_lstPY_TSlice.AddItem(cuCmulf_cd(m_pTHostLoopDensity[i], fFactor));
                if (m_bMeasureAbs)
                {
                    m_lstPY_TSliceAbs.AddItem(m_pTHostLoopDensityAbs[i] * fFactor);
                }
            }
        }

        //P4-3.3: local spatial sum -> global, normalized by the GLOBAL volume.
        const cuDoubleComplex cSumY = appGetCudaHelper()->ReduceComplex(m_pTmpLoopY, _HC_Lx * _HC_Lz * _HC_Lt);
        cuDoubleComplex cLoopAvgY = cSumY;
        appGlobalSum(cLoopAvgY);
        const DOUBLE fVolY = GlobalL(0) * GlobalL(2) * GlobalL(3);
        m_lstLoopY.AddItem(cuCdivf_cd_host(cLoopAvgY, fVolY));
    }

    if (m_bMeasureLoopZ)
    {
#if _CLG_MULTI_GPU
        //I9: the Z-direction Polyakov chain is a per-site product along z;
        //with a split z direction each rank holds only a partial chain, so the
        //loop average AND the R-distribution input are silently wrong. Reject
        //here so both branches below (distribution and plain average) are
        //covered (cross-rank z chain not implemented).
        if (NULL != appGetComm() && appGetComm()->GpuGrid()[2] > 1)
        {
            appCrucial(_T("CMeasurePolyakovXY: loopZ is not supported on multi-GPU with a split z direction (cross-rank z chain not implemented). Rejected.\n"));
            return;
        }
#endif
        pAcceptGauge->PolyakovOnSpatialSite(m_pTmpLoopZ, 2);
        if (m_bMeasureDistribution)
        {
            _ZeroXYPlane(m_pXYDeviceLoopDensity);
            if (m_bMeasureAbs)
            {
                _ZeroXYPlane(m_pXYDeviceLoopDensityAbs);
            }
        }
        if (m_bMeasureXSlice)
        {
            _ZeroSlice(m_pXDeviceLoopDensity, 0);
            if (m_bMeasureAbs)
            {
                _ZeroSlice(m_pXDeviceLoopDensityAbs, 0);
            }
        }
        if (m_bMeasureYSlice)
        {
            _ZeroSlice(m_pYDeviceLoopDensity, 1);
            if (m_bMeasureAbs)
            {
                _ZeroSlice(m_pYDeviceLoopDensityAbs, 1);
            }
        }
        if (m_bMeasureTSlice)
        {
            _ZeroSlice(m_pTDeviceLoopDensity, 3);
            if (m_bMeasureAbs)
            {
                _ZeroSlice(m_pTDeviceLoopDensityAbs, 3);
            }
        }
        preparethread_Sxyt;
        _LAUNCH_KERNEL(_kernelPolyakovTraceOfSite_Z, block3dxyt, threads3dxyt, m_pTmpLoopZ,
            m_bMeasureDistribution ? m_pXYDeviceLoopDensity : NULL,
            m_bMeasureXSlice ? m_pXDeviceLoopDensity : NULL,
            m_bMeasureYSlice ? m_pYDeviceLoopDensity : NULL,
            m_bMeasureTSlice ? m_pTDeviceLoopDensity : NULL,
            (m_bMeasureDistribution && m_bMeasureAbs) ? m_pXYDeviceLoopDensityAbs : NULL,
            (m_bMeasureXSlice && m_bMeasureAbs) ? m_pXDeviceLoopDensityAbs : NULL,
            (m_bMeasureYSlice && m_bMeasureAbs) ? m_pYDeviceLoopDensityAbs : NULL,
            (m_bMeasureTSlice && m_bMeasureAbs) ? m_pTDeviceLoopDensityAbs : NULL);

        if (m_bMeasureDistribution)
        {
            checkCudaErrors(cudaMemcpy(m_pXYHostLoopDensity, m_pXYDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
            if (m_bMeasureAbs)
            {
                checkCudaErrors(cudaMemcpy(m_pXYHostLoopDensityAbs, m_pXYDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
            }
#if _CLG_MULTI_GPU
            //P4-3.3: the XY density holds per-rank partial sums over the LOCAL
            //extent; sum across ranks so the transform sees the GLOBAL density.
            //Requires x/y NOT split (each rank then covers the same (x,y) index
            //set, so per-(x,y) values can be summed); a split x/y would need a
            //global-index gather, out of P4-3.3 scope.
            if (NULL != appGetComm())
            {
                if (appGetComm()->GpuGrid()[0] > 1 || appGetComm()->GpuGrid()[1] > 1)
                {
                    appCrucial(_T("CMeasurePolyakovXY: the XY distribution is not supported on multi-GPU with a split x/y direction. Rejected.\n"));
                    return;
                }
                appGlobalSum(m_pXYHostLoopDensity, _HC_Lx * _HC_Ly);
                if (m_bMeasureAbs)
                {
                    GlobalSumRealArray(m_pXYHostLoopDensityAbs, _HC_Lx * _HC_Ly);
                }
                checkCudaErrors(cudaMemcpy(m_pXYDeviceLoopDensity, m_pXYHostLoopDensity, sizeof(cuDoubleComplex) * _HC_Lx * _HC_Ly, cudaMemcpyHostToDevice));
                if (m_bMeasureAbs)
                {
                    checkCudaErrors(cudaMemcpy(m_pXYDeviceLoopDensityAbs, m_pXYHostLoopDensityAbs, sizeof(DOUBLE) * _HC_Lx * _HC_Ly, cudaMemcpyHostToDevice));
                }
            }
#endif

            TransformFromXYDataToRDataOnce(
                m_bShiftCenter,
                m_pXYDeviceLoopDensity,
                m_pDistributionR,
                m_pDistributionP,
                m_pHostDistributionR,
                m_pHostDistributionP,
                m_uiMaxR,
                m_uiEdgeR,
                FALSE,
                GetGaugeFieldIdSingleField(),
                m_lstPZ,
                &m_lstLoopInnerZ,
                m_lstLoopZ,
                m_lstR,
                m_uiConfigurationCount,
                F(1.0) / static_cast<Real>(GlobalL(3))
            );

            if (m_bMeasureAbs)
            {
                TransformFromXYDataToRDataOnce(
                    m_bShiftCenter,
                    m_pXYDeviceLoopDensityAbs,
                    m_pDistributionR,
                    m_pDistributionPAbs,
                    m_pHostDistributionR,
                    m_pHostDistributionPAbs,
                    m_uiMaxR,
                    m_uiEdgeR,
                    FALSE,
                    GetGaugeFieldIdSingleField(),
                    m_lstPZAbs,
                    &m_lstLoopAbsInnerZ,
                    m_lstLoopAbsZ,
                    m_lstR,
                    m_uiConfigurationCount,
                    F(1.0) / static_cast<Real>(GlobalL(3))
                );
            }
        }
        else
        {
            //P4-3.3: local spatial sum -> global, normalized by the GLOBAL volume.
            //(The split-z rejection now lives at the top of the loopZ section,
            //covering both branches.)
            const cuDoubleComplex cSumZ = appGetCudaHelper()->ReduceComplex(m_pTmpLoopZ, _HC_Lx * _HC_Ly * _HC_Lt);
            cuDoubleComplex cLoopAvgZ = cSumZ;
            appGlobalSum(cLoopAvgZ);
            const DOUBLE fVolZ = GlobalL(0) * GlobalL(1) * GlobalL(3);
            m_lstLoopZ.AddItem(cuCdivf_cd_host(cLoopAvgZ, fVolZ));
        }

        if (m_bMeasureXSlice)
        {
            //P4-3.3: see main-section slices.
            const DOUBLE fFactor = 1.0 / static_cast<DOUBLE>(GlobalL(1) * GlobalL(3));
            checkCudaErrors(cudaMemcpy(m_pXHostLoopDensity, m_pXDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lx, cudaMemcpyDeviceToHost));
            if (m_bMeasureAbs)
            {
                checkCudaErrors(cudaMemcpy(m_pXHostLoopDensityAbs, m_pXDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lx, cudaMemcpyDeviceToHost));
            }
#if _CLG_MULTI_GPU
            if (NULL != appGetComm())
            {
                appGlobalSum(m_pXHostLoopDensity, _HC_Lx);
                if (m_bMeasureAbs)
                {
                    GlobalSumRealArray(m_pXHostLoopDensityAbs, _HC_Lx);
                }
            }
#endif
            for (UINT i = 0; i < _HC_Lx; ++i)
            {
                m_lstPZ_XSlice.AddItem(cuCmulf_cd(m_pXHostLoopDensity[i], fFactor));
                if (m_bMeasureAbs)
                {
                    m_lstPZ_XSliceAbs.AddItem(m_pXHostLoopDensityAbs[i] * fFactor);
                }
            }
        }
        if (m_bMeasureYSlice)
        {
            //P4-3.3: see main-section slices.
            const DOUBLE fFactor = 1.0 / static_cast<DOUBLE>(GlobalL(0) * GlobalL(3));
            checkCudaErrors(cudaMemcpy(m_pYHostLoopDensity, m_pYDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Ly, cudaMemcpyDeviceToHost));
            if (m_bMeasureAbs)
            {
                checkCudaErrors(cudaMemcpy(m_pYHostLoopDensityAbs, m_pYDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Ly, cudaMemcpyDeviceToHost));
            }
#if _CLG_MULTI_GPU
            if (NULL != appGetComm())
            {
                appGlobalSum(m_pYHostLoopDensity, _HC_Ly);
                if (m_bMeasureAbs)
                {
                    GlobalSumRealArray(m_pYHostLoopDensityAbs, _HC_Ly);
                }
            }
#endif
            for (UINT i = 0; i < _HC_Ly; ++i)
            {
                m_lstPZ_YSlice.AddItem(cuCmulf_cd(m_pYHostLoopDensity[i], fFactor));
                if (m_bMeasureAbs)
                {
                    m_lstPZ_YSliceAbs.AddItem(m_pYHostLoopDensityAbs[i] * fFactor);
                }
            }
        }
        if (m_bMeasureTSlice)
        {
            //P4-3.3: see main-section slices.
            const DOUBLE fFactor = 1.0 / static_cast<DOUBLE>(GlobalL(0) * GlobalL(1));
            checkCudaErrors(cudaMemcpy(m_pTHostLoopDensity, m_pTDeviceLoopDensity, sizeof(cuDoubleComplex) * _HC_Lt, cudaMemcpyDeviceToHost));
            if (m_bMeasureAbs)
            {
                checkCudaErrors(cudaMemcpy(m_pTHostLoopDensityAbs, m_pTDeviceLoopDensityAbs, sizeof(DOUBLE) * _HC_Lt, cudaMemcpyDeviceToHost));
            }
#if _CLG_MULTI_GPU
            if (NULL != appGetComm())
            {
                appGlobalSum(m_pTHostLoopDensity, _HC_Lt);
                if (m_bMeasureAbs)
                {
                    GlobalSumRealArray(m_pTHostLoopDensityAbs, _HC_Lt);
                }
            }
#endif
            for (UINT i = 0; i < _HC_Lt; ++i)
            {
                m_lstPZ_TSlice.AddItem(cuCmulf_cd(m_pTHostLoopDensity[i], fFactor));
                if (m_bMeasureAbs)
                {
                    m_lstPZ_TSliceAbs.AddItem(m_pTHostLoopDensityAbs[i] * fFactor);
                }
            }
        }
    }

    if (m_bShowResult)
    {
        appDetailed(_T("\n\n ==================== Polyakov Loop (%d con)============================ \n\n"), m_uiConfigurationCount);
    }

    if (m_bShowResult)
    {
        appPushLogDate(FALSE);
        appGeneral(_T("Loop is "));
        LogGeneralComplex(m_lstLoop[m_lstLoop.GetCount() - 1]);
        if (m_bMeasureAbs)
        {
            appGeneral(_T(" Abs is %f\n"), m_lstLoopAbs[m_lstLoopAbs.GetCount() - 1]);
        }
        appPopLogDate();
    }

    //if (m_bShowResult)
    //{
    //    for (UINT i = 1; i < _HC_Lx; ++i)
    //    {
    //        appDetailed(_T("{"));
    //        for (UINT j = 1; j < _HC_Ly; ++j)
    //        {
    //            appDetailed(_T("%1.12f %s %1.12f I%s"),
    //                m_pXYHostLoopDensity[i * _HC_Ly + j].x,
    //                m_pXYHostLoopDensity[i * _HC_Ly + j].y < F(0.0) ? _T("-") : _T("+"),
    //                appAbs(m_pXYHostLoopDensity[i * _HC_Ly + j].y),
    //                (j == _HC_Ly - 1) ? _T("},\n") : _T(",   ")
    //            );
    //        }
    //    }
    //}

    if (m_bShowResult)
    {
        appGeneral(_T("\n"));
    }

    if (m_bShowResult)
    {
        appDetailed(_T("\n=====================================================\n"), m_uiConfigurationCount);
    }

    if (NULL != m_pOwner)
    {
        const cuDoubleComplex& cLoopT = m_lstLoop[m_lstLoop.GetCount() - 1];
        m_pOwner->AddOneConfigurationResult(this, _T("PolyakovT"),
            _make_cuComplex(static_cast<Real>(cLoopT.x), static_cast<Real>(cLoopT.y)));

        if (m_bMeasureDistribution)
        {
            const cuDoubleComplex& cLoopInner = m_lstLoopInner[m_lstLoopInner.GetCount() - 1];
            m_pOwner->AddOneConfigurationResult(this, _T("PolyakovInner"),
                _make_cuComplex(static_cast<Real>(cLoopInner.x), static_cast<Real>(cLoopInner.y)));

            const INT iRCount = m_lstR.Num();
            TArray<CLGComplex> thisConfigurationP;
            for (INT i = 0; i < iRCount; ++i)
            {
                const cuDoubleComplex& c = m_lstP[m_lstP.Num() - iRCount + i];
                thisConfigurationP.AddItem(_make_cuComplex(static_cast<Real>(c.x), static_cast<Real>(c.y)));
            }
            m_pOwner->AddOneConfigurationResult(this, _T("PolyakovOverR"), thisConfigurationP);

            if (m_bMeasureAbs)
            {
                m_pOwner->AddOneConfigurationResult(this, _T("PolyakovTAbs"), m_lstLoopAbs[m_lstLoopAbs.GetCount() - 1]);
                m_pOwner->AddOneConfigurationResult(this, _T("PolyakovInnerAbs"), m_lstLoopAbsInner[m_lstLoopAbsInner.GetCount() - 1]);

                TArray<DOUBLE> thisConfigurationPAbs;
                for (INT i = 0; i < iRCount; ++i)
                {
                    thisConfigurationPAbs.AddItem(m_lstPAbs[m_lstPAbs.Num() - iRCount + i]);
                }
                m_pOwner->AddOneConfigurationResult(this, _T("PolyakovOverRAbs"), thisConfigurationPAbs);
            }
        }

        if (m_bMeasureLoopX)
        {
            const cuDoubleComplex& cLoopX = m_lstLoopX[m_lstLoopX.GetCount() - 1];
            m_pOwner->AddOneConfigurationResult(this, _T("PolyakovX"),
                _make_cuComplex(static_cast<Real>(cLoopX.x), static_cast<Real>(cLoopX.y)));
        }
        if (m_bMeasureLoopY)
        {
            const cuDoubleComplex& cLoopY = m_lstLoopY[m_lstLoopY.GetCount() - 1];
            m_pOwner->AddOneConfigurationResult(this, _T("PolyakovY"),
                _make_cuComplex(static_cast<Real>(cLoopY.x), static_cast<Real>(cLoopY.y)));
        }
        if (m_bMeasureLoopZ)
        {
            const cuDoubleComplex& cLoopZ = m_lstLoopZ[m_lstLoopZ.GetCount() - 1];
            m_pOwner->AddOneConfigurationResult(this, _T("PolyakovZ"),
                _make_cuComplex(static_cast<Real>(cLoopZ.x), static_cast<Real>(cLoopZ.y)));

            if (m_bMeasureDistribution)
            {
                const cuDoubleComplex& cLoopInnerZ = m_lstLoopInnerZ[m_lstLoopInnerZ.GetCount() - 1];
                m_pOwner->AddOneConfigurationResult(this, _T("PolyakovZInner"),
                    _make_cuComplex(static_cast<Real>(cLoopInnerZ.x), static_cast<Real>(cLoopInnerZ.y)));

                const INT iRCount = m_lstR.Num();
                TArray<CLGComplex> thisConfigurationPZ;
                for (INT i = 0; i < iRCount; ++i)
                {
                    const cuDoubleComplex& c = m_lstPZ[m_lstPZ.Num() - iRCount + i];
                    thisConfigurationPZ.AddItem(_make_cuComplex(static_cast<Real>(c.x), static_cast<Real>(c.y)));
                }
                m_pOwner->AddOneConfigurationResult(this, _T("PolyakovZOverR"), thisConfigurationPZ);

                if (m_bMeasureAbs)
                {
                    m_pOwner->AddOneConfigurationResult(this, _T("PolyakovZAbs"), m_lstLoopAbsZ[m_lstLoopAbsZ.GetCount() - 1]);
                    m_pOwner->AddOneConfigurationResult(this, _T("PolyakovZInnerAbs"), m_lstLoopAbsInnerZ[m_lstLoopAbsInnerZ.GetCount() - 1]);

                    TArray<DOUBLE> thisConfigurationPZAbs;
                    for (INT i = 0; i < iRCount; ++i)
                    {
                        thisConfigurationPZAbs.AddItem(m_lstPZAbs[m_lstPZAbs.Num() - iRCount + i]);
                    }
                    m_pOwner->AddOneConfigurationResult(this, _T("PolyakovZOverRAbs"), thisConfigurationPZAbs);
                }
            }
        }
    }

    // Expose the per-slice profiles of the Polyakov loops to the dictionary
    if (NULL != m_pOwner)
    {
        if (m_bMeasureXSlice) { _AddSliceToDictionary(m_pOwner, this, _T("PolyakovXSlice"), m_lstP_XSlice, m_lstP_XSliceAbs, _HC_Lx, m_bMeasureAbs); }
        if (m_bMeasureYSlice) { _AddSliceToDictionary(m_pOwner, this, _T("PolyakovYSlice"), m_lstP_YSlice, m_lstP_YSliceAbs, _HC_Ly, m_bMeasureAbs); }
        if (m_bMeasureZSlice) { _AddSliceToDictionary(m_pOwner, this, _T("PolyakovZSlice"), m_lstP_ZSlice, m_lstP_ZSliceAbs, _HC_Lz, m_bMeasureAbs); }
        if (m_bMeasureLoopX)
        {
            if (m_bMeasureYSlice) { _AddSliceToDictionary(m_pOwner, this, _T("PolyakovX_YSlice"), m_lstPX_YSlice, m_lstPX_YSliceAbs, _HC_Ly, m_bMeasureAbs); }
            if (m_bMeasureZSlice) { _AddSliceToDictionary(m_pOwner, this, _T("PolyakovX_ZSlice"), m_lstPX_ZSlice, m_lstPX_ZSliceAbs, _HC_Lz, m_bMeasureAbs); }
            if (m_bMeasureTSlice) { _AddSliceToDictionary(m_pOwner, this, _T("PolyakovX_TSlice"), m_lstPX_TSlice, m_lstPX_TSliceAbs, _HC_Lt, m_bMeasureAbs); }
        }
        if (m_bMeasureLoopY)
        {
            if (m_bMeasureXSlice) { _AddSliceToDictionary(m_pOwner, this, _T("PolyakovY_XSlice"), m_lstPY_XSlice, m_lstPY_XSliceAbs, _HC_Lx, m_bMeasureAbs); }
            if (m_bMeasureZSlice) { _AddSliceToDictionary(m_pOwner, this, _T("PolyakovY_ZSlice"), m_lstPY_ZSlice, m_lstPY_ZSliceAbs, _HC_Lz, m_bMeasureAbs); }
            if (m_bMeasureTSlice) { _AddSliceToDictionary(m_pOwner, this, _T("PolyakovY_TSlice"), m_lstPY_TSlice, m_lstPY_TSliceAbs, _HC_Lt, m_bMeasureAbs); }
        }
        if (m_bMeasureLoopZ)
        {
            if (m_bMeasureXSlice) { _AddSliceToDictionary(m_pOwner, this, _T("PolyakovZ_XSlice"), m_lstPZ_XSlice, m_lstPZ_XSliceAbs, _HC_Lx, m_bMeasureAbs); }
            if (m_bMeasureYSlice) { _AddSliceToDictionary(m_pOwner, this, _T("PolyakovZ_YSlice"), m_lstPZ_YSlice, m_lstPZ_YSliceAbs, _HC_Ly, m_bMeasureAbs); }
            if (m_bMeasureTSlice) { _AddSliceToDictionary(m_pOwner, this, _T("PolyakovZ_TSlice"), m_lstPZ_TSlice, m_lstPZ_TSliceAbs, _HC_Lt, m_bMeasureAbs); }
        }
    }

    ++m_uiConfigurationCount;
}

void CMeasurePolyakovXY::Report()
{
    //UINT uiHalf = (_HC_Lx + 1) / 2;
    appAssert(m_uiConfigurationCount == static_cast<UINT>(m_lstLoop.Num()));
    //appAssert(static_cast<UINT>(m_uiConfigurationCount * uiHalf)
    //    == static_cast<UINT>(m_lstLoopDensity.Num()));

    appPushLogDate(FALSE);
    CLGComplex tmpChargeSum = _make_cuComplex(F(0.0), F(0.0));
    //m_lstAverageLoopDensity.RemoveAll();

    appGeneral(_T("\n\n==========================================================================\n"));
    appGeneral(_T("==================== Polyakov Loop (%d con)============================\n"), m_uiConfigurationCount);

    appGeneral(_T("\n ----------- Loop ------------- \n"));

    appGeneral(_T("{"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        tmpChargeSum.x += static_cast<Real>(m_lstLoop[i].x);
        tmpChargeSum.y += static_cast<Real>(m_lstLoop[i].y);
        LogGeneralComplex(m_lstLoop[i]);
    }
    appGeneral(_T("}\n"));

    tmpChargeSum.x = tmpChargeSum.x / m_uiConfigurationCount;
    tmpChargeSum.y = tmpChargeSum.y / m_uiConfigurationCount;
    m_cAverageLoop = tmpChargeSum;
    appGeneral(_T("\n ----------- average Loop |<P>| = %2.12f arg(P) = %2.12f ------------- \n"), _cuCabsf(tmpChargeSum), __cuCargf(tmpChargeSum));

    //appGeneral(_T("\n ----------- Loop density ------------- \n"));

    //appGeneral(_T("{\n"));
    //for (UINT k = 0; k < m_uiConfigurationCount; ++k)
    //{
    //    //appGeneral(_T("{"));
    //    for (UINT i = 0; i < uiHalf; ++i)
    //    {
    //        //LogGeneralComplex(m_lstLoopDensity[k * uiHalf + i]);

    //        //if (0 == k)
    //        //{
    //        //    m_lstAverageLoopDensity.AddItem(m_lstLoopDensity[k * uiHalf + i]);
    //        //}
    //        //else
    //        //{
    //        //    m_lstAverageLoopDensity[i] = _cuCaddf(m_lstAverageLoopDensity[i], m_lstLoopDensity[k * uiHalf + i]);
    //        //}

    //        if (k == m_uiConfigurationCount - 1)
    //        {
    //            m_lstAverageLoopDensity[i].x = m_lstAverageLoopDensity[i].x / m_uiConfigurationCount;
    //            m_lstAverageLoopDensity[i].y = m_lstAverageLoopDensity[i].y / m_uiConfigurationCount;
    //        }
    //    }
    //    //appGeneral(_T("}\n"));
    //}
    //appGeneral(_T("}\n"));

    appGeneral(_T("\n==========================================================================\n"));
    appGeneral(_T("==========================================================================\n\n"));
    appPopLogDate();
}

void CMeasurePolyakovXY::Reset()
{
    CMeasure::Reset();
    m_lstLoop.RemoveAll();
    m_lstLoopX.RemoveAll();
    m_lstLoopY.RemoveAll();
    m_lstLoopZ.RemoveAll();
    m_lstLoopInner.RemoveAll();
    m_lstLoopInnerZ.RemoveAll();
    m_lstLoopAbs.RemoveAll();
    m_lstLoopAbsZ.RemoveAll();
    m_lstLoopAbsInner.RemoveAll();
    m_lstLoopAbsInnerZ.RemoveAll();

    m_lstR.RemoveAll();
    m_lstP.RemoveAll();
    m_lstPAbs.RemoveAll();
    m_lstPZ.RemoveAll();
    m_lstPZAbs.RemoveAll();

    m_lstP_XSlice.RemoveAll();
    m_lstP_YSlice.RemoveAll();
    m_lstP_ZSlice.RemoveAll();
    m_lstP_XSliceAbs.RemoveAll();
    m_lstP_YSliceAbs.RemoveAll();
    m_lstP_ZSliceAbs.RemoveAll();

    m_lstPX_YSlice.RemoveAll();
    m_lstPX_ZSlice.RemoveAll();
    m_lstPX_TSlice.RemoveAll();
    m_lstPX_YSliceAbs.RemoveAll();
    m_lstPX_ZSliceAbs.RemoveAll();
    m_lstPX_TSliceAbs.RemoveAll();

    m_lstPY_XSlice.RemoveAll();
    m_lstPY_ZSlice.RemoveAll();
    m_lstPY_TSlice.RemoveAll();
    m_lstPY_XSliceAbs.RemoveAll();
    m_lstPY_ZSliceAbs.RemoveAll();
    m_lstPY_TSliceAbs.RemoveAll();

    m_lstPZ_XSlice.RemoveAll();
    m_lstPZ_YSlice.RemoveAll();
    m_lstPZ_TSlice.RemoveAll();
    m_lstPZ_XSliceAbs.RemoveAll();
    m_lstPZ_YSliceAbs.RemoveAll();
    m_lstPZ_TSliceAbs.RemoveAll();
}

void CMeasurePolyakovXY::Export(const CCString& sCSV, UINT iStartN, UINT iEndN, const CCString& sOName, UINT uiOmega, UINT iListStart) const
{
    appAssert(static_cast<INT>(iEndN - iStartN + 1) == m_lstLoop.Num());
    CCString sFileNameWrite;
    if (m_bMeasureDistribution)
    {
        appAssert(static_cast<INT>(iEndN - iStartN + 1) * m_lstR.Num() == m_lstP.Num());
        TArray<Real> lstR;
        if (uiOmega == iListStart)
        {
            for (INT i = 0; i < m_lstR.Num(); ++i)
            {
                if (m_bShiftCenter)
                {
                    lstR.AddItem(F(0.5) * _hostsqrt(static_cast<Real>(m_lstR[i])));
                }
                else
                {
                    lstR.AddItem(_hostsqrt(static_cast<Real>(m_lstR[i])));
                }
            }
            sFileNameWrite.Format(_T("%s_polyakov_R.csv"), sCSV.c_str());
            WriteRealArray(sFileNameWrite, lstR);
        }
    }
    TArray<cuDoubleComplex> polyIn;
    TArray<DOUBLE> polyInAbs;
    TArray<cuDoubleComplex> polyOut;
    TArray<DOUBLE> polyOutAbs;
    TArray<TArray<cuDoubleComplex>> polyakovSlice;
    TArray<TArray<DOUBLE>> polyakovSliceAbs;
    for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
    {
        if (m_bMeasureDistribution)
        {
            polyIn.AddItem(m_lstLoopInner[j]);
            if (m_bMeasureAbs)
            {
                polyInAbs.AddItem(m_lstLoopAbsInner[j]);
            }
        }

        polyOut.AddItem(m_lstLoop[j]);
        if (m_bMeasureDistribution)
        {
            if (m_bMeasureAbs)
            {
                //Note here, even Abs was ture, if does not measure distribution, the abs was not measured
                polyOutAbs.AddItem(m_lstLoopAbs[j]);
            }

            TArray<cuDoubleComplex> thisConfigurationSlice;
            TArray<DOUBLE> thisConfigurationSliceAbs;
            for (INT i = 0; i < m_lstR.Num(); ++i)
            {
                thisConfigurationSlice.AddItem(m_lstP[j * m_lstR.Num() + i]);
                if (m_bMeasureAbs)
                {
                    thisConfigurationSliceAbs.AddItem(m_lstPAbs[j * m_lstR.Num() + i]);
                }
            }
            polyakovSlice.AddItem(thisConfigurationSlice);
            if (m_bMeasureAbs)
            {
                polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
            }
        }
    }
    if (m_bMeasureDistribution)
    {
        sFileNameWrite.Format(_T("%s_%s_polyakov_In.csv"), sCSV.c_str(), sOName.c_str());
        WriteComplexArray(sFileNameWrite, polyIn);
        sFileNameWrite.Format(_T("%s_%s_polyakov_OverR.csv"), sCSV.c_str(), sOName.c_str());
        WriteComplexArray2(sFileNameWrite, polyakovSlice);
        if (m_bMeasureAbs)
        {
            sFileNameWrite.Format(_T("%s_%s_polyakov_InAbs.csv"), sCSV.c_str(), sOName.c_str());
            WriteRealArray(sFileNameWrite, polyInAbs);
            sFileNameWrite.Format(_T("%s_%s_polyakov_OverRAbs.csv"), sCSV.c_str(), sOName.c_str());
            WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
            sFileNameWrite.Format(_T("%s_%s_polyakov_Abs.csv"), sCSV.c_str(), sOName.c_str());
            WriteRealArray(sFileNameWrite, polyOutAbs);
        }
    }
    sFileNameWrite.Format(_T("%s_%s_polyakov.csv"), sCSV.c_str(), sOName.c_str());
    WriteComplexArray(sFileNameWrite, polyOut);

    if (m_bMeasureXSlice)
    {
        polyakovSlice.RemoveAll();
        polyakovSliceAbs.RemoveAll();
        for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
        {
            TArray<cuDoubleComplex> thisConfigurationSlice;
            TArray<DOUBLE> thisConfigurationSliceAbs;
            for (UINT i = 0; i < _HC_Lx; ++i)
            {
                thisConfigurationSlice.AddItem(m_lstP_XSlice[j * _HC_Lx + i]);
                if (m_bMeasureAbs)
                {
                    thisConfigurationSliceAbs.AddItem(m_lstP_XSliceAbs[j * _HC_Lx + i]);
                }
            }
            polyakovSlice.AddItem(thisConfigurationSlice);
            if (m_bMeasureAbs)
            {
                polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
            }
        }
        sFileNameWrite.Format(_T("%s_%s_polyakov_XSlice.csv"), sCSV.c_str(), sOName.c_str());
        WriteComplexArray2(sFileNameWrite, polyakovSlice);
        if (m_bMeasureAbs)
        {
            sFileNameWrite.Format(_T("%s_%s_polyakov_XSliceAbs.csv"), sCSV.c_str(), sOName.c_str());
            WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
        }
    }
    if (m_bMeasureYSlice)
    {
        polyakovSlice.RemoveAll();
        polyakovSliceAbs.RemoveAll();
        for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
        {
            TArray<cuDoubleComplex> thisConfigurationSlice;
            TArray<DOUBLE> thisConfigurationSliceAbs;
            for (UINT i = 0; i < _HC_Ly; ++i)
            {
                thisConfigurationSlice.AddItem(m_lstP_YSlice[j * _HC_Ly + i]);
                if (m_bMeasureAbs)
                {
                    thisConfigurationSliceAbs.AddItem(m_lstP_YSliceAbs[j * _HC_Ly + i]);
                }
            }
            polyakovSlice.AddItem(thisConfigurationSlice);
            if (m_bMeasureAbs)
            {
                polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
            }
        }
        sFileNameWrite.Format(_T("%s_%s_polyakov_YSlice.csv"), sCSV.c_str(), sOName.c_str());
        WriteComplexArray2(sFileNameWrite, polyakovSlice);
        if (m_bMeasureAbs)
        {
            sFileNameWrite.Format(_T("%s_%s_polyakov_YSliceAbs.csv"), sCSV.c_str(), sOName.c_str());
            WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
        }
    }
    if (m_bMeasureZSlice)
    {
        polyakovSlice.RemoveAll();
        polyakovSliceAbs.RemoveAll();
        for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
        {
            TArray<cuDoubleComplex> thisConfigurationSlice;
            TArray<DOUBLE> thisConfigurationSliceAbs;
            for (UINT i = 0; i < _HC_Lz; ++i)
            {
                thisConfigurationSlice.AddItem(m_lstP_ZSlice[j * _HC_Lz + i]);
                if (m_bMeasureAbs)
                {
                    thisConfigurationSliceAbs.AddItem(m_lstP_ZSliceAbs[j * _HC_Lz + i]);
                }
            }
            polyakovSlice.AddItem(thisConfigurationSlice);
            if (m_bMeasureAbs)
            {
                polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
            }
        }
        sFileNameWrite.Format(_T("%s_%s_polyakov_ZSlice.csv"), sCSV.c_str(), sOName.c_str());
        WriteComplexArray2(sFileNameWrite, polyakovSlice);
        if (m_bMeasureAbs)
        {
            sFileNameWrite.Format(_T("%s_%s_polyakov_ZSliceAbs.csv"), sCSV.c_str(), sOName.c_str());
            WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
        }
    }

    if (m_bMeasureLoopX)
    {
        polyOut.RemoveAll();
        for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
        {
            polyOut.AddItem(m_lstLoopX[j]);
        }
        sFileNameWrite.Format(_T("%s_%s_polyakovX.csv"), sCSV.c_str(), sOName.c_str());
        WriteComplexArray(sFileNameWrite, polyOut);

        if (m_bMeasureYSlice)
        {
            polyakovSlice.RemoveAll();
            if (m_bMeasureAbs)
            {
                polyakovSliceAbs.RemoveAll();
            }
            for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
            {
                TArray<cuDoubleComplex> thisConfigurationSlice;
                TArray<DOUBLE> thisConfigurationSliceAbs;
                for (UINT i = 0; i < _HC_Ly; ++i)
                {
                    thisConfigurationSlice.AddItem(m_lstPX_YSlice[j * _HC_Ly + i]);
                    if (m_bMeasureAbs)
                    {
                        thisConfigurationSliceAbs.AddItem(m_lstPX_YSliceAbs[j * _HC_Ly + i]);
                    }
                }
                polyakovSlice.AddItem(thisConfigurationSlice);
                if (m_bMeasureAbs)
                {
                    polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
                }
            }
            sFileNameWrite.Format(_T("%s_%s_polyakovX_YSlice.csv"), sCSV.c_str(), sOName.c_str());
            WriteComplexArray2(sFileNameWrite, polyakovSlice);
            if (m_bMeasureAbs)
            {
                sFileNameWrite.Format(_T("%s_%s_polyakovX_YSliceAbs.csv"), sCSV.c_str(), sOName.c_str());
                WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
            }
        }
        if (m_bMeasureZSlice)
        {
            polyakovSlice.RemoveAll();
            if (m_bMeasureAbs)
            {
                polyakovSliceAbs.RemoveAll();
            }
            for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
            {
                TArray<cuDoubleComplex> thisConfigurationSlice;
                TArray<DOUBLE> thisConfigurationSliceAbs;
                for (UINT i = 0; i < _HC_Lz; ++i)
                {
                    thisConfigurationSlice.AddItem(m_lstPX_ZSlice[j * _HC_Lz + i]);
                    if (m_bMeasureAbs)
                    {
                        thisConfigurationSliceAbs.AddItem(m_lstPX_ZSliceAbs[j * _HC_Lz + i]);
                    }
                }
                polyakovSlice.AddItem(thisConfigurationSlice);
                if (m_bMeasureAbs)
                {
                    polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
                }
            }
            sFileNameWrite.Format(_T("%s_%s_polyakovX_ZSlice.csv"), sCSV.c_str(), sOName.c_str());
            WriteComplexArray2(sFileNameWrite, polyakovSlice);
            if (m_bMeasureAbs)
            {
                sFileNameWrite.Format(_T("%s_%s_polyakovX_ZSliceAbs.csv"), sCSV.c_str(), sOName.c_str());
                WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
            }
        }
        if (m_bMeasureTSlice)
        {
            polyakovSlice.RemoveAll();
            if (m_bMeasureAbs)
            {
                polyakovSliceAbs.RemoveAll();
            }
            for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
            {
                TArray<cuDoubleComplex> thisConfigurationSlice;
                TArray<DOUBLE> thisConfigurationSliceAbs;
                for (UINT i = 0; i < _HC_Lt; ++i)
                {
                    thisConfigurationSlice.AddItem(m_lstPX_TSlice[j * _HC_Lt + i]);
                    if (m_bMeasureAbs)
                    {
                        thisConfigurationSliceAbs.AddItem(m_lstPX_TSliceAbs[j * _HC_Lt + i]);
                    }
                }
                polyakovSlice.AddItem(thisConfigurationSlice);
                if (m_bMeasureAbs)
                {
                    polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
                }
            }
            sFileNameWrite.Format(_T("%s_%s_polyakovX_TSlice.csv"), sCSV.c_str(), sOName.c_str());
            WriteComplexArray2(sFileNameWrite, polyakovSlice);
            if (m_bMeasureAbs)
            {
                sFileNameWrite.Format(_T("%s_%s_polyakovX_TSliceAbs.csv"), sCSV.c_str(), sOName.c_str());
                WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
            }
        }
    }

    if (m_bMeasureLoopY)
    {
        polyOut.RemoveAll();
        for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
        {
            polyOut.AddItem(m_lstLoopY[j]);
        }
        sFileNameWrite.Format(_T("%s_%s_polyakovY.csv"), sCSV.c_str(), sOName.c_str());
        WriteComplexArray(sFileNameWrite, polyOut);

        if (m_bMeasureXSlice)
        {
            polyakovSlice.RemoveAll();
            if (m_bMeasureAbs)
            {
                polyakovSliceAbs.RemoveAll();
            }
            for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
            {
                TArray<cuDoubleComplex> thisConfigurationSlice;
                TArray<DOUBLE> thisConfigurationSliceAbs;
                for (UINT i = 0; i < _HC_Lx; ++i)
                {
                    thisConfigurationSlice.AddItem(m_lstPY_XSlice[j * _HC_Lx + i]);
                    if (m_bMeasureAbs)
                    {
                        thisConfigurationSliceAbs.AddItem(m_lstPY_XSliceAbs[j * _HC_Lx + i]);
                    }
                }
                polyakovSlice.AddItem(thisConfigurationSlice);
                if (m_bMeasureAbs)
                {
                    polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
                }
            }
            sFileNameWrite.Format(_T("%s_%s_polyakovY_XSlice.csv"), sCSV.c_str(), sOName.c_str());
            WriteComplexArray2(sFileNameWrite, polyakovSlice);
            if (m_bMeasureAbs)
            {
                sFileNameWrite.Format(_T("%s_%s_polyakovY_XSliceAbs.csv"), sCSV.c_str(), sOName.c_str());
                WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
            }
        }
        if (m_bMeasureZSlice)
        {
            polyakovSlice.RemoveAll();
            if (m_bMeasureAbs)
            {
                polyakovSliceAbs.RemoveAll();
            }
            for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
            {
                TArray<cuDoubleComplex> thisConfigurationSlice;
                TArray<DOUBLE> thisConfigurationSliceAbs;
                for (UINT i = 0; i < _HC_Lz; ++i)
                {
                    thisConfigurationSlice.AddItem(m_lstPY_ZSlice[j * _HC_Lz + i]);
                    if (m_bMeasureAbs)
                    {
                        thisConfigurationSliceAbs.AddItem(m_lstPY_ZSliceAbs[j * _HC_Lz + i]);
                    }
                }
                polyakovSlice.AddItem(thisConfigurationSlice);
                if (m_bMeasureAbs)
                {
                    polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
                }
            }
            sFileNameWrite.Format(_T("%s_%s_polyakovY_ZSlice.csv"), sCSV.c_str(), sOName.c_str());
            WriteComplexArray2(sFileNameWrite, polyakovSlice);
            if (m_bMeasureAbs)
            {
                sFileNameWrite.Format(_T("%s_%s_polyakovY_ZSliceAbs.csv"), sCSV.c_str(), sOName.c_str());
                WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
            }
        }
        if (m_bMeasureTSlice)
        {
            polyakovSlice.RemoveAll();
            if (m_bMeasureAbs)
            {
                polyakovSliceAbs.RemoveAll();
            }
            for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
            {
                TArray<cuDoubleComplex> thisConfigurationSlice;
                TArray<DOUBLE> thisConfigurationSliceAbs;
                for (UINT i = 0; i < _HC_Lt; ++i)
                {
                    thisConfigurationSlice.AddItem(m_lstPY_TSlice[j * _HC_Lt + i]);
                    if (m_bMeasureAbs)
                    {
                        thisConfigurationSliceAbs.AddItem(m_lstPY_TSliceAbs[j * _HC_Lt + i]);
                    }
                }
                polyakovSlice.AddItem(thisConfigurationSlice);
                if (m_bMeasureAbs)
                {
                    polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
                }
            }
            sFileNameWrite.Format(_T("%s_%s_polyakovY_TSlice.csv"), sCSV.c_str(), sOName.c_str());
            WriteComplexArray2(sFileNameWrite, polyakovSlice);
            if (m_bMeasureAbs)
            {
                sFileNameWrite.Format(_T("%s_%s_polyakovY_TSliceAbs.csv"), sCSV.c_str(), sOName.c_str());
                WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
            }
        }
    }

    if (m_bMeasureLoopZ)
    {
        polyOut.RemoveAll();
        if (m_bMeasureDistribution)
        {
            polyIn.RemoveAll();
            polyakovSlice.RemoveAll();
            if (m_bMeasureAbs)
            {
                polyOutAbs.RemoveAll();
                polyInAbs.RemoveAll();
                polyakovSliceAbs.RemoveAll();
            }
        }

        for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
        {
            polyOut.AddItem(m_lstLoopZ[j]);
            if (m_bMeasureDistribution)
            {
                polyIn.AddItem(m_lstLoopInnerZ[j]);
                
                TArray<cuDoubleComplex> thisConfigurationSlice;
                TArray<DOUBLE> thisConfigurationSliceAbs;
                for (INT i = 0; i < m_lstR.Num(); ++i)
                {
                    thisConfigurationSlice.AddItem(m_lstPZ[j * m_lstR.Num() + i]);
                    if (m_bMeasureAbs)
                    {
                        thisConfigurationSliceAbs.AddItem(m_lstPZAbs[j * m_lstR.Num() + i]);
                    }
                }
                polyakovSlice.AddItem(thisConfigurationSlice);
                if (m_bMeasureAbs)
                {
                    polyInAbs.AddItem(m_lstLoopAbsInnerZ[j]);
                    polyOutAbs.AddItem(m_lstLoopAbsZ[j]);
                    polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
                }
            }
        }
        if (m_bMeasureDistribution)
        {
            sFileNameWrite.Format(_T("%s_%s_polyakovZ_In.csv"), sCSV.c_str(), sOName.c_str());
            WriteComplexArray(sFileNameWrite, polyIn);
            sFileNameWrite.Format(_T("%s_%s_polyakovZ_OverR.csv"), sCSV.c_str(), sOName.c_str());
            WriteComplexArray2(sFileNameWrite, polyakovSlice);
            if (m_bMeasureAbs)
            {
                sFileNameWrite.Format(_T("%s_%s_polyakovZ_InAbs.csv"), sCSV.c_str(), sOName.c_str());
                WriteRealArray(sFileNameWrite, polyInAbs);
                sFileNameWrite.Format(_T("%s_%s_polyakovZ_OverRAbs.csv"), sCSV.c_str(), sOName.c_str());
                WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
                sFileNameWrite.Format(_T("%s_%s_polyakovZ_Abs.csv"), sCSV.c_str(), sOName.c_str());
                WriteRealArray(sFileNameWrite, polyOutAbs);
            }
        }
        sFileNameWrite.Format(_T("%s_%s_polyakovZ.csv"), sCSV.c_str(), sOName.c_str());
        WriteComplexArray(sFileNameWrite, polyOut);

        if (m_bMeasureXSlice)
        {
            polyakovSlice.RemoveAll();
            if (m_bMeasureAbs)
            {
                polyakovSliceAbs.RemoveAll();
            }
            for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
            {
                TArray<cuDoubleComplex> thisConfigurationSlice;
                TArray<DOUBLE> thisConfigurationSliceAbs;
                for (UINT i = 0; i < _HC_Lx; ++i)
                {
                    thisConfigurationSlice.AddItem(m_lstPZ_XSlice[j * _HC_Lx + i]);
                    if (m_bMeasureAbs)
                    {
                        thisConfigurationSliceAbs.AddItem(m_lstPZ_XSliceAbs[j * _HC_Lx + i]);
                    }
                }
                polyakovSlice.AddItem(thisConfigurationSlice);
                if (m_bMeasureAbs)
                {
                    polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
                }
            }
            sFileNameWrite.Format(_T("%s_%s_polyakovZ_XSlice.csv"), sCSV.c_str(), sOName.c_str());
            WriteComplexArray2(sFileNameWrite, polyakovSlice);
            if (m_bMeasureAbs)
            {
                sFileNameWrite.Format(_T("%s_%s_polyakovZ_XSliceAbs.csv"), sCSV.c_str(), sOName.c_str());
                WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
            }
        }
        if (m_bMeasureYSlice)
        {
            polyakovSlice.RemoveAll();
            if (m_bMeasureAbs)
            {
                polyakovSliceAbs.RemoveAll();
            }
            for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
            {
                TArray<cuDoubleComplex> thisConfigurationSlice;
                TArray<DOUBLE> thisConfigurationSliceAbs;
                for (UINT i = 0; i < _HC_Ly; ++i)
                {
                    thisConfigurationSlice.AddItem(m_lstPZ_YSlice[j * _HC_Ly + i]);
                    if (m_bMeasureAbs)
                    {
                        thisConfigurationSliceAbs.AddItem(m_lstPZ_YSliceAbs[j * _HC_Ly + i]);
                    }
                }
                polyakovSlice.AddItem(thisConfigurationSlice);
                if (m_bMeasureAbs)
                {
                    polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
                }
            }
            sFileNameWrite.Format(_T("%s_%s_polyakovZ_YSlice.csv"), sCSV.c_str(), sOName.c_str());
            WriteComplexArray2(sFileNameWrite, polyakovSlice);
            if (m_bMeasureAbs)
            {
                sFileNameWrite.Format(_T("%s_%s_polyakovZ_YSliceAbs.csv"), sCSV.c_str(), sOName.c_str());
                WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
            }
        }
        if (m_bMeasureTSlice)
        {
            polyakovSlice.RemoveAll();
            if (m_bMeasureAbs)
            {
                polyakovSliceAbs.RemoveAll();
            }
            for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
            {
                TArray<cuDoubleComplex> thisConfigurationSlice;
                TArray<DOUBLE> thisConfigurationSliceAbs;
                for (UINT i = 0; i < _HC_Lt; ++i)
                {
                    thisConfigurationSlice.AddItem(m_lstPZ_TSlice[j * _HC_Lt + i]);
                    if (m_bMeasureAbs)
                    {
                        thisConfigurationSliceAbs.AddItem(m_lstPZ_TSliceAbs[j * _HC_Lt + i]);
                    }
                }
                polyakovSlice.AddItem(thisConfigurationSlice);
                if (m_bMeasureAbs)
                {
                    polyakovSliceAbs.AddItem(thisConfigurationSliceAbs);
                }
            }
            sFileNameWrite.Format(_T("%s_%s_polyakovZ_TSlice.csv"), sCSV.c_str(), sOName.c_str());
            WriteComplexArray2(sFileNameWrite, polyakovSlice);
            if (m_bMeasureAbs)
            {
                sFileNameWrite.Format(_T("%s_%s_polyakovZ_TSliceAbs.csv"), sCSV.c_str(), sOName.c_str());
                WriteRealArray2(sFileNameWrite, polyakovSliceAbs);
            }
        }
    }
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================