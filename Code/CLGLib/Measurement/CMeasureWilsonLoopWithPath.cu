//=============================================================================
// FILENAME : CMeasureWilsonLoopWithPath.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/10/2021 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CMeasureWilsonLoopWithPath.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureWilsonLoopWithPath)

#pragma region kernles 

template<typename T>
__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelWilsonLoopWithPathPoint(
    const T* __restrict__ pSU3,
    SSmallInt4 point, 
    SCHAR* deviceDirs,
    BYTE byPathLength,
    BYTE byFieldId,
    cuDoubleComplex* pRes
)
{
    const T loop = _deviceLinkLongT(pSU3, point, byPathLength, byFieldId, deviceDirs);
    pRes[0] = _cToDouble(_tr(loop));
}

/**
 * 
 */
template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelWilsonLoopWithPath(
    const T* __restrict__ pSU3,
    SCHAR* deviceDirs,
    BYTE byPathLength,
    BYTE byFieldId,
    cuDoubleComplex* pRes)
{
    intokernalInt4;

    const T loop = _deviceLinkLongT(pSU3, sSite4, byPathLength, byFieldId, deviceDirs);
    pRes[uiSiteIndex] = _cToDouble(_tr(loop));
}

/**
* only support Torus boundary condition
*/
template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelWilsonLoopWithPathTwoPath(
    const T* __restrict__ pSU3,
    SCHAR* deviceDirs1,
    SCHAR* deviceDirs2,
    SSmallInt4 sShift,
    BYTE byPathLength,
    BYTE byFieldId,
    cuDoubleComplex* pRes,
    DOUBLE* pdoubleRes)
{
    intokernalInt4;

    const T loop1 = _deviceLinkLongT(pSU3, sSite4, byPathLength, byFieldId, deviceDirs1);

    SSmallInt4 sSite4Shifted = sSite4;

    sSite4Shifted.x += sShift.x;
    while (sSite4Shifted.x < 0)
    {
        sSite4Shifted.x += static_cast<SCHAR>(_DC_Lx);
    }
    while (sSite4Shifted.x >= _DC_Lx)
    {
        sSite4Shifted.x -= static_cast<SCHAR>(_DC_Lx);
    }

    sSite4Shifted.y += sShift.y;
    while (sSite4Shifted.y < 0)
    {
        sSite4Shifted.y += static_cast<SCHAR>(_DC_Ly);
    }
    while (sSite4Shifted.y >= _DC_Ly)
    {
        sSite4Shifted.y -= static_cast<SCHAR>(_DC_Ly);
    }

    sSite4Shifted.z += sShift.z;
    while (sSite4Shifted.z < 0)
    {
        sSite4Shifted.z += static_cast<SCHAR>(_DC_Lz);
    }
    while (sSite4Shifted.z >= _DC_Lz)
    {
        sSite4Shifted.z -= static_cast<SCHAR>(_DC_Lz);
    }

    sSite4Shifted.w += sShift.w;
    while (sSite4Shifted.w < 0)
    {
        sSite4Shifted.w += static_cast<SCHAR>(_DC_Lt);
    }
    while (sSite4Shifted.w >= _DC_Lt)
    {
        sSite4Shifted.w -= static_cast<SCHAR>(_DC_Lt);
    }

    const T loop2 = _deviceLinkLongT(pSU3, sSite4Shifted, byPathLength, byFieldId, deviceDirs2);

    pRes[uiSiteIndex] = _cToDouble(_mulC(_tr(loop1), _tr(loop2)));
    pdoubleRes[uiSiteIndex] = static_cast<DOUBLE>(_retr(loop1) * _retr(loop2));
}

#pragma endregion

CMeasureWilsonLoopWithPath::~CMeasureWilsonLoopWithPath()
{
    if (NULL != m_pTmpDeviceRes)
    {
        checkCudaErrors(__cudaFree(m_pTmpDeviceRes));
    }

    if (NULL != m_pDevicePath)
    {
        checkCudaErrors(__cudaFree(m_pDevicePath));
    }

    if (NULL != m_pDevicePath2)
    {
        checkCudaErrors(__cudaFree(m_pDevicePath2));
    }
}

void CMeasureWilsonLoopWithPath::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    TArray<SCHAR> path;
    param.FetchValueArraySCHAR(_T("Path"), path);
    if (path.Num() < 1)
    {
        appGeneral(_T("CMeasureWilsonLoopWithPath invalid Path!\n"));
    }
    else
    {
        m_lstPath.RemoveAll();
        m_lstPath.AddItem(path);
    }

    INT iValue = 0;
    param.FetchValueINT(_T("OnePoint"), iValue);
    m_bAllPoint = (0 == iValue);

    if (!m_bAllPoint)
    {
        TArray<SCHAR> thePoint;
        param.FetchValueArraySCHAR(_T("Point"), thePoint);
        m_sPoint = _HC_Center;
        if (thePoint.Num() > 3)
        {
            memcpy(m_sPoint.m_byData4, thePoint.GetData(), sizeof(SCHAR) * 4);
        }
        checkCudaErrors(__cudaMalloc((void**)&m_pTmpDeviceRes, sizeof(cuDoubleComplex)));
    }

    checkCudaErrors(__cudaMalloc((void**)&m_pDevicePath, sizeof(SCHAR) * static_cast<UINT>(_kMaxWilsonPathLength)));
    checkCudaErrors(__cudaMalloc((void**)&m_pDevicePath2, sizeof(SCHAR) * static_cast<UINT>(_kMaxWilsonPathLength)));
    Reset();
}

void CMeasureWilsonLoopWithPath::OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pAcceptGauge)
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }

    if (m_lstPath.Num() < 1)
    {
        appCrucial(_T("CMeasureWilsonLoopWithPath invalid Path!\n"));
        return;
    }

#if _CLG_MULTI_GPU
    //The halo only carries HaloWidth layers per direction, so a path whose
    //accumulated walk in any direction exceeds that depth would read stale
    //halo slots. Reject explicitly instead of reporting a wrong loop.
    if (NULL != appGetComm())
    {
        for (INT i = 0; i < m_lstPath.Num(); ++i)
        {
            INT walk[4] = { 0, 0, 0, 0 };
            INT maxSpan[4] = { 0, 0, 0, 0 };
            for (INT s = 0; s < m_lstPath[i].Num(); ++s)
            {
                const SCHAR dir = m_lstPath[i][s];
                const SCHAR sign = (dir < 0) ? -1 : 1;
                const SCHAR axis = appAbs(dir) - 1;
                if (axis >= 0 && axis < 4)
                {
                    walk[axis] += sign;
                    if (appAbs(walk[axis]) > maxSpan[axis]) { maxSpan[axis] = appAbs(walk[axis]); }
                }
            }
            for (INT a = 0; a < 4; ++a)
            {
                if (maxSpan[a] > static_cast<INT>(_HC_HaloWidth))
                {
                    appCrucial(_T("CMeasureWilsonLoopWithPath: path %d spans %d sites in direction %d, exceeding the halo width %d; cross-rank walks beyond the halo are not supported. Rejected.\n"),
                        i, maxSpan[a], a, static_cast<INT>(_HC_HaloWidth));
                    return;
                }
            }
        }
    }
#endif

    cuDoubleComplex hostRes[1];
    preparethread;

    TArray<cuDoubleComplex> res;
    TArray<DOUBLE> fres;
    if (m_bShowResult)
    {
        appGeneral(_T("loop = "));
        appPushLogDate(FALSE);
    }
    for (INT i = 0; i < m_lstPath.Num(); ++i)
    {
        if (m_lstPath[i].Num() > _kMaxWilsonPathLength)
        {
            appCrucial(_T("The loop length(%d) is larger than %d\n"), m_lstPath[i].Num(), _kMaxWilsonPathLength);
        }
        BYTE pathL = static_cast<BYTE>(appMin(m_lstPath[i].Num(), static_cast<INT>(_kMaxWilsonPathLength)));
        checkCudaErrors(cudaMemcpy(m_pDevicePath, m_lstPath[i].GetData(), sizeof(SCHAR) * pathL, cudaMemcpyHostToDevice));
        
        if (EFT_GaugeSU3 == pAcceptGauge->GetFieldType())
        {
            const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pAcceptGauge);

            if (m_bTwoPathBackForward)
            {
                checkCudaErrors(cudaMemcpy(m_pDevicePath2, m_lstPath[i].GetData(), sizeof(SCHAR) * pathL, cudaMemcpyHostToDevice));
                _LAUNCH_KERNEL(_kernelWilsonLoopWithPathTwoPath<deviceSU3>, block, threads, 
                    pGaugeSU3->m_pDeviceData, 
                    m_pDevicePath, 
                    m_pDevicePath2, 
                    m_sShift[i], 
                    pathL, 
                    pAcceptGauge->m_byFieldId, 
                    _D_ComplexThreadBuffer,
                    _D_RealThreadBuffer);

                cuDoubleComplex sum = appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
                //P4-3.9: local spatial sum -> global, normalized by the GLOBAL volume.
                appGlobalSum(sum);
                const DOUBLE fVol1 = GlobalL(0) * GlobalL(1) * GlobalL(2) * GlobalL(3);
                sum.x = sum.x / fVol1;
                sum.y = sum.y / fVol1;
                res.AddItem(sum);

                DOUBLE fSum = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
                GlobalSumReal(fSum);
                fSum = fSum / fVol1;
                fres.AddItem(fSum);

                TArray<SCHAR> reversePath = PathDagger(m_lstPath[i]);
                checkCudaErrors(cudaMemcpy(m_pDevicePath2, reversePath.GetData(), sizeof(SCHAR) * pathL, cudaMemcpyHostToDevice));
                _LAUNCH_KERNEL(_kernelWilsonLoopWithPathTwoPath<deviceSU3>, block, threads,
                    pGaugeSU3->m_pDeviceData,
                    m_pDevicePath,
                    m_pDevicePath2,
                    m_sShift[i],
                    pathL,
                    pAcceptGauge->m_byFieldId,
                    _D_ComplexThreadBuffer,
                    _D_RealThreadBuffer);

                sum = appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
                appGlobalSum(sum);
                sum.x = sum.x / fVol1;
                sum.y = sum.y / fVol1;
                res.AddItem(sum);

                fSum = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
                GlobalSumReal(fSum);
                fSum = fSum / fVol1;
                fres.AddItem(fSum);
            }
            else
            {
                if (!m_bAllPoint)
                {
                    _LAUNCH_KERNEL(_kernelWilsonLoopWithPathPoint<deviceSU3>, 1, 1, pGaugeSU3->m_pDeviceData, m_sPoint, m_pDevicePath, pathL, pAcceptGauge->m_byFieldId, m_pTmpDeviceRes);
                }
                else
                {
                    _LAUNCH_KERNEL(_kernelWilsonLoopWithPath<deviceSU3>, block, threads, pGaugeSU3->m_pDeviceData, m_pDevicePath, pathL, pAcceptGauge->m_byFieldId, _D_ComplexThreadBuffer);
                }
            }

        }
        else
        {
            appCrucial(_T("CMeasureWilsonLoopWithPath::OnConfigurationAcceptedSingleField, field type %d not supported!\n"), pAcceptGauge->GetFieldType());
            _FAIL_EXIT;
        }

        if (!m_bTwoPathBackForward)
        {
            if (!m_bAllPoint)
            {
                checkCudaErrors(cudaMemcpy(hostRes, m_pTmpDeviceRes, sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost));
                res.AddItem(hostRes[0]);
            }
            else
            {
                cuDoubleComplex sum = appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
                //P4-3.9: local spatial sum -> global, normalized by the GLOBAL volume.
                appGlobalSum(sum);
                const DOUBLE fVol2 = GlobalL(0) * GlobalL(1) * GlobalL(2) * GlobalL(3);
                sum.x = sum.x / fVol2;
                sum.y = sum.y / fVol2;
                res.AddItem(sum);
            }
        }

        if (m_bShowResult)
        {
            appGeneral(_T("%.6f %s %.6f I%s"),
                res[res.Num() - 1].x,
                res[res.Num() - 1].y < F(0.0) ? _T("-") : _T("+"),
                appAbs(res[res.Num() - 1].y),
                i != (m_lstPath.Num() - 1) ? _T(",   ") : _T(" "));
        }
    }
    if (m_bShowResult)
    {
        appGeneral(_T("\n"));
        appPopLogDate();
    }
    ++m_uiConfigurationCount;
    m_lstV.AddItem(res);
    if (m_bTwoPathBackForward)
    {
        m_lstDV.AddItem(fres);
    }

    if (NULL != m_pOwner)
    {
        TArray<CLGComplex> loopC;
        loopC.SetSize(res.Num());
        for (INT iPath = 0; iPath < res.Num(); ++iPath)
        {
            loopC[iPath] = _cToRealC(res[iPath]);
        }
        m_pOwner->AddOneConfigurationResult(this, _T("Loop"), loopC);
        if (m_bTwoPathBackForward)
        {
            m_pOwner->AddOneConfigurationResult(this, _T("LoopReal"), fres);
        }
    }
}

void CMeasureWilsonLoopWithPath::Report()
{
    //not supported
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================