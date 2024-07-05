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
#include "CMeasureWilsonLoopWithPath.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureWilsonLoopWithPath)

#pragma region kernles 

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelWilsonLoopWithPathPoint(
    const deviceSU3* __restrict__ pSU3,
    SSmallInt4 point, 
    INT* deviceDirs,
    BYTE byPathLength,
    BYTE byFieldId,
    CLGComplex* pRes
)
{
    const deviceSU3 loop = _deviceLink(pSU3, point, byPathLength, byFieldId, deviceDirs);
    pRes[0] = loop.Tr();
}

/**
 * 
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelWilsonLoopWithPath(
    const deviceSU3* __restrict__ pSU3,
    INT* deviceDirs,
    BYTE byPathLength,
    BYTE byFieldId,
#if _CLG_DOUBLEFLOAT
    CLGComplex* pRes
#else
    cuDoubleComplex* pRes
#endif
)
{
    intokernalInt4;

    const deviceSU3 loop = _deviceLink(pSU3, sSite4, byPathLength, byFieldId, deviceDirs);
#if _CLG_DOUBLEFLOAT
    pRes[uiSiteIndex] = loop.Tr();
#else
    pRes[uiSiteIndex] = _cToDouble(loop.Tr());
#endif
}

#pragma endregion

CMeasureWilsonLoopWithPath::~CMeasureWilsonLoopWithPath()
{
    if (NULL != m_pTmpDeviceRes)
    {
        checkCudaErrors(cudaFree(m_pTmpDeviceRes));
    }

    if (NULL != m_pDevicePath)
    {
        checkCudaErrors(cudaFree(m_pDevicePath));
    }
}

void CMeasureWilsonLoopWithPath::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    param.FetchValueArrayINT(_T("Path"), m_lstPath);
    if (m_lstPath.Num() < 1)
    {
        appCrucial(_T("CMeasureWilsonLoopWithPath invalid Path!\n"));
        return;
    }

    INT iValue = 0;
    param.FetchValueINT(_T("OnePoint"), iValue);
    m_bAllPoint = (0 == iValue);

    if (!m_bAllPoint)
    {
        TArray<INT> thePoint;
        param.FetchValueArrayINT(_T("Point"), thePoint);
        m_sPoint = _HC_Center;
        if (thePoint.Num() > 3)
        {
            m_sPoint.x = static_cast<SBYTE>(thePoint[0]);
            m_sPoint.y = static_cast<SBYTE>(thePoint[1]);
            m_sPoint.z = static_cast<SBYTE>(thePoint[2]);
            m_sPoint.w = static_cast<SBYTE>(thePoint[3]);
        }
    }

    INT* hostPath = (INT*)alloca(sizeof(INT)* m_lstPath.Num());
    for (INT i = 0; i < m_lstPath.Num(); ++i)
    {
        hostPath[i] = m_lstPath[i];
    }

    checkCudaErrors(cudaMalloc((void**)&m_pTmpDeviceRes, sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)&m_pDevicePath, sizeof(INT) * m_lstPath.Num()));
    checkCudaErrors(cudaMemcpy(m_pDevicePath, hostPath, sizeof(INT) * m_lstPath.Num(), cudaMemcpyHostToDevice));

    Reset();
}

void CMeasureWilsonLoopWithPath::OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pAcceptGauge || EFT_GaugeSU3 != pAcceptGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pAcceptGauge);
    if (m_lstPath.Num() < 1)
    {
        appCrucial(_T("CMeasureWilsonLoopWithPath invalid Path!\n"));
        return;
    }

    if (!m_bAllPoint)
    {
        CLGComplex hostRes[1];
        _kernelWilsonLoopWithPathPoint << <1, 1 >> > (pGaugeSU3->m_pDeviceData, m_sPoint, m_pDevicePath, 
            static_cast<BYTE>(m_lstPath.Num()), 
            pAcceptGauge->m_byFieldId, m_pTmpDeviceRes);
        checkCudaErrors(cudaMemcpy(hostRes, m_pTmpDeviceRes, sizeof(CLGComplex), cudaMemcpyDeviceToHost));
        UpdateComplexResult(hostRes[0]);

        if (m_bShowResult)
        {
            appGeneral(_T("loop = "));
            appPushLogDate(FALSE);
            LogGeneralComplex(hostRes[0], FALSE);
            appGeneral(_T("\n"));
            appPopLogDate();
        }

        return;
    }

    preparethread;

    _kernelWilsonLoopWithPath<< <block, threads>> > (
        pGaugeSU3->m_pDeviceData, m_pDevicePath, 
        static_cast<BYTE>(m_lstPath.Num()), 
        pAcceptGauge->m_byFieldId, _D_ComplexThreadBuffer);

#if _CLG_DOUBLEFLOAT
    CLGComplex sum = appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
#else
    CLGComplex sum = _cToFloat(appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer));
#endif
    sum.x = sum.x / _HC_Volume;
    sum.y = sum.y / _HC_Volume;

    UpdateComplexResult(sum);

    if (m_bShowResult)
    {
        appGeneral(_T("loop = "));
        appPushLogDate(FALSE);
        LogGeneralComplex(sum, FALSE);
        appGeneral(_T("\n"));
        appPopLogDate();
    }
}

void CMeasureWilsonLoopWithPath::Report()
{
    Average();
    appPushLogDate(FALSE);

    appGeneral(_T("{\n"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        LogGeneralComplex(CmpResAtI(i), i + 1 < m_uiConfigurationCount);
    }
    appGeneral(_T("\n}\n"));
    if (m_uiConfigurationCount > 0)
    {
        appGeneral(_T("Average = "));
        LogGeneralComplex(GetAverageCmpRes(), FALSE);
        appGeneral(_T("\n\n"));
    }

    appPopLogDate();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================