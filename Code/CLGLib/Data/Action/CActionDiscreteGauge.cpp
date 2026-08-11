//=============================================================================
// FILENAME : CActionDiscreteGauge.cpp
//
// DESCRIPTION:
// Base class for discrete gauge actions
//
// REVISION:
//  [05/16/2026 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionDiscreteGauge.h"

__BEGIN_NAMESPACE

CActionDiscreteGauge::CActionDiscreteGauge()
    : CAction()
    , m_fBeta(1.0)
    , m_pDeviceFuncPtr(NULL)
{
}

CActionDiscreteGauge::~CActionDiscreteGauge()
{
    cudaSafeFree(m_pDeviceFuncPtr);
}

void CActionDiscreteGauge::Initial(CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    if (!param.FetchValueDOUBLE(_T("Beta"), m_fBeta))
    {
        m_fBeta = 1.0;
    }

    DeviceDiscreteLinkActionFunc hFunc = GetDeviceFunc();
    if (NULL == hFunc)
    {
        appCrucial(_T("CActionDiscreteGauge: GetDeviceFunc returned NULL!\n"));
        _FAIL_EXIT;
    }

    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceFuncPtr, sizeof(DeviceDiscreteLinkActionFunc)));
    checkCudaErrors(cudaMemcpy(m_pDeviceFuncPtr, &hFunc, sizeof(DeviceDiscreteLinkActionFunc), cudaMemcpyHostToDevice));
}

CCString CActionDiscreteGauge::GetInfos(const CCString& tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("Beta : ") + appToString(m_fBeta) + _T("\n");
    return sRet;
}

DOUBLE CActionDiscreteGauge::EnergySingleField(UBOOL bBeforeEvolution, const CFieldGauge* pGauge, const CFieldGauge* pStaple)
{
    appCrucial(_T("CActionDiscreteGauge: EnergySingleField not applicable for discrete gauge groups!\n"));
    return 0.0;
}

UBOOL CActionDiscreteGauge::CalculateForceOnGaugeSingleField(const CFieldGauge* pGauge, CFieldGauge* pForce, CFieldGauge* pStaple, ESolverPhase ePhase) const
{
    appCrucial(_T("CActionDiscreteGauge: CalculateForce not applicable for discrete gauge groups!\n"));
    return FALSE;
}

void CActionDiscreteGauge::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    appCrucial(_T("CActionDiscreteGauge: PrepareForHMC not applicable for discrete gauge groups!\n"));
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
