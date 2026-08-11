//=============================================================================
// FILENAME : CHeatbath.cu
//
// DESCRIPTION:
// Heat bath updater for discrete gauge groups (Z_N and D_N)
// GPU kernels and host dispatch for Z_N and D_N heat bath updates.
// Supports multiple gauge fields with independent actions and betas --
// each field is updated via its own kernel dispatch within a single Update() call.
//
// REVISION:
//  [05/16/2026 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Update/Discrete/CHeatbath.h"
#include "Data/Action/CActionDiscreteGauge.h"
#include "Data/Field/Gauge/CFieldGaugeLink.h"
#include "Data/Field/Gauge/CFieldGaugeKernel.h"
#include "Data/Lattice/CLatticeData.h"

__BEGIN_NAMESPACE

// ============================================
// Discrete group Heat Bath GPU Kernel
// ============================================
// For each link, compute the action weight w_k = pActionFunc(staple, k, beta)
// and select k with probability p_k = exp(-w_k) / sum_j exp(-w_j).
// The multi-gauge fields action need special treatment, so we support only single gauge field action for now.
// The staple is optional(might be NULL)
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND 
_kernelDiscreteHeatbath(
    deviceGauge * pGauge, 
    const deviceGauge * __restrict__ pStaple,
    DeviceDiscreteLinkActionFunc* pActionFunc, 
    DOUBLE beta, 
    UBOOL bEven,
    const BYTE* __restrict__ pEtaTable, UINT linkCount)
{
    intokernalDirEO_NoDir;

    if (uiLinkIndex >= linkCount)
    {
        return;
    }

    DOUBLE weights[deviceGauge::elementN];
    DOUBLE totalWeight = 0.0;

    for (INT k = 0; k < deviceGauge::elementN; ++k)
    {
        //field count, field ptrs, link index, group element index, beta
        DOUBLE w = (*pActionFunc)(pGauge, pStaple, uiLinkIndex, k, beta);
        weights[k] = _exp(-w);
        totalWeight += weights[k];
    }

    DOUBLE r = _deviceRandomDOUBLE(uiLinkIndex) * totalWeight;
    DOUBLE cumulative = 0.0;
    INT selectedK = 0;

    for (INT k = 0; k < deviceGauge::elementN; ++k)
    {
        cumulative += weights[k];
        if (r <= cumulative)
        {
            selectedK = k;
            break;
        }
    }

    _SetAsK(pGauge[uiLinkIndex], selectedK);
}

//template<INT N>
//void _DiscreteHeatbathKernelDispatch(deviceZN<N>* pGauge, const deviceZN<N>* pStaple, DeviceDiscreteLinkActionFunc* pActionFunc, DOUBLE beta, UBOOL bEven, UINT linkCount)
//{
//    preparethreadDir
//    _LAUNCH_KERNEL(_kernelZNHeatbath TMPARG(N), block, threads, pGauge, pStaple, pActionFunc, beta, bEven,
//        appGetLattice()->m_pIndexCache->m_pEtaMu, linkCount);
//}
//
//template<INT N>
//void _DNHeatbathKernelDispatch(deviceDN<N>* pGauge, const deviceDN<N>* pStaple, DeviceDiscreteLinkActionFunc* pActionFunc, DOUBLE beta, UBOOL bEven, UINT linkCount)
//{
//    preparethreadDir
//    _LAUNCH_KERNEL(_kernelDNHeatbath TMPARG(N), block, threads, pGauge, pStaple, pActionFunc, beta, bEven,
//        appGetLattice()->m_pIndexCache->m_pEtaMu, linkCount);
//}
//
//template<typename deviceGauge, INT matrixN>
//void _DiscreteComputeStaple(const deviceGauge* pGauge, deviceGauge* pStaple, BYTE byFieldId)
//{
//    CFieldGaugeKernel<deviceGauge, matrixN>::CalculateOnlyStaple(pGauge, byFieldId, pStaple);
//}

CHeatbath::CHeatbath()
    : CUpdator()
{
}

CHeatbath::~CHeatbath()
{

}

void CHeatbath::Initial(CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;

    //if (0 == pOwner->m_pGaugeField.Num())
    //{
    //    appCrucial(_T("CHeatbath requires at least one gauge field!\n"));
    //    _FAIL_EXIT;
    //}

    //m_Fields.RemoveAll();

    //for (INT i = 0; i < pOwner->m_pGaugeField.Num(); ++i)
    //{
    //    CFieldGauge* pGauge = pOwner->m_pGaugeField[i];
    //    if (NULL == pGauge) continue;

    //    // Find matching discrete action for this gauge field
    //    CActionDiscreteGauge* pAction = NULL;
    //    for (INT j = 0; j < pOwner->m_pActionList.Num(); ++j)
    //    {
    //        CActionDiscreteGauge* pCandidate = dynamic_cast<CActionDiscreteGauge*>(pOwner->m_pActionList[j]);
    //        if (NULL != pCandidate
    //            && pCandidate->GetPrimaryGaugeFieldId() == pGauge->GetFieldId())
    //        {
    //            pAction = pCandidate;
    //            break;
    //        }
    //    }
    //    if (NULL == pAction)
    //    {
    //        appCrucial(_T("CHeatbath: no CActionDiscreteGauge found for gauge field id %d!\n"), (INT)pGauge->GetFieldId());
    //        _FAIL_EXIT;
    //    }

    //    SHeatbathField field;
    //    field.m_pGaugeField = pGauge;
    //    field.m_pAction = pAction;
    //    field.m_pStapleField = NULL;
    //    m_Fields.AddItem(field);
    //}

    //if (0 == m_Fields.Num())
    //{
    //    appCrucial(_T("CHeatbath: no valid gauge field / action pairs found!\n"));
    //    _FAIL_EXIT;
    //}
}

CCString CHeatbath::GetInfos(const CCString& tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CHeatbath\n");
    //for (INT i = 0; i < m_Fields.Num(); ++i)
    //{
    //    sRet += tab + _T("  Field[") + appToString(i) + _T("]")
    //         + _T(" type=") + appToString(static_cast<INT>(m_Fields[i].m_pGaugeField->GetFieldType()))
    //         + _T(", Beta=") + appToString(m_Fields[i].m_pAction->GetBeta()) + _T("\n");
    //}
    return sRet;
}

template<typename deviceGauge>
void CHeatbath::UpdateOneField(deviceGauge* pGauge, deviceGauge* pStaple, DeviceDiscreteLinkActionFunc* pSweapFuncion, DOUBLE fBeta, UBOOL bEven, BYTE byFieldId)
{
    CFieldGaugeKernel<deviceGauge, deviceGauge::matrixN>::CalculateOnlyStaple(pGauge, byFieldId, pStaple);
    UINT uiLinkCount = static_cast<UINT>(_HC_Volume * _HC_Dir);
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelDiscreteHeatbath TMPARG(deviceGauge), block, threads, pGauge, pStaple, pSweapFuncion, fBeta, bEven,
        appGetLattice()->m_pIndexCache->m_pEtaMu, 
        uiLinkCount
        );
}

//void CHeatbath::DispatchHeatbathKernel(SHeatbathField& field, UBOOL bEven)
//{
//    BYTE byFieldId = field.m_pGaugeField->GetFieldId();
//    DeviceDiscreteLinkActionFunc* pFunc = field.m_pAction->GetDeviceFuncPtr();
//    DOUBLE fBeta = field.m_pAction->GetBeta();
//    UINT linkCount = field.m_pGaugeField->GetLinkCount();
//
//    if (NULL == field.m_pStapleField)
//    {
//        field.m_pStapleField = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));
//    }
//    void* pStapleData = field.m_pStapleField->GetData();
//
//    switch (field.m_pGaugeField->GetFieldType())
//    {
//#if _CLG_Z2_GAUGE
//    case EFT_GaugeZ2:
//        _DiscreteComputeStaple<deviceZN<2>, 1>((const deviceZN<2>*)field.m_pGaugeField->GetData(), (deviceZN<2>*)pStapleData, byFieldId);
//        _DiscreteHeatbathKernelDispatch<2>((deviceZN<2>*)field.m_pGaugeField->GetData(), (const deviceZN<2>*)pStapleData, pFunc, fBeta, bEven, linkCount);
//        break;
//#endif
//#if _CLG_Z3_GAUGE
//    case EFT_GaugeZ3:
//        _DiscreteComputeStaple<deviceZN<3>, 1>((const deviceZN<3>*)field.m_pGaugeField->GetData(), (deviceZN<3>*)pStapleData, byFieldId);
//        _DiscreteHeatbathKernelDispatch<3>((deviceZN<3>*)field.m_pGaugeField->GetData(), (const deviceZN<3>*)pStapleData, pFunc, fBeta, bEven, linkCount);
//        break;
//#endif
//#if _CLG_Z4_GAUGE
//    case EFT_GaugeZ4:
//        _DiscreteComputeStaple<deviceZN<4>, 1>((const deviceZN<4>*)field.m_pGaugeField->GetData(), (deviceZN<4>*)pStapleData, byFieldId);
//        _DiscreteHeatbathKernelDispatch<4>((deviceZN<4>*)field.m_pGaugeField->GetData(), (const deviceZN<4>*)pStapleData, pFunc, fBeta, bEven, linkCount);
//        break;
//#endif
//#if _CLG_Z5_GAUGE
//    case EFT_GaugeZ5:
//        _DiscreteComputeStaple<deviceZN<5>, 1>((const deviceZN<5>*)field.m_pGaugeField->GetData(), (deviceZN<5>*)pStapleData, byFieldId);
//        _DiscreteHeatbathKernelDispatch<5>((deviceZN<5>*)field.m_pGaugeField->GetData(), (const deviceZN<5>*)pStapleData, pFunc, fBeta, bEven, linkCount);
//        break;
//#endif
//#if _CLG_Z6_GAUGE
//    case EFT_GaugeZ6:
//        _DiscreteComputeStaple<deviceZN<6>, 1>((const deviceZN<6>*)field.m_pGaugeField->GetData(), (deviceZN<6>*)pStapleData, byFieldId);
//        _DiscreteHeatbathKernelDispatch<6>((deviceZN<6>*)field.m_pGaugeField->GetData(), (const deviceZN<6>*)pStapleData, pFunc, fBeta, bEven, linkCount);
//        break;
//#endif
//#if _CLG_D3_GAUGE
//    case EFT_GaugeD3:
//        _DiscreteComputeStaple<deviceDN<3>, 2>((const deviceDN<3>*)field.m_pGaugeField->GetData(), (deviceDN<3>*)pStapleData, byFieldId);
//        _DNHeatbathKernelDispatch<3>((deviceDN<3>*)field.m_pGaugeField->GetData(), (const deviceDN<3>*)pStapleData, pFunc, fBeta, bEven, linkCount);
//        break;
//#endif
//#if _CLG_D4_GAUGE
//    case EFT_GaugeD4:
//        _DiscreteComputeStaple<deviceDN<4>, 2>((const deviceDN<4>*)field.m_pGaugeField->GetData(), (deviceDN<4>*)pStapleData, byFieldId);
//        _DNHeatbathKernelDispatch<4>((deviceDN<4>*)field.m_pGaugeField->GetData(), (const deviceDN<4>*)pStapleData, pFunc, fBeta, bEven, linkCount);
//        break;
//#endif
//#if _CLG_D8_GAUGE
//    case EFT_GaugeD8:
//        _DiscreteComputeStaple<deviceDN<8>, 2>((const deviceDN<8>*)field.m_pGaugeField->GetData(), (deviceDN<8>*)pStapleData, byFieldId);
//        _DNHeatbathKernelDispatch<8>((deviceDN<8>*)field.m_pGaugeField->GetData(), (const deviceDN<8>*)pStapleData, pFunc, fBeta, bEven, linkCount);
//        break;
//#endif
//    default:
//        appCrucial(_T("CHeatbath: unsupported gauge field type in DispatchHeatbathKernel!\n"));
//        _FAIL_EXIT;
//    }
//}

void CHeatbath::UpdateOneParity(UBOOL bEven)
{
    for (INT i = 0; i < m_pOwner->m_pActionList.Num(); ++i)
    {
        if (m_pOwner->m_pActionList[i]->IsDiscreteGauge())
        {
            CActionDiscreteGauge* pAction = dynamic_cast<CActionDiscreteGauge*>(m_pOwner->m_pActionList[i]);
            if (NULL == pAction)
            {
                appCrucial(_T("CHeatbath: action in action list is marked as discrete gauge but cannot be cast to CActionDiscreteGauge!\n"));
                continue;
            }
            BYTE byFieldId = pAction->GetPrimaryGaugeFieldId();
            CFieldGauge* pGaugeField = dynamic_cast<CFieldGauge*>(appGetLattice()->GetFieldById(byFieldId));
            if (NULL == pGaugeField)
            {
                appCrucial(_T("CHeatbath: no gauge field found for field id %d!\n"), (INT)byFieldId);
                continue;
            }

            CFieldGauge* pStapleField = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(byFieldId, _T(__FILE__), __LINE__));

            switch (pGaugeField->GetFieldType())
            {
#if _CLG_Z2_GAUGE
            case EFT_GaugeZ2:
                UpdateOneField<deviceZN<2>>((deviceZN<2>*)pGaugeField->GetData(), (deviceZN<2>*)pStapleField->GetData(), pAction->GetDeviceFuncPtr(), pAction->GetBeta(), bEven, byFieldId);
                break;
#endif
#if _CLG_Z3_GAUGE
            case EFT_GaugeZ3:
                UpdateOneField<deviceZN<3>>((deviceZN<3>*)pGaugeField->GetData(), (deviceZN<3>*)pStapleField->GetData(), pAction->GetDeviceFuncPtr(), pAction->GetBeta(), bEven, byFieldId);
                break;
#endif
#if _CLG_Z4_GAUGE
            case EFT_GaugeZ4:
                UpdateOneField<deviceZN<4>>((deviceZN<4>*)pGaugeField->GetData(), (deviceZN<4>*)pStapleField->GetData(), pAction->GetDeviceFuncPtr(), pAction->GetBeta(), bEven, byFieldId);
                break;
#endif
#if _CLG_Z5_GAUGE
            case EFT_GaugeZ5:
                UpdateOneField<deviceZN<5>>((deviceZN<5>*)pGaugeField->GetData(), (deviceZN<5>*)pStapleField->GetData(), pAction->GetDeviceFuncPtr(), pAction->GetBeta(), bEven, byFieldId);
                break;
#endif
#if _CLG_Z6_GAUGE
            case EFT_GaugeZ6:
                UpdateOneField<deviceZN<6>>((deviceZN<6>*)pGaugeField->GetData(), (deviceZN<6>*)pStapleField->GetData(), pAction->GetDeviceFuncPtr(), pAction->GetBeta(), bEven, byFieldId);
                break;
#endif
#if _CLG_D3_GAUGE
            case EFT_GaugeD3:
                UpdateOneField<deviceDN<3>>((deviceDN<3>*)pGaugeField->GetData(), (deviceDN<3>*)pStapleField->GetData(), pAction->GetDeviceFuncPtr(), pAction->GetBeta(), bEven, byFieldId);
                break;
#endif
#if _CLG_D4_GAUGE
            case EFT_GaugeD4:
                UpdateOneField<deviceDN<4>>((deviceDN<4>*)pGaugeField->GetData(), (deviceDN<4>*)pStapleField->GetData(), pAction->GetDeviceFuncPtr(), pAction->GetBeta(), bEven, byFieldId);
                break;
#endif
#if _CLG_D8_GAUGE
            case EFT_GaugeD8:
                UpdateOneField<deviceDN<8>>((deviceDN<8>*)pGaugeField->GetData(), (deviceDN<8>*)pStapleField->GetData(), pAction->GetDeviceFuncPtr(), pAction->GetBeta(), bEven, byFieldId);
                break;
#endif
            default:
                appCrucial(_T("CHeatbath: unsupported gauge field type in DispatchHeatbathKernel!\n"));
                break;
            }
            pStapleField->Return();
        }
    }
}

UINT CHeatbath::Update(UINT iSteps, UBOOL bMeasure)
{
    //if (0 == m_Fields.Num())
    //{
    //    appCrucial(_T("CHeatbath is not initialized!\n"));
    //    _FAIL_EXIT;
    //}

    for (UINT i = 0; i < iSteps; ++i)
    {
        UpdateOneParity(FALSE);
        UpdateOneParity(TRUE);
        m_pOwner->FixAllFieldBoundary();
        checkCudaErrors(cudaDeviceSynchronize());
        ++m_uiUpdateCall;

        if (bMeasure)
        {
            //TArray<const CFieldGauge*> gauges;
            //for (INT j = 0; j < m_Fields.Num(); ++j)
            //{
            //    gauges.AddItem(m_Fields[j].m_pGaugeField);
            //}
            
            m_pOwner->OnUpdatorConfigurationAccepted(m_pOwner->m_pGaugeField.Num(), 
                m_pOwner->m_pBosonField.Num(), 
                m_pOwner->m_pTensor2Field.Num(), 
                m_pOwner->m_pGaugeField.GetData(), 
                m_pOwner->m_pBosonField.GetData(), 
                m_pOwner->m_pTensor2Field.GetData(), 
                NULL);
        }
    }

    checkCudaErrors(cudaDeviceSynchronize());
    m_pOwner->OnUpdatorFinished(bMeasure, m_bReport);
    return iSteps;
}

__CLGIMPLEMENT_CLASS(CHeatbath)

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
