//=============================================================================
// FILENAME : CActionFermionKSImprove.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [mm/dd/yy]
//  [11/12/2024 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionFermionKSImprove.h"
#include "Data/Field/Staggered/CFieldFermionKSHISQ.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionFermionKSImprove)

CActionFermionKSImprove::CActionFermionKSImprove()
    : CActionFermionKS()
    , m_bHISQ(FALSE)
{

}

CActionFermionKSImprove::~CActionFermionKSImprove()
{

}

void CActionFermionKSImprove::Initial(CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CActionFermionKS::Initial(pOwner, param, byId);

    if (m_byGaugeFieldIds.Num() < 1)
    {
        m_byGaugeFieldIds.AddItem(1);
    }

    INT iValue = 0;
    if (param.FetchValueINT(_T("HISQ"), iValue))
    {
        m_bHISQ = (0 != iValue);
    }
}

void CActionFermionKSImprove::PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate)
{
    _RECORD(CActionFermionKSImprove::PrepareForHMC);
    const CGaugeSmearing* pSmearing = appGetGaugeSmearing(m_byGaugeFieldIds[0]);
    TArray<const CFieldGauge*> effectivegauge;
    for (INT i = 0; i < gaugeNum; ++i)
    {
        if (NULL != gaugeFields[i] && gaugeFields[i]->m_byFieldId == m_byGaugeFieldIds[0])
        {
            effectivegauge.AddItem(pSmearing->GetEffectiveGauge());
        }
        else
        {
            effectivegauge.AddItem(gaugeFields[i]);
        }
    }
    m_pFerimionField->PrepareForHMC(gaugeNum, bosonNum, effectivegauge.GetData(), bosonFields);
}

/**
* To make it constant, we need to build a few temp fields outside this class
*/
UBOOL CActionFermionKSImprove::CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
    CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
    CFieldGauge* const* stapleFields, ESolverPhase ePhase) const
{
    _RECORD(CActionFermionKSImprove::CalculateForce);
    const CGaugeSmearing* pSmearing = appGetGaugeSmearing(m_byGaugeFieldIds[0]);
    CFieldGauge* f0 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byGaugeFieldIds[0], _T(__FILE__), __LINE__));
    f0->Zero();
    TArray<const CFieldGauge*> effectivegauge;
    TArray<CFieldGauge*> effectivegaugeforce;
    INT uiGaugeId = -1;
    for (INT i = 0; i < gaugeNum; ++i)
    {
        if (NULL != gaugeFields[i] && gaugeFields[i]->m_byFieldId == m_byGaugeFieldIds[0])
        {
            //gaugeFields[i]->CopyTo(peffective);
            //pSmearing->GaugeSmearing(peffective, NULL);
            //effectivegauge.AddItem(peffective);
            effectivegauge.AddItem(pSmearing->GetEffectiveGauge());
            effectivegaugeforce.AddItem(f0);
            uiGaugeId = i;
        }
        else
        {
            effectivegauge.AddItem(gaugeFields[i]);
            effectivegaugeforce.AddItem(gaugeForces[i]);
        }
    }

    UBOOL bRet = FALSE;
    if (uiGaugeId >= 0)
    {
        if (m_bHISQ)
        {
            CFieldGauge* naikf0 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byGaugeFieldIds[0], _T(__FILE__), __LINE__));
            CFieldGauge* epsilonterm = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byGaugeFieldIds[0], _T(__FILE__), __LINE__));
            naikf0->Zero();
            epsilonterm->Zero();
            m_pFerimionField->CalculateF0AndNaik(effectivegauge[uiGaugeId], f0, epsilonterm, naikf0);
            CFieldGauge* naikforce = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byGaugeFieldIds[0], _T(__FILE__), __LINE__));
            naikforce->Zero();
            naikforce->AddNaikForce(naikf0);
            naikf0->Return();
            naikforce->AxpyPlus(epsilonterm);
            epsilonterm->Return();
            pSmearing->DerivateOnU(effectivegauge[uiGaugeId], gaugeFields[uiGaugeId], naikforce, f0);
            naikforce->Return();
            gaugeForces[uiGaugeId]->AxpyPlus(f0);
            bRet = TRUE;
        }
        else
        {
            bRet = m_pFerimionField->CalculateForce(gaugeNum, bosonNum, effectivegauge.GetData(), bosonFields, effectivegaugeforce.GetData(), bosonForces, ePhase);
            pSmearing->DerivateOnU(effectivegauge[uiGaugeId], gaugeFields[uiGaugeId], NULL, f0);
            gaugeForces[uiGaugeId]->AxpyPlus(f0);
        }
    }
    else
    {
        bRet = m_pFerimionField->CalculateForce(gaugeNum, bosonNum, gaugeFields, bosonFields, gaugeForces, bosonForces, ePhase);
    }
    f0->Return();

    return bRet;
}

DOUBLE CActionFermionKSImprove::Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* stapleFields)
{
    _RECORD(CActionFermionKSImprove::Energy);
    const CGaugeSmearing* pSmearing = appGetGaugeSmearing(m_byGaugeFieldIds[0]);
    TArray<const CFieldGauge*> effectivegauge;
    for (INT i = 0; i < gaugeNum; ++i)
    {
        if (NULL != gaugeFields[i] && gaugeFields[i]->m_byFieldId == m_byGaugeFieldIds[0])
        {
            effectivegauge.AddItem(pSmearing->GetEffectiveGauge());
        }
        else
        {
            effectivegauge.AddItem(gaugeFields[i]);
        }
    }

    return m_pFerimionField->Energy(gaugeNum, bosonNum, tensor2Num, effectivegauge.GetData(), bosonFields, tensor2Fields);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================