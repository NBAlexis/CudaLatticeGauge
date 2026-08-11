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
#include "CActionFermionKSImproveCombined.h"
#include "Data/Field/Staggered/CFieldFermionKSHISQ.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionFermionKSImproveCombined)

CActionFermionKSImproveCombined::CActionFermionKSImproveCombined()
    : CActionFermionKSCombined()
{

}

void CActionFermionKSImproveCombined::PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate)
{
    _RECORD(CActionFermionKSImproveCombined::PrepareForHMC);
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
    for (INT i = 0; i < m_pFerimionField.Num(); ++i)
    {
        m_pFerimionField[i]->PrepareForHMC(gaugeNum, bosonNum, effectivegauge.GetData(), bosonFields);
    }
    
}

/**
* To make it constant, we need to build a few temp fields outside this class
*/
UBOOL CActionFermionKSImproveCombined::CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
    CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
    CFieldGauge* const* stapleFields, ESolverPhase ePhase) const
{
    _RECORD(CActionFermionKSImproveCombined::CalculateForce);
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
        for (INT i = 0; i < m_pFerimionField.Num(); ++i)
        {
            bRet = m_pFerimionField[i]->CalculateForce(gaugeNum, bosonNum, effectivegauge.GetData(), bosonFields, effectivegaugeforce.GetData(), bosonForces, ePhase);
        }
        pSmearing->DerivateOnU(effectivegauge[uiGaugeId], gaugeFields[uiGaugeId], NULL, f0);
        gaugeForces[uiGaugeId]->AxpyPlus(f0);
    }
    else
    {
        for (INT i = 0; i < m_pFerimionField.Num(); ++i)
        {
            m_pFerimionField[i]->CalculateForce(gaugeNum, bosonNum, gaugeFields, bosonFields, gaugeForces, bosonForces, ePhase);
        }
    }
    f0->Return();
    return bRet;
}

DOUBLE CActionFermionKSImproveCombined::Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* stapleFields)
{
    _RECORD(CActionFermionKSImproveCombined::Energy);
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

    DOUBLE res = 0.0;
    for (INT i = 0; i < m_pFerimionField.Num(); ++i)
    {
        res += m_pFerimionField[i]->Energy(gaugeNum, bosonNum, tensor2Num, effectivegauge.GetData(), bosonFields, tensor2Fields);
    }
    return res;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================