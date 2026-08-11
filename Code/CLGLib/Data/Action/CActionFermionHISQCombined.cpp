//=============================================================================
// FILENAME : CActionFermionHISQCombined.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [mm/dd/yy]
//  [01/26/2025 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionFermionHISQCombined.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionFermionHISQCombined)

CActionFermionHISQCombined::CActionFermionHISQCombined()
    : CAction()
{
}


void CActionFermionHISQCombined::Initial(CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);
    TArray<INT> fermionIds;
    param.FetchValueArrayINT(_T("FieldIds"), fermionIds);
    //find fermion field
    appAssert(fermionIds.Num() > 0);
    for (INT i = 0; i < fermionIds.Num(); ++i)
    {
        CFieldFermion* pFerimionField = dynamic_cast<CFieldFermion*>(appGetLattice()->GetFieldById(static_cast<BYTE>(fermionIds[i])));
        if (NULL == pFerimionField)
        {
            appCrucial(_T("CActionFermionKS work with only CFieldFermionKS!\n"));
            continue;
        }
        //if (!pFerimionField->m_bEvenPseudofermion)
        //{
        //    appCrucial(_T("CActionFermionKS work with only Even CFieldFermionKS!\n"));
        //    continue;
        //}
        m_pFerimionField.AddItem(pFerimionField);
    }
}

void CActionFermionHISQCombined::PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate)
{
    _RECORD(CActionFermionHISQCombined::PrepareForHMC);
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
UBOOL CActionFermionHISQCombined::CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
    CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
    CFieldGauge* const* stapleFields, ESolverPhase ePhase) const
{
    _RECORD(CActionFermionHISQCombined::CalculateForce);
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

    //f0 part
    CFieldGauge* f0 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byGaugeFieldIds[0], _T(__FILE__), __LINE__));
    CFieldGauge* naikf0 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byGaugeFieldIds[0], _T(__FILE__), __LINE__));
    CFieldGauge* epsilonterm = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byGaugeFieldIds[0], _T(__FILE__), __LINE__));
    
    f0->Zero();
    naikf0->Zero();
    epsilonterm->Zero();
    INT idx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);

    for (INT i = 0; i < m_pFerimionField.Num(); ++i)
    {
        m_pFerimionField[i]->CalculateF0AndNaik(effectivegauge[idx], f0, epsilonterm, naikf0);
    }

    CFieldGauge* naikforce = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byGaugeFieldIds[0], _T(__FILE__), __LINE__));
    naikforce->Zero();
    naikforce->AddNaikForce(naikf0);
    naikf0->Return();
    naikforce->AxpyPlus(epsilonterm);
    epsilonterm->Return();
    pSmearing->DerivateOnU(effectivegauge[idx], gaugeFields[idx], naikforce, f0);
    naikforce->Return();

    //naik part
    //f0->LeftMul(gaugeFields[idx], FALSE, TRUE);
    //f0->TA();
    gaugeForces[idx]->AxpyPlus(f0);
    f0->Return();
    return TRUE;
}

DOUBLE CActionFermionHISQCombined::Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* stapleFields)
{
    _RECORD(CActionFermionHISQCombined::Energy);
    //[ (DD)^(-1/4) phi ]^2
    
    //CFieldFermionKSSU3* pPooled = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(static_cast<BYTE>(m_pFerimionField->m_byFieldId)));
    //appAssert(NULL != pPooled);
    //m_pFerimionField->CopyTo(pPooled);
    //pPooled->D_EN(pGauge);
    //const CLGComplex res = pPooled->Dot(pPooled);

    const CGaugeSmearing* pSmearing = appGetGaugeSmearing(m_byGaugeFieldIds[0]);
    //CFieldGauge* peffective = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byGaugeFieldIds[0]));
    TArray<const CFieldGauge*> effectivegauge;
    for (INT i = 0; i < gaugeNum; ++i)
    {
        if (NULL != gaugeFields[i] && gaugeFields[i]->m_byFieldId == m_byGaugeFieldIds[0])
        {
            //gaugeFields[i]->CopyTo(peffective);
            //pSmearing->GaugeSmearing(peffective, NULL);
            //effectivegauge.AddItem(peffective);
            effectivegauge.AddItem(pSmearing->GetEffectiveGauge());
        }
        else
        {
            effectivegauge.AddItem(gaugeFields[i]);
        }
    }
    DOUBLE fres = 0.0;
    for (INT i = 0; i < m_pFerimionField.Num(); ++i)
    {
        //CFieldFermionKS* pPooled = dynamic_cast<CFieldFermionKS*>(appGetLattice()->GetPooledFieldById(static_cast<BYTE>(m_pFerimionField[i]->m_byFieldId)));
        //m_pFerimionField[i]->CopyTo(pPooled);
        //pPooled->D_MD(gaugeNum, bosonNum, effectivegauge.GetData(), bosonFields);
        //const cuDoubleComplex res = pPooled->Dot(m_pFerimionField[i]);
        //appDetailed(_T("CActionFermionKS : Energy = %f%s%fi\n"), res.x, res.y > 0 ? "+" : " ", res.y);
        //fres += res.x;
        //pPooled->Return();
        fres += m_pFerimionField[i]->Energy(gaugeNum, bosonNum, tensor2Num, effectivegauge.GetData(), bosonFields, tensor2Fields);
    }
    return fres;
}


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
