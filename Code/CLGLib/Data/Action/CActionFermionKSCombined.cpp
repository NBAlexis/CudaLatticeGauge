//=============================================================================
// FILENAME : CActionFermionKSCombined.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [mm/dd/yy]
//  [01/26/2025 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionFermionKSCombined.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionFermionKSCombined)

CActionFermionKSCombined::CActionFermionKSCombined()
    : CAction()
{
}


void CActionFermionKSCombined::Initial(CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);
    TArray<INT> fermionIds;
    param.FetchValueArrayINT(_T("FieldIds"), fermionIds);
    //find fermion field
    appAssert(fermionIds.Num() > 0);
    for (INT i = 0; i < fermionIds.Num(); ++i)
    {
        CFieldFermion* pFerimionField = dynamic_cast<CFieldFermion*>(appGetLattice()->GetFieldById(static_cast<BYTE>(fermionIds[i])));
        //if (NULL == pFerimionField)
        //{
        //    appCrucial(_T("CActionFermionKS work with only CFieldFermionKS!\n"));
        //    continue;
        //}
        //if (!pFerimionField->m_bEvenPseudofermion)
        //{
        //    appCrucial(_T("CActionFermionKS work with only Even CFieldFermionKS!\n"));
        //    continue;
        //}
        m_pFerimionField.AddItem(pFerimionField);
    }
}

void CActionFermionKSCombined::PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate)
{
    for (INT i = 0; i < m_pFerimionField.Num(); ++i)
    {
        m_pFerimionField[i]->PrepareForHMC(gaugeNum, bosonNum, gaugeFields, bosonFields);
    }
}

/**
* To make it constant, we need to build a few temp fields outside this class
*/
UBOOL CActionFermionKSCombined::CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
    CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
    CFieldGauge* const* stapleFields, ESolverPhase ePhase) const
{
    //CFieldGauge* thisfieldforce = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byGaugeFieldIds[0]));
    //thisfieldforce->Zero();
    //INT idx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
    for (INT i = 0; i < m_pFerimionField.Num(); ++i)
    {
        m_pFerimionField[i]->CalculateForce(gaugeNum, bosonNum, gaugeFields, bosonFields, gaugeForces, bosonForces, ePhase);
    }
    //thisfieldforce->LeftMul(gaugeFields[idx], FALSE, TRUE);
    //thisfieldforce->TA();
    //gaugeForces[idx]->AxpyPlus(thisfieldforce);
    //thisfieldforce->Return();
    return TRUE;
}

DOUBLE CActionFermionKSCombined::Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* stapleFields)
{
    //[ (DD)^(-1/4) phi ]^2
    
    //CFieldFermionKSSU3* pPooled = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(static_cast<BYTE>(m_pFerimionField->m_byFieldId)));
    //appAssert(NULL != pPooled);
    //m_pFerimionField->CopyTo(pPooled);
    //pPooled->D_EN(pGauge);
    //const CLGComplex res = pPooled->Dot(pPooled);

    DOUBLE fres = 0.0;
    for (INT i = 0; i < m_pFerimionField.Num(); ++i)
    {
        fres += m_pFerimionField[i]->Energy(gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields);
    }
    return fres;
}


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================