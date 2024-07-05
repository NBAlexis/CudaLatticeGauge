//=============================================================================
// FILENAME : CActionFermionKS.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [06/30/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionFermionKS.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionFermionKS)

CActionFermionKS::CActionFermionKS()
    : CAction()
{
}


void CActionFermionKS::Initial(CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    //find fermion field
    INT iFieldId = -1;
    param.FetchValueINT(_T("FieldId"), iFieldId);
    m_pFerimionField = dynamic_cast<CFieldFermionKS*>(appGetLattice()->GetFieldById(static_cast<BYTE>(iFieldId)));
    if (NULL == m_pFerimionField)
    {
        appCrucial(_T("CActionFermionKS work with only CFieldFermionKS!\n"));
    }
}

void CActionFermionKS::PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate)
{
    m_pFerimionField->PrepareForHMC(gaugeNum, bosonNum, gaugeFields, bosonFields);
}

/**
* To make it constant, we need to build a few temp fields outside this class
*/
UBOOL CActionFermionKS::CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
    CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
    CFieldGauge* const* stapleFields, ESolverPhase ePhase) const
{
    return m_pFerimionField->CalculateForce(gaugeNum, bosonNum, gaugeFields, bosonFields, gaugeForces, bosonForces, ePhase);
}

DOUBLE CActionFermionKS::Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldGauge* const* stableFields)
{
    //[ (DD)^(-1/4) phi ]^2
    
    //CFieldFermionKSSU3* pPooled = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(static_cast<BYTE>(m_pFerimionField->m_byFieldId)));
    //assert(NULL != pPooled);
    //m_pFerimionField->CopyTo(pPooled);
    //pPooled->D_EN(pGauge);
    //const CLGComplex res = pPooled->Dot(pPooled);

    CFieldFermionKS* pPooled = dynamic_cast<CFieldFermionKS*>(appGetLattice()->GetPooledFieldById(static_cast<BYTE>(m_pFerimionField->m_byFieldId)));
    m_pFerimionField->CopyTo(pPooled);
    pPooled->D_MD(gaugeNum, bosonNum, gaugeFields, bosonFields);
    const cuDoubleComplex res = pPooled->Dot(m_pFerimionField);

    appDetailed(_T("CActionFermionKS : Energy = %f%s%fi\n"), res.x, res.y > 0 ? "+" : " ", res.y);

    pPooled->Return();
    return res.x;
}


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================