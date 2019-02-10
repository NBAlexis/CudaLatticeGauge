//=============================================================================
// FILENAME : CActionFermionWilsonNf2.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [02/06/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionFermionWilsonNf2)

CActionFermionWilsonNf2::CActionFermionWilsonNf2()
    : CAction()
{
}

void CActionFermionWilsonNf2::Initial(CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    m_pOwner = pOwner;
    m_byActionId = byId;

    //find fermion field
    INT iFieldId = -1;
    param.FetchValueINT(_T("FieldId"), iFieldId);
    m_pFerimionField = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetFieldById(static_cast<BYTE>(iFieldId)));
    if (NULL == m_pFerimionField)
    {
        appCrucial(_T("CActionFermionWilsonNf2 work with only CFieldFermionWilsonSquareSU3!\n"));
    }
}

void CActionFermionWilsonNf2::PrepareForHMC(const CFieldGauge* pGauge)
{
    m_pFerimionField->PrepareForHMC(pGauge);
}

/**
* To make it constant, we need to build a few temp fields outside this class
*/
UBOOL CActionFermionWilsonNf2::CalculateForceOnGauge(const CFieldGauge* pGauge, CFieldGauge* pForce, CFieldGauge *) const
{
    return m_pFerimionField->CalculateForce(pGauge, pForce);
}

Real CActionFermionWilsonNf2::Energy(const CFieldGauge* pGauge) const
{
    //(D^-1 phi)^2
    CFieldFermionWilsonSquareSU3* pPooled = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(static_cast<BYTE>(m_pFerimionField->m_byFieldId)));
    assert(NULL != pPooled);
    m_pFerimionField->CopyTo(pPooled);
    pPooled->InverseD(pGauge);
    _Complex res = pPooled->Dot(pPooled);
    appGeneral(_T("CActionFermionWilsonNf2 : Energy = %f%s%fi\n"), res.x, res.y > 0 ? "+" : " ", res.y);

    pPooled->Return();
    return res.x;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================