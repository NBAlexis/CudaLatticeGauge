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

void CActionFermionWilsonNf2::PrepareForHMC(const CFieldGauge* pGauge, UINT )
{
    m_pFerimionField->PrepareForHMC(pGauge);
}

/**
* To make it constant, we need to build a few temp fields outside this class
*/
UBOOL CActionFermionWilsonNf2::CalculateForceOnGauge(const CFieldGauge* pGauge, CFieldGauge* pForce, CFieldGauge * /*staple*/, ESolverPhase ePhase) const
{
    return m_pFerimionField->CalculateForce(pGauge, pForce, ePhase);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CActionFermionWilsonNf2::Energy(UBOOL, const CFieldGauge* pGauge, const CFieldGauge*)
#else
Real CActionFermionWilsonNf2::Energy(UBOOL , const CFieldGauge* pGauge, const CFieldGauge* )
#endif
{
    //(D^-1 phi)^2
    CFieldFermionWilsonSquareSU3* pPooled = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(static_cast<BYTE>(m_pFerimionField->m_byFieldId)));
    assert(NULL != pPooled);
    m_pFerimionField->CopyTo(pPooled);
    pPooled->InverseD(pGauge);
#if !_CLG_DOUBLEFLOAT
    const cuDoubleComplex res = pPooled->Dot(pPooled);
#else
    const CLGComplex res = pPooled->Dot(pPooled);
#endif
    //pPooled->InverseDDdagger(pGauge);
    //const CLGComplex res = m_pFerimionField->Dot(pPooled);
    appDetailed(_T("CActionFermionWilsonNf2 : Energy = %f%s%fi\n"), res.x, res.y > 0 ? "+" : " ", res.y);

    pPooled->Return();
    return res.x;
}

CCString CActionFermionWilsonNf2::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CActionFermionWilsonNf2\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================