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


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionFermionKS)

CActionFermionKS::CActionFermionKS()
    : CAction()
{
}


void CActionFermionKS::Initial(CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    m_pOwner = pOwner;
    m_byActionId = byId;

    //find fermion field
    INT iFieldId = -1;
    param.FetchValueINT(_T("FieldId"), iFieldId);
    m_pFerimionField = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetFieldById(static_cast<BYTE>(iFieldId)));
    if (NULL == m_pFerimionField)
    {
        appCrucial(_T("CActionFermionKS work with only CFieldFermionWilsonSquareSU3!\n"));
    }
}

void CActionFermionKS::PrepareForHMC(const CFieldGauge* pGauge, UINT )
{
    m_pFerimionField->PrepareForHMC(pGauge);
}

/**
* To make it constant, we need to build a few temp fields outside this class
*/
UBOOL CActionFermionKS::CalculateForceOnGauge(const CFieldGauge* pGauge, CFieldGauge* pForce, CFieldGauge * /*staple*/, ESolverPhase ePhase) const
{
    return m_pFerimionField->CalculateForce(pGauge, pForce, ePhase);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CActionFermionKS::Energy(UBOOL, const CFieldGauge* pGauge, const CFieldGauge*)
#else
Real CActionFermionKS::Energy(UBOOL , const CFieldGauge* pGauge, const CFieldGauge* )
#endif
{
    //[ (DD)^(-1/4) phi ]^2
    
    //CFieldFermionKSSU3* pPooled = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(static_cast<BYTE>(m_pFerimionField->m_byFieldId)));
    //assert(NULL != pPooled);
    //m_pFerimionField->CopyTo(pPooled);
    //pPooled->D_EN(pGauge);
    //const CLGComplex res = pPooled->Dot(pPooled);

    CFieldFermionKSSU3* pPooled = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(static_cast<BYTE>(m_pFerimionField->m_byFieldId)));
    m_pFerimionField->CopyTo(pPooled);
    pPooled->D_MD(pGauge);
#if !_CLG_DOUBLEFLOAT
    const cuDoubleComplex res = pPooled->Dot(m_pFerimionField);
#else
    const CLGComplex res = pPooled->Dot(m_pFerimionField);
#endif

    appDetailed(_T("CActionFermionKS : Energy = %f%s%fi\n"), res.x, res.y > 0 ? "+" : " ", res.y);

    pPooled->Return();
    return res.x;
}

CCString CActionFermionKS::GetInfos(const CCString &tab) const
{
    CCString sRet = tab + _T("Name : CFieldFermionKSSU3\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================