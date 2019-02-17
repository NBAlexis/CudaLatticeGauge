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
    , m_uiMutiStep(1)
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

    INT iMultiStep = 1;
    if (param.FetchValueINT(_T("MultiStep"), iMultiStep))
    {
        if (iMultiStep > 0)
        {
            m_uiMutiStep = static_cast<UINT>(iMultiStep);
        }
        else
        {
            appCrucial(_T("CActionFermionWilsonNf2 MultiStep must > 0\n"));
        }
    }
}

void CActionFermionWilsonNf2::PrepareForHMC(const CFieldGauge* pGauge, UINT )
{
    m_pFerimionField->PrepareForHMC(pGauge);
}

/**
* To make it constant, we need to build a few temp fields outside this class
*/
UBOOL CActionFermionWilsonNf2::CalculateForceOnGauge(UINT uiStep, const CFieldGauge* pGauge, CFieldGauge* pForce, CFieldGauge *) const
{
    if (m_uiMutiStep > 1)
    {
        CFieldGauge* pCachedForce = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pFieldCache->GetCachedField(m_byActionId + CFieldCache::CachedForceFieldStart));
        if (0 == uiStep % m_uiMutiStep || NULL == pCachedForce)
        {
            if (NULL == pCachedForce)
            {
                pCachedForce = dynamic_cast<CFieldGauge*>(pForce->GetCopy());
                appGetLattice()->m_pFieldCache->CacheField(m_byActionId + CFieldCache::CachedForceFieldStart, pCachedForce);
            }
            return m_pFerimionField->CalculateForce(pGauge, pForce, pCachedForce);
        }
        pForce->AxpyPlus(pCachedForce);
        return TRUE;
    }
    return m_pFerimionField->CalculateForce(pGauge, pForce, NULL);
}

Real CActionFermionWilsonNf2::Energy(UBOOL , const CFieldGauge* pGauge, const CFieldGauge* )
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