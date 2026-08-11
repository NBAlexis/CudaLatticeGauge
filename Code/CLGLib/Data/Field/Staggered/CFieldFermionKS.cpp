//=============================================================================
// FILENAME : CFieldFermion.cpp
// 
// DESCRIPTION:
// There are functions for common fermions
//
// REVISION:
//  [12/7/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CFieldFermionKS.h"

__BEGIN_NAMESPACE

void CFieldFermionKS::InitialOtherParameters(CParameters& params)
{
    CFieldFermion::InitialOtherParameters(params);

    params.FetchValueReal(_T("Mass"), m_f2am);
    if (m_f2am < F(0.00000001))
    {
        appWarning(_T("CFieldFermionKS: Mass is nearly 0\n"));
    }

    INT iEachEta = 0;
    params.FetchValueINT(_T("EachSiteEta"), iEachEta);
    m_bEachSiteEta = (0 != iEachEta);

    m_eRational = ER_AllRational;
}

CCString CFieldFermionKS::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldFermion::GetInfos(tab);
    sRet = sRet + tab + _T("Mass (2am) : ") + appToString(m_f2am) + _T("\n");
    sRet = sRet + tab + _T("Diagonal Mass : ") + appToString(m_bDiagonalMass) + _T("\n");
    sRet = sRet + tab + _T("Each site eta : ") + appToString(m_bEachSiteEta) + _T("\n");
    return sRet;
}

DOUBLE CFieldFermionKS::EnergyS(const CFieldGauge* pGauge) const
{
    CFieldFermion* pPooled = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(static_cast<BYTE>(m_byFieldId), _T(__FILE__), __LINE__));
    CopyTo(pPooled);
    const CFieldGauge* gaugeFields[1] = { pGauge };
    pPooled->D_MD(1, 0, 0, gaugeFields, NULL, NULL);
    const cuDoubleComplex res = pPooled->Dot(this);
    pPooled->Return();
    return res.x;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================