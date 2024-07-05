//=============================================================================
// FILENAME : CIntegratorOmelyan.cpp
// 
// DESCRIPTION:
// This is the Omelyan integrator for HMC
//
// REVISION:
//  [02/12/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CIntegratorOmelyan.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CIntegratorOmelyan)

void CIntegratorOmelyan::Initial(CHMC* pOwner, CLatticeData* pLattice, const CParameters& params)
{
    CIntegrator::Initial(pOwner, pLattice, params);
    m_f2Lambda = OmelyanLambda2;
    if (!params.FetchValueReal(_T("Omelyan2Lambda"), m_f2Lambda))
    {
        m_f2Lambda = OmelyanLambda2;
    }
}

CCString CIntegratorOmelyan::GetInfos(const CCString& sTab) const
{
    CCString sRet;
    sRet = sTab + _T("Name : Omelyan\n");
    sRet = sRet + sTab + _T("Epsilon : ") + appToString(m_fEStep) + _T("\n");
    sRet = sRet + sTab + _T("Step : ") + appToString(static_cast<INT>(m_uiStepCount)) + _T("\n");
    sRet = sRet + sTab + _T("##Tau is trajectory length = Epsilon x Step\n");
    sRet = sRet + sTab + _T("Tau : ") + appToString(m_fEStep * m_uiStepCount) + _T("\n");
    sRet = sRet + sTab + _T("##Omelyan2Lambda = 2 x lambda\n");
    sRet = sRet + sTab + _T("Omelyan2Lambda : ") + appToString(m_f2Lambda) + _T("\n");
    return sRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================