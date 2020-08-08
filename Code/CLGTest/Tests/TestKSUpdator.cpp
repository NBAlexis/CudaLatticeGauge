//=============================================================================
// FILENAME : TestUpdator.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [01/28/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestFermionUpdatorKS(CParameters& sParam)
{
    Real fExpected = F(0.392);
#if !_CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
#endif
    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    appGetLattice()->m_pUpdator->Update(3, FALSE);

    pMeasure->Reset();
#if !_CLG_DEBUG
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->Update(40, TRUE);
#else
    appGetLattice()->m_pUpdator->Update(5, TRUE);
    Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f"), fExpected, fRes);
    if (appAbs(fRes - fExpected) > F(0.02))
    {
        return 1;
    }
    return 0;
#endif

#if !_CLG_DEBUG
    const Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f\n"), fExpected, fRes);
    UINT uiError = 0;
    if (appAbs(fRes - fExpected) > F(0.01))
    {
        ++uiError;
    }

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
    appGeneral(_T("accept (%d/60) : expected >= 50. HDiff = %f : expected < 0.3\n (exp(-0.3) is 74%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());

    if (uiAccept < 50)
    {
        ++uiError;
    }

    if (fHDiff > F(0.3))
    {
        ++uiError;
    }

    return uiError;
#endif
}

__REGIST_TEST(TestFermionUpdatorKS, Updator, TestFermionUpdatorKS);
__REGIST_TEST(TestFermionUpdatorKS, Updator, TestFermionUpdatorKSNestedForceGradient);

