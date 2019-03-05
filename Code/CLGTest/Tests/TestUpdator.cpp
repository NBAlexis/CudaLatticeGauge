//=============================================================================
// FILENAME : TestUpdator.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [01/28/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestUpdator(CParameters& sParam)
{
    Real fExpected = F(0.2064);
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);

    //we calculate staple energy from beta = 1 - 6
    CActionGaugePlaquette * pAction = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
    if (NULL == pAction)
    {
        return 1;
    }
    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    //pAction->SetBeta(F(3.0));

    //Equilibration
    appGetLattice()->m_pUpdator->Update(10, FALSE);

    //Measure
    pMeasure->Reset();
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->Update(40, TRUE);

    Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f"), fExpected, fRes);
    UINT uiError = 0;
    if (appAbs(fRes - fExpected) > F(0.005))
    {
        ++uiError;
    }

    UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
    appGeneral(_T("accept (%d/50) : expected >= 45. HDiff = %f : expected < 0.01\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());

    if (uiAccept < 45)
    {
        ++uiError;
    }

    if (fHDiff > F(0.01))
    {
        ++uiError;
    }

    return uiError;
}

__REGIST_TEST(TestUpdator, Updator, TestUpdatorLeapFrog);

__REGIST_TEST(TestUpdator, Updator, TestUpdatorOmelyan);

__REGIST_TEST(TestUpdator, Updator, TestUpdatorForceGradient);


//=============================================================================
// END OF FILE
//=============================================================================
