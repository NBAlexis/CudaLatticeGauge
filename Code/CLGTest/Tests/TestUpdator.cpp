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

    const Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f"), fExpected, fRes);
    UINT uiError = 0;
    if (appAbs(fRes - fExpected) > F(0.005))
    {
        ++uiError;
    }

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real fHDiff = static_cast<Real>(appGetLattice()->m_pUpdator->GetHDiff());
    appGeneral(_T("accept (%d/50) : expected >= 45. HDiff = %f : expected < 0.1 (exp(-0.1)=90%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());

    if (uiAccept < 45)
    {
        ++uiError;
    }

    if (fHDiff > F(0.1))
    {
        ++uiError;
    }

    return uiError;
}

__REGIST_TEST(TestUpdator, Updator, TestUpdatorLeapFrog);

__REGIST_TEST(TestUpdator, Updator, TestUpdatorOmelyan);

__REGIST_TEST(TestUpdator, Updator, TestUpdatorForceGradient);


UINT TestWilsonLoop(CParameters& sParam)
{
    Real fExpected = F(0.2064);
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);

    //we calculate staple energy from beta = 1 - 6
    CActionGaugePlaquette* pAction = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
    if (NULL == pAction)
    {
        return 1;
    }
    CMeasureWilsonLoop* pMeasure = dynamic_cast<CMeasureWilsonLoop*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }
    UINT uiError = 0;
    //pAction->SetBeta(F(3.0));

    CFieldGauge* pGauge = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField->GetCopy());
    CFieldGauge* pStaple = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField->GetCopy());

    //Equilibration
    //appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(10, FALSE);
#else
    appGetLattice()->m_pUpdator->Update(20, FALSE);
#endif
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
    //Measure
    pMeasure->Reset();
    appGetLattice()->m_pUpdator->SetTestHdiff(FALSE);
    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    INT iAccepted = appGetLattice()->m_pUpdator->GetConfigurationCount();
#if _CLG_DEBUG
    while (iAccepted < 20)
#else
    while (iAccepted < 50)
#endif
    {
        const INT newCount = appGetLattice()->m_pUpdator->Update(1, FALSE);
        if (newCount != iAccepted)
        {
            appGetLattice()->m_pGaugeField->CopyTo(pGauge);
            pGauge->CalculateOnlyStaple(pStaple);
            appGetLattice()->m_pGaugeSmearing->GaugeSmearing(pGauge, pStaple);
            pMeasure->OnConfigurationAccepted(pGauge, NULL);
            iAccepted = newCount;
        }
    }
    
    pMeasure->Report();
    //const Real fRes = pMeasure->m_fLastRealResult;
    //appGeneral(_T("res : expected=%f res=%f"), fExpected, fRes);
    //UINT uiError = 0;
    //if (appAbs(fRes - fExpected) > F(0.005))
    //{
    //    ++uiError;
    //}

    //const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    //const Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
    //appGeneral(_T("accept (%d/50) : expected >= 45. HDiff = %f : expected < 0.1 (exp(-0.1)=90%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());

    //if (uiAccept < 45)
    //{
    //    ++uiError;
    //}

    //if (fHDiff > F(0.1))
    //{
    //    ++uiError;
    //}

    return uiError;
}

__REGIST_TEST(TestWilsonLoop, Updator, TestWilsonLoop);

//=============================================================================
// END OF FILE
//=============================================================================
