//=============================================================================
// FILENAME : TestUpdator.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [01/28/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestUpdateCommon(CParameters& sParam)
{
    Real fExpected = F(0.65);

#if !_CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
#else
    sParam.FetchValueReal(_T("ExpectedResDebug"), fExpected);
#endif

    INT iValue = 0;
#if !_CLG_DEBUG
    iValue = 3;
    sParam.FetchValueINT(_T("BeforeMetropolis"), iValue);
#else
    iValue = 1;
    sParam.FetchValueINT(_T("BeforeMetropolisDebug"), iValue);
#endif
    UINT uiBeforeMetropolis = static_cast<UINT>(iValue);

#if !_CLG_DEBUG
    iValue = 12;
    sParam.FetchValueINT(_T("Metropolis"), iValue);
#else
    iValue = 4;
    sParam.FetchValueINT(_T("MetropolisDebug"), iValue);
#endif
    UINT uiMetropolis = static_cast<UINT>(iValue);

#if !_CLG_DEBUG
    iValue = 2;
    sParam.FetchValueINT(_T("ExpectMiss"), iValue);
#else
    iValue = 1;
    sParam.FetchValueINT(_T("ExpectMissDebug"), iValue);
#endif
    UINT uiExpectAccept = static_cast<UINT>(uiMetropolis - iValue);


    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    //Equilibration
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(uiBeforeMetropolis, FALSE);
#else
    appGetLattice()->m_pUpdator->Update(uiBeforeMetropolis, FALSE);
#endif

    //Measure
    if (NULL != pMeasure)
    {
        pMeasure->Reset();
    }

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(uiMetropolis, TRUE);
#else
    appGetLattice()->m_pUpdator->Update(uiMetropolis, TRUE);
#endif

    pMeasure->Average();
    const Real fRes = pMeasure->GetAverageRealRes();
    appGeneral(_T("res : expected=%f res=%f "), fExpected, fRes);
    UINT uiError = 0;
    if (appAbs(fRes - fExpected) > F(0.15))
    {
        ++uiError;
    }

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount() - uiBeforeMetropolis;
    const Real fHDiff = static_cast<Real>(appGetLattice()->m_pUpdator->GetHDiff());
    const Real fLastHDiff = appAbs(appGetLattice()->m_pUpdator->GetLastHDiff());
#if _CLG_DEBUG
    appGeneral(_T("accept (%d/%d) : expected >= %d. HDiff = %f : expected < 0.3. last HDiff = %f : expected < 0.15\n"), uiAccept, uiMetropolis, uiExpectAccept, fHDiff, fLastHDiff);
#else
    appGeneral(_T("accept (%d/%d) : expected >= %d. HDiff = %f : expected < 0.1 (exp(-0.1)=90%%),  last HDiff = %f : expected < 0.08\n"), uiAccept, uiMetropolis, uiExpectAccept, fHDiff, fLastHDiff);
#endif

#if _CLG_DEBUG
    if (uiAccept < uiExpectAccept)
#else
    if (uiAccept < uiExpectAccept)
#endif
    {
        ++uiError;
    }

#if _CLG_DEBUG
    if (fHDiff > F(0.3))
#else
    if (fHDiff > F(0.1))
#endif
    {
        ++uiError;
    }

#if _CLG_DEBUG
    if (fLastHDiff > F(0.15))
#else
    if (fLastHDiff > F(0.08))
#endif
    {
        ++uiError;
    }

    return uiError;
}

/*
UINT TestUpdator3D(CParameters& sParam)
{
    Real fExpected = F(0.2064);
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);

    //we calculate staple energy from beta = 1 - 6
    //CActionGaugePlaquette* pAction = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
    //if (NULL == pAction)
    //{
    //    return 1;
    //}
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
*/

//__REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorLeapFrog);

//__REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorOmelyan);

__REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorForceGradient, ForceGradient);

___REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorForceGradient3D, ForceGradient3D, _TEST_BOUND);

__REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorTreeImproved, TreeImproved);


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

    CFieldGauge* pGauge = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField[0]->GetCopy());
    CFieldGauge* pStaple = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField[0]->GetCopy());
    TArray<CFieldGauge*> gaugeFields;
    gaugeFields.AddItem(pGauge);

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
            appGetLattice()->m_pGaugeField[0]->CopyTo(pGauge);
            pGauge->CalculateOnlyStaple(pStaple);
            appGetLattice()->m_pGaugeSmearing->GaugeSmearing(pGauge, pStaple);
            pMeasure->OnConfigurationAccepted(1, 0, gaugeFields.GetData(), NULL, NULL);
            iAccepted = newCount;
        }
    }
    
    pMeasure->Report();

    return uiError;
}

__REGIST_TEST(TestWilsonLoop, Updator, TestWilsonLoop, WilsonLoop);

//=============================================================================
// END OF FILE
//=============================================================================
