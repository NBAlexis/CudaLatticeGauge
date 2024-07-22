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
    TArray<Real> expectres;

#if !_CLG_DEBUG
    sParam.FetchValueArrayReal(_T("ExpectedRes"), expectres);
#else
    sParam.FetchValueArrayReal(_T("ExpectedResDebug"), expectres);
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

#if _CLG_DEBUG
    iValue = 0;
    sParam.FetchValueINT(_T("SkipThermalDebug"), iValue);
    UBOOL bSkipThermal = (0 != iValue);
#endif

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

    //Equilibration
#if !_CLG_DEBUG
    appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
#else
    if (bSkipThermal)
    {
        appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
    }
#endif
    appGetLattice()->m_pUpdator->Update(uiBeforeMetropolis, FALSE);
#if !_CLG_DEBUG
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
#else
    if (bSkipThermal)
    {
        appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
    }
#endif
    //Measure
    appGetLattice()->m_pMeasurements->Reset();

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->Update(uiMetropolis, TRUE);

    appGetLattice()->m_pMeasurements->AverageAll();
    const TArray<Real> fRes = appGetLattice()->m_pMeasurements->AverageReals();
    CCString sProblem;
    sProblem.Format(_T("res : expected=%s res=%s "), appToString(expectres).c_str(), appToString(fRes).c_str());
    appGeneral(sProblem);
    UINT uiError = 0;
    if (expectres.Num() != fRes.Num())
    {
        LastProbem(sProblem);
        ++uiError;
    }
    else
    {
        for (INT i = 0; i < fRes.Num(); ++i)
        {
            if (appAbs(fRes[i] - expectres[i]) > F(0.15))
            {
                LastProbem(sProblem);
                ++uiError;
                break;
            }
        }
    }


    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount() - uiBeforeMetropolis;
    const Real fHDiff = static_cast<Real>(appGetLattice()->m_pUpdator->GetHDiff());
    const Real fLastHDiff = appAbs(appGetLattice()->m_pUpdator->GetLastHDiff());
#if _CLG_DEBUG
    const Real fExpectedHDiff = F(0.3);
    const Real fExpectedLastHDiff = F(0.15);
#else
    const Real fExpectedHDiff = F(0.12);
    const Real fExpectedLastHDiff = F(0.1);
#endif

    appGeneral(_T("accept (%d/%d) : expected >= %d. HDiff = %f : expected < %f,  last HDiff = %f : expected < %f\n"), uiAccept, uiMetropolis, uiExpectAccept, fHDiff, fExpectedHDiff, fLastHDiff, fExpectedLastHDiff);

    if (uiAccept < uiExpectAccept)
    {
        sProblem.Format(_T("accept : %d < expect=%d "), uiAccept, uiExpectAccept);
        LastProbem(sProblem);
        ++uiError;
    }

    if (fHDiff > fExpectedHDiff)
    {
        sProblem.Format(_T("hdiff : %f > expect=%f "), fHDiff, fExpectedHDiff);
        LastProbem(sProblem);
        ++uiError;
    }

    if (fLastHDiff > fExpectedLastHDiff)
    {
        sProblem.Format(_T("last-hdiff : %f > expect=%f "), fLastHDiff, fExpectedLastHDiff);
        LastProbem(sProblem);
        ++uiError;
    }

    return uiError;
}

__REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorForceGradient, ForceGradient);
__REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorTreeImproved, TreeImproved);


UINT TestWilsonLoop(CParameters& sParam)
{
    Real fExpected = F(0.2064);
#if _CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedResD"), fExpected);
#else
    sParam.FetchValueReal(_T("ExpectedResR"), fExpected);
#endif

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

    Real fCheck = pMeasure->m_lstAverageC[0][0].x;
    appGeneral(_T("Re[averageC(r=1, t=1)] = %f (expected = %f)\n"), fCheck, fExpected);

    if (abs(fCheck - fExpected) > F(0.001))
    {
        ++uiError;
    }

    return uiError;
}

__REGIST_TEST(TestWilsonLoop, Updator, TestWilsonLoop, WilsonLoop);

//=============================================================================
// END OF FILE
//=============================================================================
