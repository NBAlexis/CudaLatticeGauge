//=============================================================================
// FILENAME : TestUpdator.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [mm/dd/yy]
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
    appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->Update(uiBeforeMetropolis, FALSE);
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);

    //Measure
    appGetLattice()->m_pMeasurements->Reset();

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->ClearHDiffHistory();
    appGetLattice()->m_pUpdator->Update(uiMetropolis, TRUE);

    appGetLattice()->m_pMeasurements->AverageAll();
    const TArray<Real> fRes = appGetLattice()->m_pMeasurements->AverageReals();
    CCString sProblem;
    sProblem.Format(_T("res : expected=%s res=%s \n"), appToString(expectres).c_str(), appToString(fRes).c_str());
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
#if !_CLG_DEBUG
            if (appAbs(fRes[i] - expectres[i]) > F(0.05))
#else
            if (appAbs(fRes[i] - expectres[i]) > F(0.15))
#endif
            {
                LastProbem(sProblem);
                ++uiError;
                break;
            }
        }
    }


    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real fHDiff = static_cast<Real>(appGetLattice()->m_pUpdator->GetHDiff());
    const Real fLastHDiff = appGetLattice()->m_pUpdator->GetLastHDiff();
    const Real fMaxHDiff = appGetLattice()->m_pUpdator->GetMaxHDiff();
#if _CLG_DEBUG
    const Real fExpectedHDiff = F(0.2);
    const Real fExpectedLastHDiff = F(0.12);
    const Real fExpectedMaxHDiff = F(0.3);
#else
    const Real fExpectedHDiff = F(0.1);
    const Real fExpectedLastHDiff = F(0.08);
    const Real fExpectedMaxHDiff = F(0.15);
#endif

    appGeneral(_T("\n accept (%s/%d) : expected >= %d\n HDiff = %s : expected < %f\n  Last HDiff = %s : expected < %f\n  Max HDiff = %s : expected < %f\n"), 
        uiAccept < uiExpectAccept ? appDressColor(EVC_RED, appToString(uiAccept).c_str()).c_str() : appDressColor(EVC_GREEN, appToString(uiAccept).c_str()).c_str(),
        uiMetropolis, 
        uiExpectAccept, 
        fHDiff > fExpectedHDiff ? appDressColor(EVC_RED, appToString(fHDiff).c_str()).c_str() : appDressColor(EVC_GREEN, appToString(fHDiff).c_str()).c_str(),
        fExpectedHDiff, 
        appAbs(fLastHDiff) > fExpectedLastHDiff ? appDressColor(EVC_RED, appToString(fLastHDiff).c_str()).c_str() : appDressColor(EVC_GREEN, appToString(fLastHDiff).c_str()).c_str(),
        fExpectedLastHDiff,
        appAbs(fMaxHDiff) > fExpectedMaxHDiff ? appDressColor(EVC_RED, appToString(fMaxHDiff).c_str()).c_str() : appDressColor(EVC_GREEN, appToString(fMaxHDiff).c_str()).c_str(),
        fExpectedMaxHDiff);

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

    if (appAbs(fLastHDiff) > fExpectedLastHDiff)
    {
        sProblem.Format(_T("last-hdiff : %f > expect=%f "), fLastHDiff, fExpectedLastHDiff);
        LastProbem(sProblem);
        ++uiError;
    }

    if (appAbs(fMaxHDiff) > fExpectedMaxHDiff)
    {
        sProblem.Format(_T("max-hdiff : %f > expect=%f "), fMaxHDiff, fExpectedMaxHDiff);
        LastProbem(sProblem);
        ++uiError;
    }

    return uiError;
}

___REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorLeapFrog, Basic, _TEST_MULTIGPU);
__REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorCylinder, Cylinder);
__REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorCylinderPlaq, CylinderPlaq);
___REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorForceGradient, ForceGradient, _TEST_MULTIGPU);
___REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorForceGradientSU3_12, ForceGradientSU3_12, _TEST_MULTIGPU);
___REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorTreeImproved, TreeImproved, _TEST_MULTIGPU);
__REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorTreeImprovedD, TreeImprovedD);
__REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorTreeImprovedP, TreeImprovedP);
___REGIST_TEST(TestUpdateCommon, Updator, TestUpdatorOneLoopImproved, OneLoopImproved, _TEST_MULTIGPU);


UINT TestHeatbath(CParameters& sParam)
{
    INT iEquil = 5;
    INT iMeasure = 10;
#if _CLG_DEBUG
    sParam.FetchValueINT(_T("EquilStepsD"), iEquil);
    sParam.FetchValueINT(_T("MeasureStepsD"), iMeasure);
#else
    sParam.FetchValueINT(_T("EquilSteps"), iEquil);
    sParam.FetchValueINT(_T("MeasureSteps"), iMeasure);
#endif
    UINT uiEquil = static_cast<UINT>(iEquil);
    UINT uiMeasure = static_cast<UINT>(iMeasure);
    UINT uiError = 0;

    // Equilibration
    appGetLattice()->m_pUpdator->Update(uiEquil, FALSE);

    // Measurement
    appGetLattice()->m_pMeasurements->Reset();
    appGetLattice()->m_pUpdator->Update(uiMeasure, TRUE);
    appGetLattice()->m_pMeasurements->AverageAll();

    const TArray<Real> fRes = appGetLattice()->m_pMeasurements->AverageReals();

    if (fRes.Num() < 4)
    {
        LastProbem(_T("Heatbath measurement returned fewer than 4 results\n"));
        return 1;
    }

    appGeneral(_T("Heatbath energies: Z2(beta=0.5)=%f, Z3(beta=1.0)=%f, D3(beta=1.5)=%f, D4(beta=2.0)=%f\n"),
        fRes[0], fRes[1], fRes[2], fRes[3]);

    TArray<Real> expectres;
#if !_CLG_DEBUG
    sParam.FetchValueArrayReal(_T("ExpectedRes"), expectres);
#else
    sParam.FetchValueArrayReal(_T("ExpectedResDebug"), expectres);
#endif

    if (expectres.Num() > 0)
    {
        if (expectres.Num() != fRes.Num())
        {
            CCString sProblem;
            sProblem.Format(_T("Heatbath expected %d results but got %d\n"), expectres.Num(), fRes.Num());
            LastProbem(sProblem);
            ++uiError;
        }
        else
        {
            for (INT i = 0; i < fRes.Num(); ++i)
            {
                if (appAbs(fRes[i] - expectres[i]) > F(0.2))
                {
                    CCString sProblem;
                    sProblem.Format(_T("Heatbath field %d energy %f differs from expected %f (tol=0.2)\n"), i, fRes[i], expectres[i]);
                    LastProbem(sProblem);
                    ++uiError;
                }
            }
        }
    }
    else
    {
        // Fallback range check when ExpectedRes is not configured in YAML
        for (INT i = 0; i < 4; ++i)
        {
            if (fRes[i] < F(-0.5) || fRes[i] > F(1.5))
            {
                CCString sProblem;
                sProblem.Format(_T("Heatbath field %d energy %f out of range [-0.5, 1.5]\n"), i, fRes[i]);
                LastProbem(sProblem);
                ++uiError;
            }
        }
    }

    return uiError;
}

___REGIST_TEST(TestHeatbath, Updator, TestHeatbath, Heatbath, _TEST_MULTIGPU);


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
            appGetLattice()->m_pGaugeSmearing[pGauge->m_byFieldId]->GaugeSmearing(pGauge, NULL, pStaple);
            pMeasure->OnConfigurationAccepted(1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL);
            iAccepted = newCount;
        }
    }

    pMeasure->Report();

    Real fCheck = pMeasure->m_lstAverageC[0][0].x;
    appGeneral(_T("Re[averageC(r=1, t=1)] = %f (expected = %f)\n"), fCheck, fExpected);

    if (abs(fCheck - fExpected) > F(0.005))
    {
        ++uiError;
    }

    return uiError;
}

__REGIST_TEST(TestWilsonLoop, Updator, TestWilsonLoop, WilsonLoop);

//=============================================================================
// END OF FILE
//=============================================================================
