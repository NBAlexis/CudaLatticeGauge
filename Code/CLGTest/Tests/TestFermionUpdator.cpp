//=============================================================================
// FILENAME : TestFermionUpdator.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [mm/dd/yy]
//  [02/07/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestFermionUpdator(CParameters& sParam)
{
    Real fExpected = F(0.17);
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    appGetLattice()->m_pUpdator->Update(20, FALSE);

    pMeasure->Reset();
#if !_CLG_DEBUG
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->Update(40, TRUE);
#else
    appGetLattice()->m_pUpdator->Update(10, TRUE);
    pMeasure->Average();
    Real fRes = pMeasure->GetAverageRealRes();
    appGeneral(_T("res : expected=%f res=%f"), fExpected, fRes);
    if (appAbs(fRes - fExpected) > F(0.02))
    {
        return 1;
    }
    return 0;
#endif
    
#if !_CLG_DEBUG
    const Real fRes = pMeasure->GetLastRealRes();
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

//__REGIST_TEST(TestUpdateCommon, Updator, TestFermionUpdator, FermionUpdator);

__REGIST_TEST(TestUpdateCommon, WDUpdator, TestFermionUpdatorOmelyanGCRODR, FermionOmelyanGCRODR);

__REGIST_TEST(TestUpdateCommon, WDUpdator, TestFermionUpdatorOmelyanGMRESMDR, FermionOmelyanGMRESMDR);

//__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdatorOmelyan, FermionOmelyan);

//__REGIST_TEST(TestUpdateCommon, Updator, TestFermionUpdatorForceGradient, WDForceGradient);

//__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdatorNestedLeapFrog, FermionNestedLeapFrog);

//REGIST_TEST(TestUpdateCommon, Updator, TestFermionUpdatorNestedOmelyan, FermionNestedOmelyan);

__REGIST_TEST(TestUpdateCommon, WDUpdator, TestFermionUpdatorNestedForceGradient, NestedForceGradient);
___REGIST_TEST(TestUpdateCommon, WDUpdator, TestFermionUpdatorCG, FermionCG, _TEST_MULTIGPU);
___REGIST_TEST(TestUpdateCommon, WDUpdator, TestFermionUpdatorDeflatedCG, FermionDCG, _TEST_MULTIGPU);

__REGIST_TEST(TestUpdateCommon, WDUpdator, TestFermionUpdatorNestedForceGradientRHMC, WDRHMC);

//#if !_CLG_DEBUG

__REGIST_TEST(TestUpdateCommon, WDUpdator, TestFermionUpdatorWilsonDiracGamma1, WilsonDiracExpGamma);

__REGIST_TEST(TestUpdateCommon, WDUpdator, TestFermionUpdatorWilsonDiracGamma2, WilsonDiracGamma);

__REGIST_TEST(TestUpdateCommon, WDUpdator, TestFermionUpdatorWilsonDiracEM, WilsonDiracEM);

__REGIST_TEST(TestUpdateCommon, WDUpdator, TestFermionUpdatorCloverWilsonNf2p1, CloverNf2p1);

__REGIST_TEST(TestUpdateCommon, WDUpdator, TestFermionUpdatorStoutLinkCloverWilsonNf2p1, StoutNf2p1);

__REGIST_TEST(TestUpdateCommon, WDUpdator, TestFermionUpdatorStoutLinkCloverWilsonEMNf2, StoutEMNf2);

__REGIST_TEST(TestUpdateCommon, Boundary, TestFermionUpdatorStoutLinkCloverWilsonDNf2p1, StoutDirichNf2p1);

//#endif

UINT TestFermionUpdatorWithMesonCorrelator(CParameters& sParam)
{
    CMeasureMesonCorrelator* pMeasure = dynamic_cast<CMeasureMesonCorrelator*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    UINT uiError = 0;

#if _CLG_DEBUG

    TArray<Real> lstResExpected;
    sParam.FetchValueArrayReal(_T("ExpectedRes"), lstResExpected);
    assert(static_cast<UINT>(lstResExpected.Num()) == _HC_Lt - 1);
    appGetLattice()->m_pUpdator->Update(10, FALSE);
    appGetLattice()->m_pUpdator->Update(20, TRUE);

    TArray<Real> lstRes;
    appGeneral(_T("res = expected vs test: \n"));
    
    for (UINT i = 1; i < _HC_Lt; ++i)
    {
#if _CLG_DOUBLEFLOAT
        Real fRes = _hostlog10(pMeasure->m_lstResults[0][i] / pMeasure->m_lstResults[0][0]);
#else
        Real fRes = static_cast<Real>(_hostlog10d(pMeasure->m_lstResults[0][i] / pMeasure->m_lstResults[0][0]));
#endif
        appGeneral(_T("%f : %f \n"), lstResExpected[i - 1], fRes);
        if (appAbs(fRes - lstResExpected[i - 1]) > F(0.2))
        {
            ++uiError;
        }
    }
    appGeneral(_T("\n"));

#else

    Real fExpected = F(0.625);
    sParam.FetchValueReal(_T("ExpectedResR"), fExpected);
    appGetLattice()->m_pUpdator->Update(20, FALSE);
    appGetLattice()->m_pUpdator->Update(50, TRUE);

    const DOUBLE fRes = pMeasure->m_lstResults[0][0];
    appGeneral(_T("res : expected=%f res=%f\n"), fExpected, fRes);
    if (appAbs(fRes - fExpected) > F(0.005))
    {
        return 1;
    }

#endif
    return uiError;
}

#if _CLG_DEBUG
__REGIST_TEST(TestFermionUpdatorWithMesonCorrelator, WDUpdator, TestFermionUpdatorWithMesonCorrelator, FermionMesonCorrelator);
#else
__REGIST_TEST(TestFermionUpdatorWithMesonCorrelator, WDUpdator, TestFermionUpdatorWithMesonCorrelatorRelease, FermionMesonCorrelator);
#endif

//Why I cannot find 'TestGaugeSmearingAPEProj'?

__REGIST_TEST(TestFermionUpdatorWithMesonCorrelator, WDUpdator, TestGaugeSmearingAPEProj, GaugeSmearingAPEProj);

__REGIST_TEST(TestFermionUpdatorWithMesonCorrelator, WDUpdator, TestGaugeSmearingAPEStout, GaugeSmearingAPEStout);

UINT TestFermionUpdatorL(CParameters& sParam)
{
    Real fExpected = F(0.17);
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    INT updates = 3;
    sParam.FetchValueINT(_T("Updates"), updates);

    if (updates > 10 || updates < 3)
    {
        updates = 3;
    }
    appGetLattice()->m_pUpdator->Update(1, FALSE);

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    pMeasure->Reset();
    appGetLattice()->m_pUpdator->Update(static_cast<UINT>(updates - 1), TRUE);

    pMeasure->Average();
    const Real fRes = pMeasure->GetAverageRealRes();
    const Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
    const Real fH = appGetLattice()->m_pUpdator->GetHValue();

#if _CLG_DOUBLEFLOAT
    appGeneral(_T("res : expected=%f res=%f |r-e| should be < 0.00001\n"), fExpected, fRes);
    UINT uiError = 0;
    if (appAbs(fRes - fExpected) > F(0.00001))
    {
        ++uiError;
    }

    appGeneral(_T("H = %f, HDiff = %1.8f : expected < 1E-8 H (which is %f)\n"), fH, F(0.00000001) * fH, fHDiff);

    if (fHDiff > F(0.00000001) * fH)
    {
        ++uiError;
    }
#else
    appGeneral(_T("res : expected=%f res=%f |r-e| should be < 0.01\n"), fExpected, fRes);
    UINT uiError = 0;
    if (appAbs(fRes - fExpected) > F(0.01))
    {
        ++uiError;
    }

    appGeneral(_T("H = %f, HDiff = %1.8f : expected < 3E-7 H (which is %f)\nThe error can be 1E-7 H, which > 1, so here we do NOT use 0.3 to judge.\n"), fH, F(0.0000003) * fH, fHDiff);

    if (fHDiff > F(0.0000003) * fH)
    {
        ++uiError;
    }
#endif

    return uiError;
}

___REGIST_TEST(TestFermionUpdatorL, WDUpdator, TestFermionUpdatorLargeScale, LargeScale, _TEST_RELEASE);
___REGIST_TEST(TestFermionUpdatorL, WDUpdator, TestFermionUpdatorLargeScaleFloat, LargeScaleFloat, _TEST_RELEASE | _TEST_SINGLE);

static UINT CheckStaggeredMesonMeasureShape(const CMeasureMesonCorrelatorStaggered* pMeasure)
{
    UINT uiError = 0;
    const INT nt = _HC_Lti;
    const INT typeCount = CMeasureMesonCorrelatorStaggered::_kMesonCorrelatorType;
    const INT nConf = pMeasure->m_lstW2WCombinedCorrelator.Num();

    if (NULL == pMeasure->m_pW2WPArray || NULL == pMeasure->m_pP2PPArray)
    {
        appGeneral(_T("CMeasureMesonCorrelatorStaggered p arrays are not allocated.\n"));
        ++uiError;
    }

    if (nConf <= 0)
    {
        appGeneral(_T("CMeasureMesonCorrelatorStaggered has no measured configurations.\n"));
        return uiError + 1;
    }

    if (pMeasure->m_lstP2PCombinedCorrelator.Num() != nConf
     || pMeasure->m_lstW2WCorrelator.Num() != nConf
     || pMeasure->m_lstP2PCorrelator.Num() != nConf)
    {
        appGeneral(_T("CMeasureMesonCorrelatorStaggered configuration count mismatch.\n"));
        return uiError + 1;
    }

    if (pMeasure->m_lstAverageResults.Num() != typeCount)
    {
        appGeneral(_T("CMeasureMesonCorrelatorStaggered average result type count mismatch.\n"));
        ++uiError;
    }

    for (INT conf = 0; conf < nConf; ++conf)
    {
        if (pMeasure->m_lstW2WCombinedCorrelator[conf].Num() != typeCount
         || pMeasure->m_lstP2PCombinedCorrelator[conf].Num() != typeCount
         || pMeasure->m_lstW2WCorrelator[conf].Num() != typeCount
         || pMeasure->m_lstP2PCorrelator[conf].Num() != typeCount)
        {
            appGeneral(_T("CMeasureMesonCorrelatorStaggered type count mismatch at conf %d.\n"), conf);
            ++uiError;
            continue;
        }

        for (INT ty = 0; ty < typeCount; ++ty)
        {
            if (pMeasure->m_lstW2WCombinedCorrelator[conf][ty].Num() != nt
             || pMeasure->m_lstP2PCombinedCorrelator[conf][ty].Num() != nt)
            {
                appGeneral(_T("CMeasureMesonCorrelatorStaggered combined Nt mismatch at conf %d type %d.\n"), conf, ty);
                ++uiError;
            }

            const INT subCount = pMeasure->m_nSubChannels[ty];
            if (pMeasure->m_lstW2WCorrelator[conf][ty].Num() != subCount
             || pMeasure->m_lstP2PCorrelator[conf][ty].Num() != subCount)
            {
                appGeneral(_T("CMeasureMesonCorrelatorStaggered sub-channel count mismatch at conf %d type %d.\n"), conf, ty);
                ++uiError;
                continue;
            }

            for (INT sub = 0; sub < subCount; ++sub)
            {
                if (pMeasure->m_lstW2WCorrelator[conf][ty][sub].Num() != nt
                 || pMeasure->m_lstP2PCorrelator[conf][ty][sub].Num() != nt)
                {
                    appGeneral(_T("CMeasureMesonCorrelatorStaggered sub-channel Nt mismatch at conf %d type %d sub %d.\n"), conf, ty, sub);
                    ++uiError;
                }
            }
        }
    }

    for (INT ty = 0; ty < pMeasure->m_lstAverageResults.Num(); ++ty)
    {
        if (pMeasure->m_lstAverageResults[ty].Num() != nt)
        {
            appGeneral(_T("CMeasureMesonCorrelatorStaggered average Nt mismatch at type %d.\n"), ty);
            ++uiError;
        }
    }

    return uiError;
}

UINT TestFermionUpdatorWithMesonCorrelatorStaggered(CParameters& sParam)
{
    CMeasureMesonCorrelatorStaggered* pMeasure = dynamic_cast<CMeasureMesonCorrelatorStaggered*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));

    //TODO, the CMeasureMesonCorrelatorStaggeredSimple is removed
    CMeasureMesonCorrelatorStaggeredSimple2* pMeasuresimple = dynamic_cast<CMeasureMesonCorrelatorStaggeredSimple2*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    if (NULL == pMeasure)
    {
        return 1;
    }

    if (NULL == pMeasuresimple)
    {
        return 1;
    }

    UINT uiError = 0;
    Real fExpected1 = F(9.0);
    Real fExpected2 = F(-24.0);
    Real fExpected3 = F(9.0);
    Real fExpected4 = F(-24.0);
#if _CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedRes1D"), fExpected1);
    sParam.FetchValueReal(_T("ExpectedRes2D"), fExpected2);
    sParam.FetchValueReal(_T("ExpectedRes3D"), fExpected3);
    sParam.FetchValueReal(_T("ExpectedRes4D"), fExpected4);
#else
    sParam.FetchValueReal(_T("ExpectedRes1R"), fExpected1);
    sParam.FetchValueReal(_T("ExpectedRes2R"), fExpected2);
    sParam.FetchValueReal(_T("ExpectedRes3R"), fExpected3);
    sParam.FetchValueReal(_T("ExpectedRes4R"), fExpected4);
#endif

#if _CLG_DEBUG

    //TArray<Real> lstResExpected;
    //sParam.FetchValueArrayReal(_T("ExpectedRes"), lstResExpected);
    //assert(static_cast<UINT>(lstResExpected.Num()) == _HC_Lt - 1);
    appGetLattice()->m_pUpdator->Update(1, FALSE);
    appGetLattice()->m_pUpdator->Update(8, TRUE);
    //TArray<Real> lstRes;
    //appGeneral(_T("res = expected vs test: "));

    //for (UINT i = 1; i < _HC_Lt; ++i)
    //{
    //    Real fRes = _hostlog10(pMeasure->m_lstResults[0][i] / pMeasure->m_lstResults[0][0]);
    //    appGeneral(_T("%f : %f, "), lstResExpected[i - 1], fRes);
    //    if (appAbs(fRes - lstResExpected[i - 1]) > F(0.2))
    //    {
    //        ++uiError;
    //    }
    //}
    //appGeneral(_T("\n"));

    uiError += CheckStaggeredMesonMeasureShape(pMeasure);
    if (uiError > 0)
    {
        return uiError;
    }

    appGeneral(_T("check1 pMeasure m_lstAverageResults[0][0] = %f, expected: %f\n"), pMeasure->m_lstAverageResults[0][0], fExpected1);
    appGeneral(_T("check2 pMeasure m_lstAverageResults[1][1] = %f, expected: %f\n"), pMeasure->m_lstAverageResults[1][1], fExpected2);
    appGeneral(_T("check1 pMeasuresimple m_lstAverageResults[0][0] = %f, expected: %f\n"), pMeasuresimple->m_lstAverageResults[0][0], fExpected3);
    appGeneral(_T("check2 pMeasuresimple m_lstAverageResults[1][1] = %f, expected: %f\n"), pMeasuresimple->m_lstAverageResults[1][1], fExpected4);

    if (appAbs(pMeasure->m_lstAverageResults[0][0] - fExpected1) > F(2.0))
    {
        ++uiError;
    }
    if (appAbs(pMeasure->m_lstAverageResults[1][1] - fExpected2) > F(2.0))
    {
        ++uiError;
    }
    if (appAbs(pMeasuresimple->m_lstAverageResults[0][0] - fExpected3) > F(2.0))
    {
        ++uiError;
    }
    if (appAbs(pMeasuresimple->m_lstAverageResults[1][1] - fExpected4) > F(2.0))
    {
        ++uiError;
    }

#else

    //Real fExpected = F(0.625);
    //sParam.FetchValueReal(_T("ExpectedResR"), fExpected);
    //appGetLattice()->m_pUpdator->Update(20, FALSE);
    //appGetLattice()->m_pUpdator->Update(50, TRUE);

    //const Real fRes = pMeasure->m_lstResults[0][0];
    //appGeneral(_T("res : expected=%f res=%f"), fExpected, fRes);
    //if (appAbs(fRes - fExpected) > F(0.005))
    //{
    //    return 1;
    //}
    appGetLattice()->m_pUpdator->Update(30, FALSE);
    appGetLattice()->m_pUpdator->Update(100, TRUE);
    //appGeneral(_T("res1=%f\n"), pMeasure->m_lstAverageResults[0][0]);
    //appGeneral(_T("res2=%f\n"), pMeasure->m_lstAverageResults[1][1]);
    uiError += CheckStaggeredMesonMeasureShape(pMeasure);
    if (uiError > 0)
    {
        return uiError;
    }

    appGeneral(_T("check1 pMeasure m_lstAverageResults[0][0] = %f, expected: %f\n"), pMeasure->m_lstAverageResults[0][0], fExpected1);
    appGeneral(_T("check2 pMeasure m_lstAverageResults[1][1] = %f, expected: %f\n"), pMeasure->m_lstAverageResults[1][1], fExpected2);
    appGeneral(_T("check1 pMeasuresimple m_lstAverageResults[0][0] = %f, expected: %f\n"), pMeasuresimple->m_lstAverageResults[0][0], fExpected3);
    appGeneral(_T("check2 pMeasuresimple m_lstAverageResults[1][1] = %f, expected: %f\n"), pMeasuresimple->m_lstAverageResults[1][1], fExpected4);

    if (appAbs(pMeasure->m_lstAverageResults[0][0] - fExpected1) > F(0.1))
    {
        ++uiError;
    }
    if (appAbs(pMeasure->m_lstAverageResults[1][1] - fExpected2) > F(0.1))
    {
        ++uiError;
    }
    if (appAbs(pMeasuresimple->m_lstAverageResults[0][0] - fExpected3) > F(0.1))
    {
        ++uiError;
    }
    if (appAbs(pMeasuresimple->m_lstAverageResults[1][1] - fExpected4) > F(0.1))
    {
        ++uiError;
    }

#endif
    return uiError;
}

__REGIST_TEST(TestFermionUpdatorWithMesonCorrelatorStaggered, UpdatorKS, TestFermionUpdatorWithMesonCorrelatorStaggered, FermionMesonCorrelatorStaggered);

UINT TestBerryPhase(CParameters& sParam)
{
    Real fExpected = F(0.2064);
#if _CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedResD"), fExpected);
#else
    sParam.FetchValueReal(_T("ExpectedResR"), fExpected);
#endif

    CMeasureBerryPhase* pMeasure = dynamic_cast<CMeasureBerryPhase*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    //appGetLattice()->m_pGaugeField->InitialField(EFIT_Identity);
    //pMeasure->m_bGuageFixing = FALSE;
    //pMeasure->OnConfigurationAccepted(appGetLattice()->m_pGaugeField);

    //TArray<CFieldGauge*> gauge;
    //gauge.AddItem(appGetLattice()->m_pGaugeField);

#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(5, FALSE);
#else
    appGetLattice()->m_pUpdator->Update(5, FALSE);
#endif
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
    //Measure
    pMeasure->Reset();
    appGetLattice()->m_pUpdator->SetTestHdiff(FALSE);
    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    INT iAccepted = appGetLattice()->m_pUpdator->GetConfigurationCount();
#if _CLG_DEBUG
    while (iAccepted < 5)
#else
    while (iAccepted < 10)
#endif
    {
        const INT newCount = appGetLattice()->m_pUpdator->Update(1, FALSE);
        if (newCount != iAccepted)
        {
            pMeasure->OnConfigurationAccepted(_FIELDS, NULL);
            iAccepted = newCount;
        }
    }

    pMeasure->Report();

    Real fCheck = static_cast<Real>(pMeasure->m_lstData[pMeasure->GetConfigurationCount() - 1][0]);
    appGeneral(_T("Berry phase of last configuration at t=0, expected: %f, res: %f\n"), fExpected, fCheck);
    if (abs(fCheck - fExpected) > F(0.0001))
    {
        return 1;
    }

    return 0;
}

__REGIST_TEST(TestBerryPhase, Updator, TestBerryPhase, BerryPhase);

//=============================================================================
// END OF FILE
//=============================================================================
