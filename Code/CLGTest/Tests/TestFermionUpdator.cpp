//=============================================================================
// FILENAME : TestFermionUpdator.cpp
// 
// DESCRIPTION:
//
// REVISION:
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

__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdator);

__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdatorOmelyanGCRODR);

__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdatorOmelyanGMRESMDR);

__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdatorOmelyan);

__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdatorForceGradient);

__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdatorNestedLeapFrog);

__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdatorNestedOmelyan);

__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdatorNestedForceGradient);

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
    appGeneral(_T("res = expected vs test: "));
    
    for (UINT i = 1; i < _HC_Lt; ++i)
    {
        Real fRes = _hostlog10(pMeasure->m_lstResults[0][i] / pMeasure->m_lstResults[0][0]);
        appGeneral(_T("%f : %f, "), lstResExpected[i - 1], fRes);
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

    const Real fRes = pMeasure->m_lstResults[0][0];
    appGeneral(_T("res : expected=%f res=%f"), fExpected, fRes);
    if (appAbs(fRes - fExpected) > F(0.005))
    {
        return 1;
    }

#endif
    return uiError;
}

#if _CLG_DEBUG
__REGIST_TEST(TestFermionUpdatorWithMesonCorrelator, Updator, TestFermionUpdatorWithMesonCorrelator);
#else
__REGIST_TEST(TestFermionUpdatorWithMesonCorrelator, Updator, TestFermionUpdatorWithMesonCorrelatorRelease);
#endif

//Why I cannot find 'TestGaugeSmearingAPEProj'?

__REGIST_TEST(TestFermionUpdatorWithMesonCorrelator, Updator, TestGaugeSmearingAPEProj);

__REGIST_TEST(TestFermionUpdatorWithMesonCorrelator, Updator, TestGaugeSmearingAPEStout);

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

    const Real fRes = pMeasure->m_fLastRealResult;
    const Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
    const Real fH = appGetLattice()->m_pUpdator->GetHValue();

    appGeneral(_T("res : expected=%f res=%f\n"), fExpected, fRes);
    UINT uiError = 0;
    if (appAbs(fRes - fExpected) > F(0.01))
    {
        ++uiError;
    }

    appGeneral(_T("H = %f, HDiff = %1.8f : expected < 3E-7 H\nThe error can be 1E-7 H, which > 1, so here we do NOT use 0.3 to judge.\n"), fH, fHDiff);

    if (fHDiff > F(0.0000003) * fH)
    {
        ++uiError;
    }

    return uiError;
}

#if !_CLG_DEBUG
__REGIST_TEST(TestFermionUpdatorL, Updator, TestFermionUpdatorLargeScale);
#endif

UINT TestFermionUpdatorWithMesonCorrelatorStaggered(CParameters& sParam)
{
    CMeasureMesonCorrelatorStaggered* pMeasure = dynamic_cast<CMeasureMesonCorrelatorStaggered*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    UINT uiError = 0;
    Real fExpected1 = F(9.0);
    sParam.FetchValueReal(_T("ExpectedRes1"), fExpected1);
    Real fExpected2 = F(-24.0);
    sParam.FetchValueReal(_T("ExpectedRes2"), fExpected2);

#if _CLG_DEBUG

    //TArray<Real> lstResExpected;
    //sParam.FetchValueArrayReal(_T("ExpectedRes"), lstResExpected);
    //assert(static_cast<UINT>(lstResExpected.Num()) == _HC_Lt - 1);
    appGetLattice()->m_pUpdator->Update(1, FALSE);
    appGetLattice()->m_pUpdator->Update(3, TRUE);
    appGeneral(_T("res1=%f\n"), pMeasure->m_lstAverageResults[0][0]);
    appGeneral(_T("res2=%f\n"), pMeasure->m_lstAverageResults[1][1]);
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
    if (appAbs(pMeasure->m_lstAverageResults[0][0] - fExpected1) > F(2.0))
    {
        ++uiError;
    }
    if (appAbs(pMeasure->m_lstAverageResults[1][1] - fExpected2) > F(2.0))
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
    if (appAbs(pMeasure->m_lstAverageResults[0][0] - fExpected1) > F(0.3))
    {
        ++uiError;
    }
    if (appAbs(pMeasure->m_lstAverageResults[1][1] - fExpected2) > F(0.3))
    {
        ++uiError;
    }

#endif
    return uiError;
}

__REGIST_TEST(TestFermionUpdatorWithMesonCorrelatorStaggered, Updator, TestFermionUpdatorWithMesonCorrelatorStaggered);

//=============================================================================
// END OF FILE
//=============================================================================
