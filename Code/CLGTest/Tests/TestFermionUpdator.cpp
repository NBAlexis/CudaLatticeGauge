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
    Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f"), fExpected, fRes);
    UINT uiError = 0;
    if (appAbs(fRes - fExpected) > F(0.01))
    {
        ++uiError;
    }

    UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
    appGeneral(_T("accept (%d/60) : expected >= 50. HDiff = %f : expected < 0.08\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());

    if (uiAccept < 50)
    {
        ++uiError;
    }

    if (fHDiff > F(0.08))
    {
        ++uiError;
    }

    return uiError;
#endif
}

__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdator);

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
        if (appAbs(fRes - lstResExpected[i - 1]) > F(0.25))
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

    Real fRes = pMeasure->m_lstResults[0][0];
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

__REGIST_TEST(TestFermionUpdatorWithMesonCorrelator, Updator, TestGaugeSmearingAPEProj);

__REGIST_TEST(TestFermionUpdatorWithMesonCorrelator, Updator, TestGaugeSmearingAPEStout);

UINT TestFermionUpdatorL(CParameters& sParam)
{
    appGetLattice()->m_pUpdator->Update(1, TRUE);
    return 0;
}

#if !_CLG_DEBUG
__REGIST_TEST(TestFermionUpdatorL, Updator, TestFermionUpdatorLargeScale);
#endif

//=============================================================================
// END OF FILE
//=============================================================================
