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
#if _CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedResDebug"), fExpected);
#else
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
#endif
    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    Real fHdiff = F(0.5);
    INT iAccept = 3;
#if _CLG_DEBUG
    sParam.FetchValueReal(_T("ExpectedHdiffDebug"), fHdiff);
    sParam.FetchValueINT(_T("ExpectedAcceptDebug"), iAccept);
#else
    sParam.FetchValueReal(_T("ExpectedHdiff"), fHdiff);
    sParam.FetchValueINT(_T("ExpectedAccept"), iAccept);
#endif

    appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
    appGetLattice()->m_pUpdator->Update(3, TRUE);
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);

    appGetLattice()->m_pMeasurements->Reset();
#if !_CLG_DEBUG
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->Update(40, TRUE);
#else
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->Update(5, TRUE);
    pMeasure->Average();
    Real fRes = pMeasure->GetAverageRealRes();

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real Hdiff = appGetLattice()->m_pUpdator->GetHDiff();
    appGeneral(_T("accepted : expected >= %d res=%d\n"), iAccept, uiAccept);
    appGeneral(_T("HDiff average : expected < %f res=%f\n"), fHdiff, Hdiff);

    appGeneral(_T("res : expected=%f res=%f\n"), fExpected, fRes);
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
    appGeneral(_T("accept (%d/43) : expected >= 35. HDiff = %f : expected < 0.3\n (exp(-0.3) is 74%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());

    if (uiAccept < 35)
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

__REGIST_TEST(TestUpdateCommon, UpdatorKS, TestFermionUpdatorKS, KS);

__REGIST_TEST(TestUpdateCommon, UpdatorKS, TestFermionUpdatorKSNestedForceGradient, nestedForceGrad);

//multi-level integrator not work efficiently now
//__REGIST_TEST(TestUpdateCommon, UpdatorKS, TestFermionUpdatorKSNestedForceGradientNf2p1, multilevelForceGrad);

//__REGIST_TEST(TestUpdateCommon, UpdatorKS, TestFermionUpdatorKSNestedOmelyanNf2p1, mutlevelOmelyan);
//__REGIST_TEST(TestUpdateCommon, UpdatorKS, TestFermionUpdatorKSNestedOmelyanNf2p1MultiField, NestedOmelyanNf2p1MultiField);
//__REGIST_TEST(TestUpdateCommon, UpdatorKS, TestFermionUpdatorKSNestedForceGradientNf2p1MultiField, NestedForceGradient);
//__REGIST_TEST(TestUpdateCommon, UpdatorKS, TestFermionUpdatorKSNested11StageNf2p1MultiField, Nested11StageNf2p1MultiField);

__REGIST_TEST(TestUpdateCommon, UpdatorKS, TestFermionUpdatorKSP4, P4);

___REGIST_TEST(TestUpdateCommon, UpdatorKS, TestFermionUpdatorKSGamma, KSGamma, _TEST_RELEASE);
___REGIST_TEST(TestUpdateCommon, UpdatorKS, TestFermionUpdatorKSGammaProj, KSGammaProj, _TEST_RELEASE);
___REGIST_TEST(TestUpdateCommon, UpdatorKS, TestFermionUpdatorKSGammaEM, KSGammaEM, _TEST_RELEASE);
___REGIST_TEST(TestUpdateCommon, UpdatorKS, TestFermionUpdatorKSGammaEMProj, KSGammaEMProj, _TEST_RELEASE);



