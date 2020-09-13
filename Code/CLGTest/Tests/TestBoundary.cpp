//=============================================================================
// FILENAME : TestBoundary.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [04/20/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestBoundary(CParameters& sParam)
{
    Real fExpected = F(0.565);
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
    appGetLattice()->m_pUpdator->Update(5, FALSE);

    //Measure
    pMeasure->Reset();
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->Update(20, TRUE);

    const Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f "), fExpected, fRes);
    UINT uiError = 0;
    if (appAbs(fRes - fExpected) > F(0.02))
    {
        ++uiError;
    }

    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
    const Real fHDiff = appGetLattice()->m_pUpdator->GetHDiff();
    appGeneral(_T("accept (%d/25) : expected >= 23. HDiff = %f : expected < 0.1 (exp(-0.1)=90%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());

    if (uiAccept < 23)
    {
        ++uiError;
    }

    if (fHDiff > F(0.1))
    {
        ++uiError;
    }

    return uiError;
}

UINT TestBoundaryMapping(CParameters& sParam)
{
    //CIndexData::DebugEdgeMapping(1);
    //CIndexData::DebugEdgeGlue(1);

    appGeneral(_T("action1: %2.20f\n"), appGetLattice()->m_pActionList[0]->Energy(
        FALSE, appGetLattice()->m_pGaugeField, NULL));
    appGeneral(_T("action2: %2.20f\n"), appGetLattice()->m_pActionList[1]->Energy(
        FALSE, appGetLattice()->m_pGaugeField, NULL));
    appGeneral(_T("action3: %2.20f\n"), appGetLattice()->m_pActionList[2]->Energy(
        FALSE, appGetLattice()->m_pGaugeField, NULL));

    CFieldGauge* pStape = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField->GetCopy());
    appGetLattice()->m_pGaugeField->CalculateOnlyStaple(pStape);
    appGeneral(_T("staple: %2.20f\n"), appGetLattice()->m_pGaugeField->CalculatePlaqutteEnergyUsingStable(CCommonData::m_fBeta / F(3.0), pStape));


    appSafeDelete(pStape);

    return 0;
}

__REGIST_TEST(TestBoundary, Updator, TestDirichletBoundary);
__REGIST_TEST(TestBoundaryMapping, Misc, TestBoundaryMapping);


//=============================================================================
// END OF FILE
//=============================================================================
