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
    UINT uiError = 0;

    //Real fE = appGetLattice()->m_pGaugeField->CalculatePlaqutteEnergy(F(1.0));
    //CIndexData::DebugEdgeMapping(1, SSmallInt4(-1, -1, 0, 1));
    //CIndexData::DebugEdgeGlue(1, SSmallInt4(-1, -1, 0, 1));
    //CIndexData::DebugEdgeMapping(2, SSmallInt4(-1, -1, 0, 1));

#if 1

    const Real fExpE1 = F(14656.11807458415933069773);
    const Real fExpE2 = F(6518.69394147992716170847);
    const Real fExpE3 = F(1167.37873127234047387901);
    const Real fExpE4 = F(14498.55613259701931383461);

    const Real fEnergy1 = appGetLattice()->m_pActionList[0]->Energy(
        FALSE, appGetLattice()->m_pGaugeField, NULL);
    const Real fEnergy2 = appGetLattice()->m_pActionList[1]->Energy(
        FALSE, appGetLattice()->m_pGaugeField, NULL);
    const Real fEnergy3 = appGetLattice()->m_pActionList[2]->Energy(
        FALSE, appGetLattice()->m_pGaugeField, NULL);

    CFieldGauge* pStape = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField->GetCopy());
    appGetLattice()->m_pGaugeField->CalculateOnlyStaple(pStape);
    const Real fEnergy4 = appGetLattice()->m_pGaugeField->CalculatePlaqutteEnergyUsingStable(
        CCommonData::m_fBeta / F(3.0), pStape);

    if (appAbs(fEnergy1 - fExpE1) < F(0.000000001))
    {
        appGeneral(_T("action1 delta: %2.20f\n"), fEnergy1 - fExpE1);
    }
    else
    {
        ++uiError;
        appGeneral(_T("action1: %2.20f expecting: %2.20f\n"), fEnergy1, fExpE1);
    }
    
    if (appAbs(fEnergy2 - fExpE2) < F(0.000000001))
    {
        appGeneral(_T("action2 delta: %2.20f\n"), fEnergy2 - fExpE2);
    }
    else
    {
        ++uiError;
        appGeneral(_T("action2: %2.20f expecting: %2.20f\n"), fEnergy2, fExpE2);
    }

    if (appAbs(fEnergy3 - fExpE3) < F(0.000000001))
    {
        appGeneral(_T("action3 delta: %2.20f\n"), fEnergy3 - fExpE3);
    }
    else
    {
        ++uiError;
        appGeneral(_T("action3: %2.20f expecting: %2.20f\n"), fEnergy3, fExpE3);
    }

    if (appAbs(fEnergy4 - fExpE4) < F(0.000000001))
    {
        appGeneral(_T("staple delta: %2.20f\n"), fEnergy4 - fExpE4);
    }
    else
    {
        ++uiError;
        appGeneral(_T("staple: %2.20f expecting: %2.20f\n"), fEnergy4, fExpE4);
    }


    appSafeDelete(pStape);
#endif

    return uiError;
}

__REGIST_TEST(TestBoundary, Updator, TestDirichletBoundary);

__REGIST_TEST(TestBoundary, Updator, TestProjectivePlaneBoundary);

__REGIST_TEST(TestBoundaryMapping, Misc, TestBoundaryMapping);


//=============================================================================
// END OF FILE
//=============================================================================
