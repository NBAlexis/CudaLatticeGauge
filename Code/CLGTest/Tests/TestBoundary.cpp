//=============================================================================
// FILENAME : TestBoundary.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [04/20/2019 nbale]
//=============================================================================

#include "CLGTest.h"

//UINT TestBoundary(CParameters& sParam)
//{
//    Real fExpected = F(0.565);
//    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
//
//    //we calculate staple energy from beta = 1 - 6
//    CActionGaugePlaquette * pAction = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
//    if (NULL == pAction)
//    {
//        return 1;
//    }
//    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
//    if (NULL == pMeasure)
//    {
//        return 1;
//    }
//
//    //pAction->SetBeta(F(3.0));
//
//    //Equilibration
//    appGetLattice()->m_pUpdator->Update(5, FALSE);
//
//    //Measure
//    pMeasure->Reset();
//    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
//    appGetLattice()->m_pUpdator->Update(20, TRUE);
//
//    pMeasure->Average();
//    Real fRes = pMeasure->GetAverageRealRes();
//    appGeneral(_T("res : expected=%f res=%f "), fExpected, fRes);
//    UINT uiError = 0;
//    if (appAbs(fRes - fExpected) > F(0.02))
//    {
//        ++uiError;
//    }
//
//    const UINT uiAccept = appGetLattice()->m_pUpdator->GetConfigurationCount();
//    const Real fHDiff = static_cast<Real>(appGetLattice()->m_pUpdator->GetHDiff());
//    appGeneral(_T("accept (%d/25) : expected >= 23. HDiff = %f : expected < 0.1 (exp(-0.1)=90%%)\n"), uiAccept, appGetLattice()->m_pUpdator->GetHDiff());
//
//    if (uiAccept < 23)
//    {
//        ++uiError;
//    }
//
//    if (fHDiff > F(0.1))
//    {
//        ++uiError;
//    }
//
//    return uiError;
//}

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

    //TArray<CFieldGauge*> gaugefields;
    //gaugefields.AddItem(appGetLattice()->m_pGaugeField);

    const Real fEnergy1 = static_cast<Real>(appGetLattice()->m_pActionList[0]->Energy(
        FALSE, _FIELDS, NULL));
    const Real fEnergy2 = static_cast<Real>(appGetLattice()->m_pActionList[1]->Energy(
        FALSE, _FIELDS, NULL));
    const Real fEnergy3 = static_cast<Real>(appGetLattice()->m_pActionList[2]->Energy(
        FALSE, _FIELDS, NULL));

    CFieldGauge* pStape = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField[0]->GetCopy());
    appGetLattice()->m_pGaugeField[0]->CalculateOnlyStaple(pStape);
    const Real fEnergy4 = static_cast<Real>(appGetLattice()->m_pGaugeField[0]->CalculatePlaqutteEnergyUsingStable(
        appGetLattice()->m_pActionList[0]->GetBetaOverN(), pStape));

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

UINT TestEtaShift(CParameters& sParam)
{
    CFieldFermionKSSU3* pF1 = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetFieldById(2));
    CFieldFermionKSSU3* pF2 = dynamic_cast<CFieldFermionKSSU3*>(pF1->GetCopy());
    CFieldFermionKSU1* pF3 = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetFieldById(3));
    CFieldFermionKSU1* pF4 = dynamic_cast<CFieldFermionKSU1*>(pF3->GetCopy());
    CFieldGaugeSU3* pSU3 = new CFieldGaugeSU3();
    pSU3->InitialField(EFIT_Random);
    pSU3->m_byFieldId = 1;
    TArray<CFieldGauge*> su3;
    su3.AddItem(pSU3);
    CFieldGaugeU1* pU1 = new CFieldGaugeU1();
    pU1->InitialField(EFIT_Random);
    pU1->m_byFieldId = 1;
    TArray<CFieldGauge*> u1;
    u1.AddItem(pU1);
    pF1->TestSetEtaShift(TRUE);
    pF2->TestSetEtaShift(FALSE);
    pF3->TestSetEtaShift(TRUE);
    pF4->TestSetEtaShift(FALSE);
    UINT uiError = 0;
    pF1->D(1, 0, su3.GetData(), NULL);
    pF2->D(1, 0, su3.GetData(), NULL);
    pF1->AxpyMinus(pF2);
    const DOUBLE fDiff1 = cuCabs(pF1->Dot(pF1));
    appGeneral(_T("Diff1 = %2.12f\n"), fDiff1);
    if (fDiff1 > F(0.000001))
    {
        ++uiError;
    }
    pF3->D(1, 0, u1.GetData(), NULL);
    pF4->D(1, 0, u1.GetData(), NULL);
    pF3->AxpyMinus(pF4);
    const DOUBLE fDiff2 = cuCabs(pF3->Dot(pF3));
    appGeneral(_T("Diff2 = %2.12f\n"), fDiff2);
    if (fDiff2 > F(0.000001))
    {
        ++uiError;
    }
    appSafeDelete(pF2);
    appSafeDelete(pF4);
    appSafeDelete(pSU3);
    appSafeDelete(pU1);

    return uiError;
}

//As long as rotation - Wilson Dirac Dirichlet, Projective plane works, the below should work
//__REGIST_TEST(TestBoundary, Updator, TestDirichletBoundary, Dirichlet);
//__REGIST_TEST(TestBoundary, Updator, TestProjectivePlaneBoundary, ProjectivePlane);

//It is now not supported
//__REGIST_TEST(TestBoundaryMapping, Misc, TestBoundaryMapping, BoundaryMapping);

__REGIST_TEST(TestEtaShift, Misc, TestEtaShift, EtaShift);


//=============================================================================
// END OF FILE
//=============================================================================
