//=============================================================================
// FILENAME : TestBoundary.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [04/20/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestAngularMomentumDebug(CParameters& sParam)
{
    CGaugeFixingRandom* pRandom = new CGaugeFixingRandom();
    appGetLattice()->m_pGaugeField->FixBoundary();
    pRandom->Initial(appGetLattice(), sParam);
    CMeasureAMomentumJG* pJG = dynamic_cast<CMeasureAMomentumJG*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureAMomentumStochastic* pJF = dynamic_cast<CMeasureAMomentumStochastic*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CFieldFermionWilsonSquareSU3DR* pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(appGetLattice()->GetPooledFieldById(2));
    CFieldFermionWilsonSquareSU3DR* pF2 = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(appGetLattice()->GetPooledFieldById(2));
    pF1->InitialField(EFIT_RandomGaussian);
    pF1->FixBoundary();

    Real fJG1 = F(0.0);
    Real fJG2 = F(0.0);
    Real fJGS1 = F(0.0);
    Real fJGS2 = F(0.0);
    Real fJGChen1 = F(0.0);
    Real fJGChen2 = F(0.0);
    Real fJFL1 = F(0.0);
    Real fJFL2 = F(0.0);
    Real fJFS1 = F(0.0);
    Real fJFS2 = F(0.0);
    Real fJFLPure1 = F(0.0);
    Real fJFLPure2 = F(0.0);
    Real fJFPot1 = F(0.0);
    Real fJFPot2 = F(0.0);

    //====== Step 1 ========
    //Gauge fixing to Coulomb
    appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField);
    
    //====== Step 2 ========
    //Set A phys and A pure
    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField);
    appGetLattice()->SetAPure(appGetLattice()->m_pGaugeField);
    
    //====== Step 3 ========
    //Measure
    appGetLattice()->m_pMeasurements->Reset();
    pF1->CopyTo(pF2);
    pF2->InverseD(appGetLattice()->m_pGaugeField);
    pJG->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
    pJF->OnConfigurationAcceptedZ4(appGetLattice()->m_pGaugeField, NULL, pF1, pF2, TRUE, TRUE);
    appGetLattice()->m_pMeasurements->Report();
    fJG1 = pJG->m_lstJGAll[0];
    fJGS1 = pJG->m_lstJGSAll[0];
    fJGChen1 = pJG->m_lstJGChen[0];
    fJFL1 = pJF->m_lstJLAll[0];
    fJFS1 = pJF->m_lstJSAll[0];
    fJFLPure1 = pJF->m_lstJLPureAll[0];
    fJFPot1 = pJF->m_lstJPotAll[0];

    //====== Step 4 ========
    //Gauge fixing to random Set A pure
    pRandom->GaugeFixing(appGetLattice()->m_pGaugeField);
    pRandom->AlsoFixingFermion(pF1);
    pRandom->AlsoFixingAphys(appGetLattice()->m_pAphys);
    appGetLattice()->SetAPure(appGetLattice()->m_pGaugeField);
    appGeneral(_T("After gauge fixing to random!\n"));

    appGetLattice()->m_pMeasurements->Reset();
    pF1->CopyTo(pF2);
    pF2->InverseD(appGetLattice()->m_pGaugeField);
    pJG->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
    pJF->OnConfigurationAcceptedZ4(appGetLattice()->m_pGaugeField, NULL, pF1, pF2, TRUE, TRUE);
    appGetLattice()->m_pMeasurements->Report();

    fJG2 = pJG->m_lstJGAll[0];
    fJGS2 = pJG->m_lstJGSAll[0];
    fJGChen2 = pJG->m_lstJGChen[0];
    fJFL2 = pJF->m_lstJLAll[0];
    fJFS2 = pJF->m_lstJSAll[0];
    fJFLPure2 = pJF->m_lstJLPureAll[0];
    fJFPot2 = pJF->m_lstJPotAll[0];

    //====== Step compare ========
    //Gauge fixing to random Set A pure

    pF1->Return();
    pF2->Return();

    appGeneral(_T("Ji Decompose JG: Coulomb = %f  vs  Random = %f\n    delta = %2.12f\n"), fJG1, fJG2, appAbs(fJG1 - fJG2));
    appGeneral(_T("Wa Decompose JGS: Coulomb = %f  vs  Random = %f\n    delta = %2.12f\n"), fJGS1, fJGS2, appAbs(fJGS1 - fJGS2));
    appGeneral(_T("Chen Decompose JGChen: Coulomb = %f  vs  Random = %f\n    delta = %2.12f\n"), fJGChen1, fJGChen2, appAbs(fJGChen1 - fJGChen2));

    appGeneral(_T("Ji Decompose JFL: Coulomb = %f  vs  Random = %f\n    delta = %2.12f\n"), fJFL1, fJFL2, appAbs(fJFL1 - fJFL2));
    appGeneral(_T("Ji Decompose JFS: Coulomb = %f  vs  Random = %f\n    delta = %2.12f\n"), fJFS1, fJFS2, appAbs(fJFS1 - fJFS2));
    appGeneral(_T("Chen Decompose JFPure: Coulomb = %f  vs  Random = %f\n    delta = %2.12f\n"), fJFLPure1, fJFLPure2, appAbs(fJFLPure1 - fJFLPure2));
    appGeneral(_T("Chen Decompose JFPot: Coulomb = %2.12f  vs  Random = %f\n    delta = %2.12f\n"), fJFPot1, fJFPot2, appAbs(fJFPot1 - fJFPot2));

    delete pRandom;


    return 0;
}

#if _CLG_DEBUG

__REGIST_TEST(TestAngularMomentumDebug, Misc, TestAngularMomentum);

#else

#endif


//=============================================================================
// END OF FILE
//=============================================================================
