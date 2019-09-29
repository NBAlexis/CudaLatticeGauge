//=============================================================================
// FILENAME : TestBoundary.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [04/20/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestAngularMomentum(CParameters& sParam)
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

    //test gauge invarience of E
    //CFieldGaugeSU3D* pGaugeCopy = dynamic_cast<CFieldGaugeSU3D*>(appGetLattice()->m_pGaugeField->GetCopy());
    //appGetLattice()->m_pGaugeField->CalculateE_Using_U(pGaugeCopy);
    //appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField);
    
    //pRandom->GaugeFixing(appGetLattice()->m_pGaugeField);
    //appGetLattice()->m_pGaugeField->CalculateE_Using_U(pGaugeCopy);

    //pF1->FixBoundary();
    //pF1->CopyTo(pF2);
    //pF2->FixBoundary();
    //pF2->D(appGetLattice()->m_pGaugeField);
    //pRandom->GaugeFixing(appGetLattice()->m_pGaugeField);
    //pRandom->AlsoFixingFermion(pF1);
    //pRandom->AlsoFixingFermion(pF2);
    //pF1->D(appGetLattice()->m_pGaugeField);
    //pF1->AxpyMinus(pF2);
    //appGeneral(_T("Dot Res = %2.12f\n"), pF1->Dot(pF1).x);

#if 1
    Real fJG1 = F(0.0);
    Real fJG2 = F(0.0);
    Real fJG3 = F(0.0);
    CLGComplex fJGS1 = _zeroc;
    CLGComplex fJGS2 = _zeroc;
    CLGComplex fJGS3 = _zeroc;
    CLGComplex fJGChen1 = _zeroc;
    CLGComplex fJGChen2 = _zeroc;
    CLGComplex fJGChen3 = _zeroc;
    Real fJFL1 = F(0.0);
    Real fJFL2 = F(0.0);
    Real fJFL3 = F(0.0);
    Real fJFS1 = F(0.0);
    Real fJFS2 = F(0.0);
    Real fJFS3 = F(0.0);
    Real fJFLPure1 = F(0.0);
    Real fJFLPure2 = F(0.0);
    Real fJFLPure3 = F(0.0);
    Real fJFPot1 = F(0.0);
    Real fJFPot2 = F(0.0);
    Real fJFPot3 = F(0.0);

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

    //====== Step 5 ========
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

    fJG3 = pJG->m_lstJGAll[0];
    fJGS3 = pJG->m_lstJGSAll[0];
    fJGChen3 = pJG->m_lstJGChen[0];
    fJFL3 = pJF->m_lstJLAll[0];
    fJFS3 = pJF->m_lstJSAll[0];
    fJFLPure3 = pJF->m_lstJLPureAll[0];
    fJFPot3 = pJF->m_lstJPotAll[0];

    pF1->Return();
    pF2->Return();

    appGeneral(_T("Ji Decompose JG: Coulomb = %f  vs  Random1 = %f  vs  Random2 = %f\n"), fJG1, fJG2, fJG3);
    appGeneral(_T("Wa Decompose JGS: Coulomb = %f + %f i vs  Random1 = %f + %f i vs  Random2 = %f + %f i\n"), fJGS1.x, fJGS1.y, fJGS2.x, fJGS2.y, fJGS3.x, fJGS3.y);
    appGeneral(_T("Chen Decompose JGChen: Coulomb = %f + %f i  vs  Random1 = %f + %f i vs  Random2 = %f + %f i\n"), fJGChen1.x, fJGChen1.y, fJGChen2.x, fJGChen2.y, fJGChen3.x, fJGChen3.y);

    appGeneral(_T("Ji Decompose JFL: Coulomb = %2.12f  vs  Random1 = %2.12f  vs  Random2 = %2.12f\n"), fJFL1, fJFL2, fJFL3);
    appGeneral(_T("Ji Decompose JFS: Coulomb = %2.12f  vs  Random1 = %2.12f  vs  Random2 = %2.12f\n"), fJFS1, fJFS2, fJFS3);
    appGeneral(_T("Chen Decompose JFPure: Coulomb = %2.12f  vs  Random1 = %2.12f  vs  Random2 = %2.12f\n"), fJFLPure1, fJFLPure2, fJFLPure3);
    appGeneral(_T("Chen Decompose JFPot: Coulomb = %2.12f  vs  Random1 = %2.12f  vs  Random2 = %2.12f\n"), fJFPot1, fJFPot2, fJFPot3);

    delete pRandom;

#endif

    return 0;
}

__REGIST_TEST(TestAngularMomentum, Misc, TestAngularMomentum);


//=============================================================================
// END OF FILE
//=============================================================================
