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
    UINT uiError = 0;

    Real fJG1 = F(0.0);
    Real fJG2 = F(0.0);
    Real fJGS1 = F(0.0);
    Real fJGS2 = F(0.0);
    Real fJGChen1 = F(0.0);
    Real fJGChen2 = F(0.0);
    Real fJGChenApprox1 = F(0.0);
    Real fJGChenApprox2 = F(0.0);
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
    fJGChen1 = pJG->m_lstJGChenAll[0];
    fJGChenApprox1 = pJG->m_lstJGChenApproxAll[0];
    fJFL1 = pJF->m_lstJLAll[0];
    fJFS1 = pJF->m_lstJSAll[0];
    fJFLPure1 = pJF->m_lstJLPureAll[0];
    fJFPot1 = pJF->m_lstJPotAll[0];

    //====== Step 4 ========
    //Gauge fixing to random Set A pure
    pRandom->GaugeFixing(appGetLattice()->m_pGaugeField);
    pRandom->AlsoFixingFermion(pF1);
    pRandom->AlsoFixingAphys(appGetLattice()->m_pAphys);
    appGeneral(_T("After gauge fixing to random!\n"));

    appGetLattice()->m_pMeasurements->Reset();
    pF1->CopyTo(pF2);
    pF2->InverseD(appGetLattice()->m_pGaugeField);
    pJG->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
    pJF->OnConfigurationAcceptedZ4(appGetLattice()->m_pGaugeField, NULL, pF1, pF2, TRUE, TRUE);
    appGetLattice()->m_pMeasurements->Report();

    fJG2 = pJG->m_lstJGAll[0];
    fJGS2 = pJG->m_lstJGSAll[0];
    fJGChen2 = pJG->m_lstJGChenAll[0];
    fJGChenApprox2 = pJG->m_lstJGChenApproxAll[0];
    fJFL2 = pJF->m_lstJLAll[0];
    fJFS2 = pJF->m_lstJSAll[0];
    fJFLPure2 = pJF->m_lstJLPureAll[0];
    fJFPot2 = pJF->m_lstJPotAll[0];

    //====== Step compare ========
    //Gauge fixing to random Set A pure
    Real fJGChen3 = F(0.0);
    Real fJGChen4 = F(0.0);
    Real fJGChen5 = F(0.0);
    Real fJGChenApprox3 = F(0.0);
    Real fJGChenApprox4 = F(0.0);
    Real fJGChenApprox5 = F(0.0);

    pRandom->GaugeFixing(appGetLattice()->m_pGaugeField);
    pRandom->AlsoFixingAphys(appGetLattice()->m_pAphys);
    appGeneral(_T("After gauge fixing to random!\n"));
    appGetLattice()->m_pMeasurements->Reset();
    pJG->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
    fJGChen3 = pJG->m_lstJGChenAll[0];
    fJGChenApprox3 = pJG->m_lstJGChenApproxAll[0];

    pRandom->GaugeFixing(appGetLattice()->m_pGaugeField);
    pRandom->AlsoFixingAphys(appGetLattice()->m_pAphys);
    appGeneral(_T("After gauge fixing to random!\n"));
    appGetLattice()->m_pMeasurements->Reset();
    pJG->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
    fJGChen4 = pJG->m_lstJGChenAll[0];
    fJGChenApprox4 = pJG->m_lstJGChenApproxAll[0];

    pRandom->GaugeFixing(appGetLattice()->m_pGaugeField);
    pRandom->AlsoFixingAphys(appGetLattice()->m_pAphys);
    appGeneral(_T("After gauge fixing to random!\n"));
    appGetLattice()->m_pMeasurements->Reset();
    pJG->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
    fJGChen5 = pJG->m_lstJGChenAll[0];
    fJGChenApprox5 = pJG->m_lstJGChenApproxAll[0];

    pF1->Return();
    pF2->Return();

    if (appAbs(fJG1 - fJG2) > F(0.000001))
    {
        ++uiError;
    }

    if (appAbs(fJGS1 - fJGS2) > F(0.000001))
    {
        ++uiError;
    }

    if (appAbs(fJFL1 - fJFL2) > F(0.000001))
    {
        ++uiError;
    }

    if (appAbs(fJFS1 - fJFS2) > F(0.000001))
    {
        ++uiError;
    }

    if (appAbs(fJFLPure1 - fJFLPure2) > F(0.000001))
    {
        ++uiError;
    }

    if (appAbs(fJFPot1 - fJFPot2) > F(0.000001))
    {
        ++uiError;
    }

    if (appAbs(fJFLPure1 - fJFL1 - fJFPot1) > F(0.000001))
    {
        ++uiError;
    }

    appGeneral(_T("Ji Decompose JG: Coulomb = %2.12f  vs  Random = %2.12f\n    delta = %2.12f\n"), fJG1, fJG2, appAbs(fJG1 - fJG2));
    appGeneral(_T("Wa Decompose JGS: Coulomb = %2.12f  vs  Random = %2.12f\n    delta = %2.12f\n"), fJGS1, fJGS2, appAbs(fJGS1 - fJGS2));
    appGeneral(_T("Chen Decompose JGChen: Coulomb = %2.12f  vs  Random = %2.12f\n    delta = %2.12f\n"), fJGChen1, fJGChen2, appAbs(fJGChen1 - fJGChen2));
    appGeneral(_T("Chen Decompose JGChen Approx: Coulomb = %2.12f  vs  Random = %2.12f\n    delta = %2.12f\n"), fJGChenApprox1, fJGChenApprox2, appAbs(fJGChenApprox1 - fJGChenApprox2));

    appGeneral(_T("Ji Decompose JFL: Coulomb = %2.12f  vs  Random = %2.12f\n    delta = %2.12f\n"), fJFL1, fJFL2, appAbs(fJFL1 - fJFL2));
    appGeneral(_T("Ji Decompose JFS: Coulomb = %2.12f  vs  Random = %2.12f\n    delta = %2.12f\n"), fJFS1, fJFS2, appAbs(fJFS1 - fJFS2));
    appGeneral(_T("Chen Decompose JFPure: Coulomb = %2.12f  vs  Random = %2.12f\n    delta = %2.12f\n"), fJFLPure1, fJFLPure2, appAbs(fJFLPure1 - fJFLPure2));
    appGeneral(_T("Chen Decompose JFPot: Coulomb = %2.12f  vs  Random = %2.12f\n    delta = %2.12f\n"), fJFPot1, fJFPot2, appAbs(fJFPot1 - fJFPot2));

    appGeneral(_T("JFPure - JFL - JFPot = %2.12f\n"), fJFLPure1 - fJFL1 - fJFPot1);
    appGeneral(_T("JG, JGS+JGChen-JFPot, %2.12f, %2.12f\n"), fJG1, fJGS1 + fJGChen1 - fJFPot1);
    appGeneral(_T("JGChen Approx 1,2,3,4,5 = %f, %f, %f, %f, %f\n"), fJGChenApprox1, fJGChenApprox2, fJGChenApprox3, fJGChenApprox4, fJGChenApprox5);

    delete pRandom;


    return uiError;
}

__REGIST_TEST(TestAngularMomentum, Misc, TestAngularMomentum);


//=============================================================================
// END OF FILE
//=============================================================================
