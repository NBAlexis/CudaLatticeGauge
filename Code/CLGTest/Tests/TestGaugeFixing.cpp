//=============================================================================
// FILENAME : TestZ3Symmetry.cpp
// 
// DESCRIPTION:
//
//     Test the Z3 Symmetry
//
// REVISION:
//  [06/23/2019 nbale]
//=============================================================================

#include "CLGTest.h"

#define _tfftMX 9
#define _tfftMY 10
#define _tfftMZ 11
#define _tfftMT 12

UINT TestFFT(CParameters&)
{
    CCLGFFTHelper::TestFFT();
    return 0;
}

UINT TestGaugeFixingLandau(CParameters&)
{
    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1)->GetCopy());
    CActionGaugePlaquette* pAction1 = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
    const Real fBeforeEnergy1 = pAction1->Energy(FALSE, pGauge, NULL);

    appGetLattice()->m_pGaugeFixing->GaugeFixing(pGauge);
    const Real fDivation = appGetLattice()->m_pGaugeFixing->CheckRes(pGauge);
    const Real fAfterEnergy1 = pAction1->Energy(FALSE, pGauge, NULL);

    UINT uiError = 0;
    if (fDivation > F(0.000001))
    {
        ++uiError;
    }

    if (abs(fBeforeEnergy1 - fAfterEnergy1) > 0.000001)
    {
        ++uiError;
    }
    appGeneral(_T("Gauge fixing with divation = %f, Before Energy = %f, After Energy = %f\n"), 
        fDivation, fBeforeEnergy1, fAfterEnergy1);

    appSafeDelete(pGauge);

    return uiError;
}

UINT TestGaugeFixingCoulombDR(CParameters&)
{
    UINT uiError = 0;
    CFieldGaugeSU3D* pGauge = dynamic_cast<CFieldGaugeSU3D*>(appGetLattice()->GetFieldById(1)->GetCopy());
    CFieldFermionWilsonSquareSU3DR* pFermion = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(appGetLattice()->GetFieldById(2));
    pFermion->PrepareForHMCOnlyRandomize();

    CFieldFermionWilsonSquareSU3DR* pFermion2 = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(pFermion->GetCopy());
    pFermion->PrepareForHMCNotRandomize(pGauge);
    CActionGaugePlaquetteRotating* pAction1 = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->GetActionById(1));
    CActionFermionWilsonNf2* pAction2 = dynamic_cast<CActionFermionWilsonNf2*>(appGetLattice()->GetActionById(2));
    const Real fEnergy1 = pAction1->Energy(FALSE, pGauge, NULL);
    pAction2->m_pFerimionField = pFermion;
    const Real fEnergy2 = pAction2->Energy(FALSE, pGauge, NULL);
    
    appGetLattice()->m_pGaugeFixing->GaugeFixing(pGauge);
    const Real fError = appGetLattice()->m_pGaugeFixing->CheckRes(pGauge);

    pFermion2->PrepareForHMCNotRandomize(pGauge);
    const Real fEnergy3 = pAction1->Energy(FALSE, pGauge, NULL);
    pAction2->m_pFerimionField = pFermion2;
    const Real fEnergy4 = pAction2->Energy(FALSE, pGauge, NULL);

    if (fError > F(0.000001))
    {
        ++uiError;
    }
    if (abs(fEnergy1 - fEnergy3) > 0.000001)
    {
        ++uiError;
    }
    if (abs(fEnergy2 - fEnergy4) > 0.000001)
    {
        ++uiError;
    }
    appGeneral(_T("Divation = %2.12f\n"), fError);
    appGeneral(_T("Before Energy1 = %f, Energy2 = %f\n"), fEnergy1, fEnergy2);
    appGeneral(_T("After Energy1 = %f, Energy2 = %f\n"), fEnergy3, fEnergy4);

    appSafeDelete(pGauge);
    appSafeDelete(pFermion2);

    return uiError;
}

//test whether the chiral condensation is respected by gauge fixing
//only the random gauge fixing can also gauge transform the fermion field
UINT TestGaugeFixingCoulombDRChiral(CParameters& sParam)
{
    UINT uiError = 0;
    CFieldGaugeSU3D* pGauge = dynamic_cast<CFieldGaugeSU3D*>(appGetLattice()->GetFieldById(1)->GetCopy());
    CGaugeFixingRandom* pRandom = new CGaugeFixingRandom();
    appGetLattice()->m_pGaugeField->FixBoundary();
    pRandom->Initial(appGetLattice(), sParam);

    //Calculate condensation
    CFieldFermionWilsonSquareSU3DR* pFermion = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(appGetLattice()->GetFieldById(2));
    pFermion->InitialField(EFIT_RandomGaussian);
    pFermion->FixBoundary();
    CFieldFermionWilsonSquareSU3DR* pFermion2 = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(pFermion->GetCopy());
    pFermion2->InverseD(pGauge);
    pFermion2->FixBoundary();
    CMeasureChiralCondensate* pCC = dynamic_cast<CMeasureChiralCondensate*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    pCC->Reset();
    pCC->OnConfigurationAcceptedZ4(pGauge, NULL, pFermion, pFermion2, TRUE, TRUE);

    //pGauge->DebugPrintMe();

    //Extract results
    CLGComplex oldAll[CMeasureChiralCondensate::_kCondMeasureCount];
    CLGComplex oldPosition0[CMeasureChiralCondensate::_kCondMeasureCount];
    for (UINT i = 0; i < CMeasureChiralCondensate::_kCondMeasureCount; ++i)
    {
        oldAll[i] = pCC->m_lstCondAll[i][0];
        oldPosition0[i] = pCC->m_lstCond[i][0];
    }

    for (INT i = 0; i < 5; ++i)
    {
        //===============================
        // 
        //===============================
        pRandom->GaugeFixing(pGauge); //the transform is randomized every GaugeFixing call
        pRandom->AlsoFixingFermion(pFermion);

        pFermion->CopyTo(pFermion2);
        pFermion2->InverseD(pGauge);
        pFermion2->FixBoundary();
        pCC->Reset();
        pCC->OnConfigurationAcceptedZ4(pGauge, NULL, pFermion, pFermion2, TRUE, TRUE);

        for (UINT i = 0; i < CMeasureChiralCondensate::_kCondMeasureCount; ++i)
        {
            //reset, so, the index is 0
            const CLGComplex toBeCompareAll = pCC->m_lstCondAll[i][0];
            const CLGComplex toBeComparePosition0 = pCC->m_lstCond[i][0];

            appGeneral(_T("Cond[%d]: before = %2.12f %2.12f I  after = %2.12f %2.12f I\n"), i, 
                oldAll[i].x, oldAll[i].y, toBeCompareAll.x, toBeCompareAll.y);
            if (__cuCabsSqf(_cuCsubf(oldAll[i], toBeCompareAll)) > F(0.00000001))
            {
                ++uiError;
            }

            appGeneral(_T("Cond[%d] at 0: before = %2.12f %2.12f I  after = %2.12f %2.12f I\n"), i,
                oldPosition0[i].x, oldPosition0[i].y, toBeComparePosition0.x, toBeComparePosition0.y);
            if (__cuCabsSqf(_cuCsubf(oldPosition0[i], toBeComparePosition0)) > F(0.00000001))
            {
                ++uiError;
            }
        }
    }

    appSafeDelete(pFermion2);

    return uiError;
}

__REGIST_TEST(TestFFT, Misc, TestFFT);

__REGIST_TEST(TestGaugeFixingLandau, Misc, TestGaugeFixingLandauCornell);

__REGIST_TEST(TestGaugeFixingLandau, Misc, TestGaugeFixingCoulombCornell);

__REGIST_TEST(TestGaugeFixingLandau, Misc, TestGaugeFixingLandauLosAlamos);

__REGIST_TEST(TestGaugeFixingLandau, Misc, TestGaugeFixingCoulombLosAlamos);

__REGIST_TEST(TestGaugeFixingCoulombDR, Misc, TestGaugeFixingCoulombCornellDR);

__REGIST_TEST(TestGaugeFixingCoulombDR, Misc, TestGaugeFixingCoulombLosAlamosDR);

__REGIST_TEST(TestGaugeFixingCoulombDRChiral, Misc, TestGaugeFixingCoulombDRChiral);


//=============================================================================
// END OF FILE
//=============================================================================
