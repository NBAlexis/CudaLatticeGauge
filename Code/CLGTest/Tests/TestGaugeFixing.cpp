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

#if !_CLG_DOUBLEFLOAT
#define _GAUGE_FIXING_EnergyERROR F(0.05)
#define _GAUGE_FIXING_ZeroERROR F(0.001)
#else
#define _GAUGE_FIXING_EnergyERROR 0.005
#define _GAUGE_FIXING_ZeroERROR 0.00005
#endif


UINT TestFFT(CParameters&)
{
    CCLGFFTHelper::TestFFT();
    return 0;
}

UINT TestGaugeFixingLandau(CParameters&)
{
    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1)->GetCopy());
    TArray<CFieldGauge*> gauge;
    gauge.AddItem(pGauge);
    CActionGaugePlaquette* pAction1 = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
    const Real fBeforeEnergy1 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, gauge.GetData(), NULL, NULL));

    appGetLattice()->m_pGaugeFixing->GaugeFixing(pGauge);
    const Real fDivation = static_cast<Real>(appGetLattice()->m_pGaugeFixing->CheckRes(pGauge));
    const Real fAfterEnergy1 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, gauge.GetData(), NULL, NULL));

    UINT uiError = 0;
    if (fDivation > _GAUGE_FIXING_ZeroERROR)
    {
        ++uiError;
    }

    if (appAbs(fBeforeEnergy1 - fAfterEnergy1) > _GAUGE_FIXING_EnergyERROR)
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
    TArray<CFieldGauge*> gaugefields;
    gaugefields.AddItem(pGauge);
    CFieldFermionWilsonSquareSU3DR* pFermion = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(appGetLattice()->GetFieldById(2));
    pFermion->PrepareForHMCOnlyRandomize();

    CFieldFermionWilsonSquareSU3DR* pFermion2 = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(pFermion->GetCopy());
    pFermion->PrepareForHMCNotRandomize(pGauge);
    CActionGaugePlaquetteRotating* pAction1 = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->GetActionById(1));
    CActionFermionWilsonNf2* pAction2 = dynamic_cast<CActionFermionWilsonNf2*>(appGetLattice()->GetActionById(2));
    const Real fEnergy1 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, gaugefields.GetData(), NULL, NULL));
    pAction2->SetFermionFieldTest(pFermion);
    const Real fEnergy2 = static_cast<Real>(pAction2->Energy(FALSE, 1, 0, gaugefields.GetData(), NULL, NULL));
    
    appGetLattice()->m_pGaugeFixing->GaugeFixing(pGauge);
    const Real fError = static_cast<Real>(appGetLattice()->m_pGaugeFixing->CheckRes(pGauge));

    pFermion2->PrepareForHMCNotRandomize(pGauge);
    const Real fEnergy3 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, gaugefields.GetData(), NULL, NULL));
    pAction2->SetFermionFieldTest(pFermion2);
    const Real fEnergy4 = static_cast<Real>(pAction2->Energy(FALSE, 1, 0, gaugefields.GetData(), NULL, NULL));

    if (fError > _GAUGE_FIXING_ZeroERROR)
    {
        ++uiError;
    }
    if (appAbs(fEnergy1 - fEnergy3) > _GAUGE_FIXING_EnergyERROR)
    {
        ++uiError;
    }
    if (appAbs(fEnergy2 - fEnergy4) > _GAUGE_FIXING_EnergyERROR)
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
    TArray<CFieldGauge*> gaugeFields;
    gaugeFields.AddItem(pGauge);
    CGaugeFixingRandom* pRandom = new CGaugeFixingRandom();
    appGetLattice()->m_pGaugeField[0]->FixBoundary();
    pRandom->Initial(appGetLattice(), sParam);

    //Calculate condensation
    CFieldFermionWilsonSquareSU3DR* pFermion = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(appGetLattice()->GetFieldById(2));
    pFermion->InitialField(EFIT_RandomGaussian);
    pFermion->FixBoundary();
    CFieldFermionWilsonSquareSU3DR* pFermion2 = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(pFermion->GetCopy());
    pFermion2->InverseD(1, 0, gaugeFields.GetData(), NULL);
    pFermion2->FixBoundary();
    CMeasureChiralCondensate* pCC = dynamic_cast<CMeasureChiralCondensate*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    pCC->Reset();
    pCC->OnConfigurationAcceptedZ4(1, 0, gaugeFields.GetData(), NULL, NULL, pFermion, pFermion2, TRUE, TRUE);

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
        pFermion2->InverseD(1, 0, gaugeFields.GetData(), NULL);
        pFermion2->FixBoundary();
        pCC->Reset();
        pCC->OnConfigurationAcceptedZ4(1, 0, gaugeFields.GetData(), NULL, NULL, pFermion, pFermion2, TRUE, TRUE);

        for (UINT i1 = 0; i1 < CMeasureChiralCondensate::_kCondMeasureCount; ++i1)
        {
            //reset, so, the index is 0
            const CLGComplex toBeCompareAll = pCC->m_lstCondAll[i1][0];
            const CLGComplex toBeComparePosition0 = pCC->m_lstCond[i1][0];

            appGeneral(_T("Cond[%d]: before = %2.12f %2.12f I  after = %2.12f %2.12f I\n"), i1, 
                oldAll[i1].x, oldAll[i1].y, toBeCompareAll.x, toBeCompareAll.y);
            if (__cuCabsSqf(_cuCsubf(oldAll[i1], toBeCompareAll)) > _GAUGE_FIXING_ZeroERROR)
            {
                ++uiError;
            }

            appGeneral(_T("Cond[%d] at 0: before = %2.12f %2.12f I  after = %2.12f %2.12f I\n"), i1,
                oldPosition0[i1].x, oldPosition0[i1].y, toBeComparePosition0.x, toBeComparePosition0.y);
            if (__cuCabsSqf(_cuCsubf(oldPosition0[i1], toBeComparePosition0)) > _GAUGE_FIXING_ZeroERROR)
            {
                ++uiError;
            }
        }
    }

    appSafeDelete(pFermion2);

    return uiError;
}

//test the action of rotation (both gauge and KS) are gauge invarient in the case of projective plane
//test the chiral and angular momentum measurement is gauge invarient in the case of projective plane
UINT TestGaugeFixingCoulombPorjectivePlane(CParameters&)
{
    UINT uiError = 0;
    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1)->GetCopy());
    TArray<CFieldGauge*> gaugeFields;
    gaugeFields.AddItem(pGauge);
    CFieldFermionKSSU3R* pFermion = dynamic_cast<CFieldFermionKSSU3R*>(appGetLattice()->GetFieldById(2));
    pFermion->PrepareForHMCOnlyRandomize();

    CFieldFermionKSSU3R* pFermion2 = dynamic_cast<CFieldFermionKSSU3R*>(pFermion->GetCopy());
    CFieldFermionKSSU3R* pF1W = dynamic_cast<CFieldFermionKSSU3R*>(pFermion->GetCopy());
    pF1W->InitialField(EFIT_RandomZ4);
    CFieldFermionKSSU3R* pF2W = dynamic_cast<CFieldFermionKSSU3R*>(pF1W->GetCopy());
    pF2W->InverseD(1, 0, gaugeFields.GetData(), NULL);
    pFermion->PrepareForHMCNotRandomize(1, 0, gaugeFields.GetData(), NULL);

    CActionGaugePlaquetteRotating* pAction1 = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->GetActionById(1));
    CActionFermionKS* pAction2 = dynamic_cast<CActionFermionKS*>(appGetLattice()->GetActionById(2));

    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureChiralCondensateKS* pCC = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CMeasureAngularMomentumKS* pAM = dynamic_cast<CMeasureAngularMomentumKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));

    const Real fEnergy1 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, gaugeFields.GetData(), NULL, NULL));
    pAction2->SetFermionFieldTest(pFermion);
    const Real fEnergy2 = static_cast<Real>(pAction2->Energy(FALSE, 1, 0, gaugeFields.GetData(), NULL, NULL));

    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField[0]);
    pPL->OnConfigurationAccepted(1, 0, gaugeFields.GetData(), NULL, NULL);
    const Real fPolyakov1 = _cuCabsf(pPL->m_lstLoop[0]);
    pCC->OnConfigurationAcceptedZ4(1, 0, gaugeFields.GetData(), NULL, NULL, pF1W, pF2W, TRUE, TRUE);
    const Real fChiralCond1 = _cuCabsf(pCC->m_lstCondAll[0][0]);
    const Real fConectSusp1 = _cuCabsf(pCC->m_lstCondAll[1][0]);
    pAM->OnConfigurationAcceptedZ4(1, 0, gaugeFields.GetData(), NULL, NULL, pF1W, pF2W, TRUE, TRUE);
    const Real fOrbital1 = _cuCabsf(pAM->m_lstCondAll[0][0]);
    const Real fSpin1 = _cuCabsf(pAM->m_lstCondAll[1][0]);
    const Real fPotential1 = _cuCabsf(pAM->m_lstCondAll[2][0]);

    appGeneral(_T("PL: %d, CC: %d, %d, AM: %d, %d, %d\n"),
        pPL->m_lstLoop.GetCount(),
        pCC->m_lstCondAll[0].GetCount(),
        pCC->m_lstCondAll[1].GetCount(),
        pAM->m_lstCondAll[0].GetCount(),
        pAM->m_lstCondAll[1].GetCount(),
        pAM->m_lstCondAll[2].GetCount()
        );

    //appGetLattice()->m_pGaugeFixing->GaugeFixing(pGauge);
    CGaugeFixingRandom* pRandom = new CGaugeFixingRandom();
    pRandom->Initial(appGetLattice(), CParameters());

    for (INT i = 0; i < 10; ++i)
    {
        pRandom->GaugeFixing(pGauge);
        pRandom->AlsoFixingFermion(pFermion2);
        pRandom->AlsoFixingFermion(pF1W);
    }

    pF1W->CopyTo(pF2W);
    pF2W->InverseD(1, 0, gaugeFields.GetData(), NULL);

    pFermion2->PrepareForHMCNotRandomize(1, 0, gaugeFields.GetData(), NULL);
    const Real fEnergy3 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, gaugeFields.GetData(), NULL, NULL));
    pAction2->SetFermionFieldTest(pFermion2);
    const Real fEnergy4 = static_cast<Real>(pAction2->Energy(FALSE, 1, 0, gaugeFields.GetData(), NULL, NULL));


    pPL->OnConfigurationAccepted(1, 0, gaugeFields.GetData(), NULL, NULL);
    const Real fPolyakov2 = _cuCabsf(pPL->m_lstLoop[0]);
    pCC->OnConfigurationAcceptedZ4(1, 0, gaugeFields.GetData(), NULL, NULL, pF1W, pF2W, TRUE, TRUE);
    const Real fChiralCond2 = _cuCabsf(pCC->m_lstCondAll[0][0]);
    const Real fConectSusp2 = _cuCabsf(pCC->m_lstCondAll[1][0]);
    pAM->OnConfigurationAcceptedZ4(1, 0, gaugeFields.GetData(), NULL, NULL, pF1W, pF2W, TRUE, TRUE);
    const Real fOrbital2 = _cuCabsf(pAM->m_lstCondAll[0][0]);
    const Real fSpin2 = _cuCabsf(pAM->m_lstCondAll[1][0]);
    const Real fPotential2 = _cuCabsf(pAM->m_lstCondAll[2][0]);

    appGeneral(_T("PL: %d, CC: %d, %d, AM: %d, %d, %d\n"),
        pPL->m_lstLoop.GetCount(),
        pCC->m_lstCondAll[0].GetCount(),
        pCC->m_lstCondAll[1].GetCount(),
        pAM->m_lstCondAll[0].GetCount(),
        pAM->m_lstCondAll[1].GetCount(),
        pAM->m_lstCondAll[2].GetCount()
    );

    if (appAbs(fEnergy1 - fEnergy3) > _GAUGE_FIXING_EnergyERROR)
    {
        ++uiError;
    }
    if (appAbs(fEnergy2 - fEnergy4) > _GAUGE_FIXING_EnergyERROR)
    {
        ++uiError;
    }

    if (appAbs(fPolyakov1 - fPolyakov2) > _GAUGE_FIXING_ZeroERROR)
    {
        ++uiError;
    }

    if (appAbs(fChiralCond1 - fChiralCond2) > _GAUGE_FIXING_ZeroERROR)
    {
        ++uiError;
    }
    if (appAbs(fConectSusp1 - fConectSusp2) > _GAUGE_FIXING_ZeroERROR)
    {
        ++uiError;
    }

    if (appAbs(fOrbital1 - fOrbital2) > _GAUGE_FIXING_ZeroERROR)
    {
        ++uiError;
    }
    if (appAbs(fSpin1 - fSpin2) > _GAUGE_FIXING_ZeroERROR)
    {
        ++uiError;
    }
    if (appAbs(fPotential1 - fPotential2) > _GAUGE_FIXING_ZeroERROR)
    {
        ++uiError;
    }

    //appGeneral(_T("Gauge Divation = %2.12f\n"), fError);
    appGeneral(_T("Gauge Energy before = %f, after = %f\n"), fEnergy1, fEnergy3);
    appGeneral(_T("Fermion Energy before = %f, after = %f\n"), fEnergy2, fEnergy4);

    appGeneral(_T("Polyakov loop before = %f, after = %f\n"), fPolyakov1, fPolyakov2);
    appGeneral(_T("Chiral Condensation before = %f, after = %f\n"), fChiralCond1, fChiralCond2);
    appGeneral(_T("Connect Susp before = %f, after = %f\n"), fConectSusp1, fConectSusp2);
    appGeneral(_T("Fermion Orbital before = %f, after = %f\n"), fOrbital1, fOrbital2);
    appGeneral(_T("Fermion Spin before = %f, after = %f\n"), fSpin1, fSpin2);
    appGeneral(_T("Fermion Potential before = %f, after = %f\n"), fPotential1, fPotential2);

    appSafeDelete(pGauge);
    appSafeDelete(pFermion2);
    appSafeDelete(pRandom);

    return uiError;
}

//test the coulomb gauge fixing works for projective plane
UINT TestGaugeFixingCoulombPorjectivePlane2(CParameters&)
{
    UINT uiError = 0;
    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1)->GetCopy());
    TArray<CFieldGauge*> gaugeFields;
    gaugeFields.AddItem(pGauge);
    CActionGaugePlaquetteRotating* pAction1 = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->GetActionById(1));

    const Real fEnergy1 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, gaugeFields.GetData(), NULL, NULL));

    appGetLattice()->m_pGaugeFixing->GaugeFixing(pGauge);
    const Real fError = static_cast<Real>(appGetLattice()->m_pGaugeFixing->CheckRes(pGauge));

    const Real fEnergy3 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, gaugeFields.GetData(), NULL, NULL));

    if (fError > _GAUGE_FIXING_ZeroERROR)
    {
        ++uiError;
    }
    if (appAbs(fEnergy1 - fEnergy3) > _GAUGE_FIXING_EnergyERROR)
    {
        ++uiError;
    }

    appGeneral(_T("Gauge Divation = %2.12f\n"), fError);
    appGeneral(_T("Gauge Energy before = %f, after = %f\n"), fEnergy1, fEnergy3);

    appSafeDelete(pGauge);

    return uiError;
}

___REGIST_TEST(TestFFT, Verify, TestFFT, FFT, _TEST_NOCHECK);

__REGIST_TEST(TestGaugeFixingLandau, GaugeFixing, TestGaugeFixingLandauCornell, LandauCornell);

//FFT not applied using single float
___REGIST_TEST(TestGaugeFixingLandau, GaugeFixing, TestGaugeFixingCoulombCornell, CoulombCornell, _TEST_DOUBLE);

__REGIST_TEST(TestGaugeFixingLandau, GaugeFixing, TestGaugeFixingLandauLosAlamos, LandauLosAlamos);

__REGIST_TEST(TestGaugeFixingLandau, GaugeFixing, TestGaugeFixingCoulombLosAlamos, CoulombLosAlamos);

//FFT not applied using single float
___REGIST_TEST(TestGaugeFixingCoulombDR, GaugeFixing, TestGaugeFixingCoulombCornellDR, CoulombCornellDR, _TEST_DOUBLE);

__REGIST_TEST(TestGaugeFixingCoulombDR, GaugeFixing, TestGaugeFixingCoulombLosAlamosDR, CoulombLosAlamosDR);

__REGIST_TEST(TestGaugeFixingCoulombDRChiral, GaugeFixing, TestGaugeFixingCoulombDRChiral, CoulombDRChiral);

__REGIST_TEST(TestGaugeFixingCoulombPorjectivePlane, GaugeFixing, TestGaugeFixingRotationKS, RotationKS);

__REGIST_TEST(TestGaugeFixingCoulombPorjectivePlane2, GaugeFixing, TestGaugeFixingRotationKS2, RotationKS2);


//=============================================================================
// END OF FILE
//=============================================================================
