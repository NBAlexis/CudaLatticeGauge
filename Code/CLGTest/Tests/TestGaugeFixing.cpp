//=============================================================================
// FILENAME : TestZ3Symmetry.cpp
// 
// DESCRIPTION:
//
//     Test the Z3 Symmetry
//
// REVISION:
//  [mm/dd/yy]
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
    const Real fBeforeEnergy1 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, 0, gauge.GetData(), NULL, NULL, NULL));

    appGetLattice()->m_pGaugeFixing->GaugeFixing(pGauge);
    const Real fDivation = static_cast<Real>(appGetLattice()->m_pGaugeFixing->CheckRes(pGauge));
    const Real fAfterEnergy1 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, 0, gauge.GetData(), NULL, NULL, NULL));

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
    CActionFermionKS* pAction2 = dynamic_cast<CActionFermionKS*>(appGetLattice()->GetActionById(2));
    const Real fEnergy1 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, 0, gaugefields.GetData(), NULL, NULL, NULL));
    pAction2->SetFermionFieldTest(pFermion);
    const Real fEnergy2 = static_cast<Real>(pAction2->Energy(FALSE, 1, 0, 0, gaugefields.GetData(), NULL, NULL, NULL));
    
    appGetLattice()->m_pGaugeFixing->GaugeFixing(pGauge);
    const Real fError = static_cast<Real>(appGetLattice()->m_pGaugeFixing->CheckRes(pGauge));

    pFermion2->PrepareForHMCNotRandomize(pGauge);
    const Real fEnergy3 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, 0, gaugefields.GetData(), NULL, NULL, NULL));
    pAction2->SetFermionFieldTest(pFermion2);
    const Real fEnergy4 = static_cast<Real>(pAction2->Energy(FALSE, 1, 0, 0, gaugefields.GetData(), NULL, NULL, NULL));

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
    appGetLattice()->m_pGaugeField[0]->FixBoundary(EFB_Field);
    pRandom->Initial(appGetLattice(), sParam);

    //Calculate condensation
    CFieldFermionWilsonSquareSU3DR* pFermion = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(appGetLattice()->GetFieldById(2));
    pFermion->InitialField(EFIT_RandomGaussian);
    pFermion->FixBoundary(EFB_Field);
    CFieldFermionWilsonSquareSU3DR* pFermion2 = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(pFermion->GetCopy());
    pFermion2->InverseD(1, 0, 0, gaugeFields.GetData(), NULL, NULL);
    pFermion2->FixBoundary(EFB_Field);
    CMeasureChiralCondensate* pCC = dynamic_cast<CMeasureChiralCondensate*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    pCC->Reset();
    pCC->OnConfigurationAcceptedZ4(1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL, pFermion, pFermion2, TRUE, TRUE);

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
        pFermion2->InverseD(1, 0, 0, gaugeFields.GetData(), NULL, NULL);
        pFermion2->FixBoundary(EFB_Field);
        pCC->Reset();
        pCC->OnConfigurationAcceptedZ4(1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL, pFermion, pFermion2, TRUE, TRUE);

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
    pF2W->InverseD(1, 0, 0, gaugeFields.GetData(), NULL, NULL);
    pFermion->PrepareForHMCNotRandomize(1, 0, gaugeFields.GetData(), NULL);

    CActionGaugePlaquetteRotating* pAction1 = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->GetActionById(1));
    CActionFermionKS* pAction2 = dynamic_cast<CActionFermionKS*>(appGetLattice()->GetActionById(2));

    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureChiralCondensateKS* pCC = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CMeasureAngularMomentumKS* pAM = dynamic_cast<CMeasureAngularMomentumKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));

    const Real fEnergy1 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL));
    pAction2->SetFermionFieldTest(pFermion);
    const Real fEnergy2 = static_cast<Real>(pAction2->Energy(FALSE, 1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL));

    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField[0]);
    pPL->OnConfigurationAccepted(1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL);
    const Real fPolyakov1 = static_cast<Real>(cuCabs(pPL->m_lstLoop[0]));
    pCC->OnConfigurationAcceptedZ4(1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL, pF1W, pF2W, TRUE, TRUE);
    const Real fChiralCond1 = _cuCabsf(pCC->m_lstCondAll[0][0]);
    const Real fConectSusp1 = _cuCabsf(pCC->m_lstCondAll[1][0]);
    pAM->OnConfigurationAcceptedZ4(1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL, pF1W, pF2W, TRUE, TRUE);
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
    pF2W->InverseD(1, 0, 0, gaugeFields.GetData(), NULL, NULL);

    pFermion2->PrepareForHMCNotRandomize(1, 0, gaugeFields.GetData(), NULL);
    const Real fEnergy3 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL));
    pAction2->SetFermionFieldTest(pFermion2);
    const Real fEnergy4 = static_cast<Real>(pAction2->Energy(FALSE, 1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL));


    pPL->OnConfigurationAccepted(1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL);
    const Real fPolyakov2 = static_cast<Real>(cuCabs(pPL->m_lstLoop[0]));
    pCC->OnConfigurationAcceptedZ4(1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL, pF1W, pF2W, TRUE, TRUE);
    const Real fChiralCond2 = _cuCabsf(pCC->m_lstCondAll[0][0]);
    const Real fConectSusp2 = _cuCabsf(pCC->m_lstCondAll[1][0]);
    pAM->OnConfigurationAcceptedZ4(1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL, pF1W, pF2W, TRUE, TRUE);
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

UINT TestGaugeFixingProjectionRecovery(CParameters& param)
{
    CParameters testConfig = param.GetParameter(_T("TestConfig"));

    CCString sProjectionType;
    if (!testConfig.FetchStringValue(_T("ProjectionType"), sProjectionType))
    {
        appCrucial(_T("TestGaugeFixingProjectionRecovery: ProjectionType not found in parameters\n"));
        return 1;
    }

    CCString sGaugeGroup;
    if (!testConfig.FetchStringValue(_T("GaugeGroup"), sGaugeGroup))
    {
        sGaugeGroup = _T("SU3");
    }

    CFieldGauge* pGauge = dynamic_cast<CFieldGauge*>(appGetLattice()->GetFieldById(1)->GetCopy());
    TArray<CFieldGauge*> gauge;
    gauge.AddItem(pGauge);
    CActionGaugePlaquette* pAction1 = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
    UINT uiError = 0;

    // Additional random-configuration check: the gauge-fixing loss should improve.
    CFieldGauge* pLossGauge = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
    const Real fRandomLossBefore = static_cast<Real>(appGetLattice()->m_pGaugeFixing->CheckRes(pLossGauge));
    appGetLattice()->m_pGaugeFixing->GaugeFixing(pLossGauge);
    const Real fRandomLossAfter = static_cast<Real>(appGetLattice()->m_pGaugeFixing->CheckRes(pLossGauge));
    appGeneral(_T("%s (%s) random loss before = %2.12f, after = %2.12f\n"),
        sProjectionType.c_str(), sGaugeGroup.c_str(), fRandomLossBefore, fRandomLossAfter);
    if (fRandomLossAfter > fRandomLossBefore + _GAUGE_FIXING_ZeroERROR)
    {
        ++uiError;
    }
    appSafeDelete(pLossGauge);

    // Step 0: CenterRemove on a copy of the random configuration.
    // Removing the nearest center element must bring the phase of every link
    // into the fundamental center basin: |arg(Tr U)| <= pi/3 for SU(3)
    // (Tr U >= 0 for SU(2)).
    if (sProjectionType != _T("MAG"))
    {
        CFieldGauge* pRemove = dynamic_cast<CFieldGauge*>(pGauge->GetCopy());
        if (NULL != pRemove)
        {
            switch (pGauge->GetFieldType())
            {
            case EFT_GaugeSU2:
                CGaugeFixingMCGDirect::CenterRemove(dynamic_cast<CFieldGaugeSU2*>(pRemove));
                break;
            case EFT_GaugeSU3:
                CGaugeFixingMCGDirect::CenterRemove(dynamic_cast<CFieldGaugeSU3*>(pRemove));
                break;
            default:
                break;
            }
            appSynchronize();

            UINT uiSize = 0;
            BYTE* pData = pRemove->CopyDataOut(uiSize);
            const CLGComplex* pComp = (const CLGComplex*)pData;
            UINT uiOutOfBasin = 0;
            UINT uiLinkCount = 0;
            if (EFT_GaugeSU3 == pGauge->GetFieldType())
            {
                uiLinkCount = uiSize / (sizeof(CLGComplex) * 9);
                for (UINT i = 0; i < uiLinkCount; ++i)
                {
                    const CLGComplex& m00 = pComp[i * 9 + 0];
                    const CLGComplex& m11 = pComp[i * 9 + 4];
                    const CLGComplex& m22 = pComp[i * 9 + 8];
                    const double dRe = static_cast<double>(m00.x + m11.x + m22.x);
                    const double dIm = static_cast<double>(m00.y + m11.y + m22.y);
                    const double dTheta = atan2(dIm, dRe);
                    if (fabs(dTheta) > 3.14159265358979323846 / 3.0 + 1.0e-4)
                    {
                        ++uiOutOfBasin;
                    }
                }
            }
            else if (EFT_GaugeSU2 == pGauge->GetFieldType())
            {
                uiLinkCount = uiSize / (sizeof(CLGComplex) * 4);
                for (UINT i = 0; i < uiLinkCount; ++i)
                {
                    const CLGComplex& m00 = pComp[i * 4 + 0];
                    const CLGComplex& m11 = pComp[i * 4 + 3];
                    if ((static_cast<double>(m00.x + m11.x)) < -1.0e-4)
                    {
                        ++uiOutOfBasin;
                    }
                }
            }
            free(pData);
            appGeneral(_T("%s (%s) center-remove out-of-basin links = %u / %u\n"),
                sProjectionType.c_str(), sGaugeGroup.c_str(), uiOutOfBasin, uiLinkCount);
            if (0 != uiOutOfBasin)
            {
                ++uiError;
            }
            appSafeDelete(pRemove);
        }
    }

    // Step 1: Projection
    switch (pGauge->GetFieldType())
    {
        case EFT_GaugeSU2:
        {
            CFieldGaugeSU2* pSU2 = dynamic_cast<CFieldGaugeSU2*>(pGauge);
            if (sProjectionType == _T("MAG"))
            {
                CGaugeFixingMAG::MaximalAbelianProjection(pSU2);
            }
            else
            {
                CGaugeFixingMCGDirect::CenterProjection(pSU2);
            }
        }
        break;
        case EFT_GaugeSU3:
        {
            CFieldGaugeSU3* pSU3 = dynamic_cast<CFieldGaugeSU3*>(pGauge);
            if (sProjectionType == _T("MAG"))
            {
                CGaugeFixingMAG::MaximalAbelianProjection(pSU3);
            }
            else
            {
                CGaugeFixingMCGDirect::CenterProjection(pSU3);
            }
        }
        break;
        default:
            appCrucial(_T("TestGaugeFixingProjectionRecovery: unsupported gauge group\n"));
            appSafeDelete(pGauge);
            return 1;
    }

    const Real fProjRes = static_cast<Real>(appGetLattice()->m_pGaugeFixing->CheckRes(pGauge));
    const Real fProjEnergy = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, 0, gauge.GetData(), NULL, NULL, NULL));
    appGeneral(_T("%s (%s) projection residual = %2.12f, projection energy = %f\n"),
        sProjectionType.c_str(), sGaugeGroup.c_str(), fProjRes, fProjEnergy);

    if (fProjRes > _GAUGE_FIXING_ZeroERROR)
    {
        ++uiError;
    }

    // Step 2: Random gauge transform
    CGaugeFixingRandom* pRandom = new CGaugeFixingRandom();
    pRandom->Initial(appGetLattice(), CParameters());
    pRandom->GaugeFixing(pGauge);
    const Real fAfterRandomEnergy = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, 0, gauge.GetData(), NULL, NULL, NULL));
    appGeneral(_T("%s (%s) after random transform energy = %f\n"),
        sProjectionType.c_str(), sGaugeGroup.c_str(), fAfterRandomEnergy);

    if (appAbs(fProjEnergy - fAfterRandomEnergy) > _GAUGE_FIXING_EnergyERROR)
    {
        ++uiError;
    }

    // Step 3: Gauge fixing should recover the projected state
    appGetLattice()->m_pGaugeFixing->GaugeFixing(pGauge);
    const Real fDivation = static_cast<Real>(appGetLattice()->m_pGaugeFixing->CheckRes(pGauge));
    const Real fAfterEnergy1 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, 0, gauge.GetData(), NULL, NULL, NULL));

    if (fDivation > _GAUGE_FIXING_ZeroERROR)
    {
        ++uiError;
    }

    if (appAbs(fAfterRandomEnergy - fAfterEnergy1) > _GAUGE_FIXING_EnergyERROR)
    {
        ++uiError;
    }

    appGeneral(_T("%s (%s) gauge fixing with divation = %f, Proj Energy = %f, After Energy = %f\n"),
        sProjectionType.c_str(), sGaugeGroup.c_str(), fDivation, fProjEnergy, fAfterEnergy1);

    pGauge->CopyTo(appGetLattice()->GetFieldById(1));
    appSafeDelete(pGauge);
    appSafeDelete(pRandom);

    return uiError;
}

UINT TestGaugeFixingCoulombPorjectivePlane2(CParameters&)
{
    UINT uiError = 0;
    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1)->GetCopy());
    TArray<CFieldGauge*> gaugeFields;
    gaugeFields.AddItem(pGauge);
    CActionGaugePlaquetteRotating* pAction1 = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->GetActionById(1));

    const Real fEnergy1 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL));

    appGetLattice()->m_pGaugeFixing->GaugeFixing(pGauge);
    const Real fError = static_cast<Real>(appGetLattice()->m_pGaugeFixing->CheckRes(pGauge));

    const Real fEnergy3 = static_cast<Real>(pAction1->Energy(FALSE, 1, 0, 0, gaugeFields.GetData(), NULL, NULL, NULL));

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

___REGIST_TEST(TestFFT, Verify, TestFFT, FFT, _TEST_NOCHECK | _TEST_MULTIGPU);

__REGIST_TEST(TestGaugeFixingLandau, GaugeFixing, TestGaugeFixingLandauCornell, LandauCornell);

//FFT not applied using single float
___REGIST_TEST(TestGaugeFixingLandau, GaugeFixing, TestGaugeFixingCoulombCornell, CoulombCornell, _TEST_MULTIGPU);

__REGIST_TEST(TestGaugeFixingLandau, GaugeFixing, TestGaugeFixingLandauLosAlamos, LandauLosAlamos);

__REGIST_TEST(TestGaugeFixingLandau, GaugeFixing, TestGaugeFixingCoulombLosAlamos, CoulombLosAlamos);

//TestGaugeFixingCoulombCornellDR has problem with debug and single point
__REGIST_TEST(TestGaugeFixingCoulombDR, GaugeFixing, TestGaugeFixingCoulombCornellDR, CoulombCornellDR);

__REGIST_TEST(TestGaugeFixingCoulombDR, GaugeFixing, TestGaugeFixingCoulombLosAlamosDR, CoulombLosAlamosDR);

__REGIST_TEST(TestGaugeFixingCoulombDRChiral, GaugeFixing, TestGaugeFixingCoulombDRChiral, CoulombDRChiral);

__REGIST_TEST(TestGaugeFixingCoulombPorjectivePlane, GaugeFixing, TestGaugeFixingRotationKS, RotationKS);

__REGIST_TEST(TestGaugeFixingCoulombPorjectivePlane2, GaugeFixing, TestGaugeFixingRotationKS2, RotationKS2);

__REGIST_TEST(TestGaugeFixingProjectionRecovery, GaugeFixing, TestGaugeFixingMCGDirect, MCGDirect);

__REGIST_TEST(TestGaugeFixingProjectionRecovery, GaugeFixing, TestGaugeFixingMCGDirect_SU2, MCGDirect_SU2);

__REGIST_TEST(TestGaugeFixingProjectionRecovery, GaugeFixing, TestGaugeFixingMAG_SU2, MAG_SU2);

__REGIST_TEST(TestGaugeFixingProjectionRecovery, GaugeFixing, TestGaugeFixingMAG, MAG);

__REGIST_TEST(TestGaugeFixingProjectionRecovery, GaugeFixing, TestGaugeFixingMCGIndirect, MCGIndirect);

__REGIST_TEST(TestGaugeFixingProjectionRecovery, GaugeFixing, TestGaugeFixingMCGIndirect_SU2, MCGIndirect_SU2);

__REGIST_TEST(TestGaugeFixingProjectionRecovery, GaugeFixing, TestGaugeFixingMCGIndirect_Standard, MCGIndirect_Standard);

__REGIST_TEST(TestGaugeFixingProjectionRecovery, GaugeFixing, TestGaugeFixingMCGIndirect_Standard_SU2, MCGIndirect_Standard_SU2);


//=============================================================================
// END OF FILE
//=============================================================================
