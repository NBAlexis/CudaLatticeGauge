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

UINT TestGaugeFixingCoulombCornellDR(CParameters&)
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
    appGeneral(_T("Divation = %f\n"), fError);
    appGeneral(_T("Before Energy1 = %f, Energy2 = %f\n"), fEnergy1, fEnergy2);
    appGeneral(_T("After Energy1 = %f, Energy2 = %f\n"), fEnergy3, fEnergy4);

    appSafeDelete(pGauge);
    appSafeDelete(pFermion2);

    return uiError;
}

__REGIST_TEST(TestFFT, Misc, TestFFT);

__REGIST_TEST(TestGaugeFixingLandau, Misc, TestGaugeFixingLandauCornell);

__REGIST_TEST(TestGaugeFixingLandau, Misc, TestGaugeFixingCoulombCornell);

//__REGIST_TEST(TestGaugeFixingLandau, Misc, TestGaugeFixingLandauLosAlamos);

__REGIST_TEST(TestGaugeFixingCoulombCornellDR, Misc, TestGaugeFixingCoulombCornellDR);


//=============================================================================
// END OF FILE
//=============================================================================
