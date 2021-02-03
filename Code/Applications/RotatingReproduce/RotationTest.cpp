//=============================================================================
// FILENAME : RotationTest.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [02/03/2021 nbale]
//=============================================================================

#include "RotatingReproduce.h"

void TestPlaqutteEnergy()
{
    const CFieldGaugeSU3* pGuageField = dynamic_cast<const CFieldGaugeSU3 * >(appGetLattice()->m_pGaugeField);
    CFieldGaugeSU3* pCopy = dynamic_cast<CFieldGaugeSU3*>(pGuageField->GetCopy());

    pGuageField->CalculateOnlyStaple(pCopy);

    const Real fE1 = pGuageField->CalculatePlaqutteEnergy(F(5.0));
    const Real fE2 = pGuageField->CalculatePlaqutteEnergyUseClover(F(5.0));
    const Real fE3 = pGuageField->CalculatePlaqutteEnergyUsingStable(F(5.0), pCopy);

    appGeneral(_T("============ TestPlaqutteEnergy ===========\n"));
    appGeneral(_T("E1: %f, E2: %f, E3: %f\n"), fE1, fE2, fE3);
    appGeneral(_T("===========================================\n"));

    appSafeDelete(pCopy);
}

void TestXYTerm()
{
    CActionGaugePlaquetteRotating* pActionR = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->m_pActionList[0]);
    const Real fXY1 = pActionR->XYTerm1(appGetLattice()->m_pGaugeField);
    const Real fXY2 = pActionR->XYTerm2(appGetLattice()->m_pGaugeField);

    appGeneral(_T("================= TestXYTerm ==============\n"));
    appGeneral(_T("E1: %f, E2: %f\n"), fXY1, fXY2);
    appGeneral(_T("===========================================\n"));
}

void TestSigma12()
{
    CFieldFermionWilsonSquareSU3* pFermionCopy = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetFieldById(2)->GetCopy());
    CFieldFermionWilsonSquareSU3* pFermionCopy2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetFieldById(2)->GetCopy());

    pFermionCopy->ApplyGamma(GAMMA2);
    pFermionCopy->ApplyGamma(GAMMA1);
    pFermionCopy->ScalarMultply(_imgc);

    pFermionCopy2->ApplyGamma(SIGMA12E);

    const DOUBLE fAmp1 = pFermionCopy2->Dot(pFermionCopy2).x;
    pFermionCopy2->AxpyMinus(pFermionCopy);
    const DOUBLE fAmp2 = pFermionCopy2->Dot(pFermionCopy2).x;

    appGeneral(_T("================= Sigma12 Test ==============\n"));
    appGeneral(_T("Amp: %f, Substracted: %f\n"), fAmp1, fAmp2);
    appGeneral(_T("===========================================\n"));

    appSafeDelete(pFermionCopy);
    appSafeDelete(pFermionCopy2);
}


INT RotationTest(CParameters& params)
{
    appSetupLog(params);

    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    TestPlaqutteEnergy();

    TestXYTerm();

    TestSigma12();

    appQuitCLG();

    return 0;
}


