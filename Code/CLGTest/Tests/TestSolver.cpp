//=============================================================================
// FILENAME : TestSolver.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [02/04/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestSolver(CParameters& )
{
    UINT uiError = 0;
    Real fMaxError = F(0.000001);
    CField* pField = appGetLattice()->GetFieldById(2);
    CFieldFermionWilsonSquareSU3* pFermion = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pField);
    //Real fLength = _cuCabsf(pFermion->Dot(pFermion));

    CFieldFermionWilsonSquareSU3* pResult1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult1->D(appGetLattice()->m_pGaugeField);
    pResult1->ApplyOperator(EFO_F_InverseD, appGetLattice()->m_pGaugeField);
    pResult1->AxpyMinus(pFermion);
    Real fError1 = _cuCabsf(pResult1->Dot(pResult1));
    appGeneral(_T("| D^-1 D phi - phi |^2 =%8.18f\n"), fError1);
    if (appAbs(fError1) > fMaxError)
    {
        ++uiError;
    }

    CFieldFermionWilsonSquareSU3* pResult2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult2->ApplyOperator(EFO_F_InverseD, appGetLattice()->m_pGaugeField);
    pResult2->D(appGetLattice()->m_pGaugeField);
    pResult2->AxpyMinus(pFermion);
    Real fError2 = _cuCabsf(pResult2->Dot(pResult2));
    appGeneral(_T("| D D^-1 phi - phi |^2 =%8.18f\n"), fError2);
    if (appAbs(fError2) > fMaxError)
    {
        ++uiError;
    }

    CFieldFermionWilsonSquareSU3* pResult3 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult3->Ddagger(appGetLattice()->m_pGaugeField);
    pResult3->ApplyOperator(EFO_F_InverseDdagger, appGetLattice()->m_pGaugeField);
    pResult3->AxpyMinus(pFermion);
    Real fError3 = _cuCabsf(pResult3->Dot(pResult3));
    appGeneral(_T("| D+^-1 D+ phi - phi |^2 =%8.18f\n"), fError3);
    if (appAbs(fError3) > fMaxError)
    {
        ++uiError;
    }

    CFieldFermionWilsonSquareSU3* pResult4 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult4->ApplyOperator(EFO_F_InverseDdagger, appGetLattice()->m_pGaugeField);
    pResult4->Ddagger(appGetLattice()->m_pGaugeField);
    pResult4->AxpyMinus(pFermion);
    Real fError4 = _cuCabsf(pResult4->Dot(pResult4));
    appGeneral(_T("| D+ D+^-1 phi - phi |^2 =%8.18f\n"), fError4);
    if (appAbs(fError4) > fMaxError)
    {
        ++uiError;
    }

    CFieldFermionWilsonSquareSU3* pResult5 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult5->DDdagger(appGetLattice()->m_pGaugeField);
    pResult5->ApplyOperator(EFO_F_InverseDDdagger, appGetLattice()->m_pGaugeField);
    pResult5->AxpyMinus(pFermion);
    Real fError5 = _cuCabsf(pResult5->Dot(pResult5));
    appGeneral(_T("| (DD+)^-1 (DD+) phi - phi |^2 =%8.18f\n"), fError5);
    if (appAbs(fError5) > fMaxError)
    {
        ++uiError;
    }

    CFieldFermionWilsonSquareSU3* pResult6 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult6->ApplyOperator(EFO_F_InverseDDdagger, appGetLattice()->m_pGaugeField);
    pResult6->DDdagger(appGetLattice()->m_pGaugeField);
    pResult6->AxpyMinus(pFermion);
    Real fError6 = _cuCabsf(pResult6->Dot(pResult6));
    appGeneral(_T("| (DD+) (DD+)^-1 phi - phi |^2 =%8.18f\n"), fError6);
    if (appAbs(fError6) > fMaxError)
    {
        ++uiError;
    }

    return uiError;
}

__REGIST_TEST(TestSolver, Solver, TestSolverBiCGStab);

__REGIST_TEST(TestSolver, Solver, TestSolverGMRES);

__REGIST_TEST(TestSolver, Solver, TestSolverGCR);

//__REGIST_TEST(TestSolver, Solver, TestSolverBiCGStabLowMode);

//__REGIST_TEST(TestSolver, Solver, TestSolverGMRESLowMode);

//__REGIST_TEST(TestSolver, Solver, TestSolverGCRLowMode);

//=============================================================================
// END OF FILE
//=============================================================================
