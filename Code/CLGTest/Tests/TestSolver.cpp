//=============================================================================
// FILENAME : TestSolver.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [02/04/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestSolver(CParameters& params)
{
    UINT uiError = 0;
    Real fMaxError = F(0.0001);
    params.FetchValueReal(_T("ExpectedErr"), fMaxError);

    CField* pField = appGetLattice()->GetFieldById(2);
    CFieldFermionWilsonSquareSU3* pFermion = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pField);
    CFieldFermionWilsonSquareSU3* pResult1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());
    const Real fLengthOfPhi = pResult1->DotReal(pResult1).x;
    pResult1->D(_FIELDS);
    const Real fLengthOfDPhi = pResult1->DotReal(pResult1).x;
    pResult1->ApplyOperator(EFO_F_InverseD, _FIELDS);
    pResult1->AxpyMinus(pFermion);
    const Real fError1 = _cuCabsf(pResult1->DotReal(pResult1));
    appGeneral(_T("| phi |^2 = %8.18f;  | D phi |^2 = %8.18f\n"), fLengthOfPhi, fLengthOfDPhi);
    appGeneral(_T("| D^-1 D phi - phi |^2 =%8.18f\n"), fError1);
    if (appAbs(fError1) > fMaxError)
    {
        ++uiError;
    }
    
    CFieldFermionWilsonSquareSU3* pResult2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult2->ApplyOperator(EFO_F_InverseD, _FIELDS);
    pResult2->D(_FIELDS);
    pResult2->AxpyMinus(pFermion);
    const Real fError2 = _cuCabsf(pResult2->DotReal(pResult2));
    appGeneral(_T("| phi |^2 = %8.18f;  | D phi |^2 = %8.18f\n"), fLengthOfPhi, fLengthOfDPhi);
    appGeneral(_T("| D D^-1 phi - phi |^2 =%8.18f\n"), fError2);
    if (appAbs(fError2) > fMaxError)
    {
        ++uiError;
    }

    CFieldFermionWilsonSquareSU3* pResult3 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult3->Ddagger(_FIELDS);
    const Real fLengthOfDdaggerPhi = pResult1->DotReal(pResult1).x;
    pResult3->ApplyOperator(EFO_F_InverseDdagger, _FIELDS);
    pResult3->AxpyMinus(pFermion);
    const Real fError3 = _cuCabsf(pResult3->DotReal(pResult3));
    appGeneral(_T("| phi |^2 = %8.18f;  | D+ phi |^2 = %8.18f\n"), fLengthOfPhi, fLengthOfDdaggerPhi);
    appGeneral(_T("| D+^-1 D+ phi - phi |^2 =%8.18f\n"), fError3);
    if (appAbs(fError3) > fMaxError)
    {
        ++uiError;
    }

    CFieldFermionWilsonSquareSU3* pResult4 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult4->ApplyOperator(EFO_F_InverseDdagger, _FIELDS);
    pResult4->Ddagger(_FIELDS);
    pResult4->AxpyMinus(pFermion);
    const Real fError4 = _cuCabsf(pResult4->DotReal(pResult4));
    appGeneral(_T("| phi |^2 = %8.18f;  | D+ phi |^2 = %8.18f\n"), fLengthOfPhi, fLengthOfDdaggerPhi);
    appGeneral(_T("| D+ D+^-1 phi - phi |^2 =%8.18f\n"), fError4);
    if (appAbs(fError4) > fMaxError)
    {
        ++uiError;
    }

    CFieldFermionWilsonSquareSU3* pResult5 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult5->DDdagger(_FIELDS);
    const Real fDDdaggerPhi = pResult5->DotReal(pResult5).x;
    pResult5->ApplyOperator(EFO_F_InverseDDdagger, _FIELDS);
    pResult5->AxpyMinus(pFermion);
    const Real fError5 = _cuCabsf(pResult5->DotReal(pResult5));
    appGeneral(_T("| phi |^2 = %8.18f;  | DD+ phi |^2 = %8.18f\n"), fLengthOfPhi, fDDdaggerPhi);
    appGeneral(_T("| (DD+)^-1 (DD+) phi - phi |^2 =%8.18f\n"), fError5);
    if (appAbs(fError5) > fMaxError)
    {
        ++uiError;
    }

    CFieldFermionWilsonSquareSU3* pResult6 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult6->ApplyOperator(EFO_F_InverseDDdagger, _FIELDS);
    pResult6->DDdagger(_FIELDS);
    pResult6->AxpyMinus(pFermion);
    const Real fError6 = _cuCabsf(pResult6->DotReal(pResult6));
    appGeneral(_T("| phi |^2 = %8.18f;  | DD+ phi |^2 = %8.18f\n"), fLengthOfPhi, fDDdaggerPhi);
    appGeneral(_T("| (DD+) (DD+)^-1 phi - phi |^2 =%8.18f\n"), fError6);
    if (appAbs(fError6) > fMaxError)
    {
        ++uiError;
    }
    
    return uiError;
}

___REGIST_TEST(TestSolver, Solver, TestSolverBiCGStab, BiCGStab, _TEST_BOUND);
___REGIST_TEST(TestSolver, Solver, TestSolverGMRES, GMRES, _TEST_BOUND);
//__REGIST_TEST(TestSolver, Solver, TestSolverGCR); //slow solver not used
___REGIST_TEST(TestSolver, Solver, TestSolverGCRODR, GCRODR, _TEST_BOUND);
//__REGIST_TEST(TestSolver, Solver, TestSolverTFQMR); //slow solver not used
___REGIST_TEST(TestSolver, Solver, TestSolverGMRESLowMode, GMRES2, _TEST_BOUND);
___REGIST_TEST(TestSolver, Solver, TestSolverGCRODRLowMode, GCRODR2, _TEST_BOUND);
___REGIST_TEST(TestSolver, Solver, TestSolverBiCGStabLowMode, BiCGStab2, _TEST_BOUND | _TEST_DOUBLE);

__REGIST_TEST(TestSolver, Solver, TestEOSolverBiCGStab, BiCGStabEO);
__REGIST_TEST(TestSolver, Solver, TestEOSolverGMRES, GMRESEO);

//__REGIST_TEST(TestSolver, Solver, TestEOSolverBiCGStabD);
//__REGIST_TEST(TestSolver, Solver, TestEOSolverGMRESD);
//__REGIST_TEST(TestSolver, Solver, TestEOSolverBiCGStabDR);
//__REGIST_TEST(TestSolver, Solver, TestEOSolverGMRESDR);

//=============================================================================
// END OF FILE
//=============================================================================
