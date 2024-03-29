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
    pResult1->D(appGetLattice()->m_pGaugeField);
    const Real fLengthOfDPhi = pResult1->DotReal(pResult1).x;
    pResult1->ApplyOperator(EFO_F_InverseD, appGetLattice()->m_pGaugeField);
    pResult1->AxpyMinus(pFermion);
    const Real fError1 = _cuCabsf(pResult1->DotReal(pResult1));
    appGeneral(_T("| phi |^2 = %8.18f;  | D phi |^2 = %8.18f\n"), fLengthOfPhi, fLengthOfDPhi);
    appGeneral(_T("| D^-1 D phi - phi |^2 =%8.18f\n"), fError1);
    if (appAbs(fError1) > fMaxError)
    {
        ++uiError;
    }
    
    CFieldFermionWilsonSquareSU3* pResult2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult2->ApplyOperator(EFO_F_InverseD, appGetLattice()->m_pGaugeField);
    pResult2->D(appGetLattice()->m_pGaugeField);
    pResult2->AxpyMinus(pFermion);
    const Real fError2 = _cuCabsf(pResult2->DotReal(pResult2));
    appGeneral(_T("| phi |^2 = %8.18f;  | D phi |^2 = %8.18f\n"), fLengthOfPhi, fLengthOfDPhi);
    appGeneral(_T("| D D^-1 phi - phi |^2 =%8.18f\n"), fError2);
    if (appAbs(fError2) > fMaxError)
    {
        ++uiError;
    }

    CFieldFermionWilsonSquareSU3* pResult3 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult3->Ddagger(appGetLattice()->m_pGaugeField);
    const Real fLengthOfDdaggerPhi = pResult1->DotReal(pResult1).x;
    pResult3->ApplyOperator(EFO_F_InverseDdagger, appGetLattice()->m_pGaugeField);
    pResult3->AxpyMinus(pFermion);
    const Real fError3 = _cuCabsf(pResult3->DotReal(pResult3));
    appGeneral(_T("| phi |^2 = %8.18f;  | D+ phi |^2 = %8.18f\n"), fLengthOfPhi, fLengthOfDdaggerPhi);
    appGeneral(_T("| D+^-1 D+ phi - phi |^2 =%8.18f\n"), fError3);
    if (appAbs(fError3) > fMaxError)
    {
        ++uiError;
    }

    CFieldFermionWilsonSquareSU3* pResult4 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult4->ApplyOperator(EFO_F_InverseDdagger, appGetLattice()->m_pGaugeField);
    pResult4->Ddagger(appGetLattice()->m_pGaugeField);
    pResult4->AxpyMinus(pFermion);
    const Real fError4 = _cuCabsf(pResult4->DotReal(pResult4));
    appGeneral(_T("| phi |^2 = %8.18f;  | D+ phi |^2 = %8.18f\n"), fLengthOfPhi, fLengthOfDdaggerPhi);
    appGeneral(_T("| D+ D+^-1 phi - phi |^2 =%8.18f\n"), fError4);
    if (appAbs(fError4) > fMaxError)
    {
        ++uiError;
    }

    CFieldFermionWilsonSquareSU3* pResult5 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult5->DDdagger(appGetLattice()->m_pGaugeField);
    const Real fDDdaggerPhi = pResult5->DotReal(pResult5).x;
    pResult5->ApplyOperator(EFO_F_InverseDDdagger, appGetLattice()->m_pGaugeField);
    pResult5->AxpyMinus(pFermion);
    const Real fError5 = _cuCabsf(pResult5->DotReal(pResult5));
    appGeneral(_T("| phi |^2 = %8.18f;  | DD+ phi |^2 = %8.18f\n"), fLengthOfPhi, fDDdaggerPhi);
    appGeneral(_T("| (DD+)^-1 (DD+) phi - phi |^2 =%8.18f\n"), fError5);
    if (appAbs(fError5) > fMaxError)
    {
        ++uiError;
    }

    CFieldFermionWilsonSquareSU3* pResult6 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion->GetCopy());

    pResult6->ApplyOperator(EFO_F_InverseDDdagger, appGetLattice()->m_pGaugeField);
    pResult6->DDdagger(appGetLattice()->m_pGaugeField);
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

__REGIST_TEST(TestSolver, Solver, TestSolverBiCGStab);

__REGIST_TEST(TestSolver, Solver, TestSolverGMRES);

__REGIST_TEST(TestSolver, Solver, TestSolverGCR);

__REGIST_TEST(TestSolver, Solver, TestSolverGCRODR);

__REGIST_TEST(TestSolver, Solver, TestSolverTFQMR);


__REGIST_TEST(TestSolver, Solver, TestSolverGMRESLowMode);

__REGIST_TEST(TestSolver, Solver, TestSolverGCRODRLowMode);

#if _CLG_DOUBLEFLOAT

__REGIST_TEST(TestSolver, Solver, TestSolverBiCGStabLowMode);

#endif

__REGIST_TEST(TestSolver, Solver, TestEOSolverBiCGStab);

__REGIST_TEST(TestSolver, Solver, TestEOSolverGMRES);

__REGIST_TEST(TestSolver, Solver, TestEOSolverBiCGStabD);

__REGIST_TEST(TestSolver, Solver, TestEOSolverGMRESD);

__REGIST_TEST(TestSolver, Solver, TestEOSolverBiCGStabDR);

__REGIST_TEST(TestSolver, Solver, TestEOSolverGMRESDR);

//=============================================================================
// END OF FILE
//=============================================================================
