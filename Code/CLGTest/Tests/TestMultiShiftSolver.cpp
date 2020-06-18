//=============================================================================
// FILENAME : TestMultiShiftSolver.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [16/06/2020 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestMultiShiftSolver(CParameters& params)
{
    UINT uiError = 0;
    Real fMaxError = F(0.000001);
    params.FetchValueReal(_T("ExpectedErr"), fMaxError);

    CMultiShiftSolver* pSolver = appGetMultiShiftSolver(2);
    const CField* pField = appGetLattice()->GetFieldById(2);
    const Real fLengthOfPhi = pField->Dot(pField).x;
    const CFieldFermionWilsonSquareSU3* pFermion = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(pField);
    CFieldFermionWilsonSquareSU3* pTemp = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pField->GetCopy());
    TArray<CLGComplex> constants;
    constants.AddItem(_make_cuComplex(F(2.0), F(0.0)));
    constants.AddItem(_make_cuComplex(F(1.0), F(0.0)));
    constants.AddItem(_make_cuComplex(F(0.5), F(0.0)));
    constants.AddItem(_make_cuComplex(F(0.0), F(1.0)));
    constants.AddItem(_make_cuComplex(F(0.5), F(0.5)));
    TArray<CField*> resultFields;
    for (INT i = 0; i < constants.Num(); ++i)
    {
        CFieldFermionWilsonSquareSU3* pResult = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pField->GetCopy());
        resultFields.AddItem(pResult);
    }

    //Result = (D + cn)^{-1} pField
    pSolver->Solve(resultFields, constants, pFermion, appGetLattice()->m_pGaugeField, EFO_F_D);
    for (INT i = 0; i < constants.Num(); ++i)
    {
        // Result = (D + cn) Result
        resultFields[i]->CopyTo(pTemp);
        resultFields[i]->ApplyOperator(EFO_F_D, appGetLattice()->m_pGaugeField);
        pTemp->ScalarMultply(constants[i]);
        resultFields[i]->AxpyPlus(pTemp);
    }

    for (INT i = 0; i < constants.Num(); ++i)
    {
        resultFields[i]->AxpyMinus(pFermion);
        const Real fError1 = _cuCabsf(resultFields[i]->Dot(resultFields[i]));
        appGeneral(_T("| D D^-1 phi - phi |^2=%8.18f, |phi|^2=%8.18f\n"), fError1, fLengthOfPhi);
        if (appAbs(fError1) > fMaxError)
        {
            ++uiError;
        }
    }


    pSolver->Solve(resultFields, constants, pFermion, appGetLattice()->m_pGaugeField, EFO_F_Ddagger);
    for (INT i = 0; i < constants.Num(); ++i)
    {
        // Result = (D + cn) Result
        resultFields[i]->CopyTo(pTemp);
        resultFields[i]->ApplyOperator(EFO_F_Ddagger, appGetLattice()->m_pGaugeField);
        pTemp->ScalarMultply(constants[i]);
        resultFields[i]->AxpyPlus(pTemp);
    }

    for (INT i = 0; i < constants.Num(); ++i)
    {
        resultFields[i]->AxpyMinus(pFermion);
        const Real fError1 = _cuCabsf(resultFields[i]->Dot(resultFields[i]));
        appGeneral(_T("| D+ D+^-1 phi - phi |^2=%8.18f, |phi|^2=%8.18f\n"), fError1, fLengthOfPhi);
        if (appAbs(fError1) > fMaxError)
        {
            ++uiError;
        }
    }


    pSolver->Solve(resultFields, constants, pFermion, appGetLattice()->m_pGaugeField, EFO_F_DDdagger);
    for (INT i = 0; i < constants.Num(); ++i)
    {
        // Result = (D + cn) Result
        resultFields[i]->CopyTo(pTemp);
        resultFields[i]->ApplyOperator(EFO_F_DDdagger, appGetLattice()->m_pGaugeField);
        pTemp->ScalarMultply(constants[i]);
        resultFields[i]->AxpyPlus(pTemp);
    }

    for (INT i = 0; i < constants.Num(); ++i)
    {
        resultFields[i]->AxpyMinus(pFermion);
        const Real fError1 = _cuCabsf(resultFields[i]->Dot(resultFields[i]));
        appGeneral(_T("| DD+ (DD+)^-1 phi - phi |^2=%8.18f, |phi|^2=%8.18f\n"), fError1, fLengthOfPhi);
        if (appAbs(fError1) > fMaxError)
        {
            ++uiError;
        }
    }


    appSafeDelete(pTemp);
    for (INT i = 0; i < constants.Num(); ++i)
    {
        appSafeDelete(resultFields[i]);
    }

    return uiError;
}

__REGIST_TEST(TestMultiShiftSolver, Solver, TestMSSolverGMRES);
__REGIST_TEST(TestMultiShiftSolver, Solver, TestMSSolverFOM);

