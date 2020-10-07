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
    constants.AddItem(_make_cuComplex(F(1.0), F(-1.0)));
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

UINT TestMultiShiftSolverKS(CParameters& params)
{
    TArray<Real> oneOver2;
    TArray<Real> _oneOver2;
    TArray<Real> oneOver4;
    TArray<Real> _oneOver4;

    params.FetchValueArrayReal(_T("OneOver2"), oneOver2);
    params.FetchValueArrayReal(_T("MinusOneOver2"), _oneOver2);
    params.FetchValueArrayReal(_T("OneOver4"), oneOver4);
    params.FetchValueArrayReal(_T("MinusOneOver4"), _oneOver4);
    Real fMaxError = F(0.000001);
    params.FetchValueReal(_T("ExpectedErr"), fMaxError);

    CRatinalApproximation R1Over2(oneOver2);
    CRatinalApproximation R_1Over2(_oneOver2);
    CRatinalApproximation R1Over4(oneOver4);
    CRatinalApproximation R_1Over4(_oneOver4);

    const CField* pField = appGetLattice()->GetFieldById(2);
    const CField* pGauge = appGetLattice()->m_pGaugeField;
    const Real fLengthOfPhi = pField->Dot(pField).x;
    CFieldFermion* pFieldCopy = dynamic_cast<CFieldFermion*>(pField->GetCopy());

    pFieldCopy->RationalApproximation(EFO_F_D, pGauge, &R1Over2);
    Real fLength2 = pFieldCopy->Dot(pFieldCopy).x;
    pFieldCopy->RationalApproximation(EFO_F_D, pGauge, &R_1Over2);
    pFieldCopy->AxpyMinus(pField);
    Real fLength3 = pFieldCopy->Dot(pFieldCopy).x;
    UINT uiError = 0;
    if (fLength3 > fMaxError)
    {
        uiError++;
    }
    appGeneral(_T("|phi|^2 = %2.18f, |D^{1/2}phi|^2 = %2.18f, |D^{-1/2}D^{1/2}phi-phi|^2 = %2.18f\n"), fLengthOfPhi, fLength2, fLength3);

    pField->CopyTo(pFieldCopy);

    pFieldCopy->RationalApproximation(EFO_F_D, pGauge, &R1Over4);
    fLength2 = pFieldCopy->Dot(pFieldCopy).x;
    pFieldCopy->RationalApproximation(EFO_F_D, pGauge, &R_1Over4);
    pFieldCopy->AxpyMinus(pField);
    fLength3 = pFieldCopy->Dot(pFieldCopy).x;
    if (fLength3 > fMaxError)
    {
        uiError++;
    }
    appGeneral(_T("|phi|^2 = %2.18f, |D^{1/4}phi|^2 = %2.18f, |D^{-1/4}D^{1/4}phi-phi|^2 = %2.18f\n"), fLengthOfPhi, fLength2, fLength3);

    pField->CopyTo(pFieldCopy);

    pFieldCopy->RationalApproximation(EFO_F_DDdagger, pGauge, &R1Over2);
    fLength2 = pFieldCopy->Dot(pFieldCopy).x;
    pFieldCopy->RationalApproximation(EFO_F_DDdagger, pGauge, &R_1Over2);
    pFieldCopy->AxpyMinus(pField);
    fLength3 = pFieldCopy->Dot(pFieldCopy).x;
    if (fLength3 > fMaxError)
    {
        uiError++;
    }
    appGeneral(_T("|phi|^2 = %2.18f, |DD^{1/2}phi|^2 = %2.18f, |DD^{-1/2}DD^{1/2}phi-phi|^2 = %2.18f\n"), fLengthOfPhi, fLength2, fLength3);

    pField->CopyTo(pFieldCopy);

    pFieldCopy->RationalApproximation(EFO_F_DDdagger, pGauge, &R1Over4);
    fLength2 = pFieldCopy->Dot(pFieldCopy).x;
    pFieldCopy->RationalApproximation(EFO_F_DDdagger, pGauge, &R_1Over4);
    pFieldCopy->AxpyMinus(pField);
    fLength3 = pFieldCopy->Dot(pFieldCopy).x;
    if (fLength3 > fMaxError)
    {
        uiError++;
    }
    appGeneral(_T("|phi|^2 = %2.18f, |DD^{1/4}phi|^2 = %2.18f, |DD^{-1/4}DD^{1/4}phi-phi|^2 = %2.18f\n"), fLengthOfPhi, fLength2, fLength3);


    //It seems that the rational approximation is very poor
    //Only when gauge field is set to be 1(cold), the result will be close
    //Maybe it is because the matrix is not sparse enough? we test on 8x8x8x8
    //pField->CopyTo(pFieldCopy);
    //pFieldCopy->RationalApproximation(EFO_F_D, pGauge, &R_1Over2);
    //pFieldCopy->RationalApproximation(EFO_F_D, pGauge, &R_1Over2);
    //pFieldCopy->D(pGauge);
    //fLength2 = pFieldCopy->Dot(pFieldCopy).x;
    //pFieldCopy->AxpyMinus(pField);
    //fLength3 = pFieldCopy->Dot(pFieldCopy).x;
    //appGeneral(_T("|phi|^2 = %2.18f, |DD+(DD^{-1/2})^2phi|^2 = %2.18f, |DD+(DD^{-1/2})^2phi-phi|^2 = %2.18f\n"), fLengthOfPhi, fLength2, fLength3);
    //if (fLength3 > fMaxError)
    //{
    //    uiError++;
    //}

    appSafeDelete(pFieldCopy);
    return uiError;
}

__REGIST_TEST(TestMultiShiftSolver, Solver, TestMSSolverGMRES);
__REGIST_TEST(TestMultiShiftSolver, Solver, TestMSSolverFOM);
__REGIST_TEST(TestMultiShiftSolver, Solver, TestMSSolverBiCGStab);

__REGIST_TEST(TestMultiShiftSolverKS, Solver, TestMSKSSolverGMRES);
__REGIST_TEST(TestMultiShiftSolverKS, Solver, TestMSKSSolverFOM);
__REGIST_TEST(TestMultiShiftSolverKS, Solver, TestMSKSSolverBiCGStab);
__REGIST_TEST(TestMultiShiftSolverKS, Solver, TestMSKSSolverNested);


