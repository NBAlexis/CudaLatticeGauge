//=============================================================================
// FILENAME : CSolverCG.cpp
// 
// DESCRIPTION:
// Conjugate gradient solver for Hermitian positive definite operators.
//
// REVISION:
//  [mm/dd/yy]
//  [07/24/2026 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CSolverCG.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CSLASolverCG)

CSLASolverCG::CSLASolverCG()
    : CSLASolver()
    , m_uiMaxStep(5000)
{

}

void CSLASolverCG::Configurate(const CParameters& param)
{
    INT iValue;
    Real fValue;

    if (param.FetchValueINT(_T("MaxStep"), iValue))
    {
        m_uiMaxStep = static_cast<UINT>(iValue);
    }
    if (param.FetchValueINT(_T("AbsoluteAccuracy"), iValue))
    {
        m_bAbsoluteAccuracy = (0 != iValue);
    }
    if (param.FetchValueReal(_T("Accuracy"), fValue))
    {
        m_fAccuracy = fValue;
        if (m_fAccuracy < _CLG_FLT_EPSILON * F(2.0))
        {
            m_fAccuracy = _CLG_FLT_EPSILON * F(2.0);
            appGeneral(_T("Solver accuracy too small (%2.18f), set to be %2.18f\n"), fValue, m_fAccuracy);
        }
    }
}

void CSLASolverCG::AllocateBuffers(const CField*)
{

}

UBOOL CSLASolverCG::IterateCG(CField* pFieldX, CField* pR, const CField* pFieldB,
    INT gaugeNum, INT bosonNum, INT tensor2Num,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
    EFieldOperator uiM)
{
    CField* pP = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId, _T(__FILE__), __LINE__);
    CField* pW = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId, _T(__FILE__), __LINE__);

    //use it to estimate relative error
    Real fBLength = F(1.0);
    if (!m_bAbsoluteAccuracy)
    {
        fBLength = static_cast<Real>(_sqrt(pFieldB->GetLength()));
    }
    appParanoiac(_T("-- CSLASolverCG::Solve start operator: %s-- fLength = %f --\n"), __ENUM_TO_STRING(EFieldOperator, uiM).c_str(), fBLength);

    pR->CopyTo(pP);
    //gamma = (r, r)
    Real fRNorm2 = static_cast<Real>(pR->GetLength());
    UINT uiIter = 0;
    UBOOL bConverged = (static_cast<Real>(_sqrt(fRNorm2)) < m_fAccuracy * fBLength);

    while (!bConverged && uiIter < m_uiMaxStep)
    {
        //w = A p
        pP->CopyTo(pW);
        pW->ApplyOperator(uiM, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields);
        const Real fPAQ = static_cast<Real>(pP->Dot(pW).x);
        if (fPAQ <= _CLG_FLT_MIN_)
        {
            appGeneral(_T("CSolverCG breakdown: (p, Ap) = %f\n"), fPAQ);
            break;
        }
        const Real fAlpha = fRNorm2 / fPAQ;

        //x = x + alpha * p
        pFieldX->Axpy(fAlpha, pP);
        //r = r - alpha * w
        pR->Axpy(-fAlpha, pW);

        const Real fRNorm2New = static_cast<Real>(pR->GetLength());
        const Real fBeta = fRNorm2New / fRNorm2;
        //p = r + beta * p
        pP->ScalarMultply(fBeta);
        pP->AxpyPlus(pR);

        fRNorm2 = fRNorm2New;
        ++uiIter;
        const Real fError = static_cast<Real>(_sqrt(fRNorm2));
        appParanoiac(_T("CSolverCG::Solve deviation: ---- iter %d = %8.15f\n"), uiIter, fError);
        bConverged = (fError < m_fAccuracy * fBLength);
    }

    if (!bConverged)
    {
        appGeneral(_T("CSolverCG::Solve failed to converge, %d iterations, deviation = %8.15f\n"), uiIter, static_cast<Real>(_sqrt(fRNorm2)));
    }
    else
    {
        appParanoiac(_T("CSolverCG::Solve converged, %d iterations, deviation = %8.15f\n"), uiIter, static_cast<Real>(_sqrt(fRNorm2)));
        appParanoiac(_T("CSolverCG::Solve finished, %d iterations, converged = %d\n"), uiIter, bConverged ? 1 : 0);
    }
    pP->Return();
    pW->Return();
    return bConverged;
}

UBOOL CSLASolverCG::SolveHPD(CField* pFieldX, const CField* pFieldB,
    INT gaugeNum, INT bosonNum, INT tensor2Num,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
    EFieldOperator uiM, const CField* pStart)
{
    _RECORD(CSLASolverCG::SolveHPD);
    CField* pR = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId, _T(__FILE__), __LINE__);
    CField* pW = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId, _T(__FILE__), __LINE__);

    //save b first, pFieldX may be the same field as pFieldB
    pFieldB->CopyTo(pR);

    //set initial x
    if (NULL != pStart)
    {
        pStart->CopyTo(pFieldX);
    }
    else
    {
        pFieldX->InitialField(EFIT_Zero);
    }

    //r = b - A x
    pFieldX->CopyTo(pW);
    pW->ApplyOperator(uiM, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields);
    pR->AxpyMinus(pW);

    const UBOOL bRet = IterateCG(pFieldX, pR, pFieldB, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, uiM);
    pW->Return();
    pR->Return();
    return bRet;
}

UBOOL CSLASolverCG::Solve(CField* pFieldX, const CField* pFieldB,
    INT gaugeNum, INT bosonNum, INT tensor2Num,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
    EFieldOperator uiM, ESolverPhase ePhase, const CField* pStart)
{
    _RECORD(CSLASolverCG::Solve);
    switch (uiM)
    {
    case EFO_F_DDdagger:
        //plain CG, D D^+ is Hermitian positive definite
        return SolveHPD(pFieldX, pFieldB, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, EFO_F_DDdagger, pStart);
    case EFO_F_D:
        {
            //x = D^{-1} b = D^+ (D D^+)^{-1} b
            const UBOOL bRet = SolveHPD(pFieldX, pFieldB, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, EFO_F_DDdagger, pStart);
            pFieldX->ApplyOperator(EFO_F_Ddagger, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields);
            return bRet;
        }
    case EFO_F_Ddagger:
        {
            //x = (D^+)^{-1} b = (D D^+)^{-1} D b
            CField* pW = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId, _T(__FILE__), __LINE__);
            pFieldB->CopyTo(pW);
            pW->ApplyOperator(EFO_F_D, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields);
            const UBOOL bRet = SolveHPD(pFieldX, pW, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, EFO_F_DDdagger, pStart);
            pW->Return();
            return bRet;
        }
    default:
        appCrucial(_T("CSolverCG does not support operator %s, it requires Hermitian positive definite solvers!\n"), __ENUM_TO_STRING(EFieldOperator, uiM).c_str());
        return FALSE;
    }
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
