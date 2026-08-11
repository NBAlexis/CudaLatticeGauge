//=============================================================================
// FILENAME : CMultiShiftCG.cpp
// 
// DESCRIPTION:
// Multi-shift conjugate gradient solver for (A + sigma_n) x_n = b,
// where A must be Hermitian positive definite (EFO_F_DDdagger).
// All shifted systems share the same Krylov space, the shifted iterates
// are obtained from the seed (unshifted) CG recurrence by the zeta
// recurrences, see B. Jegerlehner, hep-lat/9608029.
//
// REVISION:
//  [mm/dd/yy]
//  [07/18/2026 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CMultiShiftCG.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMultiShiftCG)

CMultiShiftCG::CMultiShiftCG()
    : CMultiShiftSolver()
    , m_uiMaxStep(5000)
{

}

void CMultiShiftCG::Configurate(const CParameters& param)
{
    CMultiShiftSolver::Configurate(param);

    INT iValue = 5000;

    if (param.FetchValueINT(_T("MaxStep"), iValue))
    {
        m_uiMaxStep = static_cast<UINT>(iValue);
    }

}

void CMultiShiftCG::AllocateBuffers(const CField*)
{

}

UBOOL CMultiShiftCG::Solve(TArray<CField*>& pFieldX, const TArray<CLGComplex>& cn, const CField* pFieldB,
    INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
    EFieldOperator uiM, ESolverPhase ePhase, const CField* pStart)
{
    //Multi-shift CG only works with Hermitian positive definite operators
    if (EFO_F_DDdagger != uiM)
    {
        appCrucial(_T("CMultiShiftCG only supports EFO_F_DDdagger, it requires a Hermitian positive definite operator!\n"));
        return FALSE;
    }

    appPushLogDate(FALSE);
    _RECORD(CMultiShiftCG::Solve);
    appAssert(0 == m_lstVectors.Num());
    const INT iShiftCount = cn.Num();

    //per-shift quantities, the seed system is the unshifted one (sigma = 0)
    TArray<DOUBLE> sigma, zeta, zetaPrev;
    TArray<UBOOL> converged;
    for (INT i = 0; i < iShiftCount; ++i)
    {
        //CG requires real shifts, the rational denominators are real
        if (appAbs(cn[i].y) > _CLG_FLT_EPSILON)
        {
            appCrucial(_T("CMultiShiftCG only supports real shifts, shift %d is (%f, %f)\n"), i, cn[i].x, cn[i].y);
            appPopLogDate();
            return FALSE;
        }
        if (cn[i].x < F(0.0))
        {
            appCrucial(_T("CMultiShiftCG requires non-negative shifts, shift %d is %f\n"), i, cn[i].x);
            appPopLogDate();
            return FALSE;
        }
    }
    for (INT i = 0; i < iShiftCount; ++i)
    {
        sigma.AddItem(static_cast<DOUBLE>(cn[i].x));
        zeta.AddItem(1.0);
        zetaPrev.AddItem(1.0);
        converged.AddItem(FALSE);
        m_lstVectors.AddItem(appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId, _T(__FILE__), __LINE__));
    }

    CField* pR = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId, _T(__FILE__), __LINE__);
    CField* pP = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId, _T(__FILE__), __LINE__);
    CField* pW = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId, _T(__FILE__), __LINE__);

    //x_n = 0, r = p = b, and p_n = b
    pFieldB->CopyTo(pR);
    pFieldB->CopyTo(pP);
    for (INT i = 0; i < iShiftCount; ++i)
    {
        pFieldX[i]->InitialField(EFIT_Zero);
        pFieldB->CopyTo(m_lstVectors[i]);
    }

    //use it to estimate relative error
    DOUBLE fBLength = 1.0;
    if (!m_bAbsoluteAccuracy)
    {
        fBLength = sqrt(pFieldB->GetLength());
    }
    const DOUBLE fTarget = static_cast<DOUBLE>(m_fAccuracy) * fBLength;
    appParanoiac(_T("-- CMultiShiftCG::Solve start operator: %s-- fLength = %f --\n"), __ENUM_TO_STRING(EFieldOperator, uiM).c_str(), fBLength);

    //gamma = (r, r)
    DOUBLE fRNorm2 = pR->GetLength();
    DOUBLE fAlphaPrev = 1.0;
    DOUBLE fBetaPrev = 0.0;

    UINT uiIter = 0;
    INT iNotConverged = iShiftCount;
    //check whether already converged, for example b = 0
    for (INT i = 0; i < iShiftCount; ++i)
    {
        const DOUBLE fError = appAbs(zeta[i]) * sqrt(fRNorm2);
        if (fError < fTarget)
        {
            converged[i] = TRUE;
            --iNotConverged;
        }
    }

    while (iNotConverged > 0 && uiIter < m_uiMaxStep)
    {
        //w = A p, one operator application per iteration, shared by all shifts
        pP->CopyTo(pW);
        pW->ApplyOperator(uiM, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields);
        const DOUBLE fPAQ = pP->Dot(pW).x;
        if (fPAQ <= _CLG_FLT_MIN_)
        {
            appGeneral(_T("CMultiShiftCG breakdown: (p, Ap) = %f\n"), fPAQ);
            break;
        }
        const DOUBLE fAlpha = fRNorm2 / fPAQ;

        //seed: r = r - alpha * w
        pR->Axpy(static_cast<Real>(-fAlpha), pW);
        const DOUBLE fRNorm2New = pR->GetLength();
        const DOUBLE fBeta = fRNorm2New / fRNorm2;
        const DOUBLE fRNormNew = sqrt(fRNorm2New);

        //shifted recurrences, see Jegerlehner hep-lat/9608029
        for (INT i = 0; i < iShiftCount; ++i)
        {
            if (converged[i])
            {
                continue;
            }
            const DOUBLE fZetaNew = zeta[i] / (1.0 + fAlpha * sigma[i]
                + (fAlpha * fBetaPrev / fAlphaPrev) * (1.0 - zeta[i] / zetaPrev[i]));
            const DOUBLE fAlphaShift = fAlpha * fZetaNew / zeta[i];
            const DOUBLE fBetaShift = fBeta * fZetaNew * fZetaNew / (zeta[i] * zeta[i]);

            //x_n = x_n + alpha_n p_n
            pFieldX[i]->Axpy(static_cast<Real>(fAlphaShift), m_lstVectors[i]);
            //p_n = zeta_new * r + beta_n p_n
            m_lstVectors[i]->ScalarMultply(static_cast<Real>(fBetaShift));
            m_lstVectors[i]->Axpy(static_cast<Real>(fZetaNew), pR);

            zetaPrev[i] = zeta[i];
            zeta[i] = fZetaNew;

            //|r_n| = |zeta_n| * |r|
            const DOUBLE fError = appAbs(fZetaNew) * fRNormNew;
            appParanoiac(_T("CMultiShiftCG::Solve deviation: ---- iter %d, shift %d = %8.15f\n"), uiIter, i, fError);
            if (fError < fTarget)
            {
                converged[i] = TRUE;
                --iNotConverged;
            }
        }

        //seed: p = r + beta * p
        pP->ScalarMultply(static_cast<Real>(fBeta));
        pP->AxpyPlus(pR);

        fAlphaPrev = fAlpha;
        fBetaPrev = fBeta;
        fRNorm2 = fRNorm2New;
        ++uiIter;
    }

    for (INT i = 0; i < iShiftCount; ++i)
    {
        appParanoiac(_T("CMultiShiftCG::Solve deviation: ---- final, shift %d zeta = %8.15f\n"), i, zeta[i]);
        m_lstVectors[i]->Return();
    }
    m_lstVectors.RemoveAll();
    if (iNotConverged > 0)
    {
        appGeneral(_T("CMultiShiftCG::Solve failed to converge, %d shifts not converged, field ID = %d\n"), iNotConverged, pFieldX[0]->m_byFieldId);
    }
    else
    {
        appParanoiac(_T("CMultiShiftCG::Solve finished, %d iterations, %d shifts not converged\n"), uiIter, iNotConverged);
    }
    
    pR->Return();
    pP->Return();
    pW->Return();
    appPopLogDate();
    return (0 == iNotConverged);
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
