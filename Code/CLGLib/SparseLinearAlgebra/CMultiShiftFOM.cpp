//=============================================================================
// FILENAME : CMultiShiftFOM.cpp
// 
// DESCRIPTION:
// This is the class for GMRES Solver
//
// REVISION:
//  [18/06/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMultiShiftFOM)

CMultiShiftFOM::CMultiShiftFOM()
    : CMultiShiftSolver()
    , m_uiReStart(3)
    , m_uiMaxDim(20)
    , m_fAccuracy(F(0.000001))
    , m_bUseCudaForSmallMatrix(FALSE)
    , m_pHelper(NULL)
    , m_pDeviceY(NULL)
    , m_pDeviceHHat(NULL)
{

}

CMultiShiftFOM::~CMultiShiftFOM()
{
    CMultiShiftFOM::ReleaseBuffers();
    appSafeDelete(m_pHelper);
}

void CMultiShiftFOM::Configurate(const CParameters& param)
{
    INT iValue;
    Real fValue;

    if (param.FetchValueINT(_T("MaxDim"), iValue))
    {
        m_uiMaxDim = static_cast<UINT>(iValue);
    }

    if (m_uiMaxDim < 5 || m_uiMaxDim > _kMaxStep)
    {
        appCrucial(_T("Max Dim must >= 5 and <= 30, set to default (20)"));
        m_uiMaxDim = 20;
    }

    iValue = 0;
    if (param.FetchValueINT(_T("UseCudaForSmallMatrix"), iValue))
    {
        m_bUseCudaForSmallMatrix = (0 != iValue);
    }
    if (m_bUseCudaForSmallMatrix && m_uiMaxDim < CLinearAlgebraHelper::_kMaxSmallDim)
    {
        m_pHelper = new CLinearAlgebraHelper(m_uiMaxDim + 1, 1);
        checkCudaErrors(cudaMalloc((void**)&m_pDeviceHHat, sizeof(CLGComplex) * m_uiMaxDim * m_uiMaxDim));
        checkCudaErrors(cudaMalloc((void**)&m_pDeviceY, sizeof(CLGComplex) * m_uiMaxDim));
    }

    if (param.FetchValueINT(_T("Restart"), iValue))
    {
        m_uiReStart = static_cast<UINT>(iValue);
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

void CMultiShiftFOM::AllocateBuffers(const CField*)
{

}

void CMultiShiftFOM::ReleaseBuffers()
{
    if (NULL != m_pHelper)
    {
        checkCudaErrors(cudaFree(m_pDeviceY));
        checkCudaErrors(cudaFree(m_pDeviceHHat));
    }
}

UBOOL CMultiShiftFOM::Solve(TArray<CField*>& pFieldX, const TArray<CLGComplex>& cn, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM, ESolverPhase ePhase, const CField* pStart)
{
    appSetLogDate(FALSE);
    assert(0 == m_lstVectors.Num());
    for (UINT i = 0; i < m_uiMaxDim; ++i)
    {
        CField* pVectors = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
        m_lstVectors.AddItem(pVectors);
    }

    CField* pW = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pR = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);

    //use it to estimate relative error
#if !_CLG_DOUBLEFLOAT
    DOUBLE fBLength = 1.0;
#else
    Real fBLength = F(1.0);
#endif
    if (!m_bAbsoluteAccuracy)
    {
        fBLength = pFieldB->GetLength();//pFieldB->Dot(pFieldB).x;
    }
    TArray<CLGComplex> beta0;
    TArray<Real> beta0Length;
    for (INT i = 0; i < cn.Num(); ++i)
    {
        beta0.AddItem(_onec);
        beta0Length.AddItem(F(1.0));
    }

    appParanoiac(_T("-- CMultiShiftFOM::Solve start operator: %s-- fLength = %f --\n"), __ENUM_TO_STRING(EFieldOperator, uiM).c_str(), fBLength);

    //set initial gauss x0 = b or pStart
    for (INT i = 0; i < cn.Num(); ++i)
    {
        pFieldX[i]->InitialField(EFIT_Zero);
    }

    Real fMaxError = F(0.0);
    for (UINT i = 0; i < m_uiReStart; ++i)
    {
        //===============================================
        // Seed system

        if (0 == i)
        {
            pFieldB->CopyTo(pR);
        }
        pR->CopyTo(m_lstVectors[0]);
        if (0 == i)
        {
#if !_CLG_DOUBLEFLOAT
            const Real fBeta = static_cast<Real>(_sqrtd(m_lstVectors[0]->Dot(m_lstVectors[0]).x));
#else
            const Real fBeta = _sqrt(m_lstVectors[0]->Dot(m_lstVectors[0]).x);
#endif
            for (INT n = 0; n < cn.Num(); ++n)
            {
                beta0[n] = _make_cuComplex(fBeta, F(0.0));
                beta0Length[n] = fBeta;
            }
            m_lstVectors[0]->ScalarMultply(F(1.0) / fBeta);  //v[0] = (b - A x0).normalize 
        }

        for (UINT j = 0; j < m_uiMaxDim; ++j)
        {
            //w = A v[j]
            m_lstVectors[j]->CopyTo(pW);
            pW->ApplyOperator(uiM, pGaugeFeild);
            for (UINT k = 0; k <= j; ++k)
            {
#if !_CLG_DOUBLEFLOAT
                const CLGComplex dotc = _cToFloat(m_lstVectors[k]->Dot(pW));
#else
                const CLGComplex dotc = m_lstVectors[k]->Dot(pW);
#endif
                m_h[HIndex(k, j)] = dotc;
                //w -= h[k,j] v[k]
                pW->Axpy(_make_cuComplex(-dotc.x, -dotc.y), m_lstVectors[k]);
            }

            //h[j + 1, j] = ||w||
#if !_CLG_DOUBLEFLOAT
            const Real fWNorm = static_cast<Real>(_sqrtd(pW->Dot(pW).x));
#else
            const Real fWNorm = _sqrt(pW->Dot(pW).x);
#endif
            m_h[HIndex(j + 1, j)] = _make_cuComplex(fWNorm, F(0.0));
            //v[j + 1] = w / ||w||
            if (j < m_uiMaxDim - 1)
            {
                pW->ScalarMultply(F(1.0) / fWNorm);
                pW->CopyTo(m_lstVectors[j + 1]);
            }
            else
            {
                pW->ScalarMultply(F(1.0) / fWNorm);
                pW->CopyTo(pR);
            }
        }

        //===============================================
        // Add system
        fMaxError = F(0.0);
        for (INT n = 0; n < pFieldX.Num(); ++n)
        {
            if (beta0Length[n] < m_fAccuracy * fBLength)
            {
                continue;
            }

            //=============
            //1. For each cn, calculate HHat
            memcpy(m_hCopy, m_h, sizeof(CLGComplex) * m_uiMaxDim * m_uiMaxDim);
            for (UINT x = 0; x < m_uiMaxDim; ++x)
            {
                const UINT idxH = x + x * m_uiMaxDim;
                m_hCopy[idxH] = _cuCaddf(m_hCopy[idxH], cn[n]);
            }
            if (NULL != m_pHelper)
            {
                m_pHelper->InitialZeroHost(m_g, m_uiMaxDim, 1);
            }
            //If do not use helper, the rest of g is irrelevant, only initial the first element
            m_g[0] = beta0[n];

            RotateHSolveY(m_hCopy, m_g, m_y);

            beta0[n] = cuCmulf_cr(m_y[m_uiMaxDim - 1], -m_h[m_uiMaxDim - 1 + m_uiMaxDim * m_uiMaxDim].x);
            const Real fErrorN = _cuCabsf(beta0[n]);
            beta0Length[n] = fErrorN;
            if (fErrorN > fMaxError)
            {
                fMaxError = fErrorN;
            }
            appParanoiac(_T("CMultiShiftFOM::Solve deviation: ----  last beta0 %d = %8.15f\n"), n, fErrorN);
            for (UINT j = 0; j < m_uiMaxDim; ++j)
            {
                pFieldX[n]->Axpy(m_y[j], m_lstVectors[j]);
            }
        }

        //Stop if beta is very small
        if (fMaxError < m_fAccuracy * fBLength)
        {
            appParanoiac(_T("CMultiShiftFOM::Solve deviation: ---- finished ----. last beta = %8.15f\n"), fMaxError);
            pR->Return();
            pW->Return();
            for (UINT k = 0; k < m_uiMaxDim; ++k)
            {
                m_lstVectors[k]->Return();
            }
            m_lstVectors.RemoveAll();
            return TRUE;
        }

        appParanoiac(_T("CMultiShiftFOM::Solve deviation: ----  last beta = %8.15f\n"), fMaxError);
    }

    appGeneral(_T("CMultiShiftFOM::Solve failed: last divation = %8.15f\n"), fMaxError);
    pR->Return();
    pW->Return();
    for (UINT i = 0; i < m_uiMaxDim; ++i)
    {
        m_lstVectors[i]->Return();
    }
    m_lstVectors.RemoveAll();
    return FALSE;
}

void CMultiShiftFOM::RotateH(CLGComplex* h, CLGComplex* g, UINT uiDimX, UINT uiDimY)
{
    for (UINT i = 0; i < uiDimY - 1; ++i)
    {
        const UINT ii = i + i * uiDimX;
        const UINT i1i = i + (i + 1) * uiDimX;
        const Real denomi = F(1.0) / _sqrt(__cuCabsSqf(h[ii]) + __cuCabsSqf(h[i1i]));
        const CLGComplex cs = cuCmulf_cr(h[ii], denomi);
        const CLGComplex sn = cuCmulf_cr(h[i1i], denomi);
        const CLGComplex cs_h = _cuConjf(cs);
        const CLGComplex sn_h = _cuConjf(sn);

        for (UINT j = i; j < uiDimX; ++j)
        {
            const UINT ij = j + i * uiDimX;
            const UINT i1j = j + (i + 1) * uiDimX;

            const CLGComplex hij = h[ij];
            h[ij] = _cuCaddf(_cuCmulf(cs_h, hij), _cuCmulf(sn_h, h[i1j]));
            h[i1j] = _cuCsubf(_cuCmulf(cs, h[i1j]), _cuCmulf(sn, hij));
        }

        const CLGComplex minus_gi = _make_cuComplex(-g[i].x, -g[i].y);
        g[i] = _cuCmulf(cs_h, g[i]);
        g[i + 1] = _cuCmulf(sn, minus_gi);

        //m_g[i] = _cuCaddf(_cuCmulf(cs_h, m_g[i]), _cuCmulf(sn_h, m_g[i + 1]));
        //m_g[i + 1] = _cuCaddf(_cuCmulf(sn, minus_gi), _cuCmulf(cs, m_g[i + 1]));
    }
}

void CMultiShiftFOM::SolveY(CLGComplex* h, CLGComplex* g, CLGComplex* y, UINT uiDimY, UINT uiDimX)
{
    const INT iHeisenbergDim = static_cast<INT>(uiDimX);
    for (INT i = uiDimY - 1; i > -1; --i)
    {
        for (INT j = i + 1; j < iHeisenbergDim; ++j)
        {
            g[i] = _cuCsubf(g[i], _cuCmulf(h[j + i * uiDimX], y[j]));
        }
        y[i] = _cuCdivf(g[i], h[i + i * uiDimX]);
    }
}

void CMultiShiftFOM::RotateHSolveY(CLGComplex* h, CLGComplex* g, CLGComplex* y) const
{
    if (NULL == m_pHelper)
    {
        RotateH(h, g, m_uiMaxDim, m_uiMaxDim);
        SolveY(h, g, y, m_uiMaxDim, m_uiMaxDim);
        return;
    }

    checkCudaErrors(cudaMemcpy(m_pDeviceHHat, h, sizeof(CLGComplex) * m_uiMaxDim * m_uiMaxDim, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(m_pDeviceY, g, sizeof(CLGComplex) * m_uiMaxDim, cudaMemcpyHostToDevice));
    m_pHelper->RotateHenssenberg(m_pDeviceHHat, m_pDeviceY, m_uiMaxDim, m_uiMaxDim);
    m_pHelper->SolveY(m_pDeviceY, m_pDeviceHHat, 1, m_uiMaxDim);
    checkCudaErrors(cudaMemcpy(y, m_pDeviceY, sizeof(CLGComplex) * m_uiMaxDim, cudaMemcpyDeviceToHost));
}


__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================