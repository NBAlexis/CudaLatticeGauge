//=============================================================================
// FILENAME : CSolverBiCGStab.cpp
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [02/12/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CSLASolverGMRES)

CSLASolverGMRES::CSLASolverGMRES()
    : m_pOwner(NULL) 
    , m_uiReStart(3)
    , m_uiMaxDim(20)
    , m_fAccuracy(F(0.000001))
    , m_fBeta(F(0.0))
    , m_bAbsoluteAccuracy(FALSE)
{

}

CSLASolverGMRES::~CSLASolverGMRES()
{
    ReleaseBuffers();
}

void CSLASolverGMRES::Configurate(const CParameters& param)
{
    INT iValue;
    Real fValue;

    if (param.FetchValueINT(_T("MaxDim"), iValue))
    {
        m_uiMaxDim = static_cast<UINT>(iValue);
    }

    if (m_uiMaxDim < 5 || m_uiMaxDim > _kMaxStep)
    {
        appCrucial(_T("Max Dim must >= 5 and <= 100, set to default (20)"));
        m_uiMaxDim = 20;
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
    }
}

void CSLASolverGMRES::AllocateBuffers(const CField* )
{

}

void CSLASolverGMRES::ReleaseBuffers()
{

}

UBOOL CSLASolverGMRES::Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM, const CField* pStart)
{
    assert(0 == m_lstVectors.Num());
    for (UINT i = 0; i < m_uiMaxDim; ++i)
    {
        CField* pVectors = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
        m_lstVectors.AddItem(pVectors);
    }

    CField* pX = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pW = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pR = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);

    //use it to estimate relative error
    Real fBLength = F(1.0);
    if (!m_bAbsoluteAccuracy)
    {
        fBLength = pFieldB->Dot(pFieldB).x;
    }

    appParanoiac(_T("-- CSLASolverGMRES::Solve start operator: %s-- fLength = %f --\n"), __ENUM_TO_STRING(EFieldOperator, uiM).c_str(), fBLength);

    //set initial gauss x0 = b or pStart
    if (NULL == pStart)
    {
        pFieldB->CopyTo(pX);
    }
    else
    {
        pStart->CopyTo(pX);
    }

    Real fLastDiavation = F(0.0);
    for (UINT i = 0; i < m_uiReStart; ++i)
    {
        //r = b - A x0
        //v[0] = r.normalize
        //s = x0
        pX->CopyTo(pR); //x0 need to be preserved
        pR->ApplyOperator(uiM, pGaugeFeild, EOCT_Minus); //x0 = -A x0
        pR->AxpyPlus(pFieldB); //x0 = b-Ax0
        pR->CopyTo(m_lstVectors[0]);
        m_fBeta = _sqrt(m_lstVectors[0]->Dot(m_lstVectors[0]).x);
        m_lstVectors[0]->ScalarMultply(F(1.0) / m_fBeta);  //v[0] = (b - A x0).normalize
        for (UINT j = 0; j < m_uiMaxDim; ++j)
        {
            //w = A v[j]
            m_lstVectors[j]->CopyTo(pW);
            pW->ApplyOperator(uiM, pGaugeFeild);
            for (UINT k = 0; k <= j; ++k)
            {
                _Complex dotc = m_lstVectors[k]->Dot(pW);
                m_h[HIndex(k, j)] = dotc;
                //w -= h[k,j] v[k]
                pW->Axpy(_make_cuComplex(-dotc.x, -dotc.y), m_lstVectors[k]);
            }

            //h[j + 1, j] = ||w||
            Real fWNorm = _sqrt(pW->Dot(pW).x);
            m_h[HIndex(j + 1, j)] = _make_cuComplex(fWNorm, F(0.0));
            //v[j + 1] = w / ||w||
            if (j < m_uiMaxDim - 1)
            {
                pW->ScalarMultply(F(1.0) / fWNorm);
                pW->CopyTo(m_lstVectors[j + 1]);
            }
        }
        RotateH(m_uiMaxDim);
        fLastDiavation = __cuCabsSqf(m_g[m_uiMaxDim]);
        SolveY(m_uiMaxDim);
        for (UINT j = 0; j < m_uiMaxDim; ++j)
        {
            pX->Axpy(m_y[j], m_lstVectors[j]);
        }
        appParanoiac(_T("CSLASolverGMRES::Solve deviation: ---- restart ----. last divation = %f\n"), fLastDiavation);
        if (fLastDiavation < m_fAccuracy * fBLength)
        {
            pX->CopyTo(pFieldX);

            pX->Return();
            pW->Return();
            pR->Return();
            for (UINT i = 0; i < m_uiMaxDim; ++i)
            {
                m_lstVectors[i]->Return();
            }
            m_lstVectors.RemoveAll();
            return TRUE;
        }
    }

    appParanoiac(_T("CSLASolverGMRES::Solve failed: last divation = %f\n"), fLastDiavation);
    pX->CopyTo(pFieldX);

    pX->Return();
    pW->Return();
    pR->Return();
    for (UINT i = 0; i < m_uiMaxDim; ++i)
    {
        m_lstVectors[i]->Return();
    }
    m_lstVectors.RemoveAll();
    return FALSE;
}

void CSLASolverGMRES::RotateH(UINT uiHeisenbergDim)
{
    //======================= reset g ==================
    m_g[0] = _make_cuComplex(m_fBeta, F(0.0));
    for (UINT i = 0; i < uiHeisenbergDim; ++i)
    {
        UINT ii = HIndex(i, i);
        UINT i1i = HIndex(i + 1, i);
        Real denomi = F(1.0) / _sqrt(__cuCabsSqf(m_h[ii]) + __cuCabsSqf(m_h[i1i]));
        _Complex cs = cuCmulf_cr(m_h[ii], denomi);
        _Complex sn = cuCmulf_cr(m_h[i1i], denomi);
        _Complex cs_h = _cuConjf(cs);
        _Complex sn_h = _cuConjf(sn);

        for (UINT j = i; j < uiHeisenbergDim; ++j)
        {
            UINT ij = HIndex(i, j);
            UINT i1j = HIndex(i + 1, j);

            _Complex hij = m_h[ij];
            m_h[ij] = _cuCaddf(_cuCmulf(cs_h, hij), _cuCmulf(sn_h, m_h[i1j]));
            m_h[i1j] = _cuCsubf(_cuCmulf(cs, m_h[i1j]), _cuCmulf(sn, hij));
        }

        _Complex minus_gi = _make_cuComplex(-m_g[i].x, -m_g[i].y);
        m_g[i] = _cuCmulf(cs_h, m_g[i]);
        m_g[i + 1] = _cuCmulf(sn, minus_gi);
    }
}

void CSLASolverGMRES::SolveY(UINT uiHeisenbergDim)
{
    for (INT i = uiHeisenbergDim - 1; i > -1; --i)
    {
        for (INT j = i + 1; j < uiHeisenbergDim; ++j)
        {
            m_g[i] = _cuCsubf(m_g[i], _cuCmulf(m_h[HIndex(i, j)], m_y[j]));
        }
        m_y[i] = cuCdivf(m_g[i], m_h[HIndex(i,i)]);
    }
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================