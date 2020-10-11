//=============================================================================
// FILENAME : CMultiShiftGMRES.cpp
// 
// DESCRIPTION:
// This is the class for GMRES Solver
//
// REVISION:
//  [15/06/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMultiShiftGMRES)

CMultiShiftGMRES::CMultiShiftGMRES()
    : CMultiShiftSolver()
    , m_uiReStart(3)
    , m_uiMaxDim(20)
    , m_fAccuracy(F(0.000001))
    , m_fBeta(F(0.0))
    , m_bUseCudaForSmallMatrix(FALSE)
    , m_bCheckAddSystem(FALSE)
    , m_pHelper(NULL)
    , m_pDeviceY(NULL)
    , m_pDeviceZ(NULL)
    , m_pDeviceHHat(NULL)
    , m_pDeviceQ(NULL)
    , m_pDeviceR(NULL)
{

}

CMultiShiftGMRES::~CMultiShiftGMRES()
{
    CMultiShiftGMRES::ReleaseBuffers();
    appSafeDelete(m_pHelper);
}

void CMultiShiftGMRES::Configurate(const CParameters& param)
{
    INT iValue;
    Real fValue;

    if (param.FetchValueINT(_T("MaxDim"), iValue))
    {
        m_uiMaxDim = static_cast<UINT>(iValue);
    }

    if (m_uiMaxDim < 5 || m_uiMaxDim > CLinearAlgebraHelper::_kMaxSmallDim - 2)
    {
        appCrucial(_T("Max Dim must >= 5 and <= 30, set to default (20)"));
        m_uiMaxDim = 20;
    }

    m_pHelper = new CLinearAlgebraHelper(m_uiMaxDim + 1, 1);

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceHHat, sizeof(CLGComplex) * (m_uiMaxDim + 1) * (m_uiMaxDim + 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceQ, sizeof(CLGComplex) * (m_uiMaxDim + 1) * (m_uiMaxDim + 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceR, sizeof(CLGComplex) * (m_uiMaxDim + 1) * (m_uiMaxDim + 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceY, sizeof(CLGComplex) * (m_uiMaxDim + 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceZ, sizeof(CLGComplex) * (m_uiMaxDim + 1)));

    iValue = 0;
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
    iValue = 0;
    if (param.FetchValueINT(_T("CheckAddSystem"), iValue))
    {
        m_bCheckAddSystem = (0 != iValue);
    }

    for (UINT i = 0; i < _kMaxStep * (_kMaxStep + 1); ++i)
    {
        m_h[i] = _zeroc;
    }
}

void CMultiShiftGMRES::AllocateBuffers(const CField*)
{

}

void CMultiShiftGMRES::ReleaseBuffers()
{
    if (NULL != m_pDeviceHHat)
    {
        checkCudaErrors(cudaFree(m_pDeviceHHat));
        checkCudaErrors(cudaFree(m_pDeviceQ));
        checkCudaErrors(cudaFree(m_pDeviceR));
        checkCudaErrors(cudaFree(m_pDeviceY));
        checkCudaErrors(cudaFree(m_pDeviceZ));
    }
}

UBOOL CMultiShiftGMRES::Solve(TArray<CField*>& pFieldX, const TArray<CLGComplex>& cn, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM, ESolverPhase ePhase, const CField* pStart)
{
    appSetLogDate(FALSE);
    assert(0 == m_lstVectors.Num());
    for (UINT i = 0; i < m_uiMaxDim; ++i)
    {
        CField* pVectors = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
        m_lstVectors.AddItem(pVectors);
    }

    CField* pX = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pW = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pR = m_lstVectors[0];

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
    for (INT i = 0; i < cn.Num(); ++i)
    {
        beta0.AddItem(_onec);
    }

    appParanoiac(_T("-- CMultiShiftGMRES::Solve start operator: %s-- fLength = %f --\n"), __ENUM_TO_STRING(EFieldOperator, uiM).c_str(), fBLength);

    pX->InitialField(EFIT_Zero);
    for (INT i = 0; i < cn.Num(); ++i)
    {
        pFieldX[i]->InitialField(EFIT_Zero);
    }

    Real fMaxError = F(0.0);
    for (UINT i = 0; i < m_uiReStart; ++i)
    {
        //===============================================
        // Seed system
        //r = b - A x0
        //v[0] = r.normalize
        //s = x0
        if (0 == i)
        {
            pFieldB->CopyTo(pR);
        }
        else
        {
            pX->CopyTo(pR); //x0 need to be preserved
            pR->ApplyOperator(uiM, pGaugeFeild, EOCT_Minus); //x0 = -A x0
            pR->AxpyPlus(pFieldB); //x0 = b-Ax0
        }
#if !_CLG_DOUBLEFLOAT
        m_fBeta = static_cast<Real>(_sqrtd(m_lstVectors[0]->Dot(m_lstVectors[0]).x));
#else
        m_fBeta = _sqrt(m_lstVectors[0]->Dot(m_lstVectors[0]).x);
#endif
        m_lstVectors[0]->ScalarMultply(F(1.0) / m_fBeta);
        
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
        }
        checkCudaErrors(cudaMemcpy(m_pDeviceHHat, m_h, sizeof(CLGComplex) * (m_uiMaxDim + 1) * m_uiMaxDim, cudaMemcpyHostToDevice));
        m_pHelper->InitialZeroHost(m_g, m_uiMaxDim + 1, 1);
        m_g[0] = _make_cuComplex(m_fBeta, F(0.0));
        checkCudaErrors(cudaMemcpy(m_pDeviceY, m_g, sizeof(CLGComplex) * (m_uiMaxDim + 1), cudaMemcpyHostToDevice));
        m_pHelper->RotateHenssenberg(m_pDeviceHHat, m_pDeviceY, m_uiMaxDim);
        m_pHelper->SolveY(m_pDeviceY, m_pDeviceHHat, 1, m_uiMaxDim);
        checkCudaErrors(cudaMemcpy(m_y, m_pDeviceY, sizeof(CLGComplex) * m_uiMaxDim, cudaMemcpyDeviceToHost));
        for (UINT j = 0; j < m_uiMaxDim; ++j)
        {
            pX->Axpy(m_y[j], m_lstVectors[j]);
        }

        //===============================================
        // Add system
        // zm+1 = H * y
        checkCudaErrors(cudaMemcpy(m_pDeviceHHat, m_h, sizeof(CLGComplex) * (m_uiMaxDim + 1) * m_uiMaxDim, cudaMemcpyHostToDevice));
        //m_pHelper->PrintDeviceMatrix(m_pDeviceHHat, m_uiMaxDim + 1, m_uiMaxDim);
        //m_pHelper->PrintDeviceMatrix(m_pDeviceY, m_uiMaxDim, 1);
        m_pHelper->SmallMatrixMult(m_pDeviceZ, m_pDeviceHHat, m_pDeviceY, m_uiMaxDim + 1, m_uiMaxDim, 1, FALSE, FALSE);
        checkCudaErrors(cudaMemcpy(m_z, m_pDeviceZ, sizeof(CLGComplex) * (m_uiMaxDim + 1), cudaMemcpyDeviceToHost));
        for (UINT j = 0; j < m_uiMaxDim + 1; ++j)
        {
            if (0 == j)
            {
                m_z[j] = _cuCsubf(_make_cuComplex(m_fBeta, F(0.0)), m_z[j]);
            }
            else
            {
                m_z[j] = _make_cuComplex(-m_z[j].x, -m_z[j].y);
            }
        }

        //The unshifted system is always the last to converge, use m_fBeta may overestimate
        //fMaxError = m_fBeta;
        fMaxError = F(0.0);
        for (INT n = 0; n < pFieldX.Num(); ++n)
        {

            //=============
            //1. For each cn, calculate HHat
            for (UINT x = 0; x < m_uiMaxDim + 1; ++x)
            {
                for (UINT y = 0; y < m_uiMaxDim + 1; ++y)
                {
                    const UINT idxH = x + y * m_uiMaxDim;
                    const UINT idxHHat = x + y * (m_uiMaxDim + 1);
                    if (x == m_uiMaxDim)
                    {
                        m_hCopy[idxHHat] = m_z[y];
                    }
                    else
                    {
                        m_hCopy[idxHHat] = m_h[idxH];
                        if (x == y)
                        {
                            m_hCopy[idxHHat] = _cuCaddf(m_hCopy[idxHHat], cn[n]);
                        }
                    }
                }
            }

            //=============
            //2. QR factorization of HHat
            checkCudaErrors(cudaMemcpy(m_pDeviceHHat, m_hCopy, sizeof(CLGComplex) * (m_uiMaxDim + 1) * (m_uiMaxDim + 1), cudaMemcpyHostToDevice));
            //m_pHelper->PrintDeviceMatrix(m_pDeviceHHat, m_uiMaxDim + 1, m_uiMaxDim + 1);
            m_pHelper->QRFactorization(m_pDeviceQ, m_pDeviceR, m_pDeviceHHat, m_uiMaxDim + 1);

            //=============
            //3. Solve y hat and beta m
            m_pHelper->InitialZeroHost(m_g, m_uiMaxDim + 1, 1);
            m_g[0] = cuCmulf_cr(beta0[n], m_fBeta);

            checkCudaErrors(cudaMemcpy(m_pDeviceZ, m_g, sizeof(CLGComplex) * (m_uiMaxDim + 1), cudaMemcpyHostToDevice));
            m_pHelper->Dagger(m_pDeviceQ, m_uiMaxDim + 1, m_uiMaxDim + 1);
            m_pHelper->SmallMatrixMult(m_pDeviceY, m_pDeviceQ, m_pDeviceZ, m_uiMaxDim + 1, m_uiMaxDim + 1, 1, FALSE, FALSE);
            m_pHelper->SolveY(m_pDeviceY, m_pDeviceR, 1, m_uiMaxDim + 1);

            //Now the first m of m_pDeviceY is YHat, and the last of m_pDeviceY is new beta0
            checkCudaErrors(cudaMemcpy(m_yhat, m_pDeviceY, sizeof(CLGComplex) * (m_uiMaxDim + 1), cudaMemcpyDeviceToHost));

            beta0[n] = m_yhat[m_uiMaxDim];
            const Real fErrorN = _cuCabsf(beta0[n]);
            if (fErrorN > fMaxError)
            {
                fMaxError = fErrorN;
            }
            appParanoiac(_T("CMultiShiftGMRES::Solve deviation: ----  last beta0 %d = %8.15f\n"), n, fErrorN);

            for (UINT j = 0; j < m_uiMaxDim; ++j)
            {
                pFieldX[n]->Axpy(m_yhat[j], m_lstVectors[j]);
            }
        }

        //Stop if beta is very small
        if (fMaxError < m_fAccuracy * fBLength)
        {
            appParanoiac(_T("CMultiShiftGMRES::Solve deviation: ---- finished ----. last beta = %8.15f\n"), m_fBeta);
            pX->Return();
            pW->Return();
            for (UINT k = 0; k < m_uiMaxDim; ++k)
            {
                m_lstVectors[k]->Return();
            }
            m_lstVectors.RemoveAll();
            return TRUE;
        }

        appParanoiac(_T("CMultiShiftGMRES::Solve deviation: ----k:%d  last beta = %8.15f\n"), i, fMaxError);
    }

    appGeneral(_T("CSLASolverGMRES::Solve failed: last divation = %8.15f\n"), fMaxError);
    pX->Return();
    pW->Return();
    for (UINT i = 0; i < m_uiMaxDim; ++i)
    {
        m_lstVectors[i]->Return();
    }
    m_lstVectors.RemoveAll();
    return FALSE;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================