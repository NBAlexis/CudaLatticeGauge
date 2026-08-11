//=============================================================================
// FILENAME : CSolverDeflatedCG.cpp
// 
// DESCRIPTION:
// Deflated conjugate gradient solver for Hermitian positive definite
// operators.
//
// REVISION:
//  [mm/dd/yy]
//  [07/24/2026 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CSolverDeflatedCG.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CSLASolverDeflatedCG)

/**
 * Solve A c = b in place with partial pivoting, A is dm x dm and destroyed
 */
static void _SolveLinearSystemHost(CLGComplex* A, CLGComplex* b, UINT dm)
{
    for (UINT col = 0; col < dm; ++col)
    {
        UINT pivot = col;
        Real fMax = _cuCabsf(A[pivot * dm + col]);
        for (UINT row = col + 1; row < dm; ++row)
        {
            const Real fAbs = _cuCabsf(A[row * dm + col]);
            if (fAbs > fMax)
            {
                fMax = fAbs;
                pivot = row;
            }
        }
        if (fMax < _CLG_FLT_MIN_)
        {
            appCrucial(_T("CSolverDeflatedCG: the deflation matrix is singular!\n"));
            return;
        }
        if (pivot != col)
        {
            for (UINT j = 0; j < dm; ++j)
            {
                const CLGComplex tmp = A[pivot * dm + j];
                A[pivot * dm + j] = A[col * dm + j];
                A[col * dm + j] = tmp;
            }
            const CLGComplex tmp = b[pivot];
            b[pivot] = b[col];
            b[col] = tmp;
        }
        for (UINT row = col + 1; row < dm; ++row)
        {
            const CLGComplex f = _cuCdivf(A[row * dm + col], A[col * dm + col]);
            for (UINT j = col; j < dm; ++j)
            {
                A[row * dm + j] = _cuCsubf(A[row * dm + j], _cuCmulf(f, A[col * dm + j]));
            }
            b[row] = _cuCsubf(b[row], _cuCmulf(f, b[col]));
        }
    }
    for (INT row = static_cast<INT>(dm) - 1; row >= 0; --row)
    {
        for (UINT j = row + 1; j < dm; ++j)
        {
            b[row] = _cuCsubf(b[row], _cuCmulf(A[row * dm + j], b[j]));
        }
        b[row] = _cuCdivf(b[row], A[row * dm + row]);
    }
}

CSLASolverDeflatedCG::CSLASolverDeflatedCG()
    : CSLASolverCG()
    , m_uiDeflateDim(16)
    , m_uiDeflate(8)
    , m_uiReDeflateInterval(0)
    , m_uiSolveCount(0)
    , m_bDeflated(FALSE)
    , m_pHelper(NULL)
    , m_pDeviceHm(NULL)
    , m_pDeviceEigenValue(NULL)
    , m_pDevicePk(NULL)
    , m_pHostHm(NULL)
    , m_pHostE(NULL)
    , m_pHostQ(NULL)
    , m_pHostR(NULL)
    , m_pHostY(NULL)
{

}

CSLASolverDeflatedCG::~CSLASolverDeflatedCG()
{
    for (INT i = 0; i < m_lstDeflateU.Num(); ++i)
    {
        appSafeDelete(m_lstDeflateU[i]);
    }
    if (NULL != m_pDeviceHm)
    {
        checkCudaErrors(__cudaFree(m_pDeviceHm));
        checkCudaErrors(__cudaFree(m_pDeviceEigenValue));
        checkCudaErrors(__cudaFree(m_pDevicePk));
    }
    appSafeFree(m_pHostHm);
    appSafeFree(m_pHostE);
    appSafeFree(m_pHostQ);
    appSafeFree(m_pHostR);
    appSafeFree(m_pHostY);
    appSafeDelete(m_pHelper);
}

void CSLASolverDeflatedCG::Configurate(const CParameters& param)
{
    CSLASolverCG::Configurate(param);

    INT iValue;
    if (param.FetchValueINT(_T("DeflateDim"), iValue))
    {
        m_uiDeflateDim = static_cast<UINT>(iValue);
    }
    if (param.FetchValueINT(_T("Deflate"), iValue))
    {
        m_uiDeflate = static_cast<UINT>(iValue);
    }
    if (m_uiDeflate < 1)
    {
        m_uiDeflate = 1;
    }
    if (m_uiDeflateDim < m_uiDeflate)
    {
        m_uiDeflateDim = m_uiDeflate;
    }
    if (param.FetchValueINT(_T("ReDeflateInterval"), iValue))
    {
        m_uiReDeflateInterval = static_cast<UINT>(iValue);
    }

    m_pHelper = new CLinearAlgebraHelper(m_uiDeflateDim + 1);
}

void CSLASolverDeflatedCG::AllocateBuffers(const CField*)
{
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceHm, sizeof(CLGComplex) * m_uiDeflateDim * m_uiDeflateDim));
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceEigenValue, sizeof(CLGComplex) * m_uiDeflate));
    checkCudaErrors(__cudaMalloc((void**)&m_pDevicePk, sizeof(CLGComplex) * m_uiDeflate * m_uiDeflateDim));

    m_pHostHm = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiDeflateDim * m_uiDeflateDim);
    m_pHostE = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiDeflate * m_uiDeflate);
    m_pHostQ = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiDeflate * m_uiDeflate);
    m_pHostR = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiDeflate * m_uiDeflate);
    m_pHostY = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiDeflate);
}

void CSLASolverDeflatedCG::BuildDeflationSpace(const CField* pFieldB,
    INT gaugeNum, INT bosonNum, INT tensor2Num,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
    EFieldOperator uiM)
{
    _RECORD(CSLASolverDeflatedCG::BuildDeflationSpace);
    TArray<CField*> lstV;
    for (UINT i = 0; i < m_uiDeflateDim; ++i)
    {
        lstV.AddItem(appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId, _T(__FILE__), __LINE__));
    }
    CField* pW = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId, _T(__FILE__), __LINE__);

    //Arnoldi on b, the lower triangle of the Hessenberg must be zero
    memset(m_pHostHm, 0, sizeof(CLGComplex) * m_uiDeflateDim * m_uiDeflateDim);
    pFieldB->CopyTo(lstV[0]);
    lstV[0]->ScalarMultply(F(1.0) / static_cast<Real>(_sqrt(lstV[0]->GetLength())));
    for (UINT j = 0; j < m_uiDeflateDim; ++j)
    {
        lstV[j]->CopyTo(pW);
        pW->ApplyOperator(uiM, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields);
        for (UINT i = 0; i <= j; ++i)
        {
            m_pHostHm[i * m_uiDeflateDim + j] = _cToRealC(lstV[i]->Dot(pW));
            pW->Axpy(_make_cuComplex(-m_pHostHm[i * m_uiDeflateDim + j].x, -m_pHostHm[i * m_uiDeflateDim + j].y), lstV[i]);
        }
        if (j + 1 < m_uiDeflateDim)
        {
            const Real fLength = static_cast<Real>(_sqrt(pW->GetLength()));
            m_pHostHm[(j + 1) * m_uiDeflateDim + j] = _make_cuComplex(fLength, F(0.0));
            pW->CopyTo(lstV[j + 1]);
            lstV[j + 1]->ScalarMultply(F(1.0) / fLength);
        }
    }

    checkCudaErrors(cudaMemcpy(m_pDeviceHm, m_pHostHm, sizeof(CLGComplex) * m_uiDeflateDim * m_uiDeflateDim, cudaMemcpyHostToDevice));
    m_pHelper->EigenValueProblemHessenberg(m_pDeviceHm, m_pDeviceEigenValue, m_pDevicePk, m_uiDeflateDim, m_uiDeflate);
    m_pHelper->Transpose(m_pDevicePk, m_uiDeflate, m_uiDeflateDim);

    //U = V Pk, the m_uiDeflate smallest Ritz vectors
    //Pk[j][i] is the coefficient of the j-th Arnoldi vector in the i-th Ritz vector
    CLGComplex* pHostPk = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiDeflateDim * m_uiDeflate);
    checkCudaErrors(cudaMemcpy(pHostPk, m_pDevicePk, sizeof(CLGComplex) * m_uiDeflateDim * m_uiDeflate, cudaMemcpyDeviceToHost));
    if (!m_bDeflated)
    {
        for (UINT i = 0; i < m_uiDeflate; ++i)
        {
            m_lstDeflateU.AddItem(pFieldB->GetCopy());
        }
    }
    for (UINT i = 0; i < m_uiDeflate; ++i)
    {
        m_lstDeflateU[i]->InitialField(EFIT_Zero);
        for (UINT j = 0; j < m_uiDeflateDim; ++j)
        {
            m_lstDeflateU[i]->Axpy(pHostPk[j * m_uiDeflate + i], lstV[j]);
        }
    }
    free(pHostPk);

    for (UINT i = 0; i < m_uiDeflateDim; ++i)
    {
        lstV[i]->Return();
    }
    pW->Return();
    m_bDeflated = TRUE;
}

void CSLASolverDeflatedCG::CoarseCorrection(CField* pFieldX, CField* pR,
    INT gaugeNum, INT bosonNum, INT tensor2Num,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
    EFieldOperator uiM)
{
    _RECORD(CSLASolverDeflatedCG::CoarseCorrection);
    TArray<CField*> lstW;
    for (UINT i = 0; i < m_uiDeflate; ++i)
    {
        lstW.AddItem(appGetLattice()->GetPooledFieldById(pR->m_byFieldId, _T(__FILE__), __LINE__));
    }

    //W_i = A U_i, E_ij = <U_i, W_j>, y_i = <U_i, r>
    for (UINT i = 0; i < m_uiDeflate; ++i)
    {
        m_lstDeflateU[i]->CopyTo(lstW[i]);
        lstW[i]->ApplyOperator(uiM, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields);
    }
    for (UINT i = 0; i < m_uiDeflate; ++i)
    {
        m_pHostY[i] = _cToRealC(m_lstDeflateU[i]->Dot(pR));
        for (UINT j = 0; j < m_uiDeflate; ++j)
        {
            m_pHostE[i * m_uiDeflate + j] = _cToRealC(m_lstDeflateU[i]->Dot(lstW[j]));
        }
    }

    //solve E c = y exactly on the deflation space with partial pivoting
    _SolveLinearSystemHost(m_pHostE, m_pHostY, m_uiDeflate);

    //x = x + sum_i c_i U_i, r = r - sum_i c_i (A U_i)
    for (UINT i = 0; i < m_uiDeflate; ++i)
    {
        pFieldX->Axpy(m_pHostY[i], m_lstDeflateU[i]);
        pR->Axpy(_make_cuComplex(-m_pHostY[i].x, -m_pHostY[i].y), lstW[i]);
        lstW[i]->Return();
    }
}

UBOOL CSLASolverDeflatedCG::SolveHPD(CField* pFieldX, const CField* pFieldB,
    INT gaugeNum, INT bosonNum, INT tensor2Num,
    const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields,
    EFieldOperator uiM, const CField* pStart)
{
    _RECORD(CSLASolverDeflatedCG::SolveHPD);
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

    if (!m_bDeflated || (0 != m_uiReDeflateInterval && 0 == (m_uiSolveCount % m_uiReDeflateInterval)))
    {
        BuildDeflationSpace(pR, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, uiM);
    }
    ++m_uiSolveCount;

    CoarseCorrection(pFieldX, pR, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, uiM);

    const UBOOL bRet = IterateCG(pFieldX, pR, pFieldB, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, uiM);
    pW->Return();
    pR->Return();
    return bRet;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
