//=============================================================================
// FILENAME : CSLASolverGCRODR.cpp
// 
// DESCRIPTION:
// This is the class for GMRES Solver
//
// REVISION:
//  [03/15/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CSLASolverGCRODR)

CSLASolverGCRODR::CSLASolverGCRODR()
    : CSLASolver()
    , m_pHelper(NULL)
    , m_uiReStart(3)
    , m_uiMDim(10)
    , m_uiKDim(5)
    , m_fAccuracy(F(0.000001))
    , m_bHasYk(FALSE)
    , m_fBeta(F(0.0))
    , m_fDiviation(F(0.0))
    , m_eDeflationType(EEDT_SVD)

    , m_pDeviceHm(NULL)
    , m_pDeviceEigenValue(NULL)
    , m_pDevicePk(NULL)
    , m_pDeviceHmGm(NULL)

    , m_pDeviceALeft(NULL)
    , m_pDeviceA(NULL)
    , m_pDeviceB(NULL)

    , m_pHostHmGm(NULL)

    , m_pHostY(NULL)

    , m_pHostALeft(NULL)
    , m_pHostB(NULL)
    , m_pHostPk(NULL)
    , m_pHostTmpQ(NULL)
    , m_pHostTmpR(NULL)
    , m_pHostTmpGPk(NULL)

{

}

CSLASolverGCRODR::~CSLASolverGCRODR()
{
    ReleaseBuffers();
    appSafeDelete(m_pHelper);
}

void CSLASolverGCRODR::Configurate(const CParameters& param)
{
    INT iValue;
    Real fValue;

    if (param.FetchValueINT(_T("MDim"), iValue))
    {
        m_uiMDim = static_cast<UINT>(iValue);
    }

    if (m_uiMDim < 5 || m_uiMDim > _kMaxDimDR)
    {
        appCrucial(_T("M Dim must >= 5 and <= 31, set to default (10)"));
        m_uiMDim = 10;
    }

    if (param.FetchValueINT(_T("KDim"), iValue))
    {
        m_uiKDim = static_cast<UINT>(iValue);
    }

    if (m_uiKDim < _kMinKDim || m_uiKDim > m_uiMDim - 2)
    {
        appCrucial(_T("K Dim must >= 2 and <= MDim - 2, set to default (MDim - 2)"));
        m_uiKDim = m_uiMDim - 2;
    }

    m_pHelper = new CLinearAlgebraHelper(m_uiMDim + 1);
    
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

    CCString sDeflate = _T("EEDT_SVD");
    param.FetchStringValue(_T("DeflateSpace"), sDeflate);
    m_eDeflationType = __STRING_TO_ENUM(EEigenDeflationType, sDeflate);
}

/**
* Will get called in CreateFermionSolver
*/
void CSLASolverGCRODR::AllocateBuffers(const CField* pFieldB)
{
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceHm, sizeof(CLGComplex) * m_uiMDim * m_uiMDim));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceEigenValue, sizeof(CLGComplex) * m_uiKDim));
    checkCudaErrors(cudaMalloc((void**)&m_pDevicePk, sizeof(CLGComplex) * m_uiKDim * m_uiMDim));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceHmGm, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim));

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceALeft, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceA, sizeof(CLGComplex) * m_uiMDim * m_uiMDim));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceB, sizeof(CLGComplex) * m_uiMDim * m_uiMDim));

    m_pHostHmGm = (CLGComplex*)malloc(sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim);
    m_pHostHmGmToRotate = (CLGComplex*)malloc(sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim);
    m_pHostY = (CLGComplex*)malloc(sizeof(CLGComplex) * (m_uiMDim + 1));

    m_pHostALeft = (CLGComplex*)malloc(sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim);
    m_pHostB = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiMDim * m_uiMDim);
    m_pHostPk = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiMDim * m_uiKDim);
    m_pHostTmpQ = (CLGComplex*)malloc(sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiKDim);
    m_pHostTmpR = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiKDim * m_uiKDim);
    m_pHostTmpGPk = (CLGComplex*)malloc(sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiKDim);

    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        CField* pVectors = pFieldB->GetCopy();
        m_lstU.AddItem(pVectors);
    }
}

void CSLASolverGCRODR::ReleaseBuffers()
{
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        appSafeDelete(m_lstU[i]);
    }

    if (NULL != m_pDeviceHm)
    {
        checkCudaErrors(cudaFree(m_pDeviceHm));
        checkCudaErrors(cudaFree(m_pDeviceEigenValue));
        checkCudaErrors(cudaFree(m_pDevicePk));
        checkCudaErrors(cudaFree(m_pDeviceHmGm));

        checkCudaErrors(cudaFree(m_pDeviceALeft));
        checkCudaErrors(cudaFree(m_pDeviceA));
        checkCudaErrors(cudaFree(m_pDeviceB));

        free(m_pHostHmGm);
        free(m_pHostHmGmToRotate);
        free(m_pHostY);
        free(m_pHostALeft);
        free(m_pHostB);
        free(m_pHostPk);
        free(m_pHostTmpQ);
        free(m_pHostTmpR);
        free(m_pHostTmpGPk);
    }
}

UBOOL CSLASolverGCRODR::Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pFieldGauge, EFieldOperator uiM, const CField* pStart)
{
    //use it to estimate relative error
    Real fBLength = F(1.0);
    if (!m_bAbsoluteAccuracy)
    {
        fBLength = pFieldB->GetLength();//pFieldB->Dot(pFieldB).x;
    }
    GetPooledFields(pFieldB);

    //set initial gauss x0 = b or pStart
    CField* pX = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);

    if (NULL == pStart)
    {
        pFieldB->CopyTo(pX);
    }
    else
    {
        pStart->CopyTo(pX);
    }

    //The informations:
    //Is it the first time of trajectory? (Do we have Yk?)
    if (m_bHasYk)
    {
        GenerateCUFirstTime(pX, pFieldB, pFieldGauge, uiM);
    }
    else
    {
        FirstTimeGMERESSolve(pX, pFieldB, pFieldGauge, uiM);

        //We have m_pHostHmGm set
        //No matter whether converge, we need Yk for next time solve
        FindPk1();
        //If converge, we no need to update Ck
        GenerateCU(m_fDiviation >= m_fAccuracy * fBLength, TRUE); 
    }

    if (m_fDiviation < m_fAccuracy)
    {
        appParanoiac(_T("-- GCRODR::Finish afeter First GMRES step ----\n"));
        pX->CopyTo(pFieldX);
        pX->Return();
        ReleasePooledFields();
        return TRUE;
    }

    //======================================
    //Real Iteration start here
    //Here we have Uk, Ck set
    
    for (UINT i = 0; i < m_uiReStart; ++i)
    {
        CLGComplex lastDivation = m_pHostY[m_uiMDim];
        NormUkAndSetD();
        CField* pW = m_lstTmp[0];
        //r0 = residue*(v(m+1))
        //vk+1=r0/|r0|
        GetW(m_uiMDim)->CopyTo(GetW(m_uiKDim));
        //Arnoldi
        for (UINT j = m_uiKDim; j < m_uiMDim; ++j)
        {
            CField* vj = GetW(j);
            CField* vjp1 = GetW(j + 1);
            vj->CopyTo(pW);
            //w = A v[j]
            pW->ApplyOperator(uiM, pFieldGauge);
            pW->CopyTo(vjp1);
            for (UINT k = 0; k < m_uiKDim; ++k)
            {
                //B(k, j) = C(k)^dagger AV(j)
                CLGComplex CkH_W = m_lstC[k]->Dot(pW);
                m_pHostHmGm[k * m_uiMDim + j] = CkH_W;
                //v(j+1) = (I - Ck CkH).A v(j)
                vjp1->Axpy(_make_cuComplex(-CkH_W.x, -CkH_W.y), m_lstC[k]);
            }
            for (UINT k = m_uiKDim; k <= j; ++k)
            {
                CField* vk = GetW(k);
                CLGComplex dotc = vk->Dot(vjp1);
                m_pHostHmGm[k * m_uiMDim + j] = dotc;
                //w -= h[k,j] v[k]
                vjp1->Axpy(_make_cuComplex(-dotc.x, -dotc.y), vk);
            }

            //h[j + 1, j] = ||w||
            Real fWNorm = _sqrt(vjp1->Dot(vjp1).x);
            m_pHostHmGm[(j + 1) * m_uiMDim + j] = _make_cuComplex(fWNorm, F(0.0));
            //v[j + 1] = w / ||w||

            vjp1->ScalarMultply(F(1.0) / fWNorm);
        }

        for (UINT j = 0; j <= m_uiMDim; ++j)
        {
            if (j < m_uiKDim)
            {
                m_pHostY[j] = _cuCmulf(lastDivation, m_lstC[j]->Dot(m_lstV[0]));
            }
            else if (j == m_uiKDim)
            {
                m_pHostY[j] = lastDivation;
            }
            else
            {
                m_pHostY[j] = _make_cuComplex(F(0.0), F(0.0));
            }
            //m_pHostY[j] = _cuCmulf(lastDivation, GetW(j)->Dot(m_lstV[0]));
            //appGeneral(_T("before rotation Y=%f %f\n"), m_pHostY[j].x, m_pHostY[j].y);
        }

        
        memcpy(m_pHostHmGmToRotate, m_pHostHmGm, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim);
        CLinearAlgebraHelper::PrintMatrix(m_pHostHmGmToRotate, m_uiMDim + 1, m_uiMDim);
        m_pHelper->RotateHenssenbergHost(m_pHostHmGmToRotate, m_pHostY, m_uiMDim);

        //CLinearAlgebraHelper::PrintMatrix(m_pHostHmGmToRotate, m_uiMDim + 1, m_uiMDim);
        //CLinearAlgebraHelper::PrintMatrix(m_pHostY, m_uiMDim + 1, 1);

        m_pHelper->SolveYHost(m_pHostY, m_pHostHmGmToRotate, 1, m_uiMDim);

        //CLGComplex testY[100];
        CLinearAlgebraHelper::PrintMatrix(m_pHostY, m_uiMDim + 1, 1);
        //m_pHelper->SmallMatrixMult(testY, m_pHostHmGm, m_pHostY, m_uiMDim + 1, m_uiMDim + 1, 1, FALSE, FALSE);
        //CLinearAlgebraHelper::PrintMatrix(testY, m_uiMDim + 1, 1);
        m_fDiviation = __cuCabsSqf(m_pHostY[m_uiMDim]);

        for (UINT j = 0; j < m_uiMDim; ++j)
        {
            pX->Axpy(m_pHostY[j], GetV(j));
        }
        FindPk2();
        //If converge, we no need to update Ck
        GenerateCU(m_fDiviation >= m_fAccuracy * fBLength, FALSE);

        if (m_fDiviation < m_fAccuracy * fBLength)
        {
            appParanoiac(_T("-- GCRODR::Finish afeter %d step, divation = %1.12f ----\n"), i, m_fDiviation);
            pX->CopyTo(pFieldX);
            pX->Return();
            ReleasePooledFields();
            return TRUE;
        }
        appParanoiac(_T("-- GCRODR::Solve operator: After %d step |residue|^2=%1.12f ----\n"), i, m_fDiviation);
    }

    ReleasePooledFields();
    return FALSE;
}

void CSLASolverGCRODR::FirstTimeGMERESSolve(CField* pX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM)
{
    appParanoiac(_T("-- GCRODR::Solve operator: %s-- Fisrt GMRES step ----\n"), __ENUM_TO_STRING(EFieldOperator, uiM).c_str());

    //r = b - A x0
    //v[0] = r.normalize
    //s = x0

    CField* v0 = GetW(0);
    pX->CopyTo(v0); //x0 need to be preserved
    v0->ApplyOperator(uiM, pGaugeFeild, EOCT_Minus); //x0 = -A x0
    v0->AxpyPlus(pFieldB); //x0 = b-Ax0
    m_fBeta = _hostsqrt(v0->Dot(v0).x);
    v0->ScalarMultply(F(1.0) / m_fBeta);  //v[0] = (b - A x0).normalize
    m_pHelper->InitialZeroHost(m_pHostHmGm, m_uiMDim + 1, m_uiMDim);

    //Arnoldi
    for (UINT j = 0; j < m_uiMDim; ++j)
    {
        CField* vj = GetW(j);
        CField* vjp1 = GetW(j + 1);
        vj->CopyTo(vjp1);
        //w = A v[j]
        vjp1->ApplyOperator(uiM, pGaugeFeild);
        
        for (UINT k = 0; k <= j; ++k)
        {
            CField* vk = GetW(k);
            CLGComplex dotc = vk->Dot(vjp1);
            m_pHostHmGm[k * m_uiMDim + j] = dotc;
            //w -= h[k,j] v[k]
            vjp1->Axpy(_make_cuComplex(-dotc.x, -dotc.y), vk);
        }

        //h[j + 1, j] = ||w||
        Real fWNorm = _sqrt(vjp1->Dot(vjp1).x);
        m_pHostHmGm[(j + 1) * m_uiMDim + j] = _make_cuComplex(fWNorm, F(0.0));
        //v[j + 1] = w / ||w||

        vjp1->ScalarMultply(F(1.0) / fWNorm);
    }

    m_pHelper->InitialZeroHost(m_pHostY, m_uiMDim + 1, 1);
    m_pHostY[0] = _make_cuComplex(m_fBeta, F(0.0));

    memcpy(m_pHostHmGmToRotate, m_pHostHmGm, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim);
    m_pHelper->RotateHenssenbergHost(m_pHostHmGmToRotate, m_pHostY, m_uiMDim);
    
    m_fDiviation = __cuCabsSqf(m_pHostY[m_uiMDim]);
    appParanoiac(_T("-- GCRODR::Solve operator: After Fisrt GMRES step |residue|^2=%1.12f ----\n"), m_fDiviation);
    m_pHelper->SolveYHost(m_pHostY, m_pHostHmGmToRotate, 1, m_uiMDim);

    for (UINT j = 0; j < m_uiMDim; ++j)
    {
        pX->Axpy(m_pHostY[j], GetW(j));
    }
}

/**
* x = x +- AB+ x
*/
void CSLASolverGCRODR::AxpyABdagger(CField* x, const TArray<CField*>& lstA, const TArray<CField*>& lstB, UBOOL bPlusOrMinus, UINT uiK)
{
    for (UINT i = 0; i < uiK; ++i)
    {
        CLGComplex bdaggerx = lstB[i]->Dot(x);
        if (bPlusOrMinus)
        {
            x->Axpy(bdaggerx, lstA[i]);
        }
        else 
        {
            x->Axpy(_make_cuComplex(-bdaggerx.x, -bdaggerx.y), lstA[i]);
        }
    }
}

void CSLASolverGCRODR::QRFactorAY(const CFieldGauge* pGaugeField, EFieldOperator uiM)
{
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        for (UINT j = 0; j < m_uiKDim; ++j)
        {
            m_pHostTmpR[i * m_uiKDim + j] = _make_cuComplex(F(0.0), F(0.0));
        }

        //transform Y to AY
        m_lstU[i]->CopyTo(m_lstC[i]);
        m_lstC[i]->ApplyOperator(uiM, pGaugeField);
    }

    //QR of AY
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        m_pHostTmpR[i * m_uiKDim + i] = m_lstC[i]->Dot(m_lstC[i]);
        m_lstC[i]->ScalarMultply(F(1.0) / _hostsqrt(m_pHostTmpR[i * m_uiKDim + i].x));
        for (UINT j = i + 1; j < m_uiKDim; ++j)
        {
            m_pHostTmpR[i * m_uiKDim + j] = m_lstC[i]->Dot(m_lstC[j]);
            m_lstC[j]->Axpy(
                _make_cuComplex(
                    -m_pHostTmpR[i * m_uiKDim + j].x,
                    -m_pHostTmpR[i * m_uiKDim + j].y),
                m_lstC[i]);
        }
    }
}

void CSLASolverGCRODR::FindPk1()
{
    checkCudaErrors(cudaMemcpy(m_pDeviceHmGm, m_pHostHmGm, sizeof(CLGComplex) * m_uiMDim * (m_uiMDim + 1), cudaMemcpyHostToDevice));

    switch (m_eDeflationType)
    {
    case EEDT_REV:
    {
        checkCudaErrors(cudaMemcpy(m_pDeviceHm, m_pDeviceHmGm, sizeof(CLGComplex) * m_uiMDim * m_uiMDim, cudaMemcpyDeviceToDevice));
        m_pHelper->EigenValueProblemHessenberg(m_pDeviceHm, m_pDeviceEigenValue, m_pDevicePk, m_uiMDim, m_uiKDim);
    }
    break;
    case EEDT_SVD:
    {
        m_pHelper->SmallMatrixMult(m_pDeviceHm, m_pDeviceHmGm, m_pDeviceHmGm, m_uiMDim, m_uiMDim + 1, m_uiMDim, TRUE, FALSE);
        m_pHelper->EigenValueProblem(m_pDeviceHm, m_pDeviceEigenValue, m_pDevicePk, m_uiMDim, m_uiKDim);
    }
    break;
    }

#if _CLG_DEBUG
    appParanoiac(_T("GCRO-DR: EigenValues:\n"));
    m_pHelper->PrintDeviceMatrix(m_pDeviceEigenValue, 1, m_uiKDim);
#endif

    m_pHelper->Transpose(m_pDevicePk, m_uiKDim, m_uiMDim);
    checkCudaErrors(cudaMemcpy(m_pHostPk, m_pDevicePk, sizeof(CLGComplex) * m_uiMDim * m_uiKDim, cudaMemcpyDeviceToHost));
}

void CSLASolverGCRODR::FindPk2()
{
    checkCudaErrors(cudaMemcpy(m_pDeviceHmGm, m_pHostHmGm,
        sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim, cudaMemcpyHostToDevice));
    if (EEDT_REV == m_eDeflationType)
    {
        //Aleft
        for (UINT i = 0; i < m_uiMDim; ++i)
        {
            for (UINT j = 0; j < m_uiMDim + 1; ++j)
            {
                if (i < m_uiKDim && j < m_uiKDim)
                {
                    //Ui dagger Cj
                    m_pHostALeft[i * (m_uiMDim + 1) + j] = m_lstU[i]->Dot(m_lstC[j]);
                }
                else if (i < m_uiKDim && j >= m_uiKDim)
                {
                    //Ui dagger V[j + m_uiK]
                    m_pHostALeft[i * (m_uiMDim + 1) + j] = m_lstU[i]->Dot(m_lstV[j - m_uiKDim]);
                }
                else // i >= m_uiKDim
                {
                    if (i == j)
                    {
                        m_pHostALeft[i * (m_uiMDim + 1) + j] = _make_cuComplex(F(1.0), F(0.0));
                    }
                    else
                    {
                        m_pHostALeft[i * (m_uiMDim + 1) + j] = _make_cuComplex(F(0.0), F(0.0));
                    }
                }
            }
        }

        //B
        for (UINT i = 0; i < m_uiMDim; ++i)
        {
            for (UINT j = 0; j < m_uiMDim; ++j)
            {
                if (i < m_uiKDim && j < m_uiKDim)
                {
                    if (i == j)
                    {
                        //Uk is already normalized
                        m_pHostB[i * m_uiMDim + j] = _make_cuComplex(F(1.0), F(0.0));
                    }
                    else
                    {
                        m_pHostB[i * m_uiMDim + j] = m_lstU[i]->Dot(m_lstU[j]);
                    }
                }
                else if (i < m_uiKDim && j >= m_uiKDim)
                {
                    //Ui dagger Vj
                    m_pHostB[i * m_uiMDim + j] = m_pHostALeft[i * (m_uiMDim + 1) + j];
                }
                else if (i >= m_uiKDim && j < m_uiKDim)
                {
                    //(Uj dagger Vi)^*
                    m_pHostB[i * m_uiMDim + j] = _cuConjf(m_pHostALeft[j * (m_uiMDim + 1) + i]);
                }
                else  // i >= m_uiKDim && j >= m_uiKDim
                {
                    //(Ui dagger Vj)^*
                    if (i == j)
                    {
                        m_pHostB[i * m_uiMDim + j] = _make_cuComplex(F(1.0), F(0.0));
                    }
                    else
                    {
                        m_pHostB[i * m_uiMDim + j] = _make_cuComplex(F(0.0), F(0.0));
                    }
                }
            }
        }

        //CLinearAlgebraHelper::PrintMatrix(m_pHostALeft, m_uiMDim, m_uiMDim + 1);
        //CLinearAlgebraHelper::PrintMatrix(m_pHostHmGm, m_uiMDim + 1, m_uiMDim);
        //CLinearAlgebraHelper::PrintMatrix(m_pHostB, m_uiMDim, m_uiMDim);
        //A
        checkCudaErrors(cudaMemcpy(m_pDeviceALeft, m_pHostALeft, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(m_pDeviceB, m_pHostB, sizeof(CLGComplex) * m_uiMDim * m_uiMDim, cudaMemcpyHostToDevice));
        m_pHelper->SmallMatrixMult(m_pDeviceA, m_pDeviceALeft, m_pDeviceHmGm, m_uiMDim, m_uiMDim + 1, m_uiMDim, FALSE, FALSE);

        m_pHelper->GeneralizedEigenValueProblem(m_pDeviceA, m_pDeviceB, m_pDeviceEigenValue, m_pDevicePk, m_uiMDim, m_uiKDim);
    }
    else if (EEDT_SVD == m_eDeflationType)
    {
        //B
        for (UINT i = 0; i < m_uiMDim; ++i)
        {
            for (UINT j = 0; j < m_uiMDim; ++j)
            {
                if (i < m_uiKDim && j < m_uiKDim)
                {
                    if (i == j)
                    {
                        //Uk is already normalized
                        m_pHostB[i * m_uiMDim + j] = _make_cuComplex(F(1.0), F(0.0));
                    }
                    else
                    {
                        m_pHostB[i * m_uiMDim + j] = m_lstU[i]->Dot(m_lstU[j]);
                    }
                }
                else  // i >= m_uiKDim || j >= m_uiKDim
                {
                    if (i == j)
                    {
                        m_pHostB[i * m_uiMDim + j] = _make_cuComplex(F(1.0), F(0.0));
                    }
                    else
                    {
                        m_pHostB[i * m_uiMDim + j] = _make_cuComplex(F(0.0), F(0.0));
                    }
                }
            }
        }

        checkCudaErrors(cudaMemcpy(m_pDeviceB, m_pHostB, sizeof(CLGComplex) * m_uiMDim * m_uiMDim, cudaMemcpyHostToDevice));
        m_pHelper->SmallMatrixMult(m_pDeviceA, m_pDeviceHmGm, m_pDeviceHmGm, m_uiMDim, m_uiMDim + 1, m_uiMDim, TRUE, FALSE);

        m_pHelper->GeneralizedEigenValueProblem(m_pDeviceA, m_pDeviceB, m_pDeviceEigenValue, m_pDevicePk, m_uiMDim, m_uiKDim);
    }
    m_pHelper->Transpose(m_pDevicePk, m_uiKDim, m_uiMDim);
#if _CLG_DEBUG
    appParanoiac(_T("GCRO-DR: EigenValues:\n"));
    m_pHelper->PrintDeviceMatrix(m_pDeviceEigenValue, 1, m_uiKDim);
#endif
    checkCudaErrors(cudaMemcpy(m_pHostPk, m_pDevicePk, sizeof(CLGComplex) * m_uiMDim * m_uiKDim, cudaMemcpyDeviceToHost));
}

/**
* Assume m >= k
* V=(vk[0], ... , vk[k - 1])
* W=(vk[0], ... , vk[k - 1], vmk[0], ..., vmk[m-k-1])
*
* V(v1,v2,...,vk) = W(w1,w2,...,wm) (m11, ..., m1k)
*                                   (..., ..., ...)
*                                   (mm1, ..., mmk)
* I think this is expansive... the FLOP of Ax is about 100n, but this has m x k x n
*/
void CSLASolverGCRODR::VectorMultiplyMatrix(
    //TArray<class CField*>& resultV, 
    const TArray<class CField*>& lstYorC, 
    //const TArray<class CField*>& lstVmk, 
    CLGComplex* hostMatrix, UINT uiMDim)
{
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        lstYorC[0]->CopyTo(m_lstTmp[i]);
        m_lstTmp[i]->ScalarMultply(hostMatrix[i]);

        for (UINT j = 1; j < uiMDim; ++j)
        {
            m_lstTmp[i]->Axpy(hostMatrix[j * m_uiKDim + i],
                j < m_uiKDim ? lstYorC[j] : m_lstV[j - m_uiKDim]);
        }

        //appGeneral(_T("================= %d ================\n"), i);
        //m_lstTmp[i]->DebugPrintMe();
    }
}

/**
* Assume m >= k
* V=(vk[0], ... , vk[k - 1])
* W=(vk[0], ... , vk[k - 1], vmk[0], ..., vmk[m-k-1])
*
* v1 = (m11, ..., m1m)  w1      
* ..   (..., ..., ...)  ..      
* vk   (mk1, ..., mkm)  wk      
*                       wk+1
*                       ...
*                       wm
*/
//void CSLASolverGCRODR::MatrixMultiplyVector(TArray<class CField*>& resultV, const TArray<class CField*>& lstVk, const TArray<class CField*>& lstVmk, CLGComplex* hostMatrix, 
//    UINT uiKDim, UINT uiMDim, UINT uiMDimY)
//{
//    for (UINT i = 0; i < uiKDim; ++i)
//    {
//        lstVk[0]->CopyTo(resultV[i]);
//        resultV[i]->ScalarMultply(hostMatrix[i * uiMDim]);
//
//        for (UINT j = 1; j < uiMDim; ++j)
//        {
//            resultV[i]->Axpy(hostMatrix[i * uiMDim + j],
//                j < uiKDim ? lstVk[j] : lstVmk[j - uiMDim]);
//        }
//    }
//}

// YR-1
void CSLASolverGCRODR::FieldSolveY(TArray<class CField*>& resultY, const CLGComplex* R, UINT uiDim)
{
    for (UINT i = 0; i < uiDim; ++i)
    {
        for (UINT j = 0; j < i; ++j)
        {
            //appGeneral(_T("r = %f %f\n"), R[j * uiDim + i].x, R[j * uiDim + i].y);
            resultY[i]->Axpy(
                _make_cuComplex(-R[j * uiDim + i].x, -R[j * uiDim + i].y), 
                resultY[j]);
        }
        CLGComplex divider = _cuCdivf(_make_cuComplex(F(1.0), F(0.0)), R[i * uiDim + i]);
        //appGeneral(_T("r = %f %f and divider = %f %f\n"), R[i * uiDim * i].x, R[i * uiDim * i].y, divider.x, divider.y);
        resultY[i]->ScalarMultply(divider);
    }
}

/**
* We assume m_pHostPk is obtained.
* For the last iteration, we still update Uk for the next solve, this time no need to update Ck
*/
void CSLASolverGCRODR::GenerateCU(UBOOL bUpdateCk, UBOOL bJustAfterGMRES)
{
    
    //QR = Gm(Hm) Pk
    m_pHelper->SmallMatrixMultHost(
        m_pHostTmpGPk, m_pHostHmGm, m_pHostPk, 
        m_uiMDim + 1, m_uiMDim, m_uiKDim, FALSE, FALSE);
    m_pHelper->ThinQRFactorizationHost(
        m_pHostTmpQ, m_pHostTmpR, m_pHostTmpGPk, 
        m_uiMDim + 1, m_uiKDim);

    
    //CLinearAlgebraHelper::PrintMatrix(m_pHostPk, m_uiMDim, m_uiKDim);
    //CLinearAlgebraHelper::PrintMatrix(m_pHostTmpR, m_uiKDim, m_uiKDim);
    //Vm = (Uk, V) NOTE: If it is just after GMRES, Vm = (Ck, V)
    //Uk = Vm Pk (Note this is NOT GPk)
    //Uk = Uk R-1
    //for (UINT i = 0; i < m_uiMDim - m_uiKDim + 1; ++i)
    //{
    //    m_lstV[i]->DebugPrintMe();
    //}

    VectorMultiplyMatrix(
        //m_lstTmp, 
        bJustAfterGMRES ? m_lstC : m_lstU, 
        //m_lstV, 
        m_pHostPk, m_uiMDim);

    FieldSolveY(m_lstTmp, m_pHostTmpR, m_uiKDim);
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        m_lstTmp[i]->CopyTo(m_lstU[i]);
    }

    //for (UINT i = 0; i < m_uiKDim; ++i)
    //{
    //    m_lstU[i]->DebugPrintMe();
    //}

    //CLinearAlgebraHelper::PrintMatrix(m_pHostTmpQ, m_uiMDim + 1, m_uiKDim);
    if (bUpdateCk)
    {
        //Ck = W Q, W=(Ck, V)
        //VectorMultiplyMatrix(m_lstTmp, m_lstC, m_lstV, m_pHostTmpQ, m_uiMDim + 1);
        VectorMultiplyMatrix(m_lstC, m_pHostTmpQ, m_uiMDim + 1);
        for (UINT i = 0; i < m_uiKDim; ++i)
        {
            m_lstTmp[i]->CopyTo(m_lstC[i]);
        }
    }

    //for (UINT i = 0; i < m_uiKDim; ++i)
    //{
    //    m_lstC[i]->DebugPrintMe();
    //}
}

/**
* We have Yk, (which is Uk), QR is QR of AYk
*/
void CSLASolverGCRODR::GenerateCUFirstTime(CField* pX, const CField* pFieldB, const CFieldGauge* pGaugeField, EFieldOperator uiM)
{
    QRFactorAY(pGaugeField, uiM);

    //Uk = Uk R-1
    FieldSolveY(m_lstTmp, m_pHostTmpR, m_uiKDim);
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        m_lstTmp[i]->CopyTo(m_lstU[i]);
    }

    CField* v0 = m_lstV[0];
    CField* vmp1 = m_lstV[m_uiMDim - m_uiKDim];

    pX->CopyTo(v0);
    v0->ApplyOperator(uiM, pGaugeField, EOCT_Minus); //x0 = -A x0
    v0->AxpyPlus(pFieldB); //x0 = b-Ax0
    v0->CopyTo(vmp1);
    
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        CLGComplex CkH_R0 = m_lstC[i]->Dot(v0);
        //x=x+Uk CkH r0
        pX->Axpy(CkH_R0, m_lstU[i]);
        //r=r-Ck CkH r0
        vmp1->Axpy(_make_cuComplex(-CkH_R0.x, -CkH_R0.y), m_lstC[i]);
    }
    m_pHostY[m_uiMDim] = vmp1->Dot(vmp1);
    vmp1->ScalarMultply(_hostsqrt(F(1.0) / m_pHostY[m_uiMDim].x));
}

/**
* If in here, Uk is already known
*/
void CSLASolverGCRODR::NormUkAndSetD()
{
    m_pHelper->InitialZeroHost(m_pHostHmGm, m_uiMDim + 1, m_uiMDim);

    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        CLGComplex dotres = m_lstU[i]->Dot(m_lstU[i]);
        Real fLength = _hostsqrt(dotres.x);
        m_pHostHmGm[i * m_uiMDim + i] = _make_cuComplex(fLength, F(0.0));
        m_lstU[i]->ScalarMultply(F(1.0) / fLength);
    }
}

void CSLASolverGCRODR::GetPooledFields(const CField* pFieldB)
{
    assert(0 == m_lstV.Num());
    assert(0 == m_lstC.Num());
    assert(0 == m_lstTmp.Num());
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        CField* pVectors = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
        m_lstC.AddItem(pVectors);
        pVectors = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
        m_lstTmp.AddItem(pVectors);
    }
    for (UINT i = 0; i < m_uiMDim - m_uiKDim + 1; ++i)
    {
        CField* pVectors = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
        m_lstV.AddItem(pVectors);
    }
}

void CSLASolverGCRODR::ReleasePooledFields()
{
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        m_lstC[i]->Return();
        m_lstTmp[i]->Return();
    }
    for (UINT i = 0; i < m_uiMDim - m_uiKDim + 1; ++i)
    {
        m_lstV[i]->Return();
    }
    m_lstC.RemoveAll();
    m_lstTmp.RemoveAll();
    m_lstV.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================