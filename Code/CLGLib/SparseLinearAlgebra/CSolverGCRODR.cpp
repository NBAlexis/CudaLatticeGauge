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
    , m_pHostZeroMatrix(NULL)

    , m_pFieldMatrix(NULL)
    , m_uiRecalcuateR(5)
{

}

CSLASolverGCRODR::~CSLASolverGCRODR()
{
    ReleaseBuffers();
    appSafeDelete(m_pHelper);
    appSafeDelete(m_pFieldMatrix);
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
    if (param.FetchValueINT(_T("RecalculateR"), iValue))
    {
        m_uiRecalcuateR = static_cast<UINT>(iValue);
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

    
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceTmpQ, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiKDim));

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
    m_pHostZeroMatrix = (CLGComplex*)malloc(sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim);

    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        CField* pVectors = pFieldB->GetCopy();
        m_lstU.AddItem(pVectors);
    }
    for (UINT i = 0; i < m_uiMDim * (m_uiMDim + 1); ++i)
    {
        m_pHostZeroMatrix[i] = _make_cuComplex(F(0.0), F(0.0));
    }
    m_pFieldMatrix = CFieldMatrixOperation::Create(pFieldB->GetFieldType());
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
        checkCudaErrors(cudaFree(m_pDeviceTmpQ));
        
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
        free(m_pHostZeroMatrix);
    }
}

UBOOL CSLASolverGCRODR::Solve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pFieldGauge, EFieldOperator uiM, ESolverPhase ePhase, const CField* pStart)
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
    CField* pR = m_lstV[0];//appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
    CField* pW = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);

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
    if (ESP_InTrajectory == ePhase || ESP_EndTrajectory == ePhase)
    {
        GenerateCUFirstTime(pX, pR, pFieldB, pFieldGauge, uiM);
    }
    else
    {
        FirstTimeGMERESSolve(pX, pW, pFieldB, pFieldGauge, uiM);

        //We have m_pHostHmGm set
        //No matter whether converge, we need Yk for next time solve
        FindPk1();
        //If converge, we no need to update Ck
        GenerateCU(m_fDiviation >= m_fAccuracy * fBLength, TRUE); 
        pW->CopyTo(pR);
    }

    if (m_fDiviation < m_fAccuracy * fBLength)
    {
        appParanoiac(_T("-- GCRODR::Finish afeter First GMRES(or CU step) step ----\n"));
        pX->CopyTo(pFieldX);
        pX->Return();
        pW->Return();
        ReleasePooledFields();
        return TRUE;
    }

    //======================================
    //Real Iteration start here
    //Here we have Uk, Ck set
    for (UINT i = 0; i < m_uiReStart; ++i)
    {
        OrthognalXR(pX, pR, pW);

        //m_pHostHmGm is zeroed in NormUkAndSetD
        NormUkAndSetD();
        
        //vk+1=r0/|r0|
        //pR->CopyTo(GetW(m_uiKDim));
        m_fDiviation = _hostsqrt(pR->Dot(pR).x);
        GetW(m_uiKDim)->ScalarMultply(F(1.0) / m_fDiviation);
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
            Real fWNorm = _hostsqrt(vjp1->Dot(vjp1).x);
            m_pHostHmGm[(j + 1) * m_uiMDim + j] = _make_cuComplex(fWNorm, F(0.0));
            //v[j + 1] = w / ||w||
            vjp1->ScalarMultply(F(1.0) / fWNorm);
        }

        //pW is not in use now, we use pW as pR
        pR->CopyTo(pW);
        pW->ScalarMultply(m_fDiviation);

        //=============== This is almost accurate =============
        memcpy(m_pHostY, m_pHostZeroMatrix, sizeof(CLGComplex) * (m_uiMDim + 1));
        m_pHostY[m_uiKDim] = _make_cuComplex(m_fDiviation, F(0.0));

        //for (UINT j = 0; j <= m_uiMDim; ++j)
        //{
        //    m_pHostY[j] = GetW(j)->Dot(pR);
        //}
        
        memcpy(m_pHostHmGmToRotate, m_pHostHmGm, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim);

        m_pHelper->RotateHenssenbergHost(m_pHostHmGmToRotate, m_pHostY, m_uiMDim);
        m_pHelper->SolveYHost(m_pHostY, m_pHostHmGmToRotate, 1, m_uiMDim);

        for (UINT j = 0; j < m_uiMDim; ++j)
        {
            pX->Axpy(m_pHostY[j], GetV(j));
        }

        //v[0] are still used in Eigen-Value problem, so we use pW instead
        if (0 != m_uiRecalcuateR && 0 != i && 0 == (i % m_uiRecalcuateR))
        {
            //============== This is the accurate result, though, slower and not stable =================
            pX->CopyTo(pW);
            pW->ApplyOperator(uiM, pFieldGauge, EOCT_Minus); //x0 = -A x0
            pW->AxpyPlus(pFieldB); //x0 = b-Ax0
            m_cLastDiviation = pW->Dot(pW);
            m_fDiviation = _hostsqrt(m_cLastDiviation.x);
        }
        else
        {
            //============= This is not accurate enough, but stable and faster =================
            //y(m+1) = Gm(m+1,m) y(m)
            m_pHelper->SmallMatrixMultHost(m_pHostY, m_pHostHmGm, m_pHostY, m_uiMDim + 1, m_uiMDim, 1, FALSE, FALSE);
            for (UINT j = 0; j < m_uiMDim + 1; ++j)
            {
                pW->Axpy(_make_cuComplex(-m_pHostY[j].x, -m_pHostY[j].y), GetW(j));
            }
            m_cLastDiviation = pW->Dot(pW);
            m_fDiviation = _hostsqrt(m_cLastDiviation.x);
        }

        FindPk2();
        //If converge, we no need to update Ck
        if (m_fDiviation >= m_fAccuracy * fBLength || ESP_EndTrajectory != ePhase)
        {
            //If is end trajectory, no need to update both Uk and Ck
            GenerateCU(m_fDiviation >= m_fAccuracy * fBLength, FALSE);
        }

        if (m_fDiviation < m_fAccuracy * fBLength)
        {
            appParanoiac(_T("-- GCRODR::Finish afeter %d step, divation = %1.15f ----\n"), i, m_fDiviation);
            pX->CopyTo(pFieldX);
            pX->Return();
            pW->Return();
            ReleasePooledFields();
            return TRUE;
        }
        appParanoiac(_T("-- GCRODR::Solve operator: After %d step |residue|=%1.15f ----\n"), i, m_fDiviation);
        pW->CopyTo(pR);
    }

    appParanoiac(_T("GCRODR::Solve failed: last divation |residue| = %8.15f\n"), m_fDiviation);
    pX->CopyTo(pFieldX);
    pX->Return();
    pW->Return();
    ReleasePooledFields();
    return FALSE;
}

void CSLASolverGCRODR::FirstTimeGMERESSolve(CField* pX, CField* pR, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM)
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
    memcpy(m_pHostHmGm, m_pHostZeroMatrix, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim);

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

    memcpy(m_pHostY, m_pHostZeroMatrix, sizeof(CLGComplex) * (m_uiMDim + 1));
    m_pHostY[0] = _make_cuComplex(m_fBeta, F(0.0));

    memcpy(m_pHostHmGmToRotate, m_pHostHmGm, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim);
    m_pHelper->RotateHenssenbergHost(m_pHostHmGmToRotate, m_pHostY, m_uiMDim);
    m_pHelper->SolveYHost(m_pHostY, m_pHostHmGmToRotate, 1, m_uiMDim);

    for (UINT j = 0; j < m_uiMDim; ++j)
    {
        pX->Axpy(m_pHostY[j], GetW(j));
    }

    //================ This is not accurate enough =============
    if (m_uiRecalcuateR > 1)
    {
        m_pHelper->SmallMatrixMultHost(m_pHostY, m_pHostHmGm, m_pHostY, m_uiMDim + 1, m_uiMDim, 1, FALSE, FALSE);
        m_pHostY[0].x -= m_fBeta;
        pR->Zero();
        for (UINT j = 0; j < m_uiMDim + 1; ++j)
        {
            pR->Axpy(_make_cuComplex(-m_pHostY[j].x, -m_pHostY[j].y), GetW(j));
        }
        m_cLastDiviation = pR->Dot(pR);
        m_fDiviation = _hostsqrt(m_cLastDiviation.x);
    }
    else
    {
        pX->CopyTo(pR);
        pR->ApplyOperator(uiM, pGaugeFeild, EOCT_Minus); //x0 = -A x0
        pR->AxpyPlus(pFieldB); //x0 = b-Ax0
        m_cLastDiviation = pR->Dot(pR);
        m_fDiviation = _hostsqrt(m_cLastDiviation.x);
    }

    appParanoiac(_T("-- GCRODR::Solve operator: After Fisrt GMRES step |residue|=%1.12f ----\n"), m_fDiviation);
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
        Real fLength = _hostsqrt(m_lstC[i]->Dot(m_lstC[i]).x);
        m_pHostTmpR[i * m_uiKDim + i] = _make_cuComplex(fLength, F(0.0));
        m_lstC[i]->ScalarMultply(F(1.0) / fLength);
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
    switch (m_eDeflationType)
    {
    case EEDT_REV:
    {
        checkCudaErrors(cudaMemcpy(m_pDeviceHm, m_pHostHmGm, sizeof(CLGComplex) * m_uiMDim * m_uiMDim, cudaMemcpyHostToDevice));
        m_pHelper->EigenValueProblemHessenberg(m_pDeviceHm, m_pDeviceEigenValue, m_pDevicePk, m_uiMDim, m_uiKDim);
    }
    break;
    case EEDT_HEV:
    {
        checkCudaErrors(cudaMemcpy(m_pDeviceHmGm, m_pHostHmGm, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(m_pDeviceA, m_pDeviceHmGm, sizeof(CLGComplex) * m_uiMDim * m_uiMDim, cudaMemcpyDeviceToDevice));
        m_pHelper->SmallMatrixMult(m_pDeviceB, m_pDeviceHmGm, m_pDeviceHmGm, m_uiMDim, m_uiMDim + 1, m_uiMDim, TRUE, FALSE);
        m_pHelper->GeneralizedEigenValueProblem(m_pDeviceA, m_pDeviceB, m_pDeviceEigenValue, m_pDevicePk, m_uiMDim, m_uiKDim, FALSE);
    }
    break;
    case EEDT_SVD:
    {
        checkCudaErrors(cudaMemcpy(m_pDeviceHmGm, m_pHostHmGm, sizeof(CLGComplex) * m_uiMDim * (m_uiMDim + 1), cudaMemcpyHostToDevice));
        m_pHelper->SmallMatrixMult(m_pDeviceHm, m_pDeviceHmGm, m_pDeviceHmGm, m_uiMDim, m_uiMDim + 1, m_uiMDim, TRUE, FALSE);
        m_pHelper->EigenValueProblem(m_pDeviceHm, m_pDeviceEigenValue, m_pDevicePk, m_uiMDim, m_uiKDim);
    }
    break;
    }

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

        //A
        checkCudaErrors(cudaMemcpy(m_pDeviceALeft, m_pHostALeft, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim, cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemcpy(m_pDeviceB, m_pHostB, sizeof(CLGComplex) * m_uiMDim * m_uiMDim, cudaMemcpyHostToDevice));
        m_pHelper->SmallMatrixMult(m_pDeviceA, m_pDeviceALeft, m_pDeviceHmGm, m_uiMDim, m_uiMDim + 1, m_uiMDim, FALSE, FALSE);

        m_pHelper->GeneralizedEigenValueProblem(m_pDeviceA, m_pDeviceB, m_pDeviceEigenValue, m_pDevicePk, m_uiMDim, m_uiKDim);
    }
    else if (EEDT_HEV == m_eDeflationType)
    {
        //Aleft
        for (UINT i = 0; i < m_uiMDim + 1; ++i)
        {
            for (UINT j = 0; j < m_uiMDim; ++j)
            {
                if (i < m_uiKDim && j < m_uiKDim)
                {
                    //Ui dagger Cj
                    m_pHostALeft[i * m_uiMDim + j] = m_lstC[i]->Dot(m_lstU[j]);
                }
                else if (i >= m_uiKDim && j < m_uiKDim)
                {
                    //Ui dagger V[j + m_uiK]
                    m_pHostALeft[i * m_uiMDim + j] = m_lstV[i - m_uiKDim]->Dot(m_lstU[j]);
                }
                else // j >= m_uiKDim
                {
                    if (i == j)
                    {
                        m_pHostALeft[i * m_uiMDim + j] = _make_cuComplex(F(1.0), F(0.0));
                    }
                    else
                    {
                        m_pHostALeft[i * m_uiMDim + j] = _make_cuComplex(F(0.0), F(0.0));
                    }
                }
            }
        }
        checkCudaErrors(cudaMemcpy(m_pDeviceALeft, m_pHostALeft, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim, cudaMemcpyHostToDevice));
        m_pHelper->SmallMatrixMult(m_pDeviceA, m_pDeviceHmGm, m_pDeviceALeft, m_uiMDim, m_uiMDim + 1, m_uiMDim, TRUE, FALSE);
        m_pHelper->SmallMatrixMult(m_pDeviceB, m_pDeviceHmGm, m_pDeviceHmGm, m_uiMDim, m_uiMDim + 1, m_uiMDim, TRUE, FALSE);
        m_pHelper->GeneralizedEigenValueProblem(m_pDeviceA, m_pDeviceB, m_pDeviceEigenValue, m_pDevicePk, m_uiMDim, m_uiKDim, FALSE);
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
    //appParanoiac(_T("GCRO-DR: EigenValues:\n"));
    //m_pHelper->PrintDeviceMatrix(m_pDeviceEigenValue, 1, m_uiKDim);
#endif
    checkCudaErrors(cudaMemcpy(m_pHostPk, m_pDevicePk, sizeof(CLGComplex) * m_uiMDim * m_uiKDim, cudaMemcpyDeviceToHost));
}


/**
* +
* U=YR-1
* U^T=(R^T)-1Y^T
*/
void CSLASolverGCRODR::FieldSolveY(TArray<class CField*>& resultY, const CLGComplex* R, UINT uiDim)
{
    for (UINT i = 0; i < uiDim; ++i)
    {
        for (UINT j = 0; j < i; ++j)
        {
            resultY[i]->Axpy(
                _make_cuComplex(-R[j * uiDim + i].x, -R[j * uiDim + i].y), 
                resultY[j]);
        }
        CLGComplex divider = _cuCdivf(_make_cuComplex(F(1.0), F(0.0)), R[i * uiDim + i]);
        resultY[i]->ScalarMultply(divider);
    }
}

/**
* We assume m_pDevicePk and m_pHostPk are obtained.
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

    checkCudaErrors(cudaMemcpy(m_pDeviceTmpQ, m_pHostTmpQ, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiKDim, cudaMemcpyHostToDevice));

    //Vm = (Uk, V) NOTE: If it is just after GMRES, Vm = (Ck, V)
    //Uk = Vm Pk (Note this is NOT GPk)
    //Uk = Uk R-1
    if (bJustAfterGMRES)
    {
        for (UINT i = 0; i < m_uiKDim; ++i)
        {
            m_lstC[i]->CopyTo(m_lstU[i]);
        }
    }
    m_pFieldMatrix->VectorMultiplyMatrix(m_lstU, m_lstV, m_pDevicePk, m_uiMDim, m_uiKDim);
    FieldSolveY(m_lstU, m_pHostTmpR, m_uiKDim);

    if (bUpdateCk)
    {
        //Ck = W Q, W=(Ck, V)
        m_pFieldMatrix->VectorMultiplyMatrix(m_lstC, m_lstV, m_pDeviceTmpQ, m_uiMDim + 1, m_uiKDim);
    }
}

/**
* We have Yk, (which is Uk), QR is QR of AYk
*/
void CSLASolverGCRODR::GenerateCUFirstTime(CField* pX, CField* pR, const CField* pFieldB, const CFieldGauge* pGaugeField, EFieldOperator uiM)
{
    QRFactorAY(pGaugeField, uiM);

    //Uk = Uk R-1
    FieldSolveY(m_lstU, m_pHostTmpR, m_uiKDim);

    //CField* v0 = m_lstV[0];

    pX->CopyTo(pR);
    pR->ApplyOperator(uiM, pGaugeField, EOCT_Minus); //r0 = -A x0
    pR->AxpyPlus(pFieldB); //r0 = b-Ax0
    
    m_cLastDiviation = pR->Dot(pR);
    m_fDiviation = _hostsqrt(__cuCabsSqf(m_cLastDiviation));
}

/**
* If in here, Uk is already known
*/
void CSLASolverGCRODR::NormUkAndSetD()
{
    memcpy(m_pHostHmGm, m_pHostZeroMatrix, sizeof(CLGComplex) * (m_uiMDim + 1) * m_uiMDim);

    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        CLGComplex dotres = m_lstU[i]->Dot(m_lstU[i]);
        Real fLength = F(1.0) / _hostsqrt(dotres.x);
        m_pHostHmGm[i * m_uiMDim + i] = _make_cuComplex(fLength, F(0.0));
        m_lstU[i]->ScalarMultply(fLength);
    }
}

void CSLASolverGCRODR::OrthognalXR(CField* pX, CField* pR, CField* pTmp)
{
    pR->CopyTo(pTmp);
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        CLGComplex CkH_R0 = m_lstC[i]->Dot(pTmp);
        //x=x+Uk CkH r0
        pX->Axpy(CkH_R0, m_lstU[i]);
        //r=r-Ck CkH r0
        pR->Axpy(_make_cuComplex(-CkH_R0.x, -CkH_R0.y), m_lstC[i]);
    }
}

void CSLASolverGCRODR::GetPooledFields(const CField* pFieldB)
{
    assert(0 == m_lstV.Num());
    assert(0 == m_lstC.Num());
    for (UINT i = 0; i < m_uiKDim; ++i)
    {
        CField* pVectors = appGetLattice()->GetPooledFieldById(pFieldB->m_byFieldId);
        m_lstC.AddItem(pVectors);
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
    }
    for (UINT i = 0; i < m_uiMDim - m_uiKDim + 1; ++i)
    {
        m_lstV[i]->Return();
    }
    m_lstC.RemoveAll();
    m_lstV.RemoveAll();
}


__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================