//=============================================================================
// FILENAME : CLinearAlgebraHelper.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [03/14/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

CLinearAlgebraHelper::CLinearAlgebraHelper(UINT uiDim)
    : m_uiDim(uiDim)
    , m_pDeviceIntBuffer(NULL)
    , m_pDeviceFloatBuffer(NULL)
    , m_pDeviceComplexBuffer1(NULL)
    , m_pDeviceComplexBuffer2(NULL)
{
    if (uiDim > _kMaxSmallDim)
    {
        appGeneral(_T("CLinearAlgebraHelper only support dim <= 32 !!\n"));
        m_uiDim = 0;
        return;
    }

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceIntBuffer, sizeof(INT) * m_uiDim));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceFloatBuffer, sizeof(Real) * m_uiDim));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceComplexBuffer1, sizeof(CLGComplex) * m_uiDim));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceComplexBuffer2, sizeof(CLGComplex) * m_uiDim));

    //5 is enough
    AddTempMatrix(5);
}

CLinearAlgebraHelper::~CLinearAlgebraHelper()
{
    if (NULL != m_pDeviceIntBuffer)
    {
        checkCudaErrors(cudaFree(m_pDeviceIntBuffer));
    }
    if (NULL != m_pDeviceFloatBuffer)
    {
        checkCudaErrors(cudaFree(m_pDeviceFloatBuffer));
    }
    if (NULL != m_pDeviceComplexBuffer1)
    {
        checkCudaErrors(cudaFree(m_pDeviceComplexBuffer1));
    }
    if (NULL != m_pDeviceComplexBuffer2)
    {
        checkCudaErrors(cudaFree(m_pDeviceComplexBuffer2));
    }

    for (int i = 0; i < m_lstTmpMatrix.Num(); ++i)
    {
        checkCudaErrors(cudaFree(m_lstTmpMatrix[i].m_pMatrix));
    }
}

void CLinearAlgebraHelper::TestSmallMatrix()
{
    const INT testDim1 = 30;
    const INT testDim2 = 15;
    CLinearAlgebraHelper* pHelper = new CLinearAlgebraHelper(testDim1);
    appGeneral(_T("\n(* ============= test QR Factorization ============== *)\n"));
    CLGComplex mij[testDim1 * testDim2];
    for (INT x = 0; x < testDim1; ++x)
    {
        for (INT y = 0; y < testDim2; ++y)
        {
            mij[x * testDim2 + y].x = (rand() % 101 - 50) / 50.0f;
            mij[x * testDim2 + y].y = (rand() % 101 - 50) / 50.0f;
        }
    }

    CLGComplex* deviceM = NULL;
    CLGComplex* deviceQ = NULL;
    CLGComplex* deviceR = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceM, sizeof(CLGComplex) * testDim1 * testDim2));
    checkCudaErrors(cudaMalloc((void**)&deviceQ, sizeof(CLGComplex) * testDim1 * testDim2));
    checkCudaErrors(cudaMalloc((void**)&deviceR, sizeof(CLGComplex) * testDim2 * testDim2));

    checkCudaErrors(cudaMemcpy(deviceM, mij, sizeof(CLGComplex) * testDim1 * testDim2, cudaMemcpyHostToDevice));
    appGeneral(_T("\nt=\n"));
    PrintMatrix(mij, testDim1, testDim2);
    pHelper->ThinQRFactorization(deviceQ, deviceR, deviceM, testDim1, testDim2);
    
    checkCudaErrors(cudaMemcpy(mij, deviceQ, sizeof(CLGComplex) * testDim1 * testDim2, cudaMemcpyDeviceToHost));
    appGeneral(_T("\nq=\n"));
    PrintMatrix(mij, testDim1, testDim2);
    checkCudaErrors(cudaMemcpy(mij, deviceR, sizeof(CLGComplex) * testDim2 * testDim2, cudaMemcpyDeviceToHost));
    appGeneral(_T("\nr=\n"));
    PrintMatrix(mij, testDim2, testDim2);

    appGeneral(_T("\n\n(* ============= test Eigen Vector ============== *)\n"));
    CLGComplex hij[testDim1 * testDim1];
    for (INT x = 0; x < testDim1; ++x)
    {
        for (INT y = 0; y < testDim1; ++y)
        {
            hij[x * testDim1 + y].x = (rand() % 101 - 50) / 50.0f;
            hij[x * testDim1 + y].y = (rand() % 101 - 50) / 50.0f;
        }
    }

    CLGComplex* deviceH = NULL;
    CLGComplex* deviceE = NULL;
    CLGComplex* deviceV = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceH, sizeof(CLGComplex) * testDim1 * testDim1));
    checkCudaErrors(cudaMalloc((void**)&deviceE, sizeof(CLGComplex) * testDim2));
    checkCudaErrors(cudaMalloc((void**)&deviceV, sizeof(CLGComplex) * testDim1 * testDim2));

    checkCudaErrors(cudaMemcpy(deviceH, hij, sizeof(CLGComplex) * testDim1 * testDim1, cudaMemcpyHostToDevice));
    appGeneral(_T("\nh=\n"));
    PrintMatrix(hij, testDim1, testDim1);
    pHelper->EigenValueProblem(deviceH, deviceE, deviceV, testDim1, testDim2);

    checkCudaErrors(cudaMemcpy(hij, deviceE, sizeof(CLGComplex) * testDim2, cudaMemcpyDeviceToHost));
    appGeneral(_T("\n(* Smallest Eigen Values *)\n ev1="));
    PrintMatrix(hij, 1, testDim2);

    checkCudaErrors(cudaMemcpy(hij, deviceV, sizeof(CLGComplex) * testDim1 * testDim2, cudaMemcpyDeviceToHost));
    appGeneral(_T("\n(* Corresponding Eigen Vectors *)\n vv1="));
    PrintMatrix(hij, testDim2, testDim1);

    appGeneral(_T("\n\n(* ============= test Generalized Eigen Vector ==============*)\n"));

    CLGComplex aij[testDim1 * testDim1];
    CLGComplex bij[testDim1 * testDim1];
    for (INT x = 0; x < testDim1; ++x)
    {
        for (INT y = 0; y < testDim1; ++y)
        {
            aij[x * testDim1 + y].x = (rand() % 101 - 50) / 50.0f;
            aij[x * testDim1 + y].y = (rand() % 101 - 50) / 50.0f;
            bij[x * testDim1 + y].x = (rand() % 101 - 50) / 50.0f;
            bij[x * testDim1 + y].y = (rand() % 101 - 50) / 50.0f;
        }
    }

    CLGComplex* deviceA = NULL;
    CLGComplex* deviceB = NULL;
    checkCudaErrors(cudaMalloc((void**)&deviceA, sizeof(CLGComplex) * testDim1 * testDim1));
    checkCudaErrors(cudaMalloc((void**)&deviceB, sizeof(CLGComplex) * testDim1 * testDim1));
    checkCudaErrors(cudaMemcpy(deviceA, aij, sizeof(CLGComplex) * testDim1 * testDim1, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(deviceB, bij, sizeof(CLGComplex) * testDim1 * testDim1, cudaMemcpyHostToDevice));

    appGeneral(_T("\na=\n"));
    PrintMatrix(aij, testDim1, testDim1);

    appGeneral(_T("\nb=\n"));
    PrintMatrix(bij, testDim1, testDim1);

    pHelper->GeneralizedEigenValueProblem(deviceA, deviceB, deviceE, deviceV, testDim1, testDim2);

    checkCudaErrors(cudaMemcpy(hij, deviceE, sizeof(CLGComplex) * testDim2, cudaMemcpyDeviceToHost));
    appGeneral(_T("\n(* Smallest Eigen Values*) \n ev2="));
    PrintMatrix(hij, 1, testDim2);

    checkCudaErrors(cudaMemcpy(hij, deviceV, sizeof(CLGComplex) * testDim1 * testDim2, cudaMemcpyDeviceToHost));
    appGeneral(_T("\n(* Corresponding Eigen Vectors*)\n vv2="));
    PrintMatrix(hij, testDim2, testDim1);

    appGeneral(_T("\n\n(*============= Please copy those results to Mathematica to check, all should be nearly zero ============ *)\n"));

    appGeneral(_T("\n(*============= test QR ============ *)\n"));
    appGeneral(_T("Max[t - q.r // Flatten // Abs]\n"));
    appGeneral(_T("Tr[Conjugate[Transpose[q]].q] / %d - 1\n"), testDim2);
    appGeneral(_T("\n(*============= test EV ============ *)\n"));

    appGeneral(_T("Eigenvalues[h][[Table[%d - i, {i, 1, %d}]]] - ev1[[1]] // Abs // Max\n"), testDim1 + 1, testDim2);

    for (INT i = 0; i < testDim2; ++i)
    {
        appGeneral(_T("Max[Eigensystem[h][[2]][[%d]]/Eigensystem[h][[2]][[%d]][[%d]] - vv1[[%d]] / vv1[[%d]][[%d]] //Abs]\n"), 
            testDim1 - i, testDim1 - i, testDim1, i + 1, i + 1, testDim1);
    }

    appGeneral(_T("Eigenvalues[Inverse[b].a][[Table[%d - i, {i, 1, %d}]]] - ev2[[1]] // Abs // Max\n"), testDim1 + 1, testDim2);

    for (INT i = 0; i < testDim2; ++i)
    {
        appGeneral(_T("Max[Eigensystem[Inverse[b].a][[2]][[%d]]/Eigensystem[Inverse[b].a][[2]][[%d]][[%d]] - vv2[[%d]] / vv2[[%d]][[%d]] //Abs]\n"), 
            testDim1 - i, testDim1 - i, testDim1, i + 1, i + 1, testDim1);
    }

    checkCudaErrors(cudaFree(deviceM));
    checkCudaErrors(cudaFree(deviceQ));
    checkCudaErrors(cudaFree(deviceR));
    checkCudaErrors(cudaFree(deviceH));
    checkCudaErrors(cudaFree(deviceE));
    checkCudaErrors(cudaFree(deviceV));
    checkCudaErrors(cudaFree(deviceA));
    checkCudaErrors(cudaFree(deviceB));
}

#pragma region Common

#pragma region Initial

__global__ void _CLG_LAUNCH_BOUND
_kernelInitialZero(CLGComplex* R, UINT dy)
{
    R[threadIdx.x * dy + threadIdx.y] = _make_cuComplex(F(0.0), F(0.0));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelInitialOne(CLGComplex* R, UINT dy)
{
    UINT i = threadIdx.x;
    UINT j = threadIdx.y;
    if (i == j)
    {
        R[threadIdx.x * dy + threadIdx.y] = _make_cuComplex(F(1.0), F(0.0));
    }
    else
    {
        R[threadIdx.x * dy + threadIdx.y] = _make_cuComplex(F(0.0), F(0.0));
    }
}

void CLinearAlgebraHelper::InitialZero(CLGComplex* deviceMatrix, UINT dx, UINT dy)
{
    dim3 block(1, 1, 1);
    dim3 thread(dx, dy, 1);
    _kernelInitialZero << <block, thread >> > (deviceMatrix, dy);
}

void CLinearAlgebraHelper::InitialOne(CLGComplex* deviceMatrix, UINT dx)
{
    dim3 block(1, 1, 1);
    dim3 thread(dx, dx, 1);
    _kernelInitialOne << <block, thread >> > (deviceMatrix, dx);
}

#pragma endregion

#pragma region Multiply

/**
* This is tested to be faster than for, even with launch bound
* If Left = X*Y, Right = Y*Z
* Res = X*Z, block(Y,1,1) thread(X,Z,1)
* leftDim = Y, midDim = Z
*/
__global__ void _kernelSmallMatrixMult_NN(CLGComplex* res,
    const CLGComplex* __restrict__ left,
    const CLGComplex* __restrict__ right,
    UINT leftDim, UINT midDim)
{
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;
    UINT n = blockIdx.x;

    //left is dx x n, right is n x dy matrix
    CLGComplex toAdd = _cuCmulf(left[x * leftDim + n], right[n * midDim + y]);

    atomicAdd(&res[x * midDim + y].x, toAdd.x);
    atomicAdd(&res[x * midDim + y].y, toAdd.y);
}

/**
* Left * Right^+
* If Left = X*Y, Right = Z*Y (Right^+ = Y*Z)
* Res = X*Z, block(Y,1,1) thread(X,Z,1)
* leftDim = Y, midDim = Z
*/
__global__ void _kernelSmallMatrixMult_ND(CLGComplex* res,
    const CLGComplex* __restrict__ left,
    const CLGComplex* __restrict__ right,
    UINT leftDim, UINT midDim)
{
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;
    UINT n = blockIdx.x;

    //left is dx x n, right is n x dy matrix
    //left = dx dy
    //mid = dx dz
    //right = dz dy (right dagger = dy dz)
    //n->0->dy
    //x->0->dx
    //y->0->dz
    //leftDim = dy
    //rightDim = dy
    //midDim = dz
    CLGComplex toAdd = _cuCmulf(left[x * leftDim + n], _cuConjf(right[y * leftDim + n]));

    atomicAdd(&res[x * midDim + y].x, toAdd.x);
    atomicAdd(&res[x * midDim + y].y, toAdd.y);
}

/**
* Left^+ * Right
* If Left = Y*X (Left^+ = X*Y), Right = Y*Z
* Res = X*Z, block(Y,1,1) thread(X,Z,1)
* leftDim = X, midDim = Z
*/
__global__ void _kernelSmallMatrixMult_DN(CLGComplex* res,
    const CLGComplex* __restrict__ left,
    const CLGComplex* __restrict__ right,
    UINT leftDim, UINT midDim)
{
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;
    UINT n = blockIdx.x;

    CLGComplex toAdd = _cuCmulf(_cuConjf(left[n * leftDim + x]), right[n * midDim + y]);

    atomicAdd(&res[x * midDim + y].x, toAdd.x);
    atomicAdd(&res[x * midDim + y].y, toAdd.y);
}


void CLinearAlgebraHelper::SmallMatrixMult(
    CLGComplex * deviceRes, 
    const CLGComplex* left, 
    const CLGComplex* right, 
    UINT dLeft, UINT dMid, UINT dRight, 
    UBOOL bLeftDagger, UBOOL bRightDagger)
{
    InitialZero(deviceRes, dLeft, dRight);
    dim3 block(dMid, 1, 1);
    dim3 thread(dLeft, dRight, 1);
    if (bLeftDagger)
    {
        _kernelSmallMatrixMult_DN << <block, thread >> >(deviceRes, left, right, dLeft, dRight);
    }
    else if (bRightDagger)
    {
        _kernelSmallMatrixMult_ND<<<block, thread >>>(deviceRes, left, right, dMid, dRight);
    }
    else
    {
        _kernelSmallMatrixMult_NN << <block, thread >> > (deviceRes, left, right, dMid, dRight);
    }
}

#pragma endregion

#pragma region Block Multply

/**
* left = 1 0 0
*        0 U 0
*        0 0 1
*
* right = A1 A2 A3
*         B1 B2 B3
*         C1 C2 C3
*
* res = A1  A2  A3
*       UB1 UB2 UB3
*       C1  C2  C3
*
* if U is dy x dy
* block  = dy, 1, 1
* thread = dx, dx, 1
* Assume res is zeroed
* Y Dir ----->
*/
__global__ void _kernelMatrixBlockMult_LNN(CLGComplex* res, const CLGComplex* __restrict__ left, const CLGComplex* __restrict__ right, UINT iStart, UINT iEnd, UINT dm)
{
    UINT n = blockIdx.x;
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    if (x >= iStart && x < iEnd)
    {
        UINT mid = n + iStart;
        CLGComplex toAdd = _cuCmulf(left[x * dm + mid], right[mid * dm + y]);
        atomicAdd(&res[x * dm + y].x, toAdd.x);
        atomicAdd(&res[x * dm + y].y, toAdd.y);
    }
    else
    {
        if (0 == n)
        {
            res[x * dm + y] = right[x * dm + y];
        }
    }
}

__global__ void _kernelMatrixBlockMult_LDN(CLGComplex* res, const CLGComplex* __restrict__ left, const CLGComplex* __restrict__ right, UINT iStart, UINT iEnd, UINT dm)
{
    UINT n = blockIdx.x;
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    if (x >= iStart && x < iEnd)
    {
        UINT mid = n + iStart;
        CLGComplex toAdd = _cuCmulf(_cuConjf(left[mid * dm + x]), right[mid * dm + y]);
        atomicAdd(&res[x * dm + y].x, toAdd.x);
        atomicAdd(&res[x * dm + y].y, toAdd.y);
    }
    else
    {
        if (0 == n)
        {
            res[x * dm + y] = right[x * dm + y];
        }
    }
}

__global__ void _kernelMatrixBlockMult_LND(CLGComplex* res, const CLGComplex* __restrict__ left, const CLGComplex* __restrict__ right, UINT iStart, UINT iEnd, UINT dm)
{
    UINT n = blockIdx.x;
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    if (x >= iStart && x < iEnd)
    {
        UINT mid = n + iStart;
        CLGComplex toAdd = _cuCmulf(left[x * dm + mid], _cuConjf(right[y * dm + mid]));
        atomicAdd(&res[x * dm + y].x, toAdd.x);
        atomicAdd(&res[x * dm + y].y, toAdd.y);
    }
    else
    {
        if (0 == n)
        {
            res[x * dm + y] = _cuConjf(right[y * dm + x]);
        }
    }
}

/**
* left = A1 A2 A3
*        B1 B2 B3
*        C1 C2 C3
*
* right = 1 0 0
*         0 U 0
*         0 0 1
*
* res = A1 A2U A3
*       B1 B2U B3
*       C1 C2U C3
*
* if U is dy x dy
* block  = dy, 1, 1
* thread = dx, dx, 1
* Assume res is zeroed
* Y Dir ----->
*/
__global__ void _kernelMatrixBlockMult_RNN(CLGComplex* res, const CLGComplex* __restrict__ left, const CLGComplex* __restrict__ right, UINT iStart, UINT iEnd, UINT dm)
{
    UINT n = blockIdx.x;
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    if (y >= iStart && y < iEnd)
    {
        UINT mid = n + iStart;
        CLGComplex toAdd = _cuCmulf(left[x * dm + mid], right[mid * dm + y]);
        atomicAdd(&res[x * dm + y].x, toAdd.x);
        atomicAdd(&res[x * dm + y].y, toAdd.y);
    }
    else
    {
        if (0 == n)
        {
            res[x * dm + y] = left[x * dm + y];
        }
    }
}

__global__ void _kernelMatrixBlockMult_RDN(CLGComplex* res, const CLGComplex* __restrict__ left, const CLGComplex* __restrict__ right, UINT iStart, UINT iEnd, UINT dm)
{
    UINT n = blockIdx.x;
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    if (y >= iStart && y < iEnd)
    {
        UINT mid = n + iStart;
        CLGComplex toAdd = _cuCmulf(_cuConjf(left[mid * dm + x]), right[mid * dm + y]);
        atomicAdd(&res[x * dm + y].x, toAdd.x);
        atomicAdd(&res[x * dm + y].y, toAdd.y);
    }
    else
    {
        if (0 == n)
        {
            res[x * dm + y] = _cuConjf(left[y * dm + x]);
        }
    }
}

__global__ void _kernelMatrixBlockMult_RND(CLGComplex* res, const CLGComplex* __restrict__ left, const CLGComplex* __restrict__ right, UINT iStart, UINT iEnd, UINT dm)
{
    UINT n = blockIdx.x;
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    if (y >= iStart && y < iEnd)
    {
        UINT mid = n + iStart;
        CLGComplex toAdd = _cuCmulf(left[x * dm + mid], _cuConjf(right[y * dm + mid]));
        atomicAdd(&res[x * dm + y].x, toAdd.x);
        atomicAdd(&res[x * dm + y].y, toAdd.y);
    }
    else
    {
        if (0 == n)
        {
            res[x * dm + y] = left[x * dm + y];
        }
    }
}

void CLinearAlgebraHelper::BlockMatrixMult(
    CLGComplex * deviceRes,
    const CLGComplex* left,
    const CLGComplex* right,
    UINT dDim, UINT uiStart, UINT uiEnd,
    UBOOL bLeft, UBOOL bLeftDagger, UBOOL bRightDagger)
{
    InitialZero(deviceRes, dDim, dDim);
    dim3 block(uiEnd - uiStart, 1, 1);
    dim3 thread(dDim, dDim, 1);

    if (bLeft)
    {
        if (bLeftDagger)
        {
            _kernelMatrixBlockMult_LDN << <block, thread >> >(deviceRes, left, right, uiStart, uiEnd, dDim);
        }
        else if (bRightDagger)
        {
            _kernelMatrixBlockMult_LND << <block, thread >> >(deviceRes, left, right, uiStart, uiEnd, dDim);
        }
        else
        {
            _kernelMatrixBlockMult_LNN << <block, thread >> > (deviceRes, left, right, uiStart, uiEnd, dDim);
        }
    }
    else
    {
        if (bLeftDagger)
        {
            _kernelMatrixBlockMult_RDN << <block, thread >> >(deviceRes, left, right, uiStart, uiEnd, dDim);
        }
        else if (bRightDagger)
        {
            _kernelMatrixBlockMult_RND << <block, thread >> >(deviceRes, left, right, uiStart, uiEnd, dDim);
        }
        else
        {
            _kernelMatrixBlockMult_RNN << <block, thread >> > (deviceRes, left, right, uiStart, uiEnd, dDim);
        }
    }
}

#pragma endregion

#pragma region Add, Minus etc

/**
* M=M+cI, c is device buffer
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelMatrixAddConstant(CLGComplex* m, CLGComplex* c, UINT dy)
{
    UINT i = threadIdx.x;
    m[i * dy + i] = _cuCaddf(m[i * dy + i], c[0]);
}

#pragma endregion

#pragma region Block Copy

/**
* thread.xy = lx,ly
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCopyMatrixXY(
    CLGComplex* mtr,
    const CLGComplex* __restrict__ orignal,
    UINT newdy, UINT olddy)
{
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    mtr[x * newdy + y] = orignal[x * olddy + y];
}

void CLinearAlgebraHelper::BlockCopy(CLGComplex* deviceDest, const CLGComplex* deviceSrc,
    UINT lengthX, UINT lengthY, UINT dimDest, UINT dimSrc)
{
    dim3 block(1, 1, 1);
    dim3 thread(lengthX, lengthY, 1);
    _kernelCopyMatrixXY << <block, thread >> > (deviceDest, deviceSrc, dimDest, dimSrc);
}

#pragma endregion

void CLinearAlgebraHelper::PrintMatrix(const CLGComplex* mtr, UINT dx, UINT dy)
{
    appGeneral("\n{");
    for (UINT i = 0; i < dx; ++i)
    {
        for (UINT j = 0; j < dy; ++j)
        {
            appGeneral("%s%1.5f %s %1.5f I%s ",
                0 == j ? "{" : "",
                mtr[i * dy + j].x,
                mtr[i * dy + j].y < F(0.0) ? "" : "+",
                mtr[i * dy + j].y,
                dy - 1 == j ? "}" : ",");
        }
        if (i == dx - 1)
        {
            appGeneral("};\n");
        }
        else
        {
            appGeneral(",\n");
        }
    }
}

#pragma endregion

#pragma region QR Decomposition

__global__ void _CLG_LAUNCH_BOUND
_kernelOneStepHouseHolderQR(
    CLGComplex* Q,
    const CLGComplex* __restrict__ R,
    UINT i, UINT dy)
{
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    __shared__ Real length;
    __shared__ Real lengthu;
    __shared__ CLGComplex u[CLinearAlgebraHelper::_kMaxSmallDim];

    if (0 == x && 0 == y)
    {
        length = F(0.0);
        lengthu = F(0.0);
    }

    __syncthreads();

    if (0 == y && x >= i)
    {
        atomicAdd(&length, R[x * dy + i].x * R[x * dy + i].x + R[x * dy + i].y * R[x * dy + i].y);
    }

    __syncthreads();

    if (0 == x && 0 == y)
    {
        length = _sqrt(length);
    }

    __syncthreads();

    if (0 == y && x >= i)
    {
        u[x] = R[x * dy + i];

        if (x == i)
        {
            Real fuxL = __div(F(1.0), _cuCabsf(u[x]));
            Real fCos = u[x].x * fuxL;
            Real fSin = u[x].y * fuxL;
            u[x] = _cuCaddf(u[x], _make_cuComplex(length * fCos, length * fSin));
        }
        atomicAdd(&lengthu, u[x].x * u[x].x + u[x].y * u[x].y);
    }

    __syncthreads();

    if (0 == x && 0 == y)
    {
        lengthu = __div(F(1.0), _sqrt(lengthu * F(0.5)));
    }

    __syncthreads();

    if (0 == y && x >= i)
    {
        u[x].x = u[x].x * lengthu;
        u[x].y = u[x].y * lengthu;
    }

    __syncthreads();

    //uk = A[i + 1->n, i] - |A[i+1]|
    if (x < i || y < i)
    {
        if (x == y)
        {
            Q[x * dy + y] = _make_cuComplex(F(1.0), F(0.0));
        }
        else
        {
            Q[x * dy + y] = _make_cuComplex(F(0.0), F(0.0));
        }
    }
    else
    {
        Q[x * dy + y] = _cuCmulf(_cuConjf(u[y]), u[x]);
        Q[x * dy + y].x = -Q[x * dy + y].x;
        Q[x * dy + y].y = -Q[x * dy + y].y;
        if (x == y)
        {
            Q[x * dy + y].x = Q[x * dy + y].x + F(1.0);
        }
    }
}

void CLinearAlgebraHelper::QRFactorization(
    CLGComplex* Q, 
    CLGComplex* R, 
    const CLGComplex* T, 
    UINT dy)
{
    if (dy > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }
    dim3 block(1, 1, 1);
    dim3 thread1(dy, dy, 1);
    STmpMatrix sTmpQ = GetTmpMatrix();
    CLGComplex* tmpQ = sTmpQ.m_pMatrix;
    STmpMatrix sTmpM = GetTmpMatrix();
    CLGComplex* tmpM = sTmpM.m_pMatrix;
    checkCudaErrors(cudaMemcpy(R, T, sizeof(CLGComplex) * dy * dy, cudaMemcpyDeviceToDevice));
    InitialOne(Q, dy);
    for (UINT i = 0; i < dy - 1; ++i)
    {
        _kernelOneStepHouseHolderQR << <block, thread1 >> > (tmpQ, R, i, dy);

        //left is block
        BlockMatrixMult(tmpM, tmpQ, R, dy, i, dy, TRUE, FALSE, FALSE);

        checkCudaErrors(cudaMemcpy(R, tmpM, sizeof(CLGComplex) * dy * dy, cudaMemcpyDeviceToDevice));

        //right is block and right is dagger
        BlockMatrixMult(tmpM, Q, tmpQ, dy, i, dy, FALSE, FALSE, TRUE);

        checkCudaErrors(cudaMemcpy(Q, tmpM, sizeof(CLGComplex) * dy * dy, cudaMemcpyDeviceToDevice));
    }
    sTmpQ.Free();
    sTmpM.Free();
}

void CLinearAlgebraHelper::ThinQRFactorization(
    CLGComplex* Q,
    CLGComplex* R,
    const CLGComplex* T,
    UINT dx, UINT dy)
{
    if (dy > m_uiDim || dx > m_uiDim || dy > dx)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d, or it is not thin!\n"), m_uiDim);
        return;
    }

    dim3 thread2(dx, dx, 1);

    STmpMatrix sTmpR = GetTmpMatrix();
    CLGComplex* tmpR = sTmpR.m_pMatrix;
    STmpMatrix sTmpQ = GetTmpMatrix();
    CLGComplex* tmpQ = sTmpQ.m_pMatrix;
    STmpMatrix sTmpQ2 = GetTmpMatrix();
    CLGComplex* tmpQ2 = sTmpQ2.m_pMatrix;
    STmpMatrix sTmpM = GetTmpMatrix();
    CLGComplex* tmpM = sTmpM.m_pMatrix;

    InitialZero(tmpR, dx, dx);
    BlockCopy(tmpR, T, dx, dy, dx, dy);
    InitialOne(tmpQ2, dx);
    dim3 block(1, 1, 1);
    for (UINT i = 0; i < dx - 1; ++i)
    {
        _kernelOneStepHouseHolderQR << <block, thread2 >> > (tmpQ, tmpR, i, dx);

        //m = q x r
        BlockMatrixMult(tmpM, tmpQ, tmpR, dx, i, dx, TRUE, FALSE, FALSE);

        //r = m
        checkCudaErrors(cudaMemcpy(tmpR, tmpM, sizeof(CLGComplex) * dx * dx, cudaMemcpyDeviceToDevice));

        //m = q2 * q
        BlockMatrixMult(tmpM, tmpQ2, tmpQ, dx, i, dx, FALSE, FALSE, TRUE);

        //q2 = m
        checkCudaErrors(cudaMemcpy(tmpQ2, tmpM, sizeof(CLGComplex) * dx * dx, cudaMemcpyDeviceToDevice));
    }
    
    //r=r, q=q2
    BlockCopy(R, tmpR, dy, dy, dy, dx);
    BlockCopy(Q, tmpQ2, dx, dy, dy, dx);
    sTmpQ.Free();
    sTmpM.Free();
    sTmpQ2.Free();
    sTmpR.Free();
}

#pragma endregion

#pragma region QR Hensenberg

__global__ void _CLG_LAUNCH_BOUND
_kernelOneStepHouseHolder(
    CLGComplex* U, 
    const CLGComplex* __restrict__ A, 
    UINT i, UINT dx)
{
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    __shared__ Real length;
    __shared__ Real lengthu;
    __shared__ CLGComplex u[CLinearAlgebraHelper::_kMaxSmallDim];
    if (0 == x && 0 == y)
    {
        length = F(0.0);
        lengthu = F(0.0);
    }

    __syncthreads();

    if (0 == y && x > i)
    {
        atomicAdd(&length, A[x * dx + i].x * A[x * dx + i].x + A[x * dx + i].y * A[x * dx + i].y);
    }

    __syncthreads();

    if (0 == x && 0 == y)
    {
        length = _sqrt(length);
    }

    __syncthreads();

    if (0 == y && x > i)
    {
        u[x] = A[x * dx + i];

        if (x == i + 1)
        {
            Real fuxL = __div(F(1.0), _sqrt(_cuCabsf(u[x])));
            Real fCos = u[x].x * fuxL;
            Real fSin = u[x].y * fuxL;
            u[x] = _cuCaddf(u[x], _make_cuComplex(length * fCos, length * fSin));
        }
        atomicAdd(&lengthu, u[x].x * u[x].x + u[x].y * u[x].y);
    }

    __syncthreads();

    if (0 == x && 0 == y)
    {
        lengthu = __div(F(1.0), _sqrt(lengthu * F(0.5)));
    }

    __syncthreads();

    if (0 == y && x > i)
    {
        u[x].x = u[x].x * lengthu;
        u[x].y = u[x].y * lengthu;
    }

    __syncthreads();

    //uk = A[i + 1->n, i] - |A[i+1]|
    if (x <= i || y <= i)
    {
        if (x == y)
        {
            U[x * dx + y] = _make_cuComplex(F(1.0), F(0.0));
        }
        else
        {
            U[x * dx + y] = _make_cuComplex(F(0.0), F(0.0));
        }
    }
    else
    {
        U[x * dx + y] = _cuCmulf(_cuConjf(u[y]), u[x]);
        U[x * dx + y].x = -U[x * dx + y].x;
        U[x * dx + y].y = -U[x * dx + y].y;
        if (x == y)
        {
            U[x * dx + y].x = U[x * dx + y].x + F(1.0);
        }
    }
}

void CLinearAlgebraHelper::Henssenberg(CLGComplex* T, UINT dx)
{
    if (dx > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }
    dim3 block(1, 1, 1);
    dim3 thread1(dx, dx, 1);
    STmpMatrix sTmpU = GetTmpMatrix();
    CLGComplex* tmpU = sTmpU.m_pMatrix;
    STmpMatrix sTmpM = GetTmpMatrix();
    CLGComplex* tmpM = sTmpM.m_pMatrix;

    for (UINT i = 0; i < dx - 2; ++i)
    {
        _kernelOneStepHouseHolder << <block, thread1 >> > (tmpU, T, i, dx);
        BlockMatrixMult(tmpM, tmpU, T, dx, i + 1, dx, TRUE, FALSE, FALSE);
        BlockMatrixMult(T, tmpM, tmpU, dx, i + 1, dx, FALSE, FALSE, TRUE);
    }
    sTmpU.Free();
    sTmpM.Free();
}

#pragma endregion

#pragma region Back Shift

__global__ void _CLG_LAUNCH_BOUND
_kernelOneLineReduceBS(
    CLGComplex* y, 
    const CLGComplex* __restrict__ R, 
    UINT i, UINT dk, UINT dx)
{
    UINT j = threadIdx.x + i + 1; //j=i+1 to dx
    UINT n = threadIdx.y;

    if (j < dx)
    {
        CLGComplex toAdd = _cuCmulf(R[i * dx + j], y[j * dk + n]);
        atomicAdd(&y[i * dk + n].x, -toAdd.x);
        atomicAdd(&y[i * dk + n].y, -toAdd.y);
    }

    __syncthreads();

    //i + 1 == j means threadIdx = 0
    if (i + 1 == j)
    {
        y[i * dk + n] = _cuCdivf(y[i * dk + n], R[i * dx + i]);
    }
}

void CLinearAlgebraHelper::SolveY(CLGComplex* deviceY, const CLGComplex* deviceR, UINT dk, UINT dx)
{
    dim3 block(1, 1, 1);
    for (INT i = dx - 1; i >= 0; --i)
    {
        if (i == static_cast<INT>(dx - 1))
        {
            dim3 thread(1, dk, 1);
            _kernelOneLineReduceBS << <block, thread >> > (deviceY, deviceR, i, dk, dx);
        }
        else
        {
            dim3 thread(dx - i - 1, dk, 1);
            _kernelOneLineReduceBS << <block, thread >> > (deviceY, deviceR, i, dk, dx);
        }
    }
}


#pragma endregion

#pragma region Shift QR Iteration

__global__ void _CLG_LAUNCH_BOUND
_kernelWilkinsonShift(CLGComplex* m, CLGComplex* c, UINT dim)
{
    UINT i = threadIdx.x;
    if (0 == i)
    {
        //d
        c[0] = m[dim * dim - 1];

        //bc
        CLGComplex omega = _cuCmulf(m[dim * dim - dim - 1], m[dim * dim - 2]);

        Real fOmegaSq = __cuCabsSqf(omega);
        if (fOmegaSq > F(0.0000001))
        {
            //(d-a)/2
            CLGComplex xi = _make_cuComplex(
                F(0.5) * (c[0].x - m[dim * dim - dim - 2].x),
                F(0.5) * (c[0].y - m[dim * dim - dim - 2].y));

            //sqrt(((d-a)/2)^2 + bc)
            CLGComplex eta = __cuCsqrtf(_cuCaddf(_cuCmulf(xi, xi), omega));
            if (xi.x * eta.x + xi.y * eta.y < F(0.0))
            {
                c[0] = _cuCsubf(c[0], _cuCdivf(omega, _cuCsubf(eta, xi)));
            }
            else
            {
                c[0] = _cuCaddf(c[0], _cuCdivf(omega, _cuCaddf(eta, xi)));
            }
        }
    }

    __syncthreads();

    m[i * dim + i] = _cuCsubf(m[i * dim + i], c[0]);
}

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelCheckMatrix(CLGComplex* mtr, INT* decomp, UINT dx, Real fCrit)
{
    decomp[0] = dx;
    for (INT i = dx - 2; i >= 0; --i)
    {
        if (__cuCabsSqf(mtr[(i + 1) * dx + i]) < 
            fCrit * (__cuCabsSqf(mtr[(i + 1) * dx + i + 1]) + __cuCabsSqf(mtr[i * dx + i])))
        {
            mtr[(i + 1) * dx + i].x = F(0.0);
            mtr[(i + 1) * dx + i].y = F(0.0);

            if (decomp[0] == i + 2)
            {
                decomp[0] = i + 1;
            }
        }
    }
}

void CLinearAlgebraHelper::QRIterate(CLGComplex* T, UINT dx, Real fCrit, UINT iCrit)
{
    if (dx > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpQ = GetTmpMatrix();
    CLGComplex* tmpQ = sTmpQ.m_pMatrix;
    STmpMatrix sTmpR = GetTmpMatrix();
    CLGComplex* tmpR = sTmpR.m_pMatrix;
    STmpMatrix sTmpT = GetTmpMatrix();
    CLGComplex* tmpT = sTmpT.m_pMatrix;

    CLGComplex* tmpShift = m_pDeviceComplexBuffer1;
    INT* tmpDecomp = m_pDeviceIntBuffer;
    
    dim3 block(1, 1, 1);
    dim3 thread1(dx, dx, 1);
    INT endindex[1];

    checkCudaErrors(cudaMemcpy(tmpT, T, sizeof(CLGComplex) * dx * dx, cudaMemcpyDeviceToDevice));

    UINT iLastDim = dx;
    for (UINT i = 0; i < iCrit; ++i)
    {
        //find decomp
        _kernelCheckMatrix << <1, 1 >> > (tmpT, tmpDecomp, iLastDim, fCrit);
        checkCudaErrors(cudaMemcpy(endindex, tmpDecomp, sizeof(INT), cudaMemcpyDeviceToHost));

        if (endindex[0] < static_cast<INT>(iLastDim))
        {
            //copy matrix
            BlockCopy(T, tmpT, iLastDim, iLastDim, dx, iLastDim);
            if (1 == endindex[0])
            {
                //finished
                sTmpQ.Free();
                sTmpR.Free();
                sTmpT.Free();

                return;
            }

            iLastDim = endindex[0];
            //must do the copy, because [x * dim + y] the "dim" is changed.
            BlockCopy(tmpT, T, iLastDim, iLastDim, iLastDim, dx);
        }

        //shift
        //T = T - sigma I, tmpDeviceFloat[0] = sigma
        dim3 thread(iLastDim, 1, 1);
        _kernelWilkinsonShift << <block, thread >> > (tmpT, tmpShift, iLastDim);

        //QR decompose
        QRFactorization(tmpQ, tmpR, tmpT, iLastDim);

        //Update H
        //T = R Q + sigma I
        SmallMatrixMult(tmpT, tmpR, tmpQ, iLastDim, iLastDim, iLastDim, FALSE, FALSE);
        _kernelMatrixAddConstant << <block, thread >> > (tmpT, tmpShift, iLastDim);
    }

    sTmpQ.Free();
    sTmpR.Free();
    sTmpT.Free();
}

#pragma endregion

#pragma region Eigen Problem

__global__ void _CLG_LAUNCH_BOUND
_kernelSortEigenValues(const CLGComplex* __restrict__ R,
    CLGComplex* outV, Real* tmpF, INT* tmpO, UINT k, UINT dx)
{
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    if (0 == x)
    {
        tmpF[y] = R[y * dx + y].x * R[y * dx + y].x + R[y * dx + y].y * R[y * dx + y].y;
        tmpO[y] = 0;
    }

    __syncthreads();

    if (x != y)
    {
        if (tmpF[x] < tmpF[y])
        {
            atomicAdd(&tmpO[y], 1);
        }
    }

    __syncthreads();

    if (0 == x)
    {
        if (tmpO[y] < k)
        {
            outV[tmpO[y]] = R[y * dx + y];
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDaggerVector(CLGComplex* y, const CLGComplex* __restrict__ Q, UINT dx)
{
    UINT j = threadIdx.x;
    y[j] = cuConjf(Q[j * dx]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelInverseIterateShift(CLGComplex* A, const CLGComplex* __restrict__ outV, UINT k, UINT dx)
{
    UINT x = threadIdx.x;
    A[x * dx + x] = cuCsubf(A[x * dx + x], outV[k]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelNormVectors(CLGComplex* v, UINT dx)
{
    UINT x = threadIdx.y;
    UINT y = threadIdx.x;
    __shared__ Real fAmp[CLinearAlgebraHelper::_kMaxSmallDim];
    if (0 == y)
    {
        fAmp[x] = F(0.0);
    }

    __syncthreads();

    atomicAdd(&fAmp[x], v[x * dx + y].x * v[x * dx + y].x + v[x * dx + y].y * v[x * dx + y].y);

    __syncthreads();

    if (0 == y)
    {
        fAmp[x] = __div(F(1.0), _sqrt(fAmp[x]));
    }

    v[x * dx + y].x = v[x * dx + y].x * fAmp[x];
    v[x * dx + y].y = v[x * dx + y].y * fAmp[x];
}

__global__ void _CLG_LAUNCH_BOUND
_kernelErrorCheck(Real* outE, CLGComplex* v, const CLGComplex* __restrict__ A, UINT dx)
{
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    __shared__ Real length;
    __shared__ CLGComplex afterMult[CLinearAlgebraHelper::_kMaxSmallDim];

    if (0 == x && 0 == y)
    {
        length = F(0.0);
    }

    __syncthreads();

    if (0 == x)
    {
        atomicAdd(&length, v[y].x * v[y].x + v[y].y * v[y].y);
        afterMult[y] = _make_cuComplex(F(0.0), F(0.0));
    }

    __syncthreads();

    if (0 == x && 0 == y)
    {
        length = __div(F(1.0), _sqrt(length));
    }

    __syncthreads();

    if (0 == x)
    {
        v[y].x = v[y].x * length;
        v[y].y = v[y].y * length;
    }

    __syncthreads();

    cuComplex toAdd = cuCmulf(A[x * dx + y], v[y]);
    atomicAdd(&afterMult[x].x, toAdd.x);
    atomicAdd(&afterMult[x].y, toAdd.y);

    __syncthreads();

    if (0 == x)
    {
        atomicAdd(outE, afterMult[y].x * afterMult[y].x + afterMult[y].y * afterMult[y].y);
    }
}


void CLinearAlgebraHelper::EigenValueProblem(
    CLGComplex* H,
    CLGComplex* outEigenValue,
    CLGComplex* outEigenVector,
    UINT dm, UINT dk,
    Real fEigenCrit,
    UINT iMaxEigenIterate,
    Real fCrit,
    UINT iMaxIterate)
{
    if (dm > m_uiDim || dk > dm)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d! or required eigen vector number larger than dimension!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpH = GetTmpMatrix();
    CLGComplex* tmpH = sTmpH.m_pMatrix;

    CLGComplex* tmpVector = m_pDeviceComplexBuffer1;

    //preserve H for solve eigen vectors
    checkCudaErrors(cudaMemcpy(tmpH, H, sizeof(CLGComplex) * dm * dm, cudaMemcpyDeviceToDevice));

    Henssenberg(tmpH, dm);
    QRIterate(tmpH, dm, fCrit, iMaxIterate);

    STmpMatrix sTmpQ = GetTmpMatrix();
    CLGComplex* tmpQ = sTmpQ.m_pMatrix;
    STmpMatrix sTmpR = GetTmpMatrix();
    CLGComplex* tmpR = sTmpR.m_pMatrix;

    dim3 block(1, 1, 1);
    dim3 thread1(dm, dm, 1);
    dim3 thread2(dm, 1, 1);

    _kernelSortEigenValues << <block, thread1 >> > (tmpH, outEigenValue, m_pDeviceFloatBuffer, m_pDeviceIntBuffer, dk, dm);

    for (UINT i = 0; i < dk; ++i)
    {
        //Inverse Iterate
        checkCudaErrors(cudaMemcpy(tmpH, H, sizeof(CLGComplex) * dm * dm, cudaMemcpyDeviceToDevice));
        _kernelInverseIterateShift << <block, thread2 >> > (tmpH, outEigenValue, i, dm);

        QRFactorization(tmpQ, tmpR, tmpH, dm);

        //q=tmpM2, r=tmpM3
        _kernelDaggerVector << <block, thread2 >> > (tmpVector, tmpQ, dm);
        SolveY(tmpVector, tmpR, 1, dm);

        // Sometimes One Iteration is NOT enough!
        Real fErr[1];
        Real* tmpF = m_pDeviceFloatBuffer;
        STmpMatrix sTmpM = GetTmpMatrix();
        CLGComplex* tmpM = sTmpM.m_pMatrix;
        for (UINT j = 0; j < iMaxEigenIterate; ++j)
        {
            fErr[0] = F(0.0);
            checkCudaErrors(cudaMemcpy(tmpF, fErr, sizeof(Real), cudaMemcpyHostToDevice));

            _kernelErrorCheck << <block, thread1 >> > (tmpF, tmpVector, tmpH, dm);

            checkCudaErrors(cudaMemcpy(fErr, tmpF, sizeof(float), cudaMemcpyDeviceToHost));

            if (j == iMaxEigenIterate - 1)
            {
                appParanoiac(_T("(* Eigen vector(%d-iterate %d) error now = %1.12f *)\n"), i + 1, j, fErr[0]);
            }

            if (fErr[0] < fEigenCrit)
            {
                break;
            }
            SmallMatrixMult(tmpM, tmpQ, tmpVector, dm, dm, 1, TRUE, FALSE);
            SolveY(tmpM, tmpR, 1, dm);
            checkCudaErrors(cudaMemcpy(tmpVector, tmpM, sizeof(CLGComplex) * dm, cudaMemcpyDeviceToDevice));
        }

        checkCudaErrors(cudaMemcpy(outEigenVector + dm * i, tmpVector, sizeof(CLGComplex) * dm, cudaMemcpyDeviceToDevice));
    }

    //It is normalized in _kernelErrorCheck
    //dim3 thread3(dm, dk, 1);
    //_kernelNormVectors << <block, thread3 >> > (outEigenVector, dm);

    sTmpH.Free();
    sTmpQ.Free();
    sTmpR.Free();
}

#pragma endregion

#pragma region Generalized Eigen Problem

void CLinearAlgebraHelper::GeneralizedEigenValueProblem(
    CLGComplex* A,
    CLGComplex* B,
    CLGComplex* outEigenValue,
    CLGComplex* outEigenVector,
    UINT dm, UINT dk,
    Real fEigenCrit,
    UINT iMaxEigenIterate,
    Real fCrit,
    UINT iMaxIterate)
{
    if (dm > m_uiDim || dk > dm)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d! or required eigen vector number larger than dimension!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpQ = GetTmpMatrix();
    CLGComplex* tmpQ = sTmpQ.m_pMatrix;
    STmpMatrix sTmpR = GetTmpMatrix();
    CLGComplex* tmpR = sTmpR.m_pMatrix;

    QRFactorization(tmpQ, tmpR, B, dm);
    SmallMatrixMult(B, tmpQ, A, dm, dm, dm, TRUE, FALSE);
    SolveY(B, tmpR, dm, dm);

    sTmpQ.Free();
    sTmpR.Free();

    EigenValueProblem(B, outEigenValue, outEigenVector, dm, dk, fEigenCrit, iMaxEigenIterate, fCrit, iMaxIterate);
}

#pragma endregion

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
