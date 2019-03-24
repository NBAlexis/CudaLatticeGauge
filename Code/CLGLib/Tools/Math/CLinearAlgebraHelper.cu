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
    , m_pOneDeviceC(NULL)
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
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceTmpHouseHolder, sizeof(CLGComplex) * 9));
    checkCudaErrors(cudaMalloc((void**)&m_pOneDeviceC, sizeof(CLGComplex)));

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
    if (NULL != m_pDeviceTmpHouseHolder)
    {
        checkCudaErrors(cudaFree(m_pDeviceTmpHouseHolder));
    }
    if (NULL != m_pOneDeviceC)
    {
        checkCudaErrors(cudaFree(m_pOneDeviceC));
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
    
    CLGComplex mij[testDim1 * testDim2]; //Test thin QR, mij dagger mij
    CLGComplex hij[testDim1 * testDim1]; //Test EigenValue hij block mult iij, gev
    CLGComplex iij[testDim1 * testDim1]; //block multiply
    CLGComplex aij[testDim1 * testDim1]; //test gev
    CLGComplex bij[testDim1 * testDim1]; //test gev
    CLGComplex tij[testDim1 * testDim1]; //test triangular
    CLGComplex hesij[testDim1 * testDim1]; //test hessenberg eigen value
    CLGComplex heij[(testDim1 - 1) * testDim1];
    CLGComplex yij[testDim1];
    CLGComplex res1[testDim1 * testDim1];
    CLGComplex res2[testDim1 * testDim1];
    for (INT x = 0; x < testDim1; ++x)
    {
        for (INT y = 0; y < testDim2; ++y)
        {
            mij[x * testDim2 + y].x = (rand() % 101 - 50) / 50.0f;
            mij[x * testDim2 + y].y = (rand() % 101 - 50) / 50.0f;
        }

        for (INT y = 0; y < testDim1; ++y)
        {
            
            hij[x * testDim1 + y].x = GetRandomReal() * F(2.0) - F(1.0);
            hij[x * testDim1 + y].y = GetRandomReal() * F(2.0) - F(1.0);
            aij[x * testDim1 + y].x = GetRandomReal() * F(2.0) - F(1.0);
            aij[x * testDim1 + y].y = GetRandomReal() * F(2.0) - F(1.0);
            bij[x * testDim1 + y].x = GetRandomReal() * F(2.0) - F(1.0);
            bij[x * testDim1 + y].y = GetRandomReal() * F(2.0) - F(1.0);
            if (y >= x)
            {
                tij[x * testDim1 + y].x = GetRandomReal() * F(2.0) - F(1.0);
                tij[x * testDim1 + y].y = GetRandomReal() * F(2.0) - F(1.0);
            }
            else
            {
                tij[x * testDim1 + y] = _make_cuComplex(F(0.0), F(0.0));
            }

            if (y >= x - 1)
            {
                hesij[x * testDim1 + y].x = GetRandomReal() * F(2.0) - F(1.0);
                hesij[x * testDim1 + y].y = GetRandomReal() * F(2.0) - F(1.0);
            }
            else
            {
                hesij[x * testDim1 + y] = _make_cuComplex(F(0.0), F(0.0));
            }

            if (x < testDim2 && y < testDim2)
            {
                iij[x * testDim1 + y].x = GetRandomReal() * F(2.0) - F(1.0);
                iij[x * testDim1 + y].y = GetRandomReal() * F(2.0) - F(1.0);
            }
            else
            {
                if (y == x)
                {
                    iij[x * testDim1 + y] = _make_cuComplex(F(1.0), F(0.0));
                }
                else
                {
                    iij[x * testDim1 + y] = _make_cuComplex(F(0.0), F(0.0));
                }
            }
            if (y < testDim1 - 1)
            {
                if (y >= x - 1)
                {
                    heij[x * (testDim1 - 1) + y].x = GetRandomReal() * F(2.0) - F(1.0);
                    heij[x * (testDim1 - 1) + y].y = GetRandomReal() * F(2.0) - F(1.0);
                }
                else
                {
                    heij[x * (testDim1 - 1) + y] = _make_cuComplex(F(0.0), F(0.0));
                }
            }
        }

        //if (x < (testDim1 - 1))
        if (0 == x)
        {
            yij[x] = _make_cuComplex(
                GetRandomReal() * F(2.0) - F(1.0), 
                GetRandomReal() * F(2.0) - F(1.0));
        }
        else
        {
            yij[x] = _make_cuComplex(F(0.0), F(0.0));
        }
    }

    appGeneral(_T("\n(* ============= Copy the output below to Mathematica ============== *)\n"));

    appGeneral(_T("\nm=\n"));
    PrintMatrix(mij, testDim1, testDim2);
    appGeneral(_T("\nh=\n"));
    PrintMatrix(hij, testDim1, testDim1);
    appGeneral(_T("\nb=\n"));
    PrintMatrix(bij, testDim1, testDim1);
    appGeneral(_T("\na=\n"));
    PrintMatrix(aij, testDim1, testDim1);
    appGeneral(_T("\nii=\n"));
    PrintMatrix(iij, testDim1, testDim1);
    appGeneral(_T("\nt=\n"));
    PrintMatrix(tij, testDim1, testDim1);
    appGeneral(_T("\nhe=\n"));
    PrintMatrix(heij, testDim1, testDim1 - 1);
    appGeneral(_T("\nyy=\n"));
    PrintMatrix(yij, testDim1, 1);
    appGeneral(_T("\nhes=\n"));
    PrintMatrix(hesij, testDim1, testDim1);

    CLGComplex norm2 = pHelper->Normal2Host(yij, hij, testDim1);

    //transpose
    pHelper->TransposeHost(mij, testDim1, testDim2);
    appGeneral(_T("\ntm=\n"));
    PrintMatrix(mij, testDim2, testDim1);
    pHelper->TransposeHost(mij, testDim2, testDim1);

    //matrix multiply
    pHelper->SmallMatrixMultHost(res1, mij, mij, testDim2, testDim1, testDim2, TRUE, FALSE);
    appGeneral(_T("\nmm=\n"));
    PrintMatrix(res1, testDim2, testDim2);

    //block multiply
    pHelper->BlockMatrixMultHost(res1, hij, iij, testDim1, 0, testDim2, FALSE, TRUE, FALSE);
    appGeneral(_T("\nhdaggeri=\n"));
    PrintMatrix(res1, testDim1, testDim1);

    //norm2
    appGeneral(_T("\nnorm2res=%f %s %f I;\n"), 
        norm2.x, norm2.y < F(0.0) ? _T("") : _T("+"), norm2.y);

    //Solve Y
    pHelper->ThinQRFactorizationHost(res1, res2, mij, testDim1, testDim2);
    pHelper->SolveYHost(mij, tij, testDim2, testDim1);
    appGeneral(_T("\ninversetm=\n"));
    PrintMatrix(mij, testDim1, testDim2);

    //QR factor m is changed in Solve Y, so put this in fromt of Solve Y
    //pHelper->ThinQRFactorizationHost(res1, res2, mij, testDim1, testDim2);
    appGeneral(_T("\nq=\n"));
    PrintMatrix(res1, testDim1, testDim2);
    appGeneral(_T("\nr=\n"));
    PrintMatrix(res2, testDim2, testDim2);

    //Henssenberg
    pHelper->RotateHenssenbergHost(heij, yij, testDim1 - 1);
    appGeneral(_T("\nrhe=\n"));
    PrintMatrix(heij, testDim1 - 1, testDim1 - 1);
    appGeneral(_T("\nryy=\n"));
    PrintMatrix(yij, testDim1 - 1, 1);
    appGeneral(_T("\nrheresidue=%f %s %f I;\n"), 
        yij[testDim1 - 1].x,
        yij[testDim1 - 1].y < 0 ? _T("") : _T("+"),
        yij[testDim1 - 1].y);

    //EV
    pHelper->EigenValueProblemHost(hij, res1, res2, testDim1, testDim2);
    appGeneral(_T("\nev1=\n"));
    PrintMatrix(res1, 1, testDim2);
    appGeneral(_T("\nvv1=\n"));
    PrintMatrix(res2, testDim2, testDim1);

    //GEV
    pHelper->GeneralizedEigenValueProblemHost(aij, bij, res1, res2, testDim1, testDim2);
    appGeneral(_T("\nev2=\n"));
    PrintMatrix(res1, 1, testDim2);
    appGeneral(_T("\nvv2=\n"));
    PrintMatrix(res2, testDim2, testDim1);

    //TEV
    pHelper->UpperTriangularEigenVectorsHost(tij, res1, res2, testDim1, testDim2);
    appGeneral(_T("\nev3=\n"));
    PrintMatrix(res1, 1, testDim2);
    appGeneral(_T("\nvv3=\n"));
    PrintMatrix(res2, testDim2, testDim1);

    //HEV
    pHelper->EigenValueProblemHessenbergHost(hesij, res1, res2, testDim1, testDim2);
    appGeneral(_T("\nev4=\n"));
    PrintMatrix(res1, 1, testDim2);
    appGeneral(_T("\nvv4=\n"));
    PrintMatrix(res2, testDim2, testDim1);

    appGeneral(_T("\nPrint[\"Test matrix transpose\"]\n"));
    appGeneral(_T("\nMax[Transpose[m]-tm //Flatten//Abs]\n"));
    appGeneral(_T("\nPrint[\"Test matrix multiply\"]\n"));
    appGeneral(_T("\nMax[Conjugate[Transpose[m]].m-mm //Flatten//Abs]\n"));
    appGeneral(_T("\nMax[Conjugate[Transpose[h]].ii-hdaggeri //Flatten//Abs]\n"));
    appGeneral(_T("\nPrint[\"Test matrix normal2\"]\n"));
    appGeneral(_T("\nAbs[Conjugate[Transpose[yy]].h.yy-norm2res]\n"));
    appGeneral(_T("\nPrint[\"Test backward substitution\"]\n"));
    appGeneral(_T("\nMax[Inverse[t].m-inversetm //Flatten//Abs]/Max[Inverse[t].m // Flatten // Abs]\n"));

    appGeneral(_T("\nPrint[\"Test QR factorization\"]\n"));
    appGeneral(_T("\nMax[q.r-m //Flatten//Abs]\n"));
    appGeneral(_T("Tr[Conjugate[Transpose[q]].q] / %d - 1\n"), testDim2);

    appGeneral(_T("\nPrint[\"Test Henssenberg rotation\"]\n"));
    appGeneral(_T("\nNorm[he.(Inverse[rhe].ryy) - yy] - Abs[rheresidue]\n"));

    appGeneral(_T("\nPrint[\"Test Eigen Value\"]\n"));

    appGeneral(_T("Eigenvalues[h][[Table[%d - i, {i, 1, %d}]]] - ev1[[1]] // Abs // Max\n"), testDim1 + 1, testDim2);
    for (INT i = 0; i < testDim2; ++i)
    {
        appGeneral(_T("Max[Eigensystem[h][[2]][[%d]]/Eigensystem[h][[2]][[%d]][[%d]] - vv1[[%d]] / vv1[[%d]][[%d]] //Abs]\n"),
            testDim1 - i, testDim1 - i, testDim1, i + 1, i + 1, testDim1);
    }

    appGeneral(_T("\nPrint[\"Test Generalized Eigen Value\"]\n"));

    appGeneral(_T("Eigenvalues[Inverse[b].a][[Table[%d - i, {i, 1, %d}]]] - ev2[[1]] // Abs // Max\n"), testDim1 + 1, testDim2);
    for (INT i = 0; i < testDim2; ++i)
    {
        appGeneral(_T("Max[Eigensystem[Inverse[b].a][[2]][[%d]]/Eigensystem[Inverse[b].a][[2]][[%d]][[%d]] - vv2[[%d]] / vv2[[%d]][[%d]] //Abs]\n"),
            testDim1 - i, testDim1 - i, testDim1, i + 1, i + 1, testDim1);
    }

    appGeneral(_T("\nPrint[\"Test Upper triangular Eigen Value\"]\n"));

    appGeneral(_T("Eigenvalues[t][[Table[%d - i, {i, 1, %d}]]] - ev3[[1]] // Abs // Max\n"), testDim1 + 1, testDim2);
    for (INT i = 0; i < testDim2; ++i)
    {
        appGeneral(_T("Abs[Conjugate[vv3[[%d]]].Normalize[Eigensystem[t][[2]][[%d]]]] - 1\n"),
            i + 1, testDim1 - i);
    }

    appGeneral(_T("\nPrint[\"Test Hessenberg Eigen Value\"]\n"));

    appGeneral(_T("Eigenvalues[hes][[Table[%d - i, {i, 1, %d}]]] - ev4[[1]] // Abs // Max\n"), testDim1 + 1, testDim2);
    for (INT i = 0; i < testDim2; ++i)
    {
        appGeneral(_T("Abs[Conjugate[vv4[[%d]]].Normalize[Eigensystem[hes][[2]][[%d]]]] - 1\n"),
            i + 1, testDim1 - i);
    }

    appGeneral(_T("\n\n(*============= Please copy those results to Mathematica to check, all should be nearly zero ============ *)\n"));
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
__global__ void 
_kernelSmallMatrixMult_NN(CLGComplex* res,
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
__global__ void 
_kernelSmallMatrixMult_ND(CLGComplex* res,
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
__global__ void 
_kernelSmallMatrixMult_DN(CLGComplex* res,
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

__global__ void _CLG_LAUNCH_BOUND
_kernelNormal2(
    const CLGComplex* __restrict__ v,
    const CLGComplex* __restrict__ m,
    CLGComplex* res,
    UINT dm)
{
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    //assume res is zeroed
    CLGComplex toAdd = _cuCmulf(_cuConjf(v[x]), _cuCmulf(m[x * dm + y], v[y]));
    atomicAdd(&(res->x), toAdd.x);
    atomicAdd(&(res->y), toAdd.y);
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

void CLinearAlgebraHelper::Normal2(const CLGComplex* v, const CLGComplex* matrix, CLGComplex* devicees, UINT dm)
{
    dim3 block(1, 1, 1);
    dim3 thread(dm, dm, 1);
    _kernelNormal2 << <block, thread >> > (v, matrix, devicees, dm);
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
__global__ void
_kernelMatrixBlockMult_LNN(CLGComplex* res, const CLGComplex* __restrict__ left, const CLGComplex* __restrict__ right, UINT iStart, UINT iEnd, UINT dm)
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

__global__ void
_kernelMatrixBlockMult_LDN(CLGComplex* res, const CLGComplex* __restrict__ left, const CLGComplex* __restrict__ right, UINT iStart, UINT iEnd, UINT dm)
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

__global__ void
_kernelMatrixBlockMult_LND(CLGComplex* res, const CLGComplex* __restrict__ left, const CLGComplex* __restrict__ right, UINT iStart, UINT iEnd, UINT dm)
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
__global__ void
_kernelMatrixBlockMult_RNN(CLGComplex* res, const CLGComplex* __restrict__ left, const CLGComplex* __restrict__ right, UINT iStart, UINT iEnd, UINT dm)
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

__global__ void
_kernelMatrixBlockMult_RDN(CLGComplex* res, const CLGComplex* __restrict__ left, const CLGComplex* __restrict__ right, UINT iStart, UINT iEnd, UINT dm)
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

__global__ void
_kernelMatrixBlockMult_RND(CLGComplex* res, const CLGComplex* __restrict__ left, const CLGComplex* __restrict__ right, UINT iStart, UINT iEnd, UINT dm)
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

__global__ void _CLG_LAUNCH_BOUND
_kernelMatrixTranspose(const CLGComplex* __restrict__ m, CLGComplex* tmpM, UINT dx, UINT dy)
{
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    tmpM[y * dx + x] = m[x * dy + y];
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMatrixDagger(const CLGComplex* __restrict__ m, CLGComplex* tmpM, UINT dx, UINT dy)
{
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    tmpM[y * dx + x] = _cuConjf(m[x * dy + y]);
}

void CLinearAlgebraHelper::Transpose(CLGComplex* deviceMatrix, UINT dx, UINT dy)
{
    if (dx > m_uiDim || dy > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpMRes = GetTmpMatrix();
    CLGComplex* tmpMRes = sTmpMRes.m_pMatrix;

    dim3 block(1, 1, 1);
    dim3 thread(dx, dy, 1);
    _kernelMatrixTranspose << <block, thread >> > (deviceMatrix, tmpMRes, dx, dy);
    checkCudaErrors(cudaMemcpy(deviceMatrix, tmpMRes, sizeof(CLGComplex) * dx * dy, cudaMemcpyDeviceToDevice));

    sTmpMRes.Free();
}

void CLinearAlgebraHelper::Dagger(CLGComplex* deviceMatrix, UINT dx, UINT dy)
{
    if (dx > m_uiDim || dy > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpMRes = GetTmpMatrix();
    CLGComplex* tmpMRes = sTmpMRes.m_pMatrix;

    dim3 block(1, 1, 1);
    dim3 thread(dx, dy, 1);
    _kernelMatrixDagger << <block, thread >> > (deviceMatrix, tmpMRes, dx, dy);
    checkCudaErrors(cudaMemcpy(deviceMatrix, tmpMRes, sizeof(CLGComplex) * dx * dy, cudaMemcpyDeviceToDevice));

    sTmpMRes.Free();
}

void CLinearAlgebraHelper::TransposeHost(CLGComplex* hostMatrix, UINT dx, UINT dy)
{
    if (dx > m_uiDim || dy > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpMRes = GetTmpMatrix();
    CLGComplex* tmpMRes = sTmpMRes.m_pMatrix;
    checkCudaErrors(cudaMemcpy(tmpMRes, hostMatrix, sizeof(CLGComplex) * dx * dy, cudaMemcpyHostToDevice));
    Transpose(tmpMRes, dx, dy);
    checkCudaErrors(cudaMemcpy(hostMatrix, tmpMRes, sizeof(CLGComplex) * dx * dy, cudaMemcpyDeviceToHost));

    sTmpMRes.Free();
}

void CLinearAlgebraHelper::DaggerHost(CLGComplex* hostMatrix, UINT dx, UINT dy)
{
    if (dx > m_uiDim || dy > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpMRes = GetTmpMatrix();
    CLGComplex* tmpMRes = sTmpMRes.m_pMatrix;
    checkCudaErrors(cudaMemcpy(tmpMRes, hostMatrix, sizeof(CLGComplex) * dx * dy, cudaMemcpyHostToDevice));
    Dagger(tmpMRes, dx, dy);
    checkCudaErrors(cudaMemcpy(hostMatrix, tmpMRes, sizeof(CLGComplex) * dx * dy, cudaMemcpyDeviceToHost));

    sTmpMRes.Free();
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

__global__ void _CLG_LAUNCH_BOUND
_kernelMatrixBlockCopy(
    CLGComplex* dest, 
    const CLGComplex* __restrict__ src, 
    UINT srcX, UINT srcY, UINT destX, UINT destY, UINT srcDim, UINT destDim)
{
    UINT x = threadIdx.x;
    UINT y = threadIdx.y;

    dest[(x + destX) * destDim + (y + destY)] = src[(x + srcX) * srcDim + (y + srcY)];
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
            appGeneral("%s%1.8f %s %1.8f I%s ",
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
    __shared__ UBOOL bValidLength;

    if (0 == x && 0 == y)
    {
        length = F(0.0);
        lengthu = F(0.0);
        bValidLength = TRUE;
    }

    __syncthreads();

    if (0 == y && x >= i)
    {
        atomicAdd(&length, R[x * dy + i].x * R[x * dy + i].x + R[x * dy + i].y * R[x * dy + i].y);
    }

    __syncthreads();

    if (0 == x && 0 == y)
    {
        if (length < _CLG_FLT_MIN_)
        {
            bValidLength = FALSE;
        }
        else
        {
            length = _sqrt(length);
        }
    }

    __syncthreads();

    if (0 == y && x >= i && bValidLength)
    {
        u[x] = R[x * dy + i];

        if (x == i)
        {
            Real fLengthuxL = __cuCabsSqf(u[x]);
            Real fCos = F(0.0);
            Real fSin = F(0.0);
            if (fLengthuxL > _CLG_FLT_MIN_)
            { 
                Real fuxL = __div(F(1.0), _sqrt(fLengthuxL));
                fCos = u[x].x * fuxL;
                fSin = u[x].y * fuxL;
            }
            else
            {
                Real fArg = _atan2(u[x].y, u[x].x);
                fCos = _cos(fArg);
                fSin = _sin(fArg);
            }
            u[x] = _cuCaddf(u[x], _make_cuComplex(length * fCos, length * fSin));
        }
        atomicAdd(&lengthu, u[x].x * u[x].x + u[x].y * u[x].y);
    }

    __syncthreads();

    if (0 == x && 0 == y && bValidLength)
    {
        lengthu = lengthu * F(0.5);
        if (lengthu < _CLG_FLT_MIN_)
        {
            bValidLength = FALSE;
        }
        else
        {
            lengthu = __div(F(1.0), _sqrt(lengthu));
        }
    }

    __syncthreads();

    if (0 == y && x >= i)
    {
        if (bValidLength)
        {
            u[x].x = u[x].x * lengthu;
            u[x].y = u[x].y * lengthu;
        }
        else
        {
            u[x] = _make_cuComplex(F(0.0), F(0.0));
        }
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
    __shared__ UBOOL bValidLength;
    __shared__ CLGComplex u[CLinearAlgebraHelper::_kMaxSmallDim];
    if (0 == x && 0 == y)
    {
        length = F(0.0);
        lengthu = F(0.0);
        bValidLength = TRUE;
    }

    __syncthreads();

    if (0 == y && x > i)
    {
        atomicAdd(&length, A[x * dx + i].x * A[x * dx + i].x + A[x * dx + i].y * A[x * dx + i].y);
    }

    __syncthreads();

    if (0 == x && 0 == y)
    {
        if (length < _CLG_FLT_MIN_)
        {
            bValidLength = FALSE;
        }
        else
        {
            length = _sqrt(length);
        }
    }

    __syncthreads();

    if (0 == y && x > i && bValidLength)
    {
        u[x] = A[x * dx + i];

        if (x == i + 1)
        {
            Real fuxL = __cuCabsSqf(u[x]);
            Real fCos = F(0.0);
            Real fSin = F(0.0);
            if (fuxL > _CLG_FLT_MIN_)
            {
                fuxL = __div(F(1.0), _sqrt(fuxL));
                fCos = u[x].x * fuxL;
                fSin = u[x].y * fuxL;
            }
            else
            {
                Real fArg = _atan2(u[x].y, u[x].x);
                fCos = _cos(fArg);
                fSin = _sin(fArg);
            }
            u[x] = _cuCaddf(u[x], _make_cuComplex(length * fCos, length * fSin));
        }
        atomicAdd(&lengthu, u[x].x * u[x].x + u[x].y * u[x].y);
    }

    __syncthreads();

    if (0 == x && 0 == y && bValidLength)
    {
        lengthu = lengthu * F(0.5);
        if (lengthu < _CLG_FLT_MIN_)
        {
            bValidLength = FALSE;
        }
        else
        {
            lengthu = __div(F(1.0), _sqrt(lengthu));
        }
    }

    __syncthreads();

    if (0 == y && x > i)
    {
        if (bValidLength)
        {
            u[x].x = u[x].x * lengthu;
            u[x].y = u[x].y * lengthu;
        }
        else
        {
            u[x] = _make_cuComplex(F(0.0), F(0.0));
        }
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

#pragma region Backward Substitude

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
_kernelWilkinsonShift(CLGComplex* m, CLGComplex* c, Real fCrit, UINT dim)
{
    UINT i = threadIdx.x;
    if (0 == i)
    {
        //d
        c[0] = m[dim * dim - 1];

        //bc
        CLGComplex omega = _cuCmulf(m[dim * dim - dim - 1], m[dim * dim - 2]);

        Real fOmegaSq = __cuCabsSqf(omega);
        if (fOmegaSq > fCrit)
        {
            //(d-a)/2
            CLGComplex xi = _make_cuComplex(
                F(0.5) * (c[0].x - m[dim * dim - dim - 2].x),
                F(0.5) * (c[0].y - m[dim * dim - dim - 2].y));

            //sqrt(((d-a)/2)^2 + bc)
            CLGComplex eta = _cuCaddf(_cuCmulf(xi, xi), omega);
            if (__cuCabsSqf(eta) > _CLG_FLT_MIN_)
            {
                eta = __cuCsqrtf(eta);
            }
            else
            {
                eta = _make_cuComplex(F(0.0), F(0.0));
            }
            if (xi.x * eta.x + xi.y * eta.y < F(0.0))
            {
                CLGComplex divider = _cuCsubf(eta, xi);
                if (__cuCabsSqf(divider) > _CLG_FLT_MIN_)
                {
                    c[0] = _cuCsubf(c[0], _cuCdivf(omega, divider));
                }
            }
            else
            {
                CLGComplex divider = _cuCaddf(eta, xi);
                if (__cuCabsSqf(divider) > _CLG_FLT_MIN_)
                {
                    c[0] = _cuCsubf(c[0], _cuCdivf(omega, divider));
                }
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



__device__ __inline__ static void _deviceCalculateEigenValueTwo(
    CLGComplex& h00,
    CLGComplex& h01,
    CLGComplex& h10,
    CLGComplex& h11,
    Real fCrit)
{
    //bc
    CLGComplex omega = _cuCmulf(h10, h01);

    Real fOmegaSq = omega.x * omega.x + omega.y * omega.y;
    if (fOmegaSq > fCrit)
    {
        //(d-a)/2
        CLGComplex xi = _make_cuComplex(
            F(0.5) * (h11.x - h00.x),
            F(0.5) * (h11.y - h00.y));
        //sqrt(((d-a)/2)^2 + bc)
        //when bc > fCrit, eta=sqrt(xi^2+bc) cannot be close to xi
        CLGComplex eta = _cuCaddf(_cuCmulf(xi, xi), omega);
        if (__cuCabsSqf(eta) > _CLG_FLT_MIN_)
        {
            eta = __cuCsqrtf(eta);
        }
        else
        {
            eta = _make_cuComplex(F(0.0), F(0.0));
        }
        CLGComplex divider1 = _cuCaddf(eta, xi);
        CLGComplex divider2 = _cuCsubf(eta, xi);
        if (__cuCabsSqf(divider1) > _CLG_FLT_MIN_ && __cuCabsSqf(divider2) > _CLG_FLT_MIN_)
        {
            if (xi.x * eta.x + xi.y * eta.y < F(0.0))
            {
                h00 = _cuCaddf(h11, _cuCdivf(omega, _cuCaddf(eta, xi)));
                h11 = _cuCsubf(h11, _cuCdivf(omega, _cuCsubf(eta, xi)));
            }
            else
            {
                h00 = _cuCsubf(h11, _cuCdivf(omega, _cuCsubf(eta, xi)));
                h11 = _cuCaddf(h11, _cuCdivf(omega, _cuCaddf(eta, xi)));
            }
        }
    }
    h10 = _make_cuComplex(F(0.0), F(0.0));
}

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernel2By2Eigen(CLGComplex* matrix, INT* decomp, UINT dx, Real fCrit)
{
    _deviceCalculateEigenValueTwo(
        matrix[decomp[0] * dx + decomp[0]],
        matrix[decomp[0] * dx + decomp[0] + 1],
        matrix[(decomp[0] + 1) * dx + decomp[0]],
        matrix[(decomp[0] + 1) * dx + decomp[0] + 1],
        fCrit
    );
}

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelCheckMatrixDoubleShift(CLGComplex* mtr, INT* decomp, INT dx, Real fCrit)
{
    decomp[0] = 0;
    decomp[1] = dx;

    for (INT i = dx - 2; i >= 0; --i)
    {
        if (__cuCabsSqf(mtr[(i + 1) * dx + i]) <
            fCrit * (__cuCabsSqf(mtr[(i + 1) * dx + i + 1]) + __cuCabsSqf(mtr[i * dx + i])))
        {
            mtr[(i + 1) * dx + i] = _make_cuComplex(F(0.0), F(0.0));

            if (decomp[1] == i + 2)
            {
                decomp[1] = i + 1;
            }

            if (i + 1 > decomp[0] && i + 1 < decomp[1])
            {
                decomp[0] = i + 1;
            }
        }
    }

    //printf("decomp %d to %d \n", decomp[0], decomp[1]);
}

__device__ __inline__ static void _deviceTwoHouseHolder(CLGComplex& a, CLGComplex& b)
{
    Real len = __cuCabsSqf(a) + __cuCabsSqf(b);
    if (len < _CLG_FLT_MIN_)
    {
        a = _make_cuComplex(F(0.0), F(0.0));
        b = _make_cuComplex(F(0.0), F(0.0));
        return;
    }
    len = _sqrt(len);
    Real lena = a.x * a.x + a.y * a.y;
    Real fCos = F(0.0);
    Real fSin = F(0.0);
    if (lena < _CLG_FLT_MIN_)
    {
        Real fArg = _atan2(a.y, a.x);
        fCos = _cos(fArg);
        fSin = _sin(fArg);
    }
    else
    {
        lena = __div(F(1.0), _sqrt(lena));
        fCos = a.x * lena;
        fSin = a.y * lena;
    }

    a = _cuCaddf(a, _make_cuComplex(len * fCos, len * fSin));

    Real len2 = F(0.5) * (a.x * a.x + a.y * a.y + b.x * b.x + b.y * b.y);
    if (len2 < _CLG_FLT_MIN_)
    {
        a = _make_cuComplex(F(0.0), F(0.0));
        b = _make_cuComplex(F(0.0), F(0.0));
        return;
    }
    len2 = __div(F(1.0), _sqrt(len2));
    a.x = a.x * len2;
    a.y = a.y * len2;
    b.x = b.x * len2;
    b.y = b.y * len2;
}

__device__ __inline__ static void _deviceThreeHouseHolder(CLGComplex& a, CLGComplex& b, CLGComplex& c)
{
    Real len = __cuCabsSqf(a) + __cuCabsSqf(b) + __cuCabsSqf(c);
    if (len < _CLG_FLT_MIN_)
    {
        a = _make_cuComplex(F(0.0), F(0.0));
        b = _make_cuComplex(F(0.0), F(0.0));
        c = _make_cuComplex(F(0.0), F(0.0));
        return;
    }
    len = _sqrt(len);
    Real lena = a.x * a.x + a.y * a.y;
    Real fCos = F(0.0);
    Real fSin = F(0.0);
    if (lena < _CLG_FLT_MIN_)
    {
        Real fArg = _atan2(a.y, a.x);
        fCos = _cos(fArg);
        fSin = _sin(fArg);
    }
    else
    {
        lena = __div(F(1.0), _sqrt(lena));
        fCos = a.x * lena;
        fSin = a.y * lena;
    }
    a = _cuCaddf(a, _make_cuComplex(len * fCos, len * fSin));

    Real len2 = F(0.5) * (a.x * a.x + a.y * a.y + b.x * b.x + b.y * b.y + c.x * c.x + c.y * c.y);
    if (len2 < _CLG_FLT_MIN_)
    {
        a = _make_cuComplex(F(0.0), F(0.0));
        b = _make_cuComplex(F(0.0), F(0.0));
        c = _make_cuComplex(F(0.0), F(0.0));
        return;
    }
    len2 = __div(F(1.0), _sqrt(len2));
    a.x = a.x * len2;
    a.y = a.y * len2;
    b.x = b.x * len2;
    b.y = b.y * len2;
    c.x = c.x * len2;
    c.y = c.y * len2;
}

/**
* H = h00 h01
*     h10 h11
* a1, a2 be the eigenvalues
* s = a1 + a2
* t = a1.a2
*/
__device__ __inline__ static void _deviceDoubleEigen(CLGComplex& s, CLGComplex& t,
    const CLGComplex& h00, const CLGComplex& h01,
    const CLGComplex& h10, const CLGComplex& h11)
{
    s = _cuCaddf(h11, h00);
    t = _cuCsubf(_cuCmulf(h00, h11), _cuCmulf(h10, h01));
}

/**
* H = h00 h01
*     h10 h11
* a1, a2 be the eigenvalues
* s = a1 + a2
* t = a1.a2
*/
__device__ __inline__ static void _deviceDoubleShift(
    CLGComplex& a, CLGComplex& b, CLGComplex& c,
    const CLGComplex& s, const CLGComplex& t,
    const CLGComplex& h00, const CLGComplex& h01,
    const CLGComplex& h10, const CLGComplex& h11,
    const CLGComplex& h21)
{
    a = _cuCaddf(
        _cuCmulf(h00, _cuCsubf(h00, s)),
        _cuCaddf(_cuCmulf(h10, h01), t));
    b = _cuCmulf(h10,
        _cuCsubf(_cuCaddf(h00, h11), s)
    );
    c = _cuCmulf(h10, h21);
}

/**
* thread = 1,1,1
*/
__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelStartStep(const CLGComplex* __restrict__ H, CLGComplex* xyz, UINT dm)
{
    CLGComplex s, t;
    _deviceDoubleEigen(s, t,
        H[dm * dm - dm - 2/*(dm - 2) * dm + dm - 2*/],
        H[dm * dm - dm - 1/*(dm - 2) * dm + dm - 1*/],
        H[dm * dm - 2/*(dm - 1) * dm + dm - 2*/],
        H[dm * dm - 1/*(dm - 1) * dm + dm - 1*/]);

    _deviceDoubleShift(xyz[0], xyz[1], xyz[2], s, t,
        H[0],
        H[1],
        H[dm],
        H[dm + 1],
        H[2 * dm + 1]);

    //printf("s=%f %f, t=%f %f, x=%f %f, y=%f %f, z=%f %f\n",
    //    s.x, s.y,
    //    t.x, t.y,
    //    xyz[0].x, xyz[0].y,
    //    xyz[1].x, xyz[1].y,
    //    xyz[2].x, xyz[2].y);
}

/**
* thread(tx, ty, 1)
* k = 0, 1, tx = dm
* k = 2,... tx = dm - k + 1
* ty = 9
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelStepK_1(
    CLGComplex* H,
    CLGComplex* um,
    const CLGComplex* __restrict__ xyz,
    UINT k, UINT dm)
{
    //k=0, 1,  q = 0 to dm - 1
    //k=2,..., q = k-1 to dm - 1
    __shared__ CLGComplex lines[3][CLinearAlgebraHelper::_kMaxSmallDim];
    __shared__ CLGComplex house[3];

    UINT y = threadIdx.x;
    if (k >= 2)
    {
        y = y + k - 1;
    }
    UINT ux = threadIdx.y / 3;
    UINT x = ux + k;
    UINT uy = threadIdx.y % 3;

    //we are going to change all H(x, y)
    if (0 == threadIdx.x && 0 == ux && 1 == uy)
    {
        house[0] = xyz[0];
        house[1] = xyz[1];
        house[2] = xyz[2];
        _deviceThreeHouseHolder(house[0], house[1], house[2]);
        //u [i * 3 + j] = h[j]* h[i]
        um[0] = _cuCmulf(_cuConjf(house[0]), house[0]);
        um[1] = _cuCmulf(_cuConjf(house[1]), house[0]);
        um[2] = _cuCmulf(_cuConjf(house[2]), house[0]);
        um[3] = _cuCmulf(_cuConjf(house[0]), house[1]);
        um[4] = _cuCmulf(_cuConjf(house[1]), house[1]);
        um[5] = _cuCmulf(_cuConjf(house[2]), house[1]);
        um[6] = _cuCmulf(_cuConjf(house[0]), house[2]);
        um[7] = _cuCmulf(_cuConjf(house[1]), house[2]);
        um[8] = _cuCmulf(_cuConjf(house[2]), house[2]);

        //u[i * 3 + i] -= 1
        um[0].x = um[0].x - F(1.0);
        um[4].x = um[4].x - F(1.0);
        um[8].x = um[8].x - F(1.0);
    }

    if (0 == uy)
    {
        lines[ux][y] = H[x * dm + y];
        H[x * dm + y] = _make_cuComplex(F(0.0), F(0.0));
    }

    __syncthreads();

    //U = 1 - vvT
    CLGComplex u = _cuCmulf(um[ux * 3 + uy], lines[uy][y]);

    atomicAdd(&H[x * dm + y].x, -u.x);
    atomicAdd(&H[x * dm + y].y, -u.y);

}

/**
* thread(tx, ty, 1)
* k = 0,..,n-4 tx = k + 4
* k = n-3      tx = dm
* ty = 9
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelStepK_2(
    CLGComplex* H,
    const CLGComplex* __restrict__ um,
    CLGComplex* xyz,
    UINT k, UINT dm)
{
    UINT x = threadIdx.x;
    UINT ux = threadIdx.y / 3;
    UINT uy = threadIdx.y % 3;

    UINT y = uy + k;

    __shared__ CLGComplex lines[CLinearAlgebraHelper::_kMaxSmallDim][3];
    if (0 == ux)
    {
        lines[x][uy] = H[x * dm + y];
        H[x * dm + y] = _make_cuComplex(F(0.0), F(0.0));
    }

    __syncthreads();

    CLGComplex u = _cuCmulf(lines[x][ux], um[ux * 3 + uy]);

    atomicAdd(&H[x * dm + y].x, -u.x);
    atomicAdd(&H[x * dm + y].y, -u.y);

    __syncthreads();

    if (0 == x && 0 == ux)
    {
        int nextrow = k + uy + 1;
        if (nextrow < dm)
        {
            xyz[uy] = H[nextrow * dm + k];
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDoubleShiftFinal(CLGComplex* H, CLGComplex* xyz, UINT dm)
{
    UINT x = threadIdx.x;
    __shared__ CLGComplex house[2];
    __shared__ CLGComplex um[4];
    __shared__ CLGComplex lines[2][CLinearAlgebraHelper::_kMaxSmallDim];

    if (0 == x)
    {
        house[0] = xyz[0];
        house[1] = xyz[1];
        _deviceTwoHouseHolder(house[0], house[1]);
        //u [i * 2 + j] = h[j]* h[i]
        um[0] = _cuCmulf(_cuConjf(house[0]), house[0]);
        um[1] = _cuCmulf(_cuConjf(house[1]), house[0]);
        um[2] = _cuCmulf(_cuConjf(house[0]), house[1]);
        um[3] = _cuCmulf(_cuConjf(house[1]), house[1]);

        //u[i * 2 + i] -= 1
        um[0].x = um[0].x - F(1.0);
        um[3].x = um[3].x - F(1.0);
    }

    if (x >= dm - 3)
    {
        lines[0][x] = H[(dm - 2) * dm + x];
        lines[1][x] = H[(dm - 1) * dm + x];
    }

    __syncthreads();

    if (x >= dm - 3)
    {
        CLGComplex res1 = _cuCmulf(um[0 * 2 + 0], lines[0][x]);
        CLGComplex res2 = _cuCmulf(um[0 * 2 + 1], lines[1][x]);
        H[(dm - 2) * dm + x] = _make_cuComplex(-res1.x - res2.x, -res1.y - res2.y);
        res1 = _cuCmulf(um[1 * 2 + 0], lines[0][x]);
        res2 = _cuCmulf(um[1 * 2 + 1], lines[1][x]);
        H[(dm - 1) * dm + x] = _make_cuComplex(-res1.x - res2.x, -res1.y - res2.y);
    }

    __syncthreads();

    lines[0][x] = H[x * dm + (dm - 2)];
    lines[1][x] = H[x * dm + (dm - 1)];

    __syncthreads();

    CLGComplex res1 = _cuCmulf(um[0 * 2 + 0], lines[0][x]);
    CLGComplex res2 = _cuCmulf(um[1 * 2 + 0], lines[1][x]);
    H[x * dm + (dm - 2)] = _make_cuComplex(-res1.x - res2.x, -res1.y - res2.y);
    res1 = _cuCmulf(um[0 * 2 + 1], lines[0][x]);
    res2 = _cuCmulf(um[1 * 2 + 1], lines[1][x]);
    H[x * dm + (dm - 1)] = _make_cuComplex(-res1.x - res2.x, -res1.y - res2.y);
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

#if _CLG_QR_Left_Top_Deflation

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
                //appParanoiac(_T("(* total iteration = %d *)\n"), i + 1);
                return;
            }

            iLastDim = endindex[0];
            //must do the copy, because [x * dim + y] the "dim" is changed.
            BlockCopy(tmpT, T, iLastDim, iLastDim, iLastDim, dx);
        }

        //shift
        //T = T - sigma I, tmpDeviceFloat[0] = sigma
        dim3 thread(iLastDim, 1, 1);
        _kernelWilkinsonShift << <block, thread >> > (tmpT, tmpShift, fCrit, iLastDim);

        //QR decompose
        QRFactorization(tmpQ, tmpR, tmpT, iLastDim);

        //Update H
        //T = R Q + sigma I
        SmallMatrixMult(tmpT, tmpR, tmpQ, iLastDim, iLastDim, iLastDim, FALSE, FALSE);
        _kernelMatrixAddConstant << <block, thread >> > (tmpT, tmpShift, iLastDim);
    }
    BlockCopy(tmpT, T, iLastDim, iLastDim, iLastDim, dx);
#else
    INT endindex[2];
    for (UINT i = 0; i < iCrit; ++i)
    {
        //find decomp
        _kernelCheckMatrixDoubleShift << <1, 1 >> > (T, tmpDecomp, dx, fCrit);

        checkCudaErrors(cudaMemcpy(endindex, tmpDecomp, sizeof(INT) * 2, cudaMemcpyDeviceToHost));
        INT iLength = endindex[1] - endindex[0];
        if (iLength < 2)
        {
            sTmpT.Free();
            //appParanoiac(_T("(* total iteration = %d *)\n"), i + 1);
            //finished
            return;
        }
        else if (2 == iLength)
        {
            _kernel2By2Eigen << <1, 1 >> > (T, tmpDecomp, dx, fCrit);
        }
        else
        {
            dim3 threadCopy(iLength, iLength, 1);
            _kernelMatrixBlockCopy << <block, threadCopy >> > (tmpT, T,
                endindex[0], endindex[0], 0, 0, dx, iLength);

            //shift
            //T = T - sigma I, tmpDeviceFloat[0] = sigma
            dim3 thread(iLength, 1, 1);
            _kernelWilkinsonShift << <block, thread >> > (tmpT, tmpShift, fCrit, iLength);

            //QR decompose
            QRFactorization(tmpQ, tmpR, tmpT, iLength);

            //Update H
            //T = R Q + sigma I
            SmallMatrixMult(tmpT, tmpR, tmpQ, iLength, iLength, iLength, FALSE, FALSE);
            _kernelMatrixAddConstant << <block, thread >> > (tmpT, tmpShift, iLength);

            _kernelMatrixBlockCopy << <block, threadCopy >> > (T, tmpT,
                0, 0, endindex[0], endindex[0], iLength, dx);
        }
    }
#endif

    appCrucial(_T("QRITeration not converge!!\n"));
    sTmpQ.Free();
    sTmpR.Free();
    sTmpT.Free();
}

void CLinearAlgebraHelper::FrancisQRIterateBlock(CLGComplex* T, CLGComplex * tmpXYZ, UINT uiBlockDim)
{
    dim3 block(1, 1, 1);
    _kernelStartStep << <1, 1 >> > (T, tmpXYZ, uiBlockDim);

    for (UINT k = 0; k <= uiBlockDim - 3; ++k)
    {
        dim3 thread1(k < 2 ? uiBlockDim : (uiBlockDim - k + 1), 9, 1);
        _kernelStepK_1 << <block, thread1 >> > (T, 
            m_pDeviceTmpHouseHolder, tmpXYZ, k, uiBlockDim);

        dim3 thread2((uiBlockDim - 3 == k) ? uiBlockDim : (k + 4), 9, 1);
        _kernelStepK_2 << <block, thread2 >> > (T, 
            m_pDeviceTmpHouseHolder, tmpXYZ, k, uiBlockDim);
    }

    dim3 thread3(uiBlockDim, 1, 1);
    _kernelDoubleShiftFinal << <block, thread3 >> > (T, tmpXYZ, uiBlockDim);
}

void CLinearAlgebraHelper::FrancisQRIterate(CLGComplex* T, UINT dx, Real fCrit, UINT iCrit)
{
    if (dx > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    dim3 block1(1, 1, 1);
    INT endindex[2];

    STmpMatrix sTmpT = GetTmpMatrix();
    CLGComplex* tmpT = sTmpT.m_pMatrix;
    INT* tmpDecomp = m_pDeviceIntBuffer;
    CLGComplex * tmpXYZ = m_pDeviceComplexBuffer2;

    for (UINT i = 0; i < iCrit; ++i)
    {
        //find decomp
        _kernelCheckMatrixDoubleShift << <1, 1 >> > (T, tmpDecomp, dx, fCrit);

        checkCudaErrors(cudaMemcpy(endindex, tmpDecomp, sizeof(INT) * 2, cudaMemcpyDeviceToHost));
        INT iLength = endindex[1] - endindex[0];
        if (iLength < 2)
        {
            sTmpT.Free();
            //appParanoiac(_T("(* total iteration = %d *)\n"), i + 1);
            //finished
            return;
        }
        else if (2 == iLength)
        {
            _kernel2By2Eigen << <1, 1 >> > (T, tmpDecomp, dx, fCrit);
        }
        else
        {
            dim3 threadCopy(iLength, iLength, 1);
            _kernelMatrixBlockCopy << <block1, threadCopy >> > (tmpT, T,
                endindex[0], endindex[0], 0, 0, dx, iLength);
            FrancisQRIterateBlock(tmpT, tmpXYZ, iLength);
            _kernelMatrixBlockCopy << <block1, threadCopy >> > (T, tmpT,
                0, 0, endindex[0], endindex[0], iLength, dx);
        }
    }

    appCrucial(_T("FrancisQRIterate not converge!!\n"));
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
_kernelSortEigenValuesBig(const CLGComplex* __restrict__ R,
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
        if (tmpF[x] > tmpF[y])
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
    y[j] = _cuConjf(Q[j * dx]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelInverseIterateShift(CLGComplex* A, const CLGComplex* __restrict__ outV, UINT k, UINT dx)
{
    UINT x = threadIdx.x;
    A[x * dx + x] = _cuCsubf(A[x * dx + x], outV[k]);
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

    if (fAmp[x] < _CLG_FLT_MIN)
    {
        return;
    }

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

    CLGComplex toAdd = _cuCmulf(A[x * dx + y], v[y]);
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
    UINT dm, UINT dk, UBOOL bSmall,
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
#if _CLG_IMPLICITE_QR_SHIFT
    FrancisQRIterate(tmpH, dm, fCrit, iMaxIterate);
#else
    QRIterate(tmpH, dm, fCrit, iMaxIterate);
#endif

    STmpMatrix sTmpQ = GetTmpMatrix();
    CLGComplex* tmpQ = sTmpQ.m_pMatrix;
    STmpMatrix sTmpR = GetTmpMatrix();
    CLGComplex* tmpR = sTmpR.m_pMatrix;

    dim3 block(1, 1, 1);
    dim3 thread1(dm, dm, 1);
    dim3 thread2(dm, 1, 1);

    if (bSmall)
    {
        _kernelSortEigenValues << <block, thread1 >> > (tmpH, outEigenValue, m_pDeviceFloatBuffer, m_pDeviceIntBuffer, dk, dm);
    }
    else
    {
        _kernelSortEigenValuesBig << <block, thread1 >> > (tmpH, outEigenValue, m_pDeviceFloatBuffer, m_pDeviceIntBuffer, dk, dm);
    }
#if _CLG_QRIterate_Update_EigenValue
    CLGComplex res[1];
    res[0] = _make_cuComplex(F(0.0), F(0.0));
#endif
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

            //It is normalized in _kernelErrorCheck
            _kernelErrorCheck << <block, thread1 >> > (tmpF, tmpVector, tmpH, dm);

            checkCudaErrors(cudaMemcpy(fErr, tmpF, sizeof(float), cudaMemcpyDeviceToHost));

            if (j == iMaxEigenIterate - 1)
            {
                appParanoiac(_T("(* Eigen vector(%d-iterate %d) error now = %1.12f *)\n"), i + 1, j, fErr[0]);
            }

            if (fErr[0] < fEigenCrit)
            {
                //appParanoiac(_T("(* ev converge after %d iteration *)\n"), j + 1);
                break;
            }

#if _CLG_QRIterate_Update_EigenValue
            checkCudaErrors(cudaMemcpy(m_pOneDeviceC, res, sizeof(CLGComplex), cudaMemcpyHostToDevice));
            Normal2(tmpVector, H, m_pOneDeviceC, dm);
            checkCudaErrors(cudaMemcpy(tmpH, H, sizeof(CLGComplex) * dm * dm, cudaMemcpyDeviceToDevice));
            _kernelInverseIterateShift << <block, thread2 >> > (tmpH, m_pOneDeviceC, 0, dm);
            QRFactorization(tmpQ, tmpR, tmpH, dm);
#endif
            //Do the inverse here
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

void CLinearAlgebraHelper::EigenValueProblemHessenberg(
    CLGComplex* H, CLGComplex* outEigenValue, CLGComplex* outEigenVector,
    UINT dm, UINT dk, UBOOL bSmall, Real fEigenCrit, UINT iMaxEigenIterate,
    Real fQRCrit, UINT iMaxIterate)
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

#if _CLG_IMPLICITE_QR_SHIFT
    FrancisQRIterate(tmpH, dm, fQRCrit, iMaxIterate);
#else
    QRIterate(tmpH, dm, fQRCrit, iMaxIterate);
#endif

    STmpMatrix sTmpQ = GetTmpMatrix();
    CLGComplex* tmpQ = sTmpQ.m_pMatrix;
    STmpMatrix sTmpR = GetTmpMatrix();
    CLGComplex* tmpR = sTmpR.m_pMatrix;

    dim3 block(1, 1, 1);
    dim3 thread1(dm, dm, 1);
    dim3 thread2(dm, 1, 1);

    if (bSmall)
    {
        _kernelSortEigenValues << <block, thread1 >> > (tmpH, outEigenValue, m_pDeviceFloatBuffer, m_pDeviceIntBuffer, dk, dm);
    }
    else
    {
        _kernelSortEigenValuesBig << <block, thread1 >> > (tmpH, outEigenValue, m_pDeviceFloatBuffer, m_pDeviceIntBuffer, dk, dm);
    }
#if _CLG_QRIterate_Update_EigenValue
    CLGComplex res[1];
    res[0] = _make_cuComplex(F(0.0), F(0.0));
#endif
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
                //appParanoiac(_T("(* ev converge after %d iteration *)\n"), j + 1);
                break;
            }

#if _CLG_QRIterate_Update_EigenValue
            checkCudaErrors(cudaMemcpy(m_pOneDeviceC, res, sizeof(CLGComplex), cudaMemcpyHostToDevice));
            Normal2(tmpVector, H, m_pOneDeviceC, dm);
            checkCudaErrors(cudaMemcpy(tmpH, H, sizeof(CLGComplex) * dm * dm, cudaMemcpyDeviceToDevice));
            _kernelInverseIterateShift << <block, thread2 >> > (tmpH, m_pOneDeviceC, 0, dm);
            QRFactorization(tmpQ, tmpR, tmpH, dm);
#endif

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

__global__ void _CLG_LAUNCH_BOUND
_kernelExchangeOrders(INT* orders)
{
    UINT x = threadIdx.x; //0 to dm
    UINT y = threadIdx.y; //0 to dk
    __shared__ INT outOrders[CLinearAlgebraHelper::_kMaxSmallDim];

    if (y == orders[x])
    {
        //if 0 == orders[3], outOrders[3] = 0
        outOrders[y] = x;
    }

    __syncthreads();

    if (0 == x)
    {
        orders[y] = outOrders[y];
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelInitialE1Vector(CLGComplex* v, UINT dm, UINT x)
{
    UINT y = threadIdx.x;
    if (0 == y)
    {
        v[x * dm + y] = _make_cuComplex(F(1.0), F(0.0));
    }
    else
    {
        v[x * dm + y] = _make_cuComplex(F(0.0), F(0.0));
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCreateBackshiftProblem(
    const CLGComplex* __restrict__ triangular, const CLGComplex* __restrict__ eigenValue,
    UINT i, UINT iOrder, //the k row, iOrder = k - 1
    UINT dm,
    CLGComplex* resultMatrixR, CLGComplex* resultVecotr)
{
    UINT x = threadIdx.x; //0 -> iOrder - 1
    UINT y = threadIdx.y; //0 -> iOrder - 1

    if (y < x)
    {
        resultMatrixR[x * iOrder + y] = _make_cuComplex(F(0.0), F(0.0));
    }
    else if (y == x)
    {
        resultMatrixR[x * iOrder + y] = _cuCsubf(triangular[x * dm + y], eigenValue[i]);
    }
    else
    {
        resultMatrixR[x * iOrder + y] = triangular[x * dm + y];
    }

    if (0 == y)
    {
        resultVecotr[x].x = -triangular[x * dm + iOrder].x;
        resultVecotr[x].y = -triangular[x * dm + iOrder].y;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFinalNorm(
    const CLGComplex* __restrict__ triangular,
    UINT iOrder, //the k row, iOrder = k - 1
    UINT dm,
    CLGComplex* resultVecotr)
{
    UINT x = threadIdx.x;
    __shared__ Real fLength;
    if (1 == iOrder)
    {
        if (0 == x)
        {
            //norm of
            //- r[0,1]/(r[0,0] - r[1,1])
            CLGComplex a = _cuCdivf(triangular[1], _cuCsubf(triangular[dm + 1], triangular[0]));
            fLength = __div(F(1.0), _sqrt(__cuCabsSqf(a) + F(1.0)));
            resultVecotr[0] = a;
            resultVecotr[0].x = resultVecotr[0].x * fLength;
            resultVecotr[0].y = resultVecotr[0].y * fLength;
            resultVecotr[1] = _make_cuComplex(fLength, F(0.0));
        }
        else if (x > 1)
        {
            resultVecotr[x] = _make_cuComplex(F(0.0), F(0.0));
        }
    }
    else
    {
        if (0 == x)
        {
            fLength = F(0.0);
        }

        __syncthreads();

        if (x < iOrder)
        {
            atomicAdd(&fLength, __cuCabsSqf(resultVecotr[x]));
        }

        __syncthreads();

        if (0 == x)
        {
            fLength = __div(F(1.0), _sqrt(fLength + F(1.0)));
        }

        __syncthreads();

        if (x < iOrder)
        {
            resultVecotr[x].x = resultVecotr[x].x * fLength;
            resultVecotr[x].y = resultVecotr[x].y * fLength;
        }
        else if (x == iOrder)
        {
            resultVecotr[x] = _make_cuComplex(fLength, F(0.0));
        }
        else
        {
            resultVecotr[x] = _make_cuComplex(F(0.0), F(0.0));
        }
    }
}

void CLinearAlgebraHelper::UpperTriangularEigenVectors(
    const CLGComplex* upperTriangular, CLGComplex* outEigenValue, CLGComplex* outEigenVector,
    UINT dm, UINT dk, UBOOL bSmall)
{
    if (dm > m_uiDim || dk > dm)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d! or required eigen vector number larger than dimension!\n"), m_uiDim);
        return;
    }

    dim3 block(1, 1, 1);
    dim3 thread1(dm, dm, 1);
    dim3 thread2(dm, 1, 1);
    dim3 thread3(dm, dk, 1);

    if (bSmall)
    {
        _kernelSortEigenValues << <block, thread1 >> > (upperTriangular, outEigenValue, m_pDeviceFloatBuffer, m_pDeviceIntBuffer, dk, dm);
    }
    else
    {
        _kernelSortEigenValuesBig << <block, thread1 >> > (upperTriangular, outEigenValue, m_pDeviceFloatBuffer, m_pDeviceIntBuffer, dk, dm);
    }
    _kernelExchangeOrders << <block, thread3 >>> (m_pDeviceIntBuffer);
    INT orders[_kMaxSmallDim];
    checkCudaErrors(cudaMemcpy(orders, m_pDeviceIntBuffer, sizeof(INT) * dk, cudaMemcpyDeviceToHost));

    CLGComplex* tmpVector = m_pDeviceComplexBuffer1;
    STmpMatrix sTmpR = GetTmpMatrix();
    CLGComplex* tmpR = sTmpR.m_pMatrix;

    for (UINT i = 0; i < dk; ++i)
    {
        if (0 == orders[i])
        {
            //usually it is not 0
            _kernelInitialE1Vector << <block, thread2 >> > (outEigenVector, dm, i);
        }
        else if (1 == orders[i])
        {
            _kernelFinalNorm << <block, thread2 >> > (upperTriangular, 1, dm, tmpVector);
            checkCudaErrors(cudaMemcpy(outEigenVector + i * dm, tmpVector, sizeof(CLGComplex) * dm, cudaMemcpyDeviceToDevice));
        }
        else
        {
            //create back shift problem
            UINT toSolveDim = orders[i];
            //if it is the number 2 eigen-value, it is the 3rd eigen-value
            //when 3rd eigen-value, we need a 2x2 matrix.
            dim3 thread4(toSolveDim, toSolveDim, 1);
            _kernelCreateBackshiftProblem << <block, thread4 >> > (upperTriangular, outEigenValue, 
                i, toSolveDim, dm, tmpR, tmpVector);

            //solve back shift
            SolveY(tmpVector, tmpR, 1, toSolveDim);

            _kernelFinalNorm << <block, thread2 >> > (upperTriangular, toSolveDim, dm, tmpVector);

            checkCudaErrors(cudaMemcpy(outEigenVector + i * dm, tmpVector, sizeof(CLGComplex) * dm, cudaMemcpyDeviceToDevice));
        }
    }

    sTmpR.Free();
}

#pragma endregion

#pragma region Generalized Eigen Problem

void CLinearAlgebraHelper::GeneralizedEigenValueProblem(
    CLGComplex* A,
    CLGComplex* B,
    CLGComplex* outEigenValue,
    CLGComplex* outEigenVector,
    UINT dm, UINT dk, UBOOL bSmall,
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

    EigenValueProblem(B, outEigenValue, outEigenVector, dm, dk, bSmall, fEigenCrit, iMaxEigenIterate, fCrit, iMaxIterate);
}

#pragma endregion

#pragma region Givens

/**
* Left Given
*
* AX=B
* A'X=GAX=GB
* where A'[i-1, i] is zeroed.
* A is assumed to be a Henssenberg
*
* j from 0 to n-3
* i from n-1 to j+1
*
*
* left:
*
* h00* h10*   h00 h01  =  +  +
* -h10 h00    h10 h11     0  +
*
* right:
* h00 h01   h11 h10*   = +  +
* h10 h11  -h10 h11*     0  +
*
* A = dm x dmA, B = dm x 1
*
* thread.x = dmA - j
* thread.y = dm
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelLeftGivenHessenberg(UINT i, UINT j,
    CLGComplex* A, CLGComplex* g, UINT dm)
{
    //A is Henssenberg, so no need to calculate all 
    //UINT affectedDimA = dm - j;
    UINT x = threadIdx.x;
    
    __shared__ CLGComplex lineAi[CLinearAlgebraHelper::_kMaxSmallDim];
    __shared__ CLGComplex lineAj[CLinearAlgebraHelper::_kMaxSmallDim];
    __shared__ CLGComplex c0;
    __shared__ CLGComplex s0;
    __shared__ CLGComplex c0h;
    __shared__ CLGComplex s0h;

    lineAi[x] = A[(i - 1) * dm + j + x];
    lineAj[x] = A[i * dm + j + x];

    if (0 == x)
    {
        CLGComplex h00 = A[(i - 1) * dm + j];
        CLGComplex h10 = A[i * dm + j];
        Real fh00Addh10 = __cuCabsSqf(h00) + __cuCabsSqf(h10);
        if (fh00Addh10 > _CLG_FLT_MIN_)
        {
            Real fDemon = __div(F(1.0), _sqrt(fh00Addh10));
            c0 = cuCmulf_cr(h00, fDemon);
            s0 = cuCmulf_cr(h10, fDemon);
            c0h = _cuConjf(c0);
            s0h = _cuConjf(s0);
        }
        else
        {
            c0 = _make_cuComplex(F(1.0), F(0.0));
            s0 = _make_cuComplex(F(0.0), F(0.0));
            c0h = _make_cuComplex(F(1.0), F(0.0));
            s0h = _make_cuComplex(F(0.0), F(0.0));
        }

        //  c0h s0h
        //  -s0 c0
        CLGComplex g_im1 = g[i - 1];
        g[i - 1] = _cuCaddf(_cuCmulf(c0h, g[i - 1]), _cuCmulf(s0h, g[i]));
        g[i] = _cuCsubf(_cuCmulf(c0, g[i]), _cuCmulf(s0, g_im1));
    }

    __syncthreads();

    A[(i - 1) * dm + x + j] = _cuCaddf(_cuCmulf(c0h, lineAi[x]), _cuCmulf(s0h, lineAj[x]));
    A[i * dm + x + j] = _cuCsubf(_cuCmulf(c0, lineAj[x]), _cuCmulf(s0, lineAi[x]));
}

void CLinearAlgebraHelper::RotateHenssenberg(CLGComplex* H, CLGComplex* Y, UINT dm)
{
    dim3 block(1, 1, 1);
    for (UINT i = 0; i < dm; ++i)
    {
        dim3 thread(dm - i, 1, 1);
        _kernelLeftGivenHessenberg << <block, thread >> > (i + 1, i, H, Y, dm);
    }
}

#pragma endregion

#pragma region Host Functions

void CLinearAlgebraHelper::InitialZeroHost(CLGComplex* hostMatrix, UINT dx, UINT dy)
{
    if (dx > m_uiDim || dy > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpM = GetTmpMatrix();
    CLGComplex* tmpM = sTmpM.m_pMatrix;
    checkCudaErrors(cudaMemcpy(tmpM, hostMatrix, sizeof(CLGComplex) * dx * dy, cudaMemcpyHostToDevice));
    InitialZero(tmpM, dx, dy);
    checkCudaErrors(cudaMemcpy(hostMatrix, tmpM, sizeof(CLGComplex) * dx * dy, cudaMemcpyDeviceToHost));
    sTmpM.Free();
}

void CLinearAlgebraHelper::InitialOneHost(CLGComplex* hostMatrix, UINT dx)
{
    if (dx > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpM = GetTmpMatrix();
    CLGComplex* tmpM = sTmpM.m_pMatrix;
    checkCudaErrors(cudaMemcpy(tmpM, hostMatrix, sizeof(CLGComplex) * dx * dx, cudaMemcpyHostToDevice));
    InitialOne(tmpM, dx);
    checkCudaErrors(cudaMemcpy(hostMatrix, tmpM, sizeof(CLGComplex) * dx * dx, cudaMemcpyDeviceToHost));
    sTmpM.Free();
}

CLGComplex CLinearAlgebraHelper::Normal2Host(const CLGComplex* v, const CLGComplex* matrix, UINT dm)
{
    if (dm > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return _make_cuComplex(F(0.0), F(0.0));
    }

    STmpMatrix sTmpM = GetTmpMatrix();
    CLGComplex* tmpM = sTmpM.m_pMatrix;

    STmpMatrix sTmpV = GetTmpMatrix();
    CLGComplex* tmpV = sTmpV.m_pMatrix;

    CLGComplex res[1];
    res[0] = _make_cuComplex(F(0.0), F(0.0));
    checkCudaErrors(cudaMemcpy(m_pOneDeviceC, res, sizeof(CLGComplex), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(tmpM, matrix, sizeof(CLGComplex) * dm * dm, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(tmpV, v, sizeof(CLGComplex) * dm, cudaMemcpyHostToDevice));

    Normal2(tmpV, tmpM, m_pOneDeviceC, dm);

    checkCudaErrors(cudaMemcpy(res, m_pOneDeviceC, sizeof(CLGComplex), cudaMemcpyDeviceToHost));
    sTmpM.Free();
    sTmpV.Free();

    return res[0];
}

void CLinearAlgebraHelper::SmallMatrixMultHost(
    CLGComplex * hostRes,
    const CLGComplex* left,
    const CLGComplex* right,
    UINT dLeft, UINT dMid, UINT dRight,
    UBOOL bLeftDagger, UBOOL bRightDagger)
{
    if (dLeft > m_uiDim || dMid > m_uiDim || dRight > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpMRes = GetTmpMatrix();
    CLGComplex* tmpMRes = sTmpMRes.m_pMatrix;
    STmpMatrix sTmpMLeft = GetTmpMatrix();
    CLGComplex* tmpMLeft = sTmpMLeft.m_pMatrix;
    STmpMatrix sTmpMRight = GetTmpMatrix();
    CLGComplex* tmpMRight = sTmpMRight.m_pMatrix;

    checkCudaErrors(cudaMemcpy(tmpMLeft, left, sizeof(CLGComplex) * dLeft * dMid, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(tmpMRight, right, sizeof(CLGComplex) * dMid * dRight, cudaMemcpyHostToDevice));

    SmallMatrixMult(
        tmpMRes, tmpMLeft, tmpMRight,
        dLeft, dMid, dRight,
        bLeftDagger, bRightDagger);

    checkCudaErrors(cudaMemcpy(hostRes, tmpMRes, sizeof(CLGComplex) * dLeft * dRight, cudaMemcpyDeviceToHost));
    sTmpMRes.Free();
    sTmpMLeft.Free();
    sTmpMRight.Free();
}

void CLinearAlgebraHelper::BlockMatrixMultHost(
    CLGComplex * hostRes,
    const CLGComplex* left,
    const CLGComplex* right,
    UINT dDim, UINT uiStart, UINT uiEnd,
    UBOOL bLeft, UBOOL bLeftDagger, UBOOL bRightDagger)
{
    if (dDim > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpMRes = GetTmpMatrix();
    CLGComplex* tmpMRes = sTmpMRes.m_pMatrix;
    STmpMatrix sTmpMLeft = GetTmpMatrix();
    CLGComplex* tmpMLeft = sTmpMLeft.m_pMatrix;
    STmpMatrix sTmpMRight = GetTmpMatrix();
    CLGComplex* tmpMRight = sTmpMRight.m_pMatrix;

    checkCudaErrors(cudaMemcpy(tmpMLeft, left, sizeof(CLGComplex) * dDim * dDim, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(tmpMRight, right, sizeof(CLGComplex) * dDim * dDim, cudaMemcpyHostToDevice));

    BlockMatrixMult(
        tmpMRes, tmpMLeft, tmpMRight,
        dDim, uiStart, uiEnd,
        bLeft, bLeftDagger, bRightDagger);

    checkCudaErrors(cudaMemcpy(hostRes, tmpMRes, sizeof(CLGComplex) * dDim * dDim, cudaMemcpyDeviceToHost));
    sTmpMRes.Free();
    sTmpMLeft.Free();
    sTmpMRight.Free();
}

void CLinearAlgebraHelper::BlockCopyHost(CLGComplex* hostDest, const CLGComplex* hostSrc,
    UINT lengthX, UINT lengthY, UINT dimDest, UINT dimSrc)
{
    if (lengthX > m_uiDim || lengthY > m_uiDim || dimDest > m_uiDim || dimSrc > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpMRes = GetTmpMatrix();
    CLGComplex* tmpMRes = sTmpMRes.m_pMatrix;
    STmpMatrix sTmpMSrc = GetTmpMatrix();
    CLGComplex* tmpMSrc = sTmpMSrc.m_pMatrix;

    checkCudaErrors(cudaMemcpy(tmpMSrc, hostSrc, sizeof(CLGComplex) * dimSrc * lengthY, cudaMemcpyHostToDevice));
    //copy the unchanged elements
    checkCudaErrors(cudaMemcpy(tmpMRes, hostDest, sizeof(CLGComplex) * dimDest * lengthY, cudaMemcpyHostToDevice));

    BlockCopy(tmpMRes, tmpMSrc, lengthX, lengthY, dimDest, dimSrc);

    checkCudaErrors(cudaMemcpy(hostDest, tmpMRes, sizeof(CLGComplex) * dimDest * lengthY, cudaMemcpyDeviceToHost));

    sTmpMRes.Free();
    sTmpMSrc.Free();
}

void CLinearAlgebraHelper::QRFactorizationHost(CLGComplex* Q, CLGComplex* R, const CLGComplex* T, UINT uiDim)
{
    if (uiDim > m_uiDim)
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

    checkCudaErrors(cudaMemcpy(tmpT, T, sizeof(CLGComplex) * uiDim * uiDim, cudaMemcpyHostToDevice));

    QRFactorization(tmpQ, tmpR, tmpT, uiDim);

    checkCudaErrors(cudaMemcpy(Q, tmpQ, sizeof(CLGComplex) * uiDim * uiDim, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(R, tmpR, sizeof(CLGComplex) * uiDim * uiDim, cudaMemcpyDeviceToHost));

    sTmpQ.Free();
    sTmpR.Free();
    sTmpT.Free();
}

void CLinearAlgebraHelper::ThinQRFactorizationHost(CLGComplex* Q, CLGComplex* R, const CLGComplex* T, UINT dx, UINT dy)
{
    if (dx > m_uiDim || dy > m_uiDim)
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

    checkCudaErrors(cudaMemcpy(tmpT, T, sizeof(CLGComplex) * dx * dy, cudaMemcpyHostToDevice));

    ThinQRFactorization(tmpQ, tmpR, tmpT, dx, dy);

    checkCudaErrors(cudaMemcpy(Q, tmpQ, sizeof(CLGComplex) * dx * dy, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(R, tmpR, sizeof(CLGComplex) * dy * dy, cudaMemcpyDeviceToHost));

    sTmpQ.Free();
    sTmpR.Free();
    sTmpT.Free();
}

void CLinearAlgebraHelper::SolveYHost(CLGComplex* Y, const CLGComplex* R, UINT dk, UINT dx)
{
    if (dx > m_uiDim || dk > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpY = GetTmpMatrix();
    CLGComplex* tmpY = sTmpY.m_pMatrix;
    STmpMatrix sTmpR = GetTmpMatrix();
    CLGComplex* tmpR = sTmpR.m_pMatrix;

    checkCudaErrors(cudaMemcpy(tmpR, R, sizeof(CLGComplex) * dx * dx, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(tmpY, Y, sizeof(CLGComplex) * dx * dk, cudaMemcpyHostToDevice));

    SolveY(tmpY, tmpR, dk, dx);

    checkCudaErrors(cudaMemcpy(Y, tmpY, sizeof(CLGComplex) * dx * dk, cudaMemcpyDeviceToHost));

    sTmpY.Free();
    sTmpR.Free();
}

void CLinearAlgebraHelper::UpperTriangularEigenVectorsHost(
    const CLGComplex* upperTriangular, CLGComplex* outEigenValue, CLGComplex* outEigenVector,
    UINT dm, UINT dk, UBOOL bSmall)
{
    if (dm > m_uiDim || dk > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpT = GetTmpMatrix();
    CLGComplex* tmpT = sTmpT.m_pMatrix;
    STmpMatrix sTmpE = GetTmpMatrix();
    CLGComplex* tmpE = sTmpE.m_pMatrix;
    STmpMatrix sTmpV = GetTmpMatrix();
    CLGComplex* tmpV = sTmpV.m_pMatrix;

    checkCudaErrors(cudaMemcpy(tmpT, upperTriangular, sizeof(CLGComplex) * dm * dm, cudaMemcpyHostToDevice));

    UpperTriangularEigenVectors(tmpT, tmpE, tmpV, dm, dk, bSmall);

    checkCudaErrors(cudaMemcpy(outEigenValue, tmpE, sizeof(CLGComplex) * dk, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(outEigenVector, tmpV, sizeof(CLGComplex) * dm * dk, cudaMemcpyDeviceToHost));

    sTmpT.Free();
    sTmpE.Free();
    sTmpV.Free();
}

void CLinearAlgebraHelper::EigenValueProblemHost(CLGComplex* H, CLGComplex* outEigenValue, CLGComplex* outEigenVector,
    UINT dm, UINT dk, UBOOL bSmall, Real fEigenCrit, UINT iMaxEigenIter, Real fQRCrit, UINT iMaxIterate)
{
    if (dm > m_uiDim || dk > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpT = GetTmpMatrix();
    CLGComplex* tmpT = sTmpT.m_pMatrix;
    STmpMatrix sTmpE = GetTmpMatrix();
    CLGComplex* tmpE = sTmpE.m_pMatrix;
    STmpMatrix sTmpV = GetTmpMatrix();
    CLGComplex* tmpV = sTmpV.m_pMatrix;

    checkCudaErrors(cudaMemcpy(tmpT, H, sizeof(CLGComplex) * dm * dm, cudaMemcpyHostToDevice));

    EigenValueProblem(tmpT, tmpE, tmpV, dm, dk, bSmall, fEigenCrit, iMaxEigenIter, fQRCrit, iMaxIterate);

    checkCudaErrors(cudaMemcpy(outEigenValue, tmpE, sizeof(CLGComplex) * dk, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(outEigenVector, tmpV, sizeof(CLGComplex) * dm * dk, cudaMemcpyDeviceToHost));

    sTmpT.Free();
    sTmpE.Free();
    sTmpV.Free();
}

void CLinearAlgebraHelper::EigenValueProblemHessenbergHost(CLGComplex* H, CLGComplex* outEigenValue, CLGComplex* outEigenVector,
    UINT dm, UINT dk, UBOOL bSmall, Real fEigenCrit, UINT iMaxEigenIter, Real fQRCrit, UINT iMaxIterate)
{
    if (dm > m_uiDim || dk > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpT = GetTmpMatrix();
    CLGComplex* tmpT = sTmpT.m_pMatrix;
    STmpMatrix sTmpE = GetTmpMatrix();
    CLGComplex* tmpE = sTmpE.m_pMatrix;
    STmpMatrix sTmpV = GetTmpMatrix();
    CLGComplex* tmpV = sTmpV.m_pMatrix;

    checkCudaErrors(cudaMemcpy(tmpT, H, sizeof(CLGComplex) * dm * dm, cudaMemcpyHostToDevice));

    EigenValueProblemHessenberg(tmpT, tmpE, tmpV, dm, dk, bSmall, fEigenCrit, iMaxEigenIter, fQRCrit, iMaxIterate);

    checkCudaErrors(cudaMemcpy(outEigenValue, tmpE, sizeof(CLGComplex) * dk, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(outEigenVector, tmpV, sizeof(CLGComplex) * dm * dk, cudaMemcpyDeviceToHost));

    sTmpT.Free();
    sTmpE.Free();
    sTmpV.Free();
}

void CLinearAlgebraHelper::GeneralizedEigenValueProblemHost(
    CLGComplex* A, CLGComplex* B,
    CLGComplex* outEigenValue,
    CLGComplex* outEigenVector,
    UINT dm, UINT dk, UBOOL bSmall, Real fEigenCrit, UINT iMaxEigenIter, Real fQRCrit, UINT iMaxIterate)
{
    if (dm > m_uiDim || dk > m_uiDim)
    {
        appCrucial(_T("Cannot deal with matrix larger than %d!\n"), m_uiDim);
        return;
    }

    STmpMatrix sTmpA = GetTmpMatrix();
    CLGComplex* tmpA = sTmpA.m_pMatrix;
    STmpMatrix sTmpB = GetTmpMatrix();
    CLGComplex* tmpB = sTmpB.m_pMatrix;
    STmpMatrix sTmpE = GetTmpMatrix();
    CLGComplex* tmpE = sTmpE.m_pMatrix;
    STmpMatrix sTmpV = GetTmpMatrix();
    CLGComplex* tmpV = sTmpV.m_pMatrix;

    checkCudaErrors(cudaMemcpy(tmpA, A, sizeof(CLGComplex) * dm * dm, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(tmpB, B, sizeof(CLGComplex) * dm * dm, cudaMemcpyHostToDevice));

    GeneralizedEigenValueProblem(tmpA, tmpB, tmpE, tmpV, dm, dk, bSmall, fEigenCrit, iMaxEigenIter, fQRCrit, iMaxIterate);

    checkCudaErrors(cudaMemcpy(outEigenValue, tmpE, sizeof(CLGComplex) * dk, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(outEigenVector, tmpV, sizeof(CLGComplex) * dm * dk, cudaMemcpyDeviceToHost));

    sTmpA.Free();
    sTmpB.Free();
    sTmpE.Free();
    sTmpV.Free();
}

void CLinearAlgebraHelper::RotateHenssenbergHost(CLGComplex* H, CLGComplex* Ye1, UINT dmH)
{
    STmpMatrix sTmpH = GetTmpMatrix();
    CLGComplex* tmpH = sTmpH.m_pMatrix;
    STmpMatrix sTmpY = GetTmpMatrix();
    CLGComplex* tmpY = sTmpY.m_pMatrix;

    checkCudaErrors(cudaMemcpy(tmpH, H, sizeof(CLGComplex) * (dmH + 1) * dmH, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(tmpY, Ye1, sizeof(CLGComplex) * (dmH + 1), cudaMemcpyHostToDevice));

    RotateHenssenberg(tmpH, tmpY, dmH);

    checkCudaErrors(cudaMemcpy(H, tmpH, sizeof(CLGComplex) * (dmH + 1) * dmH, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(Ye1, tmpY, sizeof(CLGComplex) * (dmH + 1), cudaMemcpyDeviceToHost));

    sTmpH.Free();
    sTmpY.Free();
}

#pragma endregion

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
