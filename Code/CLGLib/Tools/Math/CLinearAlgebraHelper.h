//=============================================================================
// FILENAME : CLinearAlgebraHelper.h
// 
// DESCRIPTION:
// This is the helper class for small matrix
//
// REVISION:
//  [14/03/2019 nbale]
//=============================================================================

#ifndef _LINEARALGEBRAHELPER_H_
#define _LINEARALGEBRAHELPER_H_

#define _CLG_QR_Left_Top_Deflation 0
#define _CLG_IMPLICITE_QR_SHIFT 1
//Update eigen-value will not improve the iteration of eigen-vector
#define _CLG_QRIterate_Update_EigenValue 0
#define _CLG_QRIterate_MaxIterate (200)

#if _CLG_IMPLICITE_QR_SHIFT
#define _CLG_DEFAULT_QR_CRIT F(0.00000001)
#define _CLG_DEFAULT_QR_VECTOR_CRIT F(0.000001)
//If after 10 iteration, it is not converged, it is not about to converge
#if _CLG_DOUBLEFLOAT
#define _CLG_DEFAULT_QR_VECTOR_ITE 20
#else
#define _CLG_DEFAULT_QR_VECTOR_ITE 10
#endif
#else
#define _CLG_DEFAULT_QR_CRIT F(0.0000000001)
#define _CLG_DEFAULT_QR_VECTOR_CRIT F(0.000001)
#if _CLG_DOUBLEFLOAT
#define _CLG_DEFAULT_QR_VECTOR_ITE 20
#else
#define _CLG_DEFAULT_QR_VECTOR_ITE 10
#endif
#endif

__BEGIN_NAMESPACE

class CLGAPI CLinearAlgebraHelper
{
public:
    enum 
    {
#if _CLG_DEBUG
        _kMaxSmallDim = 22,
#else
        _kMaxSmallDime = 32,
#endif
        _kAllocateMatrixNumber = 7,
    };

    CLinearAlgebraHelper(UINT uiDim, UINT uiPreAllocate = _kAllocateMatrixNumber);
    ~CLinearAlgebraHelper();

    static void TestSmallMatrix();

#pragma region Matrix Mult

    /**
    * Print the matrix which has same format as Mathematica
    */
    static void PrintMatrix(const CLGComplex* hostMatrix, UINT dx, UINT dy);
    static void PrintDeviceMatrix(const CLGComplex* device, UINT dx, UINT dy)
    {
        CLGComplex tmp[_kMaxSmallDim * _kMaxSmallDim];
        checkCudaErrors(cudaMemcpy(tmp, device, sizeof(CLGComplex) * dx * dy, cudaMemcpyDeviceToHost));
        PrintMatrix(tmp, dx, dy);
    }

    static void InitialZero(CLGComplex* deviceMatrix, UINT dx, UINT dy);
    void InitialZeroHost(CLGComplex* hostMatrix, UINT dx, UINT dy);

    static void InitialOne(CLGComplex* deviceMatrix, UINT dx);
    void InitialOneHost(CLGComplex* hostMatrix, UINT dx);

    /**
    * Note: the dimmesions of left, right, middle, is after dagger!
    * After dagger Left = X*Y, Right = Y*Z
    * dleft = X, dmid = Y, dright = Z
    */
    static void SmallMatrixMult(
        CLGComplex * deviceRes, 
        const CLGComplex* left, 
        const CLGComplex* right, 
        UINT dLeft, UINT dMid, UINT dRight, 
        UBOOL bLeftDagger, UBOOL bRightDagger);

    void SmallMatrixMultHost(
        CLGComplex * hostRes,
        const CLGComplex* left,
        const CLGComplex* right,
        UINT dLeft, UINT dMid, UINT dRight,
        UBOOL bLeftDagger, UBOOL bRightDagger);

    /**
    * if left is 1 0 0
    *            0 U 0
    *            0 0 1
    *  bLeft = true, otherwise if right is block, bLeft = false
    *  Square matrix is assumed, so only one "dDim" instead of dLeft, dMid, dRight
    */
    static void BlockMatrixMult(
        CLGComplex * deviceRes, 
        const CLGComplex* left, 
        const CLGComplex* right, 
        UINT dDim, UINT uiStart, UINT uiEnd,
        UBOOL bLeft, UBOOL bLeftDagger, UBOOL bRightDagger);

    void BlockMatrixMultHost(
        CLGComplex * hostRes,
        const CLGComplex* left,
        const CLGComplex* right,
        UINT dDim, UINT uiStart, UINT uiEnd,
        UBOOL bLeft, UBOOL bLeftDagger, UBOOL bRightDagger);

    /**
    * Only implement the left upper corner copy.
    */
    static void BlockCopy(CLGComplex* deviceDest, const CLGComplex* deviceSrc,
        UINT lengthX, UINT lengthY, UINT dimDest, UINT dimSrc);

    void BlockCopyHost(CLGComplex* hostDest, const CLGComplex* hostSrc,
        UINT lengthX, UINT lengthY, UINT dimDest, UINT dimSrc);

    void Transpose(CLGComplex* deviceMatrix, UINT dx, UINT dy);
    void TransposeHost(CLGComplex* hostMatrix, UINT dx, UINT dy);
    void Dagger(CLGComplex* deviceMatrix, UINT dx, UINT dy);
    void DaggerHost(CLGComplex* hostMatrix, UINT dx, UINT dy);

    /**
    * v^+ Av
    * NOTE: deviceRes must be zeroed
    */
    void Normal2(const CLGComplex* v, const CLGComplex* matrix, CLGComplex* deviceRes, UINT dm) const;
    CLGComplex Normal2Host(const CLGComplex* v, const CLGComplex* matrix, UINT dm);

#pragma endregion

    /**
     * Square matrix QR Factorization
     * House holder
     */
    void QRFactorization(CLGComplex* Q, CLGComplex* R, const CLGComplex* T, UINT uiDim);

    /**
     * Not square matrix, QR Factorization
     */
    void ThinQRFactorization(CLGComplex* Q, CLGComplex* R, const CLGComplex* T, UINT dx, UINT dy);

    /**
     * Assume R is a up-traingular dx * dx matrix (length of row, length of xdim is dx, it can be dx * dx + n)
     * Calculate R^{-1} on Y, which is a dk * dx matrix (batched solve, put y1,y2,y3 as a 3 * dx matrix)
     */
    static void SolveY(CLGComplex* Y, const CLGComplex* R, UINT dk, UINT dx);

    void RotateHenssenberg(CLGComplex* H, CLGComplex* B, UINT dmH) const;
    void RotateHenssenberg(CLGComplex* H, CLGComplex* B, UINT dmX, UINT dmY) const;
    void RotateHenssenbergHost(CLGComplex* H, CLGComplex* B, UINT dmH);

    void QRFactorizationHost(CLGComplex* Q, CLGComplex* R, const CLGComplex* T, UINT uiDim);
    void ThinQRFactorizationHost(CLGComplex* Q, CLGComplex* R, const CLGComplex* T, UINT dx, UINT dy);
    void SolveYHost(CLGComplex* Y, const CLGComplex* R, UINT dk, UINT dx);

    /**
    * For deflation preconditioner, after the first step after GMRES, we already have a upper triangular.
    * Note that, for all eigen solver, do NOT support identical eigen-values
    */
    void UpperTriangularEigenVectors(
        const CLGComplex* upperTriangular, CLGComplex* outEigenValue, CLGComplex* outEigenVector, 
        UINT dm, UINT dk, UBOOL bSmall = TRUE);

    void UpperTriangularEigenVectorsHost(
        const CLGComplex* upperTriangular, CLGComplex* outEigenValue, CLGComplex* outEigenVector,
        UINT dm, UINT dk, UBOOL bSmall = TRUE);

    void EigenValueProblem(CLGComplex* H, CLGComplex* outEigenValue, CLGComplex* outEigenVector, 
        UINT dm, UINT dk, UBOOL bSmall = TRUE, Real fEigenCrit = _CLG_DEFAULT_QR_VECTOR_CRIT, UINT iMaxEigenIter = _CLG_DEFAULT_QR_VECTOR_ITE, Real fQRCrit = _CLG_DEFAULT_QR_CRIT, UINT iMaxIterate = _CLG_QRIterate_MaxIterate);

    /**
    * For H is already Hessenberg
    */
    void EigenValueProblemHessenberg(CLGComplex* H, CLGComplex* outEigenValue, CLGComplex* outEigenVector,
        UINT dm, UINT dk, UBOOL bSmall = TRUE, Real fEigenCrit = _CLG_DEFAULT_QR_VECTOR_CRIT, UINT iMaxEigenIter = _CLG_DEFAULT_QR_VECTOR_ITE, Real fQRCrit = _CLG_DEFAULT_QR_CRIT, UINT iMaxIterate = _CLG_QRIterate_MaxIterate);

    void EigenValueProblemHost(CLGComplex* H, CLGComplex* outEigenValue, CLGComplex* outEigenVector,
        UINT dm, UINT dk, UBOOL bSmall = TRUE, Real fEigenCrit = _CLG_DEFAULT_QR_VECTOR_CRIT, UINT iMaxEigenIter = _CLG_DEFAULT_QR_VECTOR_ITE, Real fQRCrit = _CLG_DEFAULT_QR_CRIT, UINT iMaxIterate = _CLG_QRIterate_MaxIterate);

    void EigenValueProblemHessenbergHost(CLGComplex* H, CLGComplex* outEigenValue, CLGComplex* outEigenVector,
        UINT dm, UINT dk, UBOOL bSmall = TRUE, Real fEigenCrit = _CLG_DEFAULT_QR_VECTOR_CRIT, UINT iMaxEigenIter = _CLG_DEFAULT_QR_VECTOR_ITE, Real fQRCrit = _CLG_DEFAULT_QR_CRIT, UINT iMaxIterate = _CLG_QRIterate_MaxIterate);

    void GeneralizedEigenValueProblem(
        CLGComplex* A, CLGComplex* B,
        CLGComplex* outEigenValue,
        CLGComplex* outEigenVector,
        UINT dm, UINT dk, UBOOL bSmall = TRUE, Real fEigenCrit = _CLG_DEFAULT_QR_VECTOR_CRIT, UINT iMaxEigenIter = _CLG_DEFAULT_QR_VECTOR_ITE, Real fQRCrit = _CLG_DEFAULT_QR_CRIT, UINT iMaxIterate = _CLG_QRIterate_MaxIterate);

    void GeneralizedEigenValueProblemHost(
        CLGComplex* A, CLGComplex* B,
        CLGComplex* outEigenValue,
        CLGComplex* outEigenVector,
        UINT dm, UINT dk, UBOOL bSmall = TRUE, Real fEigenCrit = _CLG_DEFAULT_QR_VECTOR_CRIT, UINT iMaxEigenIter = _CLG_DEFAULT_QR_VECTOR_ITE, Real fQRCrit = _CLG_DEFAULT_QR_CRIT, UINT iMaxIterate = _CLG_QRIterate_MaxIterate);

    UINT GetMaxDim() const { return m_uiDim; }

protected:

    UINT m_uiDim;

    /**
     * Transform matrix T (dx * dx) matrix, to Henssenberg form, using house holding
     */
    void Henssenberg(CLGComplex* T, UINT dx);

    /**
    * Because of round off error and divide by a small number
    * It is not the smaller the better
    * About 1E-8 sometimes produce better result than 1E-10 
    *
    * This QRIterate is deprecated, use Double shift QRIterate instead.
    *
    * T is assumed to be Henssenberg
    */
    void QRIterate(CLGComplex* T, UINT dx, Real fCrit = F(0.0000000001), UINT iCrit = _CLG_QRIterate_MaxIterate);

    /**
    * Double shift QRIterate
    * T is assumed to be Henssenberg
    */
    void FrancisQRIterate(CLGComplex* T, UINT dx, Real fCrit = F(0.00000001), UINT iCrit = _CLG_QRIterate_MaxIterate);
    void FrancisQRIterateBlock(CLGComplex* T, CLGComplex* tmpXYZ, UINT uiBlockDim);

    struct CLGAPI STmpMatrix
    {
        STmpMatrix() : m_pOwner(NULL), m_uiIndex(0), m_pMatrix(NULL), m_bUsing(FALSE) {}
        STmpMatrix(CLinearAlgebraHelper* pOwner, UINT uiIndex, CLGComplex *pMatrix) 
            : m_pOwner(pOwner), m_uiIndex(uiIndex), m_pMatrix(pMatrix), m_bUsing(FALSE) {}

        CLinearAlgebraHelper* m_pOwner;
        UINT m_uiIndex;
        CLGComplex * m_pMatrix;
        UBOOL m_bUsing;

        void Free() const
        {
            m_pOwner->ReleaseTmpMatrix(m_uiIndex);
        }
    };

    STmpMatrix GetTmpMatrix()
    {
        for (INT i = 0; i < m_lstTmpMatrix.Num(); ++i)
        {
            if (!m_lstTmpMatrix[i].m_bUsing)
            {
                m_lstTmpMatrix[i].m_bUsing = TRUE;
                return m_lstTmpMatrix[i];
            }
        }

        CLGComplex * newM = NULL;
        checkCudaErrors(cudaMalloc((void**)&newM, sizeof(CLGComplex) * m_uiDim * m_uiDim));
        STmpMatrix newone(this, m_lstTmpMatrix.Num(), newM);
        newone.m_bUsing = TRUE;
        m_lstTmpMatrix.AddItem(newone);
        appDetailed(_T("Note: Allocating new matrix!!\n"));
        return newone;
    }

    void AddTempMatrix(UINT uiNum)
    {
        for (UINT i = 0; i < uiNum; ++i)
        {
            CLGComplex * newM = NULL;
            checkCudaErrors(cudaMalloc((void**)&newM, sizeof(CLGComplex) * m_uiDim * m_uiDim));
            STmpMatrix newone(this, m_lstTmpMatrix.Num(), newM);
            newone.m_bUsing = FALSE;
            m_lstTmpMatrix.AddItem(newone);
        }
    }

    void ReleaseTmpMatrix(UINT uiIndex) { m_lstTmpMatrix[uiIndex].m_bUsing = FALSE; }

    TArray<STmpMatrix> m_lstTmpMatrix;
    //length is uiDim
    INT* m_pDeviceIntBuffer;
    Real* m_pDeviceFloatBuffer;
    CLGComplex* m_pDeviceTmpHouseHolder;
    CLGComplex* m_pDeviceComplexBuffer1;
    CLGComplex* m_pDeviceComplexBuffer2;
    CLGComplex* m_pOneDeviceC;
};

__END_NAMESPACE

#endif //#ifndef _LINEARALGEBRAHELPER_H_


//=============================================================================
// END OF FILE
//=============================================================================
