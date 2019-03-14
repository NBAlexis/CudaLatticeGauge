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

__BEGIN_NAMESPACE

class CLGAPI CLinearAlgebraHelper
{
public:
    enum { _kMaxSmallDim = 32, };

    CLinearAlgebraHelper(UINT uiDim);
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
    static void InitialOne(CLGComplex* deviceMatrix, UINT dx);

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

    /**
    * Only implement the left upper corner copy.
    */
    static void BlockCopy(CLGComplex* deviceDest, const CLGComplex* deviceSrc,
        UINT lengthX, UINT lengthY, UINT dimDest, UINT dimSrc);

#pragma endregion

    void QRFactorization(CLGComplex* Q, CLGComplex* R, const CLGComplex* T, UINT uiDim);
    void ThinQRFactorization(CLGComplex* Q, CLGComplex* R, const CLGComplex* T, UINT dx, UINT dy);
    static void SolveY(CLGComplex* deviceY, const CLGComplex* deviceR, UINT dk, UINT dx);
    void EigenValueProblem(CLGComplex* H, CLGComplex* outEigenValue, CLGComplex* outEigenVector, 
        UINT dm, UINT dk, Real fEigenCrit = F(0.000000001), UINT iMaxEigenIter = 20, Real fQRCrit = F(0.0000000000001), UINT iMaxIterate = 200);

    void GeneralizedEigenValueProblem(
        CLGComplex* A, CLGComplex* B,
        CLGComplex* outEigenValue,
        CLGComplex* outEigenVector,
        UINT dm, UINT dk, Real fEigenCrit = F(0.000000001), UINT iMaxEigenIter = 20, Real fQRCrit = F(0.0000000000001), UINT iMaxIterate = 200);

    UINT m_uiDim;

protected:

    void Henssenberg(CLGComplex* T, UINT dx);
    void QRIterate(CLGComplex* T, UINT dx, Real fCrit = F(0.00000000001), UINT iCrit = 100);

    struct CLGAPI STmpMatrix
    {
        STmpMatrix() : m_pOwner(NULL), m_uiIndex(0), m_pMatrix(NULL), m_bUsing(FALSE) {}
        STmpMatrix(CLinearAlgebraHelper* pOwner, UINT uiIndex, CLGComplex *pMatrix) 
            : m_pOwner(pOwner), m_uiIndex(uiIndex), m_pMatrix(pMatrix), m_bUsing(FALSE) {}

        CLinearAlgebraHelper* m_pOwner;
        UINT m_uiIndex;
        CLGComplex * m_pMatrix;
        UBOOL m_bUsing;

        void Free()
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
    CLGComplex* m_pDeviceComplexBuffer1;
    CLGComplex* m_pDeviceComplexBuffer2;
};

__END_NAMESPACE

#endif //#ifndef _LINEARALGEBRAHELPER_H_


//=============================================================================
// END OF FILE
//=============================================================================
