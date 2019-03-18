//=============================================================================
// FILENAME : CSolverGCRODR.h
// 
// DESCRIPTION:
// This is the class for Sparse Linear Algebra solves.
//
// REVISION:
//  [03/15/2019 nbale]
//=============================================================================

#ifndef _CSOLVERGCRODR_H_
#define _CSOLVERGCRODR_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(EEigenDeflationType,
    EEDT_REV,
    EEDT_SVD,
    );

__CLG_REGISTER_HELPER_HEADER(CSLASolverGCRODR)

class CLGAPI CSLASolverGCRODR : public CSLASolver
{
    __CLGDECLARE_CLASS(CSLASolverGCRODR)

public:

    enum { _kMinKDim = 2, _kMaxDimDR = CLinearAlgebraHelper::_kMaxSmallDim - 1, };

    CSLASolverGCRODR();
    ~CSLASolverGCRODR();

    virtual void Configurate(const CParameters& param);
    virtual void AllocateBuffers(const CField* pField);
    virtual void ReleaseBuffers();

    virtual void OnNewTrajectory() { m_bHasYk = FALSE; }
    virtual UBOOL Solve(
        CField* pFieldX, 
        const CField* pFieldB, 
        const CFieldGauge* pGaugeFeild, 
        EFieldOperator uiM, 
        const CField* pStart = NULL);

protected:

    /**
    * x = x +- AB+ x
    */
    void AxpyABdagger(class CField* x, const TArray<class CField*>& lstA, const TArray<class CField*>& lstB, UBOOL bPlusOrMinus, UINT uiK);

    /**
    * Assume m >= k
    * V=(vk[0], ... , vk[k - 1])
    * W=(vk[0], ... , vk[k - 1], vmk[0], ..., vmk[m-k-1])
    *
    * V(v1,v2,...,vk) = W(w1,w2,...,wm) (m11, ..., m1k)
    *                                   (..., ..., ...)
    *                                   (mm1, ..., mmk)
    *
    *  VectorMultiplyMatrix is for Yk = V Pk, uiMDim = m and Ck = W Q, uiMDim = m + 1
    */
    void VectorMultiplyMatrix(
        //TArray<class CField*>& resultV, 
        const TArray<class CField*>& lstYorC, 
        //const TArray<class CField*>& lstVmk, 
        CLGComplex* hostMatrix, UINT uiMDim);
    //void MatrixMultiplyVector(TArray<class CField*>& resultV, const TArray<class CField*>& lstVk, const TArray<class CField*>& lstVmk, CLGComplex* hostMatrix);

    // YR-1
    void FieldSolveY(TArray<class CField*>& resultY, const CLGComplex* R, UINT uiDim);

    void QRFactorAY(const class CFieldGauge* pGaugeField, EFieldOperator uiM);

    void FindPk1();

    void FindPk2();

    void GenerateCUFirstTime(CField* pX, const CField* pFieldB, const CFieldGauge* pGaugeField, EFieldOperator uiM);
    void GenerateCU(UBOOL bUpdateCk, UBOOL bJustAfterGMRES);
    void NormUkAndSetD();

    void GetPooledFields(const class CField* pFieldB);
    void ReleasePooledFields();

    /**
    * m_pHostHmGm set
    */
    void FirstTimeGMERESSolve(CField* pFieldX, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM);

    class CLinearAlgebraHelper* m_pHelper;

    UINT m_uiReStart;
    UINT m_uiMDim;
    UINT m_uiKDim;
    Real m_fAccuracy;
    UBOOL m_bHasYk;

    EEigenDeflationType m_eDeflationType;

    //numer = m-k+1
    TArray<class CField*> m_lstV;
    //number = k
    TArray<class CField*> m_lstC;
    //number = k
    TArray<class CField*> m_lstU;
    //number = k
    TArray<class CField*> m_lstTmp;

    class CField* GetV(UINT uiIndex) 
    {
        if (uiIndex < m_uiKDim)
        {
            return m_lstU[uiIndex];
        }
        return m_lstV[uiIndex - m_uiKDim];
    }

    class CField* GetW(UINT uiIndex)
    {
        if (uiIndex < m_uiKDim)
        {
            return m_lstC[uiIndex];
        }
        return m_lstV[uiIndex - m_uiKDim];
    }

    CLGComplex* m_pDeviceHm;
    CLGComplex* m_pDeviceEigenValue;
    CLGComplex* m_pDevicePk;
    CLGComplex* m_pDeviceHmGm;
    
    CLGComplex* m_pDeviceALeft;
    CLGComplex* m_pDeviceA;
    CLGComplex* m_pDeviceB;

    //(m+1)xm
    CLGComplex* m_pHostHmGm;
    CLGComplex* m_pHostHmGmToRotate;
    //(m+1)x1
    CLGComplex* m_pHostY;

    CLGComplex* m_pHostALeft;
    CLGComplex* m_pHostB;
    CLGComplex* m_pHostPk;
    CLGComplex* m_pHostTmpQ;
    CLGComplex* m_pHostTmpR;
    //(m+1)xk
    CLGComplex* m_pHostTmpGPk;

    Real m_fBeta;
    Real m_fDiviation;
};

__END_NAMESPACE

#endif //#ifndef _CSOLVERGCRODR_H_

//=============================================================================
// END OF FILE
//=============================================================================