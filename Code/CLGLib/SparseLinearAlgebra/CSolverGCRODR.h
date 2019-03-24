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
    EEDT_HEV,
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

    virtual UBOOL Solve(
        CField* pFieldX, 
        const CField* pFieldB, 
        const CFieldGauge* pGaugeFeild, 
        EFieldOperator uiM, 
        ESolverPhase ePhase = ESP_Once,
        const CField* pStart = NULL);

protected:

    void FieldSolveY(TArray<class CField*>& resultY, const CLGComplex* R, UINT uiDim);

    void QRFactorAY(const class CFieldGauge* pGaugeField, EFieldOperator uiM);

    void FindPk1();

    void FindPk2();

    virtual void GenerateCUFirstTime(CField* pX, CField* pR, const CField* pFieldB, const CFieldGauge* pGaugeField, EFieldOperator uiM);
    void GenerateCU(UBOOL bUpdateCk, UBOOL bJustAfterGMRES);
    void NormUkAndSetD();
    void OrthognalXR(CField* pX, CField* pR, CField* pTmp);

    void GetPooledFields(const class CField* pFieldB);
    void ReleasePooledFields();

    /**
    * m_pHostHmGm set
    */
    void FirstTimeGMERESSolve(CField* pFieldX, CField* pR, const CField* pFieldB, const CFieldGauge* pGaugeFeild, EFieldOperator uiM);

    class CLinearAlgebraHelper* m_pHelper;

    UINT m_uiReStart;
    UINT m_uiMDim;
    UINT m_uiKDim;
    UINT m_uiRecalcuateR;
    Real m_fAccuracy;

    EEigenDeflationType m_eDeflationType;

    //numer = m-k+1
    TArray<class CField*> m_lstV;
    //number = k
    TArray<class CField*> m_lstC;
    //number = k
    TArray<class CField*> m_lstU;

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
    CLGComplex* m_pDeviceTmpQ;

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

    CLGComplex* m_pHostZeroMatrix;

    Real m_fBeta;
    Real m_fDiviation;
    CLGComplex m_cLastDiviation;

    class CFieldMatrixOperation* m_pFieldMatrix;
};

__END_NAMESPACE

#endif //#ifndef _CSOLVERGCRODR_H_

//=============================================================================
// END OF FILE
//=============================================================================