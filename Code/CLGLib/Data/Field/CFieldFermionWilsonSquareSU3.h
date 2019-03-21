//=============================================================================
// FILENAME : CFieldFermionWilson.h
// 
// DESCRIPTION:
// This is the class for pseudofermion for Wilson fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [12/25/2018 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONWILSONSQUARESU3_H_
#define _CFIELDFERMIONWILSONSQUARESU3_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareSU3)

class CLGAPI CFieldFermionWilsonSquareSU3 : public CFieldFermion
{
    __CLGDECLARE_FIELD(CFieldFermionWilsonSquareSU3)

public:

    CFieldFermionWilsonSquareSU3();
    ~CFieldFermionWilsonSquareSU3();

    virtual EFieldType GetFieldType() const 
    {
        return EFT_FermionWilsonSquareSU3;
    }

    virtual void InitialField(EFieldInitialType eInitialType);
    virtual void InitialFieldWithFile(const CCString&, EFieldFileType);
    virtual void InitialWithByte(BYTE* byData);
    virtual void InitialOtherParameters(CParameters& params);
    virtual void DebugPrintMe() const;

    virtual void Zero() { InitialField(EFIT_Zero); }
    virtual void Indentity() { appCrucial(_T("Not supported for CFermionWilsonSquareSU3!")); }

    //This is Axpy(1.0f, x)
    virtual void AxpyPlus(const CField* x);
    virtual void AxpyMinus(const CField* x);
    virtual void Axpy(Real a, const CField* x);
    virtual void Axpy(const CLGComplex& a, const CField* x);
    virtual void ScalarMultply(const CLGComplex& a);
    virtual void ScalarMultply(Real a);
    virtual CLGComplex Dot(const CField* other) const;

    //=================================
    //It is tested, although, the DEBUG Mode, this is faster
    //But the RELEASE Mode, the above version is faster.
    void AxpyPlus1(const CField* x);
    void AxpyMinus1(const CField* x);
    void Axpy1(Real a, const CField* x);
    void Axpy1(const CLGComplex& a, const CField* x);
    void ScalarMultply1(const CLGComplex& a);
    void ScalarMultply1(Real a);
    CLGComplex Dot1(const CField* other) const;

    //pGauge must be gauge SU3
    virtual void D(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0));
    virtual void Ddagger(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0));
    virtual void DDdagger(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0));
    virtual UBOOL InverseD(const CField* pGauge);
    virtual UBOOL InverseDdagger(const CField* pGauge);
    virtual UBOOL InverseDDdagger(const CField* pGauge);
    virtual void ApplyGamma(EGammaMatrix eGamma);
    virtual void PrepareForHMC(const CFieldGauge* pGauge);

    virtual UBOOL CalculateForce(const CFieldGauge* pGauge, CFieldGauge* pForce, CFieldGauge* pCachedForce) const;
    virtual void InitialAsSource(const SFermionSource& sourceData);
    virtual void SaveToFile(const CCString &fileName) const;
    virtual BYTE* CopyDataOut(UINT &uiSize) const;
    virtual CCString GetInfos(const CCString &tab) const;

    void SetKai(Real fKai);

    deviceWilsonVectorSU3 * m_pDeviceData;

protected:

    Real m_fKai;
    Real* m_tmpBuffer2;

};

class CLGAPI CFieldMatrixOperationWilsonSquareSU3 : public CFieldMatrixOperation
{
public:
    CFieldMatrixOperationWilsonSquareSU3();
    ~CFieldMatrixOperationWilsonSquareSU3();
    
    //real left = (res,left)
    virtual void VectorMultiplyMatrix(TArray<CField*>& res, const TArray<CField*>& left, const CLGComplex* deviceMatrix, UINT uiDimX, UINT uiDimY);

    deviceWilsonVectorSU3** m_pResBuffer;
    deviceWilsonVectorSU3** m_pLeftBuffer;
    deviceWilsonVectorSU3** m_pHostResBuffer;
    deviceWilsonVectorSU3** m_pHostLeftBuffer;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3_H_

//=============================================================================
// END OF FILE
//=============================================================================