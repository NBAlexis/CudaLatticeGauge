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

    EFieldType GetFieldType() const override
    {
        return EFT_FermionWilsonSquareSU3;
    }

    void InitialField(EFieldInitialType eInitialType) override;
    void InitialFieldWithFile(const CCString&, EFieldFileType) override;
    void InitialWithByte(BYTE* byData) override;
    void InitialOtherParameters(CParameters& params) override;
    void DebugPrintMe() const override;

    void Zero() override { InitialField(EFIT_Zero); }

    void Identity() override
    { appCrucial(_T("Not supported for CFermionWilsonSquareSU3!")); }

    void Dagger() override;

    //This is Axpy(1.0f, x)
    void AxpyPlus(const CField* x) override;
    void AxpyMinus(const CField* x) override;
    void Axpy(Real a, const CField* x) override;
    void Axpy(const CLGComplex& a, const CField* x) override;
    void Mul(const CField* other, UBOOL bDagger = TRUE) override;
    void ScalarMultply(const CLGComplex& a) override;
    void ScalarMultply(Real a) override;
#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex Dot(const CField* other) const override;
#else
    CLGComplex Dot(const CField* other) const override;
#endif

    //=================================
    //It is tested, although, the DEBUG Mode, this is faster
    //But the RELEASE Mode, the above version is faster.
    //void AxpyPlus1(const CField* x);
    //void AxpyMinus1(const CField* x);
    //void Axpy1(Real a, const CField* x);
    //void Axpy1(const CLGComplex& a, const CField* x);
    //void ScalarMultply1(const CLGComplex& a);
    //void ScalarMultply1(Real a);
    //CLGComplex Dot1(const CField* other) const;


protected:

    //pGauge must be gauge SU3
    void DS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DDS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DDdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;

    void DWithMassS(const CField* , Real , EOperatorCoefficientType , Real , Real ) override
    {
        appCrucial(_T("Not supported for Wilson direct fermion!\n"));
    }

    void DdaggerWithMassS(const CField* , Real , EOperatorCoefficientType , Real , Real ) override
    {
        appCrucial(_T("Not supported for Wilson direct fermion!\n"));
    }

    void DDdaggerWithMassS(const CField* , Real , EOperatorCoefficientType , Real , Real ) override
    {
        appCrucial(_T("Not supported for Wilson direct fermion!\n"));
    }

    void DDWithMassS(const CField*, Real, EOperatorCoefficientType, Real, Real) override
    {
        appCrucial(_T("Not supported for Wilson direct fermion!\n"));
    }

    void PrepareForHMCS(const CFieldGauge* pGauge) override;
    UBOOL CalculateForceS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const override;

public:

    void ApplyGamma(EGammaMatrix eGamma) override;
    TArray<CFieldFermion*> GetSourcesAtSiteFromPool(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const SSmallInt4& site) const override;

    void InitialAsSource(const SFermionSource& sourceData) override;
    BYTE* CopyDataOut(UINT &uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;
    
    CCString GetInfos(const CCString &tab) const override;

    void SetKai(Real fKai);
    UINT TestGamma5Hermitian(const CFieldGauge* pGauge, UBOOL bTestGamma5 = FALSE) const;

    deviceWilsonVectorSU3* m_pDeviceData;

    _GetData

protected:

    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, 
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const override;

    Real m_fKai;

    //Not using, this is used in "Dot1" which create a thread for each element of a Wilson vector
    //In Debug, it is faster, in Release, it is slower, so not using.
    Real* m_tmpBuffer2;

};


class CLGAPI CFieldMatrixOperationWilsonSquareSU3 : public CFieldMatrixOperation
{
public:
    CFieldMatrixOperationWilsonSquareSU3();
    ~CFieldMatrixOperationWilsonSquareSU3();
    
    //real left = (res,left)
    void VectorMultiplyMatrix(TArray<CField*>& res, const TArray<CField*>& left, const CLGComplex* deviceMatrix, UINT uiDimX, UINT uiDimY) override;

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