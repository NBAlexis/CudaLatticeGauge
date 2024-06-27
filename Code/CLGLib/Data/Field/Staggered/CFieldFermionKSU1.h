//=============================================================================
// FILENAME : CFieldFermionKSU1.h
// 
// DESCRIPTION:
// This is the class for Kogut-Susskind staggered fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [10/01/2021 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONKSU1_H_
#define _CFIELDFERMIONKSU1_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSU1)

class CLGAPI CFieldFermionKSU1 : public CFieldFermionKS
{
    __CLGDECLARE_FIELD(CFieldFermionKSU1)

public:

    CFieldFermionKSU1();
    ~CFieldFermionKSU1();

    EFieldType GetFieldType() const override
    {
        return EFT_FermionStaggeredU1;
    }

    void InitialField(EFieldInitialType eInitialType) override;
    void InitialFieldWithFile(const CCString&, EFieldFileType) override;
    void InitialWithByte(BYTE* byData) override;
    void InitialOtherParameters(CParameters& params) override;
    void DebugPrintMe() const override;

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

protected:

    //pGauge must be gauge SU3
    //These are for Sparse linear algebra
    void DS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DDS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DDdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DdaggerWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DDWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DDdaggerWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;

public:

    void ApplyGamma(EGammaMatrix eGamma) override;
    TArray<CFieldFermion*> GetSourcesAtSiteFromPool(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const SSmallInt4& site) const override;
    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) override;

protected:

    void ApplyGammaKSS(const CFieldGauge* pGauge, EGammaMatrix eGamma) override;

    //================= test anti-hermitian =========
    UINT TestAntiHermitianS(const CFieldGauge* pGauge) const override;

    //These are truely D or InverseD etc.

    /**
     * Use to calculate action, it is (D^+D)^{-1/4}
     */
    void D0S(const CField* pGauge) override;
    UBOOL CalculateForceS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const override;

public:

    //For test only
    void PrepareForHMCOnlyRandomize() override;
    void PrepareForHMCNotRandomize(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson) override;

    void InitialAsSource(const SFermionSource& sourceData) override;
    BYTE* CopyDataOut(UINT& uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;

protected:

    //============================
    //Override these two functions for KS
    void DerivateD0(void* pForce, const void* pGaugeBuffer) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    //============================

public:

    /**
     * Do not override me
     */
    CLGComplex* m_pDeviceData;

    _GetData

    #pragma region Help functions to implement higher orders

    void OnlyMass(void* pTarget, Real f2am, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) override;

protected:

    void OneLinkS(const void* pGuage, BYTE byGaugeFieldId, void* pTarget, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx, 
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) override;
    void OneLinkForceS(const void* pGuage, BYTE byGaugeFieldId, void* pForce, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx) const override;

    #pragma endregion

protected:

    //phi _i and Dst0 phi _i
    CLGComplex** m_pRationalFieldPointers;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSU1_H_

//=============================================================================
// END OF FILE
//=============================================================================