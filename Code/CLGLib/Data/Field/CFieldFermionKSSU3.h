//=============================================================================
// FILENAME : CFieldFermionKSSU3.h
// 
// DESCRIPTION:
// This is the class for Kogut-Susskind staggered fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [12/08/2019 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONKSSU3_H_
#define _CFIELDFERMIONKSSU3_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3)

class CLGAPI CFieldFermionKSSU3 : public CFieldFermion
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3)

public:

    CFieldFermionKSSU3();
    ~CFieldFermionKSSU3();

    EFieldType GetFieldType() const override
    {
        return EFT_FermionStaggeredSU3;
    }

    void InitialField(EFieldInitialType eInitialType) override;
    void InitialFieldWithFile(const CCString&, EFieldFileType) override;
    void InitialWithByte(BYTE* byData) override;
    void InitialOtherParameters(CParameters& params) override;
    void DebugPrintMe() const override;

    void Zero() override { InitialField(EFIT_Zero); }
    void Dagger() override;

    void Identity() override
    {
        appCrucial(_T("Not supported for CFermionWilsonSquareSU3!"));
    }

    //This is Axpy(1.0f, x)
    void AxpyPlus(const CField* x) override;
    void AxpyMinus(const CField* x) override;
    void Axpy(Real a, const CField* x) override;
    void Axpy(const CLGComplex& a, const CField* x) override;
    void ScalarMultply(const CLGComplex& a) override;
    void ScalarMultply(Real a) override;
    CLGComplex Dot(const CField* other) const override;

    //pGauge must be gauge SU3
    //These are for Sparse linear algebra
    void D(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void Ddagger(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DDdagger(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    UBOOL InverseD(const CField* pGauge) override;
    UBOOL InverseDdagger(const CField* pGauge) override;
    UBOOL InverseDDdagger(const CField* pGauge) override;
    void ApplyGamma(EGammaMatrix eGamma) override;

    //These are truely D or InverseD etc.

    /**
     * Use to calculate action, it is (D^+D)^{-1/4}
     */
    virtual void D_MD(const CField* pGauge);
    virtual void D0(const CField* pGauge);
    virtual void D_MC(const CField* pGauge);

    void PrepareForHMC(const CFieldGauge* pGauge) override;
    UBOOL CalculateForce(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const override;

    void InitialAsSource(const SFermionSource& sourceData) override;
    void SaveToFile(const CCString& fileName) const override;
    BYTE* CopyDataOut(UINT& uiSize) const override;
    TArray<CFieldFermion*> GetSourcesAtSiteFromPool(const class CFieldGauge* pGauge, const SSmallInt4& site) const override;
    CCString GetInfos(const CCString& tab) const override;

    void SetMass(Real f2am);

    //============================
    //Override these two functions for KS
    virtual void DerivateD0(void* pForce, const void* pGaugeBuffer) const;
    virtual void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const;
    //============================

    /**
     * Do not override me
     */
    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override
    {
        DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, m_f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);
    }

    /**
     * DerivateDOperator not implemented for KS
     */
    void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const override;

    deviceSU3Vector* m_pDeviceData;

protected:

    Real m_f2am;

    // r(x) = x^{1/4} use to prepare for Nf=2
    // r(x) = x^{3/8} use as s quark for Nf=2+1
    // r(x) = (x+dm/x)^{-1/4} use as u,d quark for Nf=2+1
    CRatinalApproximation m_rMC;

    // r(x) = x^{-1/2} use to calculate force and action for Nf=2
    // r(x) = x^{-3/4} use to s quark for Nf=2+1 for Nf=2
    // r(x) = (x+dm/x)^{1/2} use as u,d quark for Nf=2+1
    CRatinalApproximation m_rMD;


    //phi _i and Dst0 phi _i
    deviceSU3Vector** m_pRationalFieldPointers;

    Real* m_pMDNumerator;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3_H_

//=============================================================================
// END OF FILE
//=============================================================================