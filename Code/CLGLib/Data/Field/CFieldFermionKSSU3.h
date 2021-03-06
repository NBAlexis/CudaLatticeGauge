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

    enum { _kKSLinkLength = 3 };

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
#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex Dot(const CField* other) const override;
#else
    CLGComplex Dot(const CField* other) const override;
#endif

    //pGauge must be gauge SU3
    //These are for Sparse linear algebra
    void D(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void Ddagger(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DD(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DDdagger(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DWithMass(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DdaggerWithMass(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DDWithMass(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DDdaggerWithMass(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    UBOOL InverseD(const CField* pGauge) override;
    UBOOL InverseDdagger(const CField* pGauge) override;
    UBOOL InverseDD(const CField* pGauge) override;
    UBOOL InverseDDdagger(const CField* pGauge) override;
    void ApplyGamma(EGammaMatrix eGamma) override;
    virtual void ApplyGammaKS(const CFieldGauge* pGauge, EGammaMatrix eGamma);

    //================= test anti-hermitian =========
    UINT TestAntiHermitian(const CFieldGauge* pGauge) const;

    //These are truely D or InverseD etc.

    /**
     * Use to calculate action, it is (D^+D)^{-1/4}
     */
    virtual void D_MD(const CField* pGauge);
    virtual void D0(const CField* pGauge);
    virtual void D_MC(const CField* pGauge);

    void PrepareForHMC(const CFieldGauge* pGauge) override;
    UBOOL CalculateForce(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const override;

    //For test only
    void PrepareForHMCOnlyRandomize();
    void PrepareForHMCNotRandomize(const CFieldGauge* pGauge);

    void InitialAsSource(const SFermionSource& sourceData) override;
    BYTE* CopyDataOut(UINT& uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;
    TArray<CFieldFermion*> GetSourcesAtSiteFromPool(const class CFieldGauge* pGauge, const SSmallInt4& site) const override;
    CCString GetInfos(const CCString& tab) const override;

    void SetMass(Real f2am);
    Real GetMass() const { return m_f2am; }

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

    #pragma region Help functions to implement higher orders

    void OnlyMass(deviceSU3Vector* pTarget, Real f2am, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff);
    void OneLink(const deviceSU3* pGuage, BYTE byGaugeFieldId, deviceSU3Vector* pTarget, Real fCoefficient, 
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx, 
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff);
    void OneLinkForce(const deviceSU3* pGuage, BYTE byGaugeFieldId, deviceSU3* pForce, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx) const;

    #pragma endregion

protected:

    //For some strange boundary condition
    //Normally, eta_{\mu}(n+\mu)=eta_{\mu}, so set this = FALSE
    UBOOL m_bEachSiteEta;

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

#pragma region Some help functions to implement higher orders

/**
 * Same as CFieldFermionKSSU3R
 * full is a list of path directions with length = iLength
 * it will be divided into two list, where l is full[0, iSep], r is (full[iSep, iLength])^dagger
 * l, r should be allocated on device
 */
static __device__ __inline__ void _deviceSeperate(const INT* __restrict__ full, INT iSep, UINT iLength, INT* l, INT* r, BYTE& LL, BYTE& RL)
{
    LL = static_cast<BYTE>(iSep);
    RL = static_cast<BYTE>(iLength - iSep);

    for (INT i = 0; i < LL; ++i)
    {
        l[i] = -full[iSep - i - 1];
    }

    for (INT i = 0; i < RL; ++i)
    {
        r[i] = full[iSep + i];
    }
}

static __device__ __inline__ void _devicePathDagger(const INT* __restrict__ path, INT* res, UINT iLength)
{
    for (UINT i = 0; i < iLength; ++i)
    {
        res[i] = -path[iLength - i - 1];
    }
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3_H_

//=============================================================================
// END OF FILE
//=============================================================================