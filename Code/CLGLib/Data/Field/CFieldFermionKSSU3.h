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

class CLGAPI CFieldFermionKSSU3 : public CFieldFermionKS
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
    void DebugPrintRed() const;

    void Dagger() override;

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
    void ApplyGammaKS(const CFieldGauge* pGauge, EGammaMatrix eGamma) override;

    //================= test anti-hermitian =========
    UINT TestAntiHermitian(const CFieldGauge* pGauge) const override;

    //These are truely D or InverseD etc.

    /**
     * Use to calculate action, it is (D^+D)^{-1/4}
     */
    void D_MD(const CField* pGauge) override;
    void D0(const CField* pGauge) override;
    void D_MC(const CField* pGauge) override;

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

    //============================
    //Override these two functions for KS
    void DerivateD0(void* pForce, const void* pGaugeBuffer) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

    //============================
    deviceSU3Vector* m_pDeviceData;

    #pragma region Help functions to implement higher orders

    void OnlyMass(void* pTarget, Real f2am, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) override;
    void OneLink(const void* pGuage, BYTE byGaugeFieldId, void* pTarget, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx, 
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) override;
    void OneLinkForce(const void* pGuage, BYTE byGaugeFieldId, void* pForce, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx) const override;

    #pragma endregion

protected:

    static void Seperate(INT* full, INT iSep, INT* l, INT* r, BYTE& LL, BYTE& RL)
    {
        LL = static_cast<BYTE>(iSep);
        RL = static_cast<BYTE>(3 - iSep);

        for (INT i = 0; i < LL; ++i)
        {
            //trace back
            l[i] = -full[iSep - i - 1];

            //If iSep = 0, This loop will not enter
            //If iSep = 1, This is -full[0]
            //If iSep = 2, This is -full[1], -full[0]
        }

        for (INT i = 0; i < RL; ++i)
        {
            r[i] = full[iSep + i];
        }
    }

    //phi _i and Dst0 phi _i
    deviceSU3Vector** m_pRationalFieldPointers;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3_H_

//=============================================================================
// END OF FILE
//=============================================================================