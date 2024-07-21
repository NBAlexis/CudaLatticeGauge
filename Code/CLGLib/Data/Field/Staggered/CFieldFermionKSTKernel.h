//=============================================================================
// FILENAME : CFieldFermionKSTKernel.h
// 
// DESCRIPTION:
// This is the class for Kogut-Susskind staggered fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [07/21/2024 nbale]
//=============================================================================
#ifndef _CFIELDFERMIONKST_KERNEL_H_
#define _CFIELDFERMIONKST_KERNEL_H_

__BEGIN_NAMESPACE

#if 0

template<typename deviceVector, typename deviceGauge, INT vectorN>
class __DLL_EXPORT CFieldFermionKST : public CFieldFermionKS
{

public:

    CFieldFermionKST();
    ~CFieldFermionKST();
    EFieldType GetFieldType() const override { return EFT_Max; }
    CField* GetCopy() const override
    {
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* ret = new CFieldFermionKST<deviceVector, deviceGauge, vectorN>();
        CopyTo(ret);
        return ret;
    }

    void InitialField(EFieldInitialType eInitialType) override;
    void InitialFieldWithFile(const CCString&, EFieldFileType) override;
    void InitialWithByte(BYTE* byData) override;
    void InitialOtherParameters(CParameters& params) override;
    void DebugPrintMe() const override;
    void CopyTo(CField* target) const override;

    void Dagger() override;

    //This is Axpy(1.0f, x)
    void AxpyPlus(const CField* x) override;
    void AxpyMinus(const CField* x) override;
    void Axpy(Real a, const CField* x) override;
    void Axpy(const CLGComplex& a, const CField* x) override;
    void Mul(const CField* other, UBOOL bDagger = TRUE) override;
    void ScalarMultply(const CLGComplex& a) override;
    void ScalarMultply(Real a) override;
    cuDoubleComplex Dot(const CField* other) const override;

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

    void InitialAsSource(const SFermionBosonSource& sourceData) override;
    BYTE* CopyDataOut(UINT& uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;

protected:

    //============================
    //Override these two functions for KS
    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;


public:

    #pragma region Help functions to implement higher orders

    void OnlyMass(void* pTarget, Real f2am, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) override;

    #pragma endregion

    deviceVector* m_pDeviceData;

    _GetData

protected:

    void OneLinkS(const void* pGuage, BYTE byGaugeFieldId, void* pTarget, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) override;
    void OneLinkForceS(const void* pGuage, BYTE byGaugeFieldId, void* pForce, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx) const override;


    //phi _i and Dst0 phi _i
    deviceVector** m_pRationalFieldPointers;

public:

    /**
     * This is for simulation, 2a is already multiplied.
     * 2a qbar Gamma q
     * for example, gamma_i  -> 1 x chichi
     *              sigma ij -> 1/2 x chichi
     *              gamma 5i -> 1/4 x chichi
     *              gamma 5  -> 1/8 x chichi
     *
     * Note: 2a is multiplied, therefore when measuring, one should use half coefficient
     * Note: Gamma_mu, and Sigma _ ij, the "i" is already multiplied so that no sign problem when simulating, it should be "-i" if recover the sign problem
     * Note: SIGMA31 is SIGMA13
     *       SIGMA41 is SIGMA14
     *       SIGMA42 is SIGMA24
     *       SIGMA43 is SIGMA34
     *
     */
    static void appApplyGammaKS(
        void* pTargetBuffer,
        const void* pBuffer,
        const void* pGaugeBuffer,
        EGammaMatrix eGamma,
        UBOOL bShiftCenter,
        UBOOL bDagger,
        Real fGammaCoeff,
        EOperatorCoefficientType eOCT,
        Real fRealCoeff,
        CLGComplex cCmpCoeff,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    /**
     * devicePathBuffer must be larger than 4
     */
    static void GammaKSForce(
        void* pForce,
        const void* pGaugeBuffer,
        const deviceVector* const* pRationalFields,
        const Real* pRationalNumerator,
        UINT uiRationalDegree,
        Real fCoeff,
        EGammaMatrix eGamma,
        INT* devicePathBuffer,
        BYTE byFieldID,
        BYTE byGaugeFieldID);
};

template<typename deviceVector, typename deviceGauge, INT vectorN>
class __DLL_EXPORT CFieldMatrixOperationKST : public CFieldMatrixOperation
{
public:
    CFieldMatrixOperationKST();
    ~CFieldMatrixOperationKST();

    //real left = (res,left)
    void VectorMultiplyMatrix(TArray<CField*>& res, const TArray<CField*>& left, const CLGComplex* deviceMatrix, UINT uiDimX, UINT uiDimY) override;

    deviceVector** m_pResBuffer;
    deviceVector** m_pLeftBuffer;
    deviceVector** m_pHostResBuffer;
    deviceVector** m_pHostLeftBuffer;
};

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSU1)
class CLGAPI CFieldFermionKSU1 : public CFieldFermionKST<CLGComplex, CLGComplex, 1>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionKSU1)
public:
    EFieldType GetFieldType() const override { return EFT_FermionStaggeredU1; }
};

#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKST_KERNEL_H_

//=============================================================================
// END OF FILE
//=============================================================================