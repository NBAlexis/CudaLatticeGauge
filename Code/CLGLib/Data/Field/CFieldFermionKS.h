//=============================================================================
// FILENAME : CFieldFermionKS.h
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

#ifndef _CFIELDFERMIONKS_H_
#define _CFIELDFERMIONKS_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKS)

class CLGAPI CFieldFermionKS : public CFieldFermion
{
    __CLGDECLARE_FIELD(CFieldFermionKS)

public:

    CFieldFermionKS();
    ~CFieldFermionKS();

    EFieldType GetFieldType() const override
    {
        return EFT_FermionStaggered;
    }

    void InitialField(EFieldInitialType eInitialType) override;
    void InitialFieldWithFile(const CCString&, EFieldFileType) override;
    void InitialWithByte(BYTE* byData) override;
    void InitialOtherParameters(CParameters& params) override;
    void DebugPrintMe() const override;

    void Zero() override { InitialField(EFIT_Zero); }

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
    void D(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void Ddagger(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void DDdagger(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    UBOOL InverseD(const CField* pGauge) override;
    UBOOL InverseDdagger(const CField* pGauge) override;
    UBOOL InverseDDdagger(const CField* pGauge) override;
    void ApplyGamma(EGammaMatrix eGamma) override;
    void PrepareForHMC(const CFieldGauge* pGauge) override;

    UBOOL CalculateForce(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const override;
    void InitialAsSource(const SFermionSource& sourceData) override;
    void SaveToFile(const CCString& fileName) const override;
    BYTE* CopyDataOut(UINT& uiSize) const override;
    TArray<CFieldFermion*> GetSourcesAtSiteFromPool(const class CFieldGauge* pGauge, const SSmallInt4& site) const override;
    CCString GetInfos(const CCString& tab) const override;

    void SetKai(Real fKai);

    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const override;

    deviceSU3Vector* m_pDeviceData;

protected:

    Real m_fKai;
};

class CLGAPI CFieldMatrixOperationKS : public CFieldMatrixOperation
{
public:
    CFieldMatrixOperationKS();
    ~CFieldMatrixOperationKS();

    //real left = (res,left)
    void VectorMultiplyMatrix(TArray<CField*>& res, const TArray<CField*>& left, const CLGComplex* deviceMatrix, UINT uiDimX, UINT uiDimY) override;

    deviceSU3Vector** m_pResBuffer;
    deviceSU3Vector** m_pLeftBuffer;
    deviceSU3Vector** m_pHostResBuffer;
    deviceSU3Vector** m_pHostLeftBuffer;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3_H_

//=============================================================================
// END OF FILE
//=============================================================================