//=============================================================================
// FILENAME : CFieldBosonVN.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [3/31/2024 nbale]
//=============================================================================

#ifndef _CFIELDBOSONVN_H_
#define _CFIELDBOSONVN_H_

__BEGIN_NAMESPACE

template<typename deviceDataBoson, typename deviceDataGauge>
class __DLL_EXPORT CFieldBosonVN : public CFieldBoson
{
public:
    CFieldBosonVN();
    ~CFieldBosonVN();

    void CopyTo(CField* U) const override;

    /**
    * This should be momentum field
    */
    void MakeRandomMomentum() override;

    void D(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void ForceOnGauge(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, CFieldGauge* const* pGaugeForce, const CFieldBoson* const* pBoson) const override;

    void InitialField(EFieldInitialType eInitialType) override;
    void InitialFieldWithFile(const CCString&, EFieldFileType) override;
    void InitialWithByte(BYTE* byData) override;
    //void InitialOtherParameters(CParameters& params) override;
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
    cuDoubleComplex Dot(const CField* other) const override;

    BYTE* CopyDataOut(UINT& uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;

    virtual UINT VectorN() const = 0;
    virtual UINT FloatN() const = 0;

    deviceDataBoson* m_pDeviceData;

    _GetData

};

__CLG_REGISTER_HELPER_HEADER(CFieldBosonU1)
class CLGAPI CFieldBosonU1 : public CFieldBosonVN<CLGComplex, CLGComplex>
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldBosonU1)
public:
    EFieldType GetFieldType() const override { return EFT_BosonU1; }
    UINT VectorN() const { return 1; }
    UINT FloatN() const { return 2; }
};

__END_NAMESPACE

#endif //#ifndef _CFIELDBOSONVN_H_

//=============================================================================
// END OF FILE
//=============================================================================