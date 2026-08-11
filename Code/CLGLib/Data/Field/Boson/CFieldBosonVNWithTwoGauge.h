//=============================================================================
// FILENAME : CFieldBosonVNTwoGauge.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [3/31/2024 nbale]
//=============================================================================

#ifndef _CFIELDBOSONVNTWOGAUGE_H_
#define _CFIELDBOSONVNTWOGAUGE_H_

#define __DEFINE_BOSON_FIELD(FIELD_NAME, TYPE_BOSON, TYPE_GAUGE, VECTOR_N, FLOAT_N, ELEMENT_TYPE) \
__CLG_REGISTER_HELPER_HEADER(FIELD_NAME) \
class CLGAPI FIELD_NAME : public CFieldBosonVN<TYPE_BOSON, TYPE_GAUGE> \
{ \
    __CLGDECLARE_FIELDWITHOUTCOPYTO(FIELD_NAME) \
public: \
    EFieldType GetFieldType() const override { return ELEMENT_TYPE; } \
    UINT VectorN() const override { return VECTOR_N; } \
    UINT FloatN() const override { return FLOAT_N; } \
};


__BEGIN_NAMESPACE

template<typename deviceDataBoson, typename deviceDataGauge1, typename deviceDataGauge2>
class __DLL_EXPORT CFieldBosonVNTwoGauge : public CFieldBoson
{
public:
    CFieldBosonVNTwoGauge();
    ~CFieldBosonVNTwoGauge();

    void CopyTo(CField* U) const override;

    /**
    * This should be momentum field
    */
    void MakeRandomMomentum() override;

    void D(INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void ForceOnGauge(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, CFieldGauge* const* pGaugeForce, const CFieldBoson* const* pBoson) const override;

    void InitialField(EFieldInitialType eInitialType) override;
    void InitialFieldWithFile(const CCString&, EFieldFileType) override;
    void InitialWithByte(BYTE* byData) override;
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
    TArray<DOUBLE> Sum() const override;

    BYTE* CopyDataOut(UINT& uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;

    deviceDataBoson* m_pDeviceData;

    _GetData

};


__END_NAMESPACE

#endif //#ifndef _CFIELDBOSONVN_H_

//=============================================================================
// END OF FILE
//=============================================================================