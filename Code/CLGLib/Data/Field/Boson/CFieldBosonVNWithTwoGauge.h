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
    TArray<DOUBLE> Sum() const override;

    BYTE* CopyDataOut(UINT& uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;

    deviceDataBoson* m_pDeviceData;

    _GetData

};

__DEFINE_BOSON_FIELD(CFieldBosonU1, CLGComplex, CLGComplex, 1, 2, EFT_BosonComplex)
__DEFINE_BOSON_FIELD(CFieldBosonSU2, deviceSU2Vector, deviceSU2, 2, 4, EFT_BosonComplexVector2)
__DEFINE_BOSON_FIELD(CFieldBosonSU3, deviceSU3Vector, deviceSU3, 3, 6, EFT_BosonComplexVector3)
__DEFINE_BOSON_FIELD(CFieldBosonSU4, deviceSU4Vector, deviceSU4, 4, 8, EFT_BosonComplexVector4)
__DEFINE_BOSON_FIELD(CFieldBosonSU5, deviceSU5Vector, deviceSU5, 5, 10, EFT_BosonComplexVector5)
__DEFINE_BOSON_FIELD(CFieldBosonSU6, deviceSU6Vector, deviceSU6, 6, 12, EFT_BosonComplexVector6)
__DEFINE_BOSON_FIELD(CFieldBosonSU7, deviceSU7Vector, deviceSU7, 7, 14, EFT_BosonComplexVector7)
__DEFINE_BOSON_FIELD(CFieldBosonSU8, deviceSU8Vector, deviceSU8, 8, 16, EFT_BosonComplexVector8)

__END_NAMESPACE

#endif //#ifndef _CFIELDBOSONVN_H_

//=============================================================================
// END OF FILE
//=============================================================================