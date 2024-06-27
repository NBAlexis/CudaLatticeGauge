//=============================================================================
// FILENAME : CFieldBosonU1NoGauge.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [3/31/2024 nbale]
//=============================================================================

#ifndef _CFIELDBOSONU1_H_
#define _CFIELDBOSONU1_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldBosonU1)

class CLGAPI CFieldBosonU1 : public CFieldBoson
{
    __CLGDECLARE_FIELD(CFieldBosonU1)

public:
    CFieldBosonU1();
    ~CFieldBosonU1();

    /**
    * This should be momentum field
    */
    void MakeRandomMomentum() override;

    void D(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldBoson* const* pBoson, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    void ForceOnGauge(INT gaugeNum, INT bosonNum, const CFieldGauge* const* pGauge, const CFieldGauge** pGaugeForce, const CFieldBoson* const* pBoson) override;
    //void DD(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    //void Square(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;

    EFieldType GetFieldType() const override
    {
        return EFT_BosonU1;
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
    void FieldMultply(const CFieldBoson* x, UBOOL bConj = TRUE) override;
    cuDoubleComplex Dot(const CField* other) const override;

    BYTE* CopyDataOut(UINT& uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;

    CLGComplex* m_pDeviceData;

    _GetData

};


__END_NAMESPACE

#endif //#ifndef _CFIELDBOSONU1_H_

//=============================================================================
// END OF FILE
//=============================================================================