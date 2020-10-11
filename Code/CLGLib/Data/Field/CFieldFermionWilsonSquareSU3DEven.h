//=============================================================================
// FILENAME : CFieldFermionWilsonSU3DEven.h
// 
// DESCRIPTION:
//   
//   Do not directly use it, copy a CFieldFermionWilsonSquareSU3 or CFieldFermionWilsonSquareSU3D instead
//   Assume Nx * Ny is even, for convinient to decompse threads
//
// REVISION:
//  [05/18/2019 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONWILSONSU3DEVEN_H_
#define _CFIELDFERMIONWILSONSU3DEVEN_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSU3DEven)

class CLGAPI CFieldFermionWilsonSU3DEven : public CFieldFermionWilsonSquareSU3D
{
    __CLGDECLARE_FIELD(CFieldFermionWilsonSU3DEven)

public:

    CFieldFermionWilsonSU3DEven() : CFieldFermionWilsonSquareSU3D() {}

    //1 - kappa^2 Doo^(-1) Deo Dee^(-1)Doe
    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

    void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const override;

    CCString GetInfos(const CCString& tab) const override;

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
    void Dagger() override;

    void WriteEvenSites(const CFieldFermion* pParentField, const CFieldGauge* pGauge, UBOOL bDdagger) override;
    void WriteBackEvenSites(CFieldFermion* pParentField, const CFieldGauge* pGauge, UBOOL bDdagger) const override;
    UBOOL IsEvenField() const override { return TRUE; }

    static void SetOddZero(deviceWilsonVectorSU3* pBuffer);
    BYTE m_byParentId;
};


__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSU3DEVEN_H_

//=============================================================================
// END OF FILE
//=============================================================================