//=============================================================================
// FILENAME : CFieldFermionWilson.h
// 
// DESCRIPTION:
// This is the class for pseudofermion for Wilson fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [12/25/2018 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONWILSONSQUARESU3_H_
#define _CFIELDFERMIONWILSONSQUARESU3_H_

__BEGIN_NAMESPACE

class CLGAPI CFieldFermionWilsonSquareSU3 : public CFieldFermion
{
    __CLGDECLARE_FIELD(CFieldFermionWilsonSquareSU3)

public:

    CFieldFermionWilsonSquareSU3();
    ~CFieldFermionWilsonSquareSU3();

    virtual EFieldType GetFieldType() const 
    {
        return EFT_FermionWilsonSquareSU3;
    }

    virtual void InitialField(EFieldInitialType eInitialType);
    virtual void InitialFieldWithFile(const CCString&, EFieldFileType);
    virtual void InitialOtherParameters(CParameters& params);
    virtual void DebugPrintMe() const;

    virtual void Zero() { InitialField(EFIT_Zero); }
    virtual void Indentity() { appCrucial(_T("Not supported for CFermionWilsonSquareSU3!")); }

    //This is Axpy(1.0f, x)
    virtual void AxpyPlus(const CField* x);
    virtual void AxpyMinus(const CField* x);
    virtual void Axpy(Real a, const CField* x);
    virtual void Axpy(const _Complex& a, const CField* x);
    virtual void ScalarMultply(const _Complex& a);
    virtual void ScalarMultply(Real a);

    virtual _Complex Dot(const CField* other) const;

    //pGauge must be gauge SU3
    virtual void D(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0));
    virtual void Ddagger(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0));
    virtual void DDdagger(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0));
    virtual UBOOL InverseD(const CField* pGauge);
    virtual UBOOL InverseDdagger(const CField* pGauge);
    virtual UBOOL InverseDDdagger(const CField* pGauge);
    virtual void ApplyGamma(EGammaMatrix eGamma);
    virtual void PrepareForHMC(const CFieldGauge* pGauge);

    virtual UBOOL CalculateForce(const CFieldGauge* pGauge, CFieldGauge* pForce, CFieldGauge* pCachedForce) const;
    virtual void SaveToFile(const CCString &fileName) const;
    virtual CCString GetInfos(const CCString &tab) const;
    void SetKai(Real fKai);

protected:

    Real m_fKai;

    deviceWilsonVectorSU3* m_pDeviceData;

};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3_H_

//=============================================================================
// END OF FILE
//=============================================================================