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

    virtual void ExpMult(const _Complex& , CField* ) const { appCrucial(_T("Not supported for CFermionWilsonSquareSU3!")); }

    //pGauge must be gauge SU3
    virtual void D(const CField* pGauge);
    virtual void Ddagger(const CField* pGauge);
    virtual void DDdagger(const CField* pGauge);
    virtual UBOOL InverseD(const CField* pGauge);
    virtual UBOOL InverseDdagger(const CField* pGauge);
    virtual UBOOL InverseDDdagger(const CField* pGauge);

    virtual void PrepareForHMC(const CFieldGauge* pGauge);

    virtual void CalculateForce(const CFieldGauge* pGauge, CFieldGauge* pForce);

protected:

    Real m_fKai;


    deviceWilsonVectorSU3* m_pDeviceData;
    //When calculating D phi etc, the neigbour will change, so we need to copy it before calculate
    deviceWilsonVectorSU3* m_pDeviceDataCopy;

    //For SU3, this is 8 x sizeof(m_pDeviceData), where 8 is the number of generators of SU3
    deviceWilsonVectorSU3* m_pForceRightVector;
    deviceWilsonVectorSU3* m_pForceRightVectorCopy;
    deviceWilsonVectorSU3* m_pForceLeftVector;
    deviceWilsonVectorSU3* m_pForceLeftVectorCopy;

};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3_H_

//=============================================================================
// END OF FILE
//=============================================================================