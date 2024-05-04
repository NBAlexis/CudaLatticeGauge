//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3DRigidAcc.h
// 
// DESCRIPTION:
//
// Dirichlet and Rigid acceleration
//
//  It doesn't work!!!
//  It has a sign problem similar as the finite chemical potential!
//
// REVISION:
//  [08/01/2020 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONWILSONSQUARESU3DRIGIDACC_H_
#define _CFIELDFERMIONWILSONSQUARESU3DRIGIDACC_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareSU3DRigidAcc)

class CLGAPI CFieldFermionWilsonSquareSU3DRigidAcc : public CFieldFermionWilsonSquareSU3D
{
    __CLGDECLARE_FIELD(CFieldFermionWilsonSquareSU3DRigidAcc)

public:

    CFieldFermionWilsonSquareSU3DRigidAcc()
        : CFieldFermionWilsonSquareSU3D()
        , m_fG2(F(0.1))
    {
    }

protected:

    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, 
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const override;

public:

    void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString &tab) const override;
    void SetG2(Real fG2) { m_fG2 = fG2; }
    Real m_fG2;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3DRIGIDACC_H_

//=============================================================================
// END OF FILE
//=============================================================================