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
        , m_bNaive(TRUE)
        , m_bExponential(FALSE)
    {
    }

    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, 
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const override;

    //void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString &tab) const override;

    UBOOL m_bNaive;
    UBOOL m_bExponential;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3DRIGIDACC_H_

//=============================================================================
// END OF FILE
//=============================================================================