//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3DR.h
// 
// DESCRIPTION:
//
// Dirichlet and rotation
// 
// It also supports torus and projective plane boundary condition
// NOTE: In the case of projective plane, the Omega GA4 (x Dy - y Dx) term does not satisfy Gamma5-Hermiticity on the boundary, so we made a modification
//
// REVISION:
//  [05/19/2019 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONWILSONSQUARESU3DR_H_
#define _CFIELDFERMIONWILSONSQUARESU3DR_H_

//Not sure this is faster, need to test
#define _CLG_ROTATING_NEW_IMP 1

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareSU3DR)

class CLGAPI CFieldFermionWilsonSquareSU3DR : public CFieldFermionWilsonSquareSU3D
{
    __CLGDECLARE_FIELD(CFieldFermionWilsonSquareSU3DR)

public:

    CFieldFermionWilsonSquareSU3DR() 
        : CFieldFermionWilsonSquareSU3D()
        , m_bNaive(TRUE)
        , m_bExponential(FALSE)
        , m_bShiftCenter(FALSE)
    {
    }

protected:

    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;

public:

    void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString &tab) const override;

    UBOOL m_bNaive;
    UBOOL m_bExponential;
    UBOOL m_bShiftCenter;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3DR_H_

//=============================================================================
// END OF FILE
//=============================================================================