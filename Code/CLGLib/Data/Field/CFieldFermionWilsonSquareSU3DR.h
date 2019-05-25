//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3DR.h
// 
// DESCRIPTION:
//
// Dirichlet and rotation
//
// REVISION:
//  [05/19/2019 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONWILSONSQUARESU3DR_H_
#define _CFIELDFERMIONWILSONSQUARESU3DR_H_

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
    {
    }

    virtual void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, 
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const;
    virtual void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const;

    virtual void InitialOtherParameters(CParameters& params);
    virtual CCString GetInfos(const CCString &tab) const;

    UBOOL m_bNaive;
    UBOOL m_bExponential;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3DR_H_

//=============================================================================
// END OF FILE
//=============================================================================