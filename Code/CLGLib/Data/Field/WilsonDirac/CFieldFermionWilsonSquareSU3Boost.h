//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3Boost.h
// 
// DESCRIPTION:
//
// Dirichlet and rotation
//
// REVISION:
//  [08/03/2020 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONWILSONSQUARESU3BOOST_H_
#define _CFIELDFERMIONWILSONSQUARESU3BOOST_H_

//Not sure this is faster, need to test

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareSU3Boost)

class CLGAPI CFieldFermionWilsonSquareSU3Boost : public CFieldFermionWilsonSquareSU3
{
    __CLGDECLARE_FIELD(CFieldFermionWilsonSquareSU3Boost)

public:

    CFieldFermionWilsonSquareSU3Boost()
        : CFieldFermionWilsonSquareSU3()
    {
    }

protected:

    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, 
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const override;

public:

    CCString GetInfos(const CCString &tab) const override;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3ACC_H_

//=============================================================================
// END OF FILE
//=============================================================================