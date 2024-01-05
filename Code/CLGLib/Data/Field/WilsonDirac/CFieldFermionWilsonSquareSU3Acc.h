//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3Acc.h
// 
// DESCRIPTION:
//
// Dirichlet and rotation
//
// REVISION:
//  [07/27/2020 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONWILSONSQUARESU3ACC_H_
#define _CFIELDFERMIONWILSONSQUARESU3ACC_H_

//Not sure this is faster, need to test

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareSU3Acc)

class CLGAPI CFieldFermionWilsonSquareSU3Acc : public CFieldFermionWilsonSquareSU3
{
    __CLGDECLARE_FIELD(CFieldFermionWilsonSquareSU3Acc)

public:

    CFieldFermionWilsonSquareSU3Acc()
        : CFieldFermionWilsonSquareSU3()
    {
    }

    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, 
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const override;

    CCString GetInfos(const CCString &tab) const override;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3ACC_H_

//=============================================================================
// END OF FILE
//=============================================================================