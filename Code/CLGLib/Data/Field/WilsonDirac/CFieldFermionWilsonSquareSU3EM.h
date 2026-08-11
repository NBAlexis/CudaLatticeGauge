//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3Gamma.h
// 
// DESCRIPTION:
// If m_bExpGamma is turned on, it is applied as
// gamma4.Exp(gamma)
//
// REVISION:
//  [mm/dd/yy]
//  [05/01/2023 nbale]
//=============================================================================
#include "CFieldFermionWilsonSquareSU3.h"

#ifndef _CFIELDFERMIONWILSONSQUARESU3EM_H_
#define _CFIELDFERMIONWILSONSQUARESU3EM_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionWilsonSquareSU3EM)

class CLGAPI CFieldFermionWilsonSquareSU3EM : public CFieldFermionWilsonSquareSU3
{
    __CLGDECLARE_FIELD(CFieldFermionWilsonSquareSU3EM)

public:

    CFieldFermionWilsonSquareSU3EM();
    ~CFieldFermionWilsonSquareSU3EM();

    void InitialOtherParameters(CParameters& params) override;

protected:

    void DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;


    void DerivateDOperator(DOUBLE fCoeff, void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;
    
    //void DDdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType = EOCT_None, Real fCoeffReal = F(1.0), Real fCoeffImg = F(0.0)) override;
    //void PrepareForHMCS(const CFieldGauge* pGauge) override;

public:

    CCString GetInfos(const CCString& tab) const override;
    Real m_fCharge;
    BYTE m_byU1FieldId;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSONSQUARESU3EM_H_

//=============================================================================
// END OF FILE
//=============================================================================