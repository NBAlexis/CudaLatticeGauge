//=============================================================================
// FILENAME : CFieldFermionKSSU3P4.h
// 
// DESCRIPTION:
// The HISQ seems too expansive for use now, so we use P4 as a first improvement
//
// REVISION:
//  [10/08/2020 nbale]
//=============================================================================
#include "CFieldFermionKSSU3.h"

#ifndef _CFIELDFERMIONKSSU3P4_H_
#define _CFIELDFERMIONKSSU3P4_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3P4)

class CLGAPI CFieldFermionKSSU3P4 : public CFieldFermionKSSU3
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3P4)

public:

    CFieldFermionKSSU3P4();
    ~CFieldFermionKSSU3P4();

protected:

    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

public:

    void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString& tab) const override;

    Real m_fomega;
    Real m_fc10;
    Real m_fc12;

    INT* m_pDevicePathBuffer;
};


__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3P4_H_

//=============================================================================
// END OF FILE
//=============================================================================