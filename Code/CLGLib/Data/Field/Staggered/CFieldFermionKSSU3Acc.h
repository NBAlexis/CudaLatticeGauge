//=============================================================================
// FILENAME : CFieldFermionKSSU3Acc.h
// 
// DESCRIPTION:
// 
// This is only suitable when gt << 1, this uses Galilean transform
// 
// gammai (pi + iAi) + gt gammaz (pt+iAt) + m
// 
//
// REVISION:
//  [11/21/2023 nbale]
//=============================================================================
#include "CFieldFermionKSSU3.h"

#ifndef _CFIELDFERMIONKSSU3ACC_H_
#define _CFIELDFERMIONKSSU3ACC_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3Acc)

class CLGAPI CFieldFermionKSSU3Acc : public CFieldFermionKSSU3
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3Acc)

protected:

    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

public:

    CCString GetInfos(const CCString& tab) const override;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3DACC_H_

//=============================================================================
// END OF FILE
//=============================================================================