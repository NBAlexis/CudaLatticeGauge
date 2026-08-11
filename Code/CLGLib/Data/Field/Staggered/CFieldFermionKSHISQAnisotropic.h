//=============================================================================
// FILENAME : CFieldFermionKSHISQAnisotropic.h
//
// DESCRIPTION:
// The anisotropic HISQ (aHISQ) staggered fermion field.
// The bare fermion anisotropy XiF multiplies only the temporal (dir == 3)
// hopping contributions (one-link, Naik, epsilon, and the corresponding
// RHMC connections). The smearing construction is unchanged.
// All anisotropy awareness lives in this subclass, the parent classes
// (CFieldFermionHISQT and below) are untouched.
//
// REVISION:
//  [mm/dd/yy]
//  [07/25/2026 nbale]
//=============================================================================
#pragma once

#ifndef _CFIELDFERMIONKSHISQANISOTROPIC_H_
#define _CFIELDFERMIONKSHISQANISOTROPIC_H_

#include "CFieldFermionKSHISQ.h"

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionHISQSU3Anisotropic)

class CLGAPI CFieldFermionHISQSU3Anisotropic : public CFieldFermionHISQSU3
{
    __CLGDECLARE_FIELDWITHOUTCOPYTO(CFieldFermionHISQSU3Anisotropic)
public:
    CFieldFermionHISQSU3Anisotropic()
        : CFieldFermionHISQSU3()
        , m_fXiF(F(1.0))
    {

    }

    void InitialOtherParameters(CParameters& params) override;
    void CopyParamTo(CField* f) const override;
    CCString GetInfos(const CCString& tab) const override;

    //the bare fermion anisotropy
    Real m_fXiF;

protected:

    //the operator paths of the parent, with the temporal contribution weighted by m_fXiF
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void DOperatorKSOnEvenOrOdd(void* pTargetBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;

public:

    //the one-field connections, a pooled shifted field uses its own m_fXiF through these overrides
    void ConnectionSelf(void* res, Real fCoeff) const override;
    void AddConnectionSelf(void* res, Real fCoeff) const override;

    //same as the parent, except the temporal weights in NaikConnection/NaikConnection2/AddConnection2
    void CalculateF0AndNaik(const CFieldGauge* pGauge, CFieldGauge* f0, CFieldGauge* pepsilonTerm, CFieldGauge* naik) const override;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSHISQANISOTROPIC_H_

//=============================================================================
// END OF FILE
//=============================================================================
