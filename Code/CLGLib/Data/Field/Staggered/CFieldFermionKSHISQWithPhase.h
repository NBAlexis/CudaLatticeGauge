//=============================================================================
// FILENAME : CFieldFermionKSHISQWithPhase.h
// 
// DESCRIPTION:
//  Only support one gauge
//  Only support even-odd in evaluation
//  
// level2:
// one-link: (1+epsi)/8
// Lepage: -1/8
// Naik:     -(1+epsi)/24
// see: 0710.0737
// default use epsi=0
// for epsi, see:
// hep-lat/0610092
// 
// REVISION:
//  [12/30/2022 nbale]
//=============================================================================
#pragma once

#include "CFieldFermionKSHISQ.h"

#ifndef _CFIELDFERMIONKSHISQWITHPHASE_H_
#define _CFIELDFERMIONKSHISQWITHPHASE_H_

__BEGIN_NAMESPACE


__CLG_REGISTER_HELPER_HEADER(CFieldFermionHISQWithPhaseSU3)
class CLGAPI CFieldFermionHISQWithPhaseSU3 : public CFieldFermionHISQSU3
{
    __CLGDECLARE_FIELD(CFieldFermionHISQWithPhaseSU3)
public:
    CFieldFermionHISQWithPhaseSU3()
        : CFieldFermionHISQSU3()
        , m_fCharge(F(0.0))
        , m_byU1FieldId(0)
    {

    }

protected:

    /**
    * This is only used for inverse, for example the calculation of meson correlator
    */
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

    void DOperatorKSOnEvenOrOdd(void* pTargetBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

    void CalculateForceEvenOddS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const override;

public:

    void InitialOtherParameters(CParameters& params) override;
    void CalculateF0AndNaik(const CFieldGauge* pGauge, CFieldGauge* f0, CFieldGauge* pepsilonTerm, CFieldGauge* naik) const override;

    CCString GetInfos(const CCString& tab) const override;

    Real m_fCharge;
    BYTE m_byU1FieldId;
};


__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSHISQWITHPHASE_H_

//=============================================================================
// END OF FILE
//=============================================================================