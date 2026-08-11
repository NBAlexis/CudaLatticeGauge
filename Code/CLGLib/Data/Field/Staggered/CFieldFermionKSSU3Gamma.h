//=============================================================================
// FILENAME : CFieldFermionKSSU3Gamma.h
// 
// DESCRIPTION:
// This is a helper to implement the condensations
// Do not use this to simulate unless you know what this is
//
// REVISION:
//  [09/10/2022 nbale]
//=============================================================================
#pragma once

#include "CFieldFermionKST.h"

#ifndef _CFIELDFERMIONKSSU3GAMMA_H_
#define _CFIELDFERMIONKSSU3GAMMA_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3Gamma)

class CLGAPI CFieldFermionKSSU3Gamma : public CFieldFermionKSSU3
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3Gamma)

public:

    CFieldFermionKSSU3Gamma();
    ~CFieldFermionKSSU3Gamma();

protected:

    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

    void DOperatorKSOnEvenOrOdd(void* pTargetBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

    void CalculateForceEvenOddS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const override;

    void ApplyGammaS(const CFieldGauge* pGauge, EGammaMatrix eGamma) override;

public:

    void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString& tab) const override;
    

    //whether gamma1,2,3,4 are applied as imaginary
    UBOOL m_bImagine;
    Real m_fCoeffGamma1;
    Real m_fCoeffGamma2;
    Real m_fCoeffGamma3;
    Real m_fCoeffGamma4;
    Real m_fCoeffGamma5;
    Real m_fCoeffGamma51;
    Real m_fCoeffGamma52;
    Real m_fCoeffGamma53;
    Real m_fCoeffGamma54;
    Real m_fCoeffSigma12;
    Real m_fCoeffSigma13;
    Real m_fCoeffSigma14;
    Real m_fCoeffSigma23;
    Real m_fCoeffSigma24;
    Real m_fCoeffSigma34;

    SCHAR* m_pDevicePathBuffer;

};


__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3GAMMA_H_

//=============================================================================
// END OF FILE
//=============================================================================