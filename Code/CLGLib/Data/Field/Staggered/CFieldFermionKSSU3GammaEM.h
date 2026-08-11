//=============================================================================
// FILENAME : CFieldFermionKSSU3GammaEM.h
// 
// DESCRIPTION:
// This is a helper to implement the condensations
// Do not use this to simulate unless you know what this is
//
// REVISION:
//  [09/28/2022 nbale]
//=============================================================================
#pragma once

#include "CFieldFermionKSSU3Gamma.h"

#ifndef _CFIELDFERMIONKSSU3GAMMAEM_H_
#define _CFIELDFERMIONKSSU3GAMMAEM_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3GammaEM)

class CLGAPI CFieldFermionKSSU3GammaEM : public CFieldFermionKSSU3Gamma
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3GammaEM)

public:

    CFieldFermionKSSU3GammaEM();

protected:

    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void ApplyGammaS(const CFieldGauge* pGauge, EGammaMatrix eGamma) override;

    void DOperatorKSOnEvenOrOdd(void* pTargetBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

    void CalculateForceEvenOddS(const CFieldGauge* pGauge, CFieldGauge* pForce, ESolverPhase ePhase) const override;

    //void CalculateF0Additive(const CFieldGauge* pGauge, CFieldGauge* f0) const override
    //{
    //    CFieldGauge* thisf0 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(f0->m_byFieldId));
    //    thisf0->Zero();
    //    CFieldFermionKSSU3Gamma::CalculateF0Additive(pGauge, thisf0);
    //    const CFieldGaugeU1Real* pU1 = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byEMFieldID));
    //    thisf0->ApplyPhaseR(pU1, -m_fCharge);
    //    f0->AxpyPlus(thisf0);
    //    thisf0->Return();
    //}

public:

    void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString& tab) const override;

    Real m_fCharge;
    BYTE m_byEMFieldID;

    /**
     * This is for simulation, 2a is already multiplied.
     * 2a qbar Gamma q
     * for example, gamma_i  -> 1 x chichi
     *              sigma ij -> 1/2 x chichi
     *              gamma 5i -> 1/4 x chichi
     *              gamma 5  -> 1/8 x chichi
     */
    //static void appApplyGammaKSEM(
    //    void* pTargetBuffer,
    //    const void* pBuffer,
    //    const void* pGaugeBuffer,
    //    const void* pEMFieldBuffer,
    //    Real fCharge,
    //    EGammaMatrix eGamma,
    //    UBOOL bShiftCenter,
    //    UBOOL bDagger,
    //    Real fGammaCoeff,
    //    EOperatorCoefficientType eOCT,
    //    Real fRealCoeff,
    //    CLGComplex cCmpCoeff,
    //    BYTE byFieldID,
    //    BYTE byGaugeFieldID);

    //static void appApplyGammaKSEMEvenOdd(
    //    void* pTargetBuffer,
    //    const void* pGaugeBuffer,
    //    const void* pEMFieldBuffer,
    //    Real fCharge,
    //    EGammaMatrix eGamma,
    //    UBOOL bShiftCenter,
    //    UBOOL bEven,
    //    UBOOL bDagger,
    //    Real fGammaCoeff,
    //    EOperatorCoefficientType eOCT,
    //    Real fRealCoeff,
    //    CLGComplex cCmpCoeff,
    //    BYTE byFieldID,
    //    BYTE byGaugeFieldID);

    ///**
    // * devicePathBuffer must be larger than 4
    // */
    //static void GammaKSForceEM(
    //    void* pForce,
    //    const void* pGaugeBuffer,
    //    const void* pEMFieldBuffer,
    //    Real fCharge,
    //    const deviceSU3Vector* const* pRationalFields,
    //    const Real* pRationalNumerator,
    //    UINT uiRationalDegree,
    //    Real fCoeff,
    //    EGammaMatrix eGamma,
    //    SCHAR* devicePathBuffer,
    //    BYTE byFieldID,
    //    BYTE byGaugeFieldID);

    //static void DOperatorEM(
    //    void* pTargetBuffer,
    //    const void* pBuffer,
    //    const void* pGaugeBuffer,
    //    const void* pEMFieldBuffer,
    //    Real f2am,
    //    Real fCharge,
    //    UBOOL bShiftCenter,
    //    UBOOL bDagger,
    //    EOperatorCoefficientType eOCT,
    //    Real fRealCoeff,
    //    CLGComplex cCmpCoeff,
    //    BYTE byFieldID,
    //    BYTE byGaugeFieldID);

    //static void DOperatorEMEvenOdd(
    //    void* pTargetBuffer,
    //    const void* pGaugeBuffer,
    //    const void* pEMFieldBuffer,
    //    Real f2am,
    //    Real fCharge,
    //    UBOOL bShiftCenter,
    //    UBOOL bEven,
    //    UBOOL bDagger,
    //    EOperatorCoefficientType eOCT,
    //    Real fRealCoeff,
    //    CLGComplex cCmpCoeff,
    //    BYTE byFieldID,
    //    BYTE byGaugeFieldID);

    //static void KSForceEM(
    //    void* pForce,
    //    const void* pGaugeBuffer,
    //    const void* pEMFieldBuffer,
    //    Real fCharge,
    //    const deviceSU3Vector* const* pRationalFields,
    //    const Real* pRationalNumerator,
    //    UINT uiRationalDegree,
    //    BYTE byFieldID);

};


__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3GAMMAEM_H_

//=============================================================================
// END OF FILE
//=============================================================================