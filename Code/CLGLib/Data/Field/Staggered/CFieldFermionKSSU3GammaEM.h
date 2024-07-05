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
    void ApplyGammaKSS(const CFieldGauge* pGauge, EGammaMatrix eGamma) override;

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
    static void appApplyGammaKSEM(
        void* pTargetBuffer,
        const void* pBuffer,
        const void* pGaugeBuffer,
        const void* pEMFieldBuffer,
        Real fCharge,
        EGammaMatrix eGamma,
        UBOOL bShiftCenter,
        UBOOL bDagger,
        Real fGammaCoeff,
        EOperatorCoefficientType eOCT,
        Real fRealCoeff,
        CLGComplex cCmpCoeff,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    /**
     * devicePathBuffer must be larger than 4
     */
    static void GammaKSForceEM(
        void* pForce,
        const void* pGaugeBuffer,
        const void* pEMFieldBuffer,
        Real fCharge,
        const deviceSU3Vector* const* pRationalFields,
        const Real* pRationalNumerator,
        UINT uiRationalDegree,
        Real fCoeff,
        EGammaMatrix eGamma,
        INT* devicePathBuffer,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    static void DOperatorEM(
        void* pTargetBuffer,
        const void* pBuffer,
        const void* pGaugeBuffer,
        const void* pEMFieldBuffer,
        Real f2am,
        Real fCharge,
        UBOOL bShiftCenter,
        UBOOL bDagger,
        EOperatorCoefficientType eOCT,
        Real fRealCoeff,
        CLGComplex cCmpCoeff,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    static void KSForceEM(
        void* pForce,
        const void* pGaugeBuffer,
        const void* pEMFieldBuffer,
        Real fCharge,
        const deviceSU3Vector* const* pRationalFields,
        const Real* pRationalNumerator,
        UINT uiRationalDegree,
        BYTE byFieldID);

};


__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3GAMMAEM_H_

//=============================================================================
// END OF FILE
//=============================================================================