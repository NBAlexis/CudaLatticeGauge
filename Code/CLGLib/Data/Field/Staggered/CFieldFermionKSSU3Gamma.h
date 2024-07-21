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
#include "CFieldFermionKSSU3.h"

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

    void ApplyGammaKSS(const CFieldGauge* pGauge, EGammaMatrix eGamma) override;

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

    INT* m_pDevicePathBuffer;

    /**
     * This is for simulation, 2a is already multiplied.
     * 2a qbar Gamma q
     * for example, gamma_i  -> 1 x chichi
     *              sigma ij -> 1/2 x chichi
     *              gamma 5i -> 1/4 x chichi
     *              gamma 5  -> 1/8 x chichi
     * 
     * Note: 2a is multiplied, therefore when measuring, one should use half coefficient
     * Note: Gamma_mu, and Sigma _ ij, the "i" is already multiplied so that no sign problem when simulating, it should be "-i" if recover the sign problem
     * Note: SIGMA31 is SIGMA13
     *       SIGMA41 is SIGMA14
     *       SIGMA42 is SIGMA24
     *       SIGMA43 is SIGMA34
     * 
     */
    static void appApplyGammaKS(
        void* pTargetBuffer,
        const void* pBuffer,
        const void* pGaugeBuffer,
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
    static void GammaKSForce(
        void* pForce,
        const void* pGaugeBuffer,
        const deviceSU3Vector* const* pRationalFields,
        const Real* pRationalNumerator,
        UINT uiRationalDegree,
        Real fCoeff,
        EGammaMatrix eGamma,
        INT* devicePathBuffer,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    /**
     * every time set gamma coefficient, update the parameters of pooled
     */
    //void UpdatePooledParamters() const;

};


__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3GAMMA_H_

//=============================================================================
// END OF FILE
//=============================================================================