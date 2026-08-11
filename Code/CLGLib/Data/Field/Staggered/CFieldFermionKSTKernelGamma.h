//=============================================================================
// FILENAME : CFieldFermionKSTKernelGamma.h
// 
// DESCRIPTION:
// Split the kernels to improve build speed
//
// REVISION:
//  [mm/dd/yy]
//  [05/12/2025 nbale]
//=============================================================================
#pragma once

#ifndef _CFIELDFERMIONKST_KERNELGAMMA_H_
#define _CFIELDFERMIONKST_KERNELGAMMA_H_

__BEGIN_NAMESPACE

template<typename deviceVector, typename deviceGauge, INT vectorN>
#if _CLG_WIN
class __DLL_EXPORT CFieldFermionKSTKernelGamma
#else
class CFieldFermionKSTKernelGamma
#endif
{
public:

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

    static void appApplyGammaKSEvenOdd(
        void* pTargetBuffer,
        const void* pGaugeBuffer,
        EGammaMatrix eGamma,
        UBOOL bShiftCenter,
        UBOOL bEven,
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
        const deviceVector* const* pRationalFields,
        const Real* pRationalNumerator,
        UINT uiRationalDegree,
        Real fCoeff,
        EGammaMatrix eGamma,
        SCHAR* devicePathBuffer,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    static void GammaKSForceEvenOdd(
        void* pForce,
        Real fCoeff,
        EGammaMatrix eGamma,
        BYTE byFieldID);

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

    static void appApplyGammaKSEMEvenOdd(
        void* pTargetBuffer,
        const void* pGaugeBuffer,
        const void* pEMFieldBuffer,
        Real fCharge,
        EGammaMatrix eGamma,
        UBOOL bShiftCenter,
        UBOOL bEven,
        UBOOL bDagger,
        Real fGammaCoeff,
        EOperatorCoefficientType eOCT,
        Real fRealCoeff,
        CLGComplex cCmpCoeff,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    /**
     * devicePathBuffer must be larger than 4
     * except for the phase, so apply the phase for yourself
     */
    static void GammaKSForceEM(
        void* pForce,
        const void* pGaugeBuffer,
        const void* pEMFieldBuffer,
        Real fCharge,
        const deviceVector* const* pRationalFields,
        const Real* pRationalNumerator,
        UINT uiRationalDegree,
        Real fCoeff,
        EGammaMatrix eGamma,
        SCHAR* devicePathBuffer,
        BYTE byFieldID,
        BYTE byGaugeFieldID);


};

#if !_CLG_WIN
extern template class CFieldFermionKSTKernelGamma<CLGComplex, CLGComplex, 1>;
extern template class CFieldFermionKSTKernelGamma<deviceSU2Vector, deviceSU2, 2>;
extern template class CFieldFermionKSTKernelGamma<deviceSU3Vector, deviceSU3, 3>;

extern template class CFieldFermionKSTKernelGamma<deviceSU4Vector, deviceSU4, 4>;
#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKST_KERNELGAMMA_H_

//=============================================================================
// END OF FILE
//=============================================================================