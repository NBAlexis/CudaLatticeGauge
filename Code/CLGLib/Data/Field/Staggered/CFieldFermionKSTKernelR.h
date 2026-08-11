//=============================================================================
// FILENAME : CFieldFermionKSTKernelR.h
// 
// DESCRIPTION:
// Split the kernels to improve build speed
//
// REVISION:
//  [mm/dd/yy]
//  [05/12/2025 nbale]
//=============================================================================
#pragma once

#ifndef _CFIELDFERMIONKST_KERNELR_H_
#define _CFIELDFERMIONKST_KERNELR_H_

__BEGIN_NAMESPACE

template<typename deviceVector, typename deviceGauge, INT vectorN>
#if _CLG_WIN
class __DLL_EXPORT CFieldFermionKSTKernelR
#else
class CFieldFermionKSTKernelR
#endif
{
public:

    static void DOperatorKS_R_RealRotation(
        DOUBLE fOmega,
        UBOOL bShiftHalfCoord, deviceVector* pTargetBuffer, const deviceVector* pBuffer,
        const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void DOperatorKS_R_ImaginaryRotation(
        DOUBLE fOmega,
        UBOOL bShiftHalfCoord, deviceVector* pTargetBuffer, const deviceVector* pBuffer,
        const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void DOperatorKS_R_ImaginaryRotation_EM(
        DOUBLE fOmega, Real fCharge,
        UBOOL bShiftHalfCoord, deviceVector* pTargetBuffer, const deviceVector* pBuffer,
        const deviceGauge* pGaugeBuffer, const Real* pPhase, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void DOperatorKS_R_ImaginaryRotation_Cached(
        DOUBLE fOmega,
        UBOOL bShiftHalfCoord, deviceVector* pTargetBuffer, const deviceVector* pBuffer,
        const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

    //Even-odd version of DOperatorKS_R_ImaginaryRotation, in-place on one parity
    //bEven is passed to intokernalEOHalf directly, same convention as CFieldFermionKSTKernel::DOperatorKSOnEvenOrOdd
    static void DOperatorKSOnEvenOrOdd_R_ImaginaryRotation(
        DOUBLE fOmega,
        UBOOL bShiftHalfCoord, deviceVector* pData,
        const deviceGauge* pGauge, BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEven,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

    //Even-odd version of DOperatorKS_R_ImaginaryRotation with EM (U(1) phase + charge)
    static void DOperatorKSOnEvenOdd_R_ImaginaryRotation_EM(
        DOUBLE fOmega,
        UBOOL bShiftHalfCoord, deviceVector* pData,
        const deviceGauge* pGauge, const Real* pPhase, Real fCharge,
        BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEven,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

    //Even-odd version of DOperatorKS_R_ImaginaryRotation_Cached
    static void DOperatorKSOnEvenOrOdd_R_ImaginaryRotation_Cached(
        DOUBLE fOmega,
        UBOOL bShiftHalfCoord, deviceVector* pData,
        const deviceGauge* pCachedGauge, BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEven,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void DerivateD0_R(
        DOUBLE fOmega,
        UBOOL bShiftCenter,
        const deviceVector* pFermion,
        BYTE byFieldId,
        deviceGauge* pForce,
        const deviceGauge* pGaugeBuffer,
        BYTE byGaugeFieldId,
        const deviceVector* const* pRationalFields,
        const Real* pNumerator,
        UINT uiRationApproxOrder);

    //Even-odd version of DerivateD0_R
    //phi_i is non-zero only on even sites and phi_id only on odd sites, so for
    //any (n1, n2) pair exactly one of the two contracts is zero. The _EO kernels
    //only compute the surviving contract, saving half of the FLOPS.
    static void DerivateD0_ROnEvenOdd(
        DOUBLE fOmega,
        UBOOL bShiftCenter,
        const deviceVector* pFermion,
        BYTE byFieldId,
        deviceGauge* pForce,
        const deviceGauge* pGaugeBuffer,
        BYTE byGaugeFieldId,
        const deviceVector* const* pRationalFields,
        const Real* pNumerator,
        UINT uiRationApproxOrder);

    //Even-odd split version of DerivateD0_REM
    static void DerivateD0_REMOnEvenOdd(
        DOUBLE fOmega,
        UBOOL bShiftCenter,
        Real fCharge,
        const deviceVector* pFermion,
        BYTE byFieldId,
        deviceGauge* pForce,
        const deviceGauge* pGaugeBuffer,
        const Real* pPhase,
        BYTE byGaugeFieldId,
        const deviceVector* const* pRationalFields,
        const Real* pNumerator,
        UINT uiRationApproxOrder);

    static void DerivateD0_REM(
        DOUBLE fOmega,
        UBOOL bShiftCenter,
        Real fCharge,
        const deviceVector* pFermion,
        BYTE byFieldId,
        deviceGauge* pForce,
        const deviceGauge* pGaugeBuffer,
        const Real* pPhaseBuffer,
        BYTE byGaugeFieldId,
        const deviceVector* const* pRationalFields,
        const Real* pNumerator,
        UINT uiRationApproxOrder);

};

#if !_CLG_WIN
extern template class CFieldFermionKSTKernelR<CLGComplex, CLGComplex, 1>;
extern template class CFieldFermionKSTKernelR<deviceSU2Vector, deviceSU2, 2>;
extern template class CFieldFermionKSTKernelR<deviceSU3Vector, deviceSU3, 3>;
extern template class CFieldFermionKSTKernelR<deviceSU4Vector, deviceSU4, 4>;
#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKST_KERNELR_H_

//=============================================================================
// END OF FILE
//=============================================================================
