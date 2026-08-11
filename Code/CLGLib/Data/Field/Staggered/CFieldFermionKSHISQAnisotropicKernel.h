//=============================================================================
// FILENAME : CFieldFermionKSHISQAnisotropicKernel.h
// 
// DESCRIPTION:
// The kernels for the anisotropic HISQ (aHISQ) staggered fermion.
// These are the fXiF-weighted (temporal, dir == 3) versions of the KS/HISQ
// kernels, split out of CFieldFermionKSTKernel and CFieldCommonKernel so that
// the shared (performance critical) kernels stay untouched.
//
// REVISION:
//  [mm/dd/yy]
//  [07/26/2026 nbale]
//=============================================================================
#pragma once

#ifndef _CFIELDFERMIONKS_HISQANISOTROPIC_KERNEL_H_
#define _CFIELDFERMIONKS_HISQANISOTROPIC_KERNEL_H_

__BEGIN_NAMESPACE

template<typename deviceVector, typename deviceGauge, INT vectorN>
#if _CLG_WIN
class __DLL_EXPORT CFieldFermionKSHISQAnisotropicKernel
#else
class CFieldFermionKSHISQAnisotropicKernel
#endif
{
public:

    /**
    * Same as CFieldFermionKSTKernel::DOperatorKS, but the dir == 3 (temporal) contribution is multiplied by fXiF
    */
    static void DOperatorKS(UBOOL bEachSiteEta, deviceVector* pTargetBuffer, const deviceVector* pBuffer,
        const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff, Real fXiF);

    /**
    * Same as CFieldFermionKSTKernel::DOperatorKSOnEvenOrOdd, but the dir == 3 (temporal) contribution is multiplied by fXiF
    */
    static void DOperatorKSOnEvenOrOdd(deviceVector* pTargetBuffer,
        const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff, Real fXiF);

    /**
    * Same as CFieldFermionKSTKernel::DerivateD0, but the dir == 3 (temporal) contribution is multiplied by fXiF
    */
    static void DerivateD0(
        const deviceVector* pFermion,
        BYTE byFieldId,
        deviceGauge* pForce,
        const deviceGauge* pGaugeBuffer,
        BYTE byGaugeFieldId,
        const deviceVector* const* pRationalFields,
        const Real* pNumerator,
        UINT uiRationApproxOrder,
        Real fXiF);

#pragma region HISQ

    static void NaikConnection(Real fNumerator, const deviceVector* rfield, deviceGauge* naikforce, BYTE byFieldId, Real fXiF);
    static void NaikConnection2(Real fNumerator, const deviceVector* rfield, const deviceVector* rfield_p_mu, deviceGauge* naikforce, BYTE byFieldId, Real fXiF);

    /**
    * suppose to work as DerivateD0, not fully tested yet
    */
    static void AddConnection2(Real fNumerator, const deviceVector* rfield, const deviceVector* rfield_p_mu, deviceGauge* naikforce, BYTE byFieldId, Real fXiF);

    /**
    * Same as CFieldFermionKSTKernel::DOperatorNaik, but the dir == 3 (temporal) contribution is multiplied by fXiF
    */
    static void DOperatorNaik(UBOOL bEachSiteEta, deviceVector* pTargetBuffer, const deviceVector* pBuffer,
        BYTE byFieldId, BYTE byGaugeFieldId, Real fNaik, Real fEpsilon,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff, Real fXiF);

    /**
    * Same as CFieldFermionKSTKernel::DOperatorNaikOnEvenOrOdd, but the dir == 3 (temporal) contribution is multiplied by fXiF
    */
    static void DOperatorNaikOnEvenOrOdd(deviceVector* pTargetBuffer,
        BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEven, Real fNaik, Real fEpsilon,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff, Real fXiF);

#pragma endregion

    /**
    * Same as CCommonKernelMV::ConnectionOneFieldStaggered, but the dir == 3 (temporal) contribution is multiplied by fXiF
    */
    static void ConnectionOneFieldStaggered(const deviceVector* v, deviceGauge* res, BYTE byFieldId, Real fCoeff, Real fXiF);

    /**
    * Same as CCommonKernelMV::AddConnectionOneFieldStaggered, but the dir == 3 (temporal) contribution is multiplied by fXiF
    */
    static void AddConnectionOneFieldStaggered(const deviceVector* v, deviceGauge* res, BYTE byFieldId, Real fCoeff, Real fXiF);

};

#if !_CLG_WIN
extern template class CFieldFermionKSHISQAnisotropicKernel<CLGComplex, CLGComplex, 1>;
extern template class CFieldFermionKSHISQAnisotropicKernel<deviceSU2Vector, deviceSU2, 2>;
extern template class CFieldFermionKSHISQAnisotropicKernel<deviceSU3Vector, deviceSU3, 3>;
extern template class CFieldFermionKSHISQAnisotropicKernel<deviceSU4Vector, deviceSU4, 4>;
#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKS_HISQANISOTROPIC_KERNEL_H_

//=============================================================================
// END OF FILE
//=============================================================================
