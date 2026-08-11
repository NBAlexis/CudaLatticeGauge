//=============================================================================
// FILENAME : CFieldFermionKSTKernel.h
// 
// DESCRIPTION:
// This is the class for Kogut-Susskind staggered fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [mm/dd/yy]
//  [07/21/2024 nbale]
//=============================================================================
#pragma once

#ifndef _CFIELDFERMIONKST_KERNEL_H_
#define _CFIELDFERMIONKST_KERNEL_H_

__BEGIN_NAMESPACE

template<typename deviceVector, typename deviceGauge, INT vectorN>
#if _CLG_WIN
class __DLL_EXPORT CFieldFermionKSTKernel
#else
class CFieldFermionKSTKernel
#endif
{
public:

    static UINT TestAntiHermitianS(BYTE byFieldId, const CFieldGauge* pGauge);

    static void DOperatorKS(UBOOL bEachSiteEta, deviceVector* pTargetBuffer, const deviceVector* pBuffer,
        const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void DOperatorKSOnEvenOrOdd(deviceVector* pTargetBuffer,
        const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void DerivateD0(
        const deviceVector* pFermion,
        BYTE byFieldId,
        deviceGauge* pForce,
        const deviceGauge* pGaugeBuffer,
        BYTE byGaugeFieldId,
        const deviceVector* const* pRationalFields,
        const Real* pNumerator,
        UINT uiRationApproxOrder);

    static void DOperatorKS_D(UBOOL bEachSiteEta, deviceVector* pTargetBuffer, const deviceVector* pBuffer,
        const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void DerivateD0_D(
        const deviceVector* pFermion,
        BYTE byFieldId,
        deviceGauge* pForce,
        const deviceGauge* pGaugeBuffer,
        BYTE byGaugeFieldId,
        const deviceVector* const* pRationalFields,
        const Real* pNumerator,
        UINT uiRationApproxOrder);

    static void DOperatorKSOnEvenOrOdd_D(deviceVector* pTargetBuffer,
        const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEven, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void OnlyMass(const deviceVector* pSource, deviceVector* pTarget, Real f2am, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void OneLinkS(const deviceVector* pSource, BYTE byFieldId, const deviceGauge* pGauge, BYTE byGaugeFieldId, deviceVector* pTarget, Real fCoefficient,
        const SCHAR* pDevicePath, BYTE pathLength, BYTE byEtaIdx,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void OneLinkForceS(const deviceVector* pFermion, BYTE byFieldId, const deviceGauge* pGauge, BYTE byGaugeFieldId, deviceGauge* pForce, Real fCoefficient,
        const SCHAR* pDevicePath, BYTE pathLength, BYTE byEtaIdx, const deviceVector* const* pRationalFields, const Real* pNumerator, UINT uiRationApproxOrder);

    /**
    * this was never tested
    */
    static void VectorMultiplyMatrix(deviceVector** hostResBuffer, deviceVector** hostLeftBuffer, deviceVector** resBuffer, deviceVector** leftBuffer, TArray<CField*>& res, const TArray<CField*>& left, const CLGComplex* deviceMatrix, UINT uiDimX, UINT uiDimY);


#pragma region Rotation


#pragma endregion

#pragma region EM

    static void DOperatorEM(
        deviceVector* pTargetBuffer,
        const deviceVector* pBuffer,
        const deviceGauge* pGaugeBuffer,
        const Real* pEMFieldBuffer,
        Real f2am,
        Real fCharge,
        UBOOL bShiftCenter,
        UBOOL bDagger,
        EOperatorCoefficientType eOCT,
        Real fRealCoeff,
        CLGComplex cCmpCoeff,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    /**
    * For even odd operator
    * the force is the usual force
    * than apply the phase to the force
    */
    static void DOperatorEMEvenOdd(
        deviceVector* pTargetBuffer,
        const deviceGauge* pGaugeBuffer,
        const Real* pEMFieldBuffer,
        Real f2am,
        Real fCharge,
        UBOOL bShiftCenter,
        UBOOL bEven,
        UBOOL bDagger,
        EOperatorCoefficientType eOCT,
        Real fRealCoeff,
        CLGComplex cCmpCoeff,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    static void KSForceEM(
        deviceGauge* pForce,
        const deviceGauge* pGaugeBuffer,
        const Real* pEMFieldBuffer,
        Real fCharge,
        const deviceVector* const* pRationalFields,
        const Real* pRationalNumerator,
        UINT uiRationalDegree,
        BYTE byFieldID);

#pragma endregion

#pragma region HISQ

    static void NaikConnection(Real fNumerator, const deviceVector* rfield, deviceGauge* naikforce, BYTE byFieldId);
    static void NaikConnection2(Real fNumerator, const deviceVector* rfield, const deviceVector* rfield_p_mu, deviceGauge* naikforce, BYTE byFieldId);

    /**
    * suppose to work as DerivateD0, not fully tested yet
    */
    static void AddConnection2(Real fNumerator, const deviceVector* rfield, const deviceVector* rfield_p_mu, deviceGauge* naikforce, BYTE byFieldId);

    static void DOperatorNaik(UBOOL bEachSiteEta, deviceVector* pTargetBuffer, const deviceVector* pBuffer,
        BYTE byFieldId, BYTE byGaugeFieldId, Real fNaik, Real fEpsilon,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void DOperatorNaikOnEvenOrOdd(deviceVector* pTargetBuffer,
        BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEven, Real fNaik, Real fEpsilon,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

#pragma endregion

};

#if !_CLG_WIN
extern template class CFieldFermionKSTKernel<CLGComplex, CLGComplex, 1>;
extern template class CFieldFermionKSTKernel<deviceSU2Vector, deviceSU2, 2>;
extern template class CFieldFermionKSTKernel<deviceSU3Vector, deviceSU3, 3>;

extern template class CFieldFermionKSTKernel<deviceSU4Vector, deviceSU4, 4>;
#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKST_KERNEL_H_

//=============================================================================
// END OF FILE
//=============================================================================