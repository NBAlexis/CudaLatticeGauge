//=============================================================================
// FILENAME : CFieldBosonVNKernel.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [mm/dd/yy]
//  [07/20/2024 nbale]
//=============================================================================
#pragma once

#include "Tools/Math/DeviceTemplates/DeviceInlineUseNoTemplateFunction.h"

#ifndef _CFIELDBOSONVN_KERNEL_H_
#define _CFIELDBOSONVN_KERNEL_H_

__BEGIN_NAMESPACE

template<typename deviceDataBoson, typename deviceDataGauge>
#if _CLG_WIN
class __DLL_EXPORT CFieldBosonVNKernel
#else
class CFieldBosonVNKernel
#endif
{
public:
    static UINT CheckHermitian(const CFieldBoson* data, UINT uiSiteCount, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const CFieldTensor2* const* tensor2Fields);
    static void ForceOnGauge(const deviceDataBoson* data, BYTE byFieldId, BYTE byGaugeFieldId, const deviceDataGauge* gaugedata, deviceDataGauge* force);
    static void DFromSource(const deviceDataBoson* source, deviceDataBoson* target, BYTE byFieldId, BYTE byGaugeFieldId, const deviceDataGauge* gaugedata, EOperatorCoefficientType eCoeffType, Real fCoeffReal, const CLGComplex& cCompCoeff);

    static void OneLink(
        deviceDataBoson* pTarget,
        BYTE byFieldId,
        const deviceDataBoson* pSource,
        const deviceDataGauge* pGauge,
        BYTE byGaugeFieldId,
        DOUBLE fCoefficient,
        _deviceCoeffFunctionPointerTwoSites fpCoeff,
        const SCHAR* pDevicePath,
        BYTE pathLength,
        EOperatorCoefficientType eOCT,
        Real fRealCoeff,
        const CLGComplex& cCmpCoeff);

    static void OneLinkForceGauge(
        const deviceDataBoson* pBoson,
        BYTE byFieldId,
        const deviceDataGauge* pGauge,
        BYTE byGaugeFieldId,
        deviceDataGauge* pForce,
        DOUBLE fCoefficient,
        _deviceCoeffFunctionPointerTwoSites fpCoeff,
        const SCHAR* pDevicePath,
        BYTE pathLength);

    static void PartialSq(
        deviceDataBoson* pTarget,
        BYTE byFieldId,
        const deviceDataBoson* pSource,
        const deviceDataGauge* pGauge,
        BYTE byGaugeFieldId,
        DOUBLE fCoefficient,
        _deviceCoeffFunctionPointer fpCoeff,
        BYTE idir,
        EOperatorCoefficientType eOCT,
        Real fRealCoeff,
        const CLGComplex& cCmpCoeff);

    static void PartialSqForceGauge(
        const deviceDataBoson* pBoson,
        BYTE byFieldId,
        const deviceDataGauge* pGauge,
        BYTE byGaugeFieldId,
        deviceDataGauge* pForce,
        DOUBLE fCoefficient,
        _deviceCoeffFunctionPointer fpCoeff,
        BYTE idir);

    static void AllocatePathBuffer(SCHAR** pathbuffer);
    static void FreePathBuffer(SCHAR* pathbuffer);
    static void CopyPathBuffer(SCHAR* devicepathbuffer, const SCHAR* hostpathbuffer, BYTE length);

    #pragma region rotation

    static void CopyFunctionPointCx(_deviceCoeffFunctionPointerTwoSites* target);
    static void CopyFunctionPointCy(_deviceCoeffFunctionPointerTwoSites* target);
    static void CopyFunctionPointCxy(_deviceCoeffFunctionPointerTwoSites* target);
    static void CopyFunctionPointCxShift(_deviceCoeffFunctionPointerTwoSites* target);
    static void CopyFunctionPointCyShift(_deviceCoeffFunctionPointerTwoSites* target);
    static void CopyFunctionPointCxyShiftXYPP(_deviceCoeffFunctionPointerTwoSites* target);
    static void CopyFunctionPointCxyShiftYXPP(_deviceCoeffFunctionPointerTwoSites* target);
    static void CopyFunctionPointCxyShiftYXPM(_deviceCoeffFunctionPointerTwoSites* target);
    static void CopyFunctionPointCxyShiftXYMP(_deviceCoeffFunctionPointerTwoSites* target);

    static void CopyFunctionPointCxCySq(_deviceCoeffFunctionPointer* target);
    static void CopyFunctionPointCxCySqShift(_deviceCoeffFunctionPointer* target);
    static void CopyFunctionPointCxSq(_deviceCoeffFunctionPointer* target);
    static void CopyFunctionPointCySq(_deviceCoeffFunctionPointer* target);
    static void CopyFunctionPointCxSqShift(_deviceCoeffFunctionPointer* target);
    static void CopyFunctionPointCySqShift(_deviceCoeffFunctionPointer* target);

    #pragma endregion
};

#if !_CLG_WIN
extern template class CFieldBosonVNKernel<Real, Real>;
extern template class CFieldBosonVNKernel<CLGComplex, CLGComplex>;
extern template class CFieldBosonVNKernel<deviceSU2Vector, deviceSU2>;
extern template class CFieldBosonVNKernel<deviceSU3Vector, deviceSU3>;

#if _CLG_SU4_BOSON
extern template class CFieldBosonVNKernel<deviceSU4Vector, deviceSU4>;
#endif
#if _CLG_SU5_BOSON
extern template class CFieldBosonVNKernel<deviceSU5Vector, deviceSU5>;
#endif
#if _CLG_SU6_BOSON
extern template class CFieldBosonVNKernel<deviceSU6Vector, deviceSU6>;
#endif
#if _CLG_SU7_BOSON
extern template class CFieldBosonVNKernel<deviceSU7Vector, deviceSU7>;
#endif
#if _CLG_SU8_BOSON
extern template class CFieldBosonVNKernel<deviceSU8Vector, deviceSU8>;
#endif
#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDBOSONVN_KERNEL_H_

//=============================================================================
// END OF FILE
//=============================================================================