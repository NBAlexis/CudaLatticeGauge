//=============================================================================
// FILENAME : CFieldBosonVNKernel.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [07/20/2024 nbale]
//=============================================================================
#include "Tools/Math/DeviceTemplates/DeviceInlineUseNoTemplateFunction.h"

#ifndef _CFIELDBOSONVN_KERNEL_H_
#define _CFIELDBOSONVN_KERNEL_H_

__BEGIN_NAMESPACE

template<typename deviceDataBoson, typename deviceDataGauge>
class __DLL_EXPORT CFieldBosonVNKernel
{
public:
    static UINT CheckHermitian(const CFieldBoson* data, UINT uiSiteCount, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson);
    static void ForceOnGauge(const deviceDataBoson* data, BYTE byFieldId, BYTE byGaugeFieldId, const deviceDataGauge* gaugedata, deviceDataGauge* force);
    static void DFromSource(const deviceDataBoson* source, deviceDataBoson* target, BYTE byFieldId, BYTE byGaugeFieldId, const deviceDataGauge* gaugedata, EOperatorCoefficientType eCoeffType, Real fCoeffReal, const CLGComplex& cCompCoeff);

    static void DiagnalTerm(
        deviceDataBoson* pTarget,
        BYTE byFieldId,
        const deviceDataBoson* pSource, Real fCoeffiecient, _deviceCoeffFunctionPointer fpCoeff,
        EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void OneLink(
        deviceDataBoson* pTarget,
        BYTE byFieldId,
        const deviceDataBoson* pSource,
        const deviceDataGauge* pGauge,
        BYTE byGaugeFieldId,
        Real fCoefficient,
        _deviceCoeffFunctionPointerTwoSites fpCoeff,
        const INT* pDevicePath,
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
        Real fCoefficient,
        _deviceCoeffFunctionPointerTwoSites fpCoeff,
        const INT* pDevicePath,
        BYTE pathLength);

    static void PartialSq(
        deviceDataBoson* pTarget,
        BYTE byFieldId,
        const deviceDataBoson* pSource,
        const deviceDataGauge* pGauge,
        BYTE byGaugeFieldId,
        Real fCoefficient,
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
        Real fCoefficient,
        _deviceCoeffFunctionPointer fpCoeff,
        BYTE idir);

    static void AllocatePathBuffer(INT** pathbuffer);
    static void FreePathBuffer(INT* pathbuffer);
    static void CopyPathBuffer(INT* devicepathbuffer, const INT* hostpathbuffer, BYTE length);

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

__END_NAMESPACE

#endif //#ifndef _CFIELDBOSONVN_KERNEL_H_

//=============================================================================
// END OF FILE
//=============================================================================