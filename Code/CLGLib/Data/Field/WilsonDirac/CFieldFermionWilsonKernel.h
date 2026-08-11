//=============================================================================
// FILENAME : CFieldFermionWilsonKernel.h
// 
// DESCRIPTION:
//
// REVISION:
//  [mm/dd/yy]
//  [08/04/2025 nbale]
//=============================================================================
#pragma once

#ifndef _CFIELDFERMIONWILSON_KERNEL_H_
#define _CFIELDFERMIONWILSON_KERNEL_H_

__BEGIN_NAMESPACE

class CLGAPI CFieldFermionWilsonKernel
{
public:

#pragma region Clover

    static void DOperatorClover(deviceWilsonVectorSU3* pTargetBuffer, const deviceWilsonVectorSU3* pSource,
        const deviceSU3* pGauge, const deviceSU3* pFmunu, BYTE byFieldId, BYTE byGaugeFieldId, DOUBLE fCoef,
        UBOOL bDagger, EOperatorCoefficientType eOCT,
        Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void PrepareSigmaMunu(const deviceWilsonVectorSU3* phi, const deviceWilsonVectorSU3* phid, deviceSU3* res);

    static void CloverForce(const deviceSU3* phidphi, const deviceSU3* pGauge, deviceSU3* force, BYTE byFermionId, BYTE byGaugeId, DOUBLE fCoef);

#pragma endregion

};


__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONWILSON_KERNEL_H_

//=============================================================================
// END OF FILE
//=============================================================================