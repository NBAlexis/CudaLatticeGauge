//=============================================================================
// FILENAME : CFieldGaugeKernel.h
// 
// DESCRIPTION:
// This is the common class for all gauge fields
//
// REVISION:
//  [07/24/2024 nbale]
//=============================================================================

#ifndef _CFIELDGAUGE_KERNEL_H_
#define _CFIELDGAUGE_KERNEL_H_

__BEGIN_NAMESPACE

template<typename deviceGauge, INT matrixN>
class __DLL_EXPORT CFieldGaugeKernel
{
public:

    static void CalculateForceAndStaple(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, deviceGauge* pStaple, Real betaOverN);
    static void CalculateOnlyStaple(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pStaple);
    static DOUBLE CalculatePlaqutteEnergy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN);
    static DOUBLE CalculatePlaqutteEnergyUseClover(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN);
    static DOUBLE CalculatePlaqutteEnergyUsingStable(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN, const deviceGauge* pStaple);

    static void CalculateForceAndStaple_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, deviceGauge* pStaple, Real betaOverN);
    static void CalculateOnlyStaple_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pStaple);
    static DOUBLE CalculatePlaqutteEnergy_D(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN);
};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_KERNEL_H_

//=============================================================================
// END OF FILE
//=============================================================================