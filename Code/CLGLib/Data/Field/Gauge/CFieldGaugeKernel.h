//=============================================================================
// FILENAME : CFieldGaugeKernel.h
// 
// DESCRIPTION:
// This is the common class for all gauge fields
//
// REVISION:
//  [mm/dd/yy]
//  [07/24/2024 nbale]
//=============================================================================
#pragma once

#ifndef _CFIELDGAUGE_KERNEL_H_
#define _CFIELDGAUGE_KERNEL_H_

__BEGIN_NAMESPACE

template<typename deviceGauge, INT matrixN>
#if _CLG_WIN
class __DLL_EXPORT CFieldGaugeKernel
#else
class CFieldGaugeKernel
#endif
{
public:

    static void CalculateForceAndStaple(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, deviceGauge* pStaple, Real betaOverN);
    static void CalculateForceAndStaple(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, deviceGauge* pStaple, CLGComplex betaOverN);
    static void CalculateForceAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverN, DOUBLE xi);
    static void CalculateOnlyStaple(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pStaple);
    static void CalculateAllStaples(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* const * ppStaples);
    static void CalculateAllPlaquttes(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* const* ppPlaquettes);

    /**
    * Note that, this actually calculate:
    * (Pmunu - Pmunu^+)
    * But, Fmunu = (Pmunu - Pmunu^+) / 8i
    */
    static void CalculateFmunu(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pFmunu);

    static void CacheKSRotationGaugeBuffer(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pRes, const Real* pPhase, Real fCharge, UBOOL bHasCharge);

    static DOUBLE CalculatePlaqutteEnergy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN);
    static DOUBLE CalculatePlaqutteEnergyAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN, DOUBLE xi);
    static DOUBLE CalculatePlaqutteEnergyUseClover(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN);
    static DOUBLE CalculatePlaqutteEnergyUseCloverAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN, DOUBLE xi);
    static DOUBLE CalculatePlaqutteEnergyUsingStaple(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN, const deviceGauge* pStaple);

    static void CalculateForceAndStaple_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, deviceGauge* pStaple, Real betaOverN);
    static void CalculateForceAnisotropy_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverN, DOUBLE xi);
    static void CalculateForceAndStapleClover_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, deviceGauge* pStaple, Real betaOverN);
    static void CalculateForceCloverAnisotropy_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverN, DOUBLE xi);
    static void CalculateOnlyStaple_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pStaple);
    static DOUBLE CalculatePlaqutteEnergy_D(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN);
    static DOUBLE CalculatePlaqutteEnergyAnisotropy_D(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverN, DOUBLE xi);

    /**
    * Dirichlet boundary already considered
    */
    static DOUBLE CalculateRectangularEnergy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverNtimeRect);

    /**
    * Dirichlet boundary already considered
    */
    static DOUBLE CalculateRectangularEnergyUseClover(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverNtimeRect);

    /**
    * anisotropic rectangular energy, every loop is weighted by the plane weight:
    * spatial planes are multiplied by xi, planes containing the temporal direction by 1/xi
    */
    static DOUBLE CalculateRectangularEnergyAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverNtimeRect, DOUBLE xi);

    /**
    * Dirichlet boundary already considered
    */
    static DOUBLE CalculateRectangularEnergyUseCloverAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverNtimeRect, DOUBLE xi);

    static DOUBLE CalculateTwistedLoopEnergy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverNtimeCT);
    static void CalculateForceTwistedLoop(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeCT);

    /**
    * anisotropic twisted loop, every loop spans a triple of directions,
    * triples containing the temporal direction are multiplied by 1/xi, spatial triples by xi
    */
    static DOUBLE CalculateTwistedLoopEnergyAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, DOUBLE betaOverNtimeCT, DOUBLE xi);
    static void CalculateForceTwistedLoopAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeCT, DOUBLE xi);

    static void CalculateForceRectangular(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeRect);
    static void CalculateForceRectangular_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeRect, SSmallInt4 sPeriod);
    static void CalculateForceRectangularClover_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeRect, SSmallInt4 sPeriod);

    /**
    * anisotropic rectangular force, every staple contribution is weighted by the plane weight:
    * spatial planes are multiplied by xi, planes containing the temporal direction by 1/xi
    */
    static void CalculateForceRectangularAnisotropy(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeRect, DOUBLE xi);
    static void CalculateForceRectangularAnisotropy_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeRect, DOUBLE xi, SSmallInt4 sPeriod);
    static void CalculateForceRectangularCloverAnisotropy_D(const deviceGauge* deviceData, BYTE byFieldId, deviceGauge* pForce, DOUBLE betaOverNtimeRect, DOUBLE xi, SSmallInt4 sPeriod);

    static void NaikForce(const CFieldGauge* naikf0, CFieldGauge* naikforce);
    
};

#define _DEF_GAUGE_KERNEL_ERROR_STUBS(deviceGauge, matrixN, groupName) \
static void CalculateForceAndStaple(const deviceGauge*, BYTE, deviceGauge*, deviceGauge*, Real) \
{ \
    appCrucial(_T("CFieldGaugeKernel: CalculateForceAndStaple is not supported for " groupName " gauge group\n")); \
} \
static void CalculateForceAnisotropy(const deviceGauge*, BYTE, deviceGauge*, DOUBLE, DOUBLE) \
{ \
    appCrucial(_T("CFieldGaugeKernel: CalculateForceAnisotropy is not supported for " groupName " gauge group\n")); \
} \
static void CalculateAllStaples(const deviceGauge*, BYTE, deviceGauge* const *) \
{ \
    appCrucial(_T("CFieldGaugeKernel: CalculateAllStaples is not supported for " groupName " gauge group\n")); \
} \
static void CalculateAllPlaquttes(const deviceGauge*, BYTE, deviceGauge* const*) \
{ \
    appCrucial(_T("CFieldGaugeKernel: CalculateAllPlaquttes is not supported for " groupName " gauge group\n")); \
} \
static void CalculateFmunu(const deviceGauge*, BYTE, deviceGauge*) \
{ \
    appCrucial(_T("CFieldGaugeKernel: CalculateFmunu is not supported for " groupName " gauge group\n")); \
} \
static void CacheKSRotationGaugeBuffer(const deviceGauge*, BYTE, deviceGauge*, const Real*, Real, UBOOL) \
{ \
    appCrucial(_T("CFieldGaugeKernel: CacheKSRotationGaugeBuffer is not supported for " groupName " gauge group\n")); \
} \
static DOUBLE CalculatePlaqutteEnergyAnisotropy(const deviceGauge*, BYTE, DOUBLE, DOUBLE) \
{ \
    appCrucial(_T("CFieldGaugeKernel: CalculatePlaqutteEnergyAnisotropy is not supported for " groupName " gauge group\n")); \
    return F(0.0); \
} \
static DOUBLE CalculatePlaqutteEnergyUseClover(const deviceGauge*, BYTE, DOUBLE) \
{ \
    appCrucial(_T("CFieldGaugeKernel: CalculatePlaqutteEnergyUseClover is not supported for " groupName " gauge group\n")); \
    return F(0.0); \
} \
static DOUBLE CalculatePlaqutteEnergyUseCloverAnisotropy(const deviceGauge*, BYTE, DOUBLE, DOUBLE) \
{ \
    appCrucial(_T("CFieldGaugeKernel: CalculatePlaqutteEnergyUseCloverAnisotropy is not supported for " groupName " gauge group\n")); \
    return F(0.0); \
} \
static DOUBLE CalculatePlaqutteEnergyUsingStaple(const deviceGauge*, BYTE, DOUBLE, const deviceGauge*) \
{ \
    appCrucial(_T("CFieldGaugeKernel: CalculatePlaqutteEnergyUsingStaple is not supported for " groupName " gauge group\n")); \
    return F(0.0); \
} \
static void NaikForce(const CFieldGauge*, CFieldGauge*) \
{ \
    appCrucial(_T("CFieldGaugeKernel: NaikForce is not supported for " groupName " gauge group\n")); \
}

// Partial specialization for Z_N: only staple and energy are valid.
// All other methods (force, Fmunu, clover, anisotropy, rotation, etc.) are
// mathematical nonsense for discrete groups — call appCrucial to error out.
template<INT N>
class CFieldGaugeKernel<deviceZN<N>, 1>
{
public:
    // --- valid operations ---
    static void CalculateOnlyStaple(const deviceZN<N>* deviceData, BYTE byFieldId, deviceZN<N>* pStaple);
    static DOUBLE CalculatePlaqutteEnergy(const deviceZN<N>* deviceData, BYTE byFieldId, DOUBLE betaOverN);

    _DEF_GAUGE_KERNEL_ERROR_STUBS(deviceZN<N>, 1, "Z_N")
};

// Partial specialization for D_N (dihedral group): same matrixN=2 as SU2
// but only staple and energy are supported (discrete group).
template<INT N>
class CFieldGaugeKernel<deviceDN<N>, 2>
{
public:
    static void CalculateOnlyStaple(const deviceDN<N>* deviceData, BYTE byFieldId, deviceDN<N>* pStaple);
    static DOUBLE CalculatePlaqutteEnergy(const deviceDN<N>* deviceData, BYTE byFieldId, DOUBLE betaOverN);

    _DEF_GAUGE_KERNEL_ERROR_STUBS(deviceDN<N>, 2, "D_N")
};

#if !_CLG_WIN
//Deepseek said put them also in .h will speed up linker
extern template class CFieldGaugeKernel<CLGComplex, 1>;
extern template class CFieldGaugeKernel<deviceSU2, 2>;
extern template class CFieldGaugeKernel<deviceSU3, 3>;

#if _CLG_SU4_GAUGE
extern template class CFieldGaugeKernel<deviceSU4, 4>;
#endif
#if _CLG_SU5_GAUGE
extern template class CFieldGaugeKernel<deviceSU5, 5>;
#endif
#if _CLG_SU6_GAUGE
extern template class CFieldGaugeKernel<deviceSU6, 6>;
#endif
#if _CLG_SU7_GAUGE
extern template class CFieldGaugeKernel<deviceSU7, 7>;
#endif
#if _CLG_SU8_GAUGE
extern template class CFieldGaugeKernel<deviceSU8, 8>;
#endif
#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_KERNEL_H_

//=============================================================================
// END OF FILE
//=============================================================================