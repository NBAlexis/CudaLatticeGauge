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
//  [07/21/2024 nbale]
//=============================================================================
#ifndef _CFIELDFERMIONKST_KERNEL_H_
#define _CFIELDFERMIONKST_KERNEL_H_

__BEGIN_NAMESPACE

template<typename deviceVector, typename deviceGauge, INT vectorN>
class __DLL_EXPORT CFieldFermionKSTKernel
{
public:

    static UINT TestAntiHermitianS(BYTE byFieldId, const CFieldGauge* pGauge);

    static void DOperatorKS(UBOOL bEachSiteEta, deviceVector* pTargetBuffer, const deviceVector* pBuffer,
        const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
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

    static void OnlyMass(const deviceVector* pSource, deviceVector* pTarget, Real f2am, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void OneLinkS(const deviceVector* pSource, BYTE byFieldId, const deviceGauge* pGauge, BYTE byGaugeFieldId, deviceVector* pTarget, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff);

    static void OneLinkForceS(const deviceVector* pFermion, BYTE byFieldId, const deviceGauge* pGauge, BYTE byGaugeFieldId, deviceGauge* pForce, Real fCoefficient,
        const INT* pDevicePath, BYTE pathLength, BYTE byEtaIdx, const deviceVector* const* pRationalFields, const Real* pNumerator, UINT uiRationApproxOrder);

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
        const deviceVector* const* pRationalFields,
        const Real* pRationalNumerator,
        UINT uiRationalDegree,
        Real fCoeff,
        EGammaMatrix eGamma,
        INT* devicePathBuffer,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    static void VectorMultiplyMatrix(deviceVector** hostResBuffer, deviceVector** hostLeftBuffer, deviceVector** resBuffer, deviceVector** leftBuffer, TArray<CField*>& res, const TArray<CField*>& left, const CLGComplex* deviceMatrix, UINT uiDimX, UINT uiDimY);
};

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKST_KERNEL_H_

//=============================================================================
// END OF FILE
//=============================================================================