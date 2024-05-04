//=============================================================================
// FILENAME : CFieldFermionKSSU3GammaEM.h
// 
// DESCRIPTION:
// This is a helper to implement the condensations
// Do not use this to simulate unless you know what this is
//
// REVISION:
//  [09/28/2022 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONKSSU3GAMMAEM_H_
#define _CFIELDFERMIONKSSU3GAMMAEM_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3GammaEM)

class CLGAPI CFieldFermionKSSU3GammaEM : public CFieldFermionKSSU3Gamma
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3GammaEM)

public:

    CFieldFermionKSSU3GammaEM();

protected:

    void DerivateD0(void* pForce, const void* pGaugeBuffer) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void ApplyGammaKSS(const CFieldGauge* pGauge, EGammaMatrix eGamma) override;

public:

    void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString& tab) const override;

    Real m_fCharge;
    BYTE m_byEMFieldID;

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

    /**
     * devicePathBuffer must be larger than 4
     */
    static void GammaKSForceEM(
        void* pForce,
        const void* pGaugeBuffer,
        const void* pEMFieldBuffer,
        Real fCharge,
        const deviceSU3Vector* const* pRationalFields,
        const Real* pRationalNumerator,
        UINT uiRationalDegree,
        Real fCoeff,
        EGammaMatrix eGamma,
        INT* devicePathBuffer,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    static void DOperatorEM(
        void* pTargetBuffer,
        const void* pBuffer,
        const void* pGaugeBuffer,
        const void* pEMFieldBuffer,
        Real f2am,
        Real fCharge,
        UBOOL bShiftCenter,
        UBOOL bDagger,
        EOperatorCoefficientType eOCT,
        Real fRealCoeff,
        CLGComplex cCmpCoeff,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    static void KSForceEM(
        void* pForce,
        const void* pGaugeBuffer,
        const void* pEMFieldBuffer,
        Real fCharge,
        const deviceSU3Vector* const* pRationalFields,
        const Real* pRationalNumerator,
        UINT uiRationalDegree,
        BYTE byFieldID);

};

#pragma region device functions


static __device__ __inline__ deviceSU3 _devicePlaneDiagonalEM(
    const deviceSU3* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataReal,
    Real fCharge,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    SBYTE dim1, SBYTE dim2)
{
    INT dir1[2];

    dir1[0] = dim1; dir1[1] = dim2;
    deviceSU3 sRet(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 2, byFieldId, dir1));

    dir1[0] = dim2; dir1[1] = dim1;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 2, byFieldId, dir1));

    sRet.MulReal(F(0.5));
    return sRet;
}

/**
 * dim1, 2, 3 =
 * 1: x, -1: -x
 * 2: y, -2: -y
 * 3: z, -3: -z
 * 4: t, -4: -t
 */
static __device__ __inline__ deviceSU3 _deviceCubicDiagonalEM(
    const deviceSU3* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataReal,
    Real fCharge,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    SBYTE dim1, SBYTE dim2, SBYTE dim3)
{
    INT dir1[3];

    dir1[0] = dim1; dir1[1] = dim2; dir1[2] = dim3;
    deviceSU3 sRet(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim1; dir1[1] = dim3; dir1[2] = dim2;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim2; dir1[1] = dim1; dir1[2] = dim3;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim2; dir1[1] = dim3; dir1[2] = dim1;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim3; dir1[1] = dim1; dir1[2] = dim2;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim3; dir1[1] = dim2; dir1[2] = dim1;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 3, byFieldId, dir1));

    sRet.MulReal(OneOver6);
    return sRet;
}

static __device__ __inline__ deviceSU3 _deviceHyperCubicDiagonalEM(
    const deviceSU3* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataReal,
    Real fCharge,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    SBYTE dim1, SBYTE dim2, SBYTE dim3, SBYTE dim4)
{
    deviceSU3 sRet = deviceSU3::makeSU3Zero();
    const SBYTE dim1234[4] = { dim1, dim2, dim3, dim4 };
    INT dir1[4];
    SBYTE dim234[3];
    for (BYTE k = 0; k < 4; ++k)
    {
        dir1[0] = dim1234[k];
        for (BYTE k2 = 0; k2 < 3; ++k2)
        {
            BYTE idx = k2 + 1 + k;
            idx = idx > 3 ? (idx - 4) : idx;
            dim234[k2] = dim1234[idx];
        }

        dir1[1] = dim234[0]; dir1[2] = dim234[1]; dir1[3] = dim234[2];
        sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[0]; dir1[2] = dim234[2]; dir1[3] = dim234[1];
        sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[1]; dir1[2] = dim234[0]; dir1[3] = dim234[2];
        sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[1]; dir1[2] = dim234[2]; dir1[3] = dim234[0];
        sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[2]; dir1[2] = dim234[0]; dir1[3] = dim234[1];
        sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[2]; dir1[2] = dim234[1]; dir1[3] = dim234[0];
        sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 4, byFieldId, dir1));
    }

    sRet.MulReal(OneOver24);
    return sRet;
}


#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3GAMMAEM_H_

//=============================================================================
// END OF FILE
//=============================================================================