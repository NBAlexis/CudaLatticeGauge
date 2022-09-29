//=============================================================================
// FILENAME : CFieldFermionKSSU3Gamma.h
// 
// DESCRIPTION:
// This is a helper to implement the condensations
// Do not use this to simulate unless you know what this is
//
// REVISION:
//  [09/10/2022 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONKSSU3GAMMA_H_
#define _CFIELDFERMIONKSSU3GAMMA_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3Gamma)

class CLGAPI CFieldFermionKSSU3Gamma : public CFieldFermionKSSU3
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3Gamma)

public:

    CFieldFermionKSSU3Gamma();
    ~CFieldFermionKSSU3Gamma();

    void DerivateD0(void* pForce, const void* pGaugeBuffer) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

    void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString& tab) const override;

    Real m_fCoeffGamma1;
    Real m_fCoeffGamma2;
    Real m_fCoeffGamma3;
    Real m_fCoeffGamma4;
    Real m_fCoeffGamma5;
    Real m_fCoeffGamma51;
    Real m_fCoeffGamma52;
    Real m_fCoeffGamma53;
    Real m_fCoeffGamma54;
    Real m_fCoeffSigma12;
    Real m_fCoeffSigma13;
    Real m_fCoeffSigma14;
    Real m_fCoeffSigma23;
    Real m_fCoeffSigma24;
    Real m_fCoeffSigma34;

    INT* m_pDevicePathBuffer;

    /**
     * This is for simulation, 2a is already multiplied.
     * 2a qbar Gamma q
     * for example, gamma_i  -> 1 x chichi
     *              sigma ij -> 1/2 x chichi
     *              gamma 5i -> 1/4 x chichi
     *              gamma 5  -> 1/8 x chichi
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
        const deviceSU3Vector* const* pRationalFields,
        const Real* pRationalNumerator,
        UINT uiRationalDegree,
        Real fCoeff,
        EGammaMatrix eGamma,
        INT* devicePathBuffer,
        BYTE byFieldID,
        BYTE byGaugeFieldID);

    /**
     * every time set gamma coefficient, update the parameters of pooled
     */
    void UpdatePooledParamters() const;

};

#pragma region device functions

static __device__ __inline__ SBYTE _deviceEta2(UINT uiEta, BYTE i, BYTE j)
{
    return ((uiEta >> i) + (uiEta >> j)) & 1;
}

/**
 * eta xyz, eta yzt, eta xyt, ...
 * for 1, 3 there is a minus sign
 * missingDir:
 * 3 - xyz x:1  y:(-1)^x z:(-1)^(x+y)             res: (-1)^y
 * 2 - xyt x:1  y:(-1)^x t:(-1)^(x+y+z)           res: (-1)^(y+z)
 * 0 - yzt y:(-1)^x z:(-1)^(x+y) t:(-1)^(x+y+z)   res: (-1)^(x+z)
 * 1 - xzt x:1  z:(-1)^(x+y) t:(-1)^(x+y+z)       res: (-1)^z
 * 
 */
static __device__ __inline__ SBYTE _deviceEta3(const SSmallInt4& sSite, BYTE missingDir)
{
    switch (missingDir)
    {
    case 3:
        return (sSite.y + 1) & 1;
    case 2:
        return (sSite.y + sSite.z) & 1;
    case 0:
        return (sSite.x + sSite.z) & 1;
    default:
        return (sSite.z + 1) & 1;
    }
}

static __device__ __inline__ deviceSU3 _devicePlaneDiagonal(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    SBYTE dim1, SBYTE dim2)
{
    INT dir1[2];

    dir1[0] = dim1; dir1[1] = dim2;
    deviceSU3 sRet(_deviceLink(pDeviceData, sStartSite, 2, byFieldId, dir1));

    dir1[0] = dim2; dir1[1] = dim1;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 2, byFieldId, dir1));

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
static __device__ __inline__ deviceSU3 _deviceCubicDiagonal(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    SBYTE dim1, SBYTE dim2, SBYTE dim3)
{
    INT dir1[3];

    dir1[0] = dim1; dir1[1] = dim2; dir1[2] = dim3;
    deviceSU3 sRet(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim1; dir1[1] = dim3; dir1[2] = dim2;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim2; dir1[1] = dim1; dir1[2] = dim3;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim2; dir1[1] = dim3; dir1[2] = dim1;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim3; dir1[1] = dim1; dir1[2] = dim2;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim3; dir1[1] = dim2; dir1[2] = dim1;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    sRet.MulReal(OneOver6);
    return sRet;
}

static __device__ __inline__ deviceSU3 _deviceHyperCubicDiagonal(
    const deviceSU3* __restrict__ pDeviceData,
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
        sRet.Add(_deviceLink(pDeviceData, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[0]; dir1[2] = dim234[2]; dir1[3] = dim234[1];
        sRet.Add(_deviceLink(pDeviceData, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[1]; dir1[2] = dim234[0]; dir1[3] = dim234[2];
        sRet.Add(_deviceLink(pDeviceData, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[1]; dir1[2] = dim234[2]; dir1[3] = dim234[0];
        sRet.Add(_deviceLink(pDeviceData, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[2]; dir1[2] = dim234[0]; dir1[3] = dim234[1];
        sRet.Add(_deviceLink(pDeviceData, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[2]; dir1[2] = dim234[1]; dir1[3] = dim234[0];
        sRet.Add(_deviceLink(pDeviceData, sStartSite, 4, byFieldId, dir1));
    }

    sRet.MulReal(OneOver24);
    return sRet;
}


#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3GAMMA_H_

//=============================================================================
// END OF FILE
//=============================================================================