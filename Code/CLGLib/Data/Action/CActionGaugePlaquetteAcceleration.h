//=============================================================================
// FILENAME : CActionGaugePlaquetteAcceleration.h
// 
// DESCRIPTION:
// 
// Periodic boundary is assumed
//
// REVISION:
//  [07/27/2020 nbale]
//=============================================================================

#ifndef _CACTIONGAUGEPLAQUETTE_ACCELERATION_H_
#define _CACTIONGAUGEPLAQUETTE_ACCELERATION_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteAcceleration)

class CLGAPI CActionGaugePlaquetteAcceleration : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteAcceleration)
public:

    CActionGaugePlaquetteAcceleration();

    Real Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable) override;
    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;

    UBOOL CalculateForceOnGauge(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;
    void OnFinishTrajectory(UBOOL bAccepted) override;
    CCString GetInfos(const CCString &tab) const override;

    void SetBeta(Real fBeta);
    static void SetG(Real fOmega);
    //void SetCenter(const SSmallInt4 &newCenter);
    //Real GetEnergyPerPlaqutte() const;

protected:

    Real m_fLastEnergy;
    Real m_fNewEnergy;
    Real m_fBetaOverN;
    UINT m_uiPlaqutteCount;
};

#pragma region Device functions

//================= Put those device functions to header file because we will use them ==============

/**
 * This is 0.5 g^2 (t^2 + (t+1)^2), note that when t+1=Lt, it is 0.5 g^2 t^2
 */
static __device__ __inline__ Real _deviceGnAcc(const SSmallInt4& sSite, Real fGsq)
{
    if (sSite.w == static_cast<SBYTE>(_DC_Lt) - 1)
    {
        return F(0.5) * sSite.w * sSite.w  * fGsq;
    }
    return F(0.5) * (F(2.0) * sSite.w * sSite.w + F(2.0) * sSite.w + F(1.0)) * fGsq;
}

/**
 * If mu = 0, or 1, it is (1+g(n)) staple
 * If mu = 2, it is usual staple
 */
static __device__ __inline__ deviceSU3 _deviceStapleTerm_Acc_XY(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite,
    Real fGsq,
    BYTE byFieldId,
    UINT uiBigIndex, BYTE mu, BYTE nu)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;
    const UINT uiBIDir2 = uiBigIndex * uiDir2;
    const UINT uiN_p_mu = __idx->m_pWalkingTable[uiBIDir2 + mu + uiDir1];
    const UINT uiN_p_nu = __idx->m_pWalkingTable[uiBIDir2 + nu + uiDir1];
    const UINT uiN_m_nu = __idx->m_pWalkingTable[uiBIDir2 + nu];
    const UINT uiN_p_mu_m_nu = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + mu + uiDir1];
    const SSmallInt4 sSiteN_m_nu = __deviceSiteIndexToInt4(
        __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_m_nu].m_uiSiteIndex
    );
    deviceSU3 left(
        _deviceGetSTTerm(
            pDeviceData, uiBigIndex, uiN_p_nu, uiN_p_mu, nu, mu, nu, 0, 0, 1
        ));
    left.MulReal(_deviceGnAcc(sSite, fGsq));

    left.Add(_deviceGetSTTerm(
            pDeviceData, uiN_m_nu, uiN_m_nu, uiN_p_mu_m_nu, nu, mu, nu, 1, 0, 0
        ).MulRealC(_deviceGnAcc(sSiteN_m_nu, fGsq)));

    return left;
}

static __device__ __inline__ deviceSU3 _deviceStapleTerm_Acc_T(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite,
    Real fGsq,
    UINT uiBigIndex, BYTE mu, BYTE nu)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;
    const UINT uiBIDir2 = uiBigIndex * uiDir2;
    const UINT uiN_p_mu = __idx->m_pWalkingTable[uiBIDir2 + mu + uiDir1];
    const UINT uiN_p_nu = __idx->m_pWalkingTable[uiBIDir2 + nu + uiDir1];
    const UINT uiN_m_nu = __idx->m_pWalkingTable[uiBIDir2 + nu];
    const UINT uiN_p_mu_m_nu = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + mu + uiDir1];

    deviceSU3 left(
        _deviceGetSTTerm(
            pDeviceData, uiBigIndex, uiN_p_nu, uiN_p_mu, nu, mu, nu, 0, 0, 1
        ));
    left.Add(_deviceGetSTTerm(
        pDeviceData, uiN_m_nu, uiN_m_nu, uiN_p_mu_m_nu, nu, mu, nu, 1, 0, 0
    ));

    left.MulReal(_deviceGnAcc(sSite, fGsq));

    return left;
}

/**
* g t
*/
static __device__ __inline__ Real _deviceHi_Acc(
    const SSmallInt4& site, const SSmallInt4& site2)
{
    return (site.w + site2.w);
}

/**
* [hi(n)+hi(n+nu)]S1  U(N+nu) U(N+mu)
* mu nu
* - +,
*/
static __device__ __inline__ deviceSU3 _deviceStapleS1_Acc(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;
    const UINT uiBIDir2 = uiBigIndex * uiDir2;

    const UINT uiN_p_mu = __idx->m_pWalkingTable[uiBIDir2 + mu + uiDir1];
    const UINT uiN_p_nu = __idx->m_pWalkingTable[uiBIDir2 + nu + uiDir1];

    const UINT uiSiteN_p_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceS1(pDeviceData, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_p_nu, mu));
    ret.MulDagger(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_p_mu, nu));

    ret.MulReal(_deviceHi_Acc(sSite, __deviceSiteIndexToInt4(uiSiteN_p_nu)));

    return ret;
}

/**
* [h(N) + h(n-nu)] S2 U(n-nu)U(n+mu-nu)
* mu nu
* + +
*/
static __device__ __inline__ deviceSU3 _deviceStapleS2_Acc(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;

    const UINT uiN_m_nu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + nu];
    const UINT uiN_m_nu_p_mu = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + mu + uiDir1];

    const UINT uiSiteN_m_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceS2(pDeviceData, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_nu, mu));
    ret.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_nu_p_mu, nu));

    ret.MulReal(_deviceHi_Acc(sSite, __deviceSiteIndexToInt4(uiSiteN_m_nu)));

    return ret;
}

/**
* [h(N+mu) + h(N+mu+nu)]U(n) U(n+nu) S3
* nu mu
* + +
*/
static __device__ __inline__ deviceSU3 _deviceStapleS3_Acc(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;

    const UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu + uiDir1];
    const UINT uiN_p_mu_p_nu = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + nu + uiDir1];
    const UINT uiN_p_nu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + nu + uiDir1];

    const UINT uiSiteN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu].m_uiSiteIndex;
    const UINT uiSiteN_p_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu_p_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceGetGaugeBCSU3Dir(pDeviceData, uiBigIndex, nu));
    ret.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_p_nu, mu));
    ret.Mul(_deviceS3(pDeviceData, uiBigIndex, mu, nu, rho));

    ret.MulReal(_deviceHi_Acc(
        __deviceSiteIndexToInt4(uiSiteN_p_mu),
        __deviceSiteIndexToInt4(uiSiteN_p_mu_p_nu)));

    return ret;

}

/**
* [h(N+mu) + h(N+mu-nu)] U(n-nu) U(n-nu) S4
* nu mu
* - +
*/
static __device__ __inline__ deviceSU3 _deviceStapleS4_Acc(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;

    const UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu + uiDir1];
    const UINT uiN_p_mu_m_nu = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + nu];
    const UINT uiN_m_nu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + nu];

    const UINT uiSiteN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu].m_uiSiteIndex;
    const UINT uiSiteN_p_mu_m_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu_m_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_nu, nu));
    ret.DaggerMul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_nu, mu));
    ret.Mul(_deviceS4(pDeviceData, uiBigIndex, mu, nu, rho));

    ret.MulReal(_deviceHi_Acc(__deviceSiteIndexToInt4(uiSiteN_p_mu),
        __deviceSiteIndexToInt4(uiSiteN_p_mu_m_nu)));

    return ret;
}

/**
* [h(n+mu) + h(n+mu+nu)] U(n) T1 U(n+nu)
* mu mu, + -
*
*/
static __device__ __inline__ deviceSU3 _deviceStapleT1_Acc(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;

    const UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu + uiDir1];
    const UINT uiN_p_nu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + nu + uiDir1];
    const UINT uiN_p_mu_p_nu = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + nu + uiDir1];

    const UINT uiSiteN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu].m_uiSiteIndex;
    const UINT uiSiteN_p_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu_p_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceGetGaugeBCSU3Dir(pDeviceData, uiBigIndex, mu));
    ret.Mul(_deviceT1(pDeviceData, uiBigIndex, mu, nu, rho));
    ret.MulDagger(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_p_nu, mu));

    ret.MulReal(_deviceHi_Acc(
        __deviceSiteIndexToInt4(uiSiteN_p_mu),
        __deviceSiteIndexToInt4(uiSiteN_p_mu_p_nu)));

    return ret;
}

/**
* [h(n-mu) + h(n-mu+nu)] U(n-mu) T2 U(n+nu-mu)
* mu mu, - +
*
*/
static __device__ __inline__ deviceSU3 _deviceStapleT2_Acc(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;

    const UINT uiN_m_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu];
    const UINT uiN_m_mu_p_nu = __idx->m_pWalkingTable[uiN_m_mu * uiDir2 + nu + uiDir1];

    const UINT uiSiteN_m_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_mu].m_uiSiteIndex;
    const UINT uiSiteN_m_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_mu_p_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_mu, mu));
    ret.DaggerMul(_deviceT2(pDeviceData, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_mu_p_nu, mu));

    ret.MulReal(_deviceHi_Acc(
        __deviceSiteIndexToInt4(uiSiteN_m_mu),
        __deviceSiteIndexToInt4(uiSiteN_m_mu_p_nu)));

    return ret;
}

/**
* i = 0, 1, 2 for coefficient
* _deviceChairTerm1,2,3 for partial mu, nu, rho
* For partial mu, the staple is (1/8)(s1+s2+s3+s4)
*/
static __device__ __inline__ deviceSU3 _deviceStapleChairTerm1_Acc(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    deviceSU3 ret(_deviceStapleS1_Acc(pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleS2_Acc(pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleS3_Acc(pDeviceData, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleS4_Acc(pDeviceData, uiBigIndex, mu, nu, rho));
    return ret;
}

/**
* i = 0, 1, 2 for coefficient
* _deviceChairTerm1,2,3 for partial mu, nu, rho
* It is (1/8) * (T1+T2 + T1(mu<->rho) + T2(mu<->rho))
*/
static __device__ __inline__ deviceSU3 _deviceStapleChairTerm2_Acc(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    deviceSU3 ret(_deviceStapleT1_Acc(pDeviceData, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleT2_Acc(pDeviceData,  uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleT1_Acc(pDeviceData, uiBigIndex, rho, nu, mu));
    ret.Add(_deviceStapleT2_Acc(pDeviceData, uiBigIndex, rho, nu, mu));
    return ret;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_ACCELERATION_H_

//=============================================================================
// END OF FILE
//=============================================================================