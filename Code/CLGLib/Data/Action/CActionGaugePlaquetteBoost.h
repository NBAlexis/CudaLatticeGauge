//=============================================================================
// FILENAME : CActionGaugePlaquetteBoost.h
// 
// DESCRIPTION:
// 
// Periodic boundary is assumed
//
// REVISION:
//  [08/03/2020 nbale]
//=============================================================================

#ifndef _CACTIONGAUGEPLAQUETTE_BOOST_H_
#define _CACTIONGAUGEPLAQUETTE_BOOST_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteBoost)

class CLGAPI CActionGaugePlaquetteBoost : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteBoost)
public:

    CActionGaugePlaquetteBoost();

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;

    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;
    CCString GetInfos(const CCString &tab) const override;

    void SetBeta(Real fBeta);
    static void SetG(Real fOmega);
    //void SetCenter(const SSmallInt4 &newCenter);
    //Real GetEnergyPerPlaqutte() const;

protected:


    UINT m_uiPlaqutteCount;
};

#pragma region Device functions

//================= Put those device functions to header file because we will use them ==============


/**
 * 
 */
static __device__ __inline__ deviceSU3 _deviceStapleTerm_Boost(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));
    const SSmallInt4 n_p_mu_m_nu = _deviceSmallInt4OffsetC(n_m_nu, __fwd(mu));

    const SIndex& n__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIndex * _DC_Dir + nu];
    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + mu];
    SIndex n_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + nu];
    n_p_mu__nu.m_byTag = n_p_mu__nu.m_byTag ^ _kDaggerOrOpposite;

    const UINT n_m_nu_bi4 = __bi4(n_m_nu);
    SIndex n_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + nu];
    n_m_nu__nu.m_byTag = n_m_nu__nu.m_byTag ^ _kDaggerOrOpposite;
    const SIndex& n_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + mu];
    const SIndex& n_p_mu_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu_m_nu) + nu];

    deviceSU3 left(
        //_deviceGetSTTerm(
        //    pDeviceData, uiBigIndex, uiN_p_nu, uiN_p_mu, nu, mu, nu, 0, 0, 1
        //)
        _deviceGetSTTerm(byFieldId, pDeviceData, n__nu, n_p_nu__mu, n_p_mu__nu)
    );

    left.Add(
    //    _deviceGetSTTerm(
    //    pDeviceData, uiN_m_nu, uiN_m_nu, uiN_p_mu_m_nu, nu, mu, nu, 1, 0, 0
    //)
        _deviceGetSTTerm(byFieldId, pDeviceData, n_m_nu__nu, n_m_nu__mu, n_p_mu_m_nu__nu)
    );
    return left;
}

/**
* [hi(n)+hi(n+nu)]S1  U(N+nu) U(N+mu)
* mu nu
* - +,
*/
static __device__ __inline__ deviceSU3 _deviceStapleS1_Boost(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));

    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + mu];
    const SIndex& n_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + nu];

    deviceSU3 ret(_deviceS1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));
    ret.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_mu__nu, byFieldId));

    return ret;
}

/**
* [h(N) + h(n-nu)] S2 U(n-nu)U(n+mu-nu)
* mu nu
* + +
*/
static __device__ __inline__ deviceSU3 _deviceStapleS2_Boost(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));

    const SIndex& n_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_nu) + mu];
    const SIndex& n_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_nu_p_mu) + nu];

    deviceSU3 ret(_deviceS2(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu__mu, byFieldId));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu_p_mu__nu, byFieldId));

    return ret;
}

/**
* [h(N+mu) + h(N+mu+nu)]U(n) U(n+nu) S3
* nu mu
* + +
*/
static __device__ __inline__ deviceSU3 _deviceStapleS3_Boost(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));

    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + mu];
    const SIndex& n__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIndex * _DC_Dir + nu];

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n__nu, byFieldId));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));
    ret.Mul(_deviceS3(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

    return ret;

}

/**
* [h(N+mu) + h(N+mu-nu)] U(n-nu) U(n-nu) S4
* nu mu
* - +
*/
static __device__ __inline__ deviceSU3 _deviceStapleS4_Boost(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));
    const UINT n_m_nu_bi4 = __bi4(n_m_nu);

    const SIndex& n_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + mu];
    const SIndex& n_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + nu];

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu__nu, byFieldId));
    ret.DaggerMul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu__mu, byFieldId));
    ret.Mul(_deviceS4(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

    return ret;
}

/**
* [h(n+mu) + h(n+mu+nu)] U(n) T1 U(n+nu)
* mu mu, + -
*
*/
static __device__ __inline__ deviceSU3 _deviceStapleT1_Boost(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));

    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + mu];
    const SIndex& n__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIndex * _DC_Dir + mu];

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n__mu, byFieldId));
    ret.Mul(_deviceT1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));

    return ret;
}

/**
* [h(n-mu) + h(n-mu+nu)] U(n-mu) T2 U(n+nu-mu)
* mu mu, - +
*
*/
static __device__ __inline__ deviceSU3 _deviceStapleT2_Boost(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite, __bck(mu));
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, __fwd(nu));

    const SIndex& n_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_mu) + mu];
    const SIndex& n_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_mu_p_nu) + mu];

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_mu__mu, byFieldId));
    ret.DaggerMul(_deviceT2(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_mu_p_nu__mu, byFieldId));

    return ret;
}

/**
* i = 0, 1, 2 for coefficient
* _deviceChairTerm1,2,3 for partial mu, nu, rho
* For partial mu, the staple is (1/8)(s1+s2+s3+s4)
*/
static __device__ __inline__ deviceSU3 _deviceStapleChairTerm1_Boost(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    deviceSU3 ret(_deviceStapleS1_Boost(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleS2_Boost(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleS3_Boost(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleS4_Boost(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    return ret;
}

/**
* i = 0, 1, 2 for coefficient
* _deviceChairTerm1,2,3 for partial mu, nu, rho
* It is (1/8) * (T1+T2 + T1(mu<->rho) + T2(mu<->rho))
*/
static __device__ __inline__ deviceSU3 _deviceStapleChairTerm2_Boost(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    deviceSU3 ret(_deviceStapleT1_Boost(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleT2_Boost(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleT1_Boost(byFieldId, pDeviceData, sSite, uiBigIndex, rho, nu, mu));
    ret.Add(_deviceStapleT2_Boost(byFieldId, pDeviceData, sSite, uiBigIndex, rho, nu, mu));
    return ret;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_ACCELERATION_H_

//=============================================================================
// END OF FILE
//=============================================================================