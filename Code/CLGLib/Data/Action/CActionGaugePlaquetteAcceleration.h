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

#if !_CLG_DOUBLEFLOAT
    DOUBLE Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;
#else
    Real Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;
#endif
    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;

    UBOOL CalculateForceOnGauge(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;
    void OnFinishTrajectory(UBOOL bAccepted) override;
    CCString GetInfos(const CCString &tab) const override;

    void SetBeta(Real fBeta);
    void SetG(Real fOmega);
    //void SetCenter(const SSmallInt4 &newCenter);
    //Real GetEnergyPerPlaqutte() const;

    Real m_fG;

protected:

#if !_CLG_DOUBLEFLOAT
    DOUBLE m_fLastEnergy;
    DOUBLE m_fNewEnergy;
#else
    Real m_fLastEnergy;
    Real m_fNewEnergy;
#endif
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
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));
    const SSmallInt4 n_p_mu_m_nu = _deviceSmallInt4OffsetC(n_m_nu, __fwd(mu));

    //n->nu
    //n_p_nu->mu
    //n_p_mu->nu +
    const SIndex& n__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIndex * _DC_Dir + nu];
    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + mu];
    SIndex n_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + nu];
    n_p_mu__nu.m_byTag = n_p_mu__nu.m_byTag ^ _kDaggerOrOpposite;

    deviceSU3 left(
        //_deviceGetSTTerm(
        //    byFieldId, pDeviceData, uiBigIndex, uiN_p_nu, uiN_p_mu, nu, mu, nu, 0, 0, 1
        //)
        _deviceGetSTTerm(byFieldId, pDeviceData, n__nu, n_p_nu__mu, n_p_mu__nu)
    );
    left.MulReal(_deviceGnAcc(sSite, fGsq));


    const UINT n_m_nu_bi = __bi(n_m_nu);
    const UINT n_m_nu_bi4 = n_m_nu_bi * _DC_Dir;
    SIndex n_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + nu];
    n_m_nu__nu.m_byTag = n_m_nu__nu.m_byTag ^ _kDaggerOrOpposite;
    const SIndex& n_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + mu];
    const SIndex& n_p_mu_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu_m_nu) + nu];

    const SIndex& site_n_m_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][n_m_nu_bi];
    
    left.Add(
        //_deviceGetSTTerm(
        //    byFieldId, pDeviceData, uiN_m_nu, uiN_m_nu, uiN_p_mu_m_nu, nu, mu, nu, 1, 0, 0
        //)
        _deviceGetSTTerm(byFieldId, pDeviceData, n_m_nu__nu, n_m_nu__mu, n_p_mu_m_nu__nu)        
        .MulRealC(_deviceGnAcc(__deviceSiteIndexToInt4(site_n_m_nu.m_uiSiteIndex), fGsq)));

    return left;
}

static __device__ __inline__ deviceSU3 _deviceStapleTerm_Acc_T(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite,
    Real fGsq,
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

    left.MulReal(_deviceGnAcc(sSite, fGsq));

    return left;
}

/**
* g t
*/
static __device__ __inline__ Real _deviceHi_Acc(
    const SSmallInt4& site, const SSmallInt4& site2)
{
    return (site.w + site2.w) * F(0.5);
}

/**
* [hi(n)+hi(n+nu)]S1  U(N+nu) U(N+mu)
* mu nu
* - +,
*/
static __device__ __inline__ deviceSU3 _deviceStapleS1_Acc(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));

    const UINT uiN_p_nu = __bi(n_p_nu);

    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_p_nu * _DC_Dir + mu];
    const SIndex& n_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + nu];

    const UINT uiSiteN_p_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceS1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));
    ret.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_mu__nu, byFieldId));

    ret.MulReal(_deviceHi_Acc(sSite, __deviceSiteIndexToInt4(uiSiteN_p_nu)));

    return ret;
}

/**
* [h(N) + h(n-nu)] S2 U(n-nu)U(n+mu-nu)
* mu nu
* + +
*/
static __device__ __inline__ deviceSU3 _deviceStapleS2_Acc(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite, __bck(mu));
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, __fwd(nu));

    const UINT uiN_m_nu = __bi(n_m_mu);

    const SIndex& n_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_m_nu * _DC_Dir + mu];
    const SIndex& n_m_mu_p_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_mu_p_nu) + nu];

    const UINT uiSiteN_m_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_m_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceS2(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_mu__mu, byFieldId));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_mu_p_nu__nu, byFieldId));

    ret.MulReal(_deviceHi_Acc(sSite, __deviceSiteIndexToInt4(uiSiteN_m_nu)));

    return ret;
}

/**
* [h(N+mu) + h(N+mu+nu)]U(n) U(n+nu) S3
* nu mu
* + +
*/
static __device__ __inline__ deviceSU3 _deviceStapleS3_Acc(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_mu_p_nu = _deviceSmallInt4OffsetC(n_p_mu, __fwd(nu));
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));

    const UINT uiN_p_mu = __bi(n_p_mu);
    const UINT uiN_p_mu_p_nu = __bi(n_p_mu_p_nu);

    const SIndex& n__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIndex * _DC_Dir + mu];
    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + mu];

    const UINT uiSiteN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_mu].m_uiSiteIndex;
    const UINT uiSiteN_p_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_mu_p_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n__nu, byFieldId));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));
    ret.Mul(_deviceS3(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

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
static __device__ __inline__ deviceSU3 _deviceStapleS4_Acc(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));
    const SSmallInt4 n_p_mu_m_nu = _deviceSmallInt4OffsetC(n_m_nu, __fwd(mu));
    
    const UINT uiN_p_mu = __bi(n_p_mu);
    const UINT uiN_p_mu_m_nu = __bi(n_p_mu_m_nu);
    const UINT n_m_nu_bi4 = __bi4(n_m_nu);

    const SIndex& n_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + mu];
    const SIndex& n_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + nu];

    const UINT uiSiteN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_mu].m_uiSiteIndex;
    const UINT uiSiteN_p_mu_m_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_mu_m_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu__nu, byFieldId));
    ret.DaggerMul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu__mu, byFieldId));
    ret.Mul(_deviceS4(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

    ret.MulReal(_deviceHi_Acc(__deviceSiteIndexToInt4(uiSiteN_p_mu),
        __deviceSiteIndexToInt4(uiSiteN_p_mu_m_nu)));

    return ret;
}

/**
* [h(n+mu) + h(n+mu+nu)] U(n) T1 U(n+nu)
* mu mu, + -
*
*/
static __device__ __inline__ deviceSU3 _deviceStapleT1_Acc(BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));
    const SSmallInt4 n_p_mu_p_nu = _deviceSmallInt4OffsetC(n_p_mu, __fwd(nu));

    const UINT uiN_p_mu = __bi(n_p_mu);
    const UINT uiN_p_mu_p_nu = __bi(n_p_mu_p_nu);

    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + mu];
    const SIndex& n__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIndex * _DC_Dir + mu];

    const UINT uiSiteN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_mu].m_uiSiteIndex;
    const UINT uiSiteN_p_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_mu_p_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n__mu, byFieldId));
    ret.Mul(_deviceT1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));

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
static __device__ __inline__ deviceSU3 _deviceStapleT2_Acc(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite, __bck(mu));
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, __fwd(nu));

    const UINT uiN_m_mu = __bi(n_m_mu);
    const UINT uiN_m_mu_p_nu = __bi(n_m_mu_p_nu);

    const SIndex& n_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_m_mu * _DC_Dir + mu];
    const SIndex& n_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_m_mu_p_nu * _DC_Dir + mu];

    const UINT uiSiteN_m_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_m_mu].m_uiSiteIndex;
    const UINT uiSiteN_m_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_m_mu_p_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_mu__mu, byFieldId));
    ret.DaggerMul(_deviceT2(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_mu_p_nu__mu, byFieldId));

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
static __device__ __inline__ deviceSU3 _deviceStapleChairTerm1_Acc(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    deviceSU3 ret(_deviceStapleS1_Acc(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleS2_Acc(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleS3_Acc(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleS4_Acc(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    return ret;
}

/**
* i = 0, 1, 2 for coefficient
* _deviceChairTerm1,2,3 for partial mu, nu, rho
* It is (1/8) * (T1+T2 + T1(mu<->rho) + T2(mu<->rho))
*/
static __device__ __inline__ deviceSU3 _deviceStapleChairTerm2_Acc(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    deviceSU3 ret(_deviceStapleT1_Acc(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleT2_Acc(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Add(_deviceStapleT1_Acc(byFieldId, pDeviceData, sSite, uiBigIndex, rho, nu, mu));
    ret.Add(_deviceStapleT2_Acc(byFieldId, pDeviceData, sSite, uiBigIndex, rho, nu, mu));
    return ret;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_ACCELERATION_H_

//=============================================================================
// END OF FILE
//=============================================================================