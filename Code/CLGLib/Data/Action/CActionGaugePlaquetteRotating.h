//=============================================================================
// FILENAME : CActionGaugePlaquetteRotating.h
// 
// DESCRIPTION:
// This is the class for rotating guage action
// Open boundary condition (identity Dirichlet boundary condition) is assumed 
// 
//
// REVISION:
//  [05/07/2019 nbale]
//=============================================================================

#ifndef _CACTIONGAUGEPLAQUETTE_ROTATING_H_
#define _CACTIONGAUGEPLAQUETTE_ROTATING_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionGaugePlaquetteRotating)

class CLGAPI CActionGaugePlaquetteRotating : public CAction
{
    __CLGDECLARE_CLASS(CActionGaugePlaquetteRotating)
public:

    CActionGaugePlaquetteRotating();

    Real Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable) override;
    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;

    UBOOL CalculateForceOnGauge(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;
    void OnFinishTrajectory(UBOOL bAccepted) override;
    CCString GetInfos(const CCString &tab) const override;

    void SetBeta(Real fBeta);
    void SetOmega(Real fOmega);
    void SetCenter(const SSmallInt4 &newCenter);
    //Real GetEnergyPerPlaqutte() const;

    Real m_fOmega;

protected:

    Real m_fLastEnergy;
    Real m_fNewEnergy;
    Real m_fBetaOverN;
    UINT m_uiPlaqutteCount;
};

//================= Put those device functions to header file because we will use them ==============

#pragma region device function

#pragma region Energy

#pragma region Plaqutte term

/**
* big index is the index of walking table.
* The plaqutte index may not be cached because n may out of boundary, so we calculate every one
* n, n+mu, n+nu, n
*/
static __device__ __inline__ deviceSU3 _device1PlaqutteTermPP(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;
    const UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + byMu + uiDir1];
    const UINT uiN_p_nu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + byNu + uiDir1];

    deviceSU3 u(_deviceGetGaugeBCSU3Dir(pDeviceData, uiBigIdx, byMu));
    u.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_p_mu, byNu));
    u.MulDagger(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_p_nu, byMu));
    u.MulDagger(_deviceGetGaugeBCSU3Dir(pDeviceData, uiBigIdx, byNu));

    return u;
}

/**
* U(-mu,nu) = U^+_{mu}(N-mu) U_{nu}(N-mu) U_{mu}(N-mu+nu) U^+_{nu}(N)
*/
static __device__ __inline__ deviceSU3 _device1PlaqutteTermMP(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;
    const UINT uiN_m_mu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + byMu];
    const UINT uiN_m_mu_p_nu = __idx->m_pWalkingTable[uiN_m_mu * uiDir2 + byNu + uiDir1];

    deviceSU3 u(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_mu, byMu));
    u.DaggerMul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_mu, byNu));
    u.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_mu_p_nu, byMu));
    u.MulDagger(_deviceGetGaugeBCSU3Dir(pDeviceData, uiBigIdx, byNu));

    return u;
}

/**
* U(mu,-nu) = U(N) U(N+mu-nu) U(N-nu) U(N-nu)
*/
static __device__ __inline__ deviceSU3 _device1PlaqutteTermPM(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;
    const UINT uiN_m_nu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + byNu];
    const UINT uiN_m_nu_p_mu = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + byMu + uiDir1];

    deviceSU3 u(_deviceGetGaugeBCSU3Dir(pDeviceData, uiBigIdx, byMu));
    u.MulDagger(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_nu_p_mu, byNu));
    u.MulDagger(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_nu, byMu));
    u.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_nu, byNu));

    return u;
}

/**
* U(-mu,-nu) = U(N-mu) U(N-mu-nu) U(N-mu-nu) U(N-nu)
*/
static __device__ __inline__ deviceSU3 _device1PlaqutteTermMM(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;
    const UINT uiN_m_nu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + byNu];
    const UINT uiN_m_mu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + byMu];
    const UINT uiN_m_nu_m_mu = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + byMu];

    //u1^+ u2^+ u3 u4
    //= (u2 u1)^+ u3 u4
    deviceSU3 u(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_nu_m_mu, byNu));
    u.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_mu, byMu));
    u.DaggerMul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_nu_m_mu, byMu));
    u.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_nu, byNu));

    return u;
}


/**
* Product of 3 terms
*/
static __device__ __inline__ deviceSU3 _deviceGetSTTerm(
    const deviceSU3* __restrict__ pDeviceData,
    UINT uiBIa, UINT uiBIb, UINT uiBIc,
    BYTE byDira, BYTE byDirb, BYTE byDirc,
    BYTE byDaggera, BYTE byDaggerb, BYTE byDaggerc)
{
    deviceSU3 ret(_deviceGetGaugeBCSU3Dir(pDeviceData, uiBIa, byDira));
    if (1 == byDaggera)
    {
        ret.Dagger();
    }

    if (1 == byDaggerb)
    {
        ret.MulDagger(_deviceGetGaugeBCSU3Dir(pDeviceData, uiBIb, byDirb));
    }
    else
    {
        ret.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiBIb, byDirb));
    }

    if (1 == byDaggerc)
    {
        ret.MulDagger(_deviceGetGaugeBCSU3Dir(pDeviceData, uiBIc, byDirc));
    }
    else
    {
        ret.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiBIc, byDirc));
    }

    return ret;
}

static __device__ __inline__ Real _device1PlaqutteTermReTr(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx)
{
    return _device1PlaqutteTermPP(pDeviceData, byMu, byNu, uiBigIdx).ReTr();
}

/**
* 3 - 1/4 Retr[ U_{mu,nu}(n)+U_{-mu,nu}(n)+U_{mu,-nu}(n)+U_{-mu,-nu}(n) ]
* = 3 - 1/4 Retr[ U_{mu,nu}(n)+U^+_{mu,nu}(n-mu)+U^+_{mu,nu}(n-nu)+U_{mu,nu}(n-mu-nu) ]
* = 3 - 1/4 Retr[ U_{mu,nu}(n)+U_{mu,nu}(n-mu)+U_{mu,nu}(n-nu)+U_{mu,nu}(n-mu-nu) ]
*/
static __device__ __inline__ Real _device4PlaqutteTerm(const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIndex)
{
    const UINT uiDir2 = _DC_Dir * 2;
    const UINT uiN_m_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + byMu];
    return F(3.0) - F(0.25) * (
        _device1PlaqutteTermReTr(pDeviceData, byMu, byNu, uiBigIndex)
        + _device1PlaqutteTermReTr(pDeviceData, byMu, byNu, uiN_m_mu)
        + _device1PlaqutteTermReTr(pDeviceData, byMu, byNu, __idx->m_pWalkingTable[uiBigIndex * uiDir2 + byNu])
        + _device1PlaqutteTermReTr(pDeviceData, byMu, byNu, __idx->m_pWalkingTable[uiN_m_mu * uiDir2 + byNu])
        );
}

#pragma endregion

#pragma region Chair term

/**
* (1/8) * Retr[()()+()()] = (1/8) * Retr[left+right]
* left(n) = Retr[(U_{a,b}(n)-U^+_{a,b}(n-a))(U_{b,c}(n)-U^+_{b,c}(n-c))]
* right(n) = Retr[(U^+_{a,b}(n-b)-U_{a,b}(n-a-b))(U^+_{b,c}(n-b)-U_{b,c}(n-b-c))]
*          = Retr[(U_{a,b}(n-b)-U^+_{a,b}(n-a-b))(U_{b,c}(n-b)-U^+_{b,c}(n-b-c))]
*          = left(n-b)
*/
static __device__ __inline__ Real _deviceChairTerm(const deviceSU3* __restrict__ pDeviceData,
    BYTE mu, BYTE nu, BYTE rho, UINT uiBigIndex)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;
    const UINT uiBigIdxDir2 = uiBigIndex * uiDir2;
    const UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIdxDir2 + mu + uiDir1];
    const UINT uiN_m_mu = __idx->m_pWalkingTable[uiBigIdxDir2 + mu];
    const UINT uiN_p_nu = __idx->m_pWalkingTable[uiBigIdxDir2 + nu + uiDir1];
    const UINT uiN_m_nu = __idx->m_pWalkingTable[uiBigIdxDir2 + nu];
    const UINT uiN_p_rho = __idx->m_pWalkingTable[uiBigIdxDir2 + rho + uiDir1];
    const UINT uiN_m_rho = __idx->m_pWalkingTable[uiBigIdxDir2 + rho];

    const UINT uiN_p_mu_m_nu = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + nu];
    const UINT uiN_m_mu_p_nu = __idx->m_pWalkingTable[uiN_m_mu * uiDir2 + nu + uiDir1];
    const UINT uiN_m_mu_m_nu = __idx->m_pWalkingTable[uiN_m_mu * uiDir2 + nu];
    const UINT uiN_m_rho_p_nu = __idx->m_pWalkingTable[uiN_m_rho * uiDir2 + nu + uiDir1];
    const UINT uiN_m_rho_m_nu = __idx->m_pWalkingTable[uiN_m_rho * uiDir2 + nu];
    const UINT uiN_m_nu_p_rho = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + rho + uiDir1];

    // U_{mu}(N) U_{nu}(N+mu) U^+_{mu}(n+nu)
    deviceSU3 term1(_deviceGetSTTerm(pDeviceData,
        uiBigIndex, uiN_p_mu, uiN_p_nu, mu, nu, mu, 0, 0, 1));

    //U^+_{mu}(N-mu) U_{nu}(N-mu) U_{mu}(N-mu+nu)
    term1.Sub(_deviceGetSTTerm(pDeviceData,
        uiN_m_mu, uiN_m_mu, uiN_m_mu_p_nu, mu, nu, mu, 1, 0, 0));

    // U_{rho}(N+nu) U^+_{nu}(N+rho) U^+_{rho}(N)
    deviceSU3 term2(_deviceGetSTTerm(pDeviceData,
        uiN_p_nu, uiN_p_rho, uiBigIndex, rho, nu, rho, 0, 1, 1));

    // U^+_{rho}(N+nu-rho) U^+_{nu}(N-rho) U_{rho}(N-rho)
    term2.Sub(_deviceGetSTTerm(pDeviceData,
        uiN_m_rho_p_nu, uiN_m_rho, uiN_m_rho, rho, nu, rho, 1, 1, 0));

    term1.Mul(term2);

    //pm mu, nu
    //U(mu,-nu) = U(N) U(N+mu-nu) U(N-nu) U(N-nu), 0110
    deviceSU3 term3(_deviceGetSTTerm(pDeviceData,
        uiBigIndex, uiN_p_mu_m_nu, uiN_m_nu, mu, nu, mu, 0, 1, 1));

    //mm
    //U(-mu, -nu) = U(N - mu) U(N - mu - nu) U(N - mu - nu) U(N - nu) 1100
    term3.Sub(_deviceGetSTTerm(pDeviceData,
        uiN_m_mu, uiN_m_mu_m_nu, uiN_m_mu_m_nu, mu, nu, mu, 1, 1, 0));

    //mp, nu, rho
    //mp = U(-mu,nu) = U^+_{mu}(N-mu) U_{nu}(N-mu) U_{mu}(N-mu+nu) U^+_{nu}(N)
    deviceSU3 term4(_deviceGetSTTerm(pDeviceData,
        uiN_m_nu, uiN_m_nu_p_rho, uiBigIndex, rho, nu, rho, 0, 0, 1));

    //mm nu rho
    //U(-mu, -nu) = U(N - mu) U(N - mu - nu) U(N - mu - nu) U(N - nu) 1100
    term4.Sub(_deviceGetSTTerm(pDeviceData,
        uiN_m_rho_m_nu, uiN_m_rho_m_nu, uiN_m_rho, rho, nu, rho, 1, 0, 0));

    term3.Mul(term4);

    term1.Add(term3);

    return term1.ReTr();
}

#pragma endregion

#pragma endregion

#pragma region Force

#pragma region Plaqutte term



/**
* g1=O^2(x^2+(x+1)^2)/2
* g2=O^2(y^2+(y+1)^2)/2
* g3=O^2(x^2+y^2)
* For identity Dirichlet boundary, if site is out of boundary, {I}_TA = 0
* So we do not care whether site is out of boundary
* Note that, for x+1, it dose NOT always mean x+1
* For g1, g2, site offset is x+1 site and y+1 site,
* for g3, sSiteOffset is not using
*/
static __device__ __inline__ Real _deviceGi(
    const SSmallInt4& sCenter,
    const SSmallInt4& sSite,
    const SSmallInt4& sSiteOffset,
    UINT uiSiteBI,
    UINT uiSiteOffsetBI,
    BYTE i,
    Real fOmegaSq)
{
    if (0 == i)
    {
        const Real fX = __idx->m_pDeviceIndexPositionToSIndex[1][uiSiteBI].IsDirichlet() ? F(0.0)
            : static_cast<Real>(sSite.x - sCenter.x);
        const Real fXp1 = __idx->m_pDeviceIndexPositionToSIndex[1][uiSiteOffsetBI].IsDirichlet() ? F(0.0)
            : static_cast<Real>(sSiteOffset.x - sCenter.x);
        return F(0.5) * fOmegaSq * (fX * fX + fXp1 * fXp1);
    }
    else if (1 == i)
    {
        const Real fY = __idx->m_pDeviceIndexPositionToSIndex[1][uiSiteBI].IsDirichlet() ? F(0.0)
            : static_cast<Real>(sSite.y - sCenter.y);
        const Real fYp1 = __idx->m_pDeviceIndexPositionToSIndex[1][uiSiteOffsetBI].IsDirichlet() ? F(0.0)
            : static_cast<Real>(sSiteOffset.y - sCenter.y);
        return F(0.5) * fOmegaSq * (fY * fY + fYp1 * fYp1);
    }
    const Real fX = __idx->m_pDeviceIndexPositionToSIndex[1][uiSiteBI].IsDirichlet() ? F(0.0)
        : static_cast<Real>(sSite.x - sCenter.x);
    const Real fY = __idx->m_pDeviceIndexPositionToSIndex[1][uiSiteBI].IsDirichlet() ? F(0.0)
        : static_cast<Real>(sSite.y - sCenter.y);
    return fOmegaSq * (fX * fX + fY * fY);
}

/**
* Coefficient = (f(n)+f(n+mu)+f(n+nu)+f(n+mu+nu))/4
* Simplfy: nu is always t direction, so f(n) = f(n+nu), f(n+mu) = f(n+mu+nu)
* Coefficient = (f(n)+f(n+mu))/2
* For 3 == mu, f(n) = f(n+mu)
* This is also true for Dirichlet boundary condition, only Dirichlet on X-Y direction is assumed
*/
static __device__ __inline__ Real _deviceFi(
    const SSmallInt4& sCenter,
    UINT uiN, BYTE i, BYTE mu, BYTE nu)
{
    const SSmallInt4 sN = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[1][uiN].m_uiSiteIndex);
    if (2 == i)
    {
        if (__idx->m_pDeviceIndexPositionToSIndex[1][uiN].IsDirichlet())
        {
            return F(0.0);
        }
        //only one case: U_{zt}, so no boundary
        const INT x1 = sN.x - sCenter.x;
        const INT y1 = sN.y - sCenter.y;
        return static_cast<Real>(x1 * x1 + y1 * y1);
    }

    //for U_{xt} or U_{yt}, so f(n)=f(n+t)
    const UINT uiDir = _DC_Dir;
    const UINT uiN_p_mu = __idx->m_pWalkingTable[uiN * uiDir * 2 + mu + uiDir];
    const SSmallInt4 sN_p_m = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu].m_uiSiteIndex);
    const UBOOL bN_surface = __idx->m_pDeviceIndexPositionToSIndex[1][uiN].IsDirichlet();
    const UBOOL bN_p_mu_surface = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu].IsDirichlet();

    if (0 == i)
    {
        const INT x1 = bN_surface ? 0 : (sN.x - sCenter.x);
        const INT x2 = bN_p_mu_surface ? 0 : (sN_p_m.x - sCenter.x);
        return F(0.5) * static_cast<Real>(x1 * x1 + x2 * x2);
    }

    //else if (0 == i)
    //{
    const INT y1 = bN_surface ? 0 : (sN.y - sCenter.y);
    const INT y2 = bN_p_mu_surface ? 0 : (sN_p_m.y - sCenter.y);
    return F(0.5) * static_cast<Real>(y1 * y1 + y2 * y2);
    //}
}

/**
* Sigma = gi(n) U_nu(n) U_mu(n+nu) U^+_nu(n+mu) + gi(n-nu) U^+_nu(n-nu) U_mu(n-nu) U_nu(n-nu+mu)
*/
static __device__ __inline__ deviceSU3 _deviceStapleTerm4(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, Real fOmegaSq,
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
    deviceSU3 right(
        _deviceGetSTTerm(
            pDeviceData, uiN_m_nu, uiN_m_nu, uiN_p_mu_m_nu, nu, mu, nu, 1, 0, 0
        ));


    //for mu = 0, 1, 2, using _deviceStapleTerm123
    //for mu = 3, always has i == nu
    //assert(i == nu);

    const SSmallInt4 site_N_p_nu = __deviceSiteIndexToInt4(
        __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_nu].m_uiSiteIndex);
    const SSmallInt4 site_N_m_nu = __deviceSiteIndexToInt4(
        __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_nu].m_uiSiteIndex);

    left.MulReal(_deviceGi(sCenter, sSite, site_N_p_nu, uiBigIndex, uiN_p_nu, nu, fOmegaSq));
    right.MulReal(_deviceGi(sCenter, site_N_m_nu, sSite, uiN_m_nu, uiBigIndex, nu, fOmegaSq));
    left.Add(right);

    return left;
}

/**
* Simplified: for mu = 0,1,2, gi(n)=gi(n-nu)
*/
static __device__ __inline__ deviceSU3 _deviceStapleTerm123(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, Real fOmegaSq,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE i)
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
    const deviceSU3 right(
        _deviceGetSTTerm(
            pDeviceData, uiN_m_nu, uiN_m_nu, uiN_p_mu_m_nu, nu, mu, nu, 1, 0, 0
        ));

    if (2 == i)
    {
        //simplified here, for z, site_offset is not needed
        left.Add(right);
        left.MulReal(_deviceGi(sCenter, sSite, sSite, uiBigIndex, uiBigIndex, i, fOmegaSq));
    }
    else
    {
        const UINT uiN_p_i = __idx->m_pWalkingTable[uiBIDir2 + i + uiDir1];
        const SSmallInt4 sSiteN_p_i = __deviceSiteIndexToInt4(
            __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_i].m_uiSiteIndex
        );

        left.Add(right);
        left.MulReal(_deviceGi(sCenter, sSite, sSiteN_p_i, uiBigIndex, uiN_p_i, i, fOmegaSq));
    }

    return left;
}

#pragma endregion

#pragma region Chair terms

/**
* U(N) U(N+rho) U(N+nu) - U(N-rho) U(N-rho) U(N-rho+nu)
* rho nu rho
* + + -, - + +
*/
static __device__ __inline__ deviceSU3 _deviceS1(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;
    const UINT uiBIDir2 = uiBigIndex * uiDir2;

    const UINT uiN_p_rho = __idx->m_pWalkingTable[uiBIDir2 + rho + uiDir1];
    const UINT uiN_p_nu = __idx->m_pWalkingTable[uiBIDir2 + nu + uiDir1];
    const UINT uiN_m_rho = __idx->m_pWalkingTable[uiBIDir2 + rho];
    const UINT uiN_m_rho_p_nu = __idx->m_pWalkingTable[uiN_m_rho * uiDir2 + nu + uiDir1];

    deviceSU3 left(
        _deviceGetSTTerm(
            pDeviceData, uiBigIndex, uiN_p_rho, uiN_p_nu, rho, nu, rho, 0, 0, 1
        ));
    left.Sub(
        _deviceGetSTTerm(
            pDeviceData, uiN_m_rho, uiN_m_rho, uiN_m_rho_p_nu, rho, nu, rho, 1, 0, 0
        ));
    return left;
}

/**
* U(N) U(N-nu+rho) U(N-nu) - U(N-rho) U(N-rho-nu) U(N-rho-nu)
* rho nu rho
* + - -, - - +
*/
static __device__ __inline__ deviceSU3 _deviceS2(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;
    const UINT uiBIDir2 = uiBigIndex * uiDir2;

    const UINT uiN_m_nu = __idx->m_pWalkingTable[uiBIDir2 + nu];
    const UINT uiN_m_nu_p_rho = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + rho + uiDir1];
    const UINT uiN_m_rho = __idx->m_pWalkingTable[uiBIDir2 + rho];
    const UINT uiN_m_rho_m_nu = __idx->m_pWalkingTable[uiN_m_rho * uiDir2 + nu];

    deviceSU3 left(
        _deviceGetSTTerm(
            pDeviceData, uiBigIndex, uiN_m_nu_p_rho, uiN_m_nu, rho, nu, rho, 0, 1, 1
        ));
    left.Sub(
        _deviceGetSTTerm(
            pDeviceData, uiN_m_rho, uiN_m_rho_m_nu, uiN_m_rho_m_nu, rho, nu, rho, 1, 1, 0
        ));
    return left;
}

/**
* U(N+mu-rho+nu) U(N+mu-rho) U(N+mu-rho) - U(N+mu+nu) U(N+mu+rho) U(N+mu)
* rho nu rho
* - - +, + - -
*/
static __device__ __inline__ deviceSU3 _deviceS3(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;

    const UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu + uiDir1];
    const UINT uiN_p_mu_m_rho = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + rho];
    const UINT uiN_p_mu_p_rho = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + rho + uiDir1];
    const UINT uiN_p_mu_p_nu = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + nu + uiDir1];
    const UINT uiN_p_mu_m_rho_p_nu = __idx->m_pWalkingTable[uiN_p_mu_m_rho * uiDir2 + nu + uiDir1];

    deviceSU3 left(
        _deviceGetSTTerm(
            pDeviceData, uiN_p_mu_m_rho_p_nu, uiN_p_mu_m_rho, uiN_p_mu_m_rho, rho, nu, rho, 1, 1, 0
        ));
    left.Sub(
        _deviceGetSTTerm(
            pDeviceData, uiN_p_mu_p_nu, uiN_p_mu_p_rho, uiN_p_mu, rho, nu, rho, 0, 1, 1
        ));

    return left;

}

/**
* U(N+mu-rho-nu) U(N+mu-rho-nu) U(N+mu-rho) - U(N+mu-nu) U(N+mu+rho-nu) U(N+mu)
* rho nu rho
* - + +, + + -
*/
static __device__ __inline__ deviceSU3 _deviceS4(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;

    const UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu + uiDir1];
    const UINT uiN_p_mu_m_rho = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + rho];
    const UINT uiN_p_mu_m_nu = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + nu];

    const UINT uiN_p_mu_m_rho_m_nu = __idx->m_pWalkingTable[uiN_p_mu_m_nu * uiDir2 + rho];
    const UINT uiN_p_mu_p_rho_m_nu = __idx->m_pWalkingTable[uiN_p_mu_m_nu * uiDir2 + rho + uiDir1];

    deviceSU3 left(
        _deviceGetSTTerm(
            pDeviceData, uiN_p_mu_m_rho_m_nu, uiN_p_mu_m_rho_m_nu, uiN_p_mu_m_rho, rho, nu, rho, 1, 0, 0
        ));
    left.Sub(
        _deviceGetSTTerm(
            pDeviceData, uiN_p_mu_m_nu, uiN_p_mu_p_rho_m_nu, uiN_p_mu, rho, nu, rho, 0, 0, 1
        ));

    return left;
}

/**
* U(N+mu-rho) U(N+mu-rho) U(N+mu-rho+nu) - U(N+mu) U(N+mu+rho) U(N+mu+nu)
* rho nu rho
* - + +, + + -
*/
static __device__ __inline__ deviceSU3 _deviceT1(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;

    const UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu + uiDir1];

    const UINT uiNpmuDir2 = uiN_p_mu * uiDir2;
    const UINT uiN_p_mu_m_rho = __idx->m_pWalkingTable[uiNpmuDir2 + rho];
    const UINT uiN_p_mu_p_rho = __idx->m_pWalkingTable[uiNpmuDir2 + rho + uiDir1];
    const UINT uiN_p_mu_p_nu = __idx->m_pWalkingTable[uiNpmuDir2 + nu + uiDir1];
    const UINT uiN_p_mu_m_rho_p_nu = __idx->m_pWalkingTable[uiN_p_mu_m_rho * uiDir2 + nu + uiDir1];

    deviceSU3 left(
        _deviceGetSTTerm(
            pDeviceData, uiN_p_mu_m_rho, uiN_p_mu_m_rho, uiN_p_mu_m_rho_p_nu, rho, nu, rho, 1, 0, 0
        ));
    left.Sub(
        _deviceGetSTTerm(
            pDeviceData, uiN_p_mu, uiN_p_mu_p_rho, uiN_p_mu_p_nu, rho, nu, rho, 0, 0, 1
        ));

    return left;
}

/**
* U(N-mu) U(N-mu+rho) U(N-mu+nu) - U(N-mu-rho) U(N-mu-rho) U(N-mu-rho+nu)
* rho nu rho
* + + -, - + +
*/
static __device__ __inline__ deviceSU3 _deviceT2(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;

    const UINT uiN_m_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu];

    const UINT uiNmmuDir2 = uiN_m_mu * uiDir2;
    const UINT uiN_m_mu_m_rho = __idx->m_pWalkingTable[uiNmmuDir2 + rho];
    const UINT uiN_m_mu_p_rho = __idx->m_pWalkingTable[uiNmmuDir2 + rho + uiDir1];
    const UINT uiN_m_mu_p_nu = __idx->m_pWalkingTable[uiNmmuDir2 + nu + uiDir1];
    const UINT uiN_m_mu_p_nu_m_rho = __idx->m_pWalkingTable[uiN_m_mu_p_nu * uiDir2 + rho];

    deviceSU3 left(
        _deviceGetSTTerm(
            pDeviceData, uiN_m_mu, uiN_m_mu_p_rho, uiN_m_mu_p_nu, rho, nu, rho, 0, 0, 1
        ));
    left.Sub(
        _deviceGetSTTerm(
            pDeviceData, uiN_m_mu_m_rho, uiN_m_mu_m_rho, uiN_m_mu_p_nu_m_rho, rho, nu, rho, 1, 0, 0
        ));
    return left;
}

/**
* i = 0, 1, 2 correspond to x, y and xy
* h_i(N) = x or y or xy
* return h_i(N) + h_i(N + nu), where N is site, and N + nu (or N + mu or ...) is site2
*/
static __device__ __inline__ Real _deviceHi(const SSmallInt4 &center,
    const SSmallInt4 &site, const SSmallInt4 &site2,
    UINT uiSiteBI, UINT uiSite2BI, BYTE i)
{
    if (0 == i)
    {
        const Real fX1 = __idx->m_pDeviceIndexPositionToSIndex[1][uiSiteBI].IsDirichlet() ? F(0.0)
            : static_cast<Real>(site.x - center.x);
        const Real fX2 = __idx->m_pDeviceIndexPositionToSIndex[1][uiSite2BI].IsDirichlet() ? F(0.0)
            : static_cast<Real>(site2.x - center.x);
        return fX1 + fX2;
    }
    else if (1 == i)
    {
        const Real fY1 = __idx->m_pDeviceIndexPositionToSIndex[1][uiSiteBI].IsDirichlet() ? F(0.0)
            : static_cast<Real>(site.y - center.y);
        const Real fY2 = __idx->m_pDeviceIndexPositionToSIndex[1][uiSite2BI].IsDirichlet() ? F(0.0)
            : static_cast<Real>(site2.y - center.y);
        return -fY1 - fY2;
    }
    const Real fX1 = __idx->m_pDeviceIndexPositionToSIndex[1][uiSiteBI].IsDirichlet() ? F(0.0)
        : static_cast<Real>(site.x - center.x);
    const Real fX2 = __idx->m_pDeviceIndexPositionToSIndex[1][uiSite2BI].IsDirichlet() ? F(0.0)
        : static_cast<Real>(site2.x - center.x);
    const Real fY1 = __idx->m_pDeviceIndexPositionToSIndex[1][uiSiteBI].IsDirichlet() ? F(0.0)
        : static_cast<Real>(site.y - center.y);
    const Real fY2 = __idx->m_pDeviceIndexPositionToSIndex[1][uiSite2BI].IsDirichlet() ? F(0.0)
        : static_cast<Real>(site2.y - center.y);
    return fX1 * fY1 + fX2 * fY2;
}

/**
* [hi(n)+hi(n+nu)]S1  U(N+nu) U(N+mu)
* mu nu
* - +,
*/
static __device__ __inline__ deviceSU3 _deviceStapleS1(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
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

    ret.MulReal(_deviceHi(sCenter,
        sSite,
        __deviceSiteIndexToInt4(uiSiteN_p_nu),
        uiBigIndex, uiN_p_nu, i));

    return ret;
}

/**
* [h(N) + h(n-nu)] S2 U(n-nu)U(n+mu-nu)
* mu nu
* + +
*/
static __device__ __inline__ deviceSU3 _deviceStapleS2(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
{
    const UINT uiDir1 = _DC_Dir;
    const UINT uiDir2 = uiDir1 * 2;

    const UINT uiN_m_nu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + nu];
    const UINT uiN_m_nu_p_mu = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + mu + uiDir1];

    const UINT uiSiteN_m_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceS2(pDeviceData, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_nu, mu));
    ret.Mul(_deviceGetGaugeBCSU3Dir(pDeviceData, uiN_m_nu_p_mu, nu));

    ret.MulReal(_deviceHi(sCenter,
        sSite,
        __deviceSiteIndexToInt4(uiSiteN_m_nu),
        uiBigIndex, uiN_m_nu, i));

    return ret;
}

/**
* [h(N+mu) + h(N+mu+nu)]U(n) U(n+nu) S3
* nu mu
* + +
*/
static __device__ __inline__ deviceSU3 _deviceStapleS3(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
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

    ret.MulReal(_deviceHi(sCenter,
        __deviceSiteIndexToInt4(uiSiteN_p_mu),
        __deviceSiteIndexToInt4(uiSiteN_p_mu_p_nu),
        uiN_p_mu, uiN_p_mu_p_nu, i));

    return ret;

}

/**
* [h(N+mu) + h(N+mu-nu)] U(n-nu) U(n-nu) S4
* nu mu
* - +
*/
static __device__ __inline__ deviceSU3 _deviceStapleS4(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
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

    ret.MulReal(_deviceHi(sCenter,
        __deviceSiteIndexToInt4(uiSiteN_p_mu),
        __deviceSiteIndexToInt4(uiSiteN_p_mu_m_nu),
        uiN_p_mu, uiN_p_mu_m_nu, i));

    return ret;
}

/**
* [h(n+mu) + h(n+mu+nu)] U(n) T1 U(n+nu)
* mu mu, + -
*
*/
static __device__ __inline__ deviceSU3 _deviceStapleT1(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
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

    ret.MulReal(_deviceHi(sCenter,
        __deviceSiteIndexToInt4(uiSiteN_p_mu),
        __deviceSiteIndexToInt4(uiSiteN_p_mu_p_nu),
        uiN_p_mu, uiN_p_mu_p_nu, i));

    return ret;
}

/**
* [h(n-mu) + h(n-mu+nu)] U(n-mu) T2 U(n+nu-mu)
* mu mu, - +
*
*/
static __device__ __inline__ deviceSU3 _deviceStapleT2(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
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

    ret.MulReal(_deviceHi(sCenter,
        __deviceSiteIndexToInt4(uiSiteN_m_mu),
        __deviceSiteIndexToInt4(uiSiteN_m_mu_p_nu),
        uiN_m_mu, uiN_m_mu_p_nu, i));

    return ret;
}

/**
* i = 0, 1, 2 for coefficient
* _deviceChairTerm1,2,3 for partial mu, nu, rho
* For partial mu, the staple is (1/8)(s1+s2+s3+s4)
*/
static __device__ __inline__ deviceSU3 _deviceStapleChairTerm1(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
{
    deviceSU3 ret(_deviceStapleS1(pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, i));
    ret.Add(_deviceStapleS2(pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, i));
    ret.Add(_deviceStapleS3(pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, i));
    ret.Add(_deviceStapleS4(pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, i));
    return ret;
}

/**
* i = 0, 1, 2 for coefficient
* _deviceChairTerm1,2,3 for partial mu, nu, rho
* It is (1/8) * (T1+T2 + T1(mu<->rho) + T2(mu<->rho))
*/
static __device__ __inline__ deviceSU3 _deviceStapleChairTerm2(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
{
    deviceSU3 ret(_deviceStapleT1(pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, i));
    ret.Add(_deviceStapleT2(pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, i));
    ret.Add(_deviceStapleT1(pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, rho, nu, mu, i));
    ret.Add(_deviceStapleT2(pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, rho, nu, mu, i));
    return ret;
}

#pragma endregion

#pragma endregion

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CACTIONGAUGEPLAQUETTE_ROTATING_H_

//=============================================================================
// END OF FILE
//=============================================================================