//=============================================================================
// FILENAME : DeviceInlineGaugeChair.h
// 
// DESCRIPTION:
// This should be implemented using inherint machinism, but due to historical reasons, it is now templates
//
//
// REVISION:
//  [07/06/2024 nbale]
//=============================================================================
#ifndef _DEVICEINLINEGAUGECHAIR_H_
#define _DEVICEINLINEGAUGECHAIR_H_

__BEGIN_NAMESPACE

#pragma region Energy

#pragma region Plaqutte term

/**
* Product of 3 terms
* U(uiBIa)_{byDira} . U(uiBIb)_{byDirb} . U(uiBIc)_{byDirc}
* To calcuate staple
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceGetSTTermT(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex& linkA, const SIndex& linkB, const SIndex& linkC)
{
    deviceGauge ret(_deviceGetGaugeBCDirSIndexT(pDeviceData, linkA, byFieldId));
    _mul(ret, _deviceGetGaugeBCDirSIndexT(pDeviceData, linkB, byFieldId));
    _mul(ret, _deviceGetGaugeBCDirSIndexT(pDeviceData, linkC, byFieldId));
    return ret;
}

template<typename deviceGauge>
static __device__ __inline__ Real _device1PlaqutteTermReTrT(
    const deviceGauge* __restrict__ pDeviceData, BYTE byFieldId,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4)
{
    return _retr(_device1PlaqutteTermPPT(pDeviceData, byMu, byNu, uiBigIdx, sSite4, byFieldId));
}

/**
* 3 - 1/4 Retr[ U_{mu,nu}(n)+U_{-mu,nu}(n)+U_{mu,-nu}(n)+U_{-mu,-nu}(n) ]
* = 3 - 1/4 Retr[ U_{mu,nu}(n)+U^+_{mu,nu}(n-mu)+U^+_{mu,nu}(n-nu)+U_{mu,nu}(n-mu-nu) ]
* = 3 - 1/4 Retr[ U_{mu,nu}(n)+U_{mu,nu}(n-mu)+U_{mu,nu}(n-nu)+U_{mu,nu}(n-mu-nu) ]
* Hey! it is wrong but correct!
* In fact, it is
* 3 - 1/4 Retr[ U_{mu,nu}(n)+U^+_{-mu,nu}(n)+U^+_{mu,-nu}(n)+U_{-mu,-nu}(n) ]
* which is
* 3 - 1/4 Retr[ U_{mu,nu}(n)+U_{mu,nu}(n-mu)+U_{mu,nu}(n-nu)+U_{mu,nu}(n-mu-nu) ]
* Thanks to Retr!
*/
template<typename deviceGauge>
static __device__ __inline__ Real _device4PlaqutteTermT(const deviceGauge* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIndex, const SSmallInt4& sSite4, BYTE byFieldId)
{
    return static_cast<Real>(_dim<deviceGauge>()) - F(0.25) * _deviceCloverRetrT(pDeviceData, sSite4, uiBigIndex, byMu, byNu, byFieldId);
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
template<typename deviceGauge>
static __device__ __inline__ Real _deviceChairTermT(const deviceGauge* __restrict__ pDeviceData,
    BYTE byFieldId, const SSmallInt4& sSite,
    BYTE mu, BYTE nu, BYTE rho, UINT uiBigIndex)
{
    const SSmallInt4& n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4& n_m_mu = _deviceSmallInt4OffsetC(sSite, __bck(mu));
    const SSmallInt4& n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));
    const SSmallInt4& n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));
    const SSmallInt4& n_p_rho = _deviceSmallInt4OffsetC(sSite, __fwd(rho));
    const SSmallInt4& n_m_rho = _deviceSmallInt4OffsetC(sSite, __bck(rho));

    const SSmallInt4& n_p_mu_m_nu = _deviceSmallInt4OffsetC(n_p_mu, __bck(nu));
    const SSmallInt4& n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, __fwd(nu));
    const SSmallInt4& n_m_mu_m_nu = _deviceSmallInt4OffsetC(n_m_mu, __bck(nu));
    const SSmallInt4& n_m_rho_p_nu = _deviceSmallInt4OffsetC(n_m_rho, __fwd(nu));
    const SSmallInt4& n_m_rho_m_nu = _deviceSmallInt4OffsetC(n_m_rho, __bck(nu));
    const SSmallInt4& n_m_nu_p_rho = _deviceSmallInt4OffsetC(n_m_nu, __fwd(rho));

    const UINT n_bi4 = uiBigIndex * _DC_Dir;
    const UINT n_p_nu_bi4 = __bi4(n_p_nu);
    const UINT n_m_mu_bi4 = __bi4(n_m_mu);
    const UINT n_m_rho_bi4 = __bi4(n_m_rho);
    const UINT n_m_nu_bi4 = __bi4(n_m_nu);
    const UINT n_m_mu_m_nu_bi4 = __bi4(n_m_mu_m_nu);
    const UINT n_m_rho_m_nu_bi4 = __bi4(n_m_rho_m_nu);

    const SIndex& n__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_bi4 + mu];
    const SIndex& n_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + nu];
    const SIndex& n_p_nu__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_p_nu_bi4 + rho];
    const SIndex& n_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_mu_bi4 + mu];
    const SIndex& n_m_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_mu_bi4 + nu];
    const SIndex& n_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_mu_p_nu) + mu];
    const SIndex& n_m_rho__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_rho_bi4 + rho];
    const SIndex& n_m_mu_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_mu_m_nu_bi4 + mu];
    const SIndex& n_m_nu__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + rho];
    const SIndex& n_m_nu_p_rho__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_nu_p_rho) + nu];
    const SIndex& n_m_rho_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_rho_m_nu_bi4 + nu];

    const SIndex n_m_mu__mu_dag = n_m_mu__mu.DaggerC();

    SIndex n_p_nu__mu_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_p_nu_bi4 + mu];
    SIndex n_p_rho__nu_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_rho) + nu];
    SIndex n__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_bi4 + rho];
    SIndex n_m_rho_p_nu__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_rho_p_nu) + rho];
    SIndex n_m_rho__nu_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_rho_bi4 + nu];
    SIndex n_p_mu_m_nu__nu_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu_m_nu) + nu];
    SIndex n_m_nu__mu_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + mu];
    SIndex n_m_mu_m_nu__nu_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_mu_m_nu_bi4 + nu];
    SIndex n_m_rho_m_nu__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_rho_m_nu_bi4 + rho];

    n_p_nu__mu_dag.m_byTag = n_p_nu__mu_dag.m_byTag ^ _kDaggerOrOpposite;
    n_p_rho__nu_dag.m_byTag = n_p_rho__nu_dag.m_byTag ^ _kDaggerOrOpposite;
    n__rho_dag.m_byTag = n__rho_dag.m_byTag ^ _kDaggerOrOpposite;
    n_m_rho_p_nu__rho_dag.m_byTag = n_m_rho_p_nu__rho_dag.m_byTag ^ _kDaggerOrOpposite;
    n_m_rho__nu_dag.m_byTag = n_m_rho__nu_dag.m_byTag ^ _kDaggerOrOpposite;
    n_p_mu_m_nu__nu_dag.m_byTag = n_p_mu_m_nu__nu_dag.m_byTag ^ _kDaggerOrOpposite;
    n_m_nu__mu_dag.m_byTag = n_m_nu__mu_dag.m_byTag ^ _kDaggerOrOpposite;
    n_m_mu_m_nu__nu_dag.m_byTag = n_m_mu_m_nu__nu_dag.m_byTag ^ _kDaggerOrOpposite;
    n_m_rho_m_nu__rho_dag.m_byTag = n_m_rho_m_nu__rho_dag.m_byTag ^ _kDaggerOrOpposite;

    // U_{mu}(N) U_{nu}(N+mu) U^+_{mu}(n+nu)
    deviceGauge term1(_deviceGetSTTermT(byFieldId, pDeviceData,
        n__mu, n_p_mu__nu, n_p_nu__mu_dag));
    //uiBigIndex, uiN_p_mu, uiN_p_nu, mu, nu, mu, 0, 0, 1));

    //U^+_{mu}(N-mu) U_{nu}(N-mu) U_{mu}(N-mu+nu)
    _sub(term1, _deviceGetSTTermT(byFieldId, pDeviceData,
        n_m_mu__mu_dag, n_m_mu__nu, n_m_mu_p_nu__mu));
    //uiN_m_mu, uiN_m_mu, uiN_m_mu_p_nu, mu, nu, mu, 1, 0, 0));

    // U_{rho}(N+nu) U^+_{nu}(N+rho) U^+_{rho}(N)
    deviceGauge term2(_deviceGetSTTermT(byFieldId, pDeviceData,
        n_p_nu__rho, n_p_rho__nu_dag, n__rho_dag));
    //uiN_p_nu, uiN_p_rho, uiBigIndex, rho, nu, rho, 0, 1, 1));

    // U^+_{rho}(N+nu-rho) U^+_{nu}(N-rho) U_{rho}(N-rho)
    _sub(term2, _deviceGetSTTermT(byFieldId, pDeviceData,
        n_m_rho_p_nu__rho_dag, n_m_rho__nu_dag, n_m_rho__rho));
    //uiN_m_rho_p_nu, uiN_m_rho, uiN_m_rho, rho, nu, rho, 1, 1, 0));

    _mul(term1, term2);

    //pm mu, nu
    //U(mu,-nu) = U(N) U(N+mu-nu) U(N-nu) U(N-nu), 0110
    deviceGauge term3(_deviceGetSTTermT(byFieldId, pDeviceData,
        n__mu, n_p_mu_m_nu__nu_dag, n_m_nu__mu_dag));
    //uiBigIndex, uiN_p_mu_m_nu, uiN_m_nu, mu, nu, mu, 0, 1, 1));

    //mm
    //U(-mu, -nu) = U(N - mu) U(N - mu - nu) U(N - mu - nu) U(N - nu) 1100
    _sub(term3, _deviceGetSTTermT(byFieldId, pDeviceData,
        n_m_mu__mu_dag, n_m_mu_m_nu__nu_dag, n_m_mu_m_nu__mu));
    //uiN_m_mu, uiN_m_mu_m_nu, uiN_m_mu_m_nu, mu, nu, mu, 1, 1, 0));

    //mp, nu, rho
    //mp = U(-mu,nu) = U^+_{mu}(N-mu) U_{nu}(N-mu) U_{mu}(N-mu+nu) U^+_{nu}(N)
    deviceGauge term4(_deviceGetSTTermT(byFieldId, pDeviceData,
        n_m_nu__rho, n_m_nu_p_rho__nu, n__rho_dag));
    //uiN_m_nu, uiN_m_nu_p_rho, uiBigIndex, rho, nu, rho, 0, 0, 1));

    //mm nu rho
    //U(-mu, -nu) = U(N - mu) U(N - mu - nu) U(N - mu - nu) U(N - nu) 1100
    _sub(term4, _deviceGetSTTermT(byFieldId, pDeviceData,
        n_m_rho_m_nu__rho_dag, n_m_rho_m_nu__nu, n_m_rho__rho));
    //uiN_m_rho_m_nu, uiN_m_rho_m_nu, uiN_m_rho, rho, nu, rho, 1, 0, 0));

    _mul(term3, term4);

    _add(term1, term3);

    return _retr(term1);
}

#pragma endregion

#pragma endregion

#pragma region Force

#pragma region Plaqutte term

/**
 * Staple for U_mu, from U_{mu,nu}
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceStapleTermGfactorT(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sSite, Real fOmegaSq,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE i, UBOOL bShifted = FALSE)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));
    SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));
    const SSmallInt4 n_p_mu_m_nu = _deviceSmallInt4OffsetC(n_m_nu, __fwd(mu));

    const SIndex& n__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIndex * _DC_Dir + nu];
    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + mu];
    SIndex n_p_mu__nu_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + nu];
    n_p_mu__nu_dag.m_byTag = n_p_mu__nu_dag.m_byTag ^ _kDaggerOrOpposite;

    const UINT n_m_nu_bi4 = __bi4(n_m_nu);
    SIndex n_m_nu__nu_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + nu];
    n_m_nu__nu_dag.m_byTag = n_m_nu__nu_dag.m_byTag ^ _kDaggerOrOpposite;
    const SIndex& n_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + mu];
    const SIndex& n_p_mu_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu_m_nu) + nu];

    deviceGauge left(
        _deviceGetSTTermT(byFieldId, pDeviceData,
            //pDeviceData, uiBigIndex, uiN_p_nu, uiN_p_mu, nu, mu, nu, 0, 0, 1
            n__nu, n_p_nu__mu, n_p_mu__nu_dag
        ));
    deviceGauge right(
        _deviceGetSTTermT(byFieldId, pDeviceData,
            //pDeviceData, uiN_m_nu, uiN_m_nu, uiN_p_mu_m_nu, nu, mu, nu, 1, 0, 0
            n_m_nu__nu_dag, n_m_nu__mu, n_p_mu_m_nu__nu
        ));

    const Real fLFactor = bShifted
        ? _deviceFiShifted(byFieldId, sSite, i, mu, nu)
        : _deviceFi(byFieldId, sSite, uiBigIndex, i, mu, nu);

    if (bTorus)
    {
        n_m_nu = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(n_m_nu)].m_uiSiteIndex);
    }

    const Real fRFactor = bShifted
        ? _deviceFiShifted(byFieldId, n_m_nu, i, mu, nu)
        : _deviceFi(byFieldId, n_m_nu, __bi(n_m_nu), i, mu, nu);

    _mul(left, fLFactor * fOmegaSq);
    _mul(right, fRFactor * fOmegaSq);
    _add(left, right);

    return left;
}

#pragma endregion

#pragma region Chair terms

/**
* U(N) U(N+rho) U(N+nu) - U(N-rho) U(N-rho) U(N-rho+nu)
* rho nu rho
* + + -, - + +
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceS1T(BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_p_rho = _deviceSmallInt4OffsetC(sSite, __fwd(rho));
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));
    const SSmallInt4 n_m_rho = _deviceSmallInt4OffsetC(sSite, __bck(rho));
    const SSmallInt4 n_m_rho_p_nu = _deviceSmallInt4OffsetC(n_m_rho, __fwd(nu));

    const UINT n_m_rho_bi4 = __bi4(n_m_rho);

    const SIndex& n__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIndex * _DC_Dir + rho];
    const SIndex& n_p_rho__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_rho) + nu];
    const SIndex& n_m_rho__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_rho_bi4 + nu];
    const SIndex& n_m_rho_p_nu__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_rho_p_nu) + rho];

    SIndex n_p_nu__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + rho];
    n_p_nu__rho_dag.m_byTag = n_p_nu__rho_dag.m_byTag ^ _kDaggerOrOpposite;
    SIndex n_m_rho__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_rho_bi4 + rho];
    n_m_rho__rho_dag.m_byTag = n_m_rho__rho_dag.m_byTag ^ _kDaggerOrOpposite;

    deviceGauge left(
        _deviceGetSTTermT(byFieldId, pDeviceData,
            //pDeviceData, uiBigIndex, uiN_p_rho, uiN_p_nu, rho, nu, rho, 0, 0, 1
            n__rho, n_p_rho__nu, n_p_nu__rho_dag
        ));
    _sub(left,
        _deviceGetSTTermT(byFieldId, pDeviceData,
            //            pDeviceData, uiN_m_rho, uiN_m_rho, uiN_m_rho_p_nu, rho, nu, rho, 1, 0, 0
            n_m_rho__rho_dag, n_m_rho__nu, n_m_rho_p_nu__rho
        ));
    return left;
}

/**
* U(N) U(N-nu+rho) U(N-nu) - U(N-rho) U(N-rho-nu) U(N-rho-nu)
* rho nu rho
* + - -, - - +
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceS2T(BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));
    const SSmallInt4 n_m_nu_p_rho = _deviceSmallInt4OffsetC(n_m_nu, __fwd(rho));
    const SSmallInt4 n_m_rho = _deviceSmallInt4OffsetC(sSite, __bck(rho));
    const SSmallInt4 n_m_rho_m_nu = _deviceSmallInt4OffsetC(n_m_rho, __bck(nu));

    const UINT n_m_rho_m_nu_bi4 = __bi4(n_m_rho_m_nu);

    const SIndex& n__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIndex * _DC_Dir + rho];
    const SIndex& n_m_rho_m_nu__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_rho_m_nu_bi4 + rho];

    SIndex n_m_nu_p_rho__nu_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_nu_p_rho) + nu];
    SIndex n_m_nu__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_nu) + rho];
    SIndex n_m_rho__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_rho) + rho];
    SIndex n_m_rho_m_nu__nu_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_rho_m_nu_bi4 + nu];

    n_m_nu_p_rho__nu_dag.m_byTag = n_m_nu_p_rho__nu_dag.m_byTag ^ _kDaggerOrOpposite;
    n_m_nu__rho_dag.m_byTag = n_m_nu__rho_dag.m_byTag ^ _kDaggerOrOpposite;
    n_m_rho__rho_dag.m_byTag = n_m_rho__rho_dag.m_byTag ^ _kDaggerOrOpposite;
    n_m_rho_m_nu__nu_dag.m_byTag = n_m_rho_m_nu__nu_dag.m_byTag ^ _kDaggerOrOpposite;

    deviceGauge left(
        _deviceGetSTTermT(byFieldId, pDeviceData,
            n__rho, n_m_nu_p_rho__nu_dag, n_m_nu__rho_dag
            //pDeviceData, uiBigIndex, uiN_m_nu_p_rho, uiN_m_nu, rho, nu, rho, 0, 1, 1
        ));
    _sub(left,
        _deviceGetSTTermT(byFieldId, pDeviceData,
            //pDeviceData, uiN_m_rho, uiN_m_rho_m_nu, uiN_m_rho_m_nu, rho, nu, rho, 1, 1, 0
            n_m_rho__rho_dag, n_m_rho_m_nu__nu_dag, n_m_rho_m_nu__rho
        ));
    return left;
}

/**
* U(N+mu-rho+nu) U(N+mu-rho) U(N+mu-rho) - U(N+mu+nu) U(N+mu+rho) U(N+mu)
* rho nu rho
* - - +, + - -
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceS3T(BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_mu_m_rho = _deviceSmallInt4OffsetC(n_p_mu, __bck(rho));
    const SSmallInt4 n_p_mu_p_rho = _deviceSmallInt4OffsetC(n_p_mu, __fwd(rho));
    const SSmallInt4 n_p_mu_p_nu = _deviceSmallInt4OffsetC(n_p_mu, __fwd(nu));
    const SSmallInt4 n_p_mu_m_rho_p_nu = _deviceSmallInt4OffsetC(n_p_mu_m_rho, __fwd(nu));

    const UINT n_p_mu_m_rho_bi4 = __bi4(n_p_mu_m_rho);

    const SIndex& n_p_mu_m_rho__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_p_mu_m_rho_bi4 + rho];
    const SIndex& n_p_mu_p_nu__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu_p_nu) + rho];

    SIndex n_p_mu_m_rho_p_nu__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu_m_rho_p_nu) + rho];
    SIndex n_p_mu_m_rho__nu_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_p_mu_m_rho_bi4 + nu];
    SIndex n_p_mu_p_rho__nu_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu_p_rho) + nu];
    SIndex n_p_mu__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + rho];

    n_p_mu_m_rho_p_nu__rho_dag.m_byTag = n_p_mu_m_rho_p_nu__rho_dag.m_byTag ^ _kDaggerOrOpposite;
    n_p_mu_m_rho__nu_dag.m_byTag = n_p_mu_m_rho__nu_dag.m_byTag ^ _kDaggerOrOpposite;
    n_p_mu_p_rho__nu_dag.m_byTag = n_p_mu_p_rho__nu_dag.m_byTag ^ _kDaggerOrOpposite;
    n_p_mu__rho_dag.m_byTag = n_p_mu__rho_dag.m_byTag ^ _kDaggerOrOpposite;

    deviceGauge left(
        _deviceGetSTTermT(byFieldId, pDeviceData,
            //pDeviceData, uiN_p_mu_m_rho_p_nu, uiN_p_mu_m_rho, uiN_p_mu_m_rho, rho, nu, rho, 1, 1, 0
            n_p_mu_m_rho_p_nu__rho_dag, n_p_mu_m_rho__nu_dag, n_p_mu_m_rho__rho
        ));
    _sub(left,
        _deviceGetSTTermT(byFieldId, pDeviceData,
            n_p_mu_p_nu__rho, n_p_mu_p_rho__nu_dag, n_p_mu__rho_dag
            //pDeviceData, uiN_p_mu_p_nu, uiN_p_mu_p_rho, uiN_p_mu, rho, nu, rho, 0, 1, 1
        ));

    return left;

}

/**
* U(N+mu-rho-nu) U(N+mu-rho-nu) U(N+mu-rho) - U(N+mu-nu) U(N+mu+rho-nu) U(N+mu)
* rho nu rho
* - + +, + + -
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceS4T(BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_mu_m_rho = _deviceSmallInt4OffsetC(n_p_mu, __bck(rho));
    const SSmallInt4 n_p_mu_m_nu = _deviceSmallInt4OffsetC(n_p_mu, __bck(nu));
    const SSmallInt4 n_p_mu_m_rho_m_nu = _deviceSmallInt4OffsetC(n_p_mu_m_nu, __bck(rho));
    const SSmallInt4 n_p_mu_p_rho_m_nu = _deviceSmallInt4OffsetC(n_p_mu_m_nu, __fwd(rho));

    const UINT n_p_mu_m_rho_m_nu_bi4 = __bi4(n_p_mu_m_rho_m_nu);

    const SIndex& n_p_mu_m_rho_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_p_mu_m_rho_m_nu_bi4 + nu];
    const SIndex& n_p_mu_m_rho__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu_m_rho) + rho];
    const SIndex& n_p_mu_m_nu__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu_m_nu) + rho];
    const SIndex& n_p_mu_p_rho_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu_p_rho_m_nu) + nu];

    SIndex n_p_mu_m_rho_m_nu__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_p_mu_m_rho_m_nu_bi4 + rho];
    SIndex n_p_mu__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + rho];

    n_p_mu_m_rho_m_nu__rho_dag.m_byTag = n_p_mu_m_rho_m_nu__rho_dag.m_byTag ^ _kDaggerOrOpposite;
    n_p_mu__rho_dag.m_byTag = n_p_mu__rho_dag.m_byTag ^ _kDaggerOrOpposite;

    deviceGauge left(
        _deviceGetSTTermT(byFieldId, pDeviceData,
            n_p_mu_m_rho_m_nu__rho_dag, n_p_mu_m_rho_m_nu__nu, n_p_mu_m_rho__rho
            //pDeviceData, uiN_p_mu_m_rho_m_nu, uiN_p_mu_m_rho_m_nu, uiN_p_mu_m_rho, rho, nu, rho, 1, 0, 0
        ));
    _sub(left,
        _deviceGetSTTermT(byFieldId, pDeviceData,
            n_p_mu_m_nu__rho, n_p_mu_p_rho_m_nu__nu, n_p_mu__rho_dag
            //pDeviceData, uiN_p_mu_m_nu, uiN_p_mu_p_rho_m_nu, uiN_p_mu, rho, nu, rho, 0, 0, 1
        ));

    return left;
}

/**
* U(N+mu-rho) U(N+mu-rho) U(N+mu-rho+nu) - U(N+mu) U(N+mu+rho) U(N+mu+nu)
* rho nu rho
* - + +, + + -
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceT1T(BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_mu_m_rho = _deviceSmallInt4OffsetC(n_p_mu, __bck(rho));
    const SSmallInt4 n_p_mu_p_rho = _deviceSmallInt4OffsetC(n_p_mu, __fwd(rho));
    const SSmallInt4 n_p_mu_p_nu = _deviceSmallInt4OffsetC(n_p_mu, __fwd(nu));
    const SSmallInt4 n_p_mu_m_rho_p_nu = _deviceSmallInt4OffsetC(n_p_mu_m_rho, __fwd(nu));

    const UINT n_p_mu_m_rho_bi4 = __bi4(n_p_mu_m_rho);

    const SIndex& n_p_mu_m_rho__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_p_mu_m_rho_bi4 + nu];
    const SIndex& n_p_mu_m_rho_p_nu__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu_m_rho_p_nu) + rho];
    const SIndex& n_p_mu__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + rho];
    const SIndex& n_p_mu_p_rho__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu_p_rho) + nu];

    SIndex n_p_mu_m_rho__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_p_mu_m_rho_bi4 + rho];
    SIndex n_p_mu_p_nu__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu_p_nu) + rho];

    n_p_mu_m_rho__rho_dag.m_byTag = n_p_mu_m_rho__rho_dag.m_byTag ^ _kDaggerOrOpposite;
    n_p_mu_p_nu__rho_dag.m_byTag = n_p_mu_p_nu__rho_dag.m_byTag ^ _kDaggerOrOpposite;

    deviceGauge left(
        _deviceGetSTTermT(byFieldId, pDeviceData,
            n_p_mu_m_rho__rho_dag, n_p_mu_m_rho__nu, n_p_mu_m_rho_p_nu__rho
            //pDeviceData, uiN_p_mu_m_rho, uiN_p_mu_m_rho, uiN_p_mu_m_rho_p_nu, rho, nu, rho, 1, 0, 0
        ));
    _sub(left,
        _deviceGetSTTermT(byFieldId, pDeviceData,
            n_p_mu__rho, n_p_mu_p_rho__nu, n_p_mu_p_nu__rho_dag
            //pDeviceData, uiN_p_mu, uiN_p_mu_p_rho, uiN_p_mu_p_nu, rho, nu, rho, 0, 0, 1
        ));

    return left;
}

/**
* U(N-mu) U(N-mu+rho) U(N-mu+nu) - U(N-mu-rho) U(N-mu-rho) U(N-mu-rho+nu)
* rho nu rho
* + + -, - + +
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceT2T(BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData, const SSmallInt4& sSite,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite, __bck(mu));
    const SSmallInt4 n_m_mu_m_rho = _deviceSmallInt4OffsetC(n_m_mu, __bck(rho));
    const SSmallInt4 n_m_mu_p_rho = _deviceSmallInt4OffsetC(n_m_mu, __fwd(rho));
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, __fwd(nu));
    const SSmallInt4 n_m_mu_p_nu_m_rho = _deviceSmallInt4OffsetC(n_m_mu_m_rho, __fwd(nu));

    const UINT n_m_mu_m_rho_bi4 = __bi4(n_m_mu_m_rho);

    const SIndex& n_m_mu__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_mu) + rho];
    const SIndex& n_m_mu_p_rho__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_mu_p_rho) + nu];
    const SIndex& n_m_mu_m_rho__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_mu_m_rho_bi4 + nu];
    const SIndex& n_m_mu_p_nu_m_rho__rho = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_mu_p_nu_m_rho) + rho];

    SIndex n_m_mu_p_nu__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_mu_p_nu) + rho];
    SIndex n_m_mu_m_rho__rho_dag = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_mu_m_rho_bi4 + rho];

    n_m_mu_p_nu__rho_dag.m_byTag = n_m_mu_p_nu__rho_dag.m_byTag ^ _kDaggerOrOpposite;
    n_m_mu_m_rho__rho_dag.m_byTag = n_m_mu_m_rho__rho_dag.m_byTag ^ _kDaggerOrOpposite;

    deviceGauge left(
        _deviceGetSTTermT(byFieldId, pDeviceData,
            n_m_mu__rho, n_m_mu_p_rho__nu, n_m_mu_p_nu__rho_dag
            //pDeviceData, uiN_m_mu, uiN_m_mu_p_rho, uiN_m_mu_p_nu, rho, nu, rho, 0, 0, 1
        ));
    _sub(left,
        _deviceGetSTTermT(byFieldId, pDeviceData,
            n_m_mu_m_rho__rho_dag, n_m_mu_m_rho__nu, n_m_mu_p_nu_m_rho__rho
            //pDeviceData, uiN_m_mu_m_rho, uiN_m_mu_m_rho, uiN_m_mu_p_nu_m_rho, rho, nu, rho, 1, 0, 0
        ));
    return left;
}


/**
* [hi(n)+hi(n+nu)]S1  U(N+nu) U(N+mu)
* mu nu
* - +,
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceStapleS1T(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, mu + 1);
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, nu + 1);
    const UINT uiN_p_nu = __bi(n_p_nu);
    const SIndex& n_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + nu];
    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_p_nu * _DC_Dir + mu];
    const SIndex& uiSiteN_p_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_nu];

    deviceGauge ret(_deviceS1T(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    _mul(ret, _deviceGetGaugeBCDirSIndexT(pDeviceData, n_p_nu__mu, byFieldId));
    _muldag(ret, _deviceGetGaugeBCDirSIndexT(pDeviceData, n_p_mu__nu, byFieldId));

    _mul(ret,
        _deviceHi(byFieldId,
            sSite,
            n_p_nu,
            __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex], uiSiteN_p_nu, fpt)
    );

    return ret;
}

/**
* [h(N) + h(n-nu)] S2 U(n-nu)U(n+mu-nu)
* mu nu
* + +
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceStapleS2T(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(n_m_nu, __fwd(mu));
    const UINT uiN_m_nu = __bi(n_m_nu);
    const SIndex& n_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_m_nu * _DC_Dir + mu];
    const SIndex& n_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_nu_p_mu) + nu];

    const SIndex& uiSiteN_m_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_m_nu];

    deviceGauge ret(_deviceS2T(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    _mul(ret, _deviceGetGaugeBCDirSIndexT(pDeviceData, n_m_nu__mu, byFieldId));
    _mul(ret, _deviceGetGaugeBCDirSIndexT(pDeviceData, n_m_nu_p_mu__nu, byFieldId));

    _mul(ret,
        _deviceHi(byFieldId,
            sSite,
            n_m_nu,
            __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex], uiSiteN_m_nu, fpt)
    );

    return ret;
}

/**
* [h(N+mu) + h(N+mu+nu)]U(n) U(n+nu) S3
* nu mu
* + +
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceStapleS3T(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_mu_p_nu = _deviceSmallInt4OffsetC(n_p_mu, __fwd(nu));
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));

    const UINT uiN_p_mu = __bi(n_p_mu);
    const UINT uiN_p_mu_p_nu = __bi(n_p_mu_p_nu);

    const SIndex& n__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIndex * _DC_Dir + nu];
    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + mu];

    const SIndex& uiSiteN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_mu];
    const SIndex& uiSiteN_p_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_mu_p_nu];

    deviceGauge ret(_deviceGetGaugeBCDirSIndexT(pDeviceData, n__nu, byFieldId));
    _mul(ret, _deviceGetGaugeBCDirSIndexT(pDeviceData, n_p_nu__mu, byFieldId));
    _mul(ret, _deviceS3T(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

    _mul(ret,
        _deviceHi(byFieldId,
            n_p_mu,
            n_p_mu_p_nu,
            uiSiteN_p_mu, uiSiteN_p_mu_p_nu, fpt)
    );

    return ret;

}

/**
* [h(N+mu) + h(N+mu-nu)] U(n-nu) U(n-nu) S4
* nu mu
* - +
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceStapleS4T(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_mu_m_nu = _deviceSmallInt4OffsetC(n_p_mu, __bck(nu));
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));

    const UINT uiN_p_mu = __bi(n_p_mu);
    const UINT uiN_p_mu_m_nu = __bi(n_p_mu_m_nu);
    const UINT n_m_nu_bi4 = __bi4(n_m_nu);

    const SIndex& n_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + nu];
    const SIndex& n_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + mu];

    const SIndex& uiSiteN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_mu];
    const SIndex& uiSiteN_p_mu_m_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_mu_m_nu];

    deviceGauge ret(_deviceGetGaugeBCDirSIndexT(pDeviceData, n_m_nu__nu, byFieldId));
    _dagmul(ret, _deviceGetGaugeBCDirSIndexT(pDeviceData, n_m_nu__mu, byFieldId));
    _mul(ret, _deviceS4T(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

    _mul(ret,
        _deviceHi(byFieldId,
            n_p_mu,
            n_p_mu_m_nu,
            uiSiteN_p_mu, uiSiteN_p_mu_m_nu, fpt)
    );

    return ret;
}

/**
* [h(n+mu) + h(n+mu+nu)] U(n) T1 U(n+nu)
* mu mu, + -
*
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceStapleT1T(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));
    const SSmallInt4 n_p_mu_p_nu = _deviceSmallInt4OffsetC(n_p_mu, __fwd(nu));

    const UINT uiN_p_mu = __bi(n_p_mu);
    const UINT uiN_p_mu_p_nu = __bi(n_p_mu_p_nu);

    const SIndex& n__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIndex * _DC_Dir + mu];
    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + mu];

    const SIndex& uiSiteN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_mu];
    const SIndex& uiSiteN_p_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_mu_p_nu];

    deviceGauge ret(_deviceGetGaugeBCDirSIndexT(pDeviceData, n__mu, byFieldId));
    _mul(ret, _deviceT1T(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    _muldag(ret, _deviceGetGaugeBCDirSIndexT(pDeviceData, n_p_nu__mu, byFieldId));

    _mul(ret,
        _deviceHi(byFieldId,
            n_p_mu,
            n_p_mu_p_nu,
            uiSiteN_p_mu, uiSiteN_p_mu_p_nu, fpt)
    );

    return ret;
}

/**
* [h(n-mu) + h(n-mu+nu)] U(n-mu) T2 U(n+nu-mu)
* mu mu, - +
*
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceStapleT2T(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite, __bck(mu));
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, __fwd(nu));

    const UINT uiN_m_mu = __bi(n_m_mu);
    const UINT uiN_m_mu_p_nu = __bi(n_m_mu_p_nu);

    const SIndex& n_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_m_mu * _DC_Dir + mu];
    const SIndex& n_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_m_mu_p_nu * _DC_Dir + mu];

    const SIndex& uiSiteN_m_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_m_mu];
    const SIndex& uiSiteN_m_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_m_mu_p_nu];

    deviceGauge ret(_deviceGetGaugeBCDirSIndexT(pDeviceData, n_m_mu__mu, byFieldId));
    _dagmul(ret, _deviceT2T(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    _mul(ret, _deviceGetGaugeBCDirSIndexT(pDeviceData, n_m_mu_p_nu__mu, byFieldId));

    _mul(ret,
        _deviceHi(byFieldId,
            n_m_mu,
            n_m_mu_p_nu,
            uiSiteN_m_mu, uiSiteN_m_mu_p_nu, fpt)
    );

    return ret;
}

/**
* i = 0, 1, 2 for coefficient
* _deviceChairTerm1,2,3 for partial mu, nu, rho
* For partial mu, the staple is (1/8)(s1+s2+s3+s4)
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceStapleChairTerm1T(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    deviceGauge ret(_deviceStapleS1T(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    _add(ret, _deviceStapleS2T(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    _add(ret, _deviceStapleS3T(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    _add(ret, _deviceStapleS4T(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    return ret;
}

/**
* i = 0, 1, 2 for coefficient
* _deviceChairTerm1,2,3 for partial mu, nu, rho
* It is (1/8) * (T1+T2 + T1(mu<->rho) + T2(mu<->rho))
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceStapleChairTerm2T(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    deviceGauge ret(_deviceStapleT1T(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    _add(ret, _deviceStapleT2T(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    _add(ret, _deviceStapleT1T(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, rho, nu, mu, fpt));
    _add(ret, _deviceStapleT2T(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, rho, nu, mu, fpt));
    return ret;
}

#pragma endregion

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _DEVICEINLINEGAUGECHAIR_H_

//=============================================================================
// END OF FILE
//=============================================================================
