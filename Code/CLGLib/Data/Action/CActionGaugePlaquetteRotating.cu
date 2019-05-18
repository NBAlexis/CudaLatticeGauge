//=============================================================================
// FILENAME : CActionGaugePlaquetteRotating.cu
// 
// DESCRIPTION:
// This is the class for rotating su3
//
// REVISION:
//  [05/07/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotating)

#pragma region device function

#pragma region Energy

#pragma region Plaqutte term

/**
* big index is the index of walking table.
* The plaqutte index may not be cached because n may out of boundary, so we calculate every one
* n, n+mu, n+nu, n
*/
__device__ __inline__ deviceSU3 _device1PlaqutteTermPP(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;
    UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + byMu + uiDir1];
    UINT uiN_p_nu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + byNu + uiDir1];
    SIndex uiSiteIdxN = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    SIndex uiSiteIdxN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu];
    SIndex uiSiteIdxN_p_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_nu];

    if (uiSiteIdxN.IsDirichlet())
    {
        //only term is uiLinkIdx2 and uiLinkIdx3
        if (uiSiteIdxN_p_mu.IsDirichlet())
        {
            return uiSiteIdxN_p_nu.IsDirichlet() ? deviceSU3::makeSU3Id()
                : pDeviceData[_deviceGetLinkIndex(uiSiteIdxN_p_nu.m_uiSiteIndex, byMu)].DaggerC();
        }

        deviceSU3 u(pDeviceData[_deviceGetLinkIndex(uiSiteIdxN_p_mu.m_uiSiteIndex, byNu)]);
        if (!uiSiteIdxN_p_nu.IsDirichlet())
        {
            u.MulDagger(pDeviceData[_deviceGetLinkIndex(uiSiteIdxN_p_nu.m_uiSiteIndex, byMu)]);
        }
        return u;
    }

    deviceSU3 u(pDeviceData[_deviceGetLinkIndex(uiSiteIdxN.m_uiSiteIndex, byMu)]);
    if (!uiSiteIdxN_p_mu.IsDirichlet())
    {
        u.Mul(pDeviceData[_deviceGetLinkIndex(uiSiteIdxN_p_mu.m_uiSiteIndex, byNu)]);
    }
    if (!uiSiteIdxN_p_nu.IsDirichlet())
    {
        u.MulDagger(pDeviceData[_deviceGetLinkIndex(uiSiteIdxN_p_nu.m_uiSiteIndex, byMu)]);
    }
    u.MulDagger(pDeviceData[_deviceGetLinkIndex(uiSiteIdxN.m_uiSiteIndex, byNu)]);

    return u;
}

/**
* U(-mu,nu) = U^+_{mu}(N-mu) U_{nu}(N-mu) U_{mu}(N-mu+nu) U^+_{nu}(N)
*/
__device__ __inline__ deviceSU3 _device1PlaqutteTermMP(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;
    UINT uiN_m_mu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + byMu];

    SIndex siteN = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    SIndex siteN_m_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_mu];
    SIndex siteN_m_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[1]
        [__idx->m_pWalkingTable[uiN_m_mu * uiDir2 + byNu + uiDir1]];

    if (siteN_m_mu.IsDirichlet())
    {
        if (siteN_m_mu_p_nu.IsDirichlet())
        {
            return siteN.IsDirichlet() ? deviceSU3::makeSU3Id()
                : pDeviceData[_deviceGetLinkIndex(siteN.m_uiSiteIndex, byNu)].DaggerC();
        }
        deviceSU3 u(pDeviceData[_deviceGetLinkIndex(siteN_m_mu_p_nu.m_uiSiteIndex, byMu)]);
        if (!siteN.IsDirichlet())
        {
            u.MulDagger(pDeviceData[_deviceGetLinkIndex(siteN.m_uiSiteIndex, byNu)]);
        }
        return u;
    }

    deviceSU3 u(pDeviceData[_deviceGetLinkIndex(siteN_m_mu.m_uiSiteIndex, byMu)]);
    u.DaggerMul(pDeviceData[_deviceGetLinkIndex(siteN_m_mu.m_uiSiteIndex, byNu)]);
    if (!siteN_m_mu_p_nu.IsDirichlet())
    {
        u.Mul(pDeviceData[_deviceGetLinkIndex(siteN_m_mu_p_nu.m_uiSiteIndex, byMu)]);
    }
    if (!siteN.IsDirichlet())
    {
        u.MulDagger(pDeviceData[_deviceGetLinkIndex(siteN.m_uiSiteIndex, byNu)]);
    }

    return u;
}

/**
* U(mu,-nu) = U(N) U^+(N+mu-nu) U^+(N-nu) U(N-nu)
*/
__device__ __inline__ deviceSU3 _device1PlaqutteTermPM(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;
    UINT uiN_m_nu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + byNu];

    SIndex siteN = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    SIndex siteN_m_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_nu];
    SIndex siteN_m_nu_p_mu = __idx->m_pDeviceIndexPositionToSIndex[1]
        [__idx->m_pWalkingTable[uiN_m_nu * uiDir2 + byMu + uiDir1]];

    if (siteN.IsDirichlet())
    {
        if (siteN_m_nu.IsDirichlet())
        {
            return siteN_m_nu_p_mu.IsDirichlet() ? deviceSU3::makeSU3Id()
                : pDeviceData[_deviceGetLinkIndex(siteN_m_nu_p_mu.m_uiSiteIndex, byNu)].DaggerC();
        }

        if (siteN_m_nu_p_mu.IsDirichlet())
        {
            deviceSU3 u(pDeviceData[_deviceGetLinkIndex(siteN_m_nu.m_uiSiteIndex, byMu)]);
            u.DaggerMul(pDeviceData[_deviceGetLinkIndex(siteN_m_nu.m_uiSiteIndex, byNu)]);

            return u;
        }
        else
        {
            //u2^+ u3^+ u4
            //=(u3 u2)^+ u4
            deviceSU3 u(pDeviceData[_deviceGetLinkIndex(siteN_m_nu.m_uiSiteIndex, byMu)]);
            u.Mul(pDeviceData[_deviceGetLinkIndex(siteN_m_nu_p_mu.m_uiSiteIndex, byNu)]);
            u.DaggerMul(pDeviceData[_deviceGetLinkIndex(siteN_m_nu.m_uiSiteIndex, byNu)]);

            return u;
        }
    }

    deviceSU3 u(pDeviceData[_deviceGetLinkIndex(siteN.m_uiSiteIndex, byMu)]);
    if (!siteN_m_nu_p_mu.IsDirichlet())
    {
        u.MulDagger(pDeviceData[_deviceGetLinkIndex(siteN_m_nu_p_mu.m_uiSiteIndex, byNu)]);
    }
    if (!siteN_m_nu.IsDirichlet())
    {
        u.MulDagger(pDeviceData[_deviceGetLinkIndex(siteN_m_nu.m_uiSiteIndex, byMu)]);
        u.Mul(pDeviceData[_deviceGetLinkIndex(siteN_m_nu.m_uiSiteIndex, byNu)]);
    }

    return u;
}

/**
* U(-mu,-nu) = U^+(N-mu) U^+(N-mu-nu) U(N-mu-nu) U(N-nu)
*/
__device__ __inline__ deviceSU3 _device1PlaqutteTermMM(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;
    UINT uiN_m_nu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + byNu];
    UINT uiN_m_mu = __idx->m_pWalkingTable[uiBigIdx * uiDir2 + byMu];
    UINT uiN_m_nu_m_mu = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + byMu];
    SIndex sN_m_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_nu];
    SIndex sN_m_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_mu];
    SIndex sN_m_nu_m_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_nu_m_mu];

    if (sN_m_mu.IsDirichlet())
    {
        if (sN_m_nu_m_mu.IsDirichlet())
        {
            return sN_m_nu.IsDirichlet() ? deviceSU3::makeSU3Id()
                : pDeviceData[_deviceGetLinkIndex(sN_m_nu.m_uiSiteIndex, byNu)];
        }
        else
        {
            //u2^+ u3 u4
            deviceSU3 u(pDeviceData[_deviceGetLinkIndex(sN_m_nu_m_mu.m_uiSiteIndex, byNu)]);
            u.DaggerMul(pDeviceData[_deviceGetLinkIndex(sN_m_nu_m_mu.m_uiSiteIndex, byMu)]);
            if (!sN_m_nu.IsDirichlet())
            {
                u.Mul(pDeviceData[_deviceGetLinkIndex(sN_m_nu.m_uiSiteIndex, byNu)]);
            }
            return u;
        }
    }

    if (sN_m_nu_m_mu.IsDirichlet())
    {
        //u1^+u4
        deviceSU3 u(pDeviceData[_deviceGetLinkIndex(sN_m_mu.m_uiSiteIndex, byMu)]);
        if (!sN_m_nu.IsDirichlet())
        {
            u.DaggerMul(pDeviceData[_deviceGetLinkIndex(sN_m_nu.m_uiSiteIndex, byNu)]);
            return u;
        }
        return u.DaggerC();
    }

    //u1^+ u2^+ u3 u4
    //= (u2 u1)^+ u3 u4
    deviceSU3 u(pDeviceData[_deviceGetLinkIndex(sN_m_nu_m_mu.m_uiSiteIndex, byNu)]);
    u.Mul(pDeviceData[_deviceGetLinkIndex(sN_m_mu.m_uiSiteIndex, byMu)]);
    u.DaggerMul(pDeviceData[_deviceGetLinkIndex(sN_m_nu_m_mu.m_uiSiteIndex, byMu)]);
    if (!sN_m_nu.IsDirichlet())
    {
        u.Mul(pDeviceData[_deviceGetLinkIndex(sN_m_nu.m_uiSiteIndex, byNu)]);
    }

    return u;
}


/**
* Product of 3 terms
*/
__device__ __inline__ deviceSU3 _deviceGetSTTerm(
    const deviceSU3* __restrict__ pDeviceData,
    UINT uiBIa, UINT uiBIb, UINT uiBIc,
    BYTE byDira, BYTE byDirb, BYTE byDirc,
    BYTE byDaggera, BYTE byDaggerb, BYTE byDaggerc)
{
    UINT uiLinkIdx1 = _deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[1][uiBIa].m_uiSiteIndex, byDira);
    UINT uiLinkIdx2 = _deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[1][uiBIb].m_uiSiteIndex, byDirb);
    UINT uiLinkIdx3 = _deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[1][uiBIc].m_uiSiteIndex, byDirc);

    deviceSU3 ret(pDeviceData[uiLinkIdx1]);
    if (1 == byDaggera)
    {
        ret.Dagger();
    }

    if (1 == byDaggerb)
    {
        ret.MulDagger(pDeviceData[uiLinkIdx2]);
    }
    else
    {
        ret.Mul(pDeviceData[uiLinkIdx2]);
    }

    if (1 == byDaggerc)
    {
        ret.MulDagger(pDeviceData[uiLinkIdx3]);
    }
    else
    {
        ret.Mul(pDeviceData[uiLinkIdx3]);
    }

    return ret;
}

__device__ __inline__ Real _device1PlaqutteTermReTr(
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
__device__ __inline__ Real _device4PlaqutteTerm(const deviceSU3* __restrict__ pDeviceData, 
    BYTE byMu, BYTE byNu, UINT uiBigIndex)
{
    UINT uiDir2 = _DC_Dir * 2;
    UINT uiN_m_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + byMu];
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
__device__ __inline__ Real _deviceChairTerm(const deviceSU3* __restrict__ pDeviceData,
    BYTE mu, BYTE nu, BYTE rho, UINT uiBigIndex)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;
    UINT uiBigIdxDir2 = uiBigIndex * uiDir2;
    UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIdxDir2 + mu + uiDir1];
    UINT uiN_m_mu = __idx->m_pWalkingTable[uiBigIdxDir2 + mu];
    UINT uiN_p_nu = __idx->m_pWalkingTable[uiBigIdxDir2 + nu + uiDir1];
    UINT uiN_m_nu = __idx->m_pWalkingTable[uiBigIdxDir2 + nu];
    UINT uiN_p_rho = __idx->m_pWalkingTable[uiBigIdxDir2 + rho + uiDir1];
    UINT uiN_m_rho = __idx->m_pWalkingTable[uiBigIdxDir2 + rho];

    UINT uiN_p_mu_m_nu = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + nu];
    UINT uiN_m_mu_p_nu = __idx->m_pWalkingTable[uiN_m_mu * uiDir2 + nu + uiDir1];
    UINT uiN_m_mu_m_nu = __idx->m_pWalkingTable[uiN_m_mu * uiDir2 + nu];
    UINT uiN_m_rho_p_nu = __idx->m_pWalkingTable[uiN_m_rho * uiDir2 + nu + uiDir1];
    UINT uiN_m_rho_m_nu = __idx->m_pWalkingTable[uiN_m_rho * uiDir2 + nu];
    UINT uiN_m_nu_p_rho = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + rho + uiDir1];

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
__device__ __inline__ Real _deviceGi(
    const SSmallInt4& sCenter, 
    const SSmallInt4& sSite, 
    const SSmallInt4& sSiteOffset,
    BYTE i, 
    Real fOmegaSq)
{
    if (0 == i)
    {
        Real fX = static_cast<Real>(sSite.x - sCenter.x);
        Real fXp1 = static_cast<Real>(sSiteOffset.x - sCenter.x);
        return F(0.5) * fOmegaSq * (fX * fX + fXp1 * fXp1);
    }
    else if (1 == i)
    {
        Real fY = static_cast<Real>(sSite.y - sCenter.y);
        Real fYp1 = static_cast<Real>(sSiteOffset.y - sCenter.y);
        return F(0.5) * fOmegaSq * (fY * fY + fYp1 * fYp1);
    }
    Real fX = static_cast<Real>(sSite.x - sCenter.x);
    Real fY = static_cast<Real>(sSite.y - sCenter.y);
    return fOmegaSq * (fX * fX + fY * fY);
}

/**
* Coefficient = (f(n)+f(n+mu)+f(n+nu)+f(n+mu+nu))/4 
* Simplfy: nu is always t direction, so f(n) = f(n+nu), f(n+mu) = f(n+mu+nu)
* Coefficient = (f(n)+f(n+mu))/2
* For 3 == mu, f(n) = f(n+mu)
* This is also true for Dirichlet boundary condition, only Dirichlet on X-Y direction is assumed
*/
__device__ __inline__ Real _deviceFi(
    const SSmallInt4& sCenter,
    UINT uiN, BYTE i, BYTE mu, BYTE nu)
{
    SSmallInt4 sN = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[1][uiN].m_uiSiteIndex);
    if (2 == mu)
    {
        INT x1 = sN.x - sCenter.x;
        INT y1 = sN.y - sCenter.y;
        return static_cast<Real>(x1 * x1 + y1 * y1);
    }

    UINT uiDir = _DC_Dir;
    UINT uiN_p_mu = __idx->m_pWalkingTable[uiN * uiDir * 2 + mu + uiDir];
    SSmallInt4 sN_p_m = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu].m_uiSiteIndex);

    if (0 == i)
    {
        INT x1 = sN.x - sCenter.x;
        INT x2 = sN_p_m.x - sCenter.x;
        return F(0.5) * static_cast<Real>(x1 * x1 + x2 * x2);
    }

    //else if (0 == i)
    //{
        INT y1 = sN.y - sCenter.y;
        INT y2 = sN_p_m.y - sCenter.y;
        return F(0.5) * static_cast<Real>(y1 * y1 + y2 * y2);
    //}
}

/**
* Sigma = gi(n) U_nu(n) U_mu(n+nu) U^+_nu(n+mu) + gi(n-nu) U^+_nu(n-nu) U_mu(n-nu) U_nu(n-nu+mu)
*/
__device__ __inline__ deviceSU3 _deviceStapleTerm4(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, Real fOmegaSq, 
    UINT uiBigIndex, BYTE mu, BYTE nu)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;
    UINT uiBIDir2 = uiBigIndex * uiDir2;
    UINT uiN_p_mu = __idx->m_pWalkingTable[uiBIDir2 + mu + uiDir1];
    UINT uiN_p_nu = __idx->m_pWalkingTable[uiBIDir2 + nu + uiDir1];
    UINT uiN_m_nu = __idx->m_pWalkingTable[uiBIDir2 + nu];
    UINT uiN_p_mu_m_nu = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + mu + uiDir1];

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

    SSmallInt4 site_N_p_nu = __deviceSiteIndexToInt4(
        __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_nu].m_uiSiteIndex);
    SSmallInt4 site_N_m_nu = __deviceSiteIndexToInt4(
        __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_nu].m_uiSiteIndex);

    left.MulReal(_deviceGi(sCenter, sSite, site_N_p_nu, nu, fOmegaSq));
    right.MulReal(_deviceGi(sCenter, site_N_m_nu, sSite, nu, fOmegaSq));
    left.Add(right);

    return left;
}

/**
* Simplified: for mu = 0,1,2, gi(n)=gi(n-nu)
*/
__device__ __inline__ deviceSU3 _deviceStapleTerm123(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, Real fOmegaSq,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE i)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;
    UINT uiBIDir2 = uiBigIndex * uiDir2;
    UINT uiN_p_mu = __idx->m_pWalkingTable[uiBIDir2 + mu + uiDir1];
    UINT uiN_p_nu = __idx->m_pWalkingTable[uiBIDir2 + nu + uiDir1];
    UINT uiN_m_nu = __idx->m_pWalkingTable[uiBIDir2 + nu];
    UINT uiN_p_mu_m_nu = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + mu + uiDir1];

    deviceSU3 left(
        _deviceGetSTTerm(
            pDeviceData, uiBigIndex, uiN_p_nu, uiN_p_mu, nu, mu, nu, 0, 0, 1
        ));
    deviceSU3 right(
        _deviceGetSTTerm(
            pDeviceData, uiN_m_nu, uiN_m_nu, uiN_p_mu_m_nu, nu, mu, nu, 1, 0, 0
        ));

    if (2 == i)
    {
        //simplified here, for z, site_offset is not needed
        left.Add(right);
        left.MulReal(_deviceGi(sCenter, sSite, sSite, i, fOmegaSq));
    }
    else
    {
        SSmallInt4 sSiteN_p_i = __deviceSiteIndexToInt4(
            __idx->m_pDeviceIndexPositionToSIndex[1][
                __idx->m_pWalkingTable[uiBIDir2 + i + uiDir1]
            ].m_uiSiteIndex
        );

        left.Add(right);
        left.MulReal(_deviceGi(sCenter, sSite, sSiteN_p_i, i, fOmegaSq));
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
__device__ __inline__ deviceSU3 _deviceS1(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;
    UINT uiBIDir2 = uiBigIndex * uiDir2;

    UINT uiN_p_rho = __idx->m_pWalkingTable[uiBIDir2 + rho + uiDir1];
    UINT uiN_p_nu = __idx->m_pWalkingTable[uiBIDir2 + nu + uiDir1];
    UINT uiN_m_rho = __idx->m_pWalkingTable[uiBIDir2 + rho];
    UINT uiN_m_rho_p_nu = __idx->m_pWalkingTable[uiN_m_rho * uiDir2 + nu + uiDir1];

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
__device__ __inline__ deviceSU3 _deviceS2(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;
    UINT uiBIDir2 = uiBigIndex * uiDir2;

    UINT uiN_m_nu = __idx->m_pWalkingTable[uiBIDir2 + nu];
    UINT uiN_m_nu_p_rho = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + rho + uiDir1];
    UINT uiN_m_rho = __idx->m_pWalkingTable[uiBIDir2 + rho];
    UINT uiN_m_rho_m_nu = __idx->m_pWalkingTable[uiN_m_rho * uiDir2 + nu];

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
__device__ __inline__ deviceSU3 _deviceS3(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;

    UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu + uiDir1];
    UINT uiN_p_mu_m_rho = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + rho];
    UINT uiN_p_mu_p_rho = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + rho + uiDir1];
    UINT uiN_p_mu_p_nu = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + nu + uiDir1];
    UINT uiN_p_mu_m_rho_p_nu = __idx->m_pWalkingTable[uiN_p_mu_m_rho * uiDir2 + nu + uiDir1];

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
__device__ __inline__ deviceSU3 _deviceS4(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;

    UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu + uiDir1];
    UINT uiN_p_mu_m_rho = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + rho];
    UINT uiN_p_mu_m_nu = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + nu];

    UINT uiN_p_mu_m_rho_m_nu = __idx->m_pWalkingTable[uiN_p_mu_m_nu * uiDir2 + rho];
    UINT uiN_p_mu_p_rho_m_nu = __idx->m_pWalkingTable[uiN_p_mu_m_nu * uiDir2 + rho + uiDir1];

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
__device__ __inline__ deviceSU3 _deviceT1(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;

    UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu + uiDir1];

    UINT uiNpmuDir2 = uiN_p_mu * uiDir2;
    UINT uiN_p_mu_m_rho = __idx->m_pWalkingTable[uiNpmuDir2 + rho];
    UINT uiN_p_mu_p_rho = __idx->m_pWalkingTable[uiNpmuDir2 + rho + uiDir1];
    UINT uiN_p_mu_p_nu = __idx->m_pWalkingTable[uiNpmuDir2 + nu + uiDir1];
    UINT uiN_p_mu_m_rho_p_nu = __idx->m_pWalkingTable[uiN_p_mu_m_rho * uiDir2 + nu + uiDir1];

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
__device__ __inline__ deviceSU3 _deviceT2(const deviceSU3* __restrict__ pDeviceData,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;

    UINT uiN_m_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu];

    UINT uiNmmuDir2 = uiN_m_mu * uiDir2;
    UINT uiN_m_mu_m_rho = __idx->m_pWalkingTable[uiNmmuDir2 + rho];
    UINT uiN_m_mu_p_rho = __idx->m_pWalkingTable[uiNmmuDir2 + rho + uiDir1];
    UINT uiN_m_mu_p_nu = __idx->m_pWalkingTable[uiNmmuDir2 + nu + uiDir1];
    UINT uiN_m_mu_p_nu_m_rho = __idx->m_pWalkingTable[uiN_m_mu_p_nu * uiDir2  + rho];

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
__device__ __inline__ Real _deviceHi(const SSmallInt4 &center, 
    const SSmallInt4 &site, const SSmallInt4 &site2, BYTE i)
{
    if (0 == i)
    {
        return static_cast<Real>(site.x + site2.x) - F(2.0) * center.x;
    }
    else if (1 == i)
    {
        return F(2.0) * center.y - static_cast<Real>(site.y + site2.y);
    }
    Real fX1 = static_cast<Real>(site.x - center.x);
    Real fX2 = static_cast<Real>(site2.x - center.x);
    Real fY1 = static_cast<Real>(site.y - center.y);
    Real fY2 = static_cast<Real>(site2.y - center.y);
    return fX1 * fY1 + fX2 * fY2;
}

/**
* [hi(n)+hi(n+nu)]S1  U(N+nu) U(N+mu) 
* mu nu
* - +,
*/
__device__ __inline__ deviceSU3 _deviceStapleS1(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;
    UINT uiBIDir2 = uiBigIndex * uiDir2;

    UINT uiN_p_mu = __idx->m_pWalkingTable[uiBIDir2 + mu + uiDir1];
    UINT uiN_p_nu = __idx->m_pWalkingTable[uiBIDir2 + nu + uiDir1];

    UINT uiSiteN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu].m_uiSiteIndex;
    UINT uiSiteN_p_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_nu].m_uiSiteIndex;

    deviceSU3 ret(_deviceS1(pDeviceData, uiBigIndex, mu, nu, rho));
    ret.Mul(pDeviceData[_deviceGetLinkIndex(uiSiteN_p_nu, mu)]);
    ret.MulDagger(pDeviceData[_deviceGetLinkIndex(uiSiteN_p_mu, nu)]);

    ret.MulReal(_deviceHi(sCenter, 
        sSite,
        __deviceSiteIndexToInt4(uiSiteN_p_nu), i));

    return ret;
}

/**
* [h(N) + h(n-nu)] S2 U(n-nu)U(n+mu-nu)
* mu nu
* + +
*/
__device__ __inline__ deviceSU3 _deviceStapleS2(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;

    UINT uiN_m_nu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + nu];
    UINT uiN_m_nu_p_mu = __idx->m_pWalkingTable[uiN_m_nu * uiDir2 + mu + uiDir1];

    UINT uiSiteN_m_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_nu].m_uiSiteIndex;
    UINT uiSiteN_m_nu_p_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_nu_p_mu].m_uiSiteIndex;

    deviceSU3 ret(_deviceS2(pDeviceData, uiBigIndex, mu, nu, rho));
    ret.Mul(pDeviceData[_deviceGetLinkIndex(uiSiteN_m_nu, mu)]);
    ret.Mul(pDeviceData[_deviceGetLinkIndex(uiSiteN_m_nu_p_mu, nu)]);

    ret.MulReal(_deviceHi(sCenter,
        sSite,
        __deviceSiteIndexToInt4(uiSiteN_m_nu), i));

    return ret;
}

/**
* [h(N+mu) + h(N+mu+nu)]U(n) U(n+nu) S3
* nu mu
* + +
*/
__device__ __inline__ deviceSU3 _deviceStapleS3(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;

    UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu + uiDir1];
    UINT uiN_p_mu_p_nu = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + nu + uiDir1];
    UINT uiN_p_nu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + nu + uiDir1];

    UINT uiSiteN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu].m_uiSiteIndex;
    UINT uiSiteN_p_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_nu].m_uiSiteIndex;
    UINT uiSiteN_p_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu_p_nu].m_uiSiteIndex;

    deviceSU3 ret(pDeviceData[_deviceGetLinkIndex(uiSiteIndex, nu)]);
    ret.Mul(pDeviceData[_deviceGetLinkIndex(uiSiteN_p_nu, mu)]);
    ret.Mul(_deviceS3(pDeviceData, uiBigIndex, mu, nu, rho));

    ret.MulReal(_deviceHi(sCenter,
        __deviceSiteIndexToInt4(uiSiteN_p_mu),
        __deviceSiteIndexToInt4(uiSiteN_p_mu_p_nu), i));

    return ret;

}

/**
* [h(N+mu) + h(N+mu-nu)] U(n-nu) U(n-nu) S4
* nu mu
* - +
*/
__device__ __inline__ deviceSU3 _deviceStapleS4(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;

    UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu + uiDir1];
    UINT uiN_p_mu_m_nu = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + nu];
    UINT uiN_m_nu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + nu];

    UINT uiSiteN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu].m_uiSiteIndex;
    UINT uiSiteN_m_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_nu].m_uiSiteIndex;
    UINT uiSiteN_p_mu_m_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu_m_nu].m_uiSiteIndex;

    deviceSU3 ret(pDeviceData[_deviceGetLinkIndex(uiSiteN_m_nu, nu)]);
    ret.DaggerMul(pDeviceData[_deviceGetLinkIndex(uiSiteN_m_nu, mu)]);
    ret.Mul(_deviceS4(pDeviceData, uiBigIndex, mu, nu, rho));

    ret.MulReal(_deviceHi(sCenter,
        __deviceSiteIndexToInt4(uiSiteN_p_mu),
        __deviceSiteIndexToInt4(uiSiteN_p_mu_m_nu), i));

    return ret;
}

/**
* [h(n+mu) + h(n+mu+nu)] U(n) T1 U(n+nu)
* mu mu, + -
*
*/
__device__ __inline__ deviceSU3 _deviceStapleT1(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;

    UINT uiN_p_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu + uiDir1];
    UINT uiN_p_nu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + nu + uiDir1];
    UINT uiN_p_mu_p_nu = __idx->m_pWalkingTable[uiN_p_mu * uiDir2 + nu + uiDir1];

    UINT uiSiteN_p_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu].m_uiSiteIndex;
    UINT uiSiteN_p_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_mu_p_nu].m_uiSiteIndex;
    UINT uiSiteN_p_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_p_nu].m_uiSiteIndex;

    deviceSU3 ret(pDeviceData[_deviceGetLinkIndex(uiSiteIndex, mu)]);
    ret.Mul(_deviceT1(pDeviceData, uiBigIndex, mu, nu, rho));
    ret.MulDagger(pDeviceData[_deviceGetLinkIndex(uiSiteN_p_nu, mu)]);

    ret.MulReal(_deviceHi(sCenter,
        __deviceSiteIndexToInt4(uiSiteN_p_mu),
        __deviceSiteIndexToInt4(uiSiteN_p_mu_p_nu), i));

    return ret;
}

/**
* [h(n-mu) + h(n-mu+nu)] U(n-mu) T2 U(n+nu-mu)
* mu mu, - +
*
*/
__device__ __inline__ deviceSU3 _deviceStapleT2(const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
{
    UINT uiDir1 = _DC_Dir;
    UINT uiDir2 = uiDir1 * 2;

    UINT uiN_m_mu = __idx->m_pWalkingTable[uiBigIndex * uiDir2 + mu];
    UINT uiN_m_mu_p_nu = __idx->m_pWalkingTable[uiN_m_mu * uiDir2 + nu + uiDir1];

    UINT uiSiteN_m_mu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_mu].m_uiSiteIndex;
    UINT uiSiteN_m_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[1][uiN_m_mu_p_nu].m_uiSiteIndex;

    deviceSU3 ret(pDeviceData[_deviceGetLinkIndex(uiSiteN_m_mu, mu)]);
    ret.DaggerMul(_deviceT2(pDeviceData, uiBigIndex, mu, nu, rho));
    ret.Mul(pDeviceData[_deviceGetLinkIndex(uiSiteN_m_mu_p_nu, mu)]);

    ret.MulReal(_deviceHi(sCenter,
        __deviceSiteIndexToInt4(uiSiteN_m_mu),
        __deviceSiteIndexToInt4(uiSiteN_m_mu_p_nu), i));

    return ret;
}

/**
* i = 0, 1, 2 for coefficient
* _deviceChairTerm1,2,3 for partial mu, nu, rho
* For partial mu, the staple is (1/8)(s1+s2+s3+s4)
*/
__device__ __inline__ deviceSU3 _deviceStapleChairTerm1(const deviceSU3* __restrict__ pDeviceData,
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
__device__ __inline__ deviceSU3 _deviceStapleChairTerm2(const deviceSU3* __restrict__ pDeviceData,
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


#pragma region kernels

/**
* This is slower, just for testing
* directly calculate Retr[1 - \hat{U}]
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermSU3_Test(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    Real betaOverN, Real fOmegaSq, 
    Real* results)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    Real fXSq = (sSite4.x - sCenterSite.x);
    fXSq = fXSq * fXSq;
    Real fYSq = (sSite4.y - sCenterSite.y);
    fYSq = fYSq * fYSq;

    //======================================================
    //4-plaqutte terms
    //Omega^2 x^2 Retr[1 - U_1,4]
    Real fU14 = fOmegaSq * fXSq * _device4PlaqutteTerm(pDeviceData, 0, 3, uiBigIdx);

    //Omega^2 y^2 Retr[1 - U_2,4]
    Real fU24 = fOmegaSq * fYSq * _device4PlaqutteTerm(pDeviceData, 1, 3, uiBigIdx);

    //Omega^2 (x^2 + y^2) Retr[1 - U_3,4]
    Real fU34 = fOmegaSq * (fXSq + fYSq) * _device4PlaqutteTerm(pDeviceData, 2, 3, uiBigIdx);

    results[uiSiteIndex] = (fU14 + fU24 + fU34) * betaOverN;
}

/**
* Using plaqutte and (f(n)+f(n+mu)+f(n+nu)+f(n+mu+nu))/4 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermSU3(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
    SSmallInt4 sCenterSite,
    Real betaOverN, Real fOmegaSq,
    Real* results)
{
    intokernalInt4;

    UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    UINT plaqLength = __idx->m_pSmallData[CIndexData::kPlaqLengthIdx];
    UINT plaqCountAll = __idx->m_pSmallData[CIndexData::kPlaqPerSiteIdx] * plaqLength;
    
    //i=0: 12
    //  1: 13
    //  2: 14
    //  3: 23
    //  4: 24
    //  5: 34
    
    Real resThisThread = F(0.0);

    //========================================
    //find plaqutte 1-4
    BYTE idx = 2;
    SIndex first = pCachedPlaqutte[idx * plaqLength + uiSiteIndex * plaqCountAll];
    deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
    if (first.NeedToDagger())
    {
        toAdd.Dagger();
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[idx * plaqLength + j + uiSiteIndex * plaqCountAll];
        deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
        if (first.NeedToDagger())
        {
            toAdd.MulDagger(toMul);
        }
        else
        {
            toAdd.Mul(toMul);
        }
    }
    resThisThread += (F(3.0) - toAdd.ReTr()) * _deviceFi(sCenterSite, uiN, 0, 0, 3);

    //========================================
    //find plaqutte 2-4
    idx = 4;
    first = pCachedPlaqutte[idx * plaqLength + uiSiteIndex * plaqCountAll];
    toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];
    if (first.NeedToDagger())
    {
        toAdd.Dagger();
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[idx * plaqLength + j + uiSiteIndex * plaqCountAll];
        deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
        if (first.NeedToDagger())
        {
            toAdd.MulDagger(toMul);
        }
        else
        {
            toAdd.Mul(toMul);
        }
    }
    resThisThread += (F(3.0) - toAdd.ReTr()) * _deviceFi(sCenterSite, uiN, 1, 1, 3);

    //========================================
    //find plaqutte 3-4
    idx = 5;
    first = pCachedPlaqutte[idx * plaqLength + uiSiteIndex * plaqCountAll];
    toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];
    if (first.NeedToDagger())
    {
        toAdd.Dagger();
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[idx * plaqLength + j + uiSiteIndex * plaqCountAll];
        deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
        if (first.NeedToDagger())
        {
            toAdd.MulDagger(toMul);
        }
        else
        {
            toAdd.Mul(toMul);
        }
    }
    resThisThread += (F(3.0) - toAdd.ReTr()) * _deviceFi(sCenterSite, uiN, 2, 2, 3);

    results[uiSiteIndex] = resThisThread * betaOverN * fOmegaSq;
}


/**
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term12(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    Real betaOverN, Real fOmega,
    Real* results)
{
    intokernalInt4;

    UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    betaOverN = F(0.125) * betaOverN;
    Real fXOmega = (sSite4.x - sCenterSite.x) * fOmega;

    //===============
    //+x Omega V412
    Real fV412 = fXOmega * _deviceChairTerm(pDeviceData, 3, 0, 1, uiN);

    //===============
    //+x Omega V432
    Real fV432 = fXOmega * _deviceChairTerm(pDeviceData, 3, 2, 1, uiN);

    results[uiSiteIndex] = (fV412 + fV432) * betaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term34(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    Real betaOverN, Real fOmega,
    Real* results)
{
    intokernalInt4;

    UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    betaOverN = F(0.125) * betaOverN;
    Real fYOmega = -(sSite4.y - sCenterSite.y) * fOmega;

    //===============
    //-y Omega V421
    Real fV421 = fYOmega * _deviceChairTerm(pDeviceData, 3, 1, 0, uiN);

    //===============
    //-y Omega V431
    Real fV431 = fYOmega * _deviceChairTerm(pDeviceData, 3, 2, 0, uiN);

    results[uiSiteIndex] = (fV421 + fV431) * betaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term5(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    Real betaOverN, Real fOmegaSq,
    Real* results)
{
    intokernalInt4;

    UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    betaOverN = F(0.125) * betaOverN;
    Real fXYOmega2 = (sSite4.x - sCenterSite.x) * (sSite4.y - sCenterSite.y) * fOmegaSq;

    //===============
    //+Omega^2 xy V142
    Real fV142 = fXYOmega2 * _deviceChairTerm(pDeviceData, 0, 3, 1, uiN);

    results[uiSiteIndex] = fV142 * betaOverN;
}

/**
* 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermSU3(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmegaSq)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    UINT uiDir = _DC_Dir;

    betaOverN = betaOverN * F(-0.5);
    deviceSU3 plaqSum = deviceSU3::makeSU3Zero();

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        if (3 != idir)
        {
            //mu = idir, nu = 4, i = mu
            deviceSU3 stap(_deviceStapleTerm123(pDeviceData, sCenterSite, sSite4, fOmegaSq, uiBigIdx, idir, 3, idir));
            deviceSU3 force(pDeviceData[linkIndex]);
            force.MulDagger(stap);
            force.Ta();
            force.MulReal(betaOverN);
            pForceData[linkIndex].Add(force);
        }
        else// if (3 == idir)
        {
            //mu = idir, nu = i = sum _1-3
            deviceSU3 stap(_deviceStapleTerm4(pDeviceData, sCenterSite, sSite4, fOmegaSq, uiBigIdx, idir, 0));
            stap.Add(_deviceStapleTerm4(pDeviceData, sCenterSite, sSite4, fOmegaSq, uiBigIdx, idir, 1));
            stap.Add(_deviceStapleTerm123(pDeviceData, sCenterSite, sSite4, fOmegaSq, uiBigIdx, idir, 2, 2));
            deviceSU3 force(pDeviceData[linkIndex]);
            force.MulDagger(stap);
            force.Ta();
            force.MulReal(betaOverN);
            pForceData[linkIndex].Add(force);
        }
    }
}

/**
* Split to 5 functions to avoid max-regcount
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term1(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+x Omega V412
    //add force for mu=4
    UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);
    UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);
    UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    deviceSU3 staple_term1_4 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        3, 0, 1, 0);

    deviceSU3 staple_term1_2 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        1, 0, 3, 0);

    deviceSU3 staple_term1_1 = _deviceStapleChairTerm2(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        3, 0, 1, 0);

    //=====================================
    //all

    deviceSU3 force1(pDeviceData[uiLink1]);
    force1.MulDagger(staple_term1_1);
    force1.Ta();
    force1.MulReal(betaOverN);
    pForceData[uiLink1].Add(force1);

    deviceSU3 force2(pDeviceData[uiLink2]);
    force2.MulDagger(staple_term1_2);
    force2.Ta();
    force2.MulReal(betaOverN);
    pForceData[uiLink2].Add(force2);

    deviceSU3 force4(pDeviceData[uiLink4]);
    force4.MulDagger(staple_term1_4);
    force4.Ta();
    force4.MulReal(betaOverN);
    pForceData[uiLink4].Add(force4);

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term2(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+x Omega V412
    //add force for mu=4
    UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);
    UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);
    UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    //===============
    //+x Omega V432
    deviceSU3 staple_term2_4 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        3, 2, 1, 0);

    deviceSU3 staple_term2_2 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        1, 2, 3, 0);

    deviceSU3 staple_term2_3 = _deviceStapleChairTerm2(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        3, 2, 1, 0);

 
    //=====================================
    //all

    deviceSU3 force2(pDeviceData[uiLink2]);
    force2.MulDagger(staple_term2_2);
    force2.Ta();
    force2.MulReal(betaOverN);
    pForceData[uiLink2].Add(force2);

    deviceSU3 force3(pDeviceData[uiLink3]);
    force3.MulDagger(staple_term2_3);
    force3.Ta();
    force3.MulReal(betaOverN);
    pForceData[uiLink3].Add(force3);

    deviceSU3 force4(pDeviceData[uiLink4]);
    force4.MulDagger(staple_term2_4);
    force4.Ta();
    force4.MulReal(betaOverN);
    pForceData[uiLink4].Add(force4);

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term3(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+x Omega V412
    //add force for mu=4
    UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);
    UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);
    UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    //===============
    //-y Omega V421
    deviceSU3 staple_term3_4 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        3, 1, 0, 1);

    deviceSU3 staple_term3_1 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        0, 1, 3, 1);

    deviceSU3 staple_term3_2 = _deviceStapleChairTerm2(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        3, 1, 0, 1);


    //=====================================
    //all

    deviceSU3 force1(pDeviceData[uiLink1]);
    force1.MulDagger(staple_term3_1);
    force1.Ta();
    force1.MulReal(betaOverN);
    pForceData[uiLink1].Add(force1);

    deviceSU3 force2(pDeviceData[uiLink2]);
    force2.MulDagger(staple_term3_2);
    force2.Ta();
    force2.MulReal(betaOverN);
    pForceData[uiLink2].Add(force2);

    deviceSU3 force4(pDeviceData[uiLink4]);
    force4.MulDagger(staple_term3_4);
    force4.Ta();
    force4.MulReal(betaOverN);
    pForceData[uiLink4].Add(force4);

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term4(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+x Omega V412
    //add force for mu=4
    UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);
    UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);
    UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

 
    //===============
    //-y Omega V431
    deviceSU3 staple_term4_4 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        3, 2, 0, 1);

    deviceSU3 staple_term4_1 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        0, 2, 3, 1);

    deviceSU3 staple_term4_3 = _deviceStapleChairTerm2(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        3, 2, 0, 1);

    //=====================================
    //all

    deviceSU3 force1(pDeviceData[uiLink1]);
    force1.MulDagger(staple_term4_1);
    force1.Ta();
    force1.MulReal(betaOverN);
    pForceData[uiLink1].Add(force1);

    deviceSU3 force3(pDeviceData[uiLink3]);
    force3.MulDagger(staple_term4_3);
    force3.Ta();
    force3.MulReal(betaOverN);
    pForceData[uiLink3].Add(force3);

    deviceSU3 force4(pDeviceData[uiLink4]);
    force4.MulDagger(staple_term4_4);
    force4.Ta();
    force4.MulReal(betaOverN);
    pForceData[uiLink4].Add(force4);

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term5(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmegaSq)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmegaSq * F(0.125);

    UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);
    UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);
    UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    //===============
    //+Omega^2 xy V142
    deviceSU3 staple_term5_1 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        0, 3, 1, 2);

    deviceSU3 staple_term5_2 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        1, 3, 0, 2);

    deviceSU3 staple_term5_4 = _deviceStapleChairTerm2(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
        0, 3, 1, 2);


    //=====================================
    //all

    deviceSU3 force1(pDeviceData[uiLink1]);
    force1.MulDagger(staple_term5_1);
    force1.Ta();
    force1.MulReal(betaOverN);
    pForceData[uiLink1].Add(force1);

    deviceSU3 force2(pDeviceData[uiLink2]);
    force2.MulDagger(staple_term5_2);
    force2.Ta();
    force2.MulReal(betaOverN);
    pForceData[uiLink2].Add(force2);

    deviceSU3 force4(pDeviceData[uiLink4]);
    force4.MulDagger(staple_term5_4);
    force4.Ta();
    force4.MulReal(betaOverN);
    pForceData[uiLink4].Add(force4);

}

#pragma endregion


CActionGaugePlaquetteRotating::CActionGaugePlaquetteRotating()
    : CAction()
    , m_uiPlaqutteCount(0)
    , m_fLastEnergy(F(0.0))
    , m_fNewEnergy(F(0.0))
    , m_fBetaOverN(F(0.1))
    , m_fOmega(F(0.0))
{
}

void CActionGaugePlaquetteRotating::PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        m_fLastEnergy = Energy(FALSE, pGauge, NULL);
    }
}

void CActionGaugePlaquetteRotating::OnFinishTrajectory(UBOOL bAccepted)
{
    if (bAccepted)
    {
        m_fLastEnergy = m_fNewEnergy;
    }
}

void CActionGaugePlaquetteRotating::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    m_pOwner = pOwner;
    m_byActionId = byId;
    Real fBeta = 0.1f;
    param.FetchValueReal(_T("Beta"), fBeta);
    CCommonData::m_fBeta = fBeta;
    if (NULL != pOwner->m_pGaugeField && EFT_GaugeSU3 == pOwner->m_pGaugeField->GetFieldType())
    {
        fBeta = fBeta / F(3.0);
    }
    m_fBetaOverN = fBeta;
    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

    Real fOmega = 0.1f;
    param.FetchValueReal(_T("Omega"), fOmega);
    m_fOmega = fOmega;

    TArray<INT> centerArray;
    param.FetchValueArrayINT(_T("Center"), centerArray);
    if (centerArray.Num() > 3)
    {
        m_sCenter.x = static_cast<SBYTE>(centerArray[0]);
        m_sCenter.y = static_cast<SBYTE>(centerArray[1]);
        m_sCenter.z = static_cast<SBYTE>(centerArray[2]);
        m_sCenter.w = static_cast<SBYTE>(centerArray[3]);
    }
}

void CActionGaugePlaquetteRotating::SetBeta(Real fBeta)
{
    CCommonData::m_fBeta = fBeta;
    if (NULL != m_pOwner->m_pGaugeField && EFT_GaugeSU3 == m_pOwner->m_pGaugeField->GetFieldType())
    {
        fBeta = fBeta / F(3.0);
    }
    m_fBetaOverN = fBeta;
}

UBOOL CActionGaugePlaquetteRotating::CalculateForceOnGauge(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    pGauge->CalculateForceAndStaple(pForce, pStaple, m_fBetaOverN);

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    if (NULL == pGaugeSU3 || NULL == pForceSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3 now.\n"));
        return TRUE;
    }

    preparethread;

    _kernelAddForce4PlaqutteTermSU3 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

    //_kernelAddForceChairTermSU3 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
    //    pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
Real CActionGaugePlaquetteRotating::Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }
    if (NULL == pStable)
    {
        m_fNewEnergy = pGauge->CalculatePlaqutteEnergy(m_fBetaOverN);
    }
    else
    {
        m_fNewEnergy = pGauge->CalculatePlaqutteEnergyUsingStable(m_fBetaOverN, pStable);
    }

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    if (NULL == pGaugeSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3 now.\n"));
        return m_fNewEnergy;
    }

    preparethread;

    //======== this is only for test ================
    //_kernelAdd4PlaqutteTermSU3_Test << <block, threads >> > (
    //    pGaugeSU3->m_pDeviceData,
    //    m_sCenter,
    //    m_fBetaOverN,
    //    m_fOmega * m_fOmega,
    //    _D_RealThreadBuffer);

    _kernelAdd4PlaqutteTermSU3 << <block, threads >> > (
            pGaugeSU3->m_pDeviceData, 
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
            m_sCenter,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);

    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);


    //_kernelAddChairTermSU3_Term12 << <block, threads >> > (
    //    pGaugeSU3->m_pDeviceData,
    //    m_sCenter,
    //    m_fBetaOverN,
    //    m_fOmega,
    //    _D_RealThreadBuffer);

    //m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    //_kernelAddChairTermSU3_Term34 << <block, threads >> > (
    //    pGaugeSU3->m_pDeviceData,
    //    m_sCenter,
    //    m_fBetaOverN,
    //    m_fOmega,
    //    _D_RealThreadBuffer);

    //m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    //_kernelAddChairTermSU3_Term5 << <block, threads >> > (
    //    pGaugeSU3->m_pDeviceData,
    //    m_sCenter,
    //    m_fBetaOverN,
    //    m_fOmega * m_fOmega,
    //    _D_RealThreadBuffer);

    //m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    return m_fNewEnergy;
}

//Real CActionGaugePlaquetteRotating::GetEnergyPerPlaqutte() const
//{
//    return m_pOwner->m_pGaugeField->CalculatePlaqutteEnergy(m_fBetaOverN) / m_uiPlaqutteCount;
//}

CCString CActionGaugePlaquetteRotating::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CActionGaugePlaquetteRotating\n");
    sRet = sRet + tab + _T("Beta : ") + appFloatToString(CCommonData::m_fBeta) + _T("\n");
    sRet = sRet + tab + _T("Omega : ") + appFloatToString(m_fOmega) + _T("\n");
    CCString sCenter;
    sCenter.Format(_T("Center: [%d, %d, %d, %d]\n")
        , static_cast<INT>(m_sCenter.x)
        , static_cast<INT>(m_sCenter.y)
        , static_cast<INT>(m_sCenter.z)
        , static_cast<INT>(m_sCenter.w));
    sRet = sRet + tab + sCenter;
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================