//=============================================================================
// FILENAME : DeviceInlineGaugeRotationCoefficientFunction.h
// 
// DESCRIPTION:
// This should be implemented using inherint machinism, but due to historical reasons, it is now templates
//
//
// REVISION:
//  [07/03/2024 nbale]
//=============================================================================

#ifndef _DEVICEINLINEGAUGE_NOTEMPLATE_FUNCTION_H_
#define _DEVICEINLINEGAUGE_NOTEMPLATE_FUNCTION_H_

__BEGIN_NAMESPACE

#pragma region Rotation

/**
* g1=O^2(x^2)/2
* g2=O^2(y^2)/2
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
    const SIndex& uiSiteBI,
    const SIndex& uiSiteOffsetBI,
    BYTE i,
    Real fOmegaSq)
{
    if (0 == i)
    {
        const Real fX = uiSiteBI.IsDirichlet() ? F(0.0)
            : static_cast<Real>(sSite.x - sCenter.x);
        return F(0.5) * fOmegaSq * (fX * fX);
        //const Real fXp1 = uiSiteOffsetBI.IsDirichlet() ? F(0.0)
        //    : static_cast<Real>(sSiteOffset.x - sCenter.x);
        //return F(0.5) * fOmegaSq * (fX * fX + fXp1 * fXp1);
    }
    else if (1 == i)
    {
        const Real fY = uiSiteBI.IsDirichlet() ? F(0.0)
            : static_cast<Real>(sSite.y - sCenter.y);
        return F(0.5) * fOmegaSq * (fY * fY);
        //const Real fYp1 = uiSiteOffsetBI.IsDirichlet() ? F(0.0)
        //    : static_cast<Real>(sSiteOffset.y - sCenter.y);
        //return F(0.5) * fOmegaSq * (fY * fY + fYp1 * fYp1);
    }
    const Real fX = uiSiteBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(sSite.x - sCenter.x);
    const Real fY = uiSiteOffsetBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(sSite.y - sCenter.y);
    const Real fXp1 = uiSiteOffsetBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(sSiteOffset.x - sCenter.x);
    const Real fYp1 = uiSiteOffsetBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(sSiteOffset.y - sCenter.y);
    return F(0.5) * fOmegaSq * (fX * fX + fY * fY + fXp1 * fXp1 + fYp1 * fYp1);
}

/**
* Coefficient = (f(n)+f(n+mu)+f(n+nu)+f(n+mu+nu))/4
* Simplfy: nu is always t direction, so f(n) = f(n+nu), f(n+mu) = f(n+mu+nu)
* Coefficient = (f(n)+f(n+mu))/2
* For 3 == mu, f(n) = f(n+mu)
* This is also true for Dirichlet boundary condition, only Dirichlet on X-Y direction is assumed
*
* 0 -> r^2, 1 -> y^2, 2 -> x^2
*
* ==================================================
* Note for periodic boundary condition:
* For const SSmallInt4 sN_p_m = _deviceSmallInt4OffsetC(sSite4, mu + 1)
* sN_p_m.mu can be -1, which leads to a wrong (sN_p_m.y - sCenter.y)
* This '-1' should be set to L_mu - 1. If we consider add the plaquttes as clovers,
* then the coordinates of the centers of the clovers will always be in the lattice,
* so should be set to L_mu - 1
*/
static __device__ __inline__ Real _deviceFi(
    BYTE byFieldId,
    const SSmallInt4& sSite4,
    UINT uiN, BYTE i, BYTE mu, BYTE nu)
{
    //for torus, we need to calculate sSite4 first, because sSite4.x might be -1
    //const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sSite4Start)].m_uiSiteIndex);


    const SSmallInt4 sN_p_mu = _deviceSmallInt4OffsetC(sSite4, mu + 1);
    const SIndex& n_p_mu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_mu)];
    const SSmallInt4 site_N_p_mu = __deviceSiteIndexToInt4(n_p_mu__idx.m_uiSiteIndex);

    const SSmallInt4 sN_p_nu = _deviceSmallInt4OffsetC(sSite4, nu + 1);
    const SIndex& n_p_nu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_nu)];
    const SSmallInt4 site_N_p_nu = __deviceSiteIndexToInt4(n_p_nu__idx.m_uiSiteIndex);

    const SSmallInt4 sN_p_numu = _deviceSmallInt4OffsetC(sN_p_mu, nu + 1);
    const SIndex& n_p_numu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_numu)];
    const SSmallInt4 site_N_p_munu = __deviceSiteIndexToInt4(n_p_numu__idx.m_uiSiteIndex);

    const UBOOL bN_surface = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet();
    const UBOOL bN_p_mu_surface = n_p_mu__idx.IsDirichlet();
    const UBOOL bN_p_nu_surface = n_p_nu__idx.IsDirichlet();
    const UBOOL bN_p_munu_surface = n_p_numu__idx.IsDirichlet();

    const INT x1 = bN_surface ? 0 : (sSite4.x - _DC_Centerx);
    const INT y1 = bN_surface ? 0 : (sSite4.y - _DC_Centery);

    const INT x2 = bN_p_mu_surface ? 0 : (site_N_p_mu.x - _DC_Centerx);
    const INT y2 = bN_p_mu_surface ? 0 : (site_N_p_mu.y - _DC_Centery);

    const INT x3 = bN_p_nu_surface ? 0 : (site_N_p_nu.x - _DC_Centerx);
    const INT y3 = bN_p_nu_surface ? 0 : (site_N_p_nu.y - _DC_Centery);

    const INT x4 = bN_p_munu_surface ? 0 : (site_N_p_munu.x - _DC_Centerx);
    const INT y4 = bN_p_munu_surface ? 0 : (site_N_p_munu.y - _DC_Centery);

    if (0 == i)
    {
        //return F(0.0);
        return F(0.25) * static_cast<Real>(x1 * x1 + y1 * y1
            + x2 * x2 + y2 * y2
            + x3 * x3 + y3 * y3
            + x4 * x4 + y4 * y4);
    }

    if (1 == i)
    {
        return F(0.25) * static_cast<Real>(
            y1 * y1
            + y2 * y2
            + y3 * y3
            + y4 * y4);
    }
    return F(0.25) * static_cast<Real>(
        x1 * x1
        + x2 * x2
        + x3 * x3
        + x4 * x4);
}

static __device__ __inline__ Real _deviceFiShifted(
    BYTE byFieldId,
    const SSmallInt4& sSite4,
    BYTE i, BYTE mu, BYTE nu)
{
    const SIndex& n__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sSite4)];
    const SSmallInt4 site_N = __deviceSiteIndexToInt4(n__idx.m_uiSiteIndex);

    const SSmallInt4 sN_p_mu = _deviceSmallInt4OffsetC(sSite4, mu + 1);
    const SIndex& n_p_mu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_mu)];
    const SSmallInt4 site_N_p_mu = __deviceSiteIndexToInt4(n_p_mu__idx.m_uiSiteIndex);

    const SSmallInt4 sN_p_nu = _deviceSmallInt4OffsetC(sSite4, nu + 1);
    const SIndex& n_p_nu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_nu)];
    const SSmallInt4 site_N_p_nu = __deviceSiteIndexToInt4(n_p_nu__idx.m_uiSiteIndex);

    const SSmallInt4 sN_p_numu = _deviceSmallInt4OffsetC(sN_p_mu, nu + 1);
    const SIndex& n_p_numu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_numu)];
    const SSmallInt4 site_N_p_munu = __deviceSiteIndexToInt4(n_p_numu__idx.m_uiSiteIndex);

    //const UBOOL bN_surface = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet();
    //const UBOOL bN_p_mu_surface = n_p_mu__idx.IsDirichlet();
    //const UBOOL bN_p_nu_surface = n_p_nu__idx.IsDirichlet();
    //const UBOOL bN_p_munu_surface = n_p_numu__idx.IsDirichlet();

    const Real x1 = static_cast<Real>(site_N.x - _DC_Centerx + F(0.5));
    const Real y1 = static_cast<Real>(site_N.y - _DC_Centery + F(0.5));

    const Real x2 = static_cast<Real>(site_N_p_mu.x - _DC_Centerx + F(0.5));
    const Real y2 = static_cast<Real>(site_N_p_mu.y - _DC_Centery + F(0.5));

    const Real x3 = static_cast<Real>(site_N_p_nu.x - _DC_Centerx + F(0.5));
    const Real y3 = static_cast<Real>(site_N_p_nu.y - _DC_Centery + F(0.5));

    const Real x4 = static_cast<Real>(site_N_p_munu.x - _DC_Centerx + F(0.5));
    const Real y4 = static_cast<Real>(site_N_p_munu.y - _DC_Centery + F(0.5));

    if (0 == i)
    {
        //const UBOOL bCorner = (sSite4.x == site_N_p_munu.x) && (sSite4.y == site_N_p_munu.y);
        //if (bCorner)
        //{
        //    return F(0.0);
        //}
        return F(0.25) * (x1 * x1 + y1 * y1
            + x2 * x2 + y2 * y2
            + x3 * x3 + y3 * y3
            + x4 * x4 + y4 * y4);
    }

    if (1 == i)
    {
        return F(0.25) * (y1 * y1
            + y2 * y2
            + y3 * y3
            + y4 * y4);
    }
    return F(0.25) * (x1 * x1
        + x2 * x2
        + x3 * x3
        + x4 * x4);
}


typedef Real(*_deviceCoeffFunctionPointer) (
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI);

/**
* i = 0, 1, 2 correspond to x, y and xy
* h_i(N) = x or y or xy
* return h_i(N) + h_i(N + nu), where N is site, and N + nu (or N + mu or ...) is site2
*/
/*
template<typename deviceGauge>
static __device__ __inline__ Real _deviceHi(
    const SSmallInt4 &center,
    const SSmallInt4 &site, const SSmallInt4 &site2,
    const SIndex& uiSiteBI, const SIndex& uiSite2BI, BYTE i)
{
    if (0 == i)
    {
        const Real fX1 = uiSiteBI.IsDirichlet() ? F(0.0)
            : static_cast<Real>(site.x - center.x);
        const Real fX2 = uiSite2BI.IsDirichlet() ? F(0.0)
            : static_cast<Real>(site2.x - center.x);
        return fX1 + fX2;
    }
    else if (1 == i)
    {
        const Real fY1 = uiSiteBI.IsDirichlet() ? F(0.0)
            : static_cast<Real>(site.y - center.y);
        const Real fY2 = uiSite2BI.IsDirichlet() ? F(0.0)
            : static_cast<Real>(site2.y - center.y);
        return -fY1 - fY2;
    }
    const Real fX1 = uiSiteBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(site.x - center.x);
    const Real fX2 = uiSite2BI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(site2.x - center.x);
    const Real fY1 = uiSiteBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(site.y - center.y);
    const Real fY2 = uiSite2BI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(site2.y - center.y);
    return fX1 * fY1 + fX2 * fY2;
}
*/

static __device__ __inline__ Real _deviceHi(
    BYTE byFieldId,
    const SSmallInt4& site, const SSmallInt4& site2,
    const SIndex& uiSiteBI, const SIndex& uiSite2BI, _deviceCoeffFunctionPointer fpt)
{
    return (*fpt)(byFieldId, site, uiSiteBI) + (*fpt)(byFieldId, site2, uiSite2BI);
}

static __device__ __inline__ Real _deviceHi0(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    return uiSiteBI.IsDirichlet() ? F(0.0) : static_cast<Real>(site.x - _DC_Centerx);
}

static __device__ __inline__ Real _deviceHi1(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    return uiSiteBI.IsDirichlet() ? F(0.0) : static_cast<Real>(_DC_Centery - site.y);
}

static __device__ __inline__ Real _deviceHi2(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    const Real fX1 = uiSiteBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(site.x - _DC_Centerx);
    const Real fY1 = uiSiteBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(site.y - _DC_Centery);
    return fX1 * fY1;
}

static __device__ __inline__ Real _deviceHi0T(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    site = __deviceSiteIndexToInt4(uiSiteBI.m_uiSiteIndex);
    return uiSiteBI.IsDirichlet() ? F(0.0) : static_cast<Real>(site.x - _DC_Centerx);
}

static __device__ __inline__ Real _deviceHi1T(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    site = __deviceSiteIndexToInt4(uiSiteBI.m_uiSiteIndex);
    return uiSiteBI.IsDirichlet() ? F(0.0) : static_cast<Real>(_DC_Centery - site.y);
}

static __device__ __inline__ Real _deviceHi2T(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    site = __deviceSiteIndexToInt4(uiSiteBI.m_uiSiteIndex);
    const Real fX1 = uiSiteBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(site.x - _DC_Centerx);
    const Real fY1 = uiSiteBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(site.y - _DC_Centery);
    return fX1 * fY1;
}

#pragma region Projective plane related

//=============================
//The shifted coord should be conflict with Dirichlet, so we do not consider it
//This is for projective plane
//typedef Real _deviceCoeffPeriodic(
//    BYTE byFieldId,
//    const SSmallInt4& center,
//    SSmallInt4 site);

/*
template<typename deviceGauge>
static __device__ __inline__ Real _deviceSiteCoeff(
    UBOOL bTorus,
    SSmallInt4 sSite4, const SSmallInt4& sCenterSite, BYTE byFieldId, BYTE byType)
{
    if (0 == byType)
    {
        //x
        const UBOOL bOpposite = !bTorus && (sSite4.x >= static_cast<SBYTE>(_DC_Lx) || sSite4.x < 0);
        sSite4 = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sSite4)].m_uiSiteIndex);
        if (bOpposite)
        {
            return -sSite4.x + sCenterSite.x - F(0.5);
        }
        return sSite4.x - sCenterSite.x + F(0.5);
    }
    if (1 == byType)
    {
        //y
        const UBOOL bOpposite = !bTorus && (sSite4.y >= static_cast<SBYTE>(_DC_Ly) || sSite4.y < 0);
        sSite4 = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sSite4)].m_uiSiteIndex);
        if (bOpposite)
        {
            return sSite4.y - sCenterSite.y + F(0.5);
        }
        return -sSite4.y + sCenterSite.y - F(0.5);
    }
    if (3 == byType)
    {
        //There should be NO byType = 3?
        sSite4 = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sSite4)].m_uiSiteIndex);
        return -sSite4.y + sCenterSite.y - F(0.5);
    }

    //byType = 2 and this is XY
    const BYTE bOppositeX = (!bTorus && (sSite4.x >= static_cast<SBYTE>(_DC_Lx) || sSite4.x < 0)) ? 1 : 0;
    const BYTE bOppositeY = (!bTorus && (sSite4.y >= static_cast<SBYTE>(_DC_Ly) || sSite4.y < 0)) ? 1 : 0;
    sSite4 = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sSite4)].m_uiSiteIndex);
    const Real fRet = (sSite4.x - sCenterSite.x + F(0.5)) * (sSite4.y - sCenterSite.y + F(0.5));
    if (0 != (bOppositeX ^ bOppositeY))
    {
        return -fRet;
    }
    return fRet;
}
*/

//template<typename deviceGauge> static __device__ __inline__ Real _deviceHiPeriodic(
//    BYTE byFieldId,
//    const SSmallInt4& center,
//    SSmallInt4 site, SSmallInt4 site2, _deviceCoeffPeriodic fpt)
//{
//    return (*fpt)(byFieldId, center, site) + (*fpt)(byFieldId, center, site2);
//}

static __device__ __inline__ Real _deviceHiShifted0(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    const UBOOL bOpposite = (site.x >= static_cast<SBYTE>(_DC_Lx) || site.x < 0);
    site = __deviceSiteIndexToInt4(uiSiteBI.m_uiSiteIndex);
    if (bOpposite)
    {
        return -site.x + _DC_Centerx - F(0.5);
    }
    return site.x - _DC_Centerx + F(0.5);
}

static __device__ __inline__ Real _deviceHiShifted1(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    const UBOOL bOpposite = (site.y >= static_cast<SBYTE>(_DC_Ly) || site.y < 0);
    site = __deviceSiteIndexToInt4(uiSiteBI.m_uiSiteIndex);
    if (bOpposite)
    {
        return site.y - _DC_Centery + F(0.5);
    }
    return -site.y + _DC_Centery - F(0.5);
}

static __device__ __inline__ Real _deviceHiShifted2(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    const BYTE bOppositeX = (site.x >= static_cast<SBYTE>(_DC_Lx) || site.x < 0) ? 1 : 0;
    const BYTE bOppositeY = (site.y >= static_cast<SBYTE>(_DC_Ly) || site.y < 0) ? 1 : 0;
    site = __deviceSiteIndexToInt4(uiSiteBI.m_uiSiteIndex);
    const Real fRet = (site.x - _DC_Centerx + F(0.5)) * (site.y - _DC_Centery + F(0.5));
    if (0 != (bOppositeX ^ bOppositeY))
    {
        return -fRet;
    }
    return fRet;
}

static __device__ __inline__ Real _deviceHiShiftedT0(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    site = __deviceSiteIndexToInt4(uiSiteBI.m_uiSiteIndex);
    return site.x - _DC_Centerx + F(0.5);
}

static __device__ __inline__ Real _deviceHiShiftedT1(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    site = __deviceSiteIndexToInt4(uiSiteBI.m_uiSiteIndex);
    return -site.y + _DC_Centery - F(0.5);
}

static __device__ __inline__ Real _deviceHiShiftedT2(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    site = __deviceSiteIndexToInt4(uiSiteBI.m_uiSiteIndex);
    return (site.x - _DC_Centerx + F(0.5)) * (site.y - _DC_Centery + F(0.5));
}


#pragma endregion

#pragma endregion

#pragma region Gamma KS

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

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _DEVICEINLINEGAUGE_NOTEMPLATE_FUNCTION_H_

//=============================================================================
// END OF FILE
//=============================================================================
