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

enum { _kLinkMaxLength = 8 };

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
    //Multi-GPU: coordinates entering the physics formula must be GLOBAL
    //(P4-1.1/1.4-R1); the SIndex may also be halo-redirected, which
    //__deviceSiteIndexToInt4 would decode as garbage. Identity on single-GPU.
    //Improve-1 (3.7): 32-bit global coordinate, never narrowed through SCHAR.
    const SInt4 site_N_p_mu = _deviceSIndexToGlobalInt4(n_p_mu__idx);

    const SSmallInt4 sN_p_nu = _deviceSmallInt4OffsetC(sSite4, nu + 1);
    const SIndex& n_p_nu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_nu)];
    const SInt4 site_N_p_nu = _deviceSIndexToGlobalInt4(n_p_nu__idx);

    const SSmallInt4 sN_p_numu = _deviceSmallInt4OffsetC(sN_p_mu, nu + 1);
    const SIndex& n_p_numu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_numu)];
    const SInt4 site_N_p_munu = _deviceSIndexToGlobalInt4(n_p_numu__idx);

    const UBOOL bN_surface = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet();
    const UBOOL bN_p_mu_surface = n_p_mu__idx.IsDirichlet();
    const UBOOL bN_p_nu_surface = n_p_nu__idx.IsDirichlet();
    const UBOOL bN_p_munu_surface = n_p_numu__idx.IsDirichlet();

    //Multi-GPU: base site in global coordinates. On MG the base coordinate may
    //be an out-of-range neighbour (halo-redirected SIndex), so resolve it
    //through the position table like the neighbours above (P4-5.4: decoding a
    //halo slot via __deviceSiteIndexToInt4 would read m_pSiteMappingTable out
    //of bounds). Identity on single-GPU.
    const SIndex& n_base__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN];
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(n_base__idx);

    const INT x1 = bN_surface ? 0 : (sSite4G.x - _DC_Centerx);
    const INT y1 = bN_surface ? 0 : (sSite4G.y - _DC_Centery);

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
    //Multi-GPU: coordinates entering the physics formula must be GLOBAL and the
    //SIndex may be halo-redirected (P4-1.1/1.4-R1). Identity on single-GPU.
    //Improve-1 (3.7): 32-bit global coordinate, never narrowed through SCHAR.
    const SInt4 site_N = _deviceSIndexToGlobalInt4(n__idx);

    const SSmallInt4 sN_p_mu = _deviceSmallInt4OffsetC(sSite4, mu + 1);
    const SIndex& n_p_mu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_mu)];
    const SInt4 site_N_p_mu = _deviceSIndexToGlobalInt4(n_p_mu__idx);

    const SSmallInt4 sN_p_nu = _deviceSmallInt4OffsetC(sSite4, nu + 1);
    const SIndex& n_p_nu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_nu)];
    const SInt4 site_N_p_nu = _deviceSIndexToGlobalInt4(n_p_nu__idx);

    const SSmallInt4 sN_p_numu = _deviceSmallInt4OffsetC(sN_p_mu, nu + 1);
    const SIndex& n_p_numu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_numu)];
    const SInt4 site_N_p_munu = _deviceSIndexToGlobalInt4(n_p_numu__idx);

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

typedef Real(*_deviceCoeffFunctionPointerTwoSites) (
    BYTE byFieldId,
    SSmallInt4 site1,
    SSmallInt4 site2,
    const SIndex& uiSiteBI1,
    const SIndex& uiSiteBI2);

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
    //Multi-GPU: the raw coordinate is local; shift by this rank's offset to get
    //the global coordinate entering the physics formula (identity on
    //single-GPU, offset 0).
    return uiSiteBI.IsDirichlet() ? F(0.0) : static_cast<Real>(site.x + static_cast<INT>(_DC_OffsetX) - _DC_Centerx);
}

static __device__ __inline__ Real _deviceHi1(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    return uiSiteBI.IsDirichlet() ? F(0.0) : static_cast<Real>(_DC_Centery - (site.y + static_cast<INT>(_DC_OffsetY)));
}

static __device__ __inline__ Real _deviceHi2(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    const Real fX1 = uiSiteBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(site.x + static_cast<INT>(_DC_OffsetX) - _DC_Centerx);
    const Real fY1 = uiSiteBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(site.y + static_cast<INT>(_DC_OffsetY) - _DC_Centery);
    return fX1 * fY1;
}

static __device__ __inline__ Real _deviceHi0T(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    //Multi-GPU: resolve via SIndex to the GLOBAL coordinate (halo-aware);
    //identity on single-GPU. Improve-1 (3.7): 32-bit global coordinate with
    //its own name; `site` stays the local storage coordinate (unused here).
    const SInt4 siteG = _deviceSIndexToGlobalInt4(uiSiteBI);
    return uiSiteBI.IsDirichlet() ? F(0.0) : static_cast<Real>(siteG.x - _DC_Centerx);
}

static __device__ __inline__ Real _deviceHi1T(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    const SInt4 siteG = _deviceSIndexToGlobalInt4(uiSiteBI);
    return uiSiteBI.IsDirichlet() ? F(0.0) : static_cast<Real>(_DC_Centery - siteG.y);
}

static __device__ __inline__ Real _deviceHi2T(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    const SInt4 siteG = _deviceSIndexToGlobalInt4(uiSiteBI);
    const Real fX1 = uiSiteBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(siteG.x - _DC_Centerx);
    const Real fY1 = uiSiteBI.IsDirichlet() ? F(0.0)
        : static_cast<Real>(siteG.y - _DC_Centery);
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
        const UBOOL bOpposite = !bTorus && (sSite4.x >= static_cast<SCHAR>(_DC_Lx) || sSite4.x < 0);
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
        const UBOOL bOpposite = !bTorus && (sSite4.y >= static_cast<SCHAR>(_DC_Ly) || sSite4.y < 0);
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
    const BYTE bOppositeX = (!bTorus && (sSite4.x >= static_cast<SCHAR>(_DC_Lx) || sSite4.x < 0)) ? 1 : 0;
    const BYTE bOppositeY = (!bTorus && (sSite4.y >= static_cast<SCHAR>(_DC_Ly) || sSite4.y < 0)) ? 1 : 0;
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
    //Multi-GPU: the raw `site` is a LOCAL offset coordinate, so the
    //out-of-range (opposite) test must run on the GLOBAL raw coordinate
    //(raw local + this rank's offset vs the global extent), and the resolved
    //site must be the halo-aware global coordinate. Both degenerate to the old
    //behaviour on single-GPU (offset 0, GlobalL == L).
    //Improve-1 (3.7): the physics formula reads the 32-bit global coordinate
    //under its own name; `site` stays the local storage coordinate.
    const INT iRawGx = static_cast<INT>(site.x) + static_cast<INT>(_DC_OffsetX);
    const UBOOL bOpposite = (iRawGx >= static_cast<INT>(_DC_GlobalLx) || iRawGx < 0);
    const SInt4 siteG = _deviceSIndexToGlobalInt4(uiSiteBI);
    if (bOpposite)
    {
        return -siteG.x + _DC_Centerx - F(0.5);
    }
    return siteG.x - _DC_Centerx + F(0.5);
}

static __device__ __inline__ Real _deviceHiShifted1(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    const INT iRawGy = static_cast<INT>(site.y) + static_cast<INT>(_DC_OffsetY);
    const UBOOL bOpposite = (iRawGy >= static_cast<INT>(_DC_GlobalLy) || iRawGy < 0);
    const SInt4 siteG = _deviceSIndexToGlobalInt4(uiSiteBI);
    if (bOpposite)
    {
        return siteG.y - _DC_Centery + F(0.5);
    }
    return -siteG.y + _DC_Centery - F(0.5);
}

static __device__ __inline__ Real _deviceHiShifted2(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    const INT iRawGx = static_cast<INT>(site.x) + static_cast<INT>(_DC_OffsetX);
    const INT iRawGy = static_cast<INT>(site.y) + static_cast<INT>(_DC_OffsetY);
    const BYTE bOppositeX = (iRawGx >= static_cast<INT>(_DC_GlobalLx) || iRawGx < 0) ? 1 : 0;
    const BYTE bOppositeY = (iRawGy >= static_cast<INT>(_DC_GlobalLy) || iRawGy < 0) ? 1 : 0;
    const SInt4 siteG = _deviceSIndexToGlobalInt4(uiSiteBI);
    const Real fRet = (siteG.x - _DC_Centerx + F(0.5)) * (siteG.y - _DC_Centery + F(0.5));
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
    //Multi-GPU: resolve via SIndex to the halo-aware GLOBAL coordinate
    //(P4-1.3); identity on single-GPU. Improve-1 (3.7): 32-bit global
    //coordinate; `site` stays the local storage coordinate (unused here).
    const SInt4 siteG = _deviceSIndexToGlobalInt4(uiSiteBI);
    return siteG.x - _DC_Centerx + F(0.5);
}

static __device__ __inline__ Real _deviceHiShiftedT1(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    const SInt4 siteG = _deviceSIndexToGlobalInt4(uiSiteBI);
    return -siteG.y + _DC_Centery - F(0.5);
}

static __device__ __inline__ Real _deviceHiShiftedT2(
    BYTE byFieldId,
    SSmallInt4 site,
    const SIndex& uiSiteBI)
{
    const SInt4 siteG = _deviceSIndexToGlobalInt4(uiSiteBI);
    return (siteG.x - _DC_Centerx + F(0.5)) * (siteG.y - _DC_Centery + F(0.5));
}


#pragma endregion

#pragma endregion

#pragma region device functions Measure Staggered Meson Simple

/**
* For pole mass, there are four patterns:
* 1
* (-1)^n_t
* (-1)^n_i
* (-1)^(n_x+n_y+n_z)
* (-1)^(n_x+n_y), ny+nz, nx+nz
* 
* For screen mass, there are four patterns:
* 1
* (-1)^n_z
* (-1)^n_(i \neq z)
* (-1)^(n_x+n_y+n_t)
* 
* So, the following types are considered:
* 0: 1
* 1: (-1)^nx
* 2: (-1)^ny
* 3: (-1)^nz
* 4: (-1)^nt
* 5: (-1)^(nx+ny)
* 6: (-1)^(nx+nz)
* 7: (-1)^(nx+nt)
* 8: (-1)^(ny+nz)
* 9: (-1)^(ny+nt)
* 10: (-1)^(nz+nt)
* 11: (-1)^(ny+nz+nt)
* 12: (-1)^(nx+nz+nt)
* 13: (-1)^(nx+ny+nt)
* 14: (-1)^(nx+ny+nz)
* 15: (-1)^(nx+ny+nz+nt)
*/
static __device__ __inline__ SCHAR _deviceStaggeredFermionSimplePhase(const SSmallInt4& sSite, BYTE byType)
{
    /*
     * typeid | mask | 
     * -------|------|------------------------
     *   0    | 0x0  | 1
     *   1    | 0x1  | (-1)^nx
     *   2    | 0x2  | (-1)^ny
     *   3    | 0x4  | (-1)^nz
     *   4    | 0x8  | (-1)^nt
     *   5    | 0x3  | (-1)^(nx+ny)
     *   6    | 0x5  | (-1)^(nx+nz)
     *   7    | 0x9  | (-1)^(nx+nt)
     *   8    | 0x6  | (-1)^(ny+nz)
     *   9    | 0xA  | (-1)^(ny+nt)
     *  10    | 0xC  | (-1)^(nz+nt)
     *  11    | 0xE  | (-1)^(ny+nz+nt)
     *  12    | 0xD  | (-1)^(nx+nz+nt)
     *  13    | 0xB  | (-1)^(nx+ny+nt)
     *  14    | 0x7  | (-1)^(nx+ny+nz)
     *  15    | 0xF  | (-1)^(nx+ny+nz+nt)
     */
    static const SCHAR mask_table[16] = {
        static_cast<SCHAR>(0x0), static_cast<SCHAR>(0x1), static_cast<SCHAR>(0x2), static_cast<SCHAR>(0x4), 
        static_cast<SCHAR>(0x8), static_cast<SCHAR>(0x3), static_cast<SCHAR>(0x5), static_cast<SCHAR>(0x9),
        static_cast<SCHAR>(0x6), static_cast<SCHAR>(0xA), static_cast<SCHAR>(0xC), static_cast<SCHAR>(0xE), 
        static_cast<SCHAR>(0xD), static_cast<SCHAR>(0xB), static_cast<SCHAR>(0x7), static_cast<SCHAR>(0xF)
    };

    //Count the '1' in the binary string, and calculate parity
    //Odd '1' is 1, even '1' is zero
    static const SCHAR parity_lookup[16] = {
        static_cast<SCHAR>(0), static_cast<SCHAR>(1), static_cast<SCHAR>(1), static_cast<SCHAR>(0), 
        static_cast<SCHAR>(1), static_cast<SCHAR>(0), static_cast<SCHAR>(0), static_cast<SCHAR>(1),
        static_cast<SCHAR>(1), static_cast<SCHAR>(0), static_cast<SCHAR>(0), static_cast<SCHAR>(1), 
        static_cast<SCHAR>(0), static_cast<SCHAR>(1), static_cast<SCHAR>(1), static_cast<SCHAR>(0)
    };

    const INT parity = (sSite.x & 1) | ((sSite.y & 1) << 1) | ((sSite.z & 1) << 2) | ((sSite.w & 1) << 3);
    const SCHAR mask = mask_table[byType];
    return static_cast<SCHAR>(1 - (parity_lookup[parity & mask] << 1));
}

#pragma endregion

#pragma region Fermion - Boson link


/**
 * Same as CFieldFermionKSSU3R
 * full is a list of path directions with length = iLength
 * it will be divided into two list, where l is full[0, iSep], r is (full[iSep, iLength])^dagger
 * l, r should be allocated on device
 */
static __device__ __inline__ void _deviceSeperate(const SCHAR* __restrict__ full, INT iSep, UINT iLength, SCHAR* l, SCHAR* r, BYTE& LL, BYTE& RL)
{
    LL = static_cast<BYTE>(iSep);
    RL = static_cast<BYTE>(iLength - iSep);

    for (INT i = 0; i < LL; ++i)
    {
        l[i] = -full[iSep - i - 1];
    }

    for (INT i = 0; i < RL; ++i)
    {
        r[i] = full[iSep + i];
    }
}

static __device__ __inline__ void _devicePathDagger(const SCHAR* __restrict__ path, SCHAR* res, UINT iLength)
{
    for (UINT i = 0; i < iLength; ++i)
    {
        res[i] = -path[iLength - i - 1];
    }
}

static void Seperate(SCHAR* full, INT iSep, SCHAR* l, SCHAR* r, BYTE& LL, BYTE& RL)
{
    LL = static_cast<BYTE>(iSep);
    RL = static_cast<BYTE>(3 - iSep);

    for (INT i = 0; i < LL; ++i)
    {
        //trace back
        l[i] = -full[iSep - i - 1];

        //If iSep = 0, This loop will not enter
        //If iSep = 1, This is -full[0]
        //If iSep = 2, This is -full[1], -full[0]
    }

    for (INT i = 0; i < RL; ++i)
    {
        r[i] = full[iSep + i];
    }
}

static TArray<SCHAR> PathDagger(const TArray<SCHAR>& path)
{
    TArray<SCHAR> ret;
    for (INT i = 0; i < path.Num(); ++i)
    {
        ret.AddItem(- path[path.Num() - i - 1]);
    }
    return ret;
}

#pragma endregion


__END_NAMESPACE

#endif //#ifndef _DEVICEINLINEGAUGE_NOTEMPLATE_FUNCTION_H_

//=============================================================================
// END OF FILE
//=============================================================================
