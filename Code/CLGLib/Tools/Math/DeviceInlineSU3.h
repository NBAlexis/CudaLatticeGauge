//=============================================================================
// FILENAME : DeviceInlineSU3.h
// 
// DESCRIPTION:
// This should be removed in the furture
//
//
// REVISION:
//  [07/05/2024 nbale]
//=============================================================================
#include "Data/Field/BoundaryField/CFieldBoundaryOne.h"

#ifndef _DEVICEINLINESU3_H_
#define _DEVICEINLINESU3_H_

__BEGIN_NAMESPACE

#pragma region device functions

/**
* Note: for baked plaqutte index, the bond if is set to SIndex
* If it is a "new SIndex" instead, remember to set the m_byTag
*/
static __device__ __inline__ const deviceSU3& _deviceGetGaugeBCSU3(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pBuffer,
    const SIndex& idx)
{
    return idx.IsDirichlet() ?
        ((CFieldBoundaryGaugeSU3*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

/**
* If the bond is on surface, return the Dirichlet
* else, return the element
*/
static __device__ __inline__ const deviceSU3& _deviceGetGaugeBCSU3Dir(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        ((CFieldBoundaryGaugeSU3*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(site) * _DC_Dir + byDir
        ]
        : pBuffer[_deviceGetLinkIndex(site.m_uiSiteIndex, byDir)];
}

static __device__ __inline__ deviceSU3 _deviceGetGaugeBCSU3DirOne(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        deviceSU3::makeSU3Id()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
}

static __device__ __inline__ deviceSU3 _deviceGetGaugeBCSU3DirZero(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        deviceSU3::makeSU3Zero()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
}


static __device__ __inline__ deviceSU3 _deviceGetGaugeBCSU3DirSIndex(
    const deviceSU3* __restrict__ pBuffer,
    const SIndex& idx,
    BYTE byFieldId)
{
    deviceSU3 ret = idx.IsDirichlet() ?
        ((CFieldBoundaryGaugeSU3*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
    if (idx.NeedToDagger())
    {
        ret.Dagger();
    }
    return ret;
}

static __device__ __inline__ deviceSU3 _deviceGetGaugeBCSU3DirOneSIndex(
    const deviceSU3* __restrict__ pBuffer,
    const SIndex& idx)
{
    if (idx.IsDirichlet())
    {
        return deviceSU3::makeSU3Id();
    }
    if (idx.NeedToDagger())
    {
        return pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)].DaggerC();
    }

    return pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

/**
 * Note that, when get zero instead of one, it is minus not dagger
 */
static __device__ __inline__ deviceSU3 _deviceGetGaugeBCSU3DirZeroSIndex(
    const deviceSU3* __restrict__ pBuffer,
    const SIndex& idx)
{
    if (idx.IsDirichlet())
    {
        return deviceSU3::makeSU3Zero();
    }
    if (idx.NeedToDagger())
    {
        return pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)].MulRealC(F(-1.0));
    }

    return pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

/**
 * calculate D_mu A _nu = Delta _mu + [A_mu, A _nu]
 * Use U now to calculate A pure
 * me will be changed, so, if me is A phys, copy me first
 */
static __device__ __inline__ deviceSU3 _deviceDPureMu(
    const deviceSU3* __restrict__ piA,
    const deviceSU3* __restrict__ piApure,
    const SSmallInt4& sSite4,
    UINT uiBigIdx,
    BYTE byMu,
    BYTE byNu,
    BYTE byFieldId)
{
    //i a D A = (A_nu (n) - A_nu (n-mu)) + iApure _mu A _nu - i A _nu Apure _mu
    const UINT uiSiteBig_m_mu = __idx->_deviceGetBigIndex(
        _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1));

    deviceSU3 res = _deviceGetGaugeBCSU3DirZero(byFieldId, piApure, uiBigIdx, byMu); //Apure _mu
    deviceSU3 res2 = _deviceGetGaugeBCSU3DirZero(byFieldId, piA, uiBigIdx, byNu); //A _nu
    res2.Mul(res); //A _nu Apure _mu
    res.Mul(_deviceGetGaugeBCSU3DirZero(byFieldId, piA, uiBigIdx, byNu)); //Apure _mu A _nu
    res.Sub(res2); //[Apure, A]
    res.Add(_deviceGetGaugeBCSU3DirZero(byFieldId, piA, uiBigIdx, byNu));
    res.Sub(_deviceGetGaugeBCSU3DirZeroSIndex(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_m_mu * _DC_Dir + byNu]));
    return res;
}

/**
 * test using (A(N+mu)-A(N-mu))/2
 */
static __device__ __inline__ deviceSU3 _deviceDPureMu2(
    const deviceSU3* __restrict__ piA,
    const deviceSU3* __restrict__ piApure,
    const SSmallInt4& sSite4,
    UINT uiBigIdx,
    BYTE byMu,
    BYTE byNu,
    BYTE byFieldId)
{
    //i a D A = (A_nu (n+mu) - A_nu (n-mu))/2 + iApure _mu A _nu - i A _nu Apure _mu
    const UINT uiSiteBig_m_mu = __idx->_deviceGetBigIndex(
        _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1));
    const UINT uiSiteBig_p_mu = __idx->_deviceGetBigIndex(
        _deviceSmallInt4OffsetC(sSite4, byMu + 1));

    deviceSU3 res = _deviceGetGaugeBCSU3DirZero(byFieldId, piApure, uiBigIdx, byMu); //Apure _mu
    deviceSU3 res2 = _deviceGetGaugeBCSU3DirZero(byFieldId, piA, uiBigIdx, byNu); //A _nu
    res2.Mul(res); //A _nu Apure _mu
    res.Mul(_deviceGetGaugeBCSU3DirZero(byFieldId, piA, uiBigIdx, byNu)); //Apure _mu A _nu
    res.Sub(res2); //[Apure, A]
    res.Add(_deviceGetGaugeBCSU3DirZeroSIndex(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_p_mu * _DC_Dir + byNu]).MulRealC(F(0.5)));
    res.Sub(_deviceGetGaugeBCSU3DirZeroSIndex(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_m_mu * _DC_Dir + byNu]).MulRealC(F(0.5)));
    return res;
}


static __device__ __inline__ deviceSU3 _devicePlaqutte(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
    UINT uiSiteIndex,
    BYTE plaqIdx, //0-5, as 12, 13, 14, 23, 24, 34
    BYTE plaqLength, //Always 4
    BYTE plaqCountAll //Always 24
)
{
    SIndex first = pCachedPlaqutte[plaqIdx * plaqLength + uiSiteIndex * plaqCountAll];
    deviceSU3 toAdd(_deviceGetGaugeBCSU3(byFieldId, pDeviceData, first));
    if (first.NeedToDagger())
    {
        toAdd.Dagger();
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[plaqIdx * plaqLength + j + uiSiteIndex * plaqCountAll];
        deviceSU3 toMul(_deviceGetGaugeBCSU3(byFieldId, pDeviceData, first));
        if (first.NeedToDagger())
        {
            toAdd.MulDagger(toMul);
        }
        else
        {
            toAdd.Mul(toMul);
        }
    }
    return toAdd;
}

/**
 * pDir[] is dirs of path, the dir is:
 *  x,y,z,t : 1,2,3,4
 *  -x,-y,-z,-t: -1,-2,-3,-4
 *
 * NOTE: This function assumes the boundary is always unity
 */
static __device__ __inline__ deviceSU3 _deviceLink(
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceSU3 sRet = deviceSU3::makeSU3Id();
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet.Dagger();
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet.MulDagger(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    sRet.Mul(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
            }
        }

        if (pDir[i] > 0 && i < (byLength - 1)) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return sRet;
}

/**
 * After every move, it maps to inside the lattice
 * Do NOT use it in projective plane boundary condition
 */
static __device__ __inline__ deviceSU3 _deviceLinkLong(
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceSU3 sRet = deviceSU3::makeSU3Id();
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];
        sStartSite = __deviceSiteIndexToInt4(newLink.m_uiSiteIndex);

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet.Dagger();
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet.MulDagger(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    sRet.Mul(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
            }
        }

        if (pDir[i] > 0 && i < (byLength - 1)) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return sRet;
}

#if discard
/**
 * Get magnetic phase
 * type0:
 * Ay(n) = q B nx
 * if twisted Ax(Lx) = - q B ny
 *
 * type1:
 * Ax(n) = - q B ny
 * if twisted Ay(Ly) = q B nx
 *
 * type2:
 * Ax(n) = - q B ny/2
 * Ay(n) = qB nx / 2
 *
 * if twisted
 * Ax(Lx) = - q B ny
 * Ay(Ly) = q B nx
 *
 */
static __device__ __inline__ Real __deviceGetMagneticPhase(const SSmallInt4& sSite, const SSmallInt4& sCenter, BYTE byLink, BYTE byGaugeType, UBOOL bTwisted)
{
    if (0 != byLink && 1 != byLink)
    {
        return F(0.0);
    }

    if (0 == byGaugeType)
    {
        if (1 == byLink)
        {
            return static_cast<Real>(sSite.x - sCenter.x + F(0.5));
        }
        if (bTwisted && sSite.x == _DC_Lxi - 1)
        {
            return -static_cast<Real>(sSite.y - sCenter.y + F(0.5));
        }
        return F(0.0);
    }

    if (1 == byGaugeType)
    {
        if (0 == byLink)
        {
            return -static_cast<Real>(sSite.y - sCenter.y + F(0.5));
        }
        if (bTwisted && sSite.y == _DC_Lyi - 1)
        {
            return static_cast<Real>(sSite.x - sCenter.x + F(0.5));
        }
        return F(0.0);
    }

    if (0 == byLink)
    {
        if (bTwisted && sSite.x == _DC_Lxi - 1)
        {
            return -static_cast<Real>(sSite.y - sCenter.y + F(0.5));
        }
        return -static_cast<Real>(sSite.y - sCenter.y + F(0.5)) * F(0.5);
    }

    if (bTwisted && sSite.x == _DC_Lxi - 1)
    {
        return static_cast<Real>(sSite.x - sCenter.x + F(0.5));
    }
    return static_cast<Real>(sSite.x - sCenter.x + F(0.5)) * F(0.5);
}

/**
 * Same as _deviceLink, but with fixed electric-magnetic field
 * The electric-magnetic field is parameterized as Bz
 * For magnetic field between n1 and n2 is:
 * Ax(Lx) = - q B ny
 * Ay(n) = q B nx
 *
 * If no twisted boundary, then only Ay(n) = q B nx, others are zero
 */
static __device__ __inline__ deviceSU3 _deviceLinkMP(
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, const SSmallInt4& sCenterSite,
    BYTE byLength, BYTE byFieldId, Real fQBz, BYTE byGaugeType, UBOOL bTwistedBoundary,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceSU3 sRet = deviceSU3::makeSU3Id();
    Real fPhase = F(0.0);

    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move-back
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const UINT uiBi4StartSite = __bi4(sStartSite);
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBi4StartSite + byDir];
        const UBOOL bDaggerFinal =
            (newLink.NeedToDagger() && !bDagger)
            || (!newLink.NeedToDagger() && bDagger);
        if (!newLink.IsDirichlet())
        {
            if (0 == newLink.m_byDir || 1 == newLink.m_byDir)
            {
                const SSmallInt4 newSite = __deviceSiteIndexToInt4(newLink.m_uiSiteIndex);
                if (bDaggerFinal)
                {
                    fPhase = fPhase - __deviceGetMagneticPhase(newSite, sCenterSite, newLink.m_byDir, byGaugeType, bTwistedBoundary);
                }
                else
                {
                    fPhase = fPhase + __deviceGetMagneticPhase(newSite, sCenterSite, newLink.m_byDir, byGaugeType, bTwistedBoundary);
                }
            }
        }

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if (bDaggerFinal)
                {
                    sRet.Dagger();
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                if (bDaggerFinal)
                {
                    sRet.MulDagger(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    sRet.Mul(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
            }
        }

        if (pDir[i] > 0 && i < (byLength - 1)) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    fPhase = fPhase * fQBz;
    const CLGComplex cmpPhase = _make_cuComplex(_cos(fPhase), _sin(fPhase));
    sRet.MulComp(cmpPhase);
    return sRet;
}

/**
 * Phase of the link with:
 * Ax(Lx) = - q B ny
 * Ay(n) = q B nx
 */
static __device__ __inline__ Real _devicePhaseM(const SSmallInt4& sPos, BYTE byDir, const SSmallInt4& sCenterSite, Real fQBz)
{
    if (1 == byDir) //y-dir
    {
        return static_cast<Real>(sPos.x - sCenterSite.x + F(0.5)) * fQBz;
    }
    if (0 == byDir) //x-dir
    {
        if (sPos.x == static_cast<INT>(_DC_Lx) - 1)
        {
            return static_cast<Real>(sPos.y - sCenterSite.y + F(0.5)) * (-fQBz);
        }
    }

    return F(0.0);
}

/**
 * This function assumes _DC_Dir = 4
 */
static __device__ __inline__ Real _devicePhaseM(UINT uiLinkIndex, const SSmallInt4& sCenterSite, Real fQBz)
{
    return _devicePhaseM(__deviceLinkIndexToInt4(uiLinkIndex), static_cast<BYTE>(uiLinkIndex & 3), sCenterSite, fQBz);
}
#endif

/**
* big index is the index of walking table.
* The plaqutte index may not be cached because n may out of boundary, so we calculate every one
* n, n+mu, n+nu, n
*
*   <----- ^
*   |      |
*   |      |
*   V      |
* O ------->
*/
static __device__ __inline__ deviceSU3 _device1PlaqutteTermPP(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    //For any boundary condition it is always, site->mu, site_p_mu->nu, site_p_nu->mu+, site->nu+
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite4, byMu + 1);
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite4, byNu + 1);
    const UINT uiB4 = uiBigIdx * _DC_Dir;
    const SIndex& s_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiB4 + byMu];
    const SIndex& s_p_mu_nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_p_mu) * _DC_Dir + byNu];
    const SIndex& s_p_nu_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_p_nu) * _DC_Dir + byMu];
    const SIndex& s_nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiB4 + byNu];

    deviceSU3 u(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_mu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_p_mu_nu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_p_nu_mu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_nu, byFieldId));

    return u;
}

/**
* U(-mu,nu) = U^+_{mu}(N-mu) U_{nu}(N-mu) U_{mu}(N-mu+nu) U^+_{nu}(N)
*
*    ------->
*    ^      |
*    |      |
*    |      V
*    <------- O
*/
static __device__ __inline__ deviceSU3 _device1PlaqutteTermMP(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1);
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, byNu + 1);
    const UINT uin_m_mub4 = __idx->_deviceGetBigIndex(n_m_mu) * _DC_Dir;
    const SIndex& s_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byMu];
    const SIndex& s_m_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byNu];
    const SIndex& s_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_mu_p_nu) * _DC_Dir + byMu];
    const SIndex& s__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byNu];

    deviceSU3 u(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_mu__mu, byFieldId));
    u.DaggerMul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_mu__nu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_mu_p_nu__mu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s__nu, byFieldId));

    return u;
}

/**
* U(mu,-nu) = U(N) U^+(N+mu-nu) U^+(N-nu) U(N-nu)
*
* O  ------->
*    ^      |
*    |      |
*    |      V
*    <-------
*/
static __device__ __inline__ deviceSU3 _device1PlaqutteTermPM(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byNu) - 1);
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(n_m_nu, byMu + 1);
    const UINT uin_m_nub4 = __idx->_deviceGetBigIndex(n_m_nu) * _DC_Dir;
    const SIndex& s_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byMu];
    const SIndex& s_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byNu];
    const SIndex& s_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_nu_p_mu) * _DC_Dir + byNu];
    const SIndex& s__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byMu];

    deviceSU3 u(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s__mu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_nu_p_mu__nu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_nu__mu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_nu__nu, byFieldId));

    return u;
}

/**
* U(-mu,-nu) = U^+(N-mu) U^+(N-mu-nu) U(N-mu-nu) U(N-nu)
*
* <----- ^ O
* |      |
* |      |
* V      |
* ------->
*/
static __device__ __inline__ deviceSU3 _device1PlaqutteTermMM(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1);
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byNu) - 1);
    const SSmallInt4 n_m_nu_m_mu = _deviceSmallInt4OffsetC(n_m_nu, -static_cast<INT>(byMu) - 1);
    const UINT uin_m_nu_m_mub4 = __idx->_deviceGetBigIndex(n_m_nu_m_mu) * _DC_Dir;

    const SIndex& s_m_nu_m_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nu_m_mub4 + byNu];
    const SIndex& s_m_nu_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nu_m_mub4 + byMu];
    const SIndex& s_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_mu) * _DC_Dir + byMu];
    const SIndex& s_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_nu) * _DC_Dir + byNu];

    //u1^+ u2^+ u3 u4
    //= (u2 u1)^+ u3 u4
    deviceSU3 u(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_nu_m_mu__nu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_mu__mu, byFieldId));
    u.DaggerMul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_nu_m_mu__mu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_nu__nu, byFieldId));

    return u;
}

/**
 * U_{mu,nu}(n)+U^+_{-mu,nu}(n)+U^+_{mu,-nu}(n)+U_{-mu,-nu}(n)
 * or
 * U_{mu,nu}(n)+U_{nu,-mu}(n)+U_{-nu,mu}(n)+U_{-mu,-nu}(n) <--- we are using this one
 * or
 * U_{mu,nu}(n)+U_{mu,nu}(n-mu)+U_{mu,nu}(n-nu)+U_{mu,nu}(n-mu-nu)
 *
 */
static __device__ __inline__ deviceSU3 _deviceClover(const deviceSU3* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    deviceSU3 ret(_device1PlaqutteTermPP(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId));
    ret.Add(_device1PlaqutteTermMM(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId));
    ret.Add(_device1PlaqutteTermPM(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));
    ret.Add(_device1PlaqutteTermMP(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));

    return ret;
}

/**
 * Avoid the add of matrices
 */
static __device__ __inline__ Real _deviceCloverRetr(const deviceSU3* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    return _device1PlaqutteTermPP(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId).ReTr()
        + _device1PlaqutteTermMM(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId).ReTr()
        + _device1PlaqutteTermPM(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId).ReTr()
        + _device1PlaqutteTermMP(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId).ReTr();
}

#pragma endregion


#pragma region Helper device functions

static __device__ __inline__ CLGComplex _deviceLinkU1(
    const CLGComplex* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    CLGComplex sRet = _onec;
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet.y = -sRet.y;
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    const CLGComplex& toMul = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                    sRet = _make_cuComplex(sRet.x * toMul.x + sRet.y * toMul.y, sRet.y * toMul.x - sRet.x * toMul.y);
                }
                else
                {
                    const CLGComplex& toMul = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                    sRet = _make_cuComplex(sRet.x * toMul.x - sRet.y * toMul.y, sRet.y * toMul.x + sRet.x * toMul.y);
                }
            }
        }

        if (pDir[i] > 0) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return sRet;
}

static __device__ __inline__ Real _deviceLinkU1ArgSum(
    const CLGComplex* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    Real sRet = F(0.0);
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = __cuCargf(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet = -sRet;
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet -= __cuCargf(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    sRet += __cuCargf(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
            }
        }

        if (pDir[i] > 0) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return sRet;
}


static __device__ __inline__ CLGComplex _deviceGetGaugeBCU1DirSIndex(
    const CLGComplex* __restrict__ pBuffer,
    const SIndex& idx,
    BYTE byFieldId)
{
    CLGComplex ret = idx.IsDirichlet() ?
        ((CFieldBoundaryGaugeU1*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
    if (idx.NeedToDagger())
    {
        ret = _cuConjf(ret);
    }
    return ret;
}




/**
* big index is the index of walking table.
* The plaqutte index may not be cached because n may out of boundary, so we calculate every one
* n, n+mu, n+nu, n
*
*   <----- ^
*   |      |
*   |      |
*   V      |
* O ------->
*/
static __device__ __inline__ CLGComplex _device1PlaqutteTermPPU1(
    const CLGComplex* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    //For any boundary condition it is always, site->mu, site_p_mu->nu, site_p_nu->mu+, site->nu+
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite4, byMu + 1);
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite4, byNu + 1);
    const UINT uiB4 = uiBigIdx * _DC_Dir;
    const SIndex& s_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiB4 + byMu];
    const SIndex& s_p_mu_nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_p_mu) * _DC_Dir + byNu];
    const SIndex& s_p_nu_mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_p_nu) * _DC_Dir + byMu];
    const SIndex& s_nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiB4 + byNu];

    CLGComplex u = _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_mu, byFieldId);
    u = _cuCmulf(u, _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_p_mu_nu, byFieldId));
    u = _cuCmulf(u, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, s_p_nu_mu, byFieldId)));
    u = _cuCmulf(u, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, s_nu, byFieldId)));

    return u;
}

/**
* U(-mu,nu) = U^+_{mu}(N-mu) U_{nu}(N-mu) U_{mu}(N-mu+nu) U^+_{nu}(N)
*
*    ------->
*    ^      |
*    |      |
*    |      V
*    <------- O
*/
static __device__ __inline__ CLGComplex _device1PlaqutteTermMPU1(
    const CLGComplex* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1);
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, byNu + 1);
    const UINT uin_m_mub4 = __idx->_deviceGetBigIndex(n_m_mu) * _DC_Dir;
    const SIndex& s_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byMu];
    const SIndex& s_m_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byNu];
    const SIndex& s_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_mu_p_nu) * _DC_Dir + byMu];
    const SIndex& s__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byNu];

    CLGComplex u = _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_mu__mu, byFieldId);
    u = _cuCmulf(_cuConjf(u), _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_mu__nu, byFieldId));
    u = _cuCmulf(u, _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_mu_p_nu__mu, byFieldId));
    u = _cuCmulf(u, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, s__nu, byFieldId)));

    return u;
}

/**
* U(mu,-nu) = U(N) U^+(N+mu-nu) U^+(N-nu) U(N-nu)
*
* O  ------->
*    ^      |
*    |      |
*    |      V
*    <-------
*/
static __device__ __inline__ CLGComplex _device1PlaqutteTermPMU1(
    const CLGComplex* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byNu) - 1);
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(n_m_nu, byMu + 1);
    const UINT uin_m_nub4 = __idx->_deviceGetBigIndex(n_m_nu) * _DC_Dir;
    const SIndex& s_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byMu];
    const SIndex& s_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byNu];
    const SIndex& s_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_nu_p_mu) * _DC_Dir + byNu];
    const SIndex& s__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byMu];

    CLGComplex u = _deviceGetGaugeBCU1DirSIndex(pDeviceData, s__mu, byFieldId);
    u = _cuCmulf(u, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_nu_p_mu__nu, byFieldId)));
    u = _cuCmulf(u, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_nu__mu, byFieldId)));
    u = _cuCmulf(u, _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_nu__nu, byFieldId));

    return u;
}

/**
* U(-mu,-nu) = U^+(N-mu) U^+(N-mu-nu) U(N-mu-nu) U(N-nu)
*
* <----- ^ O
* |      |
* |      |
* V      |
* ------->
*/
static __device__ __inline__ CLGComplex _device1PlaqutteTermMMU1(
    const CLGComplex* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1);
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byNu) - 1);
    const SSmallInt4 n_m_nu_m_mu = _deviceSmallInt4OffsetC(n_m_nu, -static_cast<INT>(byMu) - 1);
    const UINT uin_m_nu_m_mub4 = __idx->_deviceGetBigIndex(n_m_nu_m_mu) * _DC_Dir;

    const SIndex& s_m_nu_m_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nu_m_mub4 + byNu];
    const SIndex& s_m_nu_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nu_m_mub4 + byMu];
    const SIndex& s_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_mu) * _DC_Dir + byMu];
    const SIndex& s_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_nu) * _DC_Dir + byNu];

    //u1^+ u2^+ u3 u4
    //= (u2 u1)^+ u3 u4
    CLGComplex u = _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_nu_m_mu__nu, byFieldId);
    u = _cuCmulf(u, _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_mu__mu, byFieldId));
    u = _cuCmulf(_cuConjf(u), _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_nu_m_mu__mu, byFieldId));
    u = _cuCmulf(u, _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_nu__nu, byFieldId));

    return u;
}


static __device__ __inline__ Real _deviceCloverRetrU1(const CLGComplex* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    return _device1PlaqutteTermPPU1(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId).x
        + _device1PlaqutteTermMMU1(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId).x
        + _device1PlaqutteTermPMU1(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId).x
        + _device1PlaqutteTermMPU1(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId).x;
}

#pragma endregion

#pragma region Helper device functions U1 Real

static __device__ __inline__ Real _deviceLinkU1Real(
    const Real* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    Real fRet = F(0.0);
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                fRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    fRet = -fRet;
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                const Real& toAdd = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    fRet = fRet - toAdd;
                }
                else
                {
                    fRet = fRet + toAdd;
                }
            }
        }

        if (pDir[i] > 0) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return fRet;
}

static __device__ __inline__ deviceSU3 _deviceLinkEM(
    const deviceSU3* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataReal,
    Real fCharge,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceSU3 sRet = deviceSU3::makeSU3Id();
    Real fRet = F(0.0);
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                fRet = pDeviceDataReal[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet.Dagger();
                    fRet = -fRet;
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                const Real& toAdd = pDeviceDataReal[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet.MulDagger(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                    fRet = fRet - toAdd;
                }
                else
                {
                    sRet.Mul(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                    fRet = fRet + toAdd;
                }
            }
        }

        if (pDir[i] > 0 && i < (byLength - 1)) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }
    fRet = fRet * fCharge;
    sRet.MulComp(_make_cuComplex(_cos(fRet), _sin(fRet)));
    return sRet;
}

#pragma endregion


#pragma region device functions tree improved

/**
* Rectangle clover
*
*
* It sums over:
*
* -------
* |     |
* ---x---
*
* ---x---
* |     |
* -------
*
* ----
* |  |
* x  |
* |  |
* ----
*
* ----
* |  |
* |  x
* |  |
* ----
*/
static __device__ __inline__ Real _deviceOneRectangleRetr(
    const BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite, INT iMu, INT iNu)
{
    INT path[6] = { iMu, iNu, -iMu, -iMu, -iNu, iMu };
    return _deviceLink(pDeviceData, sSite, 6, byFieldId, path).ReTr();
}

static __device__ __inline__ Real _deviceCloverRectangleRetr(
    const BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite, BYTE byMu, BYTE byNu)
{
    const INT ifwdMu = __fwd(byMu);
    const INT ifwdNu = __fwd(byNu);
    Real fRes = _deviceOneRectangleRetr(byFieldId, pDeviceData, sSite, ifwdMu, ifwdNu);
    fRes += _deviceOneRectangleRetr(byFieldId, pDeviceData, sSite, ifwdMu, -ifwdNu);
    fRes += _deviceOneRectangleRetr(byFieldId, pDeviceData, sSite, ifwdNu, ifwdMu);
    fRes += _deviceOneRectangleRetr(byFieldId, pDeviceData, sSite, ifwdNu, -ifwdMu);
    return fRes;
}

#pragma endregion


//================= Put those device functions to header file because we will use them ==============

#pragma region device function Rotation

#pragma region Energy

#pragma region Plaqutte term

/**
* Product of 3 terms
* U(uiBIa)_{byDira} . U(uiBIb)_{byDirb} . U(uiBIc)_{byDirc}
* To calcuate staple
*/
static __device__ __inline__ deviceSU3 _deviceGetSTTerm(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex& linkA, const SIndex& linkB, const SIndex& linkC)
{
    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, linkA, byFieldId));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, linkB, byFieldId));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, linkC, byFieldId));
    return ret;
}

static __device__ __inline__ Real _device1PlaqutteTermReTr(
    const deviceSU3* __restrict__ pDeviceData, BYTE byFieldId,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4)
{
    return _device1PlaqutteTermPP(pDeviceData, byMu, byNu, uiBigIdx, sSite4, byFieldId).ReTr();
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
static __device__ __inline__ Real _device4PlaqutteTerm(const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIndex, const SSmallInt4& sSite4, BYTE byFieldId)
{
    return F(3.0) - F(0.25) * _deviceCloverRetr(pDeviceData, sSite4, uiBigIndex, byMu, byNu, byFieldId);
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
    deviceSU3 term1(_deviceGetSTTerm(byFieldId, pDeviceData,
        n__mu, n_p_mu__nu, n_p_nu__mu_dag));
    //uiBigIndex, uiN_p_mu, uiN_p_nu, mu, nu, mu, 0, 0, 1));

//U^+_{mu}(N-mu) U_{nu}(N-mu) U_{mu}(N-mu+nu)
    term1.Sub(_deviceGetSTTerm(byFieldId, pDeviceData,
        n_m_mu__mu_dag, n_m_mu__nu, n_m_mu_p_nu__mu));
    //uiN_m_mu, uiN_m_mu, uiN_m_mu_p_nu, mu, nu, mu, 1, 0, 0));

// U_{rho}(N+nu) U^+_{nu}(N+rho) U^+_{rho}(N)
    deviceSU3 term2(_deviceGetSTTerm(byFieldId, pDeviceData,
        n_p_nu__rho, n_p_rho__nu_dag, n__rho_dag));
    //uiN_p_nu, uiN_p_rho, uiBigIndex, rho, nu, rho, 0, 1, 1));

// U^+_{rho}(N+nu-rho) U^+_{nu}(N-rho) U_{rho}(N-rho)
    term2.Sub(_deviceGetSTTerm(byFieldId, pDeviceData,
        n_m_rho_p_nu__rho_dag, n_m_rho__nu_dag, n_m_rho__rho));
    //uiN_m_rho_p_nu, uiN_m_rho, uiN_m_rho, rho, nu, rho, 1, 1, 0));

    term1.Mul(term2);

    //pm mu, nu
    //U(mu,-nu) = U(N) U(N+mu-nu) U(N-nu) U(N-nu), 0110
    deviceSU3 term3(_deviceGetSTTerm(byFieldId, pDeviceData,
        n__mu, n_p_mu_m_nu__nu_dag, n_m_nu__mu_dag));
    //uiBigIndex, uiN_p_mu_m_nu, uiN_m_nu, mu, nu, mu, 0, 1, 1));

//mm
//U(-mu, -nu) = U(N - mu) U(N - mu - nu) U(N - mu - nu) U(N - nu) 1100
    term3.Sub(_deviceGetSTTerm(byFieldId, pDeviceData,
        n_m_mu__mu_dag, n_m_mu_m_nu__nu_dag, n_m_mu_m_nu__mu));
    //uiN_m_mu, uiN_m_mu_m_nu, uiN_m_mu_m_nu, mu, nu, mu, 1, 1, 0));

//mp, nu, rho
//mp = U(-mu,nu) = U^+_{mu}(N-mu) U_{nu}(N-mu) U_{mu}(N-mu+nu) U^+_{nu}(N)
    deviceSU3 term4(_deviceGetSTTerm(byFieldId, pDeviceData,
        n_m_nu__rho, n_m_nu_p_rho__nu, n__rho_dag));
    //uiN_m_nu, uiN_m_nu_p_rho, uiBigIndex, rho, nu, rho, 0, 0, 1));

//mm nu rho
//U(-mu, -nu) = U(N - mu) U(N - mu - nu) U(N - mu - nu) U(N - nu) 1100
    term4.Sub(_deviceGetSTTerm(byFieldId, pDeviceData,
        n_m_rho_m_nu__rho_dag, n_m_rho_m_nu__nu, n_m_rho__rho));
    //uiN_m_rho_m_nu, uiN_m_rho_m_nu, uiN_m_rho, rho, nu, rho, 1, 0, 0));

    term3.Mul(term4);

    term1.Add(term3);

    return term1.ReTr();
}

#pragma endregion

#pragma endregion

#pragma region Force

#pragma region Plaqutte term



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

/**
 * Staple for U_mu, from U_{mu,nu}
 */
static __device__ __inline__ deviceSU3 _deviceStapleTermGfactor(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceSU3* __restrict__ pDeviceData,
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

    deviceSU3 left(
        _deviceGetSTTerm(byFieldId, pDeviceData,
            //pDeviceData, uiBigIndex, uiN_p_nu, uiN_p_mu, nu, mu, nu, 0, 0, 1
            n__nu, n_p_nu__mu, n_p_mu__nu_dag
        ));
    deviceSU3 right(
        _deviceGetSTTerm(byFieldId, pDeviceData,
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

    left.MulReal(fLFactor * fOmegaSq);
    right.MulReal(fRFactor * fOmegaSq);
    left.Add(right);

    return left;
}

#pragma endregion

#pragma region Chair terms

/**
* U(N) U(N+rho) U(N+nu) - U(N-rho) U(N-rho) U(N-rho+nu)
* rho nu rho
* + + -, - + +
*/
static __device__ __inline__ deviceSU3 _deviceS1(BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
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

    deviceSU3 left(
        _deviceGetSTTerm(byFieldId, pDeviceData,
            //pDeviceData, uiBigIndex, uiN_p_rho, uiN_p_nu, rho, nu, rho, 0, 0, 1
            n__rho, n_p_rho__nu, n_p_nu__rho_dag
        ));
    left.Sub(
        _deviceGetSTTerm(byFieldId, pDeviceData,
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
static __device__ __inline__ deviceSU3 _deviceS2(BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
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

    deviceSU3 left(
        _deviceGetSTTerm(byFieldId, pDeviceData,
            n__rho, n_m_nu_p_rho__nu_dag, n_m_nu__rho_dag
            //pDeviceData, uiBigIndex, uiN_m_nu_p_rho, uiN_m_nu, rho, nu, rho, 0, 1, 1
        ));
    left.Sub(
        _deviceGetSTTerm(byFieldId, pDeviceData,
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
static __device__ __inline__ deviceSU3 _deviceS3(BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
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

    deviceSU3 left(
        _deviceGetSTTerm(byFieldId, pDeviceData,
            //pDeviceData, uiN_p_mu_m_rho_p_nu, uiN_p_mu_m_rho, uiN_p_mu_m_rho, rho, nu, rho, 1, 1, 0
            n_p_mu_m_rho_p_nu__rho_dag, n_p_mu_m_rho__nu_dag, n_p_mu_m_rho__rho
        ));
    left.Sub(
        _deviceGetSTTerm(byFieldId, pDeviceData,
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
static __device__ __inline__ deviceSU3 _deviceS4(BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
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

    deviceSU3 left(
        _deviceGetSTTerm(byFieldId, pDeviceData,
            n_p_mu_m_rho_m_nu__rho_dag, n_p_mu_m_rho_m_nu__nu, n_p_mu_m_rho__rho
            //pDeviceData, uiN_p_mu_m_rho_m_nu, uiN_p_mu_m_rho_m_nu, uiN_p_mu_m_rho, rho, nu, rho, 1, 0, 0
        ));
    left.Sub(
        _deviceGetSTTerm(byFieldId, pDeviceData,
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
static __device__ __inline__ deviceSU3 _deviceT1(BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
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

    deviceSU3 left(
        _deviceGetSTTerm(byFieldId, pDeviceData,
            n_p_mu_m_rho__rho_dag, n_p_mu_m_rho__nu, n_p_mu_m_rho_p_nu__rho
            //pDeviceData, uiN_p_mu_m_rho, uiN_p_mu_m_rho, uiN_p_mu_m_rho_p_nu, rho, nu, rho, 1, 0, 0
        ));
    left.Sub(
        _deviceGetSTTerm(byFieldId, pDeviceData,
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
static __device__ __inline__ deviceSU3 _deviceT2(BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData, const SSmallInt4& sSite,
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

    deviceSU3 left(
        _deviceGetSTTerm(byFieldId, pDeviceData,
            n_m_mu__rho, n_m_mu_p_rho__nu, n_m_mu_p_nu__rho_dag
            //pDeviceData, uiN_m_mu, uiN_m_mu_p_rho, uiN_m_mu_p_nu, rho, nu, rho, 0, 0, 1
        ));
    left.Sub(
        _deviceGetSTTerm(byFieldId, pDeviceData,
            n_m_mu_m_rho__rho_dag, n_m_mu_m_rho__nu, n_m_mu_p_nu_m_rho__rho
            //pDeviceData, uiN_m_mu_m_rho, uiN_m_mu_m_rho, uiN_m_mu_p_nu_m_rho, rho, nu, rho, 1, 0, 0
        ));
    return left;
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

/**
* [hi(n)+hi(n+nu)]S1  U(N+nu) U(N+mu)
* mu nu
* - +,
*/
static __device__ __inline__ deviceSU3 _deviceStapleS1(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, mu + 1);
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, nu + 1);
    const UINT uiN_p_nu = __bi(n_p_nu);
    const SIndex& n_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + nu];
    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_p_nu * _DC_Dir + mu];
    const SIndex& uiSiteN_p_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_nu];

    deviceSU3 ret(_deviceS1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));
    ret.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_mu__nu, byFieldId));

    ret.MulReal(
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
static __device__ __inline__ deviceSU3 _deviceStapleS2(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(n_m_nu, __fwd(mu));
    const UINT uiN_m_nu = __bi(n_m_nu);
    const SIndex& n_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_m_nu * _DC_Dir + mu];
    const SIndex& n_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_nu_p_mu) + nu];

    const SIndex& uiSiteN_m_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_m_nu];

    deviceSU3 ret(_deviceS2(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu__mu, byFieldId));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu_p_mu__nu, byFieldId));

    ret.MulReal(
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
static __device__ __inline__ deviceSU3 _deviceStapleS3(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
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

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n__nu, byFieldId));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));
    ret.Mul(_deviceS3(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

    ret.MulReal(
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
static __device__ __inline__ deviceSU3 _deviceStapleS4(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
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

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu__nu, byFieldId));
    ret.DaggerMul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu__mu, byFieldId));
    ret.Mul(_deviceS4(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

    ret.MulReal(
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
static __device__ __inline__ deviceSU3 _deviceStapleT1(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
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

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n__mu, byFieldId));
    ret.Mul(_deviceT1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));

    ret.MulReal(
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
static __device__ __inline__ deviceSU3 _deviceStapleT2(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
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

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_mu__mu, byFieldId));
    ret.DaggerMul(_deviceT2(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_mu_p_nu__mu, byFieldId));

    ret.MulReal(
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
static __device__ __inline__ deviceSU3 _deviceStapleChairTerm1(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    deviceSU3 ret(_deviceStapleS1(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    ret.Add(_deviceStapleS2(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    ret.Add(_deviceStapleS3(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    ret.Add(_deviceStapleS4(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    return ret;
}

/**
* i = 0, 1, 2 for coefficient
* _deviceChairTerm1,2,3 for partial mu, nu, rho
* It is (1/8) * (T1+T2 + T1(mu<->rho) + T2(mu<->rho))
*/
static __device__ __inline__ deviceSU3 _deviceStapleChairTerm2(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    deviceSU3 ret(_deviceStapleT1(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    ret.Add(_deviceStapleT2(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    ret.Add(_deviceStapleT1(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, rho, nu, mu, fpt));
    ret.Add(_deviceStapleT2(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, rho, nu, mu, fpt));
    return ret;
}

#pragma endregion

#pragma endregion

#pragma region Projective plane related

#pragma region Chair

//=============================
//The shifted coord should be conflict with Dirichlet, so we do not consider it
//This is for projective plane
//typedef Real _deviceCoeffPeriodic(
//    BYTE byFieldId,
//    const SSmallInt4& center,
//    SSmallInt4 site);

/*
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

//static __device__ __inline__ Real _deviceHiPeriodic(
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


/*
static __device__ __inline__ deviceSU3 _deviceStapleS1Periodic(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffPeriodic fpt)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, mu + 1);
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, nu + 1);

    const SIndex& n_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + nu];
    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + mu];

    deviceSU3 ret(_deviceS1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));
    ret.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_mu__nu, byFieldId));

    ret.MulReal(_deviceHiPeriodic(byFieldId, sCenter, sSite, n_p_nu, fpt));

    return ret;
}

static __device__ __inline__ deviceSU3 _deviceStapleS2Periodic(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffPeriodic fpt)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(n_m_nu, __fwd(mu));
    const SIndex& n_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_nu) + mu];
    const SIndex& n_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_nu_p_mu) + nu];

    deviceSU3 ret(_deviceS2(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu__mu, byFieldId));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu_p_mu__nu, byFieldId));

    ret.MulReal(_deviceHiPeriodic(byFieldId, sCenter, sSite, n_m_nu, fpt));

    return ret;
}

static __device__ __inline__ deviceSU3 _deviceStapleS3Periodic(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffPeriodic fpt)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_mu_p_nu = _deviceSmallInt4OffsetC(n_p_mu, __fwd(nu));
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));

    const SIndex& n__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIndex * _DC_Dir + nu];
    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + mu];

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n__nu, byFieldId));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));
    ret.Mul(_deviceS3(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

    ret.MulReal(_deviceHiPeriodic(byFieldId, sCenter, n_p_mu, n_p_mu_p_nu, fpt));

    return ret;

}

static __device__ __inline__ deviceSU3 _deviceStapleS4Periodic(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffPeriodic fpt)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_mu_m_nu = _deviceSmallInt4OffsetC(n_p_mu, __bck(nu));
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));

    const UINT n_m_nu_bi4 = __bi4(n_m_nu);

    const SIndex& n_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + nu];
    const SIndex& n_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][n_m_nu_bi4 + mu];

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu__nu, byFieldId));
    ret.DaggerMul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_nu__mu, byFieldId));
    ret.Mul(_deviceS4(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

    ret.MulReal(_deviceHiPeriodic(byFieldId, sCenter, n_p_mu, n_p_mu_m_nu, fpt));

    return ret;
}

static __device__ __inline__ deviceSU3 _deviceStapleT1Periodic(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffPeriodic fpt)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, __fwd(mu));
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, __fwd(nu));
    const SSmallInt4 n_p_mu_p_nu = _deviceSmallInt4OffsetC(n_p_mu, __fwd(nu));

    const SIndex& n__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIndex * _DC_Dir + mu];
    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_nu) + mu];

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n__mu, byFieldId));
    ret.Mul(_deviceT1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));

    ret.MulReal(_deviceHiPeriodic(byFieldId, sCenter, n_p_mu, n_p_mu_p_nu, fpt));

    return ret;
}

static __device__ __inline__ deviceSU3 _deviceStapleT2Periodic(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffPeriodic fpt)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite, __bck(mu));
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, __fwd(nu));

    const SIndex& n_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_mu) + mu];
    const SIndex& n_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_mu_p_nu) + mu];

    deviceSU3 ret(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_mu__mu, byFieldId));
    ret.DaggerMul(_deviceT2(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, n_m_mu_p_nu__mu, byFieldId));

    ret.MulReal(_deviceHiPeriodic(byFieldId, sCenter, n_m_mu, n_m_mu_p_nu, fpt));

    return ret;
}

static __device__ __inline__ deviceSU3 _deviceStapleChairTerm1Periodic(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffPeriodic fpt)
{
    deviceSU3 ret(_deviceStapleS1Periodic(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    ret.Add(_deviceStapleS2Periodic(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    ret.Add(_deviceStapleS3Periodic(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    ret.Add(_deviceStapleS4Periodic(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    return ret;
}

static __device__ __inline__ deviceSU3 _deviceStapleChairTerm2Periodic(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffPeriodic fpt)
{
    deviceSU3 ret(_deviceStapleT1Periodic(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    ret.Add(_deviceStapleT2Periodic(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    ret.Add(_deviceStapleT1Periodic(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, rho, nu, mu, fpt));
    ret.Add(_deviceStapleT2Periodic(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, rho, nu, mu, fpt));
    return ret;
}
*/

#pragma endregion

#pragma endregion

#pragma endregion


//================= Put those device functions to header file because we will use them ==============

#pragma region device function rotation U1

#pragma region Energy

#pragma region Plaqutte term

/**
* Product of 3 terms
* U(uiBIa)_{byDira} . U(uiBIb)_{byDirb} . U(uiBIc)_{byDirc}
* To calcuate staple
*/
static __device__ __inline__ CLGComplex _deviceGetSTTermU1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    const SIndex& linkA, const SIndex& linkB, const SIndex& linkC)
{
    CLGComplex ret = _deviceGetGaugeBCU1DirSIndex(pDeviceData, linkA, byFieldId);
    ret = _cuCmulf(ret, _deviceGetGaugeBCU1DirSIndex(pDeviceData, linkB, byFieldId));
    ret = _cuCmulf(ret, _deviceGetGaugeBCU1DirSIndex(pDeviceData, linkC, byFieldId));
    return ret;
}

//static __device__ __inline__ Real _device1PlaqutteTermReTr(
//    const deviceSU3* __restrict__ pDeviceData, BYTE byFieldId,
//    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4)
//{
//    return _device1PlaqutteTermPP(pDeviceData, byMu, byNu, uiBigIdx, sSite4, byFieldId).ReTr();
//}

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
static __device__ __inline__ Real _device4PlaqutteTermU1(const CLGComplex* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIndex, const SSmallInt4& sSite4, BYTE byFieldId)
{
    return F(1.0) - F(0.25) * _deviceCloverRetrU1(pDeviceData, sSite4, uiBigIndex, byMu, byNu, byFieldId);
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
static __device__ __inline__ Real _deviceChairTermU1(const CLGComplex* __restrict__ pDeviceData,
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
    CLGComplex term1 = _deviceGetSTTermU1(byFieldId, pDeviceData,
        n__mu, n_p_mu__nu, n_p_nu__mu_dag);
    //uiBigIndex, uiN_p_mu, uiN_p_nu, mu, nu, mu, 0, 0, 1));

//U^+_{mu}(N-mu) U_{nu}(N-mu) U_{mu}(N-mu+nu)
    term1 = _cuCsubf(term1, _deviceGetSTTermU1(byFieldId, pDeviceData,
        n_m_mu__mu_dag, n_m_mu__nu, n_m_mu_p_nu__mu));
    //uiN_m_mu, uiN_m_mu, uiN_m_mu_p_nu, mu, nu, mu, 1, 0, 0));

// U_{rho}(N+nu) U^+_{nu}(N+rho) U^+_{rho}(N)
    CLGComplex term2 = _deviceGetSTTermU1(byFieldId, pDeviceData,
        n_p_nu__rho, n_p_rho__nu_dag, n__rho_dag);
    //uiN_p_nu, uiN_p_rho, uiBigIndex, rho, nu, rho, 0, 1, 1));

// U^+_{rho}(N+nu-rho) U^+_{nu}(N-rho) U_{rho}(N-rho)
    term2 = _cuCsubf(term2, _deviceGetSTTermU1(byFieldId, pDeviceData,
        n_m_rho_p_nu__rho_dag, n_m_rho__nu_dag, n_m_rho__rho));
    //uiN_m_rho_p_nu, uiN_m_rho, uiN_m_rho, rho, nu, rho, 1, 1, 0));

    term1 = _cuCmulf(term1, term2);

    //pm mu, nu
    //U(mu,-nu) = U(N) U(N+mu-nu) U(N-nu) U(N-nu), 0110
    CLGComplex term3 = _deviceGetSTTermU1(byFieldId, pDeviceData,
        n__mu, n_p_mu_m_nu__nu_dag, n_m_nu__mu_dag);
    //uiBigIndex, uiN_p_mu_m_nu, uiN_m_nu, mu, nu, mu, 0, 1, 1));

//mm
//U(-mu, -nu) = U(N - mu) U(N - mu - nu) U(N - mu - nu) U(N - nu) 1100
    term3 = _cuCsubf(term3, _deviceGetSTTermU1(byFieldId, pDeviceData,
        n_m_mu__mu_dag, n_m_mu_m_nu__nu_dag, n_m_mu_m_nu__mu));
    //uiN_m_mu, uiN_m_mu_m_nu, uiN_m_mu_m_nu, mu, nu, mu, 1, 1, 0));

//mp, nu, rho
//mp = U(-mu,nu) = U^+_{mu}(N-mu) U_{nu}(N-mu) U_{mu}(N-mu+nu) U^+_{nu}(N)
    CLGComplex term4 = _deviceGetSTTermU1(byFieldId, pDeviceData,
        n_m_nu__rho, n_m_nu_p_rho__nu, n__rho_dag);
    //uiN_m_nu, uiN_m_nu_p_rho, uiBigIndex, rho, nu, rho, 0, 0, 1));

//mm nu rho
//U(-mu, -nu) = U(N - mu) U(N - mu - nu) U(N - mu - nu) U(N - nu) 1100
    term4 = _cuCsubf(term4, _deviceGetSTTermU1(byFieldId, pDeviceData,
        n_m_rho_m_nu__rho_dag, n_m_rho_m_nu__nu, n_m_rho__rho));
    //uiN_m_rho_m_nu, uiN_m_rho_m_nu, uiN_m_rho, rho, nu, rho, 1, 0, 0));

    term3 = _cuCmulf(term3, term4);

    term1 = _cuCaddf(term1, term3);

    return term1.x;
}

#pragma endregion

#pragma endregion

#pragma region Force

#pragma region Plaqutte term



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
//static __device__ __inline__ Real _deviceGi(
//    const SSmallInt4& sCenter,
//    const SSmallInt4& sSite,
//    const SSmallInt4& sSiteOffset,
//    const SIndex& uiSiteBI,
//    const SIndex& uiSiteOffsetBI,
//    BYTE i,
//    Real fOmegaSq)
//{
//    if (0 == i)
//    {
//        const Real fX = uiSiteBI.IsDirichlet() ? F(0.0)
//            : static_cast<Real>(sSite.x - sCenter.x);
//        return F(0.5) * fOmegaSq * (fX * fX);
//        //const Real fXp1 = uiSiteOffsetBI.IsDirichlet() ? F(0.0)
//        //    : static_cast<Real>(sSiteOffset.x - sCenter.x);
//        //return F(0.5) * fOmegaSq * (fX * fX + fXp1 * fXp1);
//    }
//    else if (1 == i)
//    {
//        const Real fY = uiSiteBI.IsDirichlet() ? F(0.0)
//            : static_cast<Real>(sSite.y - sCenter.y);
//        return F(0.5) * fOmegaSq * (fY * fY);
//        //const Real fYp1 = uiSiteOffsetBI.IsDirichlet() ? F(0.0)
//        //    : static_cast<Real>(sSiteOffset.y - sCenter.y);
//        //return F(0.5) * fOmegaSq * (fY * fY + fYp1 * fYp1);
//    }
//    const Real fX = uiSiteBI.IsDirichlet() ? F(0.0)
//        : static_cast<Real>(sSite.x - sCenter.x);
//    const Real fY = uiSiteOffsetBI.IsDirichlet() ? F(0.0)
//        : static_cast<Real>(sSite.y - sCenter.y);
//    const Real fXp1 = uiSiteOffsetBI.IsDirichlet() ? F(0.0)
//        : static_cast<Real>(sSiteOffset.x - sCenter.x);
//    const Real fYp1 = uiSiteOffsetBI.IsDirichlet() ? F(0.0)
//        : static_cast<Real>(sSiteOffset.y - sCenter.y);
//    return F(0.5) * fOmegaSq * (fX * fX + fY * fY + fXp1 * fXp1 + fYp1 * fYp1);
//}

/**
* Coefficient = (f(n)+f(n+mu)+f(n+nu)+f(n+mu+nu))/4
* Simplfy: nu is always t direction, so f(n) = f(n+nu), f(n+mu) = f(n+mu+nu)
* Coefficient = (f(n)+f(n+mu))/2
* For 3 == mu, f(n) = f(n+mu)
* This is also true for Dirichlet boundary condition, only Dirichlet on X-Y direction is assumed
*
* ==================================================
* Note for periodic boundary condition:
* For const SSmallInt4 sN_p_m = _deviceSmallInt4OffsetC(sSite4, mu + 1)
* sN_p_m.mu can be -1, which leads to a wrong (sN_p_m.y - sCenter.y)
* This '-1' should be set to L_mu - 1. If we consider add the plaquttes as clovers,
* then the coordinates of the centers of the clovers will always be in the lattice,
* so should be set to L_mu - 1
*/
//static __device__ __inline__ Real _deviceFi(
//    BYTE byFieldId,
//    const SSmallInt4& sSite4,
//    const SSmallInt4& sCenter,
//    UINT uiN, BYTE i, BYTE mu, BYTE nu)
//{
//    const SSmallInt4 sN_p_mu = _deviceSmallInt4OffsetC(sSite4, mu + 1);
//    const SIndex& n_p_mu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_mu)];
//    const SSmallInt4 site_N_p_mu = __deviceSiteIndexToInt4(n_p_mu__idx.m_uiSiteIndex);
//
//    const SSmallInt4 sN_p_nu = _deviceSmallInt4OffsetC(sSite4, nu + 1);
//    const SIndex& n_p_nu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_nu)];
//    const SSmallInt4 site_N_p_nu = __deviceSiteIndexToInt4(n_p_nu__idx.m_uiSiteIndex);
//
//    const SSmallInt4 sN_p_numu = _deviceSmallInt4OffsetC(sN_p_mu, nu + 1);
//    const SIndex& n_p_numu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_numu)];
//    const SSmallInt4 site_N_p_munu = __deviceSiteIndexToInt4(n_p_numu__idx.m_uiSiteIndex);
//
//    const UBOOL bN_surface = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet();
//    const UBOOL bN_p_mu_surface = n_p_mu__idx.IsDirichlet();
//    const UBOOL bN_p_nu_surface = n_p_nu__idx.IsDirichlet();
//    const UBOOL bN_p_munu_surface = n_p_numu__idx.IsDirichlet();
//
//    const INT x1 = bN_surface ? 0 : (sSite4.x - sCenter.x);
//    const INT y1 = bN_surface ? 0 : (sSite4.y - sCenter.y);
//
//    const INT x2 = bN_p_mu_surface ? 0 : (site_N_p_mu.x - sCenter.x);
//    const INT y2 = bN_p_mu_surface ? 0 : (site_N_p_mu.y - sCenter.y);
//
//    const INT x3 = bN_p_nu_surface ? 0 : (site_N_p_nu.x - sCenter.x);
//    const INT y3 = bN_p_nu_surface ? 0 : (site_N_p_nu.y - sCenter.y);
//
//    const INT x4 = bN_p_munu_surface ? 0 : (site_N_p_munu.x - sCenter.x);
//    const INT y4 = bN_p_munu_surface ? 0 : (site_N_p_munu.y - sCenter.y);
//
//    if (0 == i)
//    {
//
//        return F(0.25) * static_cast<Real>(x1 * x1 + y1 * y1
//            + x2 * x2 + y2 * y2
//            + x3 * x3 + y3 * y3
//            + x4 * x4 + y4 * y4);
//    }
//
//    if (1 == i)
//    {
//        return F(0.25) * static_cast<Real>(
//              y1 * y1
//            + y2 * y2
//            + y3 * y3
//            + y4 * y4);
//    }
//    return F(0.25) * static_cast<Real>(
//          x1 * x1
//        + x2 * x2
//        + x3 * x3
//        + x4 * x4);
//}

//static __device__ __inline__ Real _deviceFiShifted(
//    BYTE byFieldId,
//    const SSmallInt4& sSite4,
//    const SSmallInt4& sCenter,
//    BYTE i, BYTE mu, BYTE nu)
//{
//    const SIndex& n__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sSite4)];
//    const SSmallInt4 site_N = __deviceSiteIndexToInt4(n__idx.m_uiSiteIndex);
//
//    const SSmallInt4 sN_p_mu = _deviceSmallInt4OffsetC(sSite4, mu + 1);
//    const SIndex& n_p_mu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_mu)];
//    const SSmallInt4 site_N_p_mu = __deviceSiteIndexToInt4(n_p_mu__idx.m_uiSiteIndex);
//
//    const SSmallInt4 sN_p_nu = _deviceSmallInt4OffsetC(sSite4, nu + 1);
//    const SIndex& n_p_nu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_nu)];
//    const SSmallInt4 site_N_p_nu = __deviceSiteIndexToInt4(n_p_nu__idx.m_uiSiteIndex);
//
//    const SSmallInt4 sN_p_numu = _deviceSmallInt4OffsetC(sN_p_mu, nu + 1);
//    const SIndex& n_p_numu__idx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sN_p_numu)];
//    const SSmallInt4 site_N_p_munu = __deviceSiteIndexToInt4(n_p_numu__idx.m_uiSiteIndex);
//
//    //const UBOOL bN_surface = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet();
//    //const UBOOL bN_p_mu_surface = n_p_mu__idx.IsDirichlet();
//    //const UBOOL bN_p_nu_surface = n_p_nu__idx.IsDirichlet();
//    //const UBOOL bN_p_munu_surface = n_p_numu__idx.IsDirichlet();
//
//    const Real x1 = static_cast<Real>(site_N.x - sCenter.x + F(0.5));
//    const Real y1 = static_cast<Real>(site_N.y - sCenter.y + F(0.5));
//
//    const Real x2 = static_cast<Real>(site_N_p_mu.x - sCenter.x + F(0.5));
//    const Real y2 = static_cast<Real>(site_N_p_mu.y - sCenter.y + F(0.5));
//
//    const Real x3 = static_cast<Real>(site_N_p_nu.x - sCenter.x + F(0.5));
//    const Real y3 = static_cast<Real>(site_N_p_nu.y - sCenter.y + F(0.5));
//
//    const Real x4 = static_cast<Real>(site_N_p_munu.x - sCenter.x + F(0.5));
//    const Real y4 = static_cast<Real>(site_N_p_munu.y - sCenter.y + F(0.5));
//
//    if (0 == i)
//    {
//        //const UBOOL bCorner = (sSite4.x == site_N_p_munu.x) && (sSite4.y == site_N_p_munu.y);
//        //if (bCorner)
//        //{
//        //    return F(0.0);
//        //}
//        return F(0.25) * (x1 * x1 + y1 * y1
//            + x2 * x2 + y2 * y2
//            + x3 * x3 + y3 * y3
//            + x4 * x4 + y4 * y4);
//    }
//
//    if (1 == i)
//    {
//        return F(0.25) * (y1 * y1
//            + y2 * y2
//            + y3 * y3
//            + y4 * y4);
//    }
//    return F(0.25) * (x1 * x1
//        + x2 * x2
//        + x3 * x3
//        + x4 * x4);
//}

/**
 * Staple for U_mu, from U_{mu,nu}
 */
static __device__ __inline__ CLGComplex _deviceStapleTermGfactorU1(
    BYTE byFieldId,
    UBOOL bTorus,
    const CLGComplex* __restrict__ pDeviceData,
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

    CLGComplex left =
        _deviceGetSTTermU1(byFieldId, pDeviceData,
            //pDeviceData, uiBigIndex, uiN_p_nu, uiN_p_mu, nu, mu, nu, 0, 0, 1
            n__nu, n_p_nu__mu, n_p_mu__nu_dag
        );
    CLGComplex right =
        _deviceGetSTTermU1(byFieldId, pDeviceData,
            //pDeviceData, uiN_m_nu, uiN_m_nu, uiN_p_mu_m_nu, nu, mu, nu, 1, 0, 0
            n_m_nu__nu_dag, n_m_nu__mu, n_p_mu_m_nu__nu
        );

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

    left.x = (left.x * fLFactor + right.x * fRFactor) * fOmegaSq;
    left.y = (left.y * fLFactor + right.y * fRFactor) * fOmegaSq;

    return left;
}

#pragma endregion

#pragma region Chair terms

/**
* U(N) U(N+rho) U(N+nu) - U(N-rho) U(N-rho) U(N-rho+nu)
* rho nu rho
* + + -, - + +
*/
static __device__ __inline__ CLGComplex _deviceS1U1(BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData, const SSmallInt4& sSite,
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

    CLGComplex left =
        _deviceGetSTTermU1(byFieldId, pDeviceData,
            //pDeviceData, uiBigIndex, uiN_p_rho, uiN_p_nu, rho, nu, rho, 0, 0, 1
            n__rho, n_p_rho__nu, n_p_nu__rho_dag
        );
    left = _cuCsubf(left,
        _deviceGetSTTermU1(byFieldId, pDeviceData,
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
static __device__ __inline__ CLGComplex _deviceS2U1(BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData, const SSmallInt4& sSite,
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

    CLGComplex left =
        _deviceGetSTTermU1(byFieldId, pDeviceData,
            n__rho, n_m_nu_p_rho__nu_dag, n_m_nu__rho_dag
            //pDeviceData, uiBigIndex, uiN_m_nu_p_rho, uiN_m_nu, rho, nu, rho, 0, 1, 1
        );
    left = _cuCsubf(left,
        _deviceGetSTTermU1(byFieldId, pDeviceData,
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
static __device__ __inline__ CLGComplex _deviceS3U1(BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData, const SSmallInt4& sSite,
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

    CLGComplex left =
        _deviceGetSTTermU1(byFieldId, pDeviceData,
            //pDeviceData, uiN_p_mu_m_rho_p_nu, uiN_p_mu_m_rho, uiN_p_mu_m_rho, rho, nu, rho, 1, 1, 0
            n_p_mu_m_rho_p_nu__rho_dag, n_p_mu_m_rho__nu_dag, n_p_mu_m_rho__rho
        );
    left = _cuCsubf(left,
        _deviceGetSTTermU1(byFieldId, pDeviceData,
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
static __device__ __inline__ CLGComplex _deviceS4U1(BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData, const SSmallInt4& sSite,
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

    CLGComplex left =
        _deviceGetSTTermU1(byFieldId, pDeviceData,
            n_p_mu_m_rho_m_nu__rho_dag, n_p_mu_m_rho_m_nu__nu, n_p_mu_m_rho__rho
            //pDeviceData, uiN_p_mu_m_rho_m_nu, uiN_p_mu_m_rho_m_nu, uiN_p_mu_m_rho, rho, nu, rho, 1, 0, 0
        );
    left = _cuCsubf(left,
        _deviceGetSTTermU1(byFieldId, pDeviceData,
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
static __device__ __inline__ CLGComplex _deviceT1U1(BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData, const SSmallInt4& sSite,
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

    CLGComplex left =
        _deviceGetSTTermU1(byFieldId, pDeviceData,
            n_p_mu_m_rho__rho_dag, n_p_mu_m_rho__nu, n_p_mu_m_rho_p_nu__rho
            //pDeviceData, uiN_p_mu_m_rho, uiN_p_mu_m_rho, uiN_p_mu_m_rho_p_nu, rho, nu, rho, 1, 0, 0
        );
    left = _cuCsubf(left,
        _deviceGetSTTermU1(byFieldId, pDeviceData,
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
static __device__ __inline__ CLGComplex _deviceT2U1(BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData, const SSmallInt4& sSite,
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

    CLGComplex left =
        _deviceGetSTTermU1(byFieldId, pDeviceData,
            n_m_mu__rho, n_m_mu_p_rho__nu, n_m_mu_p_nu__rho_dag
            //pDeviceData, uiN_m_mu, uiN_m_mu_p_rho, uiN_m_mu_p_nu, rho, nu, rho, 0, 0, 1
        );
    left = _cuCsubf(left,
        _deviceGetSTTermU1(byFieldId, pDeviceData,
            n_m_mu_m_rho__rho_dag, n_m_mu_m_rho__nu, n_m_mu_p_nu_m_rho__rho
            //pDeviceData, uiN_m_mu_m_rho, uiN_m_mu_m_rho, uiN_m_mu_p_nu_m_rho, rho, nu, rho, 1, 0, 0
        ));
    return left;
}

/**
* i = 0, 1, 2 correspond to x, y and xy
* h_i(N) = x or y or xy
* return h_i(N) + h_i(N + nu), where N is site, and N + nu (or N + mu or ...) is site2
*/
//static __device__ __inline__ Real _deviceHi(
//    const SSmallInt4 &center,
//    const SSmallInt4 &site, const SSmallInt4 &site2,
//    const SIndex& uiSiteBI, const SIndex& uiSite2BI, BYTE i)
//{
//    if (0 == i)
//    {
//        const Real fX1 = uiSiteBI.IsDirichlet() ? F(0.0)
//            : static_cast<Real>(site.x - center.x);
//        const Real fX2 = uiSite2BI.IsDirichlet() ? F(0.0)
//            : static_cast<Real>(site2.x - center.x);
//        return fX1 + fX2;
//    }
//    else if (1 == i)
//    {
//        const Real fY1 = uiSiteBI.IsDirichlet() ? F(0.0)
//            : static_cast<Real>(site.y - center.y);
//        const Real fY2 = uiSite2BI.IsDirichlet() ? F(0.0)
//            : static_cast<Real>(site2.y - center.y);
//        return -fY1 - fY2;
//    }
//    const Real fX1 = uiSiteBI.IsDirichlet() ? F(0.0)
//        : static_cast<Real>(site.x - center.x);
//    const Real fX2 = uiSite2BI.IsDirichlet() ? F(0.0)
//        : static_cast<Real>(site2.x - center.x);
//    const Real fY1 = uiSiteBI.IsDirichlet() ? F(0.0)
//        : static_cast<Real>(site.y - center.y);
//    const Real fY2 = uiSite2BI.IsDirichlet() ? F(0.0)
//        : static_cast<Real>(site2.y - center.y);
//    return fX1 * fY1 + fX2 * fY2;
//}

#if 0
/**
* [hi(n)+hi(n+nu)]S1  U(N+nu) U(N+mu)
* mu nu
* - +,
*/
static __device__ __inline__ CLGComplex _deviceStapleS1U1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, mu + 1);
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, nu + 1);
    const UINT uiN_p_nu = __bi(n_p_nu);
    const SIndex& n_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + nu];
    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_p_nu * _DC_Dir + mu];
    const SIndex& uiSiteN_p_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_nu];

    CLGComplex ret = _deviceS1U1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho);
    ret = _cuCmulf(ret, _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));
    ret = _cuCmulf(ret, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, n_p_mu__nu, byFieldId)));

    ret = cuCmulf_cr(ret, _deviceHi(sCenter,
        sSite,
        __deviceSiteIndexToInt4(uiSiteN_p_nu.m_uiSiteIndex),
        __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex], uiSiteN_p_nu, i));

    return ret;
}

/**
* [h(N) + h(n-nu)] S2 U(n-nu)U(n+mu-nu)
* mu nu
* + +
*/
static __device__ __inline__ CLGComplex _deviceStapleS2U1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(n_m_nu, __fwd(mu));
    const UINT uiN_m_nu = __bi(n_m_nu);
    const SIndex& n_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_m_nu * _DC_Dir + mu];
    const SIndex& n_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_nu_p_mu) + nu];

    const SIndex& uiSiteN_m_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_m_nu];

    CLGComplex ret = _deviceS2U1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho);
    ret = _cuCmulf(ret, _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_m_nu__mu, byFieldId));
    ret = _cuCmulf(ret, _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_m_nu_p_mu__nu, byFieldId));

    ret = cuCmulf_cr(ret, _deviceHi(sCenter,
        sSite,
        __deviceSiteIndexToInt4(uiSiteN_m_nu.m_uiSiteIndex),
        __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex], uiSiteN_m_nu, i));

    return ret;
}

/**
* [h(N+mu) + h(N+mu+nu)]U(n) U(n+nu) S3
* nu mu
* + +
*/
static __device__ __inline__ CLGComplex _deviceStapleS3U1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
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

    CLGComplex ret = _deviceGetGaugeBCU1DirSIndex(pDeviceData, n__nu, byFieldId);
    ret = _cuCmulf(ret, _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));
    ret = _cuCmulf(ret, _deviceS3U1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

    ret = cuCmulf_cr(ret, _deviceHi(sCenter,
        __deviceSiteIndexToInt4(uiSiteN_p_mu.m_uiSiteIndex),
        __deviceSiteIndexToInt4(uiSiteN_p_mu_p_nu.m_uiSiteIndex),
        uiSiteN_p_mu, uiSiteN_p_mu_p_nu, i));

    return ret;

}

/**
* [h(N+mu) + h(N+mu-nu)] U(n-nu) U(n-nu) S4
* nu mu
* - +
*/
static __device__ __inline__ CLGComplex _deviceStapleS4U1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
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

    CLGComplex ret = _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_m_nu__nu, byFieldId);
    ret = _cuCmulf(_cuConjf(ret), _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_m_nu__mu, byFieldId));
    ret = _cuCmulf(ret, _deviceS4U1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

    ret = cuCmulf_cr(ret, _deviceHi(sCenter,
        __deviceSiteIndexToInt4(uiSiteN_p_mu.m_uiSiteIndex),
        __deviceSiteIndexToInt4(uiSiteN_p_mu_m_nu.m_uiSiteIndex),
        uiSiteN_p_mu, uiSiteN_p_mu_m_nu, i));

    return ret;
}

/**
* [h(n+mu) + h(n+mu+nu)] U(n) T1 U(n+nu)
* mu mu, + -
*
*/
static __device__ __inline__ CLGComplex _deviceStapleT1U1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
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

    CLGComplex ret = _deviceGetGaugeBCU1DirSIndex(pDeviceData, n__mu, byFieldId);
    ret = _cuCmulf(ret, _deviceT1U1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret = _cuCmulf(ret, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, n_p_nu__mu, byFieldId)));

    ret = cuCmulf_cr(ret, _deviceHi(sCenter,
        __deviceSiteIndexToInt4(uiSiteN_p_mu.m_uiSiteIndex),
        __deviceSiteIndexToInt4(uiSiteN_p_mu_p_nu.m_uiSiteIndex),
        uiSiteN_p_mu, uiSiteN_p_mu_p_nu, i));

    return ret;
}

/**
* [h(n-mu) + h(n-mu+nu)] U(n-mu) T2 U(n+nu-mu)
* mu mu, - +
*
*/
static __device__ __inline__ CLGComplex _deviceStapleT2U1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite, __bck(mu));
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, __fwd(nu));

    const UINT uiN_m_mu = __bi(n_m_mu);
    const UINT uiN_m_mu_p_nu = __bi(n_m_mu_p_nu);

    const SIndex& n_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_m_mu * _DC_Dir + mu];
    const SIndex& n_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_m_mu_p_nu * _DC_Dir + mu];

    const SIndex& uiSiteN_m_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_m_mu];
    const SIndex& uiSiteN_m_mu_p_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_m_mu_p_nu];

    CLGComplex ret = _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_m_mu__mu, byFieldId);
    ret = _cuCmulf(_cuConjf(ret), _deviceT2U1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret = _cuCmulf(ret, _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_m_mu_p_nu__mu, byFieldId));

    ret = cuCmulf_cr(ret, _deviceHi(sCenter,
        __deviceSiteIndexToInt4(uiSiteN_m_mu.m_uiSiteIndex),
        __deviceSiteIndexToInt4(uiSiteN_m_mu_p_nu.m_uiSiteIndex),
        uiSiteN_m_mu, uiSiteN_m_mu_p_nu, i));

    return ret;
}

#endif

/**
* i = 0, 1, 2 for coefficient
* _deviceChairTerm1,2,3 for partial mu, nu, rho
* For partial mu, the staple is (1/8)(s1+s2+s3+s4)
*/
//static __device__ __inline__ CLGComplex _deviceStapleChairTerm1U1(
//    BYTE byFieldId,
//    const CLGComplex* __restrict__ pDeviceData,
//    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
//    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
//{
//    CLGComplex ret = _deviceStapleS1U1(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, i);
//    ret = _cuCaddf(ret, _deviceStapleS2U1(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, i));
//    ret = _cuCaddf(ret, _deviceStapleS3U1(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, i));
//    ret = _cuCaddf(ret, _deviceStapleS4U1(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, i));
//    return ret;
//}

/**
* i = 0, 1, 2 for coefficient
* _deviceChairTerm1,2,3 for partial mu, nu, rho
* It is (1/8) * (T1+T2 + T1(mu<->rho) + T2(mu<->rho))
*/
//static __device__ __inline__ CLGComplex _deviceStapleChairTerm2U1(
//    BYTE byFieldId,
//    const CLGComplex* __restrict__ pDeviceData,
//    const SSmallInt4& sCenter, const SSmallInt4& sSite, UINT uiSiteIndex,
//    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, BYTE i)
//{
//    CLGComplex ret = _deviceStapleT1U1(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, i);
//    ret = _cuCaddf(ret, _deviceStapleT1U1(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, i));
//    ret = _cuCaddf(ret, _deviceStapleT1U1(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, rho, nu, mu, i));
//    ret = _cuCaddf(ret, _deviceStapleT1U1(byFieldId, pDeviceData, sCenter, sSite, uiSiteIndex, uiBigIndex, rho, nu, mu, i));
//    return ret;
//}

#pragma endregion

#pragma endregion

#pragma region Projective plane related

#pragma region Chair

//=============================
//The shifted coord should be conflict with Dirichlet, so we do not consider it
//This is for projective plane
//static __device__ __inline__ Real _deviceSiteCoeff(
//    SSmallInt4 sSite4, const SSmallInt4& sCenterSite, BYTE byFieldId, BYTE byType)
//{
//    if (0 == byType)
//    {
//        //x
//        const UBOOL bOpposite = sSite4.x >= static_cast<SBYTE>(_DC_Lx) || sSite4.x < 0;
//        sSite4 = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sSite4)].m_uiSiteIndex);
//        if (bOpposite)
//        {
//            return -sSite4.x + sCenterSite.x - F(0.5);
//        }
//        return sSite4.x - sCenterSite.x + F(0.5);
//    }
//    if (1 == byType)
//    {
//        //y
//        const UBOOL bOpposite = sSite4.y >= static_cast<SBYTE>(_DC_Ly) || sSite4.y < 0;
//        sSite4 = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sSite4)].m_uiSiteIndex);
//        if (bOpposite)
//        {
//            return sSite4.y - sCenterSite.y + F(0.5);
//        }
//        return -sSite4.y + sCenterSite.y - F(0.5);
//    }
//    if (3 == byType)
//    {
//        //There should be NO byType = 3?
//        sSite4 = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sSite4)].m_uiSiteIndex);
//        return -sSite4.y + sCenterSite.y - F(0.5);
//    }
//
//    //byType = 2 and this is XY
//    const BYTE bOppositeX = (sSite4.x >= static_cast<SBYTE>(_DC_Lx) || sSite4.x < 0) ? 1 : 0;
//    const BYTE bOppositeY = (sSite4.y >= static_cast<SBYTE>(_DC_Ly) || sSite4.y < 0) ? 1 : 0;
//    sSite4 = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sSite4)].m_uiSiteIndex);
//    const Real fRet = (sSite4.x - sCenterSite.x + F(0.5)) * (sSite4.y - sCenterSite.y + F(0.5));
//    if (0 != (bOppositeX ^ bOppositeY))
//    {
//        return -fRet;
//    }
//    return fRet;
//}
//
//static __device__ __inline__ Real _deviceHiShifted(
//    BYTE byFieldId,
//    const SSmallInt4& center,
//    SSmallInt4 site, SSmallInt4 site2, BYTE i)
//{
//    return _deviceSiteCoeff(site, center, byFieldId, i) + _deviceSiteCoeff(site2, center, byFieldId, i);
//}


static __device__ __inline__ CLGComplex _deviceStapleS1ShiftedU1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    const SSmallInt4 n_p_mu = _deviceSmallInt4OffsetC(sSite, mu + 1);
    const SSmallInt4 n_p_nu = _deviceSmallInt4OffsetC(sSite, nu + 1);
    const UINT uiN_p_nu = __bi(n_p_nu);
    const SIndex& n_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_p_mu) + nu];
    const SIndex& n_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_p_nu * _DC_Dir + mu];
    const SIndex& uiSiteN_p_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_p_nu];

    CLGComplex ret = _deviceS1U1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho);
    ret = _cuCmulf(ret, _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));
    ret = _cuCmulf(ret, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, n_p_mu__nu, byFieldId)));

    ret = cuCmulf_cr(ret,
        _deviceHi(byFieldId,
            sSite,
            n_p_nu,
            __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex], uiSiteN_p_nu, fpt));

    return ret;
}

static __device__ __inline__ CLGComplex _deviceStapleS2ShiftedU1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite, __bck(nu));
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(n_m_nu, __fwd(mu));
    const UINT uiN_m_nu = __bi(n_m_nu);
    const SIndex& n_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiN_m_nu * _DC_Dir + mu];
    const SIndex& n_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_m_nu_p_mu) + nu];
    const SIndex& uiSiteN_m_nu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN_m_nu];

    CLGComplex ret = _deviceS2U1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho);
    ret = _cuCmulf(ret, _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_m_nu__mu, byFieldId));
    ret = _cuCmulf(ret, _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_m_nu_p_mu__nu, byFieldId));

    ret = cuCmulf_cr(ret, _deviceHi(byFieldId,
        sSite,
        n_m_nu,
        __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex], uiSiteN_m_nu, fpt));

    return ret;
}

static __device__ __inline__ CLGComplex _deviceStapleS3ShiftedU1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
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

    CLGComplex ret = _deviceGetGaugeBCU1DirSIndex(pDeviceData, n__nu, byFieldId);
    ret = _cuCmulf(ret, _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_p_nu__mu, byFieldId));
    ret = _cuCmulf(ret, _deviceS3U1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

    ret = cuCmulf_cr(ret, _deviceHi(byFieldId,
        n_p_mu,
        n_p_mu_p_nu,
        uiSiteN_p_mu, uiSiteN_p_mu_p_nu, fpt));

    return ret;

}

static __device__ __inline__ CLGComplex _deviceStapleS4ShiftedU1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
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

    CLGComplex ret = _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_m_nu__nu, byFieldId);
    ret = _cuCmulf(_cuConjf(ret), _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_m_nu__mu, byFieldId));
    ret = _cuCmulf(ret, _deviceS4U1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

    ret = cuCmulf_cr(ret, _deviceHi(byFieldId,
        n_p_mu,
        n_p_mu_m_nu,
        uiSiteN_p_mu, uiSiteN_p_mu_m_nu, fpt));

    return ret;
}

static __device__ __inline__ CLGComplex _deviceStapleT1ShiftedU1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
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

    CLGComplex ret = _deviceGetGaugeBCU1DirSIndex(pDeviceData, n__mu, byFieldId);
    ret = _cuCmulf(ret, _deviceT1U1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret = _cuCmulf(ret, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, n_p_nu__mu, byFieldId)));

    ret = cuCmulf_cr(ret, _deviceHi(byFieldId,
        n_p_mu,
        n_p_mu_p_nu,
        uiSiteN_p_mu, uiSiteN_p_mu_p_nu, fpt));

    return ret;
}

static __device__ __inline__ CLGComplex _deviceStapleT2ShiftedU1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
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

    CLGComplex ret = _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_m_mu__mu, byFieldId);
    ret = _cuCmulf(_cuConjf(ret), _deviceT2U1(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
    ret = _cuCmulf(ret, _deviceGetGaugeBCU1DirSIndex(pDeviceData, n_m_mu_p_nu__mu, byFieldId));

    ret = cuCmulf_cr(ret, _deviceHi(byFieldId,
        n_m_mu,
        n_m_mu_p_nu,
        uiSiteN_m_mu, uiSiteN_m_mu_p_nu, fpt));

    return ret;
}

static __device__ __inline__ CLGComplex _deviceStapleChairTerm1ShiftedU1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    CLGComplex ret = _deviceStapleS1ShiftedU1(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt);
    ret = _cuCaddf(ret, _deviceStapleS2ShiftedU1(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    ret = _cuCaddf(ret, _deviceStapleS3ShiftedU1(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    ret = _cuCaddf(ret, _deviceStapleS4ShiftedU1(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    return ret;
}

static __device__ __inline__ CLGComplex _deviceStapleChairTerm2ShiftedU1(
    BYTE byFieldId,
    const CLGComplex* __restrict__ pDeviceData,
    const SSmallInt4& sSite, UINT uiSiteIndex,
    UINT uiBigIndex, BYTE mu, BYTE nu, BYTE rho, _deviceCoeffFunctionPointer fpt)
{
    CLGComplex ret = _deviceStapleT1ShiftedU1(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt);
    ret = _cuCaddf(ret, _deviceStapleT2ShiftedU1(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, mu, nu, rho, fpt));
    ret = _cuCaddf(ret, _deviceStapleT1ShiftedU1(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, rho, nu, mu, fpt));
    ret = _cuCaddf(ret, _deviceStapleT2ShiftedU1(byFieldId, pDeviceData, sSite, uiSiteIndex, uiBigIndex, rho, nu, mu, fpt));
    return ret;
}

#pragma endregion

#pragma endregion

#pragma endregion

//11 sec
static __device__ __inline__ CLGComplex _deviceVXXTauOptimizedU1(
    const CLGComplex* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    UINT bXorY, UBOOL bPlusMu, UBOOL bPlusTau)
{
    const INT byMu = bXorY ? 0 : 1;
    const INT iMu = bPlusMu ? (byMu + 1) : (-byMu - 1);
    const INT iTau = bPlusTau ? 4 : -4;
    INT dir1[3];

    SSmallInt4 x_p_mu_p_tau = sStartSite;
    x_p_mu_p_tau.w = x_p_mu_p_tau.w + (bPlusTau ? 1 : -1);
    //x_p_mu_p_tau.m_byData4[byMu] = x_p_mu_p_tau.m_byData4[byMu] + (bPlusMu ? 1 : -1);
    x_p_mu_p_tau.m_byData4[byMu] = x_p_mu_p_tau.m_byData4[byMu] + (bPlusMu ? 1 : -2);

    dir1[0] = iMu; dir1[1] = iTau;
    CLGComplex sRet = _deviceLinkU1(pDeviceData, sStartSite, 2, byFieldId, dir1);
    dir1[0] = iTau; dir1[1] = iMu;
    CLGComplex toAdd = _deviceLinkU1(pDeviceData, sStartSite, 2, byFieldId, dir1);
    sRet.x += toAdd.x;
    sRet.y += toAdd.y;

    //dir1[0] = iMu;
    //sRet.Mul(_deviceLink(pDeviceData, x_p_mu_p_tau, 1, byFieldId, dir1));

    const SIndex& x_p_taumu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(x_p_mu_p_tau) + byMu];
    toAdd = pDeviceData[_deviceGetLinkIndex(x_p_taumu__mu.m_uiSiteIndex, x_p_taumu__mu.m_byDir)];
    if ((x_p_taumu__mu.NeedToDagger() && bPlusMu)
        || (!x_p_taumu__mu.NeedToDagger() && !bPlusMu))
    {
        toAdd.y = -toAdd.y;
        sRet = _cuCmulf(sRet, toAdd);
    }
    else
    {
        sRet = _cuCmulf(sRet, toAdd);
    }

    dir1[0] = iMu;
    dir1[1] = iMu;
    dir1[2] = iTau;

    toAdd = _deviceLinkU1(pDeviceData, sStartSite, 3, byFieldId, dir1);
    sRet.x = (sRet.x + toAdd.x) * OneOver3;
    sRet.y = (sRet.y + toAdd.y) * OneOver3;
    return sRet;
}

static __device__ __inline__ CLGComplex _deviceVXYTOptimizedU1(
    const CLGComplex* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    UBOOL bPlusX, UBOOL bPlusY, UBOOL bPlusTau)
{
    const INT iX = bPlusX ? 1 : -1;
    const INT iY = bPlusY ? 2 : -2;
    const INT iT = bPlusTau ? 4 : -4;

    SSmallInt4 n_xy = sStartSite;
    SSmallInt4 n_xt = sStartSite;
    SSmallInt4 n_yt = sStartSite;
    if (bPlusX)
    {
        n_xy.x = n_xy.x + 1;
        n_xt.x = n_xt.x + 1;
    }
    else
    {
        n_xy.x = n_xy.x - 1;
        n_xt.x = n_xt.x - 1;
        n_yt.x = n_yt.x - 1;
    }

    if (bPlusY)
    {
        n_xy.y = n_xy.y + 1;
        n_yt.y = n_yt.y + 1;
    }
    else
    {
        n_xy.y = n_xy.y - 1;
        n_yt.y = n_yt.y - 1;
        n_xt.y = n_xt.y - 1;
    }

    if (bPlusTau)
    {
        n_xt.w = n_xt.w + 1;
        n_yt.w = n_yt.w + 1;
    }
    else
    {
        n_xt.w = n_xt.w - 1;
        n_yt.w = n_yt.w - 1;
        n_xy.w = n_xy.w - 1;
    }

    INT dir1[2];

    dir1[0] = iX; dir1[1] = iY;
    CLGComplex sRet1 = _deviceLinkU1(pDeviceData, sStartSite, 2, byFieldId, dir1);

    dir1[0] = iY; dir1[1] = iX;
    CLGComplex toAdd = _deviceLinkU1(pDeviceData, sStartSite, 2, byFieldId, dir1);
    sRet1.x += toAdd.x;
    sRet1.y += toAdd.y;

    //dir1[0] = iT;
    //sRet1.Mul(_deviceLink(pDeviceData, n_xy, 1, byFieldId, dir1));
    const SIndex& n_xy__t = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_xy) + 3];
    toAdd = pDeviceData[_deviceGetLinkIndex(n_xy__t.m_uiSiteIndex, n_xy__t.m_byDir)];
    if ((n_xy__t.NeedToDagger() && bPlusTau)
        || (!n_xy__t.NeedToDagger() && !bPlusTau))
    {
        toAdd.y = -toAdd.y;
        sRet1 = _cuCmulf(sRet1, toAdd);
    }
    else
    {
        sRet1 = _cuCmulf(sRet1, toAdd);
    }

    dir1[0] = iX; dir1[1] = iT;
    CLGComplex sRet2 = _deviceLinkU1(pDeviceData, sStartSite, 2, byFieldId, dir1);

    dir1[0] = iT; dir1[1] = iX;
    toAdd = _deviceLinkU1(pDeviceData, sStartSite, 2, byFieldId, dir1);
    sRet2.x += toAdd.x;
    sRet2.y += toAdd.y;

    //dir1[0] = iY;
    //sRet2.Mul(_deviceLink(pDeviceData, n_xt, 1, byFieldId, dir1));
    const SIndex& n_xt__y = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_xt) + 1];
    toAdd = pDeviceData[_deviceGetLinkIndex(n_xt__y.m_uiSiteIndex, n_xt__y.m_byDir)];
    if ((n_xt__y.NeedToDagger() && bPlusY)
        || (!n_xt__y.NeedToDagger() && !bPlusY))
    {
        toAdd.y = -toAdd.y;
        sRet2 = _cuCmulf(sRet2, toAdd);
    }
    else
    {
        sRet2 = _cuCmulf(sRet2, toAdd);
    }

    sRet1.x += sRet2.x;
    sRet1.y += sRet2.y;

    dir1[0] = iY; dir1[1] = iT;
    sRet2 = _deviceLinkU1(pDeviceData, sStartSite, 2, byFieldId, dir1);

    dir1[0] = iT; dir1[1] = iY;
    toAdd = _deviceLinkU1(pDeviceData, sStartSite, 2, byFieldId, dir1);
    sRet2.x += toAdd.x;
    sRet2.y += toAdd.y;

    //dir1[0] = iX;
    //sRet2.Mul(_deviceLink(pDeviceData, n_yt, 1, byFieldId, dir1));
    const SIndex& n_yt__x = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_yt)];
    toAdd = pDeviceData[_deviceGetLinkIndex(n_yt__x.m_uiSiteIndex, n_yt__x.m_byDir)];
    if ((n_yt__x.NeedToDagger() && bPlusX)
        || (!n_yt__x.NeedToDagger() && !bPlusX))
    {
        toAdd.y = -toAdd.y;
        sRet2 = _cuCmulf(sRet2, toAdd);
    }
    else
    {
        sRet2 = _cuCmulf(sRet2, toAdd);
    }

    sRet1.x = (sRet1.x + sRet2.x) * OneOver6;
    sRet1.y = (sRet1.y + sRet2.y) * OneOver6;
    return sRet1;
}

#pragma region device functions GAMMA

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


#pragma region device functions GAMMA EM


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


#pragma region device functions



/**
 * V(n1,n2)
 */
static __device__ __inline__ deviceSU3 _deviceVXXTauEM(
    const deviceSU3* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataU1,
    const SSmallInt4& sStartSite, Real fCharge, BYTE byFieldId,
    UINT bXorY, UBOOL bPlusMu, UBOOL bPlusTau)
{
    const INT iMu = bXorY ? (bPlusMu ? 1 : -1) : (bPlusMu ? 2 : -2);
    const INT iTau = bPlusTau ? 4 : -4;
    INT dir1[3];

    dir1[0] = iMu;
    dir1[1] = iMu;
    dir1[2] = iTau;
    deviceSU3 sRet = _deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 3, byFieldId, dir1);

    dir1[0] = iMu;
    dir1[1] = iTau;
    dir1[2] = iMu;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iTau;
    dir1[1] = iMu;
    dir1[2] = iMu;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 3, byFieldId, dir1));

    sRet.MulReal(OneOver3);
    return sRet;
}

static __device__ __inline__ deviceSU3 _deviceVXYTEM(
    const deviceSU3* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataU1,
    const SSmallInt4& sStartSite,
    Real fCharge, BYTE byFieldId,
    UBOOL bPlusX, UBOOL bPlusY, UBOOL bPlusTau)
{
    const INT iX = bPlusX ? 1 : -1;
    const INT iY = bPlusY ? 2 : -2;
    const INT iT = bPlusTau ? 4 : -4;
    INT dir1[3];

    dir1[0] = iX; dir1[1] = iY; dir1[2] = iT;
    deviceSU3 sRet(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iX; dir1[1] = iT; dir1[2] = iY;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iY; dir1[1] = iX; dir1[2] = iT;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iY; dir1[1] = iT; dir1[2] = iX;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iT; dir1[1] = iX; dir1[2] = iY;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iT; dir1[1] = iY; dir1[2] = iX;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 3, byFieldId, dir1));

    sRet.MulReal(OneOver6);
    return sRet;
}


//11 sec
static __device__ __inline__ deviceSU3 _deviceVXXTauOptimizedEM(
    const deviceSU3* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataU1,
    const SSmallInt4& sStartSite,
    Real fCharge, BYTE byFieldId,
    UINT bXorY, UBOOL bPlusMu, UBOOL bPlusTau)
{
    const INT byMu = bXorY ? 0 : 1;
    const INT iMu = bPlusMu ? (byMu + 1) : (-byMu - 1);
    const INT iTau = bPlusTau ? 4 : -4;
    INT dir1[3];

    SSmallInt4 x_p_mu_p_tau = sStartSite;
    x_p_mu_p_tau.w = x_p_mu_p_tau.w + (bPlusTau ? 1 : -1);
    //x_p_mu_p_tau.m_byData4[byMu] = x_p_mu_p_tau.m_byData4[byMu] + (bPlusMu ? 1 : -1);
    x_p_mu_p_tau.m_byData4[byMu] = x_p_mu_p_tau.m_byData4[byMu] + (bPlusMu ? 1 : -2);

    //Sum of the left-squre
    dir1[0] = iMu; dir1[1] = iTau;
    deviceSU3 sRet = _deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1);
    dir1[0] = iTau; dir1[1] = iMu;
    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iMu;
    //sRet.Mul(_deviceLink(pDeviceData, x_p_mu_p_tau, 1, byFieldId, dir1));

    const SIndex& x_p_taumu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(x_p_mu_p_tau) + byMu];
    const UBOOL bLastLinkDagger = (x_p_taumu__mu.NeedToDagger() && bPlusMu)
        || (!x_p_taumu__mu.NeedToDagger() && !bPlusMu);
    const UINT x_p_taumu__mu__link = _deviceGetLinkIndex(x_p_taumu__mu.m_uiSiteIndex, x_p_taumu__mu.m_byDir);
    const Real x_p_taumu__mu__phase = pDeviceDataU1[x_p_taumu__mu__link] * fCharge;
    if (bLastLinkDagger)
    {
        sRet.MulDagger(pDeviceData[x_p_taumu__mu__link]);
        sRet.MulComp(_make_cuComplex(_cos(x_p_taumu__mu__phase), -_sin(x_p_taumu__mu__phase)));
    }
    else
    {
        sRet.Mul(pDeviceData[x_p_taumu__mu__link]);
        sRet.MulComp(_make_cuComplex(_cos(x_p_taumu__mu__phase), _sin(x_p_taumu__mu__phase)));
    }

    //The final line
    dir1[0] = iMu;
    dir1[1] = iMu;
    dir1[2] = iTau;

    sRet.Add(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 3, byFieldId, dir1));

    sRet.MulReal(OneOver3);
    return sRet;
}

static __device__ __inline__ deviceSU3 _deviceVXYTOptimizedEM(
    const deviceSU3* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataU1,
    const SSmallInt4& sStartSite,
    Real fCharge, BYTE byFieldId,
    UBOOL bPlusX, UBOOL bPlusY, UBOOL bPlusTau)
{
    const INT iX = bPlusX ? 1 : -1;
    const INT iY = bPlusY ? 2 : -2;
    const INT iT = bPlusTau ? 4 : -4;

    SSmallInt4 n_xy = sStartSite;
    SSmallInt4 n_xt = sStartSite;
    SSmallInt4 n_yt = sStartSite;
    if (bPlusX)
    {
        n_xy.x = n_xy.x + 1;
        n_xt.x = n_xt.x + 1;
    }
    else
    {
        n_xy.x = n_xy.x - 1;
        n_xt.x = n_xt.x - 1;
        n_yt.x = n_yt.x - 1;
    }

    if (bPlusY)
    {
        n_xy.y = n_xy.y + 1;
        n_yt.y = n_yt.y + 1;
    }
    else
    {
        n_xy.y = n_xy.y - 1;
        n_yt.y = n_yt.y - 1;
        n_xt.y = n_xt.y - 1;
    }

    if (bPlusTau)
    {
        n_xt.w = n_xt.w + 1;
        n_yt.w = n_yt.w + 1;
    }
    else
    {
        n_xt.w = n_xt.w - 1;
        n_yt.w = n_yt.w - 1;
        n_xy.w = n_xy.w - 1;
    }

    INT dir1[2];

    dir1[0] = iX; dir1[1] = iY;
    deviceSU3 sRet1(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1));

    dir1[0] = iY; dir1[1] = iX;
    sRet1.Add(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iT;
    //sRet1.Mul(_deviceLink(pDeviceData, n_xy, 1, byFieldId, dir1));
    const SIndex& n_xy__t = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_xy) + 3];
    const UINT n_xy__t_link = _deviceGetLinkIndex(n_xy__t.m_uiSiteIndex, n_xy__t.m_byDir);
    const Real n_xy__t_phase = pDeviceDataU1[n_xy__t_link] * fCharge;
    if ((n_xy__t.NeedToDagger() && bPlusTau)
        || (!n_xy__t.NeedToDagger() && !bPlusTau))
    {
        //at present, we only consider magnetic, so there is no t-dir
        sRet1.MulDagger(pDeviceData[n_xy__t_link]);
        sRet1.MulComp(_make_cuComplex(_cos(n_xy__t_phase), -_sin(n_xy__t_phase)));
    }
    else
    {
        //at present, we only consider magnetic, so there is no t-dir
        sRet1.Mul(pDeviceData[n_xy__t_link]);
        sRet1.MulComp(_make_cuComplex(_cos(n_xy__t_phase), _sin(n_xy__t_phase)));
    }

    dir1[0] = iX; dir1[1] = iT;
    deviceSU3 sRet2 = _deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1);

    dir1[0] = iT; dir1[1] = iX;
    sRet2.Add(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iY;
    //sRet2.Mul(_deviceLink(pDeviceData, n_xt, 1, byFieldId, dir1));

    //================== add magnetic phase y ========================
    const SIndex& n_xt__y = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_xt) + 1];
    const UINT n_xt__y_link = _deviceGetLinkIndex(n_xt__y.m_uiSiteIndex, n_xt__y.m_byDir);
    const Real n_xt__y_phase = pDeviceDataU1[n_xt__y_link] * fCharge;
    if ((n_xt__y.NeedToDagger() && bPlusY)
        || (!n_xt__y.NeedToDagger() && !bPlusY))
    {
        sRet2.MulDagger(pDeviceData[n_xt__y_link]);
        sRet2.MulComp(_make_cuComplex(_cos(n_xt__y_phase), -_sin(n_xt__y_phase)));
    }
    else
    {
        sRet2.Mul(pDeviceData[n_xt__y_link]);
        sRet2.MulComp(_make_cuComplex(_cos(n_xt__y_phase), _sin(n_xt__y_phase)));
    }

    sRet1.Add(sRet2);

    dir1[0] = iY; dir1[1] = iT;
    sRet2 = _deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1);

    dir1[0] = iT; dir1[1] = iY;
    sRet2.Add(_deviceLinkEM(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iX;
    //sRet2.Mul(_deviceLink(pDeviceData, n_yt, 1, byFieldId, dir1));

    //================== add magnetic phase x ========================
    const SIndex& n_yt__x = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_yt)];
    const UINT n_yt__x_link = _deviceGetLinkIndex(n_yt__x.m_uiSiteIndex, n_yt__x.m_byDir);
    const Real n_yt__x_phase = pDeviceDataU1[n_yt__x_link] * fCharge;
    if ((n_yt__x.NeedToDagger() && bPlusX)
        || (!n_yt__x.NeedToDagger() && !bPlusX))
    {
        sRet2.MulDagger(pDeviceData[n_yt__x_link]);
        sRet2.MulComp(_make_cuComplex(_cos(n_yt__x_phase), -_sin(n_yt__x_phase)));
    }
    else
    {
        sRet2.Mul(pDeviceData[n_yt__x_link]);
        sRet2.MulComp(_make_cuComplex(_cos(n_yt__x_phase), _sin(n_yt__x_phase)));
    }

    sRet1.Add(sRet2);

    sRet1.MulReal(OneOver6);
    return sRet1;
}

#pragma endregion


#pragma region device functions KS SU3 Rotation



/**
 * V(n1,n2)
 */
static __device__ __inline__ deviceSU3 _deviceVXXTau(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    UINT bXorY, UBOOL bPlusMu, UBOOL bPlusTau)
{
    const INT iMu = bXorY ? (bPlusMu ? 1 : -1) : (bPlusMu ? 2 : -2);
    const INT iTau = bPlusTau ? 4 : -4;
    INT dir1[3];

    dir1[0] = iMu;
    dir1[1] = iMu;
    dir1[2] = iTau;
    deviceSU3 sRet = _deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1);

    dir1[0] = iMu;
    dir1[1] = iTau;
    dir1[2] = iMu;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iTau;
    dir1[1] = iMu;
    dir1[2] = iMu;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    sRet.MulReal(OneOver3);
    return sRet;
}

static __device__ __inline__ Real _deviceEta124(const SSmallInt4& sSite)
{
    return (((sSite.y + sSite.z) & 1) > 0) ? (F(-1.0)) : (F(1.0));
}

static __device__ __inline__ deviceSU3 _deviceVXYT(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    UBOOL bPlusX, UBOOL bPlusY, UBOOL bPlusTau)
{
    const INT iX = bPlusX ? 1 : -1;
    const INT iY = bPlusY ? 2 : -2;
    const INT iT = bPlusTau ? 4 : -4;
    INT dir1[3];

    dir1[0] = iX; dir1[1] = iY; dir1[2] = iT;
    deviceSU3 sRet(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iX; dir1[1] = iT; dir1[2] = iY;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iY; dir1[1] = iX; dir1[2] = iT;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iY; dir1[1] = iT; dir1[2] = iX;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iT; dir1[1] = iX; dir1[2] = iY;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iT; dir1[1] = iY; dir1[2] = iX;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    sRet.MulReal(OneOver6);
    return sRet;
}


/*
static __device__ __inline__ deviceSU3 _deviceXXT_PP(
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byFieldId, BYTE byMu)
{
    const SIndex& n__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byMu];
    _deviceSmallInt4Offset(sStartSite, __fwd(byMu));
    const SIndex& n_p_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byMu];
    _deviceSmallInt4Offset(sStartSite, __fwd(byMu));
    const SIndex& n_p_2mu__t = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + 3];

    deviceSU3 ret = pDeviceData[_deviceGetLinkIndex(n__mu.m_uiSiteIndex, n__mu.m_byDir)];
    //n__mu,n_p_mu__mu,n_p_2mu__t is never daggered
    ret.Mul(pDeviceData[_deviceGetLinkIndex(n_p_mu__mu.m_uiSiteIndex, n_p_mu__mu.m_byDir)]);
    ret.Mul(pDeviceData[_deviceGetLinkIndex(n_p_2mu__t.m_uiSiteIndex, n_p_2mu__t.m_byDir)]);

    return ret;
}

static __device__ __inline__ deviceSU3 _deviceXXT_PM(
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byFieldId, BYTE byMu)
{
    const SIndex& n__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byMu];
    _deviceSmallInt4Offset(sStartSite, __fwd(byMu));
    const SIndex& n_p_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byMu];
    _deviceSmallInt4Offset(sStartSite, __fwd(byMu));
    _deviceSmallInt4Offset(sStartSite, -4);
    const SIndex& n_p_2mu_m_t__t = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + 3];

    deviceSU3 ret = pDeviceData[_deviceGetLinkIndex(n__mu.m_uiSiteIndex, n__mu.m_byDir)];
    //n__mu is never daggered
    ret.Mul(pDeviceData[_deviceGetLinkIndex(n_p_mu__mu.m_uiSiteIndex, n_p_mu__mu.m_byDir)]);
    ret.MulDagger(pDeviceData[_deviceGetLinkIndex(n_p_2mu_m_t__t.m_uiSiteIndex, n_p_2mu_m_t__t.m_byDir)]);

    return ret;
}

static __device__ __inline__ deviceSU3 _deviceXXT_MP(
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byFieldId, BYTE byMu)
{
    _deviceSmallInt4Offset(sStartSite, __bck(byMu));
    const SIndex& n_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byMu];
    _deviceSmallInt4Offset(sStartSite, __bck(byMu));
    const SIndex& n_m_2mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byMu];
    const SIndex& n_m_2mu__t = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + 3];

    deviceSU3 ret = pDeviceData[_deviceGetLinkIndex(n_m_2mu__mu.m_uiSiteIndex, n_m_2mu__mu.m_byDir)];
    ret.Mul(pDeviceData[_deviceGetLinkIndex(n_m_mu__mu.m_uiSiteIndex, n_m_mu__mu.m_byDir)]);
    ret.DaggerMul(pDeviceData[_deviceGetLinkIndex(n_m_2mu__t.m_uiSiteIndex, n_m_2mu__t.m_byDir)]);
    return ret;
}

static __device__ __inline__ deviceSU3 _deviceXXT_MM(
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byFieldId, BYTE byMu)
{
    _deviceSmallInt4Offset(sStartSite, __bck(byMu));
    const SIndex& n_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byMu];
    _deviceSmallInt4Offset(sStartSite, __bck(byMu));
    const SIndex& n_m_2mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byMu];
    _deviceSmallInt4Offset(sStartSite, -4);
    const SIndex& n_m_2mu_mu_t__t = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + 3];

    deviceSU3 ret = pDeviceData[_deviceGetLinkIndex(n_m_2mu_mu_t__t.m_uiSiteIndex, n_m_2mu_mu_t__t.m_byDir)];
    ret.Mul(pDeviceData[_deviceGetLinkIndex(n_m_2mu__mu.m_uiSiteIndex, n_m_2mu__mu.m_byDir)]);
    ret.Mul(pDeviceData[_deviceGetLinkIndex(n_m_mu__mu.m_uiSiteIndex, n_m_mu__mu.m_byDir)]);
    ret.Dagger();
    return ret;
}
*/

//11 sec
/**
* byMu = bXorY ? 0 : 1
* so that, when bXorY = 1, it is partial_X
*          when bXorY = 0, it is partial_Y
*/
static __device__ __inline__ deviceSU3 _deviceVXXTauOptimized(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    UINT bXorY, UBOOL bPlusMu, UBOOL bPlusTau)
{
    const INT byMu = bXorY ? 0 : 1;
    const INT iMu = bPlusMu ? (byMu + 1) : (-byMu - 1);
    const INT iTau = bPlusTau ? 4 : -4;
    INT dir1[3];

    SSmallInt4 x_p_mu_p_tau = sStartSite;
    x_p_mu_p_tau.w = x_p_mu_p_tau.w + (bPlusTau ? 1 : -1);
    //x_p_mu_p_tau.m_byData4[byMu] = x_p_mu_p_tau.m_byData4[byMu] + (bPlusMu ? 1 : -1);
    x_p_mu_p_tau.m_byData4[byMu] = x_p_mu_p_tau.m_byData4[byMu] + (bPlusMu ? 1 : -2);

    dir1[0] = iMu; dir1[1] = iTau;
    deviceSU3 sRet = _deviceLink(pDeviceData, sStartSite, 2, byFieldId, dir1);
    dir1[0] = iTau; dir1[1] = iMu;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iMu;
    //sRet.Mul(_deviceLink(pDeviceData, x_p_mu_p_tau, 1, byFieldId, dir1));

    const SIndex& x_p_taumu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(x_p_mu_p_tau) + byMu];
    if ((x_p_taumu__mu.NeedToDagger() && bPlusMu)
        || (!x_p_taumu__mu.NeedToDagger() && !bPlusMu))
    {
        sRet.MulDagger(pDeviceData[_deviceGetLinkIndex(x_p_taumu__mu.m_uiSiteIndex, x_p_taumu__mu.m_byDir)]);
    }
    else
    {
        sRet.Mul(pDeviceData[_deviceGetLinkIndex(x_p_taumu__mu.m_uiSiteIndex, x_p_taumu__mu.m_byDir)]);
    }

    dir1[0] = iMu;
    dir1[1] = iMu;
    dir1[2] = iTau;

    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    sRet.MulReal(OneOver3);
    return sRet;
}

static __device__ __inline__ deviceSU3 _deviceVXYTOptimized(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    UBOOL bPlusX, UBOOL bPlusY, UBOOL bPlusTau)
{
    const INT iX = bPlusX ? 1 : -1;
    const INT iY = bPlusY ? 2 : -2;
    const INT iT = bPlusTau ? 4 : -4;

    SSmallInt4 n_xy = sStartSite;
    SSmallInt4 n_xt = sStartSite;
    SSmallInt4 n_yt = sStartSite;
    if (bPlusX)
    {
        n_xy.x = n_xy.x + 1;
        n_xt.x = n_xt.x + 1;
    }
    else
    {
        n_xy.x = n_xy.x - 1;
        n_xt.x = n_xt.x - 1;
        n_yt.x = n_yt.x - 1;
    }

    if (bPlusY)
    {
        n_xy.y = n_xy.y + 1;
        n_yt.y = n_yt.y + 1;
    }
    else
    {
        n_xy.y = n_xy.y - 1;
        n_yt.y = n_yt.y - 1;
        n_xt.y = n_xt.y - 1;
    }

    if (bPlusTau)
    {
        n_xt.w = n_xt.w + 1;
        n_yt.w = n_yt.w + 1;
    }
    else
    {
        n_xt.w = n_xt.w - 1;
        n_yt.w = n_yt.w - 1;
        n_xy.w = n_xy.w - 1;
    }

    INT dir1[2];

    dir1[0] = iX; dir1[1] = iY;
    deviceSU3 sRet1(_deviceLink(pDeviceData, sStartSite, 2, byFieldId, dir1));

    dir1[0] = iY; dir1[1] = iX;
    sRet1.Add(_deviceLink(pDeviceData, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iT;
    //sRet1.Mul(_deviceLink(pDeviceData, n_xy, 1, byFieldId, dir1));
    const SIndex& n_xy__t = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_xy) + 3];
    if ((n_xy__t.NeedToDagger() && bPlusTau)
        || (!n_xy__t.NeedToDagger() && !bPlusTau))
    {
        sRet1.MulDagger(pDeviceData[_deviceGetLinkIndex(n_xy__t.m_uiSiteIndex, n_xy__t.m_byDir)]);
    }
    else
    {
        sRet1.Mul(pDeviceData[_deviceGetLinkIndex(n_xy__t.m_uiSiteIndex, n_xy__t.m_byDir)]);
    }

    dir1[0] = iX; dir1[1] = iT;
    deviceSU3 sRet2 = _deviceLink(pDeviceData, sStartSite, 2, byFieldId, dir1);

    dir1[0] = iT; dir1[1] = iX;
    sRet2.Add(_deviceLink(pDeviceData, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iY;
    //sRet2.Mul(_deviceLink(pDeviceData, n_xt, 1, byFieldId, dir1));
    const SIndex& n_xt__y = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_xt) + 1];
    if ((n_xt__y.NeedToDagger() && bPlusY)
        || (!n_xt__y.NeedToDagger() && !bPlusY))
    {
        sRet2.MulDagger(pDeviceData[_deviceGetLinkIndex(n_xt__y.m_uiSiteIndex, n_xt__y.m_byDir)]);
    }
    else
    {
        sRet2.Mul(pDeviceData[_deviceGetLinkIndex(n_xt__y.m_uiSiteIndex, n_xt__y.m_byDir)]);
    }

    sRet1.Add(sRet2);

    dir1[0] = iY; dir1[1] = iT;
    sRet2 = _deviceLink(pDeviceData, sStartSite, 2, byFieldId, dir1);

    dir1[0] = iT; dir1[1] = iY;
    sRet2.Add(_deviceLink(pDeviceData, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iX;
    //sRet2.Mul(_deviceLink(pDeviceData, n_yt, 1, byFieldId, dir1));
    const SIndex& n_yt__x = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_yt)];
    if ((n_yt__x.NeedToDagger() && bPlusX)
        || (!n_yt__x.NeedToDagger() && !bPlusX))
    {
        sRet2.MulDagger(pDeviceData[_deviceGetLinkIndex(n_yt__x.m_uiSiteIndex, n_yt__x.m_byDir)]);
    }
    else
    {
        sRet2.Mul(pDeviceData[_deviceGetLinkIndex(n_yt__x.m_uiSiteIndex, n_yt__x.m_byDir)]);
    }

    sRet1.Add(sRet2);

    sRet1.MulReal(OneOver6);
    return sRet1;
}

#pragma endregion

#pragma region KS SU3 Acc

/**
* gamma_mu partial_nu term
* used for test, when test is done, we use _deviceGMUPNUOptimized
*/
static __device__ __inline__ deviceSU3 _deviceGMUPNU(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    BYTE mu, BYTE nu, UBOOL bPlusMu, UBOOL bPlusNu)
{
    const INT iMu = bPlusMu ? (mu + 1) : (-mu - 1);
    const INT iNu = bPlusNu ? (nu + 1) : (-nu - 1);
    INT dir1[3];

    dir1[0] = iNu;
    dir1[1] = iNu;
    dir1[2] = iMu;
    deviceSU3 sRet = _deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1);

    dir1[0] = iNu;
    dir1[1] = iMu;
    dir1[2] = iNu;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iMu;
    dir1[1] = iNu;
    dir1[2] = iNu;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    sRet.MulReal(OneOver3);
    return sRet;
}

/**
* gamma_mu partial_nu term
*
*
* ------------------------> x dir is nu
*
*    +--o   +--+--o        o               o      +--o      +--+--o
*    |      |              |    =      (   |    + |   )  +  |
* x--+      x        x--+--+        x-- +--+      +         x
*
*
* x--+      x        x--+--+
*    |      |              |
*    +--o   +--+--o        o
*
*
*
*/
static __device__ __inline__ deviceSU3 _deviceGMUPNUOptimized(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    BYTE mu, BYTE nu, UBOOL bPlusMu, UBOOL bPlusNu)
{
    const INT iMu = bPlusMu ? (mu + 1) : (-mu - 1);
    const INT iNu = bPlusNu ? (nu + 1) : (-nu - 1);
    INT dir1[3];

    SSmallInt4 x_p_nu_p_mu = sStartSite;
    x_p_nu_p_mu.m_byData4[mu] = x_p_nu_p_mu.m_byData4[mu] + (bPlusMu ? 1 : -1);
    x_p_nu_p_mu.m_byData4[nu] = x_p_nu_p_mu.m_byData4[nu] + (bPlusNu ? 1 : -2);

    dir1[0] = iMu; dir1[1] = iNu;
    deviceSU3 sRet = _deviceLink(pDeviceData, sStartSite, 2, byFieldId, dir1);
    dir1[0] = iNu; dir1[1] = iMu;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 2, byFieldId, dir1));

    const SIndex& x_p_taumu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(x_p_nu_p_mu) + nu];
    if ((x_p_taumu__mu.NeedToDagger() && bPlusNu)
        || (!x_p_taumu__mu.NeedToDagger() && !bPlusNu))
    {
        sRet.MulDagger(pDeviceData[_deviceGetLinkIndex(x_p_taumu__mu.m_uiSiteIndex, x_p_taumu__mu.m_byDir)]);
    }
    else
    {
        sRet.Mul(pDeviceData[_deviceGetLinkIndex(x_p_taumu__mu.m_uiSiteIndex, x_p_taumu__mu.m_byDir)]);
    }

    dir1[0] = iNu;
    dir1[1] = iNu;
    dir1[2] = iMu;

    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    sRet.MulReal(OneOver3);
    return sRet;
}

#pragma endregion


#pragma region Device functions Gauge Plaquette Boost

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

#pragma region device functions Wilson Dirac Dirichlet

static __device__ __inline__ deviceWilsonVectorSU3 _deviceGetFermionBCWilsonSU3(
    const deviceWilsonVectorSU3* __restrict__ pBuffer,
    const SIndex& idx,
    BYTE byFieldId)
{
    return idx.IsDirichlet() ?
        deviceWilsonVectorSU3::makeZeroWilsonVectorSU3()
        : pBuffer[idx.m_uiSiteIndex];
}

static __device__ __inline__ const deviceWilsonVectorSU3& _deviceGetFermionBCWilsonSU3T(
    const deviceWilsonVectorSU3* __restrict__ pBuffer,
    const SIndex& idx,
    BYTE byFieldId)
{
    return pBuffer[idx.m_uiSiteIndex];
}

#pragma endregion

#pragma region device functions Measure Staggered Meson Simple

static __device__ __inline__ SBYTE _deviceStaggeredFermionSimplePhase(const SSmallInt4& sSite, BYTE byType)
{
    SBYTE ret = 1;
    switch (byType)
    {
    case 1:
        ret = 3 - ((sSite.x & 1) << 1)
            - ((sSite.y & 1) << 1)
            - ((sSite.z & 1) << 1);
        //printf("shift check%d = %d\n", static_cast<INT>(ret), ((sSite.x & 1) ? -1 : 1) + ((sSite.y & 1) ? -1 : 1) + ((sSite.z & 1) ? -1 : 1));
        break;
    case 2:
        ret = 3 - (((sSite.x + sSite.y) & 1) << 1)
            - (((sSite.y + sSite.z) & 1) << 1)
            - (((sSite.x + sSite.z) & 1) << 1);
        break;
    case 3:
        ret = 1 - (((sSite.x + sSite.y + sSite.z) & 1) << 1);
        break;
    case 4:
        ret = 1 - ((sSite.x & 1) << 1);
        break;
    case 5:
        ret = 1 - ((sSite.y & 1) << 1);
        break;
    case 6:
        ret = 1 - ((sSite.z & 1) << 1);
        break;
    case 7:
        ret = 1 - (((sSite.x + sSite.y) & 1) << 1);
        break;
    case 8:
        ret = 1 - (((sSite.y + sSite.z) & 1) << 1);
        break;
    case 9:
        ret = 1 - (((sSite.x + sSite.z) & 1) << 1);
        break;
    default:
        break;
    }
    return ret;
}

#pragma endregion

#pragma region device functions Measure Topological charge XY

static __device__ __inline__ Real _deviceTrImClover(const deviceSU3* __restrict__ pGaugeField, BYTE byFieldId, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE rho, BYTE sigma)
{
    return deviceSU3::TrIm(
        _deviceClover(pGaugeField, sSite4, uiBigIdx, mu, nu, byFieldId),
        _deviceClover(pGaugeField, sSite4, uiBigIdx, rho, sigma, byFieldId));
}

static __device__ __inline__ Real _deviceTopologicalCharge(const deviceSU3* __restrict__ pGaugeField, BYTE byFieldId, const SSmallInt4& sSite4, UINT uiBigIdx)
{
    Real ret = _deviceTrImClover(pGaugeField, byFieldId, sSite4, uiBigIdx, 0, 1, 2, 3);
    ret -= _deviceTrImClover(pGaugeField, byFieldId, sSite4, uiBigIdx, 0, 2, 1, 3);
    ret += _deviceTrImClover(pGaugeField, byFieldId, sSite4, uiBigIdx, 0, 3, 1, 2);
    return OneOver32PI2 * F(2.0) * ret;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _DEVICEINLINESU3_H_

//=============================================================================
// END OF FILE
//=============================================================================
