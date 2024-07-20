//=============================================================================
// FILENAME : DeviceInlineGauge.h
// 
// DESCRIPTION:
// This should be implemented using inherint machinism, but due to historical reasons, it is now templates
//
//
// REVISION:
//  [07/03/2024 nbale]
//=============================================================================

#ifndef _DEVICEINLINEGAUGE_H_
#define _DEVICEINLINEGAUGE_H_

__BEGIN_NAMESPACE

#pragma region Gauge

/**
* Note: for baked plaqutte index, the bond if is set to SIndex
* If it is a "new SIndex" instead, remember to set the m_byTag
*/
template<typename deviceGauge>
static __device__ __inline__ const deviceGauge& _deviceGetGaugeBCT(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pBuffer,
    const SIndex& idx)
{
    return idx.IsDirichlet() ?
        ((CFieldBoundary<deviceGauge>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

/**
* If the bond is on surface, return the Dirichlet
* else, return the element
*/
template<typename deviceGauge>
static __device__ __inline__ const deviceGauge& _deviceGetGaugeBCDirT(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, byDir) ?
        ((CFieldBoundary<deviceGauge>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(site) * _DC_Dir + byDir
        ]
        : pBuffer[_deviceGetLinkIndex(site.m_uiSiteIndex, byDir)];
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceGetGaugeBCDirOneT(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, byDir) ?
        _makeId<deviceGauge>()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceGetGaugeBCDirZeroT(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, byDir) ?
        _makeZero<deviceGauge>()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceGetGaugeBCDirSIndexT(
    const deviceGauge* __restrict__ pBuffer,
    const SIndex& idx,
    BYTE byFieldId)
{
    deviceGauge ret = idx.IsDirichlet() ?
        ((CFieldBoundary<deviceGauge>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
    if (idx.NeedToDagger())
    {
        _dagger(ret);
    }
    return ret;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceGetGaugeBCDirOneSIndexT(
    const deviceGauge* __restrict__ pBuffer,
    const SIndex& idx)
{
    if (idx.IsDirichlet())
    {
        return _makeId<deviceGauge>();
    }
    if (idx.NeedToDagger())
    {
        return _daggerC(pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)]);
    }

    return pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

/**
 * Note that, when get zero instead of one, it is minus not dagger
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceGetGaugeBCDirZeroSIndexT(
    const deviceGauge* __restrict__ pBuffer,
    const SIndex& idx)
{
    if (idx.IsDirichlet())
    {
        return _makeZero<deviceGauge>();
    }
    if (idx.NeedToDagger())
    {
        return _mulC(pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)], F(-1.0));
    }

    return pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

/**
 * calculate D_mu A _nu = Delta _mu + [A_mu, A _nu]
 * Use U now to calculate A pure
 * me will be changed, so, if me is A phys, copy me first
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceDPureMuT(
    const deviceGauge* __restrict__ piA,
    const deviceGauge* __restrict__ piApure,
    const SSmallInt4& sSite4,
    UINT uiBigIdx,
    BYTE byMu,
    BYTE byNu,
    BYTE byFieldId)
{
    //i a D A = (A_nu (n) - A_nu (n-mu)) + iApure _mu A _nu - i A _nu Apure _mu
    const UINT uiSiteBig_m_mu = __idx->_deviceGetBigIndex(_deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1));

    deviceGauge res = _deviceGetGaugeBCDirZeroT(byFieldId, piApure, uiBigIdx, byMu); //Apure _mu
    deviceGauge res2 = _deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu); //A _nu
    _mul(res2, res); //A _nu Apure _mu
    _mul(res, _deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu)); //Apure _mu A _nu
    _sub(res, res2); //[Apure, A]
    _add(res, _deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu));
    _sub(res, _deviceGetGaugeBCDirZeroSIndexT(piA,__idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_m_mu * _DC_Dir + byNu]));
    return res;
}

/**
 * test using (A(N+mu)-A(N-mu))/2
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceDPureMu2T(
    const deviceGauge* __restrict__ piA,
    const deviceGauge* __restrict__ piApure,
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

    deviceGauge res = _deviceGetGaugeBCDirZeroT(byFieldId, piApure, uiBigIdx, byMu); //Apure _mu
    deviceGauge res2 = _deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu); //A _nu
    _mul(res2, res); //A _nu Apure _mu
    _mul(res, _deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu)); //Apure _mu A _nu
    _sub(res, res2); //[Apure, A]
    _add(res, _mulC(_deviceGetGaugeBCDirZeroSIndexT(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_p_mu * _DC_Dir + byNu]), F(0.5)));
    _sub(res, _mulC(_deviceGetGaugeBCDirZeroSIndexT(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_m_mu * _DC_Dir + byNu]), F(0.5)));
    return res;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _devicePlaqutteT(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
    UINT uiSiteIndex,
    BYTE plaqIdx, //0-5, as 12, 13, 14, 23, 24, 34
    BYTE plaqLength, //Always 4
    BYTE plaqCountAll //Always 24
)
{
    SIndex first = pCachedPlaqutte[plaqIdx * plaqLength + uiSiteIndex * plaqCountAll];
    deviceGauge toAdd(_deviceGetGaugeBCDirOneSIndexT(pDeviceData, first));
    if (first.NeedToDagger())
    {
        _dagger(toAdd);
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[plaqIdx * plaqLength + j + uiSiteIndex * plaqCountAll];
        deviceGauge toMul(_deviceGetGaugeBCDirOneSIndexT(pDeviceData, first));
        if (first.NeedToDagger())
        {
            _muldag(toAdd, toMul);
        }
        else
        {
            _mul(toAdd, toMul);
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
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceLinkT(
    const deviceGauge* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceGauge sRet = _makeId<deviceGauge>();
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
                    _dagger(sRet);
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
                    _muldag(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    _mul(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
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
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceLinkLongT(
    const deviceGauge* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceGauge sRet = _makeId<deviceGauge>();
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
                    _dagger(sRet);
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
                    _muldag(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    _mul(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
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
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _device1PlaqutteTermPPT(
    const deviceGauge* __restrict__ pDeviceData,
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

    deviceGauge u = _deviceGetGaugeBCDirSIndexT(pDeviceData, s_mu, byFieldId);
    _mul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_p_mu_nu, byFieldId));
    _muldag(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_p_nu_mu, byFieldId));
    _muldag(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_nu, byFieldId));

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
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _device1PlaqutteTermMPT(
    const deviceGauge* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1);
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, byNu + 1);
    const UINT uin_m_mub4 = __idx->_deviceGetBigIndex(n_m_mu) * _DC_Dir;
    const SIndex& s_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byMu];
    const SIndex& s_m_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byNu];
    const SIndex& s_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_mu_p_nu) * _DC_Dir + byMu];
    const SIndex& s__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byNu];

    deviceGauge u = _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_mu__mu, byFieldId);
    _dagmul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_mu__nu, byFieldId));
    _mul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_mu_p_nu__mu, byFieldId));
    _muldag(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s__nu, byFieldId));

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
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _device1PlaqutteTermPMT(
    const deviceGauge* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byNu) - 1);
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(n_m_nu, byMu + 1);
    const UINT uin_m_nub4 = __idx->_deviceGetBigIndex(n_m_nu) * _DC_Dir;
    const SIndex& s_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byMu];
    const SIndex& s_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byNu];
    const SIndex& s_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_nu_p_mu) * _DC_Dir + byNu];
    const SIndex& s__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byMu];

    deviceGauge u = _deviceGetGaugeBCDirSIndexT(pDeviceData, s__mu, byFieldId);
    _muldag(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_nu_p_mu__nu, byFieldId));
    _muldag(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_nu__mu, byFieldId));
    _mul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_nu__nu, byFieldId));

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
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _device1PlaqutteTermMMT(
    const deviceGauge* __restrict__ pDeviceData,
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
    deviceGauge u = _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_nu_m_mu__nu, byFieldId);
    _mul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_mu__mu, byFieldId));
    _dagmul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_nu_m_mu__mu, byFieldId));
    _mul(u, _deviceGetGaugeBCDirSIndexT(pDeviceData, s_m_nu__nu, byFieldId));

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
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceCloverT(const deviceGauge* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    deviceGauge ret = _device1PlaqutteTermPPT(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId);
    _add(ret, _device1PlaqutteTermMMT(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId));
    _add(ret, _device1PlaqutteTermPMT(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));
    _add(ret, _device1PlaqutteTermMPT(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));

    return ret;
}

/**
 * Avoid the add of matrices
 */
template<typename deviceGauge>
static __device__ __inline__ Real _deviceCloverRetrT(const deviceGauge* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    return _retr(_device1PlaqutteTermPPT(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId))
         + _retr(_device1PlaqutteTermMMT(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId))
         + _retr(_device1PlaqutteTermPMT(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId))
         + _retr(_device1PlaqutteTermMPT(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _DEVICEINLINEGAUGE_H_

//=============================================================================
// END OF FILE
//=============================================================================
