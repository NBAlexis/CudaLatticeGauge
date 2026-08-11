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

#pragma region Dirichlet site

//No need for this, since the sites are eazier to map
template<typename deviceVector>
static __device__ __inline__ const deviceVector& _deviceGetVectorBCT(
    BYTE byFieldId,
    const deviceVector* __restrict__ pBuffer,
    const SIndex& idx)
{
    return idx.IsDirichlet() ?
        ((CFieldBoundary<deviceVector>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx)
        ]
        : pBuffer[idx.m_uiSiteIndex];
}

template<typename deviceVector>
static __device__ __inline__ deviceVector _deviceGetVectorBCZeroT(
    const deviceVector* __restrict__ pBuffer,
    const SIndex& idx)
{
    return idx.IsDirichlet() ? _makeZero<deviceVector>() : pBuffer[idx.m_uiSiteIndex];
}

#pragma endregion

#pragma region Gauge

/**
* Note: for baked plaqutte index, the bond if is set to SIndex
* If it is a "new SIndex" instead, remember to set the m_byTag
* Note: in old versions, _deviceGetGaugeBCT performs dagger if the SIndex need to dagger
*       however, in the furture, it dose NOT
*       do the dagger by your self!
*/
template<typename deviceGauge>
static __device__ __inline__ const deviceGauge& _deviceGetGaugeBCT(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pBuffer,
    const SIndex& idx)
{
    return idx.IsDirichlet() ?
        ((CFieldBoundary<deviceGauge>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
#if _CLG_ASSUME_SQUARE_LATTICE
            (__idx->_devcieExchangeBoundaryFieldSiteIndex(idx) << 2U) | idx.m_byDir
#else
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
#endif
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

template<typename deviceGauge>
static __device__ __inline__ const deviceGauge* _deviceGetGaugeBCTPTR(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pBuffer,
    const SIndex& idx)
{
    return idx.IsDirichlet() ?
        (((CFieldBoundary<deviceGauge>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData + 
#if _CLG_ASSUME_SQUARE_LATTICE
            ((__idx->_devcieExchangeBoundaryFieldSiteIndex(idx) << _DC_Dir) | idx.m_byDir)
#else
            (__idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir)
#endif
            )
        : (pBuffer + _deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir));
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
 * Note! Note! If you want to get a gauge force, DO NOT use this!!!
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceGetGaugeBCDirAvector_SIndexT(
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

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceGetGaugeBCDirForceSIndexT(
    const deviceGauge* __restrict__ pBuffer,
    const SIndex& idx)
{
    if (idx.IsDirichlet())
    {
        return _makeZero<deviceGauge>();
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
    _sub(res, _deviceGetGaugeBCDirAvector_SIndexT(piA,__idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_m_mu * _DC_Dir + byNu]));
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
#if _CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqIdx
#else
    BYTE plaqIdx, //0-5, as 12, 13, 14, 23, 24, 34
    BYTE plaqLength, //Always 4
    BYTE plaqCountAllSite //Always 24
#endif
)
{
    SIndex first = pCachedPlaqutte[plaqIdx * plaqLength + uiSiteIndex * plaqCountAllSite];
    deviceGauge toAdd(_deviceGetGaugeBCDirOneSIndexT(pDeviceData, first));
    if (first.NeedToDagger())
    {
        _dagger(toAdd);
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[plaqIdx * plaqLength + j + uiSiteIndex * plaqCountAllSite];
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
 * Assume there is no "0" in the path, since it is meaningless
 * but allow "0" for byLength
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceLinkT(
    const deviceGauge* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const SCHAR* __restrict__ pDir)
{
    if (0 == byLength)
    {
        return _makeId<deviceGauge>();
    }

    //+ is 0, - is -1
    SCHAR sign_mask = (pDir[0] >> (sizeof(SCHAR) * 8 - 1));
    SCHAR diridx = sign_mask ^ (pDir[0] - sign_mask - 1);
    sStartSite.m_byData4[diridx] += sign_mask;
#if _CLG_ASSUME_SQUARE_LATTICE
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
    deviceGauge sRet = newLink.IsDirichlet() ? (_makeId<deviceGauge>()) : pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];

    //BYTE debug = 0;
    //if (0 == newLink.m_uiSiteIndex && 1 == pDir[0] && 3 == pDir[1])
    //{
    //    debug = 1;
    //}
    //if (debug)
    //{
    //    printf("0 link (%d)(%d, %d, %d, %d)_ %d(%d)\n", newLink.m_uiSiteIndex, sStartSite.x, sStartSite.y, sStartSite.z, sStartSite.w, newLink.m_byDir, newLink.IsDirichlet());
    //}
    if (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask))
    {
        _dagger(sRet);
    }
    sStartSite.m_byData4[diridx] += (sign_mask + 1);
    
    for (BYTE i = 1U; i < byLength; ++i)
    {
        sign_mask = (pDir[i] >> (sizeof(SCHAR) * 8 - 1));
        diridx = sign_mask ^ (pDir[i] - sign_mask - 1);
        sStartSite.m_byData4[diridx] += sign_mask;
#if _CLG_ASSUME_SQUARE_LATTICE
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
        //if (debug)
        //{
        //    printf("%d link (%d)(%d, %d, %d, %d)_ %d(%d)\n", i, newLink.m_uiSiteIndex, sStartSite.x, sStartSite.y, sStartSite.z, sStartSite.w, newLink.m_byDir, newLink.IsDirichlet());
        //}
        if (newLink.IsDirichlet())
        {
            sStartSite.m_byData4[diridx] += (sign_mask + 1);
            continue;
        }
        if (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask))
        {
            _muldag(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
        }
        else
        {
            _mul(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
        }
        sStartSite.m_byData4[diridx] += (sign_mask + 1);
    }

    return sRet;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceLinkTTwoField(
    const deviceGauge* __restrict__ pDeviceData1,
    const deviceGauge* __restrict__ pDeviceData2,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId, BYTE replaceIndex,
    const SCHAR* __restrict__ pDir)
{
    if (0 == byLength)
    {
        return _makeId<deviceGauge>();
    }

    const deviceGauge* pDeviceDatas[2] = { pDeviceData1, pDeviceData2 };

    //+ is 0, - is -1
    SCHAR sign_mask = (pDir[0] >> (sizeof(SCHAR) * 8 - 1));
    SCHAR diridx = sign_mask ^ (pDir[0] - sign_mask - 1);
    sStartSite.m_byData4[diridx] += sign_mask;
#if _CLG_ASSUME_SQUARE_LATTICE
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
    deviceGauge sRet = newLink.IsDirichlet() ? (_makeId<deviceGauge>()) : pDeviceDatas[0U == replaceIndex][_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];

    if (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask))
    {
        _dagger(sRet);
    }
    sStartSite.m_byData4[diridx] += (sign_mask + 1);

    for (BYTE i = 1U; i < byLength; ++i)
    {
        sign_mask = (pDir[i] >> (sizeof(SCHAR) * 8 - 1));
        diridx = sign_mask ^ (pDir[i] - sign_mask - 1);
        sStartSite.m_byData4[diridx] += sign_mask;
#if _CLG_ASSUME_SQUARE_LATTICE
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
        //if (debug)
        //{
        //    printf("%d link (%d)(%d, %d, %d, %d)_ %d(%d)\n", i, newLink.m_uiSiteIndex, sStartSite.x, sStartSite.y, sStartSite.z, sStartSite.w, newLink.m_byDir, newLink.IsDirichlet());
        //}
        if (newLink.IsDirichlet())
        {
            sStartSite.m_byData4[diridx] += (sign_mask + 1);
            continue;
        }
        if (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask))
        {
            _muldag(sRet, pDeviceDatas[i == replaceIndex][_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
        }
        else
        {
            _mul(sRet, pDeviceDatas[i == replaceIndex][_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
        }
        sStartSite.m_byData4[diridx] += (sign_mask + 1);
    }

    return sRet;
}

/**
* The first move of the path is skipped
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceLinkTSkipOne(
    const deviceGauge* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const SCHAR* __restrict__ pDir)
{
    if (byLength < 2)
    {
        return _makeId<deviceGauge>();
    }

    //first move
    SCHAR sign_mask = (pDir[0] >> (sizeof(SCHAR) * 8 - 1));
    SCHAR diridx = sign_mask ^ (pDir[0] - sign_mask - 1);
    sStartSite.m_byData4[diridx] += (2 * sign_mask + 1);

    //second move
    sign_mask = (pDir[1] >> (sizeof(SCHAR) * 8 - 1));
    diridx = sign_mask ^ (pDir[1] - sign_mask - 1);
    sStartSite.m_byData4[diridx] += sign_mask;
#if _CLG_ASSUME_SQUARE_LATTICE
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
    deviceGauge sRet = newLink.IsDirichlet() ? (_makeId<deviceGauge>()) : pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
    if (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask))
    {
        _dagger(sRet);
    }
    sStartSite.m_byData4[diridx] += (sign_mask + 1);

    for (BYTE i = 2U; i < byLength; ++i)
    {
        sign_mask = (pDir[i] >> (sizeof(SCHAR) * 8 - 1));
        diridx = sign_mask ^ (pDir[i] - sign_mask - 1);
        sStartSite.m_byData4[diridx] += sign_mask;
#if _CLG_ASSUME_SQUARE_LATTICE
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
        if (newLink.IsDirichlet())
        {
            sStartSite.m_byData4[diridx] += (sign_mask + 1);
            continue;
        }
        if (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask))
        {
            _muldag(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
        }
        else
        {
            _mul(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
        }
        sStartSite.m_byData4[diridx] += (sign_mask + 1);
    }

    return sRet;
}

template<>
__device__ __inline__ Real _deviceLinkT<Real>(
    const Real* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const SCHAR* __restrict__ pDir)
{
    //length can be 0
    if (0 == byLength)
    {
        return F(0.0);
    }

    SCHAR sign_mask = static_cast<SCHAR>(pDir[0] >> (sizeof(SCHAR) * 8 - 1));
    SCHAR diridx = sign_mask ^ (pDir[0] - sign_mask - 1);
    sStartSite.m_byData4[diridx] += sign_mask;
#if _CLG_ASSUME_SQUARE_LATTICE
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
    Real fRet = newLink.IsDirichlet() ? F(0.0) : pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
    fRet = fRet * (1 - (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask) << 1));
    sStartSite.m_byData4[diridx] += sign_mask + 1;

    for (BYTE i = 1U; i < byLength; ++i)
    {
        sign_mask = (pDir[i] >> (sizeof(SCHAR) * 8 - 1));
        diridx = sign_mask ^ (pDir[i] - sign_mask - 1);
        sStartSite.m_byData4[diridx] += static_cast<SCHAR>(sign_mask);
#if _CLG_ASSUME_SQUARE_LATTICE
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
        if (newLink.IsDirichlet())
        {
            sStartSite.m_byData4[diridx] += static_cast<SCHAR>(sign_mask + 1);
            continue;
        }
        fRet += pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)] * (1 - (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask) << 1));
        sStartSite.m_byData4[diridx] += static_cast<SCHAR>(sign_mask + 1);
    }

    return fRet;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceLinkEMT(
    const deviceGauge* __restrict__ pDeviceData,
    const Real* __restrict__ pU1Real,
    Real fCharge,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const SCHAR* __restrict__ pDir)
{
    if (0 == byLength)
    {
        return _makeId<deviceGauge>();
    }

    SCHAR sign_mask = (pDir[0] >> (sizeof(SCHAR) * 8 - 1));
    SCHAR diridx = sign_mask ^ (pDir[0] - sign_mask - 1);
    sStartSite.m_byData4[diridx] += sign_mask;
#if _CLG_ASSUME_SQUARE_LATTICE
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
    UINT linkIndex = _deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir);
    deviceGauge sRet = newLink.IsDirichlet() ? (_makeId<deviceGauge>()) : pDeviceData[linkIndex];
    Real fPhase = newLink.IsDirichlet() ? F(0.0) : pU1Real[linkIndex];
    if (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask))
    {
        _dagger(sRet);
        fPhase = -fPhase;
    }
    sStartSite.m_byData4[diridx] += (sign_mask + 1);

    for (BYTE i = 1U; i < byLength; ++i)
    {
        sign_mask = (pDir[i] >> (sizeof(SCHAR) * 8 - 1));
        diridx = sign_mask ^ (pDir[i] - sign_mask - 1);
        sStartSite.m_byData4[diridx] += sign_mask;
#if _CLG_ASSUME_SQUARE_LATTICE
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
        if (newLink.IsDirichlet())
        {
            sStartSite.m_byData4[diridx] += (sign_mask + 1);
            continue;
        }

        linkIndex = _deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir);
        if (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask))
        {
            _muldag(sRet, pDeviceData[linkIndex]);
            fPhase -= pU1Real[linkIndex];
        }
        else
        {
            _mul(sRet, pDeviceData[linkIndex]);
            fPhase += pU1Real[linkIndex];
        }
        sStartSite.m_byData4[diridx] += (sign_mask + 1);
    }

    fPhase = fPhase * fCharge;
    _mul(sRet, _make_cuComplex(_cos(fPhase), _sin(fPhase)));
    return sRet;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceLinkEMTSkipOne(
    const deviceGauge* __restrict__ pDeviceData,
    const Real* __restrict__ pU1Real,
    Real fCharge,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const SCHAR* __restrict__ pDir)
{
    if (byLength < 2)
    {
        return _makeId<deviceGauge>();
    }

    //first move
    SCHAR sign_mask = (pDir[0] >> (sizeof(SCHAR) * 8 - 1));
    SCHAR diridx = sign_mask ^ (pDir[0] - sign_mask - 1);
    sStartSite.m_byData4[diridx] += (2 * sign_mask + 1);

    //second move
    sign_mask = (pDir[1] >> (sizeof(SCHAR) * 8 - 1));
    diridx = sign_mask ^ (pDir[1] - sign_mask - 1);
    sStartSite.m_byData4[diridx] += sign_mask;
#if _CLG_ASSUME_SQUARE_LATTICE
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
    UINT linkIndex = _deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir);
    deviceGauge sRet = newLink.IsDirichlet() ? (_makeId<deviceGauge>()) : pDeviceData[linkIndex];
    Real fPhase = newLink.IsDirichlet() ? F(0.0) : pU1Real[linkIndex];
    if (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask))
    {
        _dagger(sRet);
        fPhase = -fPhase;
    }
    sStartSite.m_byData4[diridx] += (sign_mask + 1);

    for (BYTE i = 2U; i < byLength; ++i)
    {
        sign_mask = (pDir[i] >> (sizeof(SCHAR) * 8 - 1));
        diridx = sign_mask ^ (pDir[i] - sign_mask - 1);
        sStartSite.m_byData4[diridx] += sign_mask;
#if _CLG_ASSUME_SQUARE_LATTICE
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
        if (newLink.IsDirichlet())
        {
            sStartSite.m_byData4[diridx] += (sign_mask + 1);
            continue;
        }

        linkIndex = _deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir);
        if (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask))
        {
            _muldag(sRet, pDeviceData[linkIndex]);
            fPhase -= pU1Real[linkIndex];
        }
        else
        {
            _mul(sRet, pDeviceData[linkIndex]);
            fPhase += pU1Real[linkIndex];
        }
        sStartSite.m_byData4[diridx] += (sign_mask + 1);
    }

    fPhase = fPhase * fCharge;
    _mul(sRet, _make_cuComplex(_cos(fPhase), _sin(fPhase)));
    return sRet;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _devicePlaneDiagonalEMT(
    const deviceGauge* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataReal,
    Real fCharge,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    SCHAR dim1, SCHAR dim2)
{
    SCHAR dir1[2];

    dir1[0] = dim1; dir1[1] = dim2;
    deviceGauge sRet(_deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 2, byFieldId, dir1));

    dir1[0] = dim2; dir1[1] = dim1;
    _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 2, byFieldId, dir1));

    _mul(sRet, F(0.5));
    return sRet;
}

/**
 * dim1, 2, 3 =
 * 1: x, -1: -x
 * 2: y, -2: -y
 * 3: z, -3: -z
 * 4: t, -4: -t
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceCubicDiagonalEMT(
    const deviceGauge* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataReal,
    Real fCharge,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    SCHAR dim1, SCHAR dim2, SCHAR dim3)
{
    SCHAR dir1[3];

    dir1[0] = dim1; dir1[1] = dim2; dir1[2] = dim3;
    deviceGauge sRet(_deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim1; dir1[1] = dim3; dir1[2] = dim2;
    _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim2; dir1[1] = dim1; dir1[2] = dim3;
    _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim2; dir1[1] = dim3; dir1[2] = dim1;
    _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim3; dir1[1] = dim1; dir1[2] = dim2;
    _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 3, byFieldId, dir1));

    dir1[0] = dim3; dir1[1] = dim2; dir1[2] = dim1;
    _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 3, byFieldId, dir1));

    _mul(sRet, OneOver6);
    return sRet;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceHyperCubicDiagonalEMT(
    const deviceGauge* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataReal,
    Real fCharge,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    SCHAR dim1, SCHAR dim2, SCHAR dim3, SCHAR dim4)
{
    deviceGauge sRet = _makeZero<deviceGauge>();
    const SCHAR dim1234[4] = { dim1, dim2, dim3, dim4 };
    SCHAR dir1[4];
    SCHAR dim234[3];
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
        _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[0]; dir1[2] = dim234[2]; dir1[3] = dim234[1];
        _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[1]; dir1[2] = dim234[0]; dir1[3] = dim234[2];
        _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[1]; dir1[2] = dim234[2]; dir1[3] = dim234[0];
        _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[2]; dir1[2] = dim234[0]; dir1[3] = dim234[1];
        _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 4, byFieldId, dir1));

        dir1[1] = dim234[2]; dir1[2] = dim234[1]; dir1[3] = dim234[0];
        _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataReal, fCharge, sStartSite, 4, byFieldId, dir1));
    }

    _mul(sRet, OneOver24);
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
    const SCHAR* __restrict__ pDir)
{
    if (0 == byLength)
    {
        return _makeId<deviceGauge>();
    }

    //UBOOL bLog = FALSE;
    //if (0 == sStartSite.x && 0 == sStartSite.y && 0 == sStartSite.z && 0 == sStartSite.w)
    //{
    //    bLog = TRUE;
    //}
    SCHAR sign_mask = (pDir[0] >> (sizeof(SCHAR) * 8 - 1));
    SCHAR diridx = sign_mask ^ (pDir[0] - sign_mask - 1);
    sStartSite.m_byData4[diridx] += sign_mask;
#if _CLG_ASSUME_SQUARE_LATTICE
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
    SIndex newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
    deviceGauge sRet = newLink.IsDirichlet() ? (_makeId<deviceGauge>()) : pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
    if (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask))
    {
        _dagger(sRet);
    }
    //if (bLog)
    //{
    //    printf("sign mask:%d ", sign_mask);
    //    newLink.DebugPrint();
    //}
    sStartSite.m_byData4[diridx] += (sign_mask + 1);

    for (BYTE i = 1U; i < byLength; ++i)
    {
        sign_mask = (pDir[i] >> (sizeof(SCHAR) * 8 - 1));
        diridx = sign_mask ^ (pDir[i] - sign_mask - 1);
        sStartSite.m_byData4[diridx] += sign_mask;
#if _CLG_ASSUME_SQUARE_LATTICE
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) | diridx];
#else
        newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + diridx];
#endif
        if (newLink.IsDirichlet())
        {
            sStartSite.m_byData4[diridx] += (sign_mask + 1);
            //only move once, since we have valid margin, we don't worry
            sStartSite = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sStartSite)].m_uiSiteIndex);
            continue;
        }
        //if (bLog)
        //{
        //    printf("sign mask:%d ", sign_mask);
        //    newLink.DebugPrint();
        //}
        if (_UBOOLXOR(newLink.NeedToDagger(), -sign_mask))
        {
            _muldag(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
        }
        else
        {
            _mul(sRet, pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
        }
        sStartSite.m_byData4[diridx] += (sign_mask + 1);
        //only move once, since we have valid margin, we don't worry
        sStartSite = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sStartSite)].m_uiSiteIndex);
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
template<typename deviceGauge>
static __device__ __inline__ Real _deviceOneRectangleRetrT(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sSite, SCHAR iMu, SCHAR iNu)
{
    SCHAR path[6] = { iMu, iNu, static_cast<SCHAR>(-iMu), static_cast<SCHAR>(-iMu), static_cast<SCHAR>(-iNu), iMu };
    return _retr(_deviceLinkT(pDeviceData, sSite, 6, byFieldId, path));
}

template<typename deviceGauge>
static __device__ __inline__ Real _deviceCloverRectangleRetrT(
    const BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sSite, BYTE byMu, BYTE byNu)
{
    const SCHAR ifwdMu = __fwd(byMu);
    const SCHAR ifwdNu = __fwd(byNu);
    Real fRes = _deviceOneRectangleRetrT(byFieldId, pDeviceData, sSite, ifwdMu, ifwdNu);
    fRes += _deviceOneRectangleRetrT(byFieldId, pDeviceData, sSite, ifwdMu, -ifwdNu);
    fRes += _deviceOneRectangleRetrT(byFieldId, pDeviceData, sSite, ifwdNu, ifwdMu);
    fRes += _deviceOneRectangleRetrT(byFieldId, pDeviceData, sSite, ifwdNu, -ifwdMu);
    return fRes;
}

#pragma endregion

#pragma endregion

#pragma region device functions Measure Topological charge XY

template<typename deviceGauge>
static __device__ __inline__ Real _deviceTrImCloverT(const deviceGauge* __restrict__ pGaugeField, BYTE byFieldId, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE rho, BYTE sigma)
{
    return _trim(
        _deviceCloverT(pGaugeField, sSite4, uiBigIdx, mu, nu, byFieldId),
        _deviceCloverT(pGaugeField, sSite4, uiBigIdx, rho, sigma, byFieldId));
}

template<typename deviceGauge>
static __device__ __inline__ Real _deviceTopologicalChargeT(const deviceGauge* __restrict__ pGaugeField, BYTE byFieldId, const SSmallInt4& sSite4, UINT uiBigIdx)
{
    Real ret = _deviceTrImCloverT(pGaugeField, byFieldId, sSite4, uiBigIdx, 0, 1, 2, 3);
    ret -= _deviceTrImCloverT(pGaugeField, byFieldId, sSite4, uiBigIdx, 0, 2, 1, 3);
    ret += _deviceTrImCloverT(pGaugeField, byFieldId, sSite4, uiBigIdx, 0, 3, 1, 2);
    return OneOver32PI2 * F(2.0) * ret;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _DEVICEINLINEGAUGE_H_

//=============================================================================
// END OF FILE
//=============================================================================
