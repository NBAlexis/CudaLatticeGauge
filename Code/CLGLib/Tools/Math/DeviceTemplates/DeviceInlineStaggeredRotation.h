//=============================================================================
// FILENAME : DeviceInlineStaggeredRotation.h
// 
// DESCRIPTION:
// This should be implemented using inherint machinism, but due to historical reasons, it is now templates
//
//
// REVISION:
//  [07/14/2024 nbale]
//=============================================================================

#ifndef _DEVICEINLINE_STAGGEREDROTATION_H_
#define _DEVICEINLINE_STAGGEREDROTATION_H_

__BEGIN_NAMESPACE

#if 0
/**
 * V(n1,n2)
 */
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceVXXTau(
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    UINT bXorY, UBOOL bPlusMu, UBOOL bPlusTau)
{
    const INT iMu = bXorY ? (bPlusMu ? 1 : -1) : (bPlusMu ? 2 : -2);
    const INT iTau = bPlusTau ? 4 : -4;
    INT dir1[3];

    dir1[0] = iMu;
    dir1[1] = iMu;
    dir1[2] = iTau;
    deviceGauge sRet = _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1);

    dir1[0] = iMu;
    dir1[1] = iTau;
    dir1[2] = iMu;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iTau;
    dir1[1] = iMu;
    dir1[2] = iMu;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    _mul(sRet, OneOver3);
    return sRet;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceVXYT(
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    UBOOL bPlusX, UBOOL bPlusY, UBOOL bPlusTau)
{
    const INT iX = bPlusX ? 1 : -1;
    const INT iY = bPlusY ? 2 : -2;
    const INT iT = bPlusTau ? 4 : -4;
    INT dir1[3];

    dir1[0] = iX; dir1[1] = iY; dir1[2] = iT;
    deviceGauge sRet(_deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iX; dir1[1] = iT; dir1[2] = iY;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iY; dir1[1] = iX; dir1[2] = iT;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iY; dir1[1] = iT; dir1[2] = iX;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iT; dir1[1] = iX; dir1[2] = iY;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iT; dir1[1] = iY; dir1[2] = iX;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    _mul(sRet, OneOver6);
    return sRet;
}

#endif

//11 sec
/**
* byMu = bXorY ? 0 : 1
* so that, when bXorY = 1, it is partial_X
*          when bXorY = 0, it is partial_Y
*/
template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceVXXTauOptimizedT(
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    UINT bXorY, UBOOL bPlusMu, UBOOL bPlusTau)
{
    const INT byMu = bXorY ? 0 : 1;
    const SCHAR iMu = bPlusMu ? (byMu + 1) : (-byMu - 1);
    const SCHAR iTau = bPlusTau ? 4 : -4;
    SCHAR dir1[3];

    SSmallInt4 x_p_mu_p_tau = sStartSite;
    x_p_mu_p_tau.w = x_p_mu_p_tau.w + (bPlusTau ? 1 : -1);
    //x_p_mu_p_tau.m_byData4[byMu] = x_p_mu_p_tau.m_byData4[byMu] + (bPlusMu ? 1 : -1);
    x_p_mu_p_tau.m_byData4[byMu] = x_p_mu_p_tau.m_byData4[byMu] + (bPlusMu ? 1 : -2);

    dir1[0] = iMu; dir1[1] = iTau;
    deviceGauge sRet = _deviceLinkT(pDeviceData, sStartSite, 2, byFieldId, dir1);
    dir1[0] = iTau; dir1[1] = iMu;
    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iMu;
    //_mul(sRet, _deviceLinkT(pDeviceData, x_p_mu_p_tau, 1, byFieldId, dir1));

    const SIndex& x_p_taumu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(x_p_mu_p_tau) + byMu];
    if ((x_p_taumu__mu.NeedToDagger() && bPlusMu)
        || (!x_p_taumu__mu.NeedToDagger() && !bPlusMu))
    {
        _muldag(sRet, pDeviceData[_deviceGetLinkIndex(x_p_taumu__mu.m_uiSiteIndex, x_p_taumu__mu.m_byDir)]);
    }
    else
    {
        _mul(sRet, pDeviceData[_deviceGetLinkIndex(x_p_taumu__mu.m_uiSiteIndex, x_p_taumu__mu.m_byDir)]);
    }

    dir1[0] = iMu;
    dir1[1] = iMu;
    dir1[2] = iTau;

    _add(sRet, _deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

    _mul(sRet, OneOver3);
    return sRet;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceVXYTOptimizedT(
    const deviceGauge* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    UBOOL bPlusX, UBOOL bPlusY, UBOOL bPlusTau)
{
    const SCHAR iX = bPlusX ? 1 : -1;
    const SCHAR iY = bPlusY ? 2 : -2;
    const SCHAR iT = bPlusTau ? 4 : -4;

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

    SCHAR dir1[2];

    dir1[0] = iX; dir1[1] = iY;
    deviceGauge sRet1(_deviceLinkT(pDeviceData, sStartSite, 2, byFieldId, dir1));

    dir1[0] = iY; dir1[1] = iX;
    _add(sRet1, _deviceLinkT(pDeviceData, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iT;
    //sRet1.Mul(_deviceLinkT(pDeviceData, n_xy, 1, byFieldId, dir1));
    const SIndex& n_xy__t = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_xy) + 3];
    if ((n_xy__t.NeedToDagger() && bPlusTau)
        || (!n_xy__t.NeedToDagger() && !bPlusTau))
    {
        _muldag(sRet1, pDeviceData[_deviceGetLinkIndex(n_xy__t.m_uiSiteIndex, n_xy__t.m_byDir)]);
    }
    else
    {
        _mul(sRet1, pDeviceData[_deviceGetLinkIndex(n_xy__t.m_uiSiteIndex, n_xy__t.m_byDir)]);
    }

    dir1[0] = iX; dir1[1] = iT;
    deviceGauge sRet2 = _deviceLinkT(pDeviceData, sStartSite, 2, byFieldId, dir1);

    dir1[0] = iT; dir1[1] = iX;
    _add(sRet2, _deviceLinkT(pDeviceData, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iY;
    //sRet2.Mul(_deviceLinkT(pDeviceData, n_xt, 1, byFieldId, dir1));
    const SIndex& n_xt__y = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_xt) + 1];
    if ((n_xt__y.NeedToDagger() && bPlusY)
        || (!n_xt__y.NeedToDagger() && !bPlusY))
    {
        _muldag(sRet2, pDeviceData[_deviceGetLinkIndex(n_xt__y.m_uiSiteIndex, n_xt__y.m_byDir)]);
    }
    else
    {
        _mul(sRet2, pDeviceData[_deviceGetLinkIndex(n_xt__y.m_uiSiteIndex, n_xt__y.m_byDir)]);
    }

    _add(sRet1, sRet2);

    dir1[0] = iY; dir1[1] = iT;
    sRet2 = _deviceLinkT(pDeviceData, sStartSite, 2, byFieldId, dir1);

    dir1[0] = iT; dir1[1] = iY;
    _add(sRet2, _deviceLinkT(pDeviceData, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iX;
    //sRet2.Mul(_deviceLinkT(pDeviceData, n_yt, 1, byFieldId, dir1));
    const SIndex& n_yt__x = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_yt)];
    if ((n_yt__x.NeedToDagger() && bPlusX)
        || (!n_yt__x.NeedToDagger() && !bPlusX))
    {
        _muldag(sRet2, pDeviceData[_deviceGetLinkIndex(n_yt__x.m_uiSiteIndex, n_yt__x.m_byDir)]);
    }
    else
    {
        _mul(sRet2, pDeviceData[_deviceGetLinkIndex(n_yt__x.m_uiSiteIndex, n_yt__x.m_byDir)]);
    }

    _add(sRet1, sRet2);

    _mul(sRet1, OneOver6);
    return sRet1;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceVXXTauOptimizedEMT(
    const deviceGauge* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataU1,
    const SSmallInt4& sStartSite,
    Real fCharge, BYTE byFieldId,
    UINT bXorY, UBOOL bPlusMu, UBOOL bPlusTau)
{
    const SCHAR byMu = bXorY ? 0 : 1;
    const SCHAR iMu = bPlusMu ? (byMu + 1) : (-byMu - 1);
    const SCHAR iTau = bPlusTau ? 4 : -4;
    SCHAR dir1[3];

    SSmallInt4 x_p_mu_p_tau = sStartSite;
    x_p_mu_p_tau.w = x_p_mu_p_tau.w + (bPlusTau ? 1 : -1);
    //x_p_mu_p_tau.m_byData4[byMu] = x_p_mu_p_tau.m_byData4[byMu] + (bPlusMu ? 1 : -1);
    x_p_mu_p_tau.m_byData4[byMu] = x_p_mu_p_tau.m_byData4[byMu] + (bPlusMu ? 1 : -2);

    //Sum of the left-squre
    dir1[0] = iMu; dir1[1] = iTau;
    deviceGauge sRet = _deviceLinkEMT(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1);
    dir1[0] = iTau; dir1[1] = iMu;
    _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iMu;
    //sRet.Mul(_deviceLink(pDeviceData, x_p_mu_p_tau, 1, byFieldId, dir1));

    const SIndex& x_p_taumu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(x_p_mu_p_tau) + byMu];
    const UBOOL bLastLinkDagger = (x_p_taumu__mu.NeedToDagger() && bPlusMu)
        || (!x_p_taumu__mu.NeedToDagger() && !bPlusMu);
    const UINT x_p_taumu__mu__link = _deviceGetLinkIndex(x_p_taumu__mu.m_uiSiteIndex, x_p_taumu__mu.m_byDir);
    const Real x_p_taumu__mu__phase = pDeviceDataU1[x_p_taumu__mu__link] * fCharge;
    if (bLastLinkDagger)
    {
        _muldag(sRet, pDeviceData[x_p_taumu__mu__link]);
        _mul(sRet, _make_cuComplex(_cos(x_p_taumu__mu__phase), -_sin(x_p_taumu__mu__phase)));
    }
    else
    {
        _mul(sRet, pDeviceData[x_p_taumu__mu__link]);
        _mul(sRet, _make_cuComplex(_cos(x_p_taumu__mu__phase), _sin(x_p_taumu__mu__phase)));
    }

    //The final line
    dir1[0] = iMu;
    dir1[1] = iMu;
    dir1[2] = iTau;

    _add(sRet, _deviceLinkEMT(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 3, byFieldId, dir1));

    _mul(sRet, OneOver3);
    return sRet;
}

template<typename deviceGauge>
static __device__ __inline__ deviceGauge _deviceVXYTOptimizedEMT(
    const deviceGauge* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataU1,
    const SSmallInt4& sStartSite,
    Real fCharge, BYTE byFieldId,
    UBOOL bPlusX, UBOOL bPlusY, UBOOL bPlusTau)
{
    const SCHAR iX = bPlusX ? 1 : -1;
    const SCHAR iY = bPlusY ? 2 : -2;
    const SCHAR iT = bPlusTau ? 4 : -4;

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

    SCHAR dir1[2];

    dir1[0] = iX; dir1[1] = iY;
    deviceGauge sRet1(_deviceLinkEMT(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1));

    dir1[0] = iY; dir1[1] = iX;
    _add(sRet1, _deviceLinkEMT(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iT;
    //sRet1.Mul(_deviceLink(pDeviceData, n_xy, 1, byFieldId, dir1));
    const SIndex& n_xy__t = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_xy) + 3];
    const UINT n_xy__t_link = _deviceGetLinkIndex(n_xy__t.m_uiSiteIndex, n_xy__t.m_byDir);
    const Real n_xy__t_phase = pDeviceDataU1[n_xy__t_link] * fCharge;
    if ((n_xy__t.NeedToDagger() && bPlusTau) || (!n_xy__t.NeedToDagger() && !bPlusTau))
    {
        //at present, we only consider magnetic, so there is no t-dir
        _muldag(sRet1, pDeviceData[n_xy__t_link]);
        _mul(sRet1, _make_cuComplex(_cos(n_xy__t_phase), -_sin(n_xy__t_phase)));
    }
    else
    {
        //at present, we only consider magnetic, so there is no t-dir
        _mul(sRet1, pDeviceData[n_xy__t_link]);
        _mul(sRet1, _make_cuComplex(_cos(n_xy__t_phase), _sin(n_xy__t_phase)));
    }

    dir1[0] = iX; dir1[1] = iT;
    deviceGauge sRet2 = _deviceLinkEMT(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1);

    dir1[0] = iT; dir1[1] = iX;
    _add(sRet2, _deviceLinkEMT(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iY;
    //sRet2.Mul(_deviceLink(pDeviceData, n_xt, 1, byFieldId, dir1));

    //================== add magnetic phase y ========================
    const SIndex& n_xt__y = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_xt) + 1];
    const UINT n_xt__y_link = _deviceGetLinkIndex(n_xt__y.m_uiSiteIndex, n_xt__y.m_byDir);
    const Real n_xt__y_phase = pDeviceDataU1[n_xt__y_link] * fCharge;
    if ((n_xt__y.NeedToDagger() && bPlusY) || (!n_xt__y.NeedToDagger() && !bPlusY))
    {
        _muldag(sRet2, pDeviceData[n_xt__y_link]);
        _mul(sRet2, _make_cuComplex(_cos(n_xt__y_phase), -_sin(n_xt__y_phase)));
    }
    else
    {
        _mul(sRet2, pDeviceData[n_xt__y_link]);
        _mul(sRet2, _make_cuComplex(_cos(n_xt__y_phase), _sin(n_xt__y_phase)));
    }

    _add(sRet1, sRet2);

    dir1[0] = iY; dir1[1] = iT;
    sRet2 = _deviceLinkEMT(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1);

    dir1[0] = iT; dir1[1] = iY;
    _add(sRet2, _deviceLinkEMT(pDeviceData, pDeviceDataU1, fCharge, sStartSite, 2, byFieldId, dir1));

    //dir1[0] = iX;
    //sRet2.Mul(_deviceLink(pDeviceData, n_yt, 1, byFieldId, dir1));

    //================== add magnetic phase x ========================
    const SIndex& n_yt__x = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(n_yt)];
    const UINT n_yt__x_link = _deviceGetLinkIndex(n_yt__x.m_uiSiteIndex, n_yt__x.m_byDir);
    const Real n_yt__x_phase = pDeviceDataU1[n_yt__x_link] * fCharge;
    if ((n_yt__x.NeedToDagger() && bPlusX) || (!n_yt__x.NeedToDagger() && !bPlusX))
    {
        _muldag(sRet2, pDeviceData[n_yt__x_link]);
        _mul(sRet2, _make_cuComplex(_cos(n_yt__x_phase), -_sin(n_yt__x_phase)));
    }
    else
    {
        _mul(sRet2, pDeviceData[n_yt__x_link]);
        _mul(sRet2, _make_cuComplex(_cos(n_yt__x_phase), _sin(n_yt__x_phase)));
    }

    _add(sRet1, sRet2);

    _mul(sRet1, OneOver6);
    return sRet1;
}


__END_NAMESPACE

#endif //#ifndef _DEVICEINLINE_STAGGEREDROTATION_H_

//=============================================================================
// END OF FILE
//=============================================================================
