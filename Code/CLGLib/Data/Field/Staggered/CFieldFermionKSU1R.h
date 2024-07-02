//=============================================================================
// FILENAME : CFieldFermionKSU1R.h
// 
// DESCRIPTION:
// This is the class for Kogut-Susskind staggered fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [10/03/2021 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONKSU1R_H_
#define _CFIELDFERMIONKSU1R_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSU1R)

class CLGAPI CFieldFermionKSU1R : public CFieldFermionKSU1
{
    __CLGDECLARE_FIELD(CFieldFermionKSU1R)

protected:

    void DerivateD0(void* pForce, const void* pGaugeBuffer, BYTE byGaugeFieldId) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

public:

    void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString& tab) const override;

    static void Seperate(INT* full, INT iSep, INT* l, INT* r, BYTE& LL, BYTE& RL)
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
};

#pragma region device functions



/**
 * V(n1,n2)
 */
//static __device__ __inline__ CLGComplex _deviceVXXTauU1(
//    const CLGComplex* __restrict__ pDeviceData,
//    const SSmallInt4& sStartSite, BYTE byFieldId,
//    UINT bXorY, UBOOL bPlusMu, UBOOL bPlusTau)
//{
//    const INT iMu = bXorY ? (bPlusMu ? 1 : -1) : (bPlusMu ? 2 : -2);
//    const INT iTau = bPlusTau ? 4 : -4;
//    INT dir1[3];
//
//    dir1[0] = iMu;
//    dir1[1] = iMu;
//    dir1[2] = iTau;
//    CLGComplex sRet = _deviceLinkU1(pDeviceData, sStartSite, 3, byFieldId, dir1);
//    
//    dir1[0] = iMu;
//    dir1[1] = iTau;
//    dir1[2] = iMu;
//    sRet = _cuCaddf(sRet, _deviceLinkU1(pDeviceData, sStartSite, 3, byFieldId, dir1));
//
//    dir1[0] = iTau;
//    dir1[1] = iMu;
//    dir1[2] = iMu;
//    sRet = _cuCaddf(sRet, _deviceLinkU1(pDeviceData, sStartSite, 3, byFieldId, dir1));
//
//    sRet = cuCmulf_cr(sRet, OneOver3);
//    return sRet;
//}
//
//static __device__ __inline__ CLGComplex _deviceVXYTU1(
//    const CLGComplex* __restrict__ pDeviceData,
//    const SSmallInt4& sStartSite, BYTE byFieldId,
//    UBOOL bPlusX, UBOOL bPlusY, UBOOL bPlusTau)
//{
//    const INT iX = bPlusX ? 1 : -1;
//    const INT iY = bPlusY ? 2 : -2;
//    const INT iT = bPlusTau ? 4 : -4;
//    INT dir1[3];
//
//    dir1[0] = iX; dir1[1] = iY; dir1[2] = iT;
//    deviceSU3 sRet(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));
//
//    dir1[0] = iX; dir1[1] = iT; dir1[2] = iY;
//    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));
//
//    dir1[0] = iY; dir1[1] = iX; dir1[2] = iT;
//    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));
//
//    dir1[0] = iY; dir1[1] = iT; dir1[2] = iX;
//    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));
//
//    dir1[0] = iT; dir1[1] = iX; dir1[2] = iY;
//    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));
//
//    dir1[0] = iT; dir1[1] = iY; dir1[2] = iX;
//    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));
//
//    sRet.MulReal(OneOver6);
//    return sRet;
//}


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

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSU1R_H_

//=============================================================================
// END OF FILE
//=============================================================================