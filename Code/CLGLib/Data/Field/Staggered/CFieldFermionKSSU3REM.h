//=============================================================================
// FILENAME : CFieldFermionKSSU3REM.h
// 
// DESCRIPTION:
// This is the class for Kogut-Susskind staggered fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// Finally, I decide to implement Magnetic field only...
//
// REVISION:
//  [12/21/2021 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONKSSU3REM_H_
#define _CFIELDFERMIONKSSU3REM_H_

#define __AlsoElectric 0

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3REM)

class CLGAPI CFieldFermionKSSU3REM : public CFieldFermionKSSU3
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3REM)

public:

    CFieldFermionKSSU3REM();
    ~CFieldFermionKSSU3REM();

protected:

    void DerivateD0(void* pForce, const void* pGaugeBuffer) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;
    void ApplyGammaKSS(const CFieldGauge* pGauge, EGammaMatrix eGamma) override;

public:

    void InitialOtherParameters(CParameters& params) override;
    CCString GetInfos(const CCString& tab) const override;

    //BYTE m_byGaugeType;
    //UBOOL m_bTwistedBoundary;
    BYTE m_byEMFieldID;
    Real m_fQ;
    INT* m_pDevicePathBuffer;
};

#pragma region device functions



/**
 * V(n1,n2)
 */
static __device__ __inline__ deviceSU3 _deviceVXXTauEM(
    const deviceSU3* __restrict__ pDeviceData,
    const Real* __restrict__ pDeviceDataU1,
    const SSmallInt4& sStartSite,  Real fCharge, BYTE byFieldId,
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

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3REM_H_

//=============================================================================
// END OF FILE
//=============================================================================