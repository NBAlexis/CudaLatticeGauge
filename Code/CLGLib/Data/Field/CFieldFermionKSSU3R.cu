//=============================================================================
// FILENAME : CFieldFermionKSSU3R.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [09/23/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionKSSU3R)

#pragma region DOperator

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
static __device__ __inline__ deviceSU3 _deviceVXXTauOptimized(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    UINT bXorY, UBOOL bPlusMu, UBOOL bPlusTau)
{
    const INT byMu = bXorY ? 0 : 1;
    const INT iMu = bPlusMu ? (byMu + 1) : (- byMu - 1);
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

#pragma region kernel

/**
* When link n and n+mu, the coordinate is stick with n
* When link n and n-mu, the coordinate is stick with n-mu
* Irrelavent with tau
* Optimization: bXorY removed, block.x *= 2 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTerm(
    const deviceSU3Vector * __restrict__ pDeviceData,
    const deviceSU3 * __restrict__ pGauge,
    const BYTE * __restrict__ pEtaTable,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();

    #pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusMu  = idx & 2;
        const UBOOL bPlusTau = idx & 4;
        const UINT bXorY = idx & 1;
        const UINT bYorX = 1 - bXorY;
        SSmallInt4 sTargetSite = sSite4;
        SSmallInt4 sMidSite = sSite4;
        sTargetSite.m_byData4[bYorX] = sTargetSite.m_byData4[bYorX] + (bPlusMu ? 2 : -2);
        sMidSite.m_byData4[bYorX] = sMidSite.m_byData4[bYorX] + (bPlusMu ? 1 : -1);
        sTargetSite.w = sTargetSite.w + (bPlusTau ? 1 : -1);
        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sTargetSite)];
        const SIndex& sMiddleBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sMidSite)];
        sMidSite = __deviceSiteIndexToInt4(sMiddleBigIndex.m_uiSiteIndex);
        
        INT eta_tau = (pEtaTable[sMiddleBigIndex.m_uiSiteIndex] >> 3) & 1;
        eta_tau = eta_tau + bXorY;       
        if (sTargetBigIndex.NeedToOpposite())
        {
            eta_tau = eta_tau + 1;
        }

        deviceSU3Vector right = _deviceVXXTauOptimized(pGauge, sSite4, byGaugeFieldId, bXorY, bPlusMu, bPlusTau).MulVector(
            pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        right.MulReal(sCenter.m_byData4[bXorY] - sMidSite.m_byData4[bXorY] - F(0.5));

        if (bPlusMu)
        {
            eta_tau = eta_tau + 1;
        }

        if (eta_tau & 1)
        {
            result.Add(right);
        }
        else
        {
            result.Sub(right);
        }
    }

    if (bDDagger)
    {
        fOmega = -fOmega;
    }
    result.MulReal(F(0.25) * fOmega);

    switch (eCoeff)
    {
    case EOCT_Real:
        result.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        result.MulComp(cCoeff);
        break;
    }

    pResultData[uiSiteIndex].Add(result);
}


__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTau_Term(
    const deviceSU3Vector* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();

    #pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusX = (0 != (idx & 1));
        const UBOOL bPlusY = (0 != (idx & 2));
        const UBOOL bPlusT = (0 != (idx & 4));

        SSmallInt4 sOffset = sSite4;
        sOffset.x = sOffset.x + (bPlusX ? 1 : -1);
        sOffset.y = sOffset.y + (bPlusY ? 1 : -1);
        sOffset.w = sOffset.w + (bPlusT ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];
        
        const deviceSU3Vector right = _deviceVXYTOptimized(pGauge, sSite4, byGaugeFieldId, bPlusX, bPlusY, bPlusT)
        .MulVector(pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta124 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta124 = bPlusT ? (sSite4.y + sSite4.z): (site_target.y + site_target.z + 1);
        if (sTargetBigIndex.NeedToOpposite())
        {
            eta124 = eta124 + 1;
        }

        if (eta124 & 1)
        {
            result.Add(right);
        }
        else
        {
            result.Sub(right);
        }
    }

    if (bDDagger)
    {
        result.MulReal(-F(0.125) * fOmega);
    }
    else
    {
        result.MulReal(F(0.125) * fOmega);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        result.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        result.MulComp(cCoeff);
        break;
    }

    pResultData[uiSiteIndex].Add(result);
}

#pragma endregion

#pragma region Derivate

/**
 * Have n, n->n1, n->n2,
 * 1. we need to obtain V_(n, n1) , V_(n, n2)
 * 2. we need phi(n1), phi(n2), phid(n1), phid(n2)
 *
 * byContribution: 0 for mu, 1 for tau, 2 for both mu and tau
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PR_XYTerm( 
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const BYTE* __restrict__ pEtaTable,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter, BYTE byMu, INT iTau,
    INT pathLdir1, INT pathLdir2, INT pathLdir3, BYTE Llength,
    INT pathRdir1, INT pathRdir2, INT pathRdir3, BYTE Rlength,
    BYTE byContribution)
{
    intokernalInt4;
    //const UINT uiBigIdx = __bi(sSite4);

    //=================================
    // 1. Find n1, n2
    INT Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    INT Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    SSmallInt4 site_n1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(_deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength))];
    //const SSmallInt4 middleSite = _deviceSmallInt4OffsetC(site_n1, byMu + 1);
    //From now on, site_n1 is smiddle
    site_n1 = _deviceSmallInt4OffsetC(site_n1, byMu + 1);
    const SIndex& smiddle = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];
    
    site_n1 = __deviceSiteIndexToInt4(smiddle.m_uiSiteIndex);
    //y Dx and -x Dy
    Real fNv = (0 == byMu)
        ? static_cast<Real>(site_n1.y - sCenter.y + F(0.5))
        : static_cast<Real>(sCenter.x - site_n1.x - F(0.5));

    //=================================
    // 2. Find V(n,n1), V(n,n2)
    const deviceSU3 vnn1 = _deviceLink(pGauge, sSite4, Llength, 1, Ldirs);
    const deviceSU3 vnn2 = _deviceLink(pGauge, sSite4, Rlength, 1, Rdirs);

    for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
    {
        const deviceSU3Vector* phi_i = pFermionPointers[rfieldId];
        const deviceSU3Vector* phi_id = pFermionPointers[rfieldId + uiRational];
        //=================================
        // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
        deviceSU3Vector phi1 = vnn1.MulVector(phi_id[sn1.m_uiSiteIndex]);
        deviceSU3Vector phi2 = vnn2.MulVector(phi_i[sn2.m_uiSiteIndex]);
        deviceSU3Vector phi3 = vnn1.MulVector(phi_i[sn1.m_uiSiteIndex]);
        deviceSU3Vector phi4 = vnn2.MulVector(phi_id[sn2.m_uiSiteIndex]);
        if (sn1.NeedToOpposite())
        {
            phi1.MulReal(F(-1.0));
            phi3.MulReal(F(-1.0));
        }
        if (sn2.NeedToOpposite())
        {
            phi2.MulReal(F(-1.0));
            phi4.MulReal(F(-1.0));
        }
        deviceSU3 res = deviceSU3::makeSU3ContractV(phi1, phi2);
        res.Add(deviceSU3::makeSU3ContractV(phi4, phi3));
        res.Ta();
        const Real eta_tau = (1 == ((pEtaTable[smiddle.m_uiSiteIndex] >> 3) & 1)) ? F(-1.0) : F(1.0);
        res.MulReal(OneOver12 * fOmega * fNv * pNumerators[rfieldId] * eta_tau);

        //For mu
        if (0 == byContribution || 2 == byContribution)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byMu);
            pForce[linkIndex].Sub(res);
        }

        //For tau
        if (1 == byContribution || 2 == byContribution)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);
            if (iTau > 0)
            {
                pForce[linkIndex].Sub(res);
            }
            else
            {
                pForce[linkIndex].Add(res);
            }
        }
    }
}

/**
 *
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PR_XYTau_Term(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    Real fOmega, 
    INT pathLdir1, INT pathLdir2, INT pathLdir3, BYTE Llength,
    INT pathRdir1, INT pathRdir2, INT pathRdir3, BYTE Rlength)
{
    intokernalInt4;
    //const UINT uiBigIdx = __bi(sSite4);

    //=================================
    // 1. Find n1, n2
    INT Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    INT Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];

    //Why use sn2? shouldn't it be sn1?
    const Real eta124 = _deviceEta124(__deviceSiteIndexToInt4(sn2.m_uiSiteIndex));
    //=================================
    // 2. Find V(n,n1), V(n,n2)
    const deviceSU3 vnn1 = _deviceLink(pGauge, sSite4, Llength, 1, Ldirs);
    const deviceSU3 vnn2 = _deviceLink(pGauge, sSite4, Rlength, 1, Rdirs);

    for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
    {
        const deviceSU3Vector* phi_i = pFermionPointers[rfieldId];
        const deviceSU3Vector* phi_id = pFermionPointers[rfieldId + uiRational];

        //=================================
        // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
        deviceSU3Vector phi1 = vnn1.MulVector(phi_id[sn1.m_uiSiteIndex]);
        deviceSU3Vector phi2 = vnn2.MulVector(phi_i[sn2.m_uiSiteIndex]);
        deviceSU3Vector phi3 = vnn1.MulVector(phi_i[sn1.m_uiSiteIndex]);
        deviceSU3Vector phi4 = vnn2.MulVector(phi_id[sn2.m_uiSiteIndex]);
        if (sn1.NeedToOpposite())
        {
            phi1.MulReal(F(-1.0));
            phi3.MulReal(F(-1.0));
        }
        if (sn2.NeedToOpposite())
        {
            phi2.MulReal(F(-1.0));
            phi4.MulReal(F(-1.0));
        }
        deviceSU3 res = deviceSU3::makeSU3ContractV(phi1, phi2);
        //This was phi2 phi1+ * eta124(n1) - phi3 phi4+ * eta124(n2)
        //The sign of the second term is because of 'dagger'
        //However, eta124(n1) = -eta124(n2), so use Add directly.
        res.Add(deviceSU3::makeSU3ContractV(phi4, phi3));
        res.Ta();
        res.MulReal(OneOver48 * fOmega * pNumerators[rfieldId] * eta124);

        if (pathLdir1 > 0)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLdir1 - 1);
            pForce[linkIndex].Add(res);
        }

        if (pathRdir1 > 0)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRdir1 - 1);
            pForce[linkIndex].Sub(res);
        }
    }

}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PR_XYTau_Term2(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    INT pathLdir1, INT pathLdir2, INT pathLdir3)
{
    intokernalInt4;
    const INT full[3] = { pathLdir1, pathLdir2, pathLdir3 };
    INT Ldirs[3];
    INT Rdirs[3];
    BYTE Llength, Rlength;

    #pragma unroll
    for (BYTE sep = 0; sep < 4; ++sep)
    {
        _deviceSeperate(full, sep, 3, Ldirs, Rdirs, Llength, Rlength);

        const UBOOL bHasLeft = Llength > 0 && Ldirs[0] > 0;
        const UBOOL bHasRight = Rlength > 0 && Rdirs[0] > 0;
        if (!bHasLeft && !bHasRight)
        {
            continue;
        }

        //=================================
        // 1. Find n1, n2
        const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
        const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength);
        const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
        const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];

        //Why use sn2? shouldn't it be sn1?
        const Real eta124 = _deviceEta124(__deviceSiteIndexToInt4(sn2.m_uiSiteIndex));
        //=================================
        // 2. Find V(n,n1), V(n,n2)
        const deviceSU3 vnn1 = _deviceLink(pGauge, sSite4, Llength, 1, Ldirs);
        const deviceSU3 vnn2 = _deviceLink(pGauge, sSite4, Rlength, 1, Rdirs);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceSU3Vector* phi_i = pFermionPointers[rfieldId];
            const deviceSU3Vector* phi_id = pFermionPointers[rfieldId + uiRational];

            //=================================
            // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
            deviceSU3Vector phi1 = vnn1.MulVector(phi_id[sn1.m_uiSiteIndex]);
            deviceSU3Vector phi2 = vnn2.MulVector(phi_i[sn2.m_uiSiteIndex]);
            deviceSU3Vector phi3 = vnn1.MulVector(phi_i[sn1.m_uiSiteIndex]);
            deviceSU3Vector phi4 = vnn2.MulVector(phi_id[sn2.m_uiSiteIndex]);
            if (sn1.NeedToOpposite())
            {
                phi1.MulReal(F(-1.0));
                phi3.MulReal(F(-1.0));
            }
            if (sn2.NeedToOpposite())
            {
                phi2.MulReal(F(-1.0));
                phi4.MulReal(F(-1.0));
            }
            deviceSU3 res = deviceSU3::makeSU3ContractV(phi1, phi2);
            //This was phi2 phi1+ * eta124(n1) - phi3 phi4+ * eta124(n2)
            //The sign of the second term is because of 'dagger'
            //However, eta124(n1) = -eta124(n2), so use Add directly.
            res.Add(deviceSU3::makeSU3ContractV(phi4, phi3));
            res.Ta();
            res.MulReal(OneOver48 * fOmega * pNumerators[rfieldId] * eta124);

            if (bHasLeft)
            {
                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, Ldirs[0] - 1);
                pForce[linkIndex].Add(res);
            }

            if (bHasRight)
            {
                const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, Rdirs[0] - 1);
                pForce[linkIndex].Sub(res);
            }
        }
    }
}

#pragma endregion


#pragma endregion

#pragma region D and derivate

void CFieldFermionKSSU3R::DOperatorKS(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    CFieldFermionKSSU3::DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);

    deviceSU3Vector* pTarget = (deviceSU3Vector*)pTargetBuffer;
    const deviceSU3Vector* pSource = (const deviceSU3Vector*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    preparethread;
    _kernelDFermionKS_PR_XYTerm << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        m_byFieldId,
        1,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
   
    _kernelDFermionKS_PR_XYTau_Term << <block, threads >> > (
        pSource,
        pGauge,
        pTarget,
        m_byFieldId,
        1,
        CCommonData::m_fOmega,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

void CFieldFermionKSSU3R::DerivateD0(
    void* pForce,
    const void* pGaugeBuffer) const
{
    CFieldFermionKSSU3::DerivateD0(pForce, pGaugeBuffer);

#if 1
    preparethread;
    #pragma region X Y Term

    INT mu[2] = { 0, 1 };
    for (INT imu = 0; imu < 2; ++imu)
    {
        INT dirs[6][3] =
        {
            {4, mu[imu] + 1, mu[imu] + 1},
            {mu[imu] + 1, 4, mu[imu] + 1},
            {mu[imu] + 1, mu[imu] + 1, 4},
            //{4, -mu[imu] - 1, -mu[imu] - 1},
            //{-mu[imu] - 1, 4, -mu[imu] - 1},
            //{-mu[imu] - 1, -mu[imu] - 1, 4},
            {mu[imu] + 1, mu[imu] + 1, -4},
            {mu[imu] + 1, -4, mu[imu] + 1},
            {-4, mu[imu] + 1, mu[imu] + 1},
        };

        INT iTau[6] = { 1, 1, 1, -1, -1, -1 };
        BYTE contributionOf[6][4] =
        {
            {1, 0, 0, 3},
            {0, 1, 0, 3},
            {0, 0, 1, 3},
            //{1, 3, 0, 0},
            //{3, 2, 3, 0},
            //{3, 0, 2, 3},
            {0, 0, 3, 1},
            {0, 3, 2, 3},
            {3, 2, 0, 3},
        };

        for (INT pathidx = 0; pathidx < 6; ++pathidx)
        {
            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                if (3 == contributionOf[pathidx][iSeperation])
                {
                    continue;
                }

                INT L[3] = { 0, 0, 0 };
                INT R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(dirs[pathidx], iSeperation, L, R, LLength, RLength);

                _kernelDFermionKSForce_PR_XYTerm << <block, threads >> > (
                    (const deviceSU3*)pGaugeBuffer,
                    (deviceSU3*)pForce,
                    appGetLattice()->m_pIndexCache->m_pEtaMu,
                    m_pRationalFieldPointers,
                    m_pMDNumerator,
                    m_rMD.m_uiDegree,
                    m_byFieldId,
                    CCommonData::m_fOmega, CCommonData::m_sCenter,
                    static_cast<BYTE>(imu), iTau[pathidx],
                    L[0], L[1], L[2], LLength,
                    R[0], R[1], R[2], RLength,
                    contributionOf[pathidx][iSeperation]
                    );
            }
        }
    }

    #pragma endregion

    #pragma region Polarization term

    //===========================
    //polarization terms
    //ilinkType is +-x +-y +t,
    //INT linkTypes[4][3] =
    //{
    //    {1, 2, 4},
    //    {1, 2, -4},
    //    {-1, 2, 4},
    //    {-1, 2, -4}
    //};
    INT linkTypes[4][3] =
    {
        {1, 2, 4},
        {1, -2, 4},
        {-1, 2, 4},
        {-1, -2, 4}
    };

    for (INT ilinkType = 0; ilinkType < 4; ++ilinkType)
    {
        INT sixlinks[6][3] =
        {
            {linkTypes[ilinkType][0], linkTypes[ilinkType][1], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][0], linkTypes[ilinkType][2], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][0], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][2], linkTypes[ilinkType][0]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][0], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][1], linkTypes[ilinkType][0]}
        };

        for (INT isixtype = 0; isixtype < 6; ++isixtype)
        {
            //bearly no change of time, because force calculation is not frequent
            _kernelDFermionKSForce_PR_XYTau_Term2 << <block, threads >> > (
                (const deviceSU3*)pGaugeBuffer,
                (deviceSU3*)pForce,
                m_pRationalFieldPointers,
                m_pMDNumerator,
                m_rMD.m_uiDegree,
                m_byFieldId,
                CCommonData::m_fOmega,
                sixlinks[isixtype][0], sixlinks[isixtype][1], sixlinks[isixtype][2]
                );
            /*
            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                INT L[3] = { 0, 0, 0 };
                INT R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(sixlinks[isixtype], iSeperation, L, R, LLength, RLength);

                const UBOOL bHasLeft = (LLength > 0) && (L[0] > 0);
                const UBOOL bHasRight = (RLength > 0) && (R[0] > 0);

                if (bHasLeft || bHasRight)
                {
                    _kernelDFermionKSForce_PR_XYTau_Term << <block, threads >> > (
                        (const deviceSU3*)pGaugeBuffer,
                        (deviceSU3*)pForce,
                        m_pRationalFieldPointers,
                        m_pMDNumerator,
                        m_rMD.m_uiDegree,
                        m_byFieldId,
                        CCommonData::m_fOmega,
                        L[0], L[1], L[2], LLength,
                        R[0], R[1], R[2], RLength
                        );
                }
            }
            */
        }
    }
    
    #pragma endregion
#endif
}

#pragma endregion

void CFieldFermionKSSU3R::InitialOtherParameters(CParameters& params)
{
    CFieldFermionKSSU3::InitialOtherParameters(params);
    m_bEachSiteEta = TRUE;
}

void CFieldFermionKSSU3R::CopyTo(CField* U) const
{
    CFieldFermionKSSU3::CopyTo(U);
}

CCString CFieldFermionKSSU3R::GetInfos(const CCString& tab) const
{
    CCString sRet = tab + _T("Name : CFieldFermionKSSU3R\n");
    sRet = sRet + tab + _T("Mass (2am) : ") + appFloatToString(m_f2am) + _T("\n");
    sRet = sRet + tab + _T("MD Rational (c) : ") + appFloatToString(m_rMD.m_fC) + _T("\n");
    sRet = sRet + tab + _T("MC Rational (c) : ") + appFloatToString(m_rMC.m_fC) + _T("\n");
    sRet = sRet + tab + _T("Omega : ") + appFloatToString(CCommonData::m_fOmega) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================