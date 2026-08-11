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
#include "DeviceTemplates/DeviceInlineUseNoTemplateFunction.h"
#include "Tools/Math/DeviceInlineTemplate.h"

#ifndef _DEVICEINLINESU3_H_
#define _DEVICEINLINESU3_H_

__BEGIN_NAMESPACE

#pragma region device functions

/**
* Note: for baked plaqutte index, the bond if is set to SIndex
* If it is a "new SIndex" instead, remember to set the m_byTag
*/
//static __device__ __inline__ const deviceSU3& _deviceGetGaugeBCSU3(
//    BYTE byFieldId,
//    const deviceSU3* __restrict__ pBuffer,
//    const SIndex& idx)
//{
//    return idx.IsDirichlet() ?
//        ((CFieldBoundaryGaugeSU3*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
//            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
//        ]
//        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
//}

/**
* If the bond is on surface, return the Dirichlet
* else, return the element
*/
//static __device__ __inline__ const deviceSU3& _deviceGetGaugeBCSU3Dir(
//    BYTE byFieldId,
//    const deviceSU3* __restrict__ pBuffer,
//    UINT uiBigIdx,
//    BYTE byDir)
//{
//    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
//    return __idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, byDir) ?
//        ((CFieldBoundaryGaugeSU3*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
//            __idx->_devcieExchangeBoundaryFieldSiteIndex(site) * _DC_Dir + byDir
//        ]
//        : pBuffer[_deviceGetLinkIndex(site.m_uiSiteIndex, byDir)];
//}

//static __device__ __inline__ deviceSU3 _deviceGetGaugeBCSU3DirOne(
//    BYTE byFieldId,
//    const deviceSU3* __restrict__ pBuffer,
//    UINT uiBigIdx,
//    BYTE byDir)
//{
//    return __idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, byDir) ?
//        deviceSU3::makeSU3Id()
//        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
//}

//static __device__ __inline__ deviceSU3 _deviceGetGaugeBCSU3DirZero(
//    BYTE byFieldId,
//    const deviceSU3* __restrict__ pBuffer,
//    UINT uiBigIdx,
//    BYTE byDir)
//{
//    return __idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, byDir) ?
//        deviceSU3::makeSU3Zero()
//        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
//}


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

    deviceSU3 res = _deviceGetGaugeBCDirZeroT(byFieldId, piApure, uiBigIdx, byMu); //Apure _mu
    deviceSU3 res2 = _deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu); //A _nu
    res2.Mul(res); //A _nu Apure _mu
    res.Mul(_deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu)); //Apure _mu A _nu
    res.Sub(res2); //[Apure, A]
    res.Add(_deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu));
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

    deviceSU3 res = _deviceGetGaugeBCDirZeroT(byFieldId, piApure, uiBigIdx, byMu); //Apure _mu
    deviceSU3 res2 = _deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu); //A _nu
    res2.Mul(res); //A _nu Apure _mu
    res.Mul(_deviceGetGaugeBCDirZeroT(byFieldId, piA, uiBigIdx, byNu)); //Apure _mu A _nu
    res.Sub(res2); //[Apure, A]
    res.Add(_deviceGetGaugeBCSU3DirZeroSIndex(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_p_mu * _DC_Dir + byNu]).MulRealC(F(0.5)));
    res.Sub(_deviceGetGaugeBCSU3DirZeroSIndex(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_m_mu * _DC_Dir + byNu]).MulRealC(F(0.5)));
    return res;
}

#pragma endregion

//================= Put those device functions to header file because we will use them ==============

#pragma region KS SU3 Acc

/**
* gamma_mu partial_nu term
* used for test, when test is done, we use _deviceGMUPNUOptimized
*/
//static __device__ __inline__ deviceSU3 _deviceGMUPNU(
//    const deviceSU3* __restrict__ pDeviceData,
//    const SSmallInt4& sStartSite, BYTE byFieldId,
//    BYTE mu, BYTE nu, UBOOL bPlusMu, UBOOL bPlusNu)
//{
//    const INT iMu = bPlusMu ? (mu + 1) : (-mu - 1);
//    const INT iNu = bPlusNu ? (nu + 1) : (-nu - 1);
//    INT dir1[3];
//
//    dir1[0] = iNu;
//    dir1[1] = iNu;
//    dir1[2] = iMu;
//    deviceSU3 sRet = _deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1);
//
//    dir1[0] = iNu;
//    dir1[1] = iMu;
//    dir1[2] = iNu;
//    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));
//
//    dir1[0] = iMu;
//    dir1[1] = iNu;
//    dir1[2] = iNu;
//    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));
//
//    sRet.MulReal(OneOver3);
//    return sRet;
//}

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
    const SCHAR iMu = bPlusMu ? (mu + 1) : (-static_cast<SCHAR>(mu) - 1);
    const SCHAR iNu = bPlusNu ? (nu + 1) : (-static_cast<SCHAR>(nu) - 1);
    SCHAR dir1[3];

    SSmallInt4 x_p_nu_p_mu = sStartSite;
    x_p_nu_p_mu.m_byData4[mu] = x_p_nu_p_mu.m_byData4[mu] + (bPlusMu ? 1 : -1);
    x_p_nu_p_mu.m_byData4[nu] = x_p_nu_p_mu.m_byData4[nu] + (bPlusNu ? 1 : -2);

    dir1[0] = iMu; dir1[1] = iNu;
    deviceSU3 sRet = _deviceLinkT(pDeviceData, sStartSite, 2, byFieldId, dir1);
    dir1[0] = iNu; dir1[1] = iMu;
    sRet.Add(_deviceLinkT(pDeviceData, sStartSite, 2, byFieldId, dir1));

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

    sRet.Add(_deviceLinkT(pDeviceData, sStartSite, 3, byFieldId, dir1));

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
        _deviceGetSTTermT(byFieldId, pDeviceData, n__nu, n_p_nu__mu, n_p_mu__nu)
    );

    left.Add(
        //    _deviceGetSTTerm(
        //    pDeviceData, uiN_m_nu, uiN_m_nu, uiN_p_mu_m_nu, nu, mu, nu, 1, 0, 0
        //)
        _deviceGetSTTermT(byFieldId, pDeviceData, n_m_nu__nu, n_m_nu__mu, n_p_mu_m_nu__nu)
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

    deviceSU3 ret(_deviceS1T(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
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

    deviceSU3 ret(_deviceS2T(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
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
    ret.Mul(_deviceS3T(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

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
    ret.Mul(_deviceS4T(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));

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
    ret.Mul(_deviceT1T(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
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
    ret.DaggerMul(_deviceT2T(byFieldId, pDeviceData, sSite, uiBigIndex, mu, nu, rho));
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

__END_NAMESPACE

#endif //#ifndef _DEVICEINLINESU3_H_

//=============================================================================
// END OF FILE
//=============================================================================
