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
#include "CFieldFermionKSU1.h"

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


#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSU1R_H_

//=============================================================================
// END OF FILE
//=============================================================================