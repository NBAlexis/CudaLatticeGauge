//=============================================================================
// FILENAME : CFieldFermionKSSU3Acc.h
// 
// DESCRIPTION:
// 
// This is only suitable when gt >> 1, this uses Galilean transform
// 
// gammai (pi + iAi) + gt gammaz (pt+iAt) + m
// 
//
// REVISION:
//  [11/21/2023 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONKSSU3ACC_H_
#define _CFIELDFERMIONKSSU3ACC_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3Acc)

class CLGAPI CFieldFermionKSSU3Acc : public CFieldFermionKSSU3
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3Acc)

public:

    void DerivateD0(void* pForce, const void* pGaugeBuffer) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

    CCString GetInfos(const CCString& tab) const override;
};

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

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3DACC_H_

//=============================================================================
// END OF FILE
//=============================================================================