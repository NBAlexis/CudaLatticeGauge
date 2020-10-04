//=============================================================================
// FILENAME : CFieldGaugeSU3TreeImproved.h
// 
// DESCRIPTION:
// Note:
// This was beta * [(5/3) plaq - (1/12) rect] with beta=6/g^2
// Use beta = 10/g^2 instead of 6/g^2,
// this is then beta * [plaq + (-1/20) rect]
//
// Note:
// Whenever staple is mentioned, we mean the staple in the sense of plaquettes only
// This is also used to define a 'fat' link
// So the rectangle terms are not included
//
// REVISION:
//  [10/03/2020 nbale]
//=============================================================================

#ifndef _CFIELDGAUGE_SU3_TREEIMPROVED_H_
#define _CFIELDGAUGE_SU3_TREEIMPROVED_H_


__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU3TreeImproved)

class CLGAPI CFieldGaugeSU3TreeImproved : public CFieldGaugeSU3
{
    __CLGDECLARE_FIELD(CFieldGaugeSU3TreeImproved)

public:
    CFieldGaugeSU3TreeImproved();
    ~CFieldGaugeSU3TreeImproved();

    void InitialOtherParameters(CParameters& param) override;

#pragma region HMC

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override;

    Real CalculatePlaqutteEnergy(Real betaOverN) const override;
    Real CalculatePlaqutteEnergyUseClover(Real betaOverN) const override;
    Real CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge *pStaple) const override;

    Real m_fRectOverPlaq;

#pragma endregion

};

#pragma region device functions

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

__END_NAMESPACE


#endif //#ifndef _CFIELDGAUGE_SU3_TREEIMPROVED_H_

//=============================================================================
// END OF FILE
//=============================================================================