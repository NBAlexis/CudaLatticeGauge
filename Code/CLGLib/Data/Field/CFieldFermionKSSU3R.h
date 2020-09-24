//=============================================================================
// FILENAME : CFieldFermionKSSU3R.h
// 
// DESCRIPTION:
// This is the class for Kogut-Susskind staggered fermions
// For pseudo fermion, this is in fact a boson field phi.
//
// Current implementation, assumes square lattice
//
// REVISION:
//  [09/23/2020 nbale]
//=============================================================================

#ifndef _CFIELDFERMIONKSSU3R_H_
#define _CFIELDFERMIONKSSU3R_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldFermionKSSU3R)

class CLGAPI CFieldFermionKSSU3R : public CFieldFermionKSSU3
{
    __CLGDECLARE_FIELD(CFieldFermionKSSU3R)

public:

    void DerivateD0(void* pForce, const void* pGaugeBuffer) const override;
    void DOperatorKS(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, Real f2am,
        UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const override;

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
static __device__ __inline__ deviceSU3 _deviceVXXTau(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    UINT bXorY, UBOOL bPlusMu, UBOOL bPlusTau)
{
    const INT iMu = bXorY ? (bPlusMu ? 1 : -1) : (bPlusMu ? 2 : -2);
    const INT iTau = bPlusTau ? 4 : -4;
    INT dir1[3];

    dir1[0] = iMu;
    dir1[1] = iMu;
    dir1[2] = iTau;
    deviceSU3 sRet = _deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1);
    
    dir1[0] = iMu;
    dir1[1] = iTau;
    dir1[2] = iMu;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iTau;
    dir1[1] = iMu;
    dir1[2] = iMu;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    sRet.MulReal(OneOver3);
    return sRet;
}

static __device__ __inline__ Real _deviceEta124(const SSmallInt4& sSite)
{
    return (((sSite.y + sSite.z) & 1) > 0) ? (F(-1.0)) : (F(1.0));
}

static __device__ __inline__ deviceSU3 _deviceVXYT(
    const deviceSU3* __restrict__ pDeviceData,
    const SSmallInt4& sStartSite, BYTE byFieldId,
    UBOOL bPlusX, UBOOL bPlusY, UBOOL bPlusTau)
{
    const INT iX = bPlusX ? 1 : -1;
    const INT iY = bPlusY ? 2 : -2;
    const INT iT = bPlusTau ? 4 : -4;
    INT dir1[3];

    dir1[0] = iX; dir1[1] = iY; dir1[2] = iT;
    deviceSU3 sRet(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iX; dir1[1] = iT; dir1[2] = iY;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iY; dir1[1] = iX; dir1[2] = iT;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iY; dir1[1] = iT; dir1[2] = iX;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iT; dir1[1] = iX; dir1[2] = iY;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    dir1[0] = iT; dir1[1] = iY; dir1[2] = iX;
    sRet.Add(_deviceLink(pDeviceData, sStartSite, 3, byFieldId, dir1));

    sRet.MulReal(OneOver6);
    return sRet;
}


#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDFERMIONKSSU3R_H_

//=============================================================================
// END OF FILE
//=============================================================================