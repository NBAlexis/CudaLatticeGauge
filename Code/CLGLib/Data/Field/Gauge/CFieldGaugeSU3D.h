//=============================================================================
// FILENAME : CFieldGaugeSU3D.h
// 
// DESCRIPTION:
// There is simplifications for periodic boundary condition
// which is invalid for Dirichlet.
// This is only for Dirichlet.
//
// REVISION:
//  [05/17/2019 nbale]
//=============================================================================

#ifndef _CFIELDGAUGE_SU3D_H_
#define _CFIELDGAUGE_SU3D_H_

#define gaugeSU3KernelFuncionStart \
    intokernaldir; \
    for (UINT idir = 0; idir < uiDir; ++idir) \
    { \
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir); 


#define gaugeSU3KernelFuncionEnd \
    } 



__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU3D)

class CLGAPI CFieldGaugeSU3D : public CFieldGaugeSU3
{
    __CLGDECLARE_FIELD(CFieldGaugeSU3D)

public:

    CFieldGaugeSU3D() : CFieldGaugeSU3() {}

#pragma region HMC

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override;
    void CalculateOnlyStaple(CFieldGauge* pStaple) const override;

    void MakeRandomGenerator() override;
#if !_CLG_DOUBLEFLOAT
    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override;
    DOUBLE CalculateKinematicEnergy() const override;
    DOUBLE CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge* pStaple) const override
    {
        return CalculatePlaqutteEnergy(betaOverN);
    }
#else
    Real CalculatePlaqutteEnergy(Real betaOverN) const override;
    Real CalculateKinematicEnergy() const override;
    Real CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge* pStaple) const override
    {
        return CalculatePlaqutteEnergy(betaOverN);
    }
#endif

#pragma endregion

#pragma region BLAS

    void FixBoundary() override;

#pragma endregion

    void ExpMult(Real a, CField* U) const override;
    CCString GetInfos(const CCString &tab) const override;

#pragma region Test Functions to test gauge invarience of angular momentum

    /**
     * iA = U.TA() / 2
     */
    void TransformToIA() override;

    /**
     * U=exp(iA)
     */
    void TransformToU() override;

#pragma endregion

};

#pragma region device functions

/**
* Note: for baked plaqutte index, the bond if is set to SIndex
* If it is a "new SIndex" instead, remember to set the m_byTag
*/
static __device__ __inline__ const deviceSU3& _deviceGetGaugeBCSU3(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pBuffer,
    const SIndex& idx)
{
    return idx.IsDirichlet() ?
        ((CFieldBoundaryGaugeSU3*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

/**
* If the bond is on surface, return the Dirichlet
* else, return the element
*/
static __device__ __inline__ const deviceSU3& _deviceGetGaugeBCSU3Dir(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        ((CFieldBoundaryGaugeSU3*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(site) * _DC_Dir + byDir
        ]
        : pBuffer[_deviceGetLinkIndex(site.m_uiSiteIndex, byDir)];
}

static __device__ __inline__ deviceSU3 _deviceGetGaugeBCSU3DirOne(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        deviceSU3::makeSU3Id()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
}

static __device__ __inline__ deviceSU3 _deviceGetGaugeBCSU3DirZero(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        deviceSU3::makeSU3Zero()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
}


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

    deviceSU3 res = _deviceGetGaugeBCSU3DirZero(byFieldId, piApure, uiBigIdx, byMu); //Apure _mu
    deviceSU3 res2 = _deviceGetGaugeBCSU3DirZero(byFieldId, piA, uiBigIdx, byNu); //A _nu
    res2.Mul(res); //A _nu Apure _mu
    res.Mul(_deviceGetGaugeBCSU3DirZero(byFieldId, piA, uiBigIdx, byNu)); //Apure _mu A _nu
    res.Sub(res2); //[Apure, A]
    res.Add(_deviceGetGaugeBCSU3DirZero(byFieldId, piA, uiBigIdx, byNu));
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

    deviceSU3 res = _deviceGetGaugeBCSU3DirZero(byFieldId, piApure, uiBigIdx, byMu); //Apure _mu
    deviceSU3 res2 = _deviceGetGaugeBCSU3DirZero(byFieldId, piA, uiBigIdx, byNu); //A _nu
    res2.Mul(res); //A _nu Apure _mu
    res.Mul(_deviceGetGaugeBCSU3DirZero(byFieldId, piA, uiBigIdx, byNu)); //Apure _mu A _nu
    res.Sub(res2); //[Apure, A]
    res.Add(_deviceGetGaugeBCSU3DirZeroSIndex(piA, 
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_p_mu * _DC_Dir + byNu]).MulRealC(F(0.5)));
    res.Sub(_deviceGetGaugeBCSU3DirZeroSIndex(piA, 
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_m_mu * _DC_Dir + byNu]).MulRealC(F(0.5)));
    return res;
}


static __device__ __inline__ deviceSU3 _devicePlaqutte(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
    UINT uiSiteIndex,
    BYTE plaqIdx, //0-5, as 12, 13, 14, 23, 24, 34
    BYTE plaqLength, //Always 4
    BYTE plaqCountAll //Always 24
    )
{
    SIndex first = pCachedPlaqutte[plaqIdx * plaqLength + uiSiteIndex * plaqCountAll];
    deviceSU3 toAdd(_deviceGetGaugeBCSU3(byFieldId, pDeviceData, first));
    if (first.NeedToDagger())
    {
        toAdd.Dagger();
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[plaqIdx * plaqLength + j + uiSiteIndex * plaqCountAll];
        deviceSU3 toMul(_deviceGetGaugeBCSU3(byFieldId, pDeviceData, first));
        if (first.NeedToDagger())
        {
            toAdd.MulDagger(toMul);
        }
        else
        {
            toAdd.Mul(toMul);
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
 */
static __device__ __inline__ deviceSU3 _deviceLink(
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceSU3 sRet = deviceSU3::makeSU3Id();
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ? 
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if ( (newLink.NeedToDagger() && !bDagger)
                 || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet.Dagger();
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                if ((newLink.NeedToDagger() && !bDagger)
                || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet.MulDagger(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    sRet.Mul(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
            }
        }

        if (pDir[i] > 0 && i < (byLength - 1)) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return sRet;
}

/**
 * After every move, it maps to inside the lattice
 * Do NOT use it in projective plane boundary condition
 */
static __device__ __inline__ deviceSU3 _deviceLinkLong(
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceSU3 sRet = deviceSU3::makeSU3Id();
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(sStartSite) + byDir];
        sStartSite = __deviceSiteIndexToInt4(newLink.m_uiSiteIndex);

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet.Dagger();
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                if ((newLink.NeedToDagger() && !bDagger)
                    || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet.MulDagger(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    sRet.Mul(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
            }
        }

        if (pDir[i] > 0 && i < (byLength - 1)) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return sRet;
}

#if discard
/**
 * Get magnetic phase
 * type0:
 * Ay(n) = q B nx
 * if twisted Ax(Lx) = - q B ny
 *
 * type1:
 * Ax(n) = - q B ny
 * if twisted Ay(Ly) = q B nx
 *
 * type2:
 * Ax(n) = - q B ny/2
 * Ay(n) = qB nx / 2
 *
 * if twisted
 * Ax(Lx) = - q B ny
 * Ay(Ly) = q B nx
 *
 */
static __device__ __inline__ Real __deviceGetMagneticPhase(const SSmallInt4& sSite, const SSmallInt4& sCenter, BYTE byLink, BYTE byGaugeType, UBOOL bTwisted)
{
    if (0 != byLink && 1 != byLink)
    {
        return F(0.0);
    }

    if (0 == byGaugeType)
    {
        if (1 == byLink)
        {
            return static_cast<Real>(sSite.x - sCenter.x + F(0.5));
        }
        if (bTwisted && sSite.x == _DC_Lxi - 1)
        {
            return -static_cast<Real>(sSite.y - sCenter.y + F(0.5));
        }
        return F(0.0);
    }

    if (1 == byGaugeType)
    {
        if (0 == byLink)
        {
            return -static_cast<Real>(sSite.y - sCenter.y + F(0.5));
        }
        if (bTwisted && sSite.y == _DC_Lyi - 1)
        {
            return static_cast<Real>(sSite.x - sCenter.x + F(0.5));
        }
        return F(0.0);
    }

    if (0 == byLink)
    {
        if (bTwisted && sSite.x == _DC_Lxi - 1)
        {
            return -static_cast<Real>(sSite.y - sCenter.y + F(0.5));
        }
        return -static_cast<Real>(sSite.y - sCenter.y + F(0.5)) * F(0.5);
    }

    if (bTwisted && sSite.x == _DC_Lxi - 1)
    {
        return static_cast<Real>(sSite.x - sCenter.x + F(0.5));
    }
    return static_cast<Real>(sSite.x - sCenter.x + F(0.5)) * F(0.5);
}

/**
 * Same as _deviceLink, but with fixed electric-magnetic field
 * The electric-magnetic field is parameterized as Bz
 * For magnetic field between n1 and n2 is:
 * Ax(Lx) = - q B ny
 * Ay(n) = q B nx
 *
 * If no twisted boundary, then only Ay(n) = q B nx, others are zero
 */
static __device__ __inline__ deviceSU3 _deviceLinkMP(
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, const SSmallInt4& sCenterSite,
    BYTE byLength, BYTE byFieldId, Real fQBz, BYTE byGaugeType, UBOOL bTwistedBoundary,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceSU3 sRet = deviceSU3::makeSU3Id();
    Real fPhase = F(0.0);
    
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == pDir[i])
        {
            continue;
        }
        UBOOL bDagger = FALSE;
        const BYTE byDir = pDir[i] > 0 ?
            static_cast<BYTE>(pDir[i] - 1) : static_cast<BYTE>(-pDir[i] - 1);

        if (pDir[i] < 0) //Move-back
        {
            bDagger = TRUE;
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
        const UINT uiBi4StartSite = __bi4(sStartSite);
        const SIndex& newLink = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBi4StartSite + byDir];
        const UBOOL bDaggerFinal = 
               (newLink.NeedToDagger() && !bDagger)
            || (!newLink.NeedToDagger() && bDagger);
        if (!newLink.IsDirichlet())
        {
            if (0 == newLink.m_byDir || 1 == newLink.m_byDir)
            {
                const SSmallInt4 newSite = __deviceSiteIndexToInt4(newLink.m_uiSiteIndex);
                if (bDaggerFinal)
                {
                    fPhase = fPhase - __deviceGetMagneticPhase(newSite, sCenterSite, newLink.m_byDir, byGaugeType, bTwistedBoundary);
                }
                else
                {
                    fPhase = fPhase + __deviceGetMagneticPhase(newSite, sCenterSite, newLink.m_byDir, byGaugeType, bTwistedBoundary);
                }
            }
        }

        if (0 == i)
        {
            if (!newLink.IsDirichlet())
            {
                sRet = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                if (bDaggerFinal)
                {
                    sRet.Dagger();
                }
            }
        }
        else
        {
            if (!newLink.IsDirichlet())
            {
                if (bDaggerFinal)
                {
                    sRet.MulDagger(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    sRet.Mul(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
            }
        }

        if (pDir[i] > 0 && i < (byLength - 1)) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    fPhase = fPhase * fQBz;
    const CLGComplex cmpPhase = _make_cuComplex(_cos(fPhase), _sin(fPhase));
    sRet.MulComp(cmpPhase);
    return sRet;
}

/**
 * Phase of the link with:
 * Ax(Lx) = - q B ny
 * Ay(n) = q B nx
 */
static __device__ __inline__ Real _devicePhaseM(const SSmallInt4& sPos, BYTE byDir, const SSmallInt4& sCenterSite, Real fQBz)
{
    if (1 == byDir) //y-dir
    {
        return static_cast<Real>(sPos.x - sCenterSite.x + F(0.5)) * fQBz;
    }
    if (0 == byDir) //x-dir
    {
        if (sPos.x == static_cast<INT>(_DC_Lx) - 1)
        {
            return static_cast<Real>(sPos.y - sCenterSite.y + F(0.5)) * (-fQBz);
        }
    }

    return F(0.0);
}

/**
 * This function assumes _DC_Dir = 4
 */
static __device__ __inline__ Real _devicePhaseM(UINT uiLinkIndex, const SSmallInt4& sCenterSite, Real fQBz)
{
    return _devicePhaseM(__deviceLinkIndexToInt4(uiLinkIndex), static_cast<BYTE>(uiLinkIndex & 3), sCenterSite, fQBz);
}
#endif

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
static __device__ __inline__ deviceSU3 _device1PlaqutteTermPP(
    const deviceSU3* __restrict__ pDeviceData,
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

    deviceSU3 u(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_mu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_p_mu_nu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_p_nu_mu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_nu, byFieldId));

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
static __device__ __inline__ deviceSU3 _device1PlaqutteTermMP(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1);
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, byNu + 1);
    const UINT uin_m_mub4 = __idx->_deviceGetBigIndex(n_m_mu) * _DC_Dir;
    const SIndex& s_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byMu];
    const SIndex& s_m_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byNu];
    const SIndex& s_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_mu_p_nu) * _DC_Dir + byMu];
    const SIndex& s__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byNu];

    deviceSU3 u(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_mu__mu, byFieldId));
    u.DaggerMul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_mu__nu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_mu_p_nu__mu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s__nu, byFieldId));

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
static __device__ __inline__ deviceSU3 _device1PlaqutteTermPM(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byNu) - 1);
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(n_m_nu, byMu + 1);
    const UINT uin_m_nub4 = __idx->_deviceGetBigIndex(n_m_nu) * _DC_Dir;
    const SIndex& s_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byMu];
    const SIndex& s_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byNu];
    const SIndex& s_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_nu_p_mu) * _DC_Dir + byNu];
    const SIndex& s__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byMu];

    deviceSU3 u(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s__mu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_nu_p_mu__nu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_nu__mu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_nu__nu, byFieldId));

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
static __device__ __inline__ deviceSU3 _device1PlaqutteTermMM(
    const deviceSU3* __restrict__ pDeviceData,
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
    deviceSU3 u(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_nu_m_mu__nu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_mu__mu, byFieldId));
    u.DaggerMul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_nu_m_mu__mu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU3DirSIndex(pDeviceData, s_m_nu__nu, byFieldId));

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
static __device__ __inline__ deviceSU3 _deviceClover(const deviceSU3* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    deviceSU3 ret(_device1PlaqutteTermPP(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId));
    ret.Add(_device1PlaqutteTermMM(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId));
    ret.Add(_device1PlaqutteTermPM(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));
    ret.Add(_device1PlaqutteTermMP(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));

    return ret;
}

/**
 * Avoid the add of matrices
 */
static __device__ __inline__ Real _deviceCloverRetr(const deviceSU3* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    return _device1PlaqutteTermPP(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId).ReTr()
         + _device1PlaqutteTermMM(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId).ReTr()
         + _device1PlaqutteTermPM(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId).ReTr()
         + _device1PlaqutteTermMP(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId).ReTr();
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_SU3D_H_

//=============================================================================
// END OF FILE
//=============================================================================