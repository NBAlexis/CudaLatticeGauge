//=============================================================================
// FILENAME : CFieldGaugeSU2.h
// 
// DESCRIPTION:
// This is the class for the gauge fields
//
// REVISION:
//  [07/02/2024 nbale]
//=============================================================================

#ifndef _CFIELDGAUGE_SU2_H_
#define _CFIELDGAUGE_SU2_H_

__BEGIN_NAMESPACE

#if 0

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU2)

class CLGAPI CFieldGaugeSU2 : public CFieldGauge
{
    __CLGDECLARE_FIELD(CFieldGaugeSU2)
public:

    CFieldGaugeSU2();
    ~CFieldGaugeSU2();
    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType) override;
    void InitialWithByte(BYTE* byData) override;
    void InitialWithByteCompressed(BYTE*) override { appCrucial(_T("Not supported by %s\n"), __ENUM_TO_STRING(EFieldType, GetFieldType()).c_str()); }
    void InitialField(EFieldInitialType eInitialType) override;
    EFieldType GetFieldType() const override { return EFT_GaugeSU2; }
    UINT MatrixN() const override { return 2; }
    void DebugPrintMe() const override;
    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override;
    void CalculateOnlyStaple(CFieldGauge* pStaple) const override;
    void MakeRandomGenerator() override;
    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override;
    DOUBLE CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const override;
    DOUBLE CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge* pStaple) const override;
    DOUBLE CalculateKinematicEnergy() const override;
    void Zero() override;
    void Identity() override;
    void Dagger() override;
    void AxpyPlus(const CField* x) override;
    void AxpyMinus(const CField* x) override;
    void Axpy(Real a, const CField* x) override;
    void Axpy(const CLGComplex& a, const CField* x) override;
    void Mul(const CField* other, UBOOL bDagger = TRUE) override;
    void ScalarMultply(const CLGComplex& a) override;
    void ScalarMultply(Real a) override;
    void SetOneDirectionUnity(BYTE byDir) override;
    void SetOneDirectionZero(BYTE byDir) override;
    void TransformToIA() override { appCrucial(_T("Not supported by %s\n"), __ENUM_TO_STRING(EFieldType, GetFieldType()).c_str()); }
    void TransformToU() override { appCrucial(_T("Not supported by %s\n"), __ENUM_TO_STRING(EFieldType, GetFieldType()).c_str()); }
    void CalculateE_Using_U(CFieldGauge* pResoult) const override { appCrucial(_T("Not supported by %s\n"), __ENUM_TO_STRING(EFieldType, GetFieldType()).c_str()); }
    void CalculateNablaE_Using_U(CFieldGauge* pResoult, UBOOL bNaive = FALSE) const override { appCrucial(_T("Not supported by %s\n"), __ENUM_TO_STRING(EFieldType, GetFieldType()).c_str()); }
    void ExpMult(Real a, CField* U) const override;
    void ElementNormalize() override;
    cuDoubleComplex Dot(const CField* other) const override;
    CCString SaveToCompressedFile(const CCString& fileName) const override { appCrucial(_T("Not supported by %s\n"), __ENUM_TO_STRING(EFieldType, GetFieldType()).c_str()); return _T(""); }
    BYTE* CopyDataOut(UINT& uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;

    void PolyakovOnSpatialSite(cuDoubleComplex* buffer) const override;

    deviceSU2* m_pDeviceData;

    _GetData
};

#pragma region device functions

/**
* Note: for baked plaqutte index, the bond if is set to SIndex
* If it is a "new SIndex" instead, remember to set the m_byTag
*/
static __device__ __inline__ const deviceSU2& _deviceGetGaugeBCSU2(
    BYTE byFieldId,
    const deviceSU2* __restrict__ pBuffer,
    const SIndex& idx)
{
    return idx.IsDirichlet() ?
        ((CFieldBoundaryGaugeSU2*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

/**
* If the bond is on surface, return the Dirichlet
* else, return the element
*/
static __device__ __inline__ const deviceSU2& _deviceGetGaugeBCSU2Dir(
    BYTE byFieldId,
    const deviceSU2* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        ((CFieldBoundaryGaugeSU2*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(site) * _DC_Dir + byDir
        ] 
        : pBuffer[_deviceGetLinkIndex(site.m_uiSiteIndex, byDir)];
}

static __device__ __inline__ deviceSU2 _deviceGetGaugeBCSU2DirOne(
    BYTE byFieldId,
    const deviceSU2* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        deviceSU2::makeSU2Id()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
}

static __device__ __inline__ deviceSU2 _deviceGetGaugeBCSU2DirZero(
    BYTE byFieldId,
    const deviceSU2* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        deviceSU2::makeSU2Zero()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
}

static __device__ __inline__ deviceSU2 _deviceGetGaugeBCSU2DirSIndex(
    const deviceSU2* __restrict__ pBuffer,
    const SIndex& idx,
    BYTE byFieldId)
{
    deviceSU2 ret = idx.IsDirichlet() ?
        ((CFieldBoundaryGaugeSU2*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
    if (idx.NeedToDagger())
    {
        ret.Dagger();
    }
    return ret;
}

static __device__ __inline__ deviceSU2 _deviceGetGaugeBCSU2DirOneSIndex(
    const deviceSU2* __restrict__ pBuffer,
    const SIndex& idx)
{
    if (idx.IsDirichlet())
    {
        return deviceSU2::makeSU2Id();
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
static __device__ __inline__ deviceSU2 _deviceGetGaugeBCSU2DirZeroSIndex(
    const deviceSU2* __restrict__ pBuffer,
    const SIndex& idx)
{
    if (idx.IsDirichlet())
    {
        return deviceSU2::makeSU2Zero();
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
static __device__ __inline__ deviceSU2 _deviceDPureMuSU2(
    const deviceSU2* __restrict__ piA,
    const deviceSU2* __restrict__ piApure,
    const SSmallInt4& sSite4,
    UINT uiBigIdx,
    BYTE byMu,
    BYTE byNu,
    BYTE byFieldId)
{
    //i a D A = (A_nu (n) - A_nu (n-mu)) + iApure _mu A _nu - i A _nu Apure _mu
    const UINT uiSiteBig_m_mu = __idx->_deviceGetBigIndex(
        _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1));

    deviceSU2 res = _deviceGetGaugeBCSU2DirZero(byFieldId, piApure, uiBigIdx, byMu); //Apure _mu
    deviceSU2 res2 = _deviceGetGaugeBCSU2DirZero(byFieldId, piA, uiBigIdx, byNu); //A _nu
    res2.Mul(res); //A _nu Apure _mu
    res.Mul(_deviceGetGaugeBCSU2DirZero(byFieldId, piA, uiBigIdx, byNu)); //Apure _mu A _nu
    res.Sub(res2); //[Apure, A]
    res.Add(_deviceGetGaugeBCSU2DirZero(byFieldId, piA, uiBigIdx, byNu));
    res.Sub(_deviceGetGaugeBCSU2DirZeroSIndex(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_m_mu * _DC_Dir + byNu]));
    return res;
}

/**
 * test using (A(N+mu)-A(N-mu))/2
 */
static __device__ __inline__ deviceSU2 _deviceDPureMu2SU2(
    const deviceSU2* __restrict__ piA,
    const deviceSU2* __restrict__ piApure,
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

    deviceSU2 res = _deviceGetGaugeBCSU2DirZero(byFieldId, piApure, uiBigIdx, byMu); //Apure _mu
    deviceSU2 res2 = _deviceGetGaugeBCSU2DirZero(byFieldId, piA, uiBigIdx, byNu); //A _nu
    res2.Mul(res); //A _nu Apure _mu
    res.Mul(_deviceGetGaugeBCSU2DirZero(byFieldId, piA, uiBigIdx, byNu)); //Apure _mu A _nu
    res.Sub(res2); //[Apure, A]
    res.Add(_deviceGetGaugeBCSU2DirZeroSIndex(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_p_mu * _DC_Dir + byNu]).MulRealC(F(0.5)));
    res.Sub(_deviceGetGaugeBCSU2DirZeroSIndex(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_m_mu * _DC_Dir + byNu]).MulRealC(F(0.5)));
    return res;
}

static __device__ __inline__ deviceSU2 _devicePlaqutteSU2(
    const deviceSU2* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
    UINT uiSiteIndex,
    BYTE plaqIdx, //0-5, as 12, 13, 14, 23, 24, 34
    BYTE plaqLength, //Always 4
    BYTE plaqCountAll //Always 24
)
{
    SIndex first = pCachedPlaqutte[plaqIdx * plaqLength + uiSiteIndex * plaqCountAll];
    deviceSU2 toAdd(_deviceGetGaugeBCSU2DirOneSIndex(pDeviceData, first));
    if (first.NeedToDagger())
    {
        toAdd.Dagger();
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[plaqIdx * plaqLength + j + uiSiteIndex * plaqCountAll];
        deviceSU2 toMul(_deviceGetGaugeBCSU2DirOneSIndex(pDeviceData, first));
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
static __device__ __inline__ deviceSU2 _deviceLinkSU2(
    const deviceSU2* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceSU2 sRet = deviceSU2::makeSU2Id();
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

/**
 * After every move, it maps to inside the lattice
 * Do NOT use it in projective plane boundary condition
 */
static __device__ __inline__ deviceSU2 _deviceLinkLongSU2(
    const deviceSU2* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceSU2 sRet = deviceSU2::makeSU2Id();
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
static __device__ __inline__ deviceSU2 _device1PlaqutteTermPPSU2(
    const deviceSU2* __restrict__ pDeviceData,
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

    deviceSU2 u(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_mu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_p_mu_nu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_p_nu_mu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_nu, byFieldId));

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
static __device__ __inline__ deviceSU2 _device1PlaqutteTermMPSU2(
    const deviceSU2* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1);
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, byNu + 1);
    const UINT uin_m_mub4 = __idx->_deviceGetBigIndex(n_m_mu) * _DC_Dir;
    const SIndex& s_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byMu];
    const SIndex& s_m_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byNu];
    const SIndex& s_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_mu_p_nu) * _DC_Dir + byMu];
    const SIndex& s__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byNu];

    deviceSU2 u(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_m_mu__mu, byFieldId));
    u.DaggerMul(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_m_mu__nu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_m_mu_p_nu__mu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s__nu, byFieldId));

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
static __device__ __inline__ deviceSU2 _device1PlaqutteTermPMSU2(
    const deviceSU2* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byNu) - 1);
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(n_m_nu, byMu + 1);
    const UINT uin_m_nub4 = __idx->_deviceGetBigIndex(n_m_nu) * _DC_Dir;
    const SIndex& s_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byMu];
    const SIndex& s_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byNu];
    const SIndex& s_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_nu_p_mu) * _DC_Dir + byNu];
    const SIndex& s__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byMu];

    deviceSU2 u(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s__mu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_m_nu_p_mu__nu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_m_nu__mu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_m_nu__nu, byFieldId));

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
static __device__ __inline__ deviceSU2 _device1PlaqutteTermMMSU2(
    const deviceSU2* __restrict__ pDeviceData,
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
    deviceSU2 u(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_m_nu_m_mu__nu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_m_mu__mu, byFieldId));
    u.DaggerMul(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_m_nu_m_mu__mu, byFieldId));
    u.Mul(_deviceGetGaugeBCSU2DirSIndex(pDeviceData, s_m_nu__nu, byFieldId));

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
static __device__ __inline__ deviceSU2 _deviceCloverSU2(const deviceSU2* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    deviceSU2 ret(_device1PlaqutteTermPPSU2(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId));
    ret.Add(_device1PlaqutteTermMMSU2(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId));
    ret.Add(_device1PlaqutteTermPMSU2(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));
    ret.Add(_device1PlaqutteTermMPSU2(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));

    return ret;
}

/**
 * Avoid the add of matrices
 */
static __device__ __inline__ Real _deviceCloverRetrSU2(const deviceSU2* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    return _device1PlaqutteTermPPSU2(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId).ReTr()
        + _device1PlaqutteTermMMSU2(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId).ReTr()
        + _device1PlaqutteTermPMSU2(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId).ReTr()
        + _device1PlaqutteTermMPSU2(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId).ReTr();
}

#pragma endregion

#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_SU2_H_

//=============================================================================
// END OF FILE
//=============================================================================