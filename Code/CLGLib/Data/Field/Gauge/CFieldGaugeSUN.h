//=============================================================================
// FILENAME : CFieldGaugeSUN.h
// 
// DESCRIPTION:
// This is the class for the gauge fields
//
// REVISION:
//  [07/01/2024 nbale]
//=============================================================================

#ifndef _CFIELDGAUGE_SUN_H_
#define _CFIELDGAUGE_SUN_H_

__BEGIN_NAMESPACE

template<INT N, INT NofE>
class __DLL_EXPORT CFieldGaugeSUN : public CFieldGauge
{
public:

    CFieldGaugeSUN();
    ~CFieldGaugeSUN();
    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType) override;
    void InitialWithByte(BYTE* byData) override;
    void InitialWithByteCompressed(BYTE*) override { appCrucial(_T("Not supported by %s\n"), __ENUM_TO_STRING(EFieldType, GetFieldType()).c_str()); }
    void InitialField(EFieldInitialType eInitialType) override;
    EFieldType GetFieldType() const override { return EFT_GaugeSUN; }
    UINT MatrixN() const override { return N; }
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
    BYTE* CopyDataOut(UINT &uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;

    void PolyakovOnSpatialSite(cuDoubleComplex* buffer) const override;

    deviceSUN<N, NofE>* m_pDeviceData;

    _GetData
};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU4)

class CLGAPI CFieldGaugeSU4 : public CFieldGaugeSUN<4, 16>
{
    __CLGDECLARE_FIELD(CFieldGaugeSU4)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU4; }
};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU5)

class CLGAPI CFieldGaugeSU5 : public CFieldGaugeSUN<5, 32>
{
    __CLGDECLARE_FIELD(CFieldGaugeSU5)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU5; }
};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU6)

class CLGAPI CFieldGaugeSU6 : public CFieldGaugeSUN<6, 64>
{
    __CLGDECLARE_FIELD(CFieldGaugeSU6)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU6; }
};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU7)

class CLGAPI CFieldGaugeSU7 : public CFieldGaugeSUN<7, 64>
{
    __CLGDECLARE_FIELD(CFieldGaugeSU7)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU7; }
};

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU8)

class CLGAPI CFieldGaugeSU8 : public CFieldGaugeSUN<8, 64>
{
    __CLGDECLARE_FIELD(CFieldGaugeSU8)
public:
    EFieldType GetFieldType() const override { return EFT_GaugeSU8; }
};


#pragma region device functions

/**
* Note: for baked plaqutte index, the bond if is set to SIndex
* If it is a "new SIndex" instead, remember to set the m_byTag
*/
template<INT N, INT NofE>
static __device__ __inline__ const deviceSUN<N, NofE>& _deviceGetGaugeBCSUN(
    BYTE byFieldId,
    const deviceSUN<N, NofE>* __restrict__ pBuffer,
    const SIndex& idx)
{
    return idx.IsDirichlet() ?
        ((CFieldBoundaryGaugeSUN<N, NofE>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
}

/**
* If the bond is on surface, return the Dirichlet
* else, return the element
*/
template<INT N, INT NofE>
static __device__ __inline__ const deviceSUN<N, NofE>& _deviceGetGaugeBCSUNDir(
    BYTE byFieldId,
    const deviceSUN<N, NofE>* __restrict__ pBuffer,
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

template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _deviceGetGaugeBCSUNDirOne(
    BYTE byFieldId,
    const deviceSUN<N, NofE>* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        deviceSUN<N, NofE>::makeSUNId()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
}

template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _deviceGetGaugeBCSUNDirZero(
    BYTE byFieldId,
    const deviceSUN<N, NofE>* __restrict__ pBuffer,
    UINT uiBigIdx,
    BYTE byDir)
{
    return __idx->_deviceIsBondOnSurface(uiBigIdx, byDir) ?
        deviceSUN<N, NofE>::makeSUNZero()
        : pBuffer[_deviceGetLinkIndex(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].m_uiSiteIndex, byDir)];
}

template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _deviceGetGaugeBCSUNDirSIndex(
    const deviceSUN<N, NofE>* __restrict__ pBuffer,
    const SIndex& idx,
    BYTE byFieldId)
{
    deviceSUN<N, NofE> ret = idx.IsDirichlet() ?
        ((CFieldBoundaryGaugeSUN<N, NofE>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];

    if (idx.NeedToDagger())
    {
        ret.Dagger();
    }
    return ret;
}

template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _deviceGetGaugeBCSUNDirOneSIndex(
    const deviceSUN<N, NofE>* __restrict__ pBuffer,
    const SIndex& idx)
{
    if (idx.IsDirichlet())
    {
        return deviceSUN<N, NofE>::makeSUNId();
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
template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _deviceGetGaugeBCSUNDirZeroSIndex(
    const deviceSUN<N, NofE>* __restrict__ pBuffer,
    const SIndex& idx)
{
    if (idx.IsDirichlet())
    {
        return deviceSUN<N, NofE>::makeSUNZero();
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
template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _deviceDPureMuSUN(
    const deviceSUN<N, NofE>* __restrict__ piA,
    const deviceSUN<N, NofE>* __restrict__ piApure,
    const SSmallInt4& sSite4,
    UINT uiBigIdx,
    BYTE byMu,
    BYTE byNu,
    BYTE byFieldId)
{
    //i a D A = (A_nu (n) - A_nu (n-mu)) + iApure _mu A _nu - i A _nu Apure _mu
    const UINT uiSiteBig_m_mu = __idx->_deviceGetBigIndex(
        _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1));

    deviceSUN<N, NofE> res = _deviceGetGaugeBCSUNDirZero(piApure, uiBigIdx, byMu); //Apure _mu
    deviceSUN<N, NofE> res2 = _deviceGetGaugeBCSUNDirZero(piA, uiBigIdx, byNu); //A _nu
    res2.Mul(res); //A _nu Apure _mu
    res.Mul(_deviceGetGaugeBCSUNDirZero(piA, uiBigIdx, byNu)); //Apure _mu A _nu
    res.Sub(res2); //[Apure, A]
    res.Add(_deviceGetGaugeBCSUNDirZero(piA, uiBigIdx, byNu));
    res.Sub(_deviceGetGaugeBCSUNDirZeroSIndex(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_m_mu * _DC_Dir + byNu]));
    return res;
}

/**
 * test using (A(N+mu)-A(N-mu))/2
 */
template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _deviceDPureMu2SUN(
    const deviceSUN<N, NofE>* __restrict__ piA,
    const deviceSUN<N, NofE>* __restrict__ piApure,
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

    deviceSUN<N, NofE> res = _deviceGetGaugeBCSUNDirZero(piApure, uiBigIdx, byMu); //Apure _mu
    deviceSUN<N, NofE> res2 = _deviceGetGaugeBCSUNDirZero(piA, uiBigIdx, byNu); //A _nu
    res2.Mul(res); //A _nu Apure _mu
    res.Mul(_deviceGetGaugeBCSUNDirZero(piA, uiBigIdx, byNu)); //Apure _mu A _nu
    res.Sub(res2); //[Apure, A]
    res.Add(_deviceGetGaugeBCSUNDirZeroSIndex(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_p_mu * _DC_Dir + byNu]).MulRealC(F(0.5)));
    res.Sub(_deviceGetGaugeBCSUNDirZeroSIndex(piA,
        __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiSiteBig_m_mu * _DC_Dir + byNu]).MulRealC(F(0.5)));
    return res;
}

template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _devicePlaqutteSUN(
    const deviceSUN<N, NofE>* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
    UINT uiSiteIndex,
    BYTE plaqIdx, //0-5, as 12, 13, 14, 23, 24, 34
    BYTE plaqLength, //Always 4
    BYTE plaqCountAll //Always 24
)
{
    SIndex first = pCachedPlaqutte[plaqIdx * plaqLength + uiSiteIndex * plaqCountAll];
    deviceSUN<N, NofE> toAdd(_deviceGetGaugeBCSUNDirOneSIndex(pDeviceData, first));
    if (first.NeedToDagger())
    {
        toAdd.Dagger();
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[plaqIdx * plaqLength + j + uiSiteIndex * plaqCountAll];
        deviceSUN<N, NofE> toMul(_deviceGetGaugeBCSUNDirOneSIndex(pDeviceData, first));
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
template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _deviceLinkSUN(
    const deviceSUN<N, NofE>* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceSUN<N, NofE> sRet = deviceSUN<N, NofE>::makeSUNId();
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
template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _deviceLinkLongSUN(
    const deviceSUN<N, NofE>* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    deviceSUN<N, NofE> sRet = deviceSUN<N, NofE>::makeSUNId();
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
template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _device1PlaqutteTermPPSUN(
    const deviceSUN<N, NofE>* __restrict__ pDeviceData,
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

    deviceSUN<N, NofE> u(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_mu, byFieldId));
    u.Mul(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_p_mu_nu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_p_nu_mu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_nu, byFieldId));

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
template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _device1PlaqutteTermMPSUN(
    const deviceSUN<N, NofE>* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1);
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, byNu + 1);
    const UINT uin_m_mub4 = __idx->_deviceGetBigIndex(n_m_mu) * _DC_Dir;
    const SIndex& s_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byMu];
    const SIndex& s_m_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byNu];
    const SIndex& s_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_mu_p_nu) * _DC_Dir + byMu];
    const SIndex& s__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byNu];

    deviceSUN<N, NofE> u(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_m_mu__mu, byFieldId));
    u.DaggerMul(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_m_mu__nu, byFieldId));
    u.Mul(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_m_mu_p_nu__mu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s__nu, byFieldId));

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
template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _device1PlaqutteTermPMSUN(
    const deviceSUN<N, NofE>* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byNu) - 1);
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(n_m_nu, byMu + 1);
    const UINT uin_m_nub4 = __idx->_deviceGetBigIndex(n_m_nu) * _DC_Dir;
    const SIndex& s_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byMu];
    const SIndex& s_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byNu];
    const SIndex& s_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_nu_p_mu) * _DC_Dir + byNu];
    const SIndex& s__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byMu];

    deviceSUN<N, NofE> u(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s__mu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_m_nu_p_mu__nu, byFieldId));
    u.MulDagger(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_m_nu__mu, byFieldId));
    u.Mul(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_m_nu__nu, byFieldId));

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
template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _device1PlaqutteTermMMSUN(
    const deviceSUN<N, NofE>* __restrict__ pDeviceData,
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
    deviceSUN<N, NofE> u(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_m_nu_m_mu__nu, byFieldId));
    u.Mul(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_m_mu__mu, byFieldId));
    u.DaggerMul(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_m_nu_m_mu__mu, byFieldId));
    u.Mul(_deviceGetGaugeBCSUNDirSIndex(pDeviceData, s_m_nu__nu, byFieldId));

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
template<INT N, INT NofE>
static __device__ __inline__ deviceSUN<N, NofE> _deviceCloverSUN(const deviceSUN<N, NofE>* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    deviceSUN<N, NofE> ret(_device1PlaqutteTermPPSUN(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId));
    ret.Add(_device1PlaqutteTermMMSUN(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId));
    ret.Add(_device1PlaqutteTermPMSUN(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));
    ret.Add(_device1PlaqutteTermMPSUN(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId));

    return ret;
}

/**
 * Avoid the add of matrices
 */
template<INT N, INT NofE>
static __device__ __inline__ Real _deviceCloverRetrSUN(const deviceSUN<N, NofE>* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    return _device1PlaqutteTermPPSUN(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId).ReTr()
        + _device1PlaqutteTermMMSUN(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId).ReTr()
        + _device1PlaqutteTermPMSUN(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId).ReTr()
        + _device1PlaqutteTermMPSUN(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId).ReTr();
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_SUN_H_

//=============================================================================
// END OF FILE
//=============================================================================