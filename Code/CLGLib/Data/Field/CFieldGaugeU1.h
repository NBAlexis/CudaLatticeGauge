//=============================================================================
// FILENAME : CFieldGaugeU1.h
// 
// DESCRIPTION:
// This is the class for the gauge fields
//
// REVISION:
//  [10/13/2020 nbale]
//=============================================================================

#ifndef _CFIELDGAUGE_U1_H_
#define _CFIELDGAUGE_U1_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CFieldGaugeU1)

class CLGAPI CFieldGaugeU1 : public CFieldGauge
{
    __CLGDECLARE_FIELD(CFieldGaugeU1)

public:
    CFieldGaugeU1();
    ~CFieldGaugeU1();

    void InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFileType) override;
    void InitialWithByte(BYTE* byData) override;
    void InitialWithByteCompressed(BYTE*) override;
    void InitialField(EFieldInitialType eInitialType) override;
    EFieldType GetFieldType() const override { return EFT_GaugeU1; }
    void DebugPrintMe() const override;

#pragma region HMC

    void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const override;
    void CalculateOnlyStaple(CFieldGauge* pStaple) const override;
    void MakeRandomGenerator() override;
#if !_CLG_DOUBLEFLOAT
    DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const override;
    DOUBLE CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const override;
    DOUBLE CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge* pStaple) const override;
    DOUBLE CalculateKinematicEnergy() const override;
#else
    Real CalculatePlaqutteEnergy(Real betaOverN) const override;
    Real CalculatePlaqutteEnergyUseClover(Real betaOverN) const override;
    Real CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge *pStaple) const override;
    Real CalculateKinematicEnergy() const override;
#endif

#pragma endregion

#pragma region BLAS

    void Zero() override;
    void Identity() override;
    void Dagger() override;

    void AxpyPlus(const CField* x) override;
    void AxpyMinus(const CField* x) override;
    void Axpy(Real a, const CField* x) override;
    void Axpy(const CLGComplex& a, const CField* x) override;
    void ScalarMultply(const CLGComplex& a) override;
    void ScalarMultply(Real a) override;

    void SetOneDirectionUnity(BYTE byDir) override;
    void SetOneDirectionZero(BYTE byDir) override;

#pragma endregion

#pragma region Test Functions to test gauge invarience of angular momentum

    /**
     * iA = U.TA() / 2
     */
    void TransformToIA() override;

    /**
     * U=exp(iA)
     */
    void TransformToU() override;

    void CalculateE_Using_U(CFieldGauge* pResoult) const override;

    void CalculateNablaE_Using_U(CFieldGauge* pResoult, UBOOL bNaive = FALSE) const override;

#pragma endregion

    void ExpMult(Real a, CField* U) const override;

    void ElementNormalize() override;
#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex Dot(const CField* other) const override;
#else
    CLGComplex Dot(const CField* other) const override;
#endif
    BYTE* CopyDataOut(UINT &uiSize) const override;
    BYTE* CopyDataOutFloat(UINT& uiSize) const override;
    BYTE* CopyDataOutDouble(UINT& uiSize) const override;
    CCString GetInfos(const CCString &tab) const override;

    CLGComplex* m_pDeviceData;

protected:

    void SetByArray(Real* array);
};

#pragma region Helper device functions

static __device__ __inline__ CLGComplex _deviceLinkU1(
    const CLGComplex* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    CLGComplex sRet = _onec;
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
                    sRet.y = -sRet.y;
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
                    const CLGComplex& toMul = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                    sRet = _make_cuComplex(sRet.x * toMul.x + sRet.y * toMul.y, sRet.y * toMul.x - sRet.x * toMul.y);
                }
                else
                {
                    const CLGComplex& toMul = pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)];
                    sRet = _make_cuComplex(sRet.x * toMul.x - sRet.y * toMul.y, sRet.y * toMul.x + sRet.x * toMul.y);
                }
            }
        }

        if (pDir[i] > 0) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return sRet;
}

static __device__ __inline__ Real _deviceLinkU1ArgSum(
    const CLGComplex* __restrict__ pDeviceData,
    SSmallInt4 sStartSite, BYTE byLength, BYTE byFieldId,
    const INT* __restrict__ pDir)
{
    //length can be 0
    Real sRet = F(0.0);
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
                sRet = __cuCargf(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                if ((newLink.NeedToDagger() && !bDagger)
                 || (!newLink.NeedToDagger() && bDagger)
                    )
                {
                    sRet = -sRet;
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
                    sRet -= __cuCargf(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
                else
                {
                    sRet += __cuCargf(pDeviceData[_deviceGetLinkIndex(newLink.m_uiSiteIndex, newLink.m_byDir)]);
                }
            }
        }

        if (pDir[i] > 0) //Move
        {
            _deviceSmallInt4Offset(sStartSite, pDir[i]);
        }
    }

    return sRet;
}


static __device__ __inline__ CLGComplex _deviceGetGaugeBCU1DirSIndex(
    const CLGComplex* __restrict__ pBuffer,
    const SIndex& idx,
    BYTE byFieldId)
{
    CLGComplex ret = idx.IsDirichlet() ?
        ((CFieldBoundaryGaugeU1*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[
            __idx->_devcieExchangeBoundaryFieldSiteIndex(idx) * _DC_Dir + idx.m_byDir
        ]
        : pBuffer[_deviceGetLinkIndex(idx.m_uiSiteIndex, idx.m_byDir)];
    if (idx.NeedToDagger())
    {
        ret = _cuConjf(ret);
    }
    return ret;
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
static __device__ __inline__ CLGComplex _device1PlaqutteTermPPU1(
    const CLGComplex* __restrict__ pDeviceData,
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

    CLGComplex u = _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_mu, byFieldId);
    u = _cuCmulf(u, _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_p_mu_nu, byFieldId));
    u = _cuCmulf(u, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, s_p_nu_mu, byFieldId)));
    u = _cuCmulf(u, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, s_nu, byFieldId)));

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
static __device__ __inline__ CLGComplex _device1PlaqutteTermMPU1(
    const CLGComplex* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_mu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byMu) - 1);
    const SSmallInt4 n_m_mu_p_nu = _deviceSmallInt4OffsetC(n_m_mu, byNu + 1);
    const UINT uin_m_mub4 = __idx->_deviceGetBigIndex(n_m_mu) * _DC_Dir;
    const SIndex& s_m_mu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byMu];
    const SIndex& s_m_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_mub4 + byNu];
    const SIndex& s_m_mu_p_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_mu_p_nu) * _DC_Dir + byMu];
    const SIndex& s__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byNu];

    CLGComplex u = _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_mu__mu, byFieldId);
    u = _cuCmulf(_cuConjf(u), _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_mu__nu, byFieldId));
    u = _cuCmulf(u, _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_mu_p_nu__mu, byFieldId));
    u = _cuCmulf(u, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, s__nu, byFieldId)));

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
static __device__ __inline__ CLGComplex _device1PlaqutteTermPMU1(
    const CLGComplex* __restrict__ pDeviceData,
    BYTE byMu, BYTE byNu, UINT uiBigIdx, const SSmallInt4& sSite4, BYTE byFieldId)
{
    const SSmallInt4 n_m_nu = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(byNu) - 1);
    const SSmallInt4 n_m_nu_p_mu = _deviceSmallInt4OffsetC(n_m_nu, byMu + 1);
    const UINT uin_m_nub4 = __idx->_deviceGetBigIndex(n_m_nu) * _DC_Dir;
    const SIndex& s_m_nu__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byMu];
    const SIndex& s_m_nu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uin_m_nub4 + byNu];
    const SIndex& s_m_nu_p_mu__nu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(n_m_nu_p_mu) * _DC_Dir + byNu];
    const SIndex& s__mu = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * _DC_Dir + byMu];

    CLGComplex u = _deviceGetGaugeBCU1DirSIndex(pDeviceData, s__mu, byFieldId);
    u = _cuCmulf(u, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_nu_p_mu__nu, byFieldId)));
    u = _cuCmulf(u, _cuConjf(_deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_nu__mu, byFieldId)));
    u = _cuCmulf(u, _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_nu__nu, byFieldId));

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
static __device__ __inline__ CLGComplex _device1PlaqutteTermMMU1(
    const CLGComplex* __restrict__ pDeviceData,
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
    CLGComplex u = _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_nu_m_mu__nu, byFieldId);
    u = _cuCmulf(u, _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_mu__mu, byFieldId));
    u = _cuCmulf(_cuConjf(u), _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_nu_m_mu__mu, byFieldId));
    u = _cuCmulf(u, _deviceGetGaugeBCU1DirSIndex(pDeviceData, s_m_nu__nu, byFieldId));

    return u;
}


static __device__ __inline__ Real _deviceCloverRetrU1(const CLGComplex* __restrict__ pGaugeField, const SSmallInt4& sSite4, UINT uiBigIdx, BYTE mu, BYTE nu, BYTE byFieldId)
{
    return _device1PlaqutteTermPPU1(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId).x
        + _device1PlaqutteTermMMU1(pGaugeField, mu, nu, uiBigIdx, sSite4, byFieldId).x
        + _device1PlaqutteTermPMU1(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId).x
        + _device1PlaqutteTermMPU1(pGaugeField, nu, mu, uiBigIdx, sSite4, byFieldId).x;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_U1_H_

//=============================================================================
// END OF FILE
//=============================================================================