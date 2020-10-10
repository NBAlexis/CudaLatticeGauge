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

__CLGIMPLEMENT_CLASS(CFieldFermionKSSU3P4)

#pragma region DOperator

#pragma region kernel

#pragma region C10 term

/**
* Note, it is am NOT 2am in this function 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_P4_Fat(
    const deviceSU3Vector * __restrict__ pDeviceData,
    const deviceSU3 * __restrict__ pGauge,
    const BYTE * __restrict__ pEtaTable,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fUFactor,
    Real fFatFactor,
    Real fam,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    //idir = mu
    #pragma unroll
    for (BYTE idir = 0; idir < 4; ++idir)
    {
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //This is in fact, -1 * eta(n + mu)
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        const Real eta_mu2 = (1 == ((pEtaTable[x_m_mu_Gauge.m_uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        deviceSU3 x_Gauge_element = pGauge[linkIndex];
        x_Gauge_element.MulReal(fUFactor);
        deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        x_m_mu_Gauge_element.MulReal(fUFactor);

        deviceSU3 toAdd = deviceSU3::makeSU3Zero();
        #pragma unroll
        for (BYTE byNu = 0; byNu < 4; ++byNu)
        {
            if (byNu != idir)
            {
                INT path[3] = { __fwd(byNu), __fwd(idir), __bck(byNu) };
                toAdd.Add(_deviceLink(pGauge, sSite4, 3, byGaugeFieldId, path));
                path[0] = path[2]; path[2] = __fwd(byNu);
                toAdd.Add(_deviceLink(pGauge, sSite4, 3, byGaugeFieldId, path));
            }
        }
        toAdd.MulReal(fFatFactor);
        x_Gauge_element.Add(toAdd);

        toAdd = deviceSU3::makeSU3Zero();
        //const SSmallInt4 sSite_m_mu = _deviceSmallInt4OffsetC(sSite4, __bck(idir));
        #pragma unroll
        for (BYTE byNu = 0; byNu < 4; ++byNu)
        {
            if (byNu != idir)
            {
                INT path[3] = { __fwd(byNu), __bck(idir), __bck(byNu) };
                toAdd.Add(_deviceLink(pGauge, sSite4, 3, byGaugeFieldId, path));
                path[0] = path[2]; path[2] = __fwd(byNu);
                toAdd.Add(_deviceLink(pGauge, sSite4, 3, byGaugeFieldId, path));
            }
        }
        toAdd.MulReal(fFatFactor);
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }
        x_m_mu_Gauge_element.Add(toAdd);

        //U(x,mu) phi(x+ mu)
        deviceSU3Vector u_phi_x_p_m = x_Gauge_element.MulVector(pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.MulReal(F(-1.0) * eta_mu);
        }
        else
        {
            u_phi_x_p_m.MulReal(eta_mu);
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceSU3Vector u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulVector(pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        u_dagger_phi_x_m_m.MulReal(eta_mu2);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.Add(u_dagger_phi_x_m_m);
        }
        else
        {
            u_phi_x_p_m.Sub(u_dagger_phi_x_m_m);
        }
        result.Add(u_phi_x_p_m);
    }

    pResultData[uiSiteIndex].MulReal(fam);
    if (bDDagger)
    {
        pResultData[uiSiteIndex].Sub(result);
    }
    else
    {
        pResultData[uiSiteIndex].Add(result);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        pResultData[uiSiteIndex].MulReal(fCoeff);
        break;
    case EOCT_Complex:
        pResultData[uiSiteIndex].MulComp(cCoeff);
        break;
    }
}


/**
 * This is a copy from CFieldFermionKSSU3.cu
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_P4_UTerm(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    Real fCoefficient,
    UINT uiRational,
    BYTE byFieldId)
{
    intokernaldir;

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];

        for (UINT uiR = 0; uiR < uiRational; ++uiR)
        {
            const deviceSU3Vector* phi_i = pFermionPointers[uiR];
            const deviceSU3Vector* phi_id = pFermionPointers[uiR + uiRational];

            deviceSU3Vector toContract = pGauge[linkIndex].MulVector(phi_i[x_p_mu_Fermion.m_uiSiteIndex]);
            deviceSU3 thisTerm = deviceSU3::makeSU3ContractV(phi_id[uiSiteIndex], toContract);

            toContract = pGauge[linkIndex].MulVector(phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
            thisTerm.Add(deviceSU3::makeSU3ContractV(toContract, phi_i[uiSiteIndex]));

            if (x_p_mu_Fermion.NeedToOpposite())
            {
                thisTerm.MulReal(eta_mu * pNumerators[uiR] * F(-1.0));
            }
            else
            {
                thisTerm.MulReal(eta_mu * pNumerators[uiR]);
            }

            thisTerm.Ta();
            thisTerm.MulReal(fCoefficient);
            pForce[linkIndex].Sub(thisTerm);
        }
    }
}

#pragma endregion

#pragma region C12 term

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_P4_C12(
    const deviceSU3Vector* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fHalfC12,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    INT path[3] = {0, 0, 0};

    //idir = mu
    #pragma unroll
    for (BYTE idir = 0; idir < 4; ++idir)
    {
        #pragma unroll
        for (BYTE byNu = 0; byNu < 4; ++byNu)
        {
            if (byNu != idir)
            {
                //==================================
                //U12
                path[0] = __fwd(idir); path[1] = __fwd(byNu); path[2] = __fwd(byNu);
                SSmallInt4 siten = _deviceSmallInt4OffsetC(sSite4, path, 3);
                SIndex sn = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten)];
                deviceSU3 vn = _deviceLink(pGauge, sSite4, 3, byGaugeFieldId, path);
                INT etamu = pEtaTable[uiSiteIndex] >> idir;
                if (sn.NeedToOpposite()) { etamu = etamu + 1; }
                if (etamu & 1)
                {
                    result.Sub(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }
                else
                {
                    result.Add(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }

                path[0] = __fwd(byNu); path[1] = __fwd(byNu); path[2] = __fwd(idir);
                vn = _deviceLink(pGauge, sSite4, 3, byGaugeFieldId, path);
                if (etamu & 1)
                {
                    result.Sub(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }
                else
                {
                    result.Add(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }

                //==================================
                //U21
                path[0] = __bck(byNu); path[1] = __bck(byNu); path[2] = __bck(idir);
                siten = _deviceSmallInt4OffsetC(sSite4, path, 3);
                sn = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten)];
                vn = _deviceLink(pGauge, sSite4, 3, byGaugeFieldId, path);
                etamu = pEtaTable[sn.m_uiSiteIndex] >> idir;
                if (sn.NeedToOpposite()) { etamu = etamu + 1; }
                if (etamu & 1)
                {
                    result.Add(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }
                else
                {
                    result.Sub(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }

                path[0] = __bck(idir); path[1] = __bck(byNu); path[2] = __bck(byNu);
                vn = _deviceLink(pGauge, sSite4, 3, byGaugeFieldId, path);
                if (etamu & 1)
                {
                    result.Add(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }
                else
                {
                    result.Sub(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }

                //==================================
                //U1,-2
                path[0] = __fwd(idir); path[1] = __bck(byNu); path[2] = __bck(byNu);
                siten = _deviceSmallInt4OffsetC(sSite4, path, 3);
                sn = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten)];
                vn = _deviceLink(pGauge, sSite4, 3, byGaugeFieldId, path);
                etamu = pEtaTable[uiSiteIndex] >> idir;
                if (sn.NeedToOpposite()) { etamu = etamu + 1; }
                if (etamu & 1)
                {
                    result.Sub(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }
                else
                {
                    result.Add(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }

                path[0] = __bck(byNu); path[1] = __bck(byNu); path[2] = __fwd(idir);
                vn = _deviceLink(pGauge, sSite4, 3, byGaugeFieldId, path);
                if (etamu & 1)
                {
                    result.Sub(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }
                else
                {
                    result.Add(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }

                //U-2,1
                path[0] = __fwd(byNu); path[1] = __fwd(byNu); path[2] = __bck(idir);
                siten = _deviceSmallInt4OffsetC(sSite4, path, 3);
                sn = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten)];
                vn = _deviceLink(pGauge, sSite4, 3, byGaugeFieldId, path);
                etamu = pEtaTable[sn.m_uiSiteIndex] >> idir;
                if (sn.NeedToOpposite()) { etamu = etamu + 1; }
                if (etamu & 1)
                {
                    result.Add(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }
                else
                {
                    result.Sub(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }

                path[0] = __bck(idir); path[1] = __fwd(byNu); path[2] = __fwd(byNu);
                vn = _deviceLink(pGauge, sSite4, 3, byGaugeFieldId, path);
                if (etamu & 1)
                {
                    result.Add(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }
                else
                {
                    result.Sub(vn.MulVector(pDeviceData[sn.m_uiSiteIndex]));
                }
            }
        }
    }

    if (bDDagger)
    {
        fHalfC12 = fHalfC12 * F(-1.0);
    }
    result.MulReal(fHalfC12);

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

#pragma endregion

#pragma region D and derivate

void CFieldFermionKSSU3P4::DOperatorKS(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceSU3Vector* pTarget = (deviceSU3Vector*)pTargetBuffer;
    const deviceSU3Vector* pSource = (const deviceSU3Vector*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    const Real fDenorm = F(1.0) + F(6.0) * m_fomega;
    const Real fUFactor = m_fc10 / fDenorm;
    const Real fFatFactor = m_fc10 * m_fomega / fDenorm;
    const Real fHalfC12 = m_fc12 * F(0.5);

    preparethread;
    _kernelDFermionKS_P4_Fat << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        pTarget,
        m_byFieldId,
        1,
        fUFactor,
        fFatFactor,
        f2am * F(0.5),
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    _kernelDFermionKS_P4_C12 << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        m_byFieldId,
        1,
        fHalfC12,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

void CFieldFermionKSSU3P4::DerivateD0(
    void* pForce,
    const void* pGaugeBuffer) const
{
    const Real fDenorm = F(1.0) + F(6.0) * m_fomega;
    const Real fUFactor = m_fc10 / fDenorm;
    const Real fFatFactor = m_fc10 * m_fomega / fDenorm;
    const Real fHalfC12 = m_fc12 * F(0.5);

    preparethread;
    _kernelDFermionKSForce_P4_UTerm << <block, threads >> > (
        (const deviceSU3*)pGaugeBuffer,
        (deviceSU3*)pForce,
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        m_pRationalFieldPointers,
        m_pMDNumerator,
        fUFactor,
        m_rMD.m_uiDegree,
        m_byFieldId);

    for (INT i = 0; i < _HC_Diri; ++i)
    {
        for (INT j = 0; j < _HC_Diri; ++j)
        {
            if (j != i)
            {
                INT path[3] = { j + 1, i + 1, -j - 1 };
                checkCudaErrors(cudaMemcpy(m_pDevicePathBuffer, path, sizeof(INT) * 3, cudaMemcpyHostToDevice));
                OneLinkForce(
                    (const deviceSU3*)pGaugeBuffer,
                    1,
                    (deviceSU3*)pForce,
                    fFatFactor,
                    m_pDevicePathBuffer,
                    3,
                    static_cast<BYTE>(i)
                );

                path[0] = - j - 1; path[2] = j + 1;
                checkCudaErrors(cudaMemcpy(m_pDevicePathBuffer, path, sizeof(INT) * 3, cudaMemcpyHostToDevice));
                OneLinkForce(
                    (const deviceSU3*)pGaugeBuffer,
                    1,
                    (deviceSU3*)pForce,
                    fFatFactor,
                    m_pDevicePathBuffer,
                    3,
                    static_cast<BYTE>(i)
                );

                //==========================================
                //c12 terms
                path[0] = i + 1; path[1] = j + 1; path[2] = j + 1;
                checkCudaErrors(cudaMemcpy(m_pDevicePathBuffer, path, sizeof(INT) * 3, cudaMemcpyHostToDevice));
                OneLinkForce(
                    (const deviceSU3*)pGaugeBuffer,
                    1,
                    (deviceSU3*)pForce,
                    fHalfC12,
                    m_pDevicePathBuffer,
                    3,
                    static_cast<BYTE>(i)
                );

                path[0] = j + 1; path[1] = j + 1; path[2] = i + 1;
                checkCudaErrors(cudaMemcpy(m_pDevicePathBuffer, path, sizeof(INT) * 3, cudaMemcpyHostToDevice));
                OneLinkForce(
                    (const deviceSU3*)pGaugeBuffer,
                    1,
                    (deviceSU3*)pForce,
                    fHalfC12,
                    m_pDevicePathBuffer,
                    3,
                    static_cast<BYTE>(i)
                );

                path[0] = i + 1; path[1] = -j - 1; path[2] = -j - 1;
                checkCudaErrors(cudaMemcpy(m_pDevicePathBuffer, path, sizeof(INT) * 3, cudaMemcpyHostToDevice));
                OneLinkForce(
                    (const deviceSU3*)pGaugeBuffer,
                    1,
                    (deviceSU3*)pForce,
                    fHalfC12,
                    m_pDevicePathBuffer,
                    3,
                    static_cast<BYTE>(i)
                );

                path[0] = -j - 1; path[1] = -j - 1; path[2] = i + 1;
                checkCudaErrors(cudaMemcpy(m_pDevicePathBuffer, path, sizeof(INT) * 3, cudaMemcpyHostToDevice));
                OneLinkForce(
                    (const deviceSU3*)pGaugeBuffer,
                    1,
                    (deviceSU3*)pForce,
                    fHalfC12,
                    m_pDevicePathBuffer,
                    3,
                    static_cast<BYTE>(i)
                );
            }
        }
    }
}

#pragma endregion

CFieldFermionKSSU3P4::CFieldFermionKSSU3P4()
    : CFieldFermionKSSU3()
    , m_fomega(F(0.25))
    , m_fc10(F(0.3915))
    , m_fc12(F(0.0180833333333333))
    , m_pDevicePathBuffer(NULL)
{
    checkCudaErrors(cudaMalloc((void**)&m_pDevicePathBuffer, sizeof(INT) * 3));
}

CFieldFermionKSSU3P4::~CFieldFermionKSSU3P4()
{
    checkCudaErrors(cudaFree(m_pDevicePathBuffer));
}

void CFieldFermionKSSU3P4::InitialOtherParameters(CParameters& params)
{
    CFieldFermionKSSU3::InitialOtherParameters(params);
    m_bEachSiteEta = TRUE;
}

void CFieldFermionKSSU3P4::CopyTo(CField* U) const
{
    CFieldFermionKSSU3::CopyTo(U);
}

CCString CFieldFermionKSSU3P4::GetInfos(const CCString& tab) const
{
    CCString sRet = tab + _T("Name : CFieldFermionKSSU3R\n");
    sRet = sRet + tab + _T("Mass (2am) : ") + appFloatToString(m_f2am) + _T("\n");
    sRet = sRet + tab + _T("MD Rational (c) : ") + appFloatToString(m_rMD.m_fC) + _T("\n");
    sRet = sRet + tab + _T("MC Rational (c) : ") + appFloatToString(m_rMC.m_fC) + _T("\n");
    sRet = sRet + tab + _T("omega : ") + appFloatToString(m_fomega) + _T("\n");
    sRet = sRet + tab + _T("c10 : ") + appFloatToString(m_fc10) + _T("\n");
    sRet = sRet + tab + _T("c12 : ") + appFloatToString(m_fc12) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================