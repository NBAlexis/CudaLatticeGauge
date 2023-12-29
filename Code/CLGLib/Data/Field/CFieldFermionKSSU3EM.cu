//=============================================================================
// FILENAME : CFieldFermionKSSU3EM.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [10/06/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionKSSU3EM)

#pragma region DOperator


#pragma region kernel

/**
 * very very strange, the eta_mu is not modified for projective plane
 * What do you mean eta_mu is not modified???
 *
 * Not support shift center
 *
 * Explain:
 * qBz: u_y = exp(i qBz x), u_x(L_x) = exp(-i qBz Lx y)
 *
 * qEz: u_t = exp(- i qEz z), u_z(L_z) = exp(i qEz Lz t)
 * 
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSEM_Simple(
    const deviceSU3Vector* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceSU3Vector* pResultData,
    Real fam,
    Real fqEz,
    Real fqBz,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    const Real fX1 = static_cast<Real>(sSite4.x - _DC_Centerx);
    const Real fY1 = static_cast<Real>(sSite4.y - _DC_Centery);
    const Real fZ1 = static_cast<Real>(sSite4.z - _DC_Centerz);
    const Real fT1 = static_cast<Real>(sSite4.w - _DC_Centert);

    Real u1pX = F(0.0);
    const Real u1pY = fX1 * fqBz;
    Real u1pZ = F(0.0);
    const Real u1pT = -fZ1 * fqEz;
    if (sSite4.x == _DC_Lx - 1)
    {
        u1pX = -fY1 * _DC_Lx * fqBz;
    }
    if (sSite4.z == _DC_Lz - 1)
    {
        u1pZ = fT1 * _DC_Lz * fqEz;
    }

    //idir = mu
    for (UINT idir = 0; idir < _DC_Dir; ++idir)
    {
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];
        const SSmallInt4 x_m_site = __deviceSiteIndexToInt4(x_m_mu_Fermion.m_uiSiteIndex);

        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        const Real eta_mu2 = (1 == ((pEtaTable[x_m_mu_Fermion.m_uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        const Real fX2 = static_cast<Real>(x_m_site.x - _DC_Centerx);
        const Real fY2 = static_cast<Real>(x_m_site.y - _DC_Centery);
        const Real fZ2 = static_cast<Real>(x_m_site.z - _DC_Centerz);
        const Real fT2 = static_cast<Real>(x_m_site.w - _DC_Centert);

        Real u1mX = F(0.0);
        const Real u1mY = fX2 * fqBz;
        Real u1mZ = F(0.0);
        const Real u1mT = - fZ2 * fqEz;
        if (x_m_site.x == _DC_Lx - 1)
        {
            u1mX = -fY2 * _DC_Lx * fqBz;
        }
        if (x_m_site.z == _DC_Lz - 1)
        {
            u1mZ = fT2 * _DC_Lz * fqEz;
        }

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        deviceSU3 x_Gauge_element = pGauge[linkIndex];
        if (0 == idir)
        {
            x_Gauge_element.MulComp(_make_cuComplex(_cos(u1pX), _sin(u1pX)));
        }
        else if (1 == idir)
        {
            x_Gauge_element.MulComp(_make_cuComplex(_cos(u1pY), _sin(u1pY)));
        }
        else if (2 == idir)
        {
            x_Gauge_element.MulComp(_make_cuComplex(_cos(u1pZ), _sin(u1pZ)));
        }
        else if (3 == idir)
        {
            x_Gauge_element.MulComp(_make_cuComplex(_cos(u1pT), _sin(u1pT)));
        }
        deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (0 == idir)
        {
            x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(u1mX), _sin(u1mX)));
        }
        else if (1 == idir)
        {
            x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(u1mY), _sin(u1mY)));
        }
        else if (2 == idir)
        {
            x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(u1mZ), _sin(u1mZ)));
        }
        else if (3 == idir)
        {
            x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(u1mT), _sin(u1mT)));
        }
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

        //U^{dagger}(x-mu) phi(x-mu)
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
        //u_dagger_phi_x_m_m.MulReal(bShiftCenter ? eta_mu2 : eta_mu);
        u_dagger_phi_x_m_m.MulReal(eta_mu2);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.Add(u_dagger_phi_x_m_m);
        }
        else
        {
            u_phi_x_p_m.Sub(u_dagger_phi_x_m_m);
        }
        //u_phi_x_p_m.MulReal(eta_mu);
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

#pragma endregion

#pragma region Derivate

/**
 * Explain:
 * qBz: u_y = exp(i qBz x), u_x(L_x) = exp(-i qBz Lx y)
 *
 * qEz: u_t = exp(- i qEz z), u_z(L_z) = exp(i qEz Lz t)
 *
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForceEM_Simple(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    Real fqEz,
    Real fqBz,
    BYTE byFieldId)
{
    intokernalInt4;

    const Real fX1 = static_cast<Real>(sSite4.x - _DC_Centerx);
    const Real fY1 = static_cast<Real>(sSite4.y - _DC_Centery);
    const Real fZ1 = static_cast<Real>(sSite4.z - _DC_Centerz);
    const Real fT1 = static_cast<Real>(sSite4.w - _DC_Centert);

    Real u1pX = F(0.0);
    const Real u1pY = fX1 * fqBz;
    Real u1pZ = F(0.0);
    const Real u1pT = -fZ1 * fqEz;
    if (sSite4.x == _DC_Lx - 1)
    {
        u1pX = -fY1 * _DC_Lx * fqBz;
    }
    if (sSite4.z == _DC_Lz - 1)
    {
        u1pZ = fT1 * _DC_Lz * fqEz;
    }

    //idir = mu
    for (UINT idir = 0; idir < _DC_Dir; ++idir)
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
            if (0 == idir)
            {
                toContract.MulComp(_make_cuComplex(_cos(u1pX), _sin(u1pX)));
            }
            else if (1 == idir)
            {
                toContract.MulComp(_make_cuComplex(_cos(u1pY), _sin(u1pY)));
            }
            else if (2 == idir)
            {
                toContract.MulComp(_make_cuComplex(_cos(u1pZ), _sin(u1pZ)));
            }
            else if (3 == idir)
            {
                toContract.MulComp(_make_cuComplex(_cos(u1pT), _sin(u1pT)));
            }
            deviceSU3 thisTerm = deviceSU3::makeSU3ContractV(phi_id[uiSiteIndex], toContract);

            toContract = pGauge[linkIndex].MulVector(phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
            if (0 == idir)
            {
                toContract.MulComp(_make_cuComplex(_cos(u1pX), _sin(u1pX)));
            }
            else if (1 == idir)
            {
                toContract.MulComp(_make_cuComplex(_cos(u1pY), _sin(u1pY)));
            }
            else if (2 == idir)
            {
                toContract.MulComp(_make_cuComplex(_cos(u1pZ), _sin(u1pZ)));
            }
            else if (3 == idir)
            {
                toContract.MulComp(_make_cuComplex(_cos(u1pT), _sin(u1pT)));
            }
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

            pForce[linkIndex].Sub(thisTerm);
        }
    }

}

__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGammaEM(
    deviceSU3Vector* pMe,
    const deviceSU3Vector* __restrict__ pOther,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    Real fqEz, 
    Real fqBz,
    BYTE byDir)
{
    intokernalInt4;
    const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> byDir) & 1)) ? F(-1.0) : F(1.0);
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byDir);
    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, byDir)];

    const Real fX1 = static_cast<Real>(sSite4.x - _DC_Centerx);
    const Real fY1 = static_cast<Real>(sSite4.y - _DC_Centery);
    const Real fZ1 = static_cast<Real>(sSite4.z - _DC_Centerz);
    const Real fT1 = static_cast<Real>(sSite4.w - _DC_Centert);

    Real u1pX = F(0.0);
    const Real u1pY = fX1 * fqBz;
    Real u1pZ = F(0.0);
    const Real u1pT = -fZ1 * fqEz;
    if (sSite4.x == _DC_Lx - 1)
    {
        u1pX = -fY1 * _DC_Lx * fqBz;
    }
    if (sSite4.z == _DC_Lz - 1)
    {
        u1pZ = fT1 * _DC_Lz * fqEz;
    }

    const SSmallInt4 x_m_site = __deviceSiteIndexToInt4(x_m_mu_Fermion.m_uiSiteIndex);
    const Real fX2 = static_cast<Real>(x_m_site.x - _DC_Centerx);
    const Real fY2 = static_cast<Real>(x_m_site.y - _DC_Centery);
    const Real fZ2 = static_cast<Real>(x_m_site.z - _DC_Centerz);
    const Real fT2 = static_cast<Real>(x_m_site.w - _DC_Centert);

    Real u1mX = F(0.0);
    const Real u1mY = fX2 * fqBz;
    Real u1mZ = F(0.0);
    const Real u1mT = -fZ2 * fqEz;
    if (x_m_site.x == _DC_Lx - 1)
    {
        u1mX = -fY2 * _DC_Lx * fqBz;
    }
    if (x_m_site.z == _DC_Lz - 1)
    {
        u1mZ = fT2 * _DC_Lz * fqEz;
    }

    if (0 == byDir)
    {
        x_Gauge_element.MulComp(_make_cuComplex(_cos(u1pX), _sin(u1pX)));
    }
    else if (1 == byDir)
    {
        x_Gauge_element.MulComp(_make_cuComplex(_cos(u1pY), _sin(u1pY)));
    }
    else if (2 == byDir)
    {
        x_Gauge_element.MulComp(_make_cuComplex(_cos(u1pZ), _sin(u1pZ)));
    }
    else if (3 == byDir)
    {
        x_Gauge_element.MulComp(_make_cuComplex(_cos(u1pT), _sin(u1pT)));
    }

    if (0 == byDir)
    {
        x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(u1mX), _sin(u1mX)));
    }
    else if (1 == byDir)
    {
        x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(u1mY), _sin(u1mY)));
    }
    else if (2 == byDir)
    {
        x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(u1mZ), _sin(u1mZ)));
    }
    else if (3 == byDir)
    {
        x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(u1mT), _sin(u1mT)));
    }

    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }
    pMe[uiSiteIndex] = x_Gauge_element.MulVector(pOther[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        pMe[uiSiteIndex].MulReal(F(-1.0));
    }
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        pMe[uiSiteIndex].Sub(x_m_mu_Gauge_element.MulVector(pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        pMe[uiSiteIndex].Add(x_m_mu_Gauge_element.MulVector(pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    pMe[uiSiteIndex].MulReal(F(0.5) * eta_mu);
}

#pragma endregion


#pragma endregion

#pragma region D and derivate

void CFieldFermionKSSU3EM::DOperatorKS(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceSU3Vector* pTarget = (deviceSU3Vector*)pTargetBuffer;
    const deviceSU3Vector* pSource = (const deviceSU3Vector*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    preparethread;
    _kernelDFermionKSEM_Simple << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        f2am,
        CCommonData::m_fEz * m_fQ,
        CCommonData::m_fBz * m_fQ,
        m_byFieldId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

void CFieldFermionKSSU3EM::DerivateD0(
    void* pForce,
    const void* pGaugeBuffer) const
{
    preparethread;
    _kernelDFermionKSForceEM_Simple << <block, threads >> > (
        (const deviceSU3*)pGaugeBuffer,
        (deviceSU3*)pForce,
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        m_pRationalFieldPointers,
        m_pMDNumerator,
        m_rMD.m_uiDegree,
        CCommonData::m_fEz * m_fQ,
        CCommonData::m_fBz * m_fQ,
        m_byFieldId);
}

#pragma endregion

void CFieldFermionKSSU3EM::InitialOtherParameters(CParameters& params)
{
    CFieldFermionKSSU3::InitialOtherParameters(params);
    Real fa2Ez = F(0.0);
    Real fa2Bz = F(0.0);
    if (params.FetchValueReal(_T("Qa2Ez"), fa2Ez))
    {
        CCommonData::m_fEz = fa2Ez;
    }
    if (params.FetchValueReal(_T("Qa2Bz"), fa2Bz))
    {
        CCommonData::m_fBz = fa2Bz;
    }
    params.FetchValueReal(_T("EMChange"), m_fQ);
}

void CFieldFermionKSSU3EM::ApplyGammaKS(const CFieldGauge* pGauge, EGammaMatrix eGamma)
{
    INT iDir = -1;
    switch (eGamma)
    {
    case UNITY:
    {
        return;
    }
    case GAMMA1:
    {
        iDir = 0;
    }
    break;
    case GAMMA2:
    {
        iDir = 1;
    }
    break;
    case GAMMA3:
    {
        iDir = 2;
    }
    break;
    case GAMMA4:
    {
        iDir = 3;
    }
    break;
    }

    if (iDir >= 0)
    {
        if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKSSU3 can only play with gauge SU3!"));
            return;
        }
        const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
        CFieldFermionKSSU3EM* pPooled = dynamic_cast<CFieldFermionKSSU3EM*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
        checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
        preparethread;
        _kernelKSApplyGammaEM << <block, threads >> > (
            m_pDeviceData,
            pPooled->m_pDeviceData,
            pFieldSU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            CCommonData::m_fEz * m_fQ,
            CCommonData::m_fBz * m_fQ,
            static_cast<BYTE>(iDir));

        pPooled->Return();
    }
    appParanoiac(_T("Calling CFieldFermionKSSU3EM apply gamma %d\n"), iDir);
}

void CFieldFermionKSSU3EM::CopyTo(CField* U) const
{
    CFieldFermionKSSU3::CopyTo(U);
    CFieldFermionKSSU3EM* target = dynamic_cast<CFieldFermionKSSU3EM*>(U);
    if (NULL != target)
    {
        target->m_fQ = m_fQ;
    }
}

CCString CFieldFermionKSSU3EM::GetInfos(const CCString& tab) const
{
    CCString sRet = tab + _T("Name : CFieldFermionKSSU3R\n");
    sRet = sRet + tab + _T("Mass (2am) : ") + appFloatToString(m_f2am) + _T("\n");
    sRet = sRet + tab + _T("MD Rational (c) : ") + appFloatToString(m_rMD.m_fC) + _T("\n");
    sRet = sRet + tab + _T("MC Rational (c) : ") + appFloatToString(m_rMC.m_fC) + _T("\n");
    sRet = sRet + tab + _T("Q x a^2Ez : ") + appFloatToString(m_fQ * CCommonData::m_fEz) + _T("\n");
    sRet = sRet + tab + _T("Q x a^2Bz : ") + appFloatToString(m_fQ * CCommonData::m_fBz) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================