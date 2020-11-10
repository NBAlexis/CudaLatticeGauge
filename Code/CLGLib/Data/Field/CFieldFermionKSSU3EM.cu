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
    UBOOL bShiftCenter,
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    Real fX = static_cast<Real>(sSite4.x - sCenter.x);
    Real fY = static_cast<Real>(sSite4.y - sCenter.y);
    Real fZ = static_cast<Real>(sSite4.z - sCenter.z);
    //SSmallInt4 site_shift = sSite4;
    //site_shift.x = site_shift.x - 1;
    //site_shift = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_shift)].m_uiSiteIndex);
    //Real fXm = static_cast<Real>(sSite4.x - sCenter.x);
    //site_shift = sSite4;
    //site_shift.y = site_shift.y - 1;
    //site_shift = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_shift)].m_uiSiteIndex);
    //Real fYm = static_cast<Real>(sSite4.y - sCenter.y);
    //site_shift = sSite4;
    //site_shift.z = site_shift.z - 1;
    //site_shift = __deviceSiteIndexToInt4(__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_shift)].m_uiSiteIndex);
    //Real fZm = static_cast<Real>(sSite4.z - sCenter.z);

    if (bShiftCenter)
    {
        fX += F(0.5);
        fY += F(0.5);
        //fXm += F(0.5);
        //fYm += F(0.5);
    }
    fX = -fX * fqBz;
    fY = fY * fqBz;
    fZ = fZ * fqEz;
    //fXm = -fXm * fqBz;
    //fYm = fYm * fqBz;
    //fZm = fZm * fqEz;

    //idir = mu
    for (UINT idir = 0; idir < _DC_Dir; ++idir)
    {
        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        deviceSU3 x_Gauge_element = pGauge[linkIndex];
        if (0 == idir)
        {
            x_Gauge_element.MulComp(_make_cuComplex(_cos(fY), _sin(fY)));
        }
        else if (1 == idir)
        {
            x_Gauge_element.MulComp(_make_cuComplex(_cos(fX), _sin(fX)));
        }
        else if (3 == idir)
        {
            x_Gauge_element.MulComp(_make_cuComplex(_cos(fZ), _sin(fZ)));
        }
        deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (0 == idir)
        {
            x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(fY), _sin(fY)));
        }
        else if (1 == idir)
        {
            x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(fX), _sin(fX)));
        }
        else if (3 == idir)
        {
            x_m_mu_Gauge_element.MulComp(_make_cuComplex(_cos(fZ), _sin(fZ)));
        }
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

        //U(x,mu) phi(x+ mu)
        deviceSU3Vector u_phi_x_p_m = x_Gauge_element.MulVector(pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.MulReal(F(-1.0));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceSU3Vector u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulVector(pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.Add(u_dagger_phi_x_m_m);
        }
        else
        {
            u_phi_x_p_m.Sub(u_dagger_phi_x_m_m);
        }
        u_phi_x_p_m.MulReal(eta_mu);
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
    UBOOL bShiftCenter,
    SSmallInt4 sCenter,
    BYTE byFieldId)
{
    intokernalInt4;

    Real fX = static_cast<Real>(sSite4.x - sCenter.x);
    Real fY = static_cast<Real>(sSite4.y - sCenter.y);
    Real fZ = static_cast<Real>(sSite4.z - sCenter.z);
    if (bShiftCenter)
    {
        fX += F(0.5);
        fY += F(0.5);
    }
    fX = -fX * fqBz;
    fY = fY * fqBz;
    fZ = fZ * fqEz;

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
                toContract.MulComp(_make_cuComplex(_cos(fY), _sin(fY)));
            }
            else if (1 == idir)
            {
                toContract.MulComp(_make_cuComplex(_cos(fX), _sin(fX)));
            }
            else if (3 == idir)
            {
                toContract.MulComp(_make_cuComplex(_cos(fZ), _sin(fZ)));
            }
            deviceSU3 thisTerm = deviceSU3::makeSU3ContractV(phi_id[uiSiteIndex], toContract);

            toContract = pGauge[linkIndex].MulVector(phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
            if (0 == idir)
            {
                toContract.MulComp(_make_cuComplex(_cos(fY), _sin(fY)));
            }
            else if (1 == idir)
            {
                toContract.MulComp(_make_cuComplex(_cos(fX), _sin(fX)));
            }
            else if (3 == idir)
            {
                toContract.MulComp(_make_cuComplex(_cos(fZ), _sin(fZ)));
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
        m_fa2Ez,
        m_fa2Bz,
        m_bEachSiteEta,
        CCommonData::m_sCenter,
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
        m_fa2Ez,
        m_fa2Bz,
        m_bEachSiteEta,
        CCommonData::m_sCenter,
        m_byFieldId);
}

#pragma endregion

void CFieldFermionKSSU3EM::InitialOtherParameters(CParameters& params)
{
    CFieldFermionKSSU3::InitialOtherParameters(params);
    params.FetchValueReal(_T("Qa2Ez"), m_fa2Ez);
    params.FetchValueReal(_T("Qa2Bz"), m_fa2Bz);
    m_bEachSiteEta = TRUE;
}

void CFieldFermionKSSU3EM::CopyTo(CField* U) const
{
    CFieldFermionKSSU3::CopyTo(U);
    CFieldFermionKSSU3EM* target = dynamic_cast<CFieldFermionKSSU3EM*>(U);
    if (NULL != target)
    {
        target->m_fa2Ez = m_fa2Ez;
        target->m_fa2Bz = m_fa2Bz;
    }
}

CCString CFieldFermionKSSU3EM::GetInfos(const CCString& tab) const
{
    CCString sRet = tab + _T("Name : CFieldFermionKSSU3R\n");
    sRet = sRet + tab + _T("Mass (2am) : ") + appFloatToString(m_f2am) + _T("\n");
    sRet = sRet + tab + _T("MD Rational (c) : ") + appFloatToString(m_rMD.m_fC) + _T("\n");
    sRet = sRet + tab + _T("MC Rational (c) : ") + appFloatToString(m_rMC.m_fC) + _T("\n");
    sRet = sRet + tab + _T("Q x a^2Ez : ") + appFloatToString(m_fa2Ez) + _T("\n");
    sRet = sRet + tab + _T("Q x a^2Bz : ") + appFloatToString(m_fa2Bz) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================