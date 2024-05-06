//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3Acc.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [07/27/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSquareSU3Acc)

#pragma region DOperator

#pragma region kernel

#pragma region hamitonian

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_Acc(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    Real fG,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    const gammaMatrix & gamma5 = __chiralGamma[GAMMA5];
    const gammaMatrix & gamma4 = __chiralGamma[GAMMA4];
    const Real fGT = fG * sSite4.w;
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];
    if (bDDagger)
    {
        pResultData[uiSiteIndex] = gamma5.MulWilsonC(pResultData[uiSiteIndex]);
    }

    //idir = mu
    #pragma unroll
    for (UINT idir = 0; idir < 4; ++idir)
    {
        //Get Gamma mu
        const gammaMatrix & gammaMu = __chiralGamma[GAMMA1 + idir];

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        deviceSU3 x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

        deviceWilsonVectorSU3 x_p_mu_Fermion_element = pDeviceData[x_p_mu_Fermion.m_uiSiteIndex];
        deviceWilsonVectorSU3 x_m_mu_Fermion_element = pDeviceData[x_m_mu_Fermion.m_uiSiteIndex];

        if (bDDagger)
        {
            x_p_mu_Fermion_element = gamma5.MulWilsonC(x_p_mu_Fermion_element);
            x_m_mu_Fermion_element = gamma5.MulWilsonC(x_m_mu_Fermion_element);
        }

        //hopping terms

        //U(x,mu) phi(x+ mu)
        deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            result.Sub(u_phi_x_p_m);

            if (2 == idir)
            {
                result.Sub(gamma4.MulWilsonC(u_phi_x_p_m.MulRealC(fGT)));
            }

            //- gammamu U(x,mu) phi(x+ mu)
            result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        }
        else
        {
            result.Add(u_phi_x_p_m);

            if (2 == idir)
            {
                result.Add(gamma4.MulWilsonC(u_phi_x_p_m.MulRealC(fGT)));
            }
            //- gammamu U(x,mu) phi(x+ mu)
            result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            result.Sub(u_dagger_phi_x_m_m);

            if (2 == idir)
            {
                result.Add(gamma4.MulWilsonC(u_phi_x_p_m.MulRealC(fGT)));
            }

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
        else
        {
            result.Add(u_dagger_phi_x_m_m);

            if (2 == idir)
            {
                result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m.MulRealC(fGT)));
            }

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
    }
    
    //result = phi(x) - kai sum _mu result
    result.MulReal(kai);
    pResultData[uiSiteIndex].Sub(result);

    if (bDDagger)
    {
        pResultData[uiSiteIndex] = gamma5.MulWilsonC(pResultData[uiSiteIndex]);
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

#pragma region force

__global__ void _CLG_LAUNCH_BOUND
_kernelDWilsonForceSU3_Acc(
    const deviceWilsonVectorSU3* __restrict__ pInverseD,
    const deviceWilsonVectorSU3* __restrict__ pInverseDDdagger,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    deviceSU3* pForce,
    Real fKai,
    Real fG,
    BYTE byFieldId)
{
    intokernalInt4;

    const gammaMatrix& gamma4 = __chiralGamma[GAMMA4];
    const Real fGT = fG * sSite4.w;

    const deviceWilsonVectorSU3 x_Left(pInverseDDdagger[uiSiteIndex]);
    const deviceWilsonVectorSU3 x_Right(pInverseD[uiSiteIndex]);

    //idir = mu
    #pragma unroll
    for (UINT idir = 0; idir < 4; ++idir)
    {
        //Get Gamma mu
        gammaMatrix gammaMu = __chiralGamma[GAMMA1 + idir];

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        SIndex x_p_mu_Fermion = pFermionMove[linkIndex * 2]; // __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, (idir + 1));

        const deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
        deviceWilsonVectorSU3 x_p_mu_Left(pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex]);

        deviceSU3 x_Gauge_element = pGauge[linkIndex];

        deviceWilsonVectorSU3 right1(x_p_mu_Right);
        right1.Sub(gammaMu.MulWilsonC(right1));
        if (2 == idir)
        {
            right1.Add(gamma4.MulWilsonC(x_p_mu_Right).MulRealC(fGT));
        }
        deviceSU3 mid = deviceSU3::makeSU3Contract(x_Left, right1);

        deviceWilsonVectorSU3 right2(x_Right);
        right2.Add(gammaMu.MulWilsonC(right2));
        if (2 == idir)
        {
            right2.Sub(gamma4.MulWilsonC(x_Right).MulRealC(fGT));
        }
        mid.Add(deviceSU3::makeSU3Contract(right2, x_p_mu_Left));

        deviceSU3 forceOfThisLink = x_Gauge_element.MulC(mid);
        forceOfThisLink.Ta();
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            forceOfThisLink.MulReal(-fKai);
        }
        else
        {
            forceOfThisLink.MulReal(fKai);
        }

        pForce[linkIndex].Add(forceOfThisLink);
    }
}

#pragma endregion

#pragma endregion

/**
* It will not produce correct result unless _CLG_USE_LAUNCH_BOUND is set to 1
* It seems max-registor count is reached without throw an error
* To avoid it, we use slower method, split it into two functions.
*/
void CFieldFermionWilsonSquareSU3Acc::DOperator(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, 
    UBOOL bDagger, EOperatorCoefficientType eOCT, 
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceWilsonVectorSU3* pTarget = (deviceWilsonVectorSU3*)pTargetBuffer;
    const deviceWilsonVectorSU3* pSource = (deviceWilsonVectorSU3*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    preparethread;
    _kernelDFermionWilsonSquareSU3_Acc << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        pTarget,
        m_fKai,
        CCommonData::m_fG,
        m_byFieldId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

void CFieldFermionWilsonSquareSU3Acc::DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const
{
    deviceSU3* pForceSU3 = (deviceSU3*)pForce;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
    const deviceWilsonVectorSU3* pDphiBuffer = (deviceWilsonVectorSU3*)pDphi;
    const deviceWilsonVectorSU3* pDDphiBuffer = (deviceWilsonVectorSU3*)pDDphi;

    preparethread;
    _kernelDWilsonForceSU3_Acc << <block, threads >> > (
        pDphiBuffer,
        pDDphiBuffer,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        pForceSU3,
        m_fKai,
        CCommonData::m_fG,
        m_byFieldId);
}

#pragma endregion

void CFieldFermionWilsonSquareSU3Acc::CopyTo(CField* U) const
{
    CFieldFermionWilsonSquareSU3::CopyTo(U);
}

CCString CFieldFermionWilsonSquareSU3Acc::GetInfos(const CCString &tab) const
{
    CCString sRet = CFieldFermionWilsonSquareSU3::GetInfos(tab);
    sRet = sRet + tab + _T("Acc : ") + appToString(CCommonData::m_fG) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================