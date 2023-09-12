//=============================================================================
// FILENAME : CFieldFermionWilsonSU3DEven.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/18/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSU3DREven)

#pragma region DOperator

#pragma region Step1

/**
* Step1 is phi_odd = Doo^{-1} Doe phi_even,
* in phi_even = Dee phi_even - Deo (Doo^{-1} Doe phi_even)
* Note: for M = c D, it neglect all the coefficient in Doo^{-1} Doe.
* Also, because gamma5 ^2 = 1, for M = Ddagger = gamma5 D gamma5,
* then it is: gamma5 (Deo Doo^{-1} Doe) gamma5, so we keep only the gamma5 on the right
*
* Note: For Naive Expoenential, Doo = 1, so it is Doe gamma5 phi_even
* Note: The first step has no coefficients
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_X_Step1(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    intokernalInt4_odd;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix gammaMu = __chiralGamma[GAMMA1];

    Real fYOmega = static_cast<Real>(sSite4.y - sCenter.y)* fOmega;

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 0);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 0);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);

    if (bDDagger)
    {
        x_p_mu_Fermion_element = gamma5.MulWilsonC(x_p_mu_Fermion_element);
        x_m_mu_Fermion_element = gamma5.MulWilsonC(x_m_mu_Fermion_element);
    }

    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);

    if (x_p_mu_Fermion.NeedToOpposite())
    {
        //printf("OppositeX x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        result.Sub(u_phi_x_p_m);
        result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fYOmega);
        result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        //==========================
        // (1 - Ga_x - y O Ga_4)
        result.Add(u_phi_x_p_m);
        //- Ga_x U(x,mu) phi(x+ mu)
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fYOmega);
        // - y O Ga_4
        result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    }


    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(u_dagger_phi_x_m_m);
        result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fYOmega);
        result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        //==========================
        // (1 + Ga_x + y O Ga_4)
        result.Add(u_dagger_phi_x_m_m);
        // + Ga_x U^{dagger}(x-mu) phi(x-mu)
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fYOmega);
        // + y O Ga_4 U^{dagger}(x-mu) phi(x-mu)
        result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }

    //this is added as pResultData[uiSiteIndex].Sub(result);
    //However, now it is pResultData[uiSiteIndex] = -result;
    //So here we multiply -kai
    //result.MulReal(kai);
    result.MulReal(-kai);

    //We only gamma5 the right vector
    //if (bDDagger)
    //{
    //    result = gamma5.MulWilsonC(result);
    //}
    pResultData[uiSiteIndex] = result;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_Y_Step1(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    intokernalInt4_odd;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix gammaMu = __chiralGamma[GAMMA2];

    Real fXOmega = static_cast<Real>(sSite4.x - sCenter.x)* fOmega;

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 1);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 1);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);

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
        //printf("OppositeY x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        result.Sub(u_phi_x_p_m);
        result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fXOmega);
        result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        result.Add(u_phi_x_p_m);
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fXOmega);
        result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(u_dagger_phi_x_m_m);
        result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fXOmega);
        result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        result.Add(u_dagger_phi_x_m_m);
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fXOmega);
        result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    result.MulReal(kai);

    pResultData[uiSiteIndex].Sub(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_Z_Step1(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    intokernalInt4_odd;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gammaMu = __chiralGamma[GAMMA3];

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 2);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 2);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);

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
        //printf("OppositeZ x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        result.Sub(u_phi_x_p_m);
        //- gammamu U(x,mu) phi(x+ mu)
        result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        result.Add(u_phi_x_p_m);
        //- gammamu U(x,mu) phi(x+ mu)
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(u_dagger_phi_x_m_m);
        //gammamu U^{dagger}(x-mu) phi(x-mu)
        result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        result.Add(u_dagger_phi_x_m_m);
        //gammamu U^{dagger}(x-mu) phi(x-mu)
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
    }

    result.MulReal(kai);
    pResultData[uiSiteIndex].Sub(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_T_Step1(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    intokernalInt4_odd;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix sigma12 = __chiralGamma[SIGMA12E];
    Real fhalfOmega = fOmega * F(0.5);
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 3);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);

    if (bDDagger)
    {
        x_p_mu_Fermion_element = gamma5.MulWilsonC(x_p_mu_Fermion_element);
        x_m_mu_Fermion_element = gamma5.MulWilsonC(x_m_mu_Fermion_element);
    }

    //hopping terms

    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
    u_phi_x_p_m.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    deviceWilsonVectorSU3 cospart = u_phi_x_p_m.MulRealC(kai * (_cos(fhalfOmega) /* - F(1.0)*/));
    //sinpart 
    u_phi_x_p_m.MulComp(_make_cuComplex(F(0.0), kai * _sin(fhalfOmega)));
    cospart.Add(sigma12.MulWilsonC(u_phi_x_p_m));

    if (x_p_mu_Fermion.NeedToOpposite())
    {
        //printf("OppositeT x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        //k(cos_m_1 - i sin sigma12)(1-gamma4)
        result.Sub(cospart);
    }
    else
    {
        //-k(cos_m_1 - i sin sigma12)(1-gamma4)
        result.Add(cospart);
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);

    u_dagger_phi_x_m_m.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    deviceWilsonVectorSU3 cospart2 = u_dagger_phi_x_m_m.MulRealC(kai * (_cos(fhalfOmega) /*- F(1.0)*/));
    //sinpart 
    u_dagger_phi_x_m_m.MulComp(_make_cuComplex(F(0.0), -kai * _sin(fhalfOmega)));
    cospart2.Add(sigma12.MulWilsonC(u_dagger_phi_x_m_m));

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        //-k(cos_m_1+isin)(1+gamma4)
        result.Sub(cospart2);
    }
    else
    {
        result.Add(cospart2);
    }

    //result = phi(x) - kai sum _mu result
    //result.MulReal(kai);
    pResultData[uiSiteIndex].Sub(result);
}

#pragma endregion

#pragma region Step2

/**
 * phi_even = Dee phi_even - Deo phi_odd
 * If we have coefficient, then, it is c (Dee phi_even - Deo phi_odd)
 * If we have dagger, then, it is c (Dee phi_even - gamma5 Deo phi_odd)
 *
 * For Naive and Exp, Dee = 1
 * it is c (phi_even - gamma5 Deo phi_odd)
 * 
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_X_Step2(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4_even;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    //For the first kernel only
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix gammaMu = __chiralGamma[GAMMA1];

    Real fYOmega = static_cast<Real>(sSite4.y - sCenter.y) * fOmega;

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 0);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 0);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pResultData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pResultData, x_m_mu_Fermion, byFieldId);

    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
    
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        //printf("OppositeX x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        result.Sub(u_phi_x_p_m);
        result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fYOmega);
        result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        //==========================
        // (1 - Ga_x - y O Ga_4)
        result.Add(u_phi_x_p_m);
        //- Ga_x U(x,mu) phi(x+ mu)
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fYOmega);
        // - y O Ga_4
        result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    }


    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(u_dagger_phi_x_m_m);
        result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fYOmega);
        result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        //==========================
        // (1 + Ga_x + y O Ga_4)
        result.Add(u_dagger_phi_x_m_m);
        // + Ga_x U^{dagger}(x-mu) phi(x-mu)
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fYOmega);
        // + y O Ga_4 U^{dagger}(x-mu) phi(x-mu)
        result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }

    result.MulReal(kai);

    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }
    //Change to Add, because it is -Deo
    //Also, change the position, so coefficients can be multiplied for all but not only result
    //Note!! The change of position is for the first kernel only!
    pResultData[uiSiteIndex].Add(result);

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

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_Y_Step2(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4_even;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix gammaMu = __chiralGamma[GAMMA2];

    Real fXOmega = static_cast<Real>(sSite4.x - sCenter.x) * fOmega;

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 1);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 1);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pResultData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pResultData, x_m_mu_Fermion, byFieldId);

    //hopping terms

    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        //printf("OppositeY x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        result.Sub(u_phi_x_p_m);
        result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fXOmega);
        result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        result.Add(u_phi_x_p_m);
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fXOmega);
        result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(u_dagger_phi_x_m_m);
        result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fXOmega);
        result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        result.Add(u_dagger_phi_x_m_m);
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fXOmega);
        result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    result.MulReal(kai);

    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        result.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        result.MulComp(cCoeff);
        break;
    }
    //Change to Add
    pResultData[uiSiteIndex].Add(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_Z_Step2(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4_even;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gammaMu = __chiralGamma[GAMMA3];

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 2);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 2);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pResultData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pResultData, x_m_mu_Fermion, byFieldId);

    //hopping terms

    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        //printf("OppositeZ x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        result.Sub(u_phi_x_p_m);
        //- gammamu U(x,mu) phi(x+ mu)
        result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        result.Add(u_phi_x_p_m);
        //- gammamu U(x,mu) phi(x+ mu)
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(u_dagger_phi_x_m_m);
        //gammamu U^{dagger}(x-mu) phi(x-mu)
        result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        result.Add(u_dagger_phi_x_m_m);
        //gammamu U^{dagger}(x-mu) phi(x-mu)
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
    }

    result.MulReal(kai);
    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }

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

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_T_Step2(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4_even;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix sigma12 = __chiralGamma[SIGMA12E];
    Real fhalfOmega = fOmega * F(0.5);
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 3);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pResultData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pResultData, x_m_mu_Fermion, byFieldId);

    //hopping terms

    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
    u_phi_x_p_m.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    deviceWilsonVectorSU3 cospart = u_phi_x_p_m.MulRealC(kai * (_cos(fhalfOmega) /* - F(1.0)*/));
    //sinpart 
    u_phi_x_p_m.MulComp(_make_cuComplex(F(0.0), kai * _sin(fhalfOmega)));
    cospart.Add(sigma12.MulWilsonC(u_phi_x_p_m));

    if (x_p_mu_Fermion.NeedToOpposite())
    {
        //printf("OppositeT x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        //k(cos_m_1 - i sin sigma12)(1-gamma4)
        result.Sub(cospart);
    }
    else
    {
        //-k(cos_m_1 - i sin sigma12)(1-gamma4)
        result.Add(cospart);
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);

    u_dagger_phi_x_m_m.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    deviceWilsonVectorSU3 cospart2 = u_dagger_phi_x_m_m.MulRealC(kai * (_cos(fhalfOmega) /*- F(1.0)*/));
    //sinpart 
    u_dagger_phi_x_m_m.MulComp(_make_cuComplex(F(0.0), -kai * _sin(fhalfOmega)));
    cospart2.Add(sigma12.MulWilsonC(u_dagger_phi_x_m_m));

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        //-k(cos_m_1+isin)(1+gamma4)
        result.Sub(cospart2);
    }
    else
    {
        result.Add(cospart2);
    }

    //result = phi(x) - kai sum _mu result
    //result.MulReal(kai);
    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }

    //==============================
    //res = [gamma5 (orig - term3 - result) gamma5] * coeff
    //    = [gamma5 orig gamma5] * coeff - [gamma5 (term3 + result) gamma5] * coeff
    //    = res0 - [gamma5 (term3 + result) gamma5] * coeff
    //term3.Add(result);
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

void CFieldFermionWilsonSU3DREven::DOperator(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceWilsonVectorSU3* pTarget = (deviceWilsonVectorSU3*)pTargetBuffer;
    const deviceWilsonVectorSU3* pSource = (deviceWilsonVectorSU3*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    
    preparethread_even;

#pragma region Step1

    _kernelDFermionWilsonSquareSU3_DR_Exponential_X_Step1 << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pTarget,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDagger);
    _kernelDFermionWilsonSquareSU3_DR_Exponential_Y_Step1 << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pTarget,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDagger);
    _kernelDFermionWilsonSquareSU3_DR_Exponential_Z_Step1 << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pTarget,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDagger);
    _kernelDFermionWilsonSquareSU3_DR_Exponential_T_Step1 << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pTarget,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDagger);

#pragma endregion

#pragma region Step1

    _kernelDFermionWilsonSquareSU3_DR_Exponential_X_Step2 << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pTarget,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
    _kernelDFermionWilsonSquareSU3_DR_Exponential_Y_Step2 << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pTarget,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
    _kernelDFermionWilsonSquareSU3_DR_Exponential_Z_Step2 << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pTarget,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
    _kernelDFermionWilsonSquareSU3_DR_Exponential_T_Step2 << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pTarget,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

#pragma endregion

    SetOddZero(pTarget);
}


#pragma region Kernel

#pragma region Write Even

/**
 * - Deo Doo^{-1} psi_odd
 * For Wilson-Dirac, Doo=1, so it is just
 * - Deo psi_odd
 *
 * For Ddagger
 *  - gamma5 Deo Doo^{-1} gamma5 psi_odd
 *  = - gamma5 Deo gamma5 psi_odd
 *
 * very similar as Step2 in D operator, but is different
 * For D operator, it is - gamma5 Deo psi_odd
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_X_WriteEven(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    intokernalInt4_even;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    //For the first kernel only
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix gammaMu = __chiralGamma[GAMMA1];

    Real fYOmega = static_cast<Real>(sSite4.y - sCenter.y)* fOmega;

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 0);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 0);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);
    if (bDDagger)
    {
        x_p_mu_Fermion_element = gamma5.MulWilsonC(x_p_mu_Fermion_element);
        x_m_mu_Fermion_element = gamma5.MulWilsonC(x_m_mu_Fermion_element);
    }

    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);

    if (x_p_mu_Fermion.NeedToOpposite())
    {
        //printf("OppositeX x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        result.Sub(u_phi_x_p_m);
        result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fYOmega);
        result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        //==========================
        // (1 - Ga_x - y O Ga_4)
        result.Add(u_phi_x_p_m);
        //- Ga_x U(x,mu) phi(x+ mu)
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fYOmega);
        // - y O Ga_4
        result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    }


    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(u_dagger_phi_x_m_m);
        result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fYOmega);
        result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        //==========================
        // (1 + Ga_x + y O Ga_4)
        result.Add(u_dagger_phi_x_m_m);
        // + Ga_x U^{dagger}(x-mu) phi(x-mu)
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fYOmega);
        // + y O Ga_4 U^{dagger}(x-mu) phi(x-mu)
        result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }

    result.MulReal(kai);

    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }
    //Change to Add, because it is -Deo
    //Also, change the position, so coefficients can be multiplied for all but not only result
    //Note!! The change of position is for the first kernel only!
    pResultData[uiSiteIndex].Add(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_Y_WriteEven(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    intokernalInt4_even;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix gammaMu = __chiralGamma[GAMMA2];

    Real fXOmega = static_cast<Real>(sSite4.x - sCenter.x)* fOmega;

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 1);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 1);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);
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
        //printf("OppositeY x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        result.Sub(u_phi_x_p_m);
        result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fXOmega);
        result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        result.Add(u_phi_x_p_m);
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fXOmega);
        result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(u_dagger_phi_x_m_m);
        result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fXOmega);
        result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        result.Add(u_dagger_phi_x_m_m);
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fXOmega);
        result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    result.MulReal(kai);

    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }

    //Change to Add
    pResultData[uiSiteIndex].Add(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_Z_WriteEven(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    intokernalInt4_even;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gammaMu = __chiralGamma[GAMMA3];

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 2);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 2);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);
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
        //printf("OppositeZ x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        result.Sub(u_phi_x_p_m);
        //- gammamu U(x,mu) phi(x+ mu)
        result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        result.Add(u_phi_x_p_m);
        //- gammamu U(x,mu) phi(x+ mu)
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(u_dagger_phi_x_m_m);
        //gammamu U^{dagger}(x-mu) phi(x-mu)
        result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        result.Add(u_dagger_phi_x_m_m);
        //gammamu U^{dagger}(x-mu) phi(x-mu)
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
    }

    result.MulReal(kai);
    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }

    pResultData[uiSiteIndex].Add(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_T_WriteEven(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    intokernalInt4_even;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix sigma12 = __chiralGamma[SIGMA12E];
    Real fhalfOmega = fOmega * F(0.5);
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 3);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);
    if (bDDagger)
    {
        x_p_mu_Fermion_element = gamma5.MulWilsonC(x_p_mu_Fermion_element);
        x_m_mu_Fermion_element = gamma5.MulWilsonC(x_m_mu_Fermion_element);
    }
    //hopping terms

    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
    u_phi_x_p_m.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    deviceWilsonVectorSU3 cospart = u_phi_x_p_m.MulRealC(kai * (_cos(fhalfOmega) /* - F(1.0)*/));
    //sinpart 
    u_phi_x_p_m.MulComp(_make_cuComplex(F(0.0), kai * _sin(fhalfOmega)));
    cospart.Add(sigma12.MulWilsonC(u_phi_x_p_m));

    if (x_p_mu_Fermion.NeedToOpposite())
    {
        //printf("OppositeT x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        //k(cos_m_1 - i sin sigma12)(1-gamma4)
        result.Sub(cospart);
    }
    else
    {
        //-k(cos_m_1 - i sin sigma12)(1-gamma4)
        result.Add(cospart);
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);

    u_dagger_phi_x_m_m.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    deviceWilsonVectorSU3 cospart2 = u_dagger_phi_x_m_m.MulRealC(kai * (_cos(fhalfOmega) /*- F(1.0)*/));
    //sinpart 
    u_dagger_phi_x_m_m.MulComp(_make_cuComplex(F(0.0), -kai * _sin(fhalfOmega)));
    cospart2.Add(sigma12.MulWilsonC(u_dagger_phi_x_m_m));

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        //-k(cos_m_1+isin)(1+gamma4)
        result.Sub(cospart2);
    }
    else
    {
        result.Add(cospart2);
    }

    //result = phi(x) - kai sum _mu result
    //result.MulReal(kai);
    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }

    //==============================
    //res = [gamma5 (orig - term3 - result) gamma5] * coeff
    //    = [gamma5 orig gamma5] * coeff - [gamma5 (term3 + result) gamma5] * coeff
    //    = res0 - [gamma5 (term3 + result) gamma5] * coeff
    //term3.Add(result);
    pResultData[uiSiteIndex].Add(result);
}

#pragma endregion

/**
 * z_even obtained, we calculate D_oo^{-1}(phi_odd - D_oe z_even)
 * Similar as Step1, but we need to keep the left gamma5
 *
 * For the first kernel, we also write z_even
 *
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_X_WriteBack(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    UINT uiEvenIndex = (((threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z)) << 1);
    UINT uiSiteIndex = uiEvenIndex + 1;
    SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);
    if (!sSite4.IsOdd())
    {
        uiSiteIndex = uiSiteIndex - 1;
        uiEvenIndex = uiEvenIndex + 1;
        sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);
    }

    //write z_even
    pResultData[uiEvenIndex] = pDeviceData[uiEvenIndex];

    //====================
    //Now we only care the odd sites
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix gammaMu = __chiralGamma[GAMMA1];

    Real fYOmega = static_cast<Real>(sSite4.y - sCenter.y)* fOmega;

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 0);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 0);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);

    if (bDDagger)
    {
        x_p_mu_Fermion_element = gamma5.MulWilsonC(x_p_mu_Fermion_element);
        x_m_mu_Fermion_element = gamma5.MulWilsonC(x_m_mu_Fermion_element);
    }

    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);

    if (x_p_mu_Fermion.NeedToOpposite())
    {
        //printf("OppositeX x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        result.Sub(u_phi_x_p_m);
        result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fYOmega);
        result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        //==========================
        // (1 - Ga_x - y O Ga_4)
        result.Add(u_phi_x_p_m);
        //- Ga_x U(x,mu) phi(x+ mu)
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fYOmega);
        // - y O Ga_4
        result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(u_dagger_phi_x_m_m);
        result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fYOmega);
        result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        //==========================
        // (1 + Ga_x + y O Ga_4)
        result.Add(u_dagger_phi_x_m_m);
        // + Ga_x U^{dagger}(x-mu) phi(x-mu)
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fYOmega);
        // + y O Ga_4 U^{dagger}(x-mu) phi(x-mu)
        result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }

    result.MulReal(kai);

    //We only gamma5 the right vector
    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }
    pResultData[uiSiteIndex].Add(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_Y_WriteBack(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    intokernalInt4_odd;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix gammaMu = __chiralGamma[GAMMA2];

    Real fXOmega = static_cast<Real>(sSite4.x - sCenter.x)* fOmega;

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 1);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 1);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);

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
        //printf("OppositeY x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        result.Sub(u_phi_x_p_m);
        result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fXOmega);
        result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        result.Add(u_phi_x_p_m);
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fXOmega);
        result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(u_dagger_phi_x_m_m);
        result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fXOmega);
        result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        result.Add(u_dagger_phi_x_m_m);
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fXOmega);
        result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    result.MulReal(kai);
    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }

    pResultData[uiSiteIndex].Add(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_Z_WriteBack(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    intokernalInt4_odd;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gammaMu = __chiralGamma[GAMMA3];

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 2);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 2);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);

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
        //printf("OppositeZ x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        result.Sub(u_phi_x_p_m);
        //- gammamu U(x,mu) phi(x+ mu)
        result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        result.Add(u_phi_x_p_m);
        //- gammamu U(x,mu) phi(x+ mu)
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(u_dagger_phi_x_m_m);
        //gammamu U^{dagger}(x-mu) phi(x-mu)
        result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        result.Add(u_dagger_phi_x_m_m);
        //gammamu U^{dagger}(x-mu) phi(x-mu)
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
    }

    result.MulReal(kai);
    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }
    pResultData[uiSiteIndex].Add(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_T_WriteBack(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    intokernalInt4_odd;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix sigma12 = __chiralGamma[SIGMA12E];
    Real fhalfOmega = fOmega * F(0.5);
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 3);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);

    if (bDDagger)
    {
        x_p_mu_Fermion_element = gamma5.MulWilsonC(x_p_mu_Fermion_element);
        x_m_mu_Fermion_element = gamma5.MulWilsonC(x_m_mu_Fermion_element);
    }

    //hopping terms

    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
    u_phi_x_p_m.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    deviceWilsonVectorSU3 cospart = u_phi_x_p_m.MulRealC(kai * (_cos(fhalfOmega) /* - F(1.0)*/));
    //sinpart 
    u_phi_x_p_m.MulComp(_make_cuComplex(F(0.0), kai * _sin(fhalfOmega)));
    cospart.Add(sigma12.MulWilsonC(u_phi_x_p_m));

    if (x_p_mu_Fermion.NeedToOpposite())
    {
        //printf("OppositeT x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        //k(cos_m_1 - i sin sigma12)(1-gamma4)
        result.Sub(cospart);
    }
    else
    {
        //-k(cos_m_1 - i sin sigma12)(1-gamma4)
        result.Add(cospart);
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);

    u_dagger_phi_x_m_m.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    deviceWilsonVectorSU3 cospart2 = u_dagger_phi_x_m_m.MulRealC(kai * (_cos(fhalfOmega) /*- F(1.0)*/));
    //sinpart 
    u_dagger_phi_x_m_m.MulComp(_make_cuComplex(F(0.0), -kai * _sin(fhalfOmega)));
    cospart2.Add(sigma12.MulWilsonC(u_dagger_phi_x_m_m));

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        //-k(cos_m_1+isin)(1+gamma4)
        result.Sub(cospart2);
    }
    else
    {
        result.Add(cospart2);
    }

    //result = phi(x) - kai sum _mu result
    //result.MulReal(kai);

    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }
    
    pResultData[uiSiteIndex].Add(result);
}

#pragma endregion


void CFieldFermionWilsonSU3DREven::CopyTo(CField* U) const
{
    CFieldFermionWilsonSU3DEven::CopyTo(U);
}

/**
 * phi_even = psi_even - Deo Doo^{-1} psi_odd
 */
void CFieldFermionWilsonSU3DREven::WriteEvenSites(const CFieldFermion* x, const CFieldGauge* pGauge, UBOOL bDagger)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3DR* pField = dynamic_cast<const CFieldFermionWilsonSquareSU3DR*>(x);
    if (NULL == pField)
    {
        appCrucial(_T("CFieldFermionWilsonSU3DREven can only play with gauge SU3DR!\n"));
        return;
    }

    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    m_byParentId = pField->m_byFieldId;

    //appGeneral(_T("----kai:%f, Omega:%f, CenterX:%d, CenterY:%d\n"), CCommonData::m_fKai, CCommonData::m_fOmega, CCommonData::m_sCenter.x, CCommonData::m_sCenter.y);

    preparethread_even;
    _kernelDFermionWilsonSquareSU3_DR_Exponential_X_WriteEven << <block, threads >> > (
        pField->m_pDeviceData,
        pFieldSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        m_pDeviceData,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDagger);
    _kernelDFermionWilsonSquareSU3_DR_Exponential_Y_WriteEven << <block, threads >> > (
        pField->m_pDeviceData,
        pFieldSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        m_pDeviceData,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDagger);
    _kernelDFermionWilsonSquareSU3_DR_Exponential_Z_WriteEven << <block, threads >> > (
        pField->m_pDeviceData,
        pFieldSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        m_pDeviceData,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDagger);
    _kernelDFermionWilsonSquareSU3_DR_Exponential_T_WriteEven << <block, threads >> > (
        pField->m_pDeviceData,
        pFieldSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        m_pDeviceData,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDagger);


    if (!pField->m_bNaive || !pField->m_bExponential)
    {
        appCrucial(_T("Error: CFieldFermionWilsonSquareSU3DREven only work with Naive and Exponential!\n"));
    }

    if (0 != (_HC_Lt & 1) || 0 != (_HC_Lt & 1))
    {
        appGeneral(_T("Warning: Even-Odd work only for even extent!\n"));
    }

    //Step2 set odd zero
    SetOddZero(m_pDeviceData);

    if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
    {
        m_fLength = Dot(this).x;
    }
}

/**
* z_even obtained, we calculate D_oo^{-1}(phi_odd - D_oe z_even)
*/
void CFieldFermionWilsonSU3DREven::WriteBackEvenSites(CFieldFermion* pParentField, const CFieldGauge* pGauge, UBOOL bDdagger) const
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    if (NULL == pParentField || EFT_FermionWilsonSquareSU3 != pParentField->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3DEven can only play with gauge CFieldFermionWilsonSquareSU3!"));
        return;
    }
    CFieldFermionWilsonSquareSU3DR* pParentFieldSU3 = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(pParentField);
    if (NULL == pParentFieldSU3)
    {
        appCrucial(_T("CFieldFermionWilsonSU3DREven can only play with gauge SU3DR!"));
        return;
    }

    preparethread_even;
    _kernelDFermionWilsonSquareSU3_DR_Exponential_X_WriteBack << <block, threads >> > (
        m_pDeviceData,
        pFieldSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pParentFieldSU3->m_pDeviceData,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDdagger);
    _kernelDFermionWilsonSquareSU3_DR_Exponential_Y_WriteBack << <block, threads >> > (
        m_pDeviceData,
        pFieldSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pParentFieldSU3->m_pDeviceData,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDdagger);
    _kernelDFermionWilsonSquareSU3_DR_Exponential_Z_WriteBack << <block, threads >> > (
        m_pDeviceData,
        pFieldSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pParentFieldSU3->m_pDeviceData,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDdagger);
    _kernelDFermionWilsonSquareSU3_DR_Exponential_T_WriteBack << <block, threads >> > (
        m_pDeviceData,
        pFieldSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pParentFieldSU3->m_pDeviceData,
        CCommonData::m_fKai,
        CCommonData::m_fOmega,
        CCommonData::m_sCenter,
        m_byParentId,
        bDdagger);
}

CCString CFieldFermionWilsonSU3DREven::GetInfos(const CCString& tab) const
{
    CCString sRet = tab + _T("Name : CFieldFermionWilsonSU3DREven\n");
    sRet = sRet + tab + _T("Hopping : ") + appFloatToString(CCommonData::m_fKai) + _T("\n");

    SSmallInt4 boundary = appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(m_byFieldId);
    sRet = sRet + tab + _T("boundary : [") + appIntToString(boundary.x) + _T(", ") + appIntToString(boundary.y) + _T(", ") + appIntToString(boundary.z) + _T(", ") + appIntToString(boundary.w) + _T("]\n");

    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================