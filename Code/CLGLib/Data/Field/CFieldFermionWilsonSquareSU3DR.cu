//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3DR.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/19/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSquareSU3DR)

#pragma region DOperator

#pragma region kernel

#pragma region hamitonian

/**
 * This is only for no exponent D operator
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR(
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
    UBOOL bNaive,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    const gammaMatrix& gamma5 = __chiralGamma[GAMMA5];
    const gammaMatrix& gamma4 = __chiralGamma[GAMMA4];
    const Real fYkOmega = static_cast<Real>(sSite4.y - sCenter.y) * fOmega * kai;
    const Real fXkOmega = static_cast<Real>(sSite4.x - sCenter.x) * fOmega * kai;

    deviceWilsonVectorSU3 right_element = bDDagger ? gamma5.MulWilsonC(pDeviceData[uiSiteIndex]) : pDeviceData[uiSiteIndex];

    deviceWilsonVectorSU3 term3(__chiralGamma[SIGMA12E].MulWilsonC(right_element));
    term3 = gamma4.MulWilsonC(term3);
    term3.MulComp(_make_cuComplex(F(0.0), -fOmega * kai));

    if (!bNaive)
    {
        //res = res + 2 (y-x) k Omega right - i k gamma4 sigma12 right
        //-2(y-x) k Omega right
        right_element.MulReal(F(2.0) * (fYkOmega - fXkOmega));
        term3.Sub(right_element);
    }

    if (bDDagger)
    {
        term3 = gamma5.MulWilsonC(term3);
    }

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //idir = mu
    #pragma unroll
    for (UINT idir = 0; idir < 2; ++idir)
    {
        //=========================
        //get things

        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        //deviceSU3 x_Gauge_element = pGauge[linkIndex];
        const deviceSU3& x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir);
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
            if (0 == idir)
            {
                //- k y Omega x_p_m
                //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
                u_phi_x_p_m.MulReal(fYkOmega);
                if (!bNaive)
                {
                    result.Sub(u_phi_x_p_m);
                }
                result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_phi_x_p_m.MulReal(fXkOmega);
                if (!bNaive)
                {
                    result.Add(u_phi_x_p_m);
                }
                result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
            }
        }
        else
        {
            if (0 == idir)
            {
                //- k y Omega x_p_m
                //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
                u_phi_x_p_m.MulReal(fYkOmega);
                if (!bNaive)
                {
                    result.Add(u_phi_x_p_m);
                }
                result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_phi_x_p_m.MulReal(fXkOmega);
                if (!bNaive)
                {
                    result.Sub(u_phi_x_p_m);
                }
                result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
            }
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            if (0 == idir)
            {
                //- k y Omega x_p_m
                //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
                u_dagger_phi_x_m_m.MulReal(fYkOmega);
                if (!bNaive)
                {
                    result.Sub(u_dagger_phi_x_m_m);
                }
                result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_dagger_phi_x_m_m.MulReal(fXkOmega);
                if (!bNaive)
                {
                    result.Add(u_dagger_phi_x_m_m);
                }
                result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
            }
        }
        else
        {
            if (0 == idir)
            {
                //- k y Omega x_p_m
                //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
                u_dagger_phi_x_m_m.MulReal(fYkOmega);
                if (!bNaive)
                {
                    result.Add(u_dagger_phi_x_m_m);
                }
                result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_dagger_phi_x_m_m.MulReal(fXkOmega);
                if (!bNaive)
                {
                    result.Sub(u_dagger_phi_x_m_m);
                }
                result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
            }
        }
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
    term3.Add(result);
    switch (eCoeff)
    {
    case EOCT_Real:
        term3.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        term3.MulComp(cCoeff);
        break;
    }
    pResultData[uiSiteIndex].Sub(term3);
}

#if _CLG_ROTATING_NEW_IMP

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_0(
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
    UBOOL bNaive,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    if (bNaive)
    {
        pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];
        switch (eCoeff)
        {
        case EOCT_Real:
            pResultData[uiSiteIndex].MulReal(fCoeff);
            break;
        case EOCT_Complex:
            pResultData[uiSiteIndex].MulComp(cCoeff);
            break;
        }
        return;
    }

    const gammaMatrix& gamma5 = __chiralGamma[GAMMA5];
    deviceWilsonVectorSU3 right_element = bDDagger ?
        gamma5.MulWilsonC(pDeviceData[uiSiteIndex]) : pDeviceData[uiSiteIndex];

    deviceWilsonVectorSU3 term3 = deviceWilsonVectorSU3(right_element);

    const Real fkOmega = fOmega * kai;
    const Real fYkOmega = static_cast<Real>(sSite4.y - sCenter.y) * fkOmega;
    const Real fXkOmega = static_cast<Real>(sSite4.x - sCenter.x) * fkOmega;

    //res = res + 2 (y-x) k Omega right - i k gamma4 sigma12 right
    //-2(y-x) k Omega right
    right_element.MulReal(F(2.0) * (fYkOmega - fXkOmega));
    //term3.Sub(right_element);
    term3.Add(right_element);

    if (bDDagger)
    {
        term3 = gamma5.MulWilsonC(term3);
    }

    //==============================
    //res = [gamma5 (orig - term3 - result) gamma5] * coeff
    //    = [gamma5 orig gamma5] * coeff - [gamma5 (term3 + result) gamma5] * coeff
    //    = res0 - [gamma5 (term3 + result) gamma5] * coeff
    //term3.Add(result);
    switch (eCoeff)
    {
    case EOCT_Real:
        term3.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        term3.MulComp(cCoeff);
        break;
    }
    pResultData[uiSiteIndex] = term3;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_X(
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
    UBOOL bNaive,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    const gammaMatrix& gamma5 = __chiralGamma[GAMMA5];
    const gammaMatrix& gamma4 = __chiralGamma[GAMMA4];
    const gammaMatrix& gammaMu = __chiralGamma[GAMMA1];

    const Real fYOmega = static_cast<Real>(sSite4.y - sCenter.y) * fOmega;

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 0);

    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    const deviceSU3& x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 0);
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
        if (!bNaive)
        {
            result.Sub(u_phi_x_p_m);
        }
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
        if (!bNaive)
        {
            result.Add(u_phi_x_p_m);
        }
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
        if (!bNaive)
        {
            result.Sub(u_dagger_phi_x_m_m);
        }
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
        if (!bNaive)
        {
            result.Add(u_dagger_phi_x_m_m);
        }
        // + y O Ga_4 U^{dagger}(x-mu) phi(x-mu)
        result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
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
    pResultData[uiSiteIndex].Sub(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_Y(
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
    UBOOL bNaive,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    const gammaMatrix& gamma5 = __chiralGamma[GAMMA5];
    const gammaMatrix& gamma4 = __chiralGamma[GAMMA4];
    const gammaMatrix& gammaMu = __chiralGamma[GAMMA2];

    const Real fXOmega = static_cast<Real>(sSite4.x - sCenter.x) * fOmega;

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 1);

    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    const deviceSU3& x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 1);
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
        if (!bNaive)
        {
            result.Add(u_phi_x_p_m);
        }
        result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        result.Add(u_phi_x_p_m);
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        u_phi_x_p_m.MulReal(fXOmega);
        if (!bNaive)
        {
            result.Sub(u_phi_x_p_m);
        }
        result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Sub(u_dagger_phi_x_m_m);
        result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fXOmega);
        if (!bNaive)
        {
            result.Add(u_dagger_phi_x_m_m);
        }
        result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        result.Add(u_dagger_phi_x_m_m);
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        u_dagger_phi_x_m_m.MulReal(fXOmega);
        if (!bNaive)
        {
            result.Sub(u_dagger_phi_x_m_m);
        }
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
    pResultData[uiSiteIndex].Sub(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_Z(
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
    UBOOL bNaive,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    const gammaMatrix& gamma5 = __chiralGamma[GAMMA5];
    const gammaMatrix& gammaMu = __chiralGamma[GAMMA3];

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 2);

    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    const deviceSU3& x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 2);
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
    const deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
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
    const deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
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
    pResultData[uiSiteIndex].Sub(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_T(
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
    UBOOL bNaive,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    const gammaMatrix& gamma5 = __chiralGamma[GAMMA5];
    const gammaMatrix& gamma4 = __chiralGamma[GAMMA4];
    const gammaMatrix& sigma12 = __chiralGamma[SIGMA12E];
    const Real fhalfOmega = fOmega * F(0.5);
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);

    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    const deviceSU3& x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 3);
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
    switch (eCoeff)
    {
    case EOCT_Real:
        result.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        result.MulComp(cCoeff);
        break;
    }
    pResultData[uiSiteIndex].Sub(result);
}

#else

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_0(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    Real fOmega,
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger,
    UBOOL bNaive,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    Real fkOmega = fOmega * kai;
    Real fYkOmega = static_cast<Real>(sSite4.y - sCenter.y) * fkOmega;
    Real fXkOmega = static_cast<Real>(sSite4.x - sCenter.x) * fkOmega;

    deviceWilsonVectorSU3 right_element = bDDagger ? gamma5.MulWilsonC(pDeviceData[uiSiteIndex]) : pDeviceData[uiSiteIndex];

    deviceWilsonVectorSU3 term3 = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    if (!bNaive)
    {
        //res = res + 2 (y-x) k Omega right - i k gamma4 sigma12 right
        //-2(y-x) k Omega right
        right_element.MulReal(F(2.0) * (fYkOmega - fXkOmega));
        term3.Sub(right_element);
    }

    if (bDDagger)
    {
        term3 = gamma5.MulWilsonC(term3);
    }

    //==============================
    //res = [gamma5 (orig - term3 - result) gamma5] * coeff
    //    = [gamma5 orig gamma5] * coeff - [gamma5 (term3 + result) gamma5] * coeff
    //    = res0 - [gamma5 (term3 + result) gamma5] * coeff
    //term3.Add(result);
    switch (eCoeff)
    {
    case EOCT_Real:
        term3.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        term3.MulComp(cCoeff);
        break;
    }
    pResultData[uiSiteIndex].Sub(term3);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_X(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    Real fOmega,
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger,
    UBOOL bNaive,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];

    Real fYkOmega = static_cast<Real>(sSite4.y - sCenter.y) * fOmega * kai;
    //Real fXkOmega = static_cast<Real>(sSite4.x - sCenter.x) * fOmega * kai;

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    BYTE idir = 0;
    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    x_m_mu_Gauge_element.Dagger();

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
        //- k y Omega x_p_m
        //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
        u_phi_x_p_m.MulReal(fYkOmega);
        if (!bNaive)
        {
            result.Sub(u_phi_x_p_m);
        }
        result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        //- k y Omega x_p_m
        //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
        u_phi_x_p_m.MulReal(fYkOmega);
        if (!bNaive)
        {
            result.Add(u_phi_x_p_m);
        }
        result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        //- k y Omega x_p_m
        //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
        u_dagger_phi_x_m_m.MulReal(fYkOmega);
        if (!bNaive)
        {
            result.Sub(u_dagger_phi_x_m_m);
        }
        result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        //- k y Omega x_p_m
        //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
        u_dagger_phi_x_m_m.MulReal(fYkOmega);
        if (!bNaive)
        {
            result.Add(u_dagger_phi_x_m_m);
        }
        result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
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
    pResultData[uiSiteIndex].Sub(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_Y(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    Real fOmega,
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger,
    UBOOL bNaive,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];

    //Real fYkOmega = static_cast<Real>(sSite4.y - sCenter.y) * fOmega * kai;
    Real fXkOmega = static_cast<Real>(sSite4.x - sCenter.x) * fOmega * kai;

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    BYTE idir = 1;
    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    x_m_mu_Gauge_element.Dagger();

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
        //+ k x Omega x_p_m
        //- k x Omega gamma4 U(x,mu) phi(x+ mu)
        u_phi_x_p_m.MulReal(fXkOmega);
        if (!bNaive)
        {
            result.Add(u_phi_x_p_m);
        }
        result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    }
    else
    {
        //+ k x Omega x_p_m
        //- k x Omega gamma4 U(x,mu) phi(x+ mu)
        u_phi_x_p_m.MulReal(fXkOmega);
        if (!bNaive)
        {
            result.Sub(u_phi_x_p_m);
        }
        result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        //+ k x Omega x_p_m
        //- k x Omega gamma4 U(x,mu) phi(x+ mu)
        u_dagger_phi_x_m_m.MulReal(fXkOmega);
        if (!bNaive)
        {
            result.Add(u_dagger_phi_x_m_m);
        }
        result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    }
    else
    {
        //+ k x Omega x_p_m
        //- k x Omega gamma4 U(x,mu) phi(x+ mu)
        u_dagger_phi_x_m_m.MulReal(fXkOmega);
        if (!bNaive)
        {
            result.Sub(u_dagger_phi_x_m_m);
        }
        result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
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
    pResultData[uiSiteIndex].Sub(result);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential_T(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    Real fOmega,
    SSmallInt4 sCenter,
    BYTE byFieldId,
    UBOOL bDDagger,
    UBOOL bNaive,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

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
    BYTE idir = 3;
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

    SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    x_m_mu_Gauge_element.Dagger();

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
    deviceWilsonVectorSU3 cospart = u_phi_x_p_m.MulRealC(kai * (_cos(fhalfOmega) - F(1.0)));
    //sinpart 
    u_phi_x_p_m.MulComp(_make_cuComplex(F(0.0), -kai * _sin(fhalfOmega)));
    cospart.Add(sigma12.MulWilsonC(u_phi_x_p_m));

    if (x_p_mu_Fermion.NeedToOpposite())
    {
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
    deviceWilsonVectorSU3 cospart2 = u_dagger_phi_x_m_m.MulRealC(kai * (_cos(fhalfOmega) - F(1.0)));
    //sinpart 
    u_dagger_phi_x_m_m.MulComp(_make_cuComplex(F(0.0), kai * _sin(fhalfOmega)));
    cospart2.Add(sigma12.MulWilsonC(u_dagger_phi_x_m_m));

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        //-(cos_m_1+isin)(1+gamma4)
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
    pResultData[uiSiteIndex].Sub(result);
}

#endif

#pragma endregion

#pragma region force

__global__ void _CLG_LAUNCH_BOUND
_kernelDWilsonForceSU3_DR_X_Naive(
    const deviceWilsonVectorSU3* __restrict__ pInverseD,
    const deviceWilsonVectorSU3* __restrict__ pInverseDDdagger,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    deviceSU3* pForce,
    SSmallInt4 sCenter,
    Real fKai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    Real fYOmega = static_cast<Real>(sSite4.y - sCenter.y) * fOmega;
    const gammaMatrix& gamma4 = __chiralGamma[GAMMA4];

    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 0);
    const SIndex& x_p_mu_Fermion = pFermionMove[linkIndex * 2];

    //If one of the sites is on surface, it has no contribution.
    //Note that, the bond on surface is equivelant to both sites on surface.

    if (x_p_mu_Fermion.IsDirichlet() || sSite.IsDirichlet())
    {
        return;
    }

    if (x_p_mu_Fermion.NeedToOpposite())
    {
        fYOmega = -fYOmega * fKai;
    }
    else
    {
        fYOmega = fYOmega * fKai;
    }

    //====================
    //Get things
    const deviceWilsonVectorSU3& x_Left = pInverseDDdagger[uiSiteIndex];
    deviceWilsonVectorSU3 x_Right(pInverseD[uiSiteIndex]);
    //all not on surface
    deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
    const deviceWilsonVectorSU3& x_p_mu_Left = pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex];

    const deviceSU3& x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

    x_p_mu_Right.MulReal(-fYOmega);
    x_p_mu_Right = gamma4.MulWilsonC(x_p_mu_Right);

    deviceSU3 mid = deviceSU3::makeSU3Contract(x_Left, x_p_mu_Right);

    x_Right.MulReal(fYOmega);
    x_Right = gamma4.MulWilsonC(x_Right);
    mid.Add(deviceSU3::makeSU3Contract(x_Right, x_p_mu_Left));

    deviceSU3 forceOfThisLink = x_Gauge_element.MulC(mid);
    forceOfThisLink.Ta();
    pForce[linkIndex].Add(forceOfThisLink);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDWilsonForceSU3_DR_Y_Naive(
    const deviceWilsonVectorSU3* __restrict__ pInverseD,
    const deviceWilsonVectorSU3* __restrict__ pInverseDDdagger,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    deviceSU3* pForce,
    SSmallInt4 sCenter,
    Real fKai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    Real fXOmega = static_cast<Real>(sSite4.x - sCenter.x) * fOmega;
    const gammaMatrix& gamma4 = __chiralGamma[GAMMA4];

    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 1);
    const SIndex& x_p_mu_Fermion = pFermionMove[linkIndex * 2];

    //If one of the sites is on surface, it has no contribution.
    //Note that, the bond on surface is equivelant to both sites on surface.

    if (x_p_mu_Fermion.IsDirichlet() || sSite.IsDirichlet())
    {
        return;
    }

    if (x_p_mu_Fermion.NeedToOpposite())
    {
        fXOmega = -fXOmega * fKai;
    }
    else
    {
        fXOmega = fXOmega * fKai;
    }

    //====================
    //Get things
    const deviceWilsonVectorSU3& x_Left = pInverseDDdagger[uiSiteIndex];
    deviceWilsonVectorSU3 x_Right(pInverseD[uiSiteIndex]);
    //all not on surface
    deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
    const deviceWilsonVectorSU3& x_p_mu_Left = pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex];

    const deviceSU3& x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

    x_p_mu_Right.MulReal(fXOmega);
    x_p_mu_Right = gamma4.MulWilsonC(x_p_mu_Right);

    deviceSU3 mid = deviceSU3::makeSU3Contract(x_Left, x_p_mu_Right);

    x_Right.MulReal(-fXOmega);
    x_Right = gamma4.MulWilsonC(x_Right);
    mid.Add(deviceSU3::makeSU3Contract(x_Right, x_p_mu_Left));

    deviceSU3 forceOfThisLink = x_Gauge_element.MulC(mid);
    forceOfThisLink.Ta();
    pForce[linkIndex].Add(forceOfThisLink);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDWilsonForceSU3_DR_X(
    const deviceWilsonVectorSU3* __restrict__ pInverseD,
    const deviceWilsonVectorSU3* __restrict__ pInverseDDdagger,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    deviceSU3* pForce,
    SSmallInt4 sCenter,
    Real fKai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    Real fYOmega = static_cast<Real>(sSite4.y - sCenter.y) * fOmega;
    const gammaMatrix& gamma4 = __chiralGamma[GAMMA4];

    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 0);
    const SIndex& x_p_mu_Fermion = pFermionMove[linkIndex * 2];

    //If one of the sites is on surface, it has no contribution.
    //Note that, the bond on surface is equivelant to both sites on surface.

    if (x_p_mu_Fermion.IsDirichlet() || sSite.IsDirichlet())
    {
        return;
    }

    if (x_p_mu_Fermion.NeedToOpposite())
    {
        fYOmega = -fYOmega * fKai;
    }
    else
    {
        fYOmega = fYOmega * fKai;
    }

    //====================
    //Get things
    const deviceWilsonVectorSU3& x_Left = pInverseDDdagger[uiSiteIndex];
    deviceWilsonVectorSU3 x_Right(pInverseD[uiSiteIndex]);
    //all not on surface
    deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
    const deviceWilsonVectorSU3& x_p_mu_Left = pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex];

    const deviceSU3& x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

    //(1-gamma4)Y
    x_p_mu_Right.MulReal(fYOmega);
    x_p_mu_Right.Sub(gamma4.MulWilsonC(x_p_mu_Right));

    deviceSU3 mid = deviceSU3::makeSU3Contract(x_Left, x_p_mu_Right);

    //(1+gamma4)Y
    x_Right.MulReal(fYOmega);
    x_Right.Add(gamma4.MulWilsonC(x_Right));
    mid.Add(deviceSU3::makeSU3Contract(x_Right, x_p_mu_Left));

    deviceSU3 forceOfThisLink = x_Gauge_element.MulC(mid);
    forceOfThisLink.Ta();
    pForce[linkIndex].Add(forceOfThisLink);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDWilsonForceSU3_DR_Y(
    const deviceWilsonVectorSU3* __restrict__ pInverseD,
    const deviceWilsonVectorSU3* __restrict__ pInverseDDdagger,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    deviceSU3* pForce,
    SSmallInt4 sCenter,
    Real fKai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    Real fXOmega = static_cast<Real>(sSite4.x - sCenter.x) * fOmega;
    const gammaMatrix& gamma4 = __chiralGamma[GAMMA4];

    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 1);
    const SIndex& x_p_mu_Fermion = pFermionMove[linkIndex * 2];

    //If one of the sites is on surface, it has no contribution.
    //Note that, the bond on surface is equivelant to both sites on surface.

    if (x_p_mu_Fermion.IsDirichlet() || sSite.IsDirichlet())
    {
        return;
    }

    if (x_p_mu_Fermion.NeedToOpposite())
    {
        fXOmega = -fXOmega * fKai;
    }
    else
    {
        fXOmega = fXOmega * fKai;
    }

    //====================
    //Get things
    const deviceWilsonVectorSU3& x_Left = pInverseDDdagger[uiSiteIndex];
    deviceWilsonVectorSU3 x_Right(pInverseD[uiSiteIndex]);
    //all not on surface
    deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
    const deviceWilsonVectorSU3& x_p_mu_Left = pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex];

    const deviceSU3& x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

    //-(1-gamam4)X
    x_p_mu_Right.MulReal(-fXOmega);
    x_p_mu_Right.Sub(gamma4.MulWilsonC(x_p_mu_Right));

    deviceSU3 mid = deviceSU3::makeSU3Contract(x_Left, x_p_mu_Right);

    //-(1+gamma4)X
    x_Right.MulReal(-fXOmega);
    x_Right.Add(gamma4.MulWilsonC(x_Right));
    mid.Add(deviceSU3::makeSU3Contract(x_Right, x_p_mu_Left));

    deviceSU3 forceOfThisLink = x_Gauge_element.MulC(mid);
    forceOfThisLink.Ta();
    pForce[linkIndex].Add(forceOfThisLink);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDWilsonForceSU3_DR_T(
    const deviceWilsonVectorSU3* __restrict__ pInverseD,
    const deviceWilsonVectorSU3* __restrict__ pInverseDDdagger,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    deviceSU3* pForce,
    SSmallInt4 sCenter,
    Real fKai,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fOmega,
#else
    Real fOmega,
#endif
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    const gammaMatrix& gamma4 = __chiralGamma[GAMMA4];
    const gammaMatrix& sigma12 = __chiralGamma[SIGMA12E];

    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);
    const SIndex& x_p_mu_Fermion = pFermionMove[linkIndex * 2];

    //If one of the sites is on surface, it has no contribution.
    //Note that, the bond on surface is equivelant to both sites on surface.

    if (x_p_mu_Fermion.IsDirichlet() || sSite.IsDirichlet())
    {
        return;
    }

    //====================
    //Get things
    const deviceWilsonVectorSU3& x_Left = pInverseDDdagger[uiSiteIndex];
    deviceWilsonVectorSU3 x_Right(pInverseD[uiSiteIndex]);
    //all not on surface
    deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
    const deviceWilsonVectorSU3& x_p_mu_Left = pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex];

    const deviceSU3& x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

    const Real fFac = (x_p_mu_Fermion.NeedToOpposite() ? F(-1.0) : F(1.0)) * fKai;
    const Real fAng = F(0.5) * fOmega;
    const Real fCos = fFac * (_cos(fAng) - F(1.0));
    const Real fSin = fFac * _sin(fAng);
;
    //first term has same sign as _Y function (-1)
    //we already use fFac = -1 or 1
    //if (x_p_mu_Fermion.NeedToOpposite())
    //{
    //    fSin = -fSin;
    //    fCos = -fCos;
    //}

    x_p_mu_Right.Sub(gamma4.MulWilsonC(x_p_mu_Right));
    deviceWilsonVectorSU3 x_p_mu_Right_real = x_p_mu_Right.MulRealC(fCos);
    x_p_mu_Right.MulComp(_make_cuComplex(F(0.0), fSin));
    x_p_mu_Right = sigma12.MulWilsonC(x_p_mu_Right);
    x_p_mu_Right_real.Add(x_p_mu_Right);
    deviceSU3 mid = deviceSU3::makeSU3Contract(x_Left, x_p_mu_Right_real);

    x_Right.Add(gamma4.MulWilsonC(x_Right));
    deviceWilsonVectorSU3 x_Right_real = x_Right.MulRealC(fCos);
    x_Right.MulComp(_make_cuComplex(F(0.0), -fSin));
    x_Right = sigma12.MulWilsonC(x_Right);
    x_Right_real.Add(x_Right);
    mid.Add(deviceSU3::makeSU3Contract(x_Right_real, x_p_mu_Left));

    deviceSU3 forceOfThisLink = x_Gauge_element.MulC(mid);
    forceOfThisLink.Ta();
    pForce[linkIndex].Add(forceOfThisLink);
}

#pragma endregion

#pragma endregion

/**
* It will not produce correct result unless _CLG_USE_LAUNCH_BOUND is set to 1
* It seems max-registor count is reached without throw an error
* To avoid it, we use slower method, split it into two functions.
*/
void CFieldFermionWilsonSquareSU3DR::DOperator(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, 
    UBOOL bDagger, EOperatorCoefficientType eOCT, 
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{

#if _CLG_ROTATING_NEW_IMP

    deviceWilsonVectorSU3* pTarget = (deviceWilsonVectorSU3*)pTargetBuffer;
    const deviceWilsonVectorSU3* pSource = (deviceWilsonVectorSU3*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
        
    preparethread;

    if (m_bExponential)
    {
        _kernelDFermionWilsonSquareSU3_DR_Exponential_0 << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pTarget,
            m_fKai,
            CCommonData::m_fOmega,
            CCommonData::m_sCenter,
            m_byFieldId,
            bDagger,
            m_bNaive,
            eOCT,
            fRealCoeff,
            cCmpCoeff);

        _kernelDFermionWilsonSquareSU3_DR_Exponential_X << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pTarget,
            m_fKai,
            CCommonData::m_fOmega,
            CCommonData::m_sCenter,
            m_byFieldId,
            bDagger,
            m_bNaive,
            eOCT,
            fRealCoeff,
            cCmpCoeff);

        _kernelDFermionWilsonSquareSU3_DR_Exponential_Y << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pTarget,
            m_fKai,
            CCommonData::m_fOmega,
            CCommonData::m_sCenter,
            m_byFieldId,
            bDagger,
            m_bNaive,
            eOCT,
            fRealCoeff,
            cCmpCoeff);

        _kernelDFermionWilsonSquareSU3_DR_Exponential_Z << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pTarget,
            m_fKai,
            CCommonData::m_fOmega,
            CCommonData::m_sCenter,
            m_byFieldId,
            bDagger,
            m_bNaive,
            eOCT,
            fRealCoeff,
            cCmpCoeff);

        _kernelDFermionWilsonSquareSU3_DR_Exponential_T << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pTarget,
            m_fKai,
            CCommonData::m_fOmega,
            CCommonData::m_sCenter,
            m_byFieldId,
            bDagger,
            m_bNaive,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }
    else
    {
        CFieldFermionWilsonSquareSU3D::DOperator(pTargetBuffer, pBuffer, pGaugeBuffer, bDagger,
            eOCT, fRealCoeff, cCmpCoeff);

        _kernelDFermionWilsonSquareSU3_DR << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pTarget,
            m_fKai,
            CCommonData::m_fOmega,
            CCommonData::m_sCenter,
            m_byFieldId,
            bDagger,
            m_bNaive,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }

#else

    CFieldFermionWilsonSquareSU3D::DOperator(pTargetBuffer, pBuffer, pGaugeBuffer, bDagger,
        eOCT, fRealCoeff, cCmpCoeff);

    deviceWilsonVectorSU3* pTarget = (deviceWilsonVectorSU3*)pTargetBuffer;
    const deviceWilsonVectorSU3* pSource = (deviceWilsonVectorSU3*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    preparethread;

    if (m_bExponential)
    {
        _kernelDFermionWilsonSquareSU3_DR_Exponential_0 << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pTarget,
            m_fKai,
            CCommonData::m_fOmega,
            CCommonData::m_sCenter,
            m_byFieldId,
            bDagger,
            m_bNaive,
            eOCT,
            fRealCoeff,
            cCmpCoeff);

        _kernelDFermionWilsonSquareSU3_DR_Exponential_X << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pTarget,
            m_fKai,
            CCommonData::m_fOmega,
            CCommonData::m_sCenter,
            m_byFieldId,
            bDagger,
            m_bNaive,
            eOCT,
            fRealCoeff,
            cCmpCoeff);

        _kernelDFermionWilsonSquareSU3_DR_Exponential_Y << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pTarget,
            m_fKai,
            CCommonData::m_fOmega,
            CCommonData::m_sCenter,
            m_byFieldId,
            bDagger,
            m_bNaive,
            eOCT,
            fRealCoeff,
            cCmpCoeff);

        _kernelDFermionWilsonSquareSU3_DR_Exponential_T << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pTarget,
            m_fKai,
            CCommonData::m_fOmega,
            CCommonData::m_sCenter,
            m_byFieldId,
            bDagger,
            m_bNaive,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }
    else
    {
        _kernelDFermionWilsonSquareSU3_DR << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pTarget,
            m_fKai,
            CCommonData::m_fOmega,
            CCommonData::m_sCenter,
            m_byFieldId,
            bDagger,
            m_bNaive,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }

#endif
}

void CFieldFermionWilsonSquareSU3DR::DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const
{
    CFieldFermionWilsonSquareSU3D::DerivateDOperator(pForce, pDphi, pDDphi, pGaugeBuffer);

    deviceSU3* pForceSU3 = (deviceSU3*)pForce;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
    const deviceWilsonVectorSU3* pDphiBuffer = (deviceWilsonVectorSU3*)pDphi;
    const deviceWilsonVectorSU3* pDDphiBuffer = (deviceWilsonVectorSU3*)pDDphi;

    preparethread;

    if (m_bNaive)
    {
        _kernelDWilsonForceSU3_DR_X_Naive << <block, threads >> > (
            pDphiBuffer,
            pDDphiBuffer,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pForceSU3,
            CCommonData::m_sCenter, 
            m_fKai, 
            CCommonData::m_fOmega, 
            m_byFieldId);

        _kernelDWilsonForceSU3_DR_Y_Naive << <block, threads >> > (
            pDphiBuffer,
            pDDphiBuffer,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pForceSU3,
            CCommonData::m_sCenter, 
            m_fKai, 
            CCommonData::m_fOmega, 
            m_byFieldId);
    }
    else
    {
        _kernelDWilsonForceSU3_DR_X << <block, threads >> > (
            pDphiBuffer,
            pDDphiBuffer,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pForceSU3,
            CCommonData::m_sCenter, 
            m_fKai, 
            CCommonData::m_fOmega, 
            m_byFieldId);

        _kernelDWilsonForceSU3_DR_Y << <block, threads >> > (
            pDphiBuffer,
            pDDphiBuffer,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pForceSU3,
            CCommonData::m_sCenter, 
            m_fKai, 
            CCommonData::m_fOmega, 
            m_byFieldId);
    }

    if (m_bExponential)
    {
        _kernelDWilsonForceSU3_DR_T << <block, threads >> > (
            pDphiBuffer,
            pDDphiBuffer,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pForceSU3,
            CCommonData::m_sCenter, 
            m_fKai, 
            CCommonData::m_fOmega, 
            m_byFieldId);
    }
}

#pragma endregion

void CFieldFermionWilsonSquareSU3DR::InitialOtherParameters(CParameters& params)
{
    CFieldFermionWilsonSquareSU3::InitialOtherParameters(params);

    INT iNaive = 1;
    params.FetchValueINT(_T("Naive"), iNaive);
    m_bNaive = 0 != iNaive;

    INT iExponential = 0;
    params.FetchValueINT(_T("Exponential"), iExponential);
    m_bExponential = 0 != iExponential;
}

void CFieldFermionWilsonSquareSU3DR::CopyTo(CField* U) const
{
    CFieldFermionWilsonSquareSU3D::CopyTo(U);
    CFieldFermionWilsonSquareSU3DR * pField = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(U);
    if (NULL != pField)
    {
        pField->m_bNaive = m_bNaive;
        pField->m_bExponential = m_bExponential;
    }
}

CCString CFieldFermionWilsonSquareSU3DR::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CFieldFermionWilsonSquareSU3DR\n");
    sRet = sRet + tab + _T("Hopping : ") + appFloatToString(CCommonData::m_fKai) + _T("\n");
    sRet = sRet + tab + _T("Naive : ") + (m_bNaive ? _T("1") : _T("0"));
    sRet = sRet + tab + _T("Exponential : ") + (m_bExponential ? _T("1") : _T("0"));
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================