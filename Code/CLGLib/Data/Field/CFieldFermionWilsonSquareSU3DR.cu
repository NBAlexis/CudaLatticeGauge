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

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR(
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
    Real fXkOmega = static_cast<Real>(sSite4.x - sCenter.x) * fOmega * kai;

    deviceWilsonVectorSU3 right_element = bDDagger ? gamma5.MulWilsonC(pDeviceData[uiSiteIndex]) : pDeviceData[uiSiteIndex];

    deviceWilsonVectorSU3 term3(__chiralGamma[SIGMA12].MulWilsonC(right_element));
    term3 = gamma4.MulWilsonC(term3);
    term3.MulComp(_make_cuComplex(F(0.0), fOmega * kai));

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
    for (UINT idir = 0; idir < 2; ++idir)
    {
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

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DR_Exponential(
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
    gammaMatrix sigma12 = __chiralGamma[SIGMA12];
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

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //idir = mu, we have x, y and t term but no z term
    for (UINT idir = 0; idir < 3; ++idir)
    {
        //=========================
        //get things
        if (2 == idir)
        {
            idir = 3;
        }
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
            else if (3 == idir)
            {
                //the sign is opppsite as "y Omega" term same as "x Omega", which is 
                //0.5 i k omega sigma12 (1 - gamma4)
                u_phi_x_p_m.MulComp(make_cuComplex(F(0.0), fkOmega * F(0.5)));
                u_phi_x_p_m.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
                result.Add(sigma12.MulWilsonC(u_phi_x_p_m));
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
            else if (3 == idir)
            {
                //the sign is opppsite as "y Omega" term same as "x Omega", which is 
                //-0.5 i k omega sigma12 (1 - gamma4)
                u_phi_x_p_m.MulComp(make_cuComplex(F(0.0), -fkOmega * F(0.5)));
                u_phi_x_p_m.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
                result.Add(sigma12.MulWilsonC(u_phi_x_p_m));
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
            else if (3 == idir)
            {
                //the sign is same as "y Omega" term opposite as "x Omega", which is 
                //-0.5 i k omega sigma12 (1 + gamma4)
                u_dagger_phi_x_m_m.MulComp(make_cuComplex(F(0.0), -fkOmega * F(0.5)));
                u_dagger_phi_x_m_m.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
                result.Add(sigma12.MulWilsonC(u_dagger_phi_x_m_m));
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
            else if (3 == idir)
            {
                //the sign is same as "y Omega" term opposite as "x Omega", which is 
                //0.5 i k omega sigma12 (1 + gamma4)
                u_dagger_phi_x_m_m.MulComp(make_cuComplex(F(0.0), fkOmega * F(0.5)));
                u_dagger_phi_x_m_m.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
                result.Add(sigma12.MulWilsonC(u_dagger_phi_x_m_m));
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
    Real fOmega,
    BYTE byFieldId)
{
    intokernalInt4;
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    Real fYOmega = static_cast<Real>(sSite4.y - sCenter.y) * fOmega;
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 0);
    SIndex x_p_mu_Fermion = pFermionMove[linkIndex * 2];

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
    deviceWilsonVectorSU3 x_Left(pInverseDDdagger[uiSiteIndex]);
    deviceWilsonVectorSU3 x_Right(pInverseD[uiSiteIndex]);
    //all not on surface
    deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
    deviceWilsonVectorSU3 x_p_mu_Left(pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex]);

    deviceSU3 x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

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
    Real fOmega,
    BYTE byFieldId)
{
    intokernalInt4;
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    Real fXOmega = static_cast<Real>(sSite4.x - sCenter.x) * fOmega;
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 1);
    SIndex x_p_mu_Fermion = pFermionMove[linkIndex * 2];

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
    deviceWilsonVectorSU3 x_Left(pInverseDDdagger[uiSiteIndex]);
    deviceWilsonVectorSU3 x_Right(pInverseD[uiSiteIndex]);
    //all not on surface
    deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
    deviceWilsonVectorSU3 x_p_mu_Left(pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex]);

    deviceSU3 x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

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
    Real fOmega,
    BYTE byFieldId)
{
    intokernalInt4;
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    Real fYOmega = static_cast<Real>(sSite4.y - sCenter.y) * fOmega;
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 0);
    SIndex x_p_mu_Fermion = pFermionMove[linkIndex * 2];

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
    deviceWilsonVectorSU3 x_Left(pInverseDDdagger[uiSiteIndex]);
    deviceWilsonVectorSU3 x_Right(pInverseD[uiSiteIndex]);
    //all not on surface
    deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
    deviceWilsonVectorSU3 x_p_mu_Left(pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex]);

    deviceSU3 x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

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
    Real fOmega,
    BYTE byFieldId)
{
    intokernalInt4;
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    Real fXOmega = static_cast<Real>(sSite4.x - sCenter.x) * fOmega;
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 1);
    SIndex x_p_mu_Fermion = pFermionMove[linkIndex * 2];

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
    deviceWilsonVectorSU3 x_Left(pInverseDDdagger[uiSiteIndex]);
    deviceWilsonVectorSU3 x_Right(pInverseD[uiSiteIndex]);
    //all not on surface
    deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
    deviceWilsonVectorSU3 x_p_mu_Left(pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex]);

    deviceSU3 x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

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
    Real fOmega,
    BYTE byFieldId)
{
    intokernalInt4;
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    CLGComplex fkOmegaP = _make_cuComplex(F(0.0), F(0.5) * fKai * fOmega);
    CLGComplex fkOmegaM = _make_cuComplex(F(0.0), -F(0.5) * fKai * fOmega);
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix sigma12 = __chiralGamma[SIGMA12];

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);
    SIndex x_p_mu_Fermion = pFermionMove[linkIndex * 2];

    //If one of the sites is on surface, it has no contribution.
    //Note that, the bond on surface is equivelant to both sites on surface.

    if (x_p_mu_Fermion.IsDirichlet() || sSite.IsDirichlet())
    {
        return;
    }

    //====================
    //Get things
    deviceWilsonVectorSU3 x_Left(pInverseDDdagger[uiSiteIndex]);
    deviceWilsonVectorSU3 x_Right(pInverseD[uiSiteIndex]);
    //all not on surface
    deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
    deviceWilsonVectorSU3 x_p_mu_Left(pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex]);

    deviceSU3 x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

    //first term has same sign as _Y function (-1)
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        x_p_mu_Right.MulComp(fkOmegaM);
    }
    else
    {
        x_p_mu_Right.MulComp(fkOmegaP);
    }
    x_p_mu_Right.Sub(gamma4.MulWilsonC(x_p_mu_Right));
    x_p_mu_Right = sigma12.MulWilsonC(x_p_mu_Right);

    deviceSU3 mid = deviceSU3::makeSU3Contract(x_Left, x_p_mu_Right);

    //second term has same sign as _X function (+1)
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        x_Right.MulComp(fkOmegaP);
    }
    else
    {
        x_Right.MulComp(fkOmegaM);
    }
    x_Right.Add(gamma4.MulWilsonC(x_Right));
    x_Right = sigma12.MulWilsonC(x_Right);
    mid.Add(deviceSU3::makeSU3Contract(x_Right, x_p_mu_Left));

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
    CFieldFermionWilsonSquareSU3D::DOperator(pTargetBuffer, pBuffer, pGaugeBuffer, bDagger,
        eOCT, fRealCoeff, cCmpCoeff);

    deviceWilsonVectorSU3* pTarget = (deviceWilsonVectorSU3*)pTargetBuffer;
    const deviceWilsonVectorSU3* pSource = (deviceWilsonVectorSU3*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
        
    preparethread;

    if (m_bExponential)
    {
        _kernelDFermionWilsonSquareSU3_DR_Exponential << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pTarget,
            m_fKai,
            m_fOmega,
            m_sCenter,
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
            m_fOmega,
            m_sCenter,
            m_byFieldId,
            bDagger,
            m_bNaive,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }

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
            m_sCenter, m_fKai, m_fOmega, m_byFieldId);

        _kernelDWilsonForceSU3_DR_Y_Naive << <block, threads >> > (
            pDphiBuffer,
            pDDphiBuffer,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pForceSU3,
            m_sCenter, m_fKai, m_fOmega, m_byFieldId);
    }
    else
    {
        _kernelDWilsonForceSU3_DR_X << <block, threads >> > (
            pDphiBuffer,
            pDDphiBuffer,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pForceSU3,
            m_sCenter, m_fKai, m_fOmega, m_byFieldId);

        _kernelDWilsonForceSU3_DR_Y << <block, threads >> > (
            pDphiBuffer,
            pDDphiBuffer,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pForceSU3,
            m_sCenter, m_fKai, m_fOmega, m_byFieldId);
    }

    if (m_bExponential)
    {
        _kernelDWilsonForceSU3_DR_T << <block, threads >> > (
            pDphiBuffer,
            pDDphiBuffer,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            pForceSU3,
            m_sCenter, m_fKai, m_fOmega, m_byFieldId);
    }
}

#pragma endregion

void CFieldFermionWilsonSquareSU3DR::InitialOtherParameters(CParameters& params)
{
    CFieldFermionWilsonSquareSU3::InitialOtherParameters(params);

    Real fOmega = 0.1f;
    params.FetchValueReal(_T("Omega"), fOmega);
    m_fOmega = fOmega;

    TArray<INT> centerArray;
    params.FetchValueArrayINT(_T("Center"), centerArray);
    if (centerArray.Num() > 3)
    {
        m_sCenter.x = static_cast<SBYTE>(centerArray[0]);
        m_sCenter.y = static_cast<SBYTE>(centerArray[1]);
        m_sCenter.z = static_cast<SBYTE>(centerArray[2]);
        m_sCenter.w = static_cast<SBYTE>(centerArray[3]);
    }

    INT iNaive = 1;
    params.FetchValueINT(_T("Naive"), iNaive);
    if (0 == iNaive)
    {
        m_bNaive = FALSE;
    }

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
        pField->m_fOmega = m_fOmega;
        pField->m_sCenter = m_sCenter;
        pField->m_bNaive = m_bNaive;
        pField->m_bExponential = m_bExponential;
    }
}

CCString CFieldFermionWilsonSquareSU3DR::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CFieldFermionWilsonSquareSU3DR\n");
    sRet = sRet + tab + _T("Hopping : ") + appFloatToString(CCommonData::m_fKai) + _T("\n");
    sRet = sRet + tab + _T("Omega : ") + appFloatToString(m_fOmega) + _T("\n");
    CCString sCenter;
    sCenter.Format(_T("Center: [%d, %d, %d, %d]\n"), m_sCenter.x, m_sCenter.y, m_sCenter.z, m_sCenter.w);
    sRet = sRet + tab + sCenter;
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================