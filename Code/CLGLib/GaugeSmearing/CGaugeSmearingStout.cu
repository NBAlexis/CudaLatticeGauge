//=============================================================================
// FILENAME : CGaugeSmearingStout.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [mm/dd/yy]
//  [07/25/2025 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CGaugeSmearingStout.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CGaugeSmearingStoutSU3)

#pragma region kernels

/**
* assert 6 == plaqCount
* assert 4 == plaqLength
*
* for dir = mu,
* for nu from 0 to 3 and skip mu
*
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStoutFat3(
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pRes,
    Real rhojk, Real rho4mu,
    const SIndex* __restrict__ pCachedStapleIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    BYTE byFieldId)
{
    intokernalDir;

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1; //3
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
#endif

    deviceGauge& resC = pRes[uiLinkIndex];
    resC = _makeZero<deviceGauge>();
    #pragma unroll
    for (BYTE elidx = 0U; elidx < 6U; ++elidx)
    {
        //for spatial link, it is: forward-backward, xy, xz, xt; yx, yz, yt; zx, zy, xt; so for the first 4 staples, use rhojk, for the last 2, use rho4mu
        //for temporal link, it is always rho4mu 
        const Real rho = (3 == dir || elidx >= 4U) ? rho4mu : rhojk;
        const UINT byIdxStart = plaqLengthm1 * elidx + uiLinkIndex * plaqCountAllLink;
        const SIndex& first = pCachedStapleIndex[byIdxStart];
        deviceGauge res(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

        if (first.NeedToDagger())
        {
            _dagger(res);
        }
        #pragma unroll
        for (BYTE j = 1U; j < plaqLengthm1; ++j)
        {
            const SIndex& nextlink = pCachedStapleIndex[byIdxStart + j];
            const deviceGauge& toMul = _deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink);

            if (nextlink.NeedToDagger())
            {
                _muldag(res, toMul);
            }
            else
            {
                _mul(res, toMul);
            }
        }

        _mul(res, rho);
        _add(resC, res);
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStoutExpStep(
    const deviceGauge* __restrict__ pOrignalGauge,
    const deviceGauge* __restrict__ pC,
    deviceGauge* pQ,
    deviceGauge* pQ2,
    deviceGauge* pExpIQ,
    deviceGauge* pResGauge,
    cuDoubleComplex* f012_exp2u_expm1u,
    DOUBLE* u_w2_cosw_xi0_denorm,
    BYTE* expCondition)
{
    intokernalDir_NoDir;

    deviceGauge& Q = pQ[uiLinkIndex];
    deviceGauge& Q2 = pQ2[uiLinkIndex];
    deviceGauge& expIQ = pExpIQ[uiLinkIndex];

    const UINT dparamstart = uiLinkIndex * 5U;
    DOUBLE& u = u_w2_cosw_xi0_denorm[dparamstart];
    DOUBLE& w2 = u_w2_cosw_xi0_denorm[dparamstart + 1];
    DOUBLE& cosw = u_w2_cosw_xi0_denorm[dparamstart + 2];
    DOUBLE& xi0 = u_w2_cosw_xi0_denorm[dparamstart + 3];
    DOUBLE& dnorm = u_w2_cosw_xi0_denorm[dparamstart + 4];
    cuDoubleComplex& f0 = f012_exp2u_expm1u[dparamstart];
    cuDoubleComplex& f1 = f012_exp2u_expm1u[dparamstart + 1];
    cuDoubleComplex& f2 = f012_exp2u_expm1u[dparamstart + 2];
    cuDoubleComplex& exp2u = f012_exp2u_expm1u[dparamstart + 3];
    cuDoubleComplex& expm1u = f012_exp2u_expm1u[dparamstart + 4];

    Q = pC[uiLinkIndex];
    _muldag(Q, pOrignalGauge[uiLinkIndex]);
    _ta(Q);
    _mul(Q, _make_cuComplex(F(0.0), -F(1.0)));

    //if (0 == uiLinkIndex)
    //{
    //    _print(Q, "Q");
    //}

    Q2 = Q;
    _mul(Q2, Q);

    //U = exp(iQ) U
#if _CLG_DOUBLEFLOAT
    DOUBLE c0 = _detv(Q).x;
#else
    DOUBLE c0 = static_cast<DOUBLE>(_detv(Q).x);
#endif
    u = _tr(Q2).x * 0.166666666666666666667;
    if (u < _CLG_FLT_MIN_)
    {
        //u is small means Q is small, exp(Q) = 1
        //pEffectiveGauge[uiLinkIndex] = pDeviceData[uiLinkIndex];
        //printf("ever here?\n");
        expCondition[uiLinkIndex] = 2;
        return;
    }

    u = sqrt(u);

    expCondition[uiLinkIndex] = 0;
    if (c0 < 0)
    {
        c0 = -c0;
        expCondition[uiLinkIndex] = 1;
    }

    //DOUBLE th = c0 / c0max;
    //This is the last time to use c0, and the only time to use c0max
    //In the following, c0 is discarded, and store th
    c0 = c0 / (2.0 * u * u * u);

    if (c0 >= 1.0)
    {
        c0 = 0.0;
    }
    else
    {
        c0 = acos(c0) * 0.33333333333333333333;
    }
    //there is no builtin acos

    const DOUBLE w = sin(c0) * u * 1.73205080756887729352744634151;
    u = u * cos(c0);
    //This is the last time to use th, so in the following, th is discarded, and store xi

    const DOUBLE u2 = u * u;
    w2 = w * w;
    if (abs(w2) < _CLG_FLT_MIN)
    {
        expCondition[uiLinkIndex] = 2;
        return;
    }
    xi0 = (abs(w) > 0.05) ? (sin(w) / w) : (
        1.0 - w2 * (
            1.0 - w2 * (
                1.0 - w2 * 0.0238095238095238095238095238095
                ) * 0.05
            ) * 0.16666666666666666666667
        );
    const DOUBLE f2u = 2.0 * u;
    const DOUBLE f8u2 = 8.0 * u2;
    const DOUBLE fu2_m_w2 = u2 - w2;
    cosw = cos(w);
    const DOUBLE f3u2 = 3.0 * u2;

    exp2u = make_cuDoubleComplex(cos(f2u), sin(f2u));
    expm1u = make_cuDoubleComplex(cos(u), -sin(u));

    if (abs(f8u2 + fu2_m_w2) < _CLG_FLT_MIN_)
    {
        expCondition[uiLinkIndex] = 2;
        return;
    }
    dnorm = 1.0 / (f8u2 + fu2_m_w2);


    expIQ = Q;
    f1 = cuCsub(cuCmulf_cd(exp2u, f2u),
        cuCmul(expm1u,
            make_cuDoubleComplex(f2u * cosw, -(f3u2 - w2) * xi0)
        ));
    f1 = cuCmulf_cd(f1, dnorm);
    //f1 = expCondition[uiLinkIndex] ? cuConj(cuCmulf_cd(f1, -dnorm)) : cuCmulf_cd(f1, dnorm);
    _mul(expIQ, _cToRealC(expCondition[uiLinkIndex] ? cuConj(cuCmulf_cd(f1, -1.0)) : f1));

    f2 = cuCsub(exp2u,
        cuCmul(expm1u,
            make_cuDoubleComplex(cosw, 3.0 * u * xi0)
        ));
    f2 = cuCmulf_cd(f2, dnorm);
    //f2 = expCondition[uiLinkIndex] ? cuConj(cuCmulf_cd(f2, dnorm)) : cuCmulf_cd(f2, dnorm);
    _add(expIQ, _mulC(Q2, _cToRealC(expCondition[uiLinkIndex] ? cuConj(f2) : f2)));

    f0 = cuCadd(cuCmulf_cd(exp2u, fu2_m_w2),
        cuCmul(expm1u,
            make_cuDoubleComplex(f8u2 * cosw, f2u * (f3u2 + w2) * xi0)
        ));
    //f0 = expCondition[uiLinkIndex] ? cuConj(cuCmulf_cd(f0, dnorm)) : cuCmulf_cd(f0, dnorm);
    f0 = cuCmulf_cd(f0, dnorm);
    _add(expIQ, _cToRealC(expCondition[uiLinkIndex] ? cuConj(f0) : f0));

    //now pEffectiveGauge is exp(iQ)
    pResGauge[uiLinkIndex] = expIQ;
    _mul(pResGauge[uiLinkIndex], pOrignalGauge[uiLinkIndex]);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStoutExpStep_NoQ(
    const deviceGauge* __restrict__ pOrignalGauge,
    const deviceGauge* __restrict__ pC,
    deviceGauge* pExpIQ,
    deviceGauge* pResGauge,
    cuDoubleComplex* f012_exp2u_expm1u,
    DOUBLE* u_w2_cosw_xi0_denorm,
    BYTE* expCondition)
{
    intokernalDir_NoDir;

    deviceGauge& expIQ = pExpIQ[uiLinkIndex];

    const UINT dparamstart = uiLinkIndex * 5U;
    DOUBLE& u = u_w2_cosw_xi0_denorm[dparamstart];
    DOUBLE& w2 = u_w2_cosw_xi0_denorm[dparamstart + 1];
    DOUBLE& cosw = u_w2_cosw_xi0_denorm[dparamstart + 2];
    DOUBLE& xi0 = u_w2_cosw_xi0_denorm[dparamstart + 3];
    DOUBLE& dnorm = u_w2_cosw_xi0_denorm[dparamstart + 4];
    cuDoubleComplex& f0 = f012_exp2u_expm1u[dparamstart];
    cuDoubleComplex& f1 = f012_exp2u_expm1u[dparamstart + 1];
    cuDoubleComplex& f2 = f012_exp2u_expm1u[dparamstart + 2];
    cuDoubleComplex& exp2u = f012_exp2u_expm1u[dparamstart + 3];
    cuDoubleComplex& expm1u = f012_exp2u_expm1u[dparamstart + 4];

    deviceGauge Q = pC[uiLinkIndex];
    _muldag(Q, pOrignalGauge[uiLinkIndex]);
    _ta(Q);
    _mul(Q, _make_cuComplex(F(0.0), -F(1.0)));

    deviceGauge Q2 = Q;
    _mul(Q2, Q);

    //U = exp(iQ) U
#if _CLG_DOUBLEFLOAT
    DOUBLE c0 = _detv(Q).x;
#else
    DOUBLE c0 = static_cast<DOUBLE>(_detv(Q).x);
#endif
    u = _tr(Q2).x * 0.166666666666666666667;
    if (u < _CLG_FLT_MIN_)
    {
        //u is small means Q is small, exp(Q) = 1
        //pEffectiveGauge[uiLinkIndex] = pDeviceData[uiLinkIndex];
        //printf("ever here?\n");
        expCondition[uiLinkIndex] = 2;
        return;
    }

    u = sqrt(u);

    expCondition[uiLinkIndex] = 0;
    if (c0 < 0)
    {
        c0 = -c0;
        expCondition[uiLinkIndex] = 1;
    }

    //DOUBLE th = c0 / c0max;
    //This is the last time to use c0, and the only time to use c0max
    //In the following, c0 is discarded, and store th
    c0 = c0 / (2.0 * u * u * u);

    if (c0 >= 1.0)
    {
        c0 = 0.0;
    }
    else
    {
        c0 = acos(c0) * 0.33333333333333333333;
    }
    //there is no builtin acos

    const DOUBLE w = sin(c0) * u * 1.73205080756887729352744634151;
    u = u * cos(c0);
    //This is the last time to use th, so in the following, th is discarded, and store xi

    const DOUBLE u2 = u * u;
    w2 = w * w;
    if (abs(w2) < _CLG_FLT_MIN)
    {
        expCondition[uiLinkIndex] = 2;
        return;
    }
    xi0 = (abs(w) > 0.05) ? (sin(w) / w) : (
        1.0 - w2 * (
            1.0 - w2 * (
                1.0 - w2 * 0.0238095238095238095238095238095
                ) * 0.05
            ) * 0.16666666666666666666667
        );
    const DOUBLE f2u = 2.0 * u;
    const DOUBLE f8u2 = 8.0 * u2;
    const DOUBLE fu2_m_w2 = u2 - w2;
    cosw = cos(w);
    const DOUBLE f3u2 = 3.0 * u2;

    exp2u = make_cuDoubleComplex(cos(f2u), sin(f2u));
    expm1u = make_cuDoubleComplex(cos(u), -sin(u));

    if (abs(f8u2 + fu2_m_w2) < _CLG_FLT_MIN)
    {
        expCondition[uiLinkIndex] = 2;
        return;
    }
    dnorm = 1.0 / (f8u2 + fu2_m_w2);

    expIQ = Q;
    f1 = cuCsub(cuCmulf_cd(exp2u, f2u),
        cuCmul(expm1u,
            make_cuDoubleComplex(f2u * cosw, -(f3u2 - w2) * xi0)
        ));
    f1 = cuCmulf_cd(f1, dnorm);
    //f1 = expCondition[uiLinkIndex] ? cuConj(cuCmulf_cd(f1, -dnorm)) : cuCmulf_cd(f1, dnorm);
    _mul(expIQ, _cToRealC(expCondition[uiLinkIndex] ? cuConj(cuCmulf_cd(f1, -1.0)) : f1));

    f2 = cuCsub(exp2u,
        cuCmul(expm1u,
            make_cuDoubleComplex(cosw, 3.0 * u * xi0)
        ));
    f2 = cuCmulf_cd(f2, dnorm);
    //f2 = expCondition[uiLinkIndex] ? cuConj(cuCmulf_cd(f2, dnorm)) : cuCmulf_cd(f2, dnorm);
    _add(expIQ, _mulC(Q2, _cToRealC(expCondition[uiLinkIndex] ? cuConj(f2) : f2)));

    f0 = cuCadd(cuCmulf_cd(exp2u, fu2_m_w2),
        cuCmul(expm1u,
            make_cuDoubleComplex(f8u2 * cosw, f2u * (f3u2 + w2) * xi0)
        ));
    //f0 = expCondition[uiLinkIndex] ? cuConj(cuCmulf_cd(f0, dnorm)) : cuCmulf_cd(f0, dnorm);
    f0 = cuCmulf_cd(f0, dnorm);
    _add(expIQ, _cToRealC(expCondition[uiLinkIndex] ? cuConj(f0) : f0));

    //now pEffectiveGauge is exp(iQ)
    pResGauge[uiLinkIndex] = expIQ;
    _mul(pResGauge[uiLinkIndex], pOrignalGauge[uiLinkIndex]);
}

/**
* calculate force also using staple, but replace one of the links with f0
* (60) to (65) with (57) and (58) of hep-lat/0311018
* with (69) (70) (73) (74)
* 
* Note that, in CLGLib, we calculate f^{dagger}, so we only need Lambda^{dagger}
* 
* This kernel has too much registers, split it
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStoutForceLambda(
    const deviceGauge* __restrict__ pOriginalGauge,
    const deviceGauge* __restrict__ pf0,
    const deviceGauge* __restrict__ pQ,
    const deviceGauge* __restrict__ pQ2,
    deviceGauge* pLambda,
    const cuDoubleComplex* __restrict__ f012_exp2u_expm1u,
    const DOUBLE* __restrict__ u_w2_cosw_xi0_denorm,
    const BYTE* __restrict__ expCondition)
{
    intokernalDir_NoDir;
    const BYTE& minus = expCondition[uiLinkIndex];

    if (2 == minus)
    {
        //when Q is zero, Lambda is also zero
        pLambda[uiLinkIndex] = _makeZero<deviceGauge>();
        return;
    }

    deviceGauge& res = pLambda[uiLinkIndex];
    const deviceGauge& Q = pQ[uiLinkIndex];
    const deviceGauge& Q2 = pQ2[uiLinkIndex];
    const deviceGauge& U = pOriginalGauge[uiLinkIndex];
    const deviceGauge& Sigma = pf0[uiLinkIndex];

    const UINT dparamstart = uiLinkIndex * 5U;
    const DOUBLE& u = u_w2_cosw_xi0_denorm[dparamstart];
    const DOUBLE& w2 = u_w2_cosw_xi0_denorm[dparamstart + 1];
    const DOUBLE& cosw = u_w2_cosw_xi0_denorm[dparamstart + 2];
    const DOUBLE& xi0 = u_w2_cosw_xi0_denorm[dparamstart + 3];
    const DOUBLE& dnorm = u_w2_cosw_xi0_denorm[dparamstart + 4];
    const cuDoubleComplex& f0 = f012_exp2u_expm1u[dparamstart];
    const cuDoubleComplex& f1 = f012_exp2u_expm1u[dparamstart + 1];
    const cuDoubleComplex& f2 = f012_exp2u_expm1u[dparamstart + 2];
    const cuDoubleComplex& exp2u = f012_exp2u_expm1u[dparamstart + 3];
    const cuDoubleComplex& expm1u = f012_exp2u_expm1u[dparamstart + 4];

    deviceGauge B = U;
    _muldag(B, Sigma); //B = U.Sigma
    res = B;
    _mul(res, Q); //res = U.Sigma.Q
    _add(res, _mulC(Q, B)); // res = U.Sigma.Q + Q.U.Sigma
    _mul(res, _cToRealC(expCondition[uiLinkIndex] ? cuConj(f2) : f2)); // res = f2 U.Sigma.Q + f2 Q.U.Sigma
    _mul(B, _cToRealC(expCondition[uiLinkIndex] ? cuConj(cuCmulf_cd(f1, -1.0)) : f1)); //B = f1 U.Sigma
    _add(res, B); //res = f1 U.Sigma + f2 U.Sigma.Q + f2 Q.U.Sigma

    //calculate six b
    const DOUBLE u2 = u * u;
    const DOUBLE xi1 = (abs(w2) < _CLG_FLT_MIN) ? (cosw - xi0) : ((cosw - xi0) / w2);
    const DOUBLE dnormb = dnorm * dnorm * 0.5;
    const DOUBLE f2u = 2.0 * u;
    const DOUBLE f3u = 3.0 * u;
    const DOUBLE f24u = 24.0 * u;
    const DOUBLE f3u2 = 3.0 * u2;
    const DOUBLE f3u2x1 = f3u2 * xi1;
    const DOUBLE f3xi0 = 3.0 * xi0;
    const DOUBLE cosw_p_xi0 = cosw + xi0;
    const DOUBLE f_3u2_m_w2 = f3u2 - w2;
    const DOUBLE f_30u2_p_2w2 = 2.0 * (15.0 * u2 + w2);

    //j = 0, r_0^(1) and r_0^(2)
    // (60)
    cuDoubleComplex rj1 = cuCadd(
        cuCmul(exp2u, make_cuDoubleComplex(f2u, 2.0 * (u2 - w2))),
        cuCmul(expm1u, make_cuDoubleComplex(16.0 * u * cosw + f2u * (f3u2 + w2) * xi0, 2.0 * (9.0 * u2 + w2) * xi0 - 8.0 * u2 * cosw))
    );
    // (63)
    cuDoubleComplex rj2 = cuCadd(
        cuCmulf_cd(exp2u, -2.0),
        cuCmul(expm1u, make_cuDoubleComplex(-8.0 * u2 * xi0, f2u * (cosw_p_xi0 + f3u2x1) ))
        );

    const cuDoubleComplex b10 = cuCmulf_cd(
        cuCsub(cuCadd(cuCmulf_cd(rj1, f2u), cuCmulf_cd(rj2, f_3u2_m_w2)),
            cuCmulf_cd(f0, f_30u2_p_2w2)),
        dnormb
    );
    const cuDoubleComplex b20 = cuCmulf_cd(
        cuCsub(cuCsub(rj1, cuCmulf_cd(rj2, f3u)),
            cuCmulf_cd(f0, f24u)),
        minus ? -dnormb : dnormb
    );

    //j = 1, r_1^(1) and r_1^(2)
    // (61)
    rj1 = cuCadd(
        cuCmul(exp2u, make_cuDoubleComplex(2.0, 4.0 * u)),
        cuCmul(expm1u, make_cuDoubleComplex(f_3u2_m_w2 * xi0 - 2.0 * cosw, f2u * (cosw + f3xi0)))
    );
    // (64)
    rj2 = cuCmul(expm1u, 
        make_cuDoubleComplex(f2u * xi0, f3u2x1 - cosw_p_xi0)
    );

    const cuDoubleComplex b11 = cuCmulf_cd(
        cuCsub(cuCadd(cuCmulf_cd(rj1, f2u), cuCmulf_cd(rj2, f_3u2_m_w2)),
            cuCmulf_cd(f1, f_30u2_p_2w2)),
        minus ? -dnormb : dnormb
    );
    const cuDoubleComplex b21 = cuCmulf_cd(
        cuCsub(cuCsub(rj1, cuCmulf_cd(rj2, f3u)),
            cuCmulf_cd(f1, f24u)),
        dnormb
    );

    //j = 2, r_2^(1) and r_2^(2)
    // (62)
    rj1 = cuCadd(
        cuCmul(exp2u, make_cuDoubleComplex(0.0, 2.0)),
        cuCmul(expm1u, make_cuDoubleComplex(-f3u * xi0, cosw - f3xi0))
    );
    // (65)
    rj2 = cuCmul(expm1u, make_cuDoubleComplex(xi0, -f3u * xi1));

    const cuDoubleComplex b12 = cuCmulf_cd(
        cuCsub(cuCadd(cuCmulf_cd(rj1, f2u), cuCmulf_cd(rj2, f_3u2_m_w2)),
            cuCmulf_cd(f2, f_30u2_p_2w2)),
        dnormb
    );
    const cuDoubleComplex b22 = cuCmulf_cd(
        cuCsub(cuCsub(rj1, cuCmulf_cd(rj2, f3u)),
            cuCmulf_cd(f2, f24u)),
        minus ? -dnormb : dnormb
    );

    //===============================================
    //B1, B2 term
    B = Q2; 
    _mul(B, _cToRealC(minus ? cuConj(b12) : b12));
    _add(B, _mulC(Q, _cToRealC(minus ? cuConj(b11) : b11)));
    _add(B, _cToRealC(minus ? cuConj(b10) : b10));

    B = _dagmulC(Sigma, B);
    _mul(B, U);
    CLGComplex traceres = _tr(B);
    B = Q;
    _mul(B, traceres);
    _add(res, B);

    B = Q2;
    _mul(B, _cToRealC(minus ? cuConj(b22) : b22));
    _add(B, _mulC(Q, _cToRealC(minus ? cuConj(b21) : b21)));
    _add(B, _cToRealC(minus ? cuConj(b20) : b20));

    B = _dagmulC(Sigma, B);
    _mul(B, U);
    traceres = _tr(B);
    B = Q2;
    _mul(B, traceres);
    _add(res, B);
    _th(res);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStoutForceLambda_NoQ(
    const deviceGauge* __restrict__ pOriginalGauge,
    const deviceGauge* __restrict__ pf0,
    const deviceGauge* __restrict__ pC,
    deviceGauge* pLambda,
    const cuDoubleComplex* __restrict__ f012_exp2u_expm1u,
    const DOUBLE* __restrict__ u_w2_cosw_xi0_denorm,
    const BYTE* __restrict__ expCondition)
{
    intokernalDir_NoDir;
    const BYTE& minus = expCondition[uiLinkIndex];

    if (2 == minus)
    {
        //when Q is zero, Lambda is also zero
        pLambda[uiLinkIndex] = _makeZero<deviceGauge>();
        return;
    }

    deviceGauge& res = pLambda[uiLinkIndex];
    deviceGauge Q = pC[uiLinkIndex];

    _muldag(Q, pOriginalGauge[uiLinkIndex]);
    _ta(Q);
    _mul(Q, _make_cuComplex(F(0.0), -F(1.0)));

    deviceGauge Q2 = Q;
    _mul(Q2, Q);

    const deviceGauge& U = pOriginalGauge[uiLinkIndex];
    const deviceGauge& Sigma = pf0[uiLinkIndex];

    const UINT dparamstart = uiLinkIndex * 5U;
    const DOUBLE& u = u_w2_cosw_xi0_denorm[dparamstart];
    const DOUBLE& w2 = u_w2_cosw_xi0_denorm[dparamstart + 1];
    const DOUBLE& cosw = u_w2_cosw_xi0_denorm[dparamstart + 2];
    const DOUBLE& xi0 = u_w2_cosw_xi0_denorm[dparamstart + 3];
    const DOUBLE& dnorm = u_w2_cosw_xi0_denorm[dparamstart + 4];
    const cuDoubleComplex& f0 = f012_exp2u_expm1u[dparamstart];
    const cuDoubleComplex& f1 = f012_exp2u_expm1u[dparamstart + 1];
    const cuDoubleComplex& f2 = f012_exp2u_expm1u[dparamstart + 2];
    const cuDoubleComplex& exp2u = f012_exp2u_expm1u[dparamstart + 3];
    const cuDoubleComplex& expm1u = f012_exp2u_expm1u[dparamstart + 4];

    deviceGauge B = U;
    _muldag(B, Sigma); //B = U.Sigma
    res = B;
    _mul(res, Q); //res = U.Sigma.Q
    _add(res, _mulC(Q, B)); // res = U.Sigma.Q + Q.U.Sigma
    _mul(res, _cToRealC(expCondition[uiLinkIndex] ? cuConj(f2) : f2)); // res = f2 U.Sigma.Q + f2 Q.U.Sigma
    _mul(B, _cToRealC(expCondition[uiLinkIndex] ? cuConj(cuCmulf_cd(f1, -1.0)) : f1)); //B = f1 U.Sigma
    _add(res, B); //res = f1 U.Sigma + f2 U.Sigma.Q + f2 Q.U.Sigma

    //calculate six b
    const DOUBLE u2 = u * u;
    const DOUBLE xi1 = (abs(w2) < _CLG_FLT_MIN) ? (cosw - xi0) : ((cosw - xi0) / w2);
    const DOUBLE dnormb = dnorm * dnorm * 0.5;
    const DOUBLE f2u = 2.0 * u;
    const DOUBLE f3u = 3.0 * u;
    const DOUBLE f24u = 24.0 * u;
    const DOUBLE f3u2 = 3.0 * u2;
    const DOUBLE f3u2x1 = f3u2 * xi1;
    const DOUBLE f3xi0 = 3.0 * xi0;
    const DOUBLE cosw_p_xi0 = cosw + xi0;
    const DOUBLE f_3u2_m_w2 = f3u2 - w2;
    const DOUBLE f_30u2_p_2w2 = 2.0 * (15.0 * u2 + w2);

    //j = 0, r_0^(1) and r_0^(2)
    // (60)
    cuDoubleComplex rj1 = cuCadd(
        cuCmul(exp2u, make_cuDoubleComplex(f2u, 2.0 * (u2 - w2))),
        cuCmul(expm1u, make_cuDoubleComplex(16.0 * u * cosw + f2u * (f3u2 + w2) * xi0, 2.0 * (9.0 * u2 + w2) * xi0 - 8.0 * u2 * cosw))
    );
    // (63)
    cuDoubleComplex rj2 = cuCadd(
        cuCmulf_cd(exp2u, -2.0),
        cuCmul(expm1u, make_cuDoubleComplex(-8.0 * u2 * xi0, f2u * (cosw_p_xi0 + f3u2x1)))
    );

    const cuDoubleComplex b10 = cuCmulf_cd(
        cuCsub(cuCadd(cuCmulf_cd(rj1, f2u), cuCmulf_cd(rj2, f_3u2_m_w2)),
            cuCmulf_cd(f0, f_30u2_p_2w2)),
        dnormb
    );
    const cuDoubleComplex b20 = cuCmulf_cd(
        cuCsub(cuCsub(rj1, cuCmulf_cd(rj2, f3u)),
            cuCmulf_cd(f0, f24u)),
        minus ? -dnormb : dnormb
    );

    //j = 1, r_1^(1) and r_1^(2)
    // (61)
    rj1 = cuCadd(
        cuCmul(exp2u, make_cuDoubleComplex(2.0, 4.0 * u)),
        cuCmul(expm1u, make_cuDoubleComplex(f_3u2_m_w2 * xi0 - 2.0 * cosw, f2u * (cosw + f3xi0)))
    );
    // (64)
    rj2 = cuCmul(expm1u,
        make_cuDoubleComplex(f2u * xi0, f3u2x1 - cosw_p_xi0)
    );

    const cuDoubleComplex b11 = cuCmulf_cd(
        cuCsub(cuCadd(cuCmulf_cd(rj1, f2u), cuCmulf_cd(rj2, f_3u2_m_w2)),
            cuCmulf_cd(f1, f_30u2_p_2w2)),
        minus ? -dnormb : dnormb
    );
    const cuDoubleComplex b21 = cuCmulf_cd(
        cuCsub(cuCsub(rj1, cuCmulf_cd(rj2, f3u)),
            cuCmulf_cd(f1, f24u)),
        dnormb
    );

    //j = 2, r_2^(1) and r_2^(2)
    // (62)
    rj1 = cuCadd(
        cuCmul(exp2u, make_cuDoubleComplex(0.0, 2.0)),
        cuCmul(expm1u, make_cuDoubleComplex(-f3u * xi0, cosw - f3xi0))
    );
    // (65)
    rj2 = cuCmul(expm1u, make_cuDoubleComplex(xi0, -f3u * xi1));

    const cuDoubleComplex b12 = cuCmulf_cd(
        cuCsub(cuCadd(cuCmulf_cd(rj1, f2u), cuCmulf_cd(rj2, f_3u2_m_w2)),
            cuCmulf_cd(f2, f_30u2_p_2w2)),
        dnormb
    );
    const cuDoubleComplex b22 = cuCmulf_cd(
        cuCsub(cuCsub(rj1, cuCmulf_cd(rj2, f3u)),
            cuCmulf_cd(f2, f24u)),
        minus ? -dnormb : dnormb
    );

    //===============================================
    //B1, B2 term
    B = Q2;
    _mul(B, _cToRealC(minus ? cuConj(b12) : b12));
    _add(B, _mulC(Q, _cToRealC(minus ? cuConj(b11) : b11)));
    _add(B, _cToRealC(minus ? cuConj(b10) : b10));

    B = _dagmulC(Sigma, B);
    _mul(B, U);
    CLGComplex traceres = _tr(B);
    B = Q;
    _mul(B, traceres);
    _add(res, B);

    B = Q2;
    _mul(B, _cToRealC(minus ? cuConj(b22) : b22));
    _add(B, _mulC(Q, _cToRealC(minus ? cuConj(b21) : b21)));
    _add(B, _cToRealC(minus ? cuConj(b20) : b20));

    B = _dagmulC(Sigma, B);
    _mul(B, U);
    traceres = _tr(B);
    B = Q2;
    _mul(B, traceres);
    _add(res, B);
    _th(res);
}

/**
* First two terms of smearing force of (75) of hep-lat/0311018
* Note that, after Lambda is calculated, f0 and C are not used anymore
* 
* Note that, in CLGLib, we calculate f^{dagger}
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStoutForceFirstTwoTerms(
    deviceGauge* pf0,
    const deviceGauge* __restrict__ pExpIQ,
    const deviceGauge* __restrict__ pC,
    const deviceGauge* __restrict__ pLambda,
    const BYTE* __restrict__ expCondition)
{
    intokernalDir_NoDir;

    if (2 == expCondition[uiLinkIndex])
    {
        //when Q is zero, force is not smeared, so keep it
        return;
    }
    //if (0 == uiLinkIndex)
    //{
    //    _print(pf0[uiLinkIndex], "f0");
    //}

    deviceGauge toadd = pLambda[uiLinkIndex];
    _mul(toadd, pC[uiLinkIndex]);
    _mul(toadd, _make_cuComplex(F(0.0), -F(1.0)));

    pf0[uiLinkIndex] = _dagmulC(pExpIQ[uiLinkIndex], pf0[uiLinkIndex]);
    _add(pf0[uiLinkIndex], toadd);
}

/**
* calculate force also using staple, last part of (75) of hep-lat/0311018
* 
* Put Lambda alone with U
* for U, it was Lambda_dagger U with plus sign
* for U_dagger, it was U_dagger, Lambda_dagger with minus sign
* 
* After all, dagger all, and then, multiply (i)
* 
* Calculate forward and backward seperately.
* Calculate 3 terms together
* 
* LU1 U2 U3
* U1 LU2 U3
* U1 U2 LU3
* 
* X = LU2
* X = X.U3 = (LU2.U3)
* Y = LU3
* Y = U2.Y = (U2 LU3)
* X = X + Y = (LU2.U3 + U2.LU3)
* X = U1.X = U1.(LU2.U3 + U2.LU3)
* Y = LU1
* Y = Y.U2 = (LU1.U2)
* Y = Y.U3 = (LU1.U2.U3)
* X = X + Y
* = U1.(LU2.U3 + U2.LU3) + LU1.U2.U3
* = LU1.U2.U3 + U1.LU2.U3 + U1.U2.LU3
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStoutFat3Force(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pLambda,
    const BYTE* __restrict__ expCondition,
    deviceGauge* pforce,
    Real rhojk, Real rho4mu,
    const SIndex* __restrict__ pCachedStapleIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    BYTE byFieldId)
{
    intokernalEDir(6U);

    if (2 == expCondition[uiLinkIndex])
    {
        //when Q is zero, force is not smeared, so keep it
        return;
    }

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1; //3
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1; //3*6 = 18
#endif
    const UINT byIdxStart = plaqLengthm1 * elementIdx + uiLinkIndex * plaqCountAllLink;
    const SIndex& link1 = pCachedStapleIndex[byIdxStart];
    const SIndex& link2 = pCachedStapleIndex[byIdxStart + 1];
    const SIndex& link3 = pCachedStapleIndex[byIdxStart + 2];
    const deviceSU3& U1 = _deviceGetGaugeBCT(byFieldId, pDeviceData, link1);
    const deviceSU3& U2 = _deviceGetGaugeBCT(byFieldId, pDeviceData, link2);
    const deviceSU3& U3 = _deviceGetGaugeBCT(byFieldId, pDeviceData, link3);
    const deviceSU3& Lam1 = _deviceGetGaugeBCDirForceSIndexT(pLambda, link1);
    const deviceSU3& Lam2 = _deviceGetGaugeBCDirForceSIndexT(pLambda, link2);
    const deviceSU3& Lam3 = _deviceGetGaugeBCDirForceSIndexT(pLambda, link3);
    //link2 is always not dagger (consider only torus boundary because hep-lat/0311018 says translational boundary is used)

    deviceGauge X = Lam1;
    _mul(X, U1);
    deviceGauge Y = U1;
    
    if (link1.NeedToDagger())
    {
        _dagger(X);
        _oppo(X);
        _dagger(Y);
    }
    
    if (link2.NeedToDagger())
    {
        _muldag(X, U2);

        _muldag(Y, U2);
        _mul(Y, Lam2);
        _oppo(Y);
    }
    else
    {
        _mul(X, U2);
        _mul(Y, Lam2);
        _mul(Y, U2);
    }
    _add(X, Y); //X = LU1 U2 + U1 LU2

    Y = U1;
    if (link1.NeedToDagger())
    {
        _dagger(Y);
    }
    if (link2.NeedToDagger())
    {
        _muldag(Y, U2);
    }
    else
    {
        _mul(Y, U2);
    }

    //Y = U1 U2
    if (link3.NeedToDagger())
    {
        _muldag(X, U3);
        _muldag(Y, U3);
        _mul(Y, Lam3);
        _oppo(Y);
    }
    else
    {
        _mul(X, U3);
        _mul(Y, Lam3);
        _mul(Y, U3);
    }

    _add(X, Y);

    if (3 == dir || elementIdx >= 4)
    {
        _mul(X, _make_cuComplex(F(0.0), rho4mu));
    }
    else
    {
        _mul(X, _make_cuComplex(F(0.0), rhojk));
    }

    #pragma unroll
    for (BYTE k = 0U; k < 6U; ++k)
    {
        if (k == elementIdx)
        {
            _add(pforce[uiLinkIndex], X);
        }
        __syncthreads();
    }
}

#pragma endregion

#pragma region Smearings

void CGaugeSmearingStoutSU3::GaugeSmearingFull(CFieldGauge* pGauge, const CFieldGauge* pOrignal, CFieldGauge* pStaple, UBOOL bProject)
{
    _RECORD(CGaugeSmearingStoutSU3::GaugeSmearingFull);
    appParanoiac(_T("CGaugeSmearingStoutSU3::GaugeSmearingFull\n"));
    if (NULL == m_pC)
    {
        m_pC = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
        m_pQ = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
        m_pQ2 = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
        m_pExpIQ = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
        m_pLambda = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
    }

    CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pGauge);
    const CFieldGaugeSU3* pGaugeOriginalSU3 = dynamic_cast<const CFieldGaugeSU3*>(pOrignal);
    CFieldGaugeSU3* pCSU3 = dynamic_cast<CFieldGaugeSU3*>(m_pC);
    CFieldGaugeSU3* pQSU3 = dynamic_cast<CFieldGaugeSU3*>(m_pQ);
    CFieldGaugeSU3* pQ2SU3 = dynamic_cast<CFieldGaugeSU3*>(m_pQ2);
    CFieldGaugeSU3* pExpSU3 = dynamic_cast<CFieldGaugeSU3*>(m_pExpIQ);

    //do the (1-rho) + rho fat3 smearing
    preparethreadDir;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelStoutFat3<deviceSU3>, block, threads,
        pGaugeOriginalSU3->m_pDeviceData,
        pCSU3->m_pDeviceData,
        m_fRhojk, m_fRho4mu,
        appGetLattice()->m_pIndexCache->m_pStappleCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        m_byFieldId
    );
#else
    _LAUNCH_KERNEL(_kernelStoutFat3<deviceSU3>, block, threads,
        pGaugeOriginalSU3->m_pDeviceData,
        pCSU3->m_pDeviceData,
        m_fRhojk, m_fRho4mu,
        appGetLattice()->m_pIndexCache->m_pStappleCache[m_byFieldId],
        m_byFieldId
    );
#endif

    _LAUNCH_KERNEL(_kernelStoutExpStep<deviceSU3>, block, threads,
        pGaugeOriginalSU3->m_pDeviceData,
        pCSU3->m_pDeviceData,
        pQSU3->m_pDeviceData,
        pQ2SU3->m_pDeviceData,
        pExpSU3->m_pDeviceData,
        pGaugeSU3->m_pDeviceData,
        m_f012_exp2u_expm1u,
        m_u_w2_cosw_xi0_denorm,
        m_byMinus
    );
    pGauge->FixBoundary(EFB_Field);
}

void CGaugeSmearingStoutSU3::DerivateOnUFull(const CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, CFieldGauge* pf0) const
{
    _RECORD(CGaugeSmearingStoutSU3::DerivateOnUFull);
    pf0->FixBoundary(EFB_Force);
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pOrignalGauge);
    const CFieldGaugeSU3* pCSU3 = dynamic_cast<const CFieldGaugeSU3*>(m_pC);
    const CFieldGaugeSU3* pQSU3 = dynamic_cast<const CFieldGaugeSU3*>(m_pQ);
    const CFieldGaugeSU3* pQ2SU3 = dynamic_cast<const CFieldGaugeSU3*>(m_pQ2);
    const CFieldGaugeSU3* pExpSU3 = dynamic_cast<const CFieldGaugeSU3*>(m_pExpIQ);
    CFieldGaugeSU3* pLambdaSU3 = dynamic_cast<CFieldGaugeSU3*>(m_pLambda);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pf0);

    {
        preparethreadDir;

        _LAUNCH_KERNEL(_kernelStoutForceLambda<deviceSU3>, block, threads,
            pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData,
            pQSU3->m_pDeviceData,
            pQ2SU3->m_pDeviceData,
            pLambdaSU3->m_pDeviceData,
            m_f012_exp2u_expm1u,
            m_u_w2_cosw_xi0_denorm,
            m_byMinus
        );

        _LAUNCH_KERNEL(_kernelStoutForceFirstTwoTerms<deviceSU3>, block, threads,
            pForceSU3->m_pDeviceData,
            pExpSU3->m_pDeviceData,
            pCSU3->m_pDeviceData,
            pLambdaSU3->m_pDeviceData,
            m_byMinus
        );
    }
    {
        preparethreadEDir(6U);

#if !_CLG_ASSUME_SQUARE_LATTICE

        _LAUNCH_KERNEL(_kernelStoutFat3Force<deviceSU3>, block, threads,
            pGaugeSU3->m_pDeviceData,
            pLambdaSU3->m_pDeviceData,
            m_byMinus,
            pForceSU3->m_pDeviceData,
            m_fRhojk, m_fRho4mu,
            appGetLattice()->m_pIndexCache->m_pStappleCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
            appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
            m_byFieldId
        );
#else
        _LAUNCH_KERNEL(_kernelStoutFat3Force<deviceSU3>, block, threads,
            pGaugeSU3->m_pDeviceData,
            pLambdaSU3->m_pDeviceData,
            m_byMinus,
            pForceSU3->m_pDeviceData,
            m_fRhojk, m_fRho4mu,
            appGetLattice()->m_pIndexCache->m_pStappleCache[m_byFieldId],
            m_byFieldId
        );
#endif
    }
    pf0->FixBoundary(EFB_Force);
}

void CGaugeSmearingStoutSU3::GaugeSmearingMedian(CFieldGauge* pGauge, const CFieldGauge* pOrignal, CFieldGauge* pStaple, UBOOL bProject)
{
    _RECORD(CGaugeSmearingStoutSU3::GaugeSmearingMedian);
    appParanoiac(_T("CGaugeSmearingStoutSU3::GaugeSmearingMedian\n"));
    if (NULL == m_pC)
    {
        m_pC = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
        m_pExpIQ = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
        m_pLambda = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(m_byFieldId, _T(__FILE__), __LINE__));
    }

    CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pGauge);
    const CFieldGaugeSU3* pGaugeOriginalSU3 = dynamic_cast<const CFieldGaugeSU3*>(pOrignal);
    CFieldGaugeSU3* pCSU3 = dynamic_cast<CFieldGaugeSU3*>(m_pC);
    CFieldGaugeSU3* pExpSU3 = dynamic_cast<CFieldGaugeSU3*>(m_pExpIQ);

    //do the (1-rho) + rho fat3 smearing
    preparethreadDir;
#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelStoutFat3<deviceSU3>, block, threads,
        pGaugeOriginalSU3->m_pDeviceData,
        pCSU3->m_pDeviceData,
        m_fRhojk, m_fRho4mu,
        appGetLattice()->m_pIndexCache->m_pStappleCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        m_byFieldId
    );
#else
    _LAUNCH_KERNEL(_kernelStoutFat3<deviceSU3>, block, threads,
        pGaugeOriginalSU3->m_pDeviceData,
        pCSU3->m_pDeviceData,
        m_fRhojk, m_fRho4mu,
        appGetLattice()->m_pIndexCache->m_pStappleCache[m_byFieldId],
        m_byFieldId
    );
#endif

    _LAUNCH_KERNEL(_kernelStoutExpStep_NoQ<deviceSU3>, block, threads,
        pGaugeOriginalSU3->m_pDeviceData,
        pCSU3->m_pDeviceData,
        pExpSU3->m_pDeviceData,
        pGaugeSU3->m_pDeviceData,
        m_f012_exp2u_expm1u,
        m_u_w2_cosw_xi0_denorm,
        m_byMinus
    );
    pGauge->FixBoundary(EFB_Field);
}

void CGaugeSmearingStoutSU3::DerivateOnUMedian(const CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, CFieldGauge* pf0) const
{
    _RECORD(CGaugeSmearingStoutSU3::DerivateOnUMedian);
    pf0->FixBoundary(EFB_Force);
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pOrignalGauge);
    const CFieldGaugeSU3* pCSU3 = dynamic_cast<const CFieldGaugeSU3*>(m_pC);
    const CFieldGaugeSU3* pExpSU3 = dynamic_cast<const CFieldGaugeSU3*>(m_pExpIQ);
    CFieldGaugeSU3* pLambdaSU3 = dynamic_cast<CFieldGaugeSU3*>(m_pLambda);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pf0);

    {
        preparethreadDir;

        _LAUNCH_KERNEL(_kernelStoutForceLambda_NoQ<deviceSU3>, block, threads,
            pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData,
            pCSU3->m_pDeviceData,
            pLambdaSU3->m_pDeviceData,
            m_f012_exp2u_expm1u,
            m_u_w2_cosw_xi0_denorm,
            m_byMinus
        );

        _LAUNCH_KERNEL(_kernelStoutForceFirstTwoTerms<deviceSU3>, block, threads,
            pForceSU3->m_pDeviceData,
            pExpSU3->m_pDeviceData,
            pCSU3->m_pDeviceData,
            pLambdaSU3->m_pDeviceData,
            m_byMinus
        );
    }
    {
        preparethreadEDir(6U);

#if !_CLG_ASSUME_SQUARE_LATTICE

        _LAUNCH_KERNEL(_kernelStoutFat3Force<deviceSU3>, block, threads,
            pGaugeSU3->m_pDeviceData,
            pLambdaSU3->m_pDeviceData,
            m_byMinus,
            pForceSU3->m_pDeviceData,
            m_fRhojk, m_fRho4mu,
            appGetLattice()->m_pIndexCache->m_pStappleCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
            appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
            m_byFieldId
        );
#else
        _LAUNCH_KERNEL(_kernelStoutFat3Force<deviceSU3>, block, threads,
            pGaugeSU3->m_pDeviceData,
            pLambdaSU3->m_pDeviceData,
            m_byMinus,
            pForceSU3->m_pDeviceData,
            m_fRhojk, m_fRho4mu,
            appGetLattice()->m_pIndexCache->m_pStappleCache[m_byFieldId],
            m_byFieldId
        );
#endif
    }
    pf0->FixBoundary(EFB_Force);
}


#pragma endregion

CGaugeSmearingStoutSU3::CGaugeSmearingStoutSU3()
    : CGaugeSmearing()
    , m_fRhojk(F(0.25))
    , m_fRho4mu(F(0.25))
    , m_pC(NULL)
    , m_pQ(NULL)
    , m_pQ2(NULL)
    , m_pExpIQ(NULL)
    , m_pLambda(NULL)
    , m_f012_exp2u_expm1u(NULL)
    , m_u_w2_cosw_xi0_denorm(NULL)
    , m_byMinus(NULL)
    , m_eCache(ESLC_Full)
{
    checkCudaErrors(__cudaMalloc((void**)&m_f012_exp2u_expm1u, sizeof(cuDoubleComplex) * _HC_Volume * _HC_Dir * 5));
    checkCudaErrors(__cudaMalloc((void**)&m_u_w2_cosw_xi0_denorm, sizeof(DOUBLE) * _HC_Volume * _HC_Dir * 5));
    checkCudaErrors(__cudaMalloc((void**)&m_byMinus, sizeof(BYTE) * _HC_Volume * _HC_Dir));
}

CGaugeSmearingStoutSU3::~CGaugeSmearingStoutSU3()
{
    checkCudaErrors(__cudaFree(m_f012_exp2u_expm1u));
    checkCudaErrors(__cudaFree(m_u_w2_cosw_xi0_denorm));
    checkCudaErrors(__cudaFree(m_byMinus));

    //if (NULL != m_pC)
    //{
    //    m_pC->Return();
    //}

    //if (NULL != m_pQ)
    //{
    //    m_pQ->Return();
    //}

    //if (NULL != m_pQ2)
    //{
    //    m_pQ2->Return();
    //}

    //if (NULL != m_pExpIQ)
    //{
    //    m_pExpIQ->Return();
    //}

    //if (NULL != m_pLambda)
    //{
    //    m_pLambda->Return();
    //}
}

void CGaugeSmearingStoutSU3::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    CGaugeSmearing::Initial(pOwner, params);
    m_bCalledWhenUpdate = TRUE;
    m_fRhojk = F(0.25);
    if (params.FetchValueReal(_T("rho"), m_fRhojk))
    {
        m_fRho4mu = m_fRhojk;
    }
    else
    {
        m_fRhojk = F(0.25);
        params.FetchValueReal(_T("rhojk"), m_fRhojk);
        m_fRho4mu = F(0.25);
        params.FetchValueReal(_T("rho4mu"), m_fRho4mu);
    }

    CCString sValue = _T("ESLC_Full");
    if (params.FetchStringValue(_T("Cache"), sValue))
    {
        m_eCache = __STRING_TO_ENUM(EStoutLinkCache, sValue);
    }
}

CCString CGaugeSmearingStoutSU3::GetInfos(const CCString &tab) const
{
    CCString sRet = CGaugeSmearing::GetInfos(tab);
    sRet = sRet + tab + _T("rhojk : ") + appToString(m_fRhojk) + _T("\n");
    sRet = sRet + tab + _T("rho4mu : ") + appToString(m_fRho4mu) + _T("\n");
    sRet = sRet + tab + _T("Cache : ") + __ENUM_TO_STRING(EStoutLinkCache, m_eCache).c_str() + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================