//=============================================================================
// FILENAME : CGaugeFixingCoulombCornell.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [09/20/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernels

#pragma region Cornell Steepest Descend

#if !_CLG_DOUBLEFLOAT
/**
 * A_mu (n) = TA(U _mu (n))/i
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateA3D(
    SBYTE uiT,
    const deviceSU3* __restrict__ pU,
    DOUBLE* pA11,
    cuDoubleComplex* pA12,
    cuDoubleComplex* pA13,
    DOUBLE* pA22,
    cuDoubleComplex* pA23)
{
    intokernalInt4_S;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //No need to calcuate A4, since we don't need it
    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        const UINT uiLinkIndex3D = uiSiteIndex3D * (uiDir - 1) + dir;
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            deviceSU3 su3A(pU[uiLinkIndex]);
            su3A.Ta();
            pA11[uiLinkIndex3D] = static_cast<DOUBLE>(su3A.m_me[0].y);
            pA12[uiLinkIndex3D] = _cToDouble(su3A.m_me[1]);
            pA13[uiLinkIndex3D] = _cToDouble(su3A.m_me[2]);
            pA22[uiLinkIndex3D] = static_cast<DOUBLE>(su3A.m_me[4].y);
            pA23[uiLinkIndex3D] = _cToDouble(su3A.m_me[5]);
        }
        else
        {
            pA11[uiLinkIndex3D] = 0.0;
            pA12[uiLinkIndex3D] = make_cuDoubleComplex(0.0, 0.0);
            pA13[uiLinkIndex3D] = make_cuDoubleComplex(0.0, 0.0);
            pA22[uiLinkIndex3D] = 0.0;
            pA23[uiLinkIndex3D] = make_cuDoubleComplex(0.0, 0.0);
        }
    }
}

#if _CLG_DEBUG
__global__ void _CLG_LAUNCH_BOUND_HALF
#else
__global__ void _CLG_LAUNCH_BOUND
#endif
_kernelCalculateA3DLog(
    SBYTE uiT,
    const deviceSU3* __restrict__ pU,
    DOUBLE* pA11,
    cuDoubleComplex* pA12,
    cuDoubleComplex* pA13,
    DOUBLE* pA22,
    cuDoubleComplex* pA23)
{
    intokernalInt4_S;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //No need to calcuate A4, since we don't need it
#pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        const UINT uiLinkIndex3D = uiSiteIndex3D * (uiDir - 1) + dir;
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            deviceSU3 su3A(pU[uiLinkIndex]);
            su3A = su3A.Log();
            pA11[uiLinkIndex3D] = static_cast<DOUBLE>(su3A.m_me[0].y);
            pA12[uiLinkIndex3D] = _cToDouble(su3A.m_me[1]);
            pA13[uiLinkIndex3D] = _cToDouble(su3A.m_me[2]);
            pA22[uiLinkIndex3D] = static_cast<DOUBLE>(su3A.m_me[4].y);
            pA23[uiLinkIndex3D] = _cToDouble(su3A.m_me[5]);
        }
        else
        {
            pA11[uiLinkIndex3D] = 0.0;
            pA12[uiLinkIndex3D] = make_cuDoubleComplex(0.0, 0.0);
            pA13[uiLinkIndex3D] = make_cuDoubleComplex(0.0, 0.0);
            pA22[uiLinkIndex3D] = 0.0;
            pA23[uiLinkIndex3D] = make_cuDoubleComplex(0.0, 0.0);
        }
    }
}

/**
 * Gamma(n) = Delta _{-mu} A(n) = \sum _mu (A_mu(n - mu) - A_mu(n))
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAGradient3D(
    BYTE byFieldId,
    SBYTE uiT,
    DOUBLE* pGamma11,
    cuDoubleComplex* pGamma12,
    cuDoubleComplex* pGamma13,
    DOUBLE* pGamma22,
    cuDoubleComplex* pGamma23,
    const DOUBLE* __restrict__ pA11,
    const cuDoubleComplex* __restrict__ pA12,
    const cuDoubleComplex* __restrict__ pA13,
    const DOUBLE* __restrict__ pA22,
    const cuDoubleComplex* __restrict__ pA23)
{
    intokernalInt4_S_Only3D;

    pGamma11[uiSiteIndex3D] = 0.0;
    pGamma12[uiSiteIndex3D] = make_cuDoubleComplex(0.0, 0.0);
    pGamma13[uiSiteIndex3D] = make_cuDoubleComplex(0.0, 0.0);
    pGamma22[uiSiteIndex3D] = 0.0;
    pGamma23[uiSiteIndex3D] = make_cuDoubleComplex(0.0, 0.0);

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    if (site.IsDirichlet())
    {
        return;
    }

    //No need to calcuate A4, since we don't need it
    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        const UINT uiLinkIndex3D = uiSiteIndex3D * (uiDir - 1) + dir;
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            pGamma11[uiSiteIndex3D] = pGamma11[uiSiteIndex3D] - pA11[uiLinkIndex3D];
            pGamma12[uiSiteIndex3D] = cuCsub(pGamma12[uiSiteIndex3D], pA12[uiLinkIndex3D]);
            pGamma13[uiSiteIndex3D] = cuCsub(pGamma13[uiSiteIndex3D], pA13[uiLinkIndex3D]);
            pGamma22[uiSiteIndex3D] = pGamma22[uiSiteIndex3D] - pA22[uiLinkIndex3D];
            pGamma23[uiSiteIndex3D] = cuCsub(pGamma23[uiSiteIndex3D], pA23[uiLinkIndex3D]);
        }

        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, __bck(dir));
        const SIndex& p_m_mu_dir = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__bi4(p_m_mu_site) + dir];
        const UINT uiLinkIndex2_3D = (p_m_mu_dir.m_uiSiteIndex / _DC_Lt) * (uiDir - 1) + p_m_mu_dir.m_byDir;

        //if (!__idx->_deviceIsBondOnSurface(p_m_mu, dir))
        if (!p_m_mu_dir.IsDirichlet())
        {
            if (p_m_mu_dir.NeedToDagger())
            {
                //dagger means A -> -A
                pGamma11[uiSiteIndex3D] = pGamma11[uiSiteIndex3D] - pA11[uiLinkIndex2_3D];
                pGamma12[uiSiteIndex3D] = cuCsub(pGamma12[uiSiteIndex3D], pA12[uiLinkIndex2_3D]);
                pGamma13[uiSiteIndex3D] = cuCsub(pGamma13[uiSiteIndex3D], pA13[uiLinkIndex2_3D]);
                pGamma22[uiSiteIndex3D] = pGamma22[uiSiteIndex3D] - pA22[uiLinkIndex2_3D];
                pGamma23[uiSiteIndex3D] = cuCsub(pGamma23[uiSiteIndex3D], pA23[uiLinkIndex2_3D]);
            }
            else
            {
                pGamma11[uiSiteIndex3D] = pGamma11[uiSiteIndex3D] + pA11[uiLinkIndex2_3D];
                pGamma12[uiSiteIndex3D] = cuCadd(pGamma12[uiSiteIndex3D], pA12[uiLinkIndex2_3D]);
                pGamma13[uiSiteIndex3D] = cuCadd(pGamma13[uiSiteIndex3D], pA13[uiLinkIndex2_3D]);
                pGamma22[uiSiteIndex3D] = pGamma22[uiSiteIndex3D] + pA22[uiLinkIndex2_3D];
                pGamma23[uiSiteIndex3D] = cuCadd(pGamma23[uiSiteIndex3D], pA23[uiLinkIndex2_3D]);
            }
        }
    }
}

/**
 * g(x)=exp(-i a Delta A)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateG3D(
    SBYTE uiT,
    deviceSU3* pG,
    const DOUBLE* __restrict__ pGamma11,
    const cuDoubleComplex* __restrict__ pGamma12,
    const cuDoubleComplex* __restrict__ pGamma13,
    const DOUBLE* __restrict__ pGamma22,
    const cuDoubleComplex* __restrict__ pGamma23,
    DOUBLE fAlpha)
{
    intokernalInt4_S_Only3D;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    if (site.IsDirichlet())
    {
        pG[uiSiteIndex3D] = deviceSU3::makeSU3Id();
    }
    else
    {
        const deviceSU3 pA = deviceSU3::makeSU3TA(
            _cToFloat(pGamma12[uiSiteIndex3D]), _cToFloat(pGamma13[uiSiteIndex3D]), _cToFloat(pGamma23[uiSiteIndex3D]),
            static_cast<Real>(pGamma11[uiSiteIndex3D]), static_cast<Real>(pGamma22[uiSiteIndex3D]));
        pG[uiSiteIndex3D] = (0 == _DC_ExpPrecision)
            ? pA.QuickExp(fAlpha)
            : pA.ExpReal(fAlpha, _DC_ExpPrecision);
    }
}
#else
/**
 * A_mu (n) = TA(U _mu (n))/i
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateA3D(
        SBYTE uiT,
        const deviceSU3* __restrict__ pU,
        Real* pA11,
        CLGComplex* pA12,
        CLGComplex* pA13,
        Real* pA22,
        CLGComplex* pA23)
{
    intokernalInt4_S;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir); 
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //No need to calcuate A4, since we don't need it
    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        const UINT uiLinkIndex3D = uiSiteIndex3D * (uiDir - 1) + dir;
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            deviceSU3 su3A(pU[uiLinkIndex]);
            su3A.Ta();
            pA11[uiLinkIndex3D] = su3A.m_me[0].y;
            pA12[uiLinkIndex3D] = su3A.m_me[1];
            pA13[uiLinkIndex3D] = su3A.m_me[2];
            pA22[uiLinkIndex3D] = su3A.m_me[4].y;
            pA23[uiLinkIndex3D] = su3A.m_me[5];
        }
        else
        {
            pA11[uiLinkIndex3D] = F(0.0);
            pA12[uiLinkIndex3D] = _zeroc;
            pA13[uiLinkIndex3D] = _zeroc;
            pA22[uiLinkIndex3D] = F(0.0);
            pA23[uiLinkIndex3D] = _zeroc;
        }
    }
}

#if _CLG_DEBUG
__global__ void _CLG_LAUNCH_BOUND_HALF
#else
__global__ void _CLG_LAUNCH_BOUND
#endif
_kernelCalculateA3DLog(
    SBYTE uiT,
    const deviceSU3* __restrict__ pU,
    Real* pA11,
    CLGComplex* pA12,
    CLGComplex* pA13,
    Real* pA22,
    CLGComplex* pA23)
{
    intokernalInt4_S;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //No need to calcuate A4, since we don't need it
#pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        const UINT uiLinkIndex3D = uiSiteIndex3D * (uiDir - 1) + dir;
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            deviceSU3 su3A(pU[uiLinkIndex]);
            su3A = su3A.Log();
            pA11[uiLinkIndex3D] = su3A.m_me[0].y;
            pA12[uiLinkIndex3D] = su3A.m_me[1];
            pA13[uiLinkIndex3D] = su3A.m_me[2];
            pA22[uiLinkIndex3D] = su3A.m_me[4].y;
            pA23[uiLinkIndex3D] = su3A.m_me[5];
        }
        else
        {
            pA11[uiLinkIndex3D] = F(0.0);
            pA12[uiLinkIndex3D] = _zeroc;
            pA13[uiLinkIndex3D] = _zeroc;
            pA22[uiLinkIndex3D] = F(0.0);
            pA23[uiLinkIndex3D] = _zeroc;
        }
    }
}

/**
 * Gamma(n) = Delta _{-mu} A(n) = \sum _mu (A_mu(n - mu) - A_mu(n))
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAGradient3D(
    BYTE byFieldId,
    SBYTE uiT,
    Real* pGamma11,
    CLGComplex* pGamma12,
    CLGComplex* pGamma13,
    Real* pGamma22,
    CLGComplex* pGamma23,
    const Real* __restrict__ pA11,
    const CLGComplex* __restrict__ pA12,
    const CLGComplex* __restrict__ pA13,
    const Real* __restrict__ pA22,
    const CLGComplex* __restrict__ pA23)
{
    intokernalInt4_S_Only3D;

    pGamma11[uiSiteIndex3D] = F(0.0);
    pGamma12[uiSiteIndex3D] = _zeroc;
    pGamma13[uiSiteIndex3D] = _zeroc;
    pGamma22[uiSiteIndex3D] = F(0.0);
    pGamma23[uiSiteIndex3D] = _zeroc;
    
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    if (site.IsDirichlet())
    {
        return;
    }

    //No need to calcuate A4, since we don't need it
    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        const UINT uiLinkIndex3D = uiSiteIndex3D * (uiDir - 1) + dir;
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            pGamma11[uiSiteIndex3D] = pGamma11[uiSiteIndex3D] - pA11[uiLinkIndex3D];
            pGamma12[uiSiteIndex3D] = _cuCsubf(pGamma12[uiSiteIndex3D], pA12[uiLinkIndex3D]);
            pGamma13[uiSiteIndex3D] = _cuCsubf(pGamma13[uiSiteIndex3D], pA13[uiLinkIndex3D]);
            pGamma22[uiSiteIndex3D] = pGamma22[uiSiteIndex3D] - pA22[uiLinkIndex3D];
            pGamma23[uiSiteIndex3D] = _cuCsubf(pGamma23[uiSiteIndex3D], pA23[uiLinkIndex3D]);
        }

        const SSmallInt4 p_m_mu_site = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(dir) - 1);
        const SIndex& p_m_mu_dir = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][__idx->_deviceGetBigIndex(p_m_mu_site) * uiDir + dir];
        const UINT uiLinkIndex2_3D = (p_m_mu_dir.m_uiSiteIndex / _DC_Lt) * (uiDir - 1) + dir;

        //if (!__idx->_deviceIsBondOnSurface(p_m_mu, dir))
        if (!p_m_mu_dir.IsDirichlet())
        {
            if (p_m_mu_dir.NeedToDagger())
            {
                pGamma11[uiSiteIndex3D] = pGamma11[uiSiteIndex3D] - pA11[uiLinkIndex2_3D];
                pGamma12[uiSiteIndex3D] = _cuCsubf(pGamma12[uiSiteIndex3D], pA12[uiLinkIndex2_3D]);
                pGamma13[uiSiteIndex3D] = _cuCsubf(pGamma13[uiSiteIndex3D], pA13[uiLinkIndex2_3D]);
                pGamma22[uiSiteIndex3D] = pGamma22[uiSiteIndex3D] - pA22[uiLinkIndex2_3D];
                pGamma23[uiSiteIndex3D] = _cuCsubf(pGamma23[uiSiteIndex3D], pA23[uiLinkIndex2_3D]);
            }
            else
            {
                pGamma11[uiSiteIndex3D] = pGamma11[uiSiteIndex3D] + pA11[uiLinkIndex2_3D];
                pGamma12[uiSiteIndex3D] = _cuCaddf(pGamma12[uiSiteIndex3D], pA12[uiLinkIndex2_3D]);
                pGamma13[uiSiteIndex3D] = _cuCaddf(pGamma13[uiSiteIndex3D], pA13[uiLinkIndex2_3D]);
                pGamma22[uiSiteIndex3D] = pGamma22[uiSiteIndex3D] + pA22[uiLinkIndex2_3D];
                pGamma23[uiSiteIndex3D] = _cuCaddf(pGamma23[uiSiteIndex3D], pA23[uiLinkIndex2_3D]);
            }
        }
    }
}

/**
 * g(x)=exp(-i a Delta A)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateG3D(
    SBYTE uiT,
    deviceSU3* pG,
    const Real* __restrict__ pGamma11,
    const CLGComplex* __restrict__ pGamma12,
    const CLGComplex* __restrict__ pGamma13,
    const Real* __restrict__ pGamma22,
    const CLGComplex* __restrict__ pGamma23,
    Real fAlpha)
{
    intokernalInt4_S_Only3D;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    if (site.IsDirichlet())
    {
        pG[uiSiteIndex3D] = deviceSU3::makeSU3Id();
    }
    else
    {
        deviceSU3 pA = deviceSU3::makeSU3TA(
            pGamma12[uiSiteIndex3D], pGamma13[uiSiteIndex3D], pGamma23[uiSiteIndex3D],
            pGamma11[uiSiteIndex3D], pGamma22[uiSiteIndex3D]);
        pG[uiSiteIndex3D] = (0 == _DC_ExpPrecision)
            ? pA.QuickExp(fAlpha)
            : pA.ExpReal(fAlpha, _DC_ExpPrecision);
    }
}
#endif

/**
 * g(n) U_mu(n) g(n+mu)^dagger
 * U_4(n) is left unchanged
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransform3D(
    BYTE byFieldId,
    SBYTE uiT,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
{
    intokernalInt4_S;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const deviceSU3 left(pGx[uiSiteIndex3D]);

    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, dir))
        {
            const UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
            deviceSU3 res(pGauge[uiLinkDir]);

            const SSmallInt4 p_p_mu_site = _deviceSmallInt4OffsetC(sSite4, __fwd(dir));
            const SIndex& site_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(p_p_mu_site)];

            if (!site_p_mu.IsDirichlet())
            {
                res.MulDagger(pGx[(site_p_mu.m_uiSiteIndex / _DC_Lt)]);
            }

            pGauge[uiLinkDir] = left.MulC(res);
        }
    }
}

/**
 * g(n) U_t(n)
 * U_t(n -t) g(n)^dagger
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransform3DT(
    BYTE byFieldId,
    SBYTE uiT,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
{
    intokernalInt4_S;

    //const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    //const BYTE uiDir2 = uiDir * 2;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
        if (!site.IsDirichlet())
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);
            pGauge[uiLinkIndex] = pGx[uiSiteIndex3D].MulC(pGauge[uiLinkIndex]);
        }
    }

    const SSmallInt4 p_m_t_site = _deviceSmallInt4OffsetC(sSite4, -4);
    const UINT p_m_t_bi = __idx->_deviceGetBigIndex(p_m_t_site);
    const SIndex& site_m_t_point = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][p_m_t_bi];
    const SIndex& site_m_t = __idx->m_pDeviceIndexLinkToSIndex[byFieldId][p_m_t_bi * _DC_Dir + 3];
    
    if (!site_m_t.IsDirichlet() && !site_m_t_point.IsDirichlet()) 
    {
        const UINT uiLinkIndex2 = _deviceGetLinkIndex(site_m_t.m_uiSiteIndex, 3);
        if (site_m_t.NeedToDagger())
        {
            //never here
            printf("ever here???\n");
            pGauge[uiLinkIndex2] = pGx[uiSiteIndex3D].MulC(pGauge[uiLinkIndex2]);
        }
        else
        {
            pGauge[uiLinkIndex2].MulDagger(pGx[uiSiteIndex3D]);
        }
    }
}


#if !_CLG_DOUBLEFLOAT
/**
* res = Tr[Delta A^2]
* If Delta A is a anti-Hermitian, Tr[Delta A^2] = 2 (|A12|^2+|A13|^2+|A23|^2 + |A11+A22|^2)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateTrAGradientSq3D(
    SBYTE uiT,
    DOUBLE* pDeviceRes,
    const DOUBLE* __restrict__ pDeltaA11,
    const cuDoubleComplex* __restrict__ pDeltaA12,
    const cuDoubleComplex* __restrict__ pDeltaA13,
    const DOUBLE* __restrict__ pDeltaA22,
    const cuDoubleComplex* __restrict__ pDeltaA23)
{
    intokernalInt4_S_Only3D;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    if (site.IsDirichlet())
    {
        pDeviceRes[uiSiteIndex3D] = 0.0;
    }

    const DOUBLE fAbs1 = cuCabs(pDeltaA12[uiSiteIndex3D]);
    const DOUBLE fAbs2 = cuCabs(pDeltaA13[uiSiteIndex3D]);
    const DOUBLE fAbs3 = cuCabs(pDeltaA23[uiSiteIndex3D]);
    const DOUBLE fM1122 = pDeltaA11[uiSiteIndex3D] + pDeltaA22[uiSiteIndex3D];
    pDeviceRes[uiSiteIndex3D] = 2.0 * (fAbs1 * fAbs1 + fAbs2 * fAbs2 + fAbs3 * fAbs3 + fM1122 * fM1122);
}
#else
/**
* res = Tr[Delta A^2]
* If Delta A is a anti-Hermitian, Tr[Delta A^2] = 2 (|A12|^2+|A13|^2+|A23|^2 + |A11+A22|^2)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateTrAGradientSq3D(
    SBYTE uiT,
    Real* pDeviceRes,
    const Real* __restrict__ pDeltaA11,
    const CLGComplex* __restrict__ pDeltaA12,
    const CLGComplex* __restrict__ pDeltaA13,
    const Real* __restrict__ pDeltaA22,
    const CLGComplex* __restrict__ pDeltaA23)
{
    intokernalInt4_S_Only3D;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];
    if (site.IsDirichlet())
    {
        pDeviceRes[uiSiteIndex3D] = F(0.0);
    }

    const Real fAbs1 = _cuCabsf(pDeltaA12[uiSiteIndex3D]);
    const Real fAbs2 = _cuCabsf(pDeltaA13[uiSiteIndex3D]);
    const Real fAbs3 = _cuCabsf(pDeltaA23[uiSiteIndex3D]);
    const Real fM1122 = pDeltaA11[uiSiteIndex3D] + pDeltaA22[uiSiteIndex3D];
    pDeviceRes[uiSiteIndex3D] = F(2.0) * (fAbs1 * fAbs1 + fAbs2 * fAbs2 + fAbs3 * fAbs3 + fM1122 * fM1122);
}
#endif

#pragma endregion

#pragma region FFT accelaration

#if !_CLG_DOUBLEFLOAT
__global__ void _CLG_LAUNCH_BOUND
_kernelBakeMomentumTable3D(DOUBLE* pP, UINT uiV)
{
    const SBYTE uiT = 0;
    intokernalInt4_S_Only3D;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    DOUBLE fDenorm = static_cast<DOUBLE>(uiDir - 1);
    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        fDenorm -= cos(2.0 * PI * sSite4.m_byData4[dir] / static_cast<DOUBLE>(_constIntegers[ECI_Lx + dir]));
    }

    if (abs(fDenorm) < _CLG_FLT_EPSILON)
    {
        fDenorm = 0.5;
    }

    //when p^2=0, p^2=1, or p^2 = 2(Nd - sum cos)
    //4 * 3 / 2(Nd - sum cos)
    pP[uiSiteIndex3D] = 6.0 / (fDenorm * uiV);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTRtoC3D(const DOUBLE* __restrict__ realBuffer, cuDoubleComplex* complexBuffer)
{
    const SBYTE uiT = 0;
    intokernalInt4_S_Only3D;
    complexBuffer[uiSiteIndex3D] = make_cuDoubleComplex(realBuffer[uiSiteIndex3D], 0.0);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTCtoR3D(const cuDoubleComplex* __restrict__ complexBuffer, DOUBLE* realBuffer)
{
    const SBYTE uiT = 0;
    intokernalInt4_S_Only3D;
    realBuffer[uiSiteIndex3D] = complexBuffer[uiSiteIndex3D].x;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTScale3D(const DOUBLE* __restrict__ pP, cuDoubleComplex* fftRes)
{
    const SBYTE uiT = 0;
    intokernalInt4_S_Only3D;
    fftRes[uiSiteIndex3D] = cuCmulf_cd(fftRes[uiSiteIndex3D], pP[uiSiteIndex3D]);
}
#else
__global__ void _CLG_LAUNCH_BOUND
_kernelBakeMomentumTable3D(Real* pP, UINT uiV)
{
    const SBYTE uiT = 0;
    intokernalInt4_S_Only3D;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    Real fDenorm = static_cast<Real>(uiDir - 1);
    for (BYTE dir = 0; dir < uiDir - 1; ++dir)
    {
        fDenorm -= _cos(F(2.0) * PI * sSite4.m_byData4[dir] / static_cast<Real>(_constIntegers[ECI_Lx + dir]));
    }

    if (abs(fDenorm) < _CLG_FLT_EPSILON)
    {
        fDenorm = F(0.5);
    }

    //when p^2=0, p^2=1, or p^2 = 2(Nd - sum cos)
    //4 * 3 / 2(Nd - sum cos)
    pP[uiSiteIndex3D] = F(6.0) / (fDenorm * uiV);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTRtoC3D(const Real* __restrict__ realBuffer, CLGComplex* complexBuffer)
{
    const SBYTE uiT = 0;
    intokernalInt4_S_Only3D;
    complexBuffer[uiSiteIndex3D] = _make_cuComplex(realBuffer[uiSiteIndex3D], F(0.0));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTCtoR3D(const CLGComplex* __restrict__ complexBuffer, Real* realBuffer)
{
    const SBYTE uiT = 0;
    intokernalInt4_S_Only3D;
    realBuffer[uiSiteIndex3D] = complexBuffer[uiSiteIndex3D].x;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFFTScale3D(const Real* __restrict__ pP, CLGComplex* fftRes)
{
    const SBYTE uiT = 0;
    intokernalInt4_S_Only3D;
    fftRes[uiSiteIndex3D] = cuCmulf_cr(fftRes[uiSiteIndex3D], pP[uiSiteIndex3D]);
}
#endif
#pragma endregion

#pragma endregion

__CLGIMPLEMENT_CLASS(CGaugeFixingCoulombCornell)

void CGaugeFixingCoulombCornell::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;

    const TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock();
    TArray<UINT> latticeDim;
    latticeDim.AddItem(_HC_Lx);
    latticeDim.AddItem(_HC_Ly);
    latticeDim.AddItem(_HC_Lz);
    TArray <UINT> decomp = _getDecompose(deviceConstraints, latticeDim);
    m_pHDecomp[0] = decomp[0];
    m_pHDecomp[1] = decomp[1];
    m_pHDecomp[2] = decomp[2];
    m_pHDecomp[3] = decomp[3];
    m_pHDecomp[4] = decomp[4];
    m_pHDecomp[5] = decomp[5];
    m_lstDims.RemoveAll();
    m_lstDims.AddItem(_HC_Lx);
    m_lstDims.AddItem(_HC_Ly);
    m_lstDims.AddItem(_HC_Lz);

    checkCudaErrors(cudaMalloc((void**)& m_pDDecomp, sizeof(UINT) * 6));

    //========== Initial Settings ==============

#if !_CLG_DOUBLEFLOAT
    if (!params.FetchValueDOUBLE(_T("Alpha"), m_fAlpha))
#else
    if (!params.FetchValueReal(_T("Alpha"), m_fAlpha))
#endif
    {
        appGeneral(_T("CGaugeFixingCoulombCornell: Alpha not set, set to 0.08 by defualt."));
    }
    
    if (!params.FetchValueReal(_T("Accuracy"), m_fAccuracy))
    {
        appGeneral(_T("CGaugeFixingCoulombCornell: Accuracy not set, set to 0.00000000001 by defualt."));
        if (m_fAccuracy < _CLG_FLT_EPSILON * F(2.0))
        {
            m_fAccuracy = _CLG_FLT_EPSILON * F(2.0);
            appGeneral(_T("Solver accuracy too small, set to be %2.18f\n"), m_fAccuracy);
        }
    }
    //m_fAccuracy = m_fAccuracy;

    INT iValue = static_cast<INT>(m_iMaxIterate);
    if (!params.FetchValueINT(_T("MaxIterate"), iValue))
    {
        appGeneral(_T("CGaugeFixingCoulombCornell: MaxIterate not set, set to 100000 by defualt."));
    }
    m_iMaxIterate = static_cast<UINT>(iValue);

    iValue = 1000;
    if (!params.FetchValueINT(_T("ShowErrorStep"), iValue))
    {
        appParanoiac(_T("CGaugeFixingCoulombCornell: ShowErrorStep not set, set to 1000 by defualt."));
    }
    m_iShowErrorStep = iValue;

    iValue = 1;
    if (!params.FetchValueINT(_T("FFT"), iValue))
    {
        appGeneral(_T("CGaugeFixingCoulombCornell: FFT not set, set to 1 by defualt."));
    }
    m_bFA = (0 != iValue);

#if !_CLG_DOUBLEFLOAT
    if (m_bFA)
    {
        appCrucial(_T("Do not use FFT for single float point\n"));
        m_bFA = FALSE;
    }
    //========== Initial Buffers ==============
    checkCudaErrors(cudaMalloc((void**)&m_pA11, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(DOUBLE)));
    checkCudaErrors(cudaMalloc((void**)&m_pA12, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(cuDoubleComplex)));
    checkCudaErrors(cudaMalloc((void**)&m_pA13, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(cuDoubleComplex)));
    checkCudaErrors(cudaMalloc((void**)&m_pA22, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(DOUBLE)));
    checkCudaErrors(cudaMalloc((void**)&m_pA23, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(cuDoubleComplex)));

    checkCudaErrors(cudaMalloc((void**)&m_pGamma11, _HC_Volume_xyz * sizeof(DOUBLE)));
    checkCudaErrors(cudaMalloc((void**)&m_pGamma12, _HC_Volume_xyz * sizeof(cuDoubleComplex)));
    checkCudaErrors(cudaMalloc((void**)&m_pGamma13, _HC_Volume_xyz * sizeof(cuDoubleComplex)));
    checkCudaErrors(cudaMalloc((void**)&m_pGamma22, _HC_Volume_xyz * sizeof(DOUBLE)));
    checkCudaErrors(cudaMalloc((void**)&m_pGamma23, _HC_Volume_xyz * sizeof(cuDoubleComplex)));
#else
    //========== Initial Buffers ==============
    checkCudaErrors(cudaMalloc((void**)& m_pA11, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(Real)));
    checkCudaErrors(cudaMalloc((void**)& m_pA12, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)& m_pA13, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)& m_pA22, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(Real)));
    checkCudaErrors(cudaMalloc((void**)& m_pA23, _HC_Volume_xyz * (_HC_Dir - 1) * sizeof(CLGComplex)));

    checkCudaErrors(cudaMalloc((void**)& m_pGamma11, _HC_Volume_xyz * sizeof(Real)));
    checkCudaErrors(cudaMalloc((void**)& m_pGamma12, _HC_Volume_xyz * sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)& m_pGamma13, _HC_Volume_xyz * sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)& m_pGamma22, _HC_Volume_xyz * sizeof(Real)));
    checkCudaErrors(cudaMalloc((void**)& m_pGamma23, _HC_Volume_xyz * sizeof(CLGComplex)));
#endif

    checkCudaErrors(cudaMalloc((void**)& m_pG, _HC_Volume_xyz * sizeof(deviceSU3)));
    if (m_bFA)
    {
#if !_CLG_DOUBLEFLOAT
        checkCudaErrors(cudaMalloc((void**)&m_pMomentumTable, _HC_Volume_xyz * sizeof(DOUBLE)));
        checkCudaErrors(cudaMalloc((void**)&m_pTempFFTBuffer, _HC_Volume_xyz * sizeof(cuDoubleComplex)));
#else
        checkCudaErrors(cudaMalloc((void**)& m_pMomentumTable, _HC_Volume_xyz * sizeof(Real)));
        checkCudaErrors(cudaMalloc((void**)& m_pTempFFTBuffer, _HC_Volume_xyz * sizeof(CLGComplex)));
#endif
        preparethread_S;
        _kernelBakeMomentumTable3D << <block, threads >> > (m_pMomentumTable, _HC_Volume_xyz);
    }
}

void CGaugeFixingCoulombCornell::GaugeFixing(CFieldGauge* pResGauge)
{
    if (NULL == pResGauge || EFT_GaugeSU3 != pResGauge->GetFieldType())
    {
        appCrucial(_T("CGaugeFixingLandauCornell only implemented with gauge SU3!\n"));
        return;
    }
    CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pResGauge);
    deviceSU3* pDeviceBufferPointer = pGaugeSU3->m_pDeviceData;
    for (SBYTE uiT = 0; uiT < static_cast<SBYTE>(_HC_Lt); ++uiT)
    {
        GaugeFixingOneTimeSlice(pDeviceBufferPointer, uiT, pGaugeSU3->m_byFieldId);
    }
}

void CGaugeFixingCoulombCornell::GaugeFixingOneTimeSlice(deviceSU3* pDeviceBufferPointer, SBYTE uiT, BYTE byFieldId)
{
    preparethread_S;
    m_iIterate = 0;
#if !_CLG_DOUBLEFLOAT
    DOUBLE fTheta = 0.0;
#else
    Real fTheta = F(0.0);
#endif

    while (m_iIterate < m_iMaxIterate)
    {
        //======= 1. Calculate Gamma    =========
        if (0 == _HC_ALog)
        {
            _kernelCalculateA3D << <block, threads >> > (
                uiT,
                pDeviceBufferPointer,
                m_pA11,
                m_pA12,
                m_pA13,
                m_pA22,
                m_pA23);
        }
        else
        {
            _kernelCalculateA3DLog << <block, threads >> > (
                uiT,
                pDeviceBufferPointer,
                m_pA11,
                m_pA12,
                m_pA13,
                m_pA22,
                m_pA23);
        }

        _kernelCalculateAGradient3D << <block, threads >> > (
            byFieldId,
            uiT,
            m_pGamma11,
            m_pGamma12,
            m_pGamma13,
            m_pGamma22,
            m_pGamma23,
            m_pA11,
            m_pA12,
            m_pA13,
            m_pA22,
            m_pA23);

        //======= 2. Calculate Theta    =========
        _kernelCalculateTrAGradientSq3D << <block, threads >> > (
            uiT,
            _D_RealThreadBuffer,
            m_pGamma11,
            m_pGamma12,
            m_pGamma13,
            m_pGamma22,
            m_pGamma23);
        fTheta = appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz) / (3 * _HC_Volume_xyz);
        if (m_iShowErrorStep > 0 && 0 == m_iIterate % m_iShowErrorStep)
        {
            appParanoiac(_T("Theta%d = %2.12f\n"), m_iIterate, fTheta);
        }

        if (fTheta < m_fAccuracy)
        {
            return;
        }

        //======= 3. FFT =========
        if (m_bFA)
        {
#if !_CLG_DOUBLEFLOAT
            CCLGFFTHelper::FFT3DWithXYZDouble(m_pGamma12, m_lstDims, TRUE);
            _kernelFFTScale3D << <block, threads >> > (m_pMomentumTable, m_pGamma12);
            CCLGFFTHelper::FFT3DWithXYZDouble(m_pGamma12, m_lstDims, FALSE);

            CCLGFFTHelper::FFT3DWithXYZDouble(m_pGamma13, m_lstDims, TRUE);
            _kernelFFTScale3D << <block, threads >> > (m_pMomentumTable, m_pGamma13);
            CCLGFFTHelper::FFT3DWithXYZDouble(m_pGamma13, m_lstDims, FALSE);

            CCLGFFTHelper::FFT3DWithXYZDouble(m_pGamma23, m_lstDims, TRUE);
            _kernelFFTScale3D << <block, threads >> > (m_pMomentumTable, m_pGamma23);
            CCLGFFTHelper::FFT3DWithXYZDouble(m_pGamma23, m_lstDims, FALSE);

            _kernelFFTRtoC3D << <block, threads >> > (m_pGamma11, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT3DWithXYZDouble(m_pTempFFTBuffer, m_lstDims, TRUE);
            _kernelFFTScale3D << <block, threads >> > (m_pMomentumTable, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT3DWithXYZDouble(m_pTempFFTBuffer, m_lstDims, FALSE);
            _kernelFFTCtoR3D << <block, threads >> > (m_pTempFFTBuffer, m_pGamma11);

            _kernelFFTRtoC3D << <block, threads >> > (m_pGamma22, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT3DWithXYZDouble(m_pTempFFTBuffer, m_lstDims, TRUE);
            _kernelFFTScale3D << <block, threads >> > (m_pMomentumTable, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT3DWithXYZDouble(m_pTempFFTBuffer, m_lstDims, FALSE);
            _kernelFFTCtoR3D << <block, threads >> > (m_pTempFFTBuffer, m_pGamma22); 
#else
            CCLGFFTHelper::FFT3DWithXYZ(m_pGamma12, m_lstDims, TRUE);
            _kernelFFTScale3D << <block, threads >> > (m_pMomentumTable, m_pGamma12);
            CCLGFFTHelper::FFT3DWithXYZ(m_pGamma12, m_lstDims, FALSE);

            CCLGFFTHelper::FFT3DWithXYZ(m_pGamma13, m_lstDims, TRUE);
            _kernelFFTScale3D << <block, threads >> > (m_pMomentumTable, m_pGamma13);
            CCLGFFTHelper::FFT3DWithXYZ(m_pGamma13, m_lstDims, FALSE);

            CCLGFFTHelper::FFT3DWithXYZ(m_pGamma23, m_lstDims, TRUE);
            _kernelFFTScale3D << <block, threads >> > (m_pMomentumTable, m_pGamma23);
            CCLGFFTHelper::FFT3DWithXYZ(m_pGamma23, m_lstDims, FALSE);

            _kernelFFTRtoC3D << <block, threads >> > (m_pGamma11, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT3DWithXYZ(m_pTempFFTBuffer, m_lstDims, TRUE);
            _kernelFFTScale3D << <block, threads >> > (m_pMomentumTable, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT3DWithXYZ(m_pTempFFTBuffer, m_lstDims, FALSE);
            _kernelFFTCtoR3D << <block, threads >> > (m_pTempFFTBuffer, m_pGamma11);

            _kernelFFTRtoC3D << <block, threads >> > (m_pGamma22, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT3DWithXYZ(m_pTempFFTBuffer, m_lstDims, TRUE);
            _kernelFFTScale3D << <block, threads >> > (m_pMomentumTable, m_pTempFFTBuffer);
            CCLGFFTHelper::FFT3DWithXYZ(m_pTempFFTBuffer, m_lstDims, FALSE);
            _kernelFFTCtoR3D << <block, threads >> > (m_pTempFFTBuffer, m_pGamma22);
#endif
        }

        //======= 4. Gauge Transform    =========
        _kernelCalculateG3D << <block, threads >> > (
            uiT,
            m_pG,
            m_pGamma11,
            m_pGamma12,
            m_pGamma13,
            m_pGamma22,
            m_pGamma23,
            m_fAlpha);
        _kernelGaugeTransform3D << <block, threads >> > (byFieldId, uiT, m_pG, pDeviceBufferPointer);
        _kernelGaugeTransform3DT << <block, threads >> > (byFieldId, uiT, m_pG, pDeviceBufferPointer);
        ++m_iIterate;
    }

    appGeneral(_T("Gauge fixing failed with last error = %2.15f\n"), fTheta);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CGaugeFixingCoulombCornell::CheckRes(const CFieldGauge* pGauge)
#else
Real CGaugeFixingCoulombCornell::CheckRes(const CFieldGauge* pGauge)
#endif
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CGaugeFixingLandauCornell only implemented with gauge SU3!\n"));
        return F(0.0);
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
#if !_CLG_DOUBLEFLOAT
    DOUBLE fRet = 2.0;
#else
    Real fRet = F(0.0);
#endif

    preparethread_S;
    for (SBYTE uiT = 0; uiT < static_cast<SBYTE>(_HC_Lt); ++uiT)
    {
        if (0 == _HC_ALog)
        {
            _kernelCalculateA3D << <block, threads >> > (
                uiT,
                pGaugeSU3->m_pDeviceData,
                m_pA11,
                m_pA12,
                m_pA13,
                m_pA22,
                m_pA23);
        }
        else
        {
            _kernelCalculateA3DLog << <block, threads >> > (
                uiT,
                pGaugeSU3->m_pDeviceData,
                m_pA11,
                m_pA12,
                m_pA13,
                m_pA22,
                m_pA23);
        }

        _kernelCalculateAGradient3D << <block, threads >> > (
            pGauge->m_byFieldId,
            uiT,
            m_pGamma11,
            m_pGamma12,
            m_pGamma13,
            m_pGamma22,
            m_pGamma23,
            m_pA11,
            m_pA12,
            m_pA13,
            m_pA22,
            m_pA23);

        _kernelCalculateTrAGradientSq3D << <block, threads >> > (
            uiT,
            _D_RealThreadBuffer,
            m_pGamma11,
            m_pGamma12,
            m_pGamma13,
            m_pGamma22,
            m_pGamma23);

        fRet += appAbs(appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz) / (3 * _HC_Volume_xyz));
    }
    return fRet / _HC_Lt;
}

CCString CGaugeFixingCoulombCornell::GetInfos(const CCString& tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeFixingCoulombCornell\n");
    sRet = sRet + tab + _T("accuray : ") + appToString(m_fAccuracy) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================