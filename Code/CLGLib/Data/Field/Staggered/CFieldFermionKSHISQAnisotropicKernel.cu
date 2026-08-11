//=============================================================================
// FILENAME : CFieldFermionKSHISQAnisotropicKernel.cu
// 
// DESCRIPTION:
// The device implementations of the anisotropic HISQ (aHISQ) staggered fermion.
// These are the fXiF-weighted (temporal, dir == 3) versions of the KS/HISQ
// kernels, split out of CFieldFermionKSTKernel and CFieldCommonKernel so that
// the shared (performance critical) kernels stay untouched.
//
// REVISION:
//  [mm/dd/yy]
//  [07/26/2026 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "Data/Field/Gauge/CFieldGaugeLink.h"
#include "CFieldFermionKSHISQAnisotropicKernel.h"

//save registor strategy uses global memory frequently, so is slower
//however, it reduces the number of registers, enabling more threads at the one time
#define _CLG_SAVE_REGISTOR_STRATEGY 0

__BEGIN_NAMESPACE

#pragma region DOperator

#pragma region kernel

/**
* Dks = 2am + \sum _{\mu} \eta_{\mu} (n) (U_{\mu}(n) \delta _{n,n+\mu} -U^+_{\mu}(n-\mu) \delta _{n,n-\mu})
* U act on su3
* gamma act on spinor
*
* If bDagger, it is just \eta 5(n) Dks \eta _5(n)
* A easier way is D_{ks}(m=0) is anti-Hermitian, so D^+ = - D_{ks,m=0} + 2am
*
*/
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    Real f2am,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    Real fXiF)
{
    intokernaldir;

    deviceVector result = _makeZero<deviceVector>();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
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
        const deviceGauge& x_Gauge_element = pGauge[linkIndex];
        deviceGauge x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            _dagger(x_m_mu_Gauge_element);
        }

        //U(x,mu) phi(x+ mu)
        deviceVector u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            _mul(u_phi_x_p_m, F(-1.0));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceVector u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            _add(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        else
        {
            _sub(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        _mul(u_phi_x_p_m, eta_mu);
        if (3 == idir)
        {
            _mul(u_phi_x_p_m, fXiF);
        }
        _add(result, u_phi_x_p_m);
    }

    _mul(pResultData[uiSiteIndex], f2am);
    if (bDDagger)
    {
        _sub(pResultData[uiSiteIndex], result);
    }
    else
    {
        _add(pResultData[uiSiteIndex], result);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(pResultData[uiSiteIndex], fCoeff);
        break;
    case EOCT_Complex:
        _mul(pResultData[uiSiteIndex], cCoeff);
        break;
    default:
        break;
    }
}

#pragma region e-o D

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2_MCR(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    Real f2am,
    UBOOL bEven,
    UBOOL bDDagger,
    Real fCoeff,
    Real fXiF)
{
    intokernalEOHalf;

    deviceVector result = pDeviceData[uiSiteIndex];
    _mul(result, f2am);

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _submulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _addmulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _adddagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _subdagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }

#pragma unroll
    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        deviceVector tmpT = _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _sub(result, tmpT);
        }
        else
        {
            _add(result, tmpT);
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        tmpT = _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _add(result, tmpT);
        }
        else
        {
            _sub(result, tmpT);
        }
    }
    _mul(result, fCoeff);
    pDeviceData[uiSiteIndex] = result;
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2_MCC(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    Real f2am,
    UBOOL bEven,
    UBOOL bDDagger,
    CLGComplex cCoeff,
    Real fXiF)
{
    intokernalEOHalf;

    deviceVector result = pDeviceData[uiSiteIndex];
    _mul(result, f2am);

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _submulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _addmulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _adddagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _subdagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }

#pragma unroll
    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        deviceVector tmpT = _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _sub(result, tmpT);
        }
        else
        {
            _add(result, tmpT);
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        tmpT = _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _add(result, tmpT);
        }
        else
        {
            _sub(result, tmpT);
        }
    }
    _mul(result, cCoeff);
    pDeviceData[uiSiteIndex] = result;
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2_M(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    Real f2am,
    UBOOL bEven,
    UBOOL bDDagger,
    Real fXiF)
{
    intokernalEOHalf;

    deviceVector result = pDeviceData[uiSiteIndex];
    _mul(result, f2am);

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _submulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _addmulP(&result, pGauge + linkIndex, pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _adddagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _subdagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }

    #pragma unroll
    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        deviceVector tmpT = _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _sub(result, tmpT);
        }
        else
        {
            _add(result, tmpT);
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        tmpT = _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _add(result, tmpT);
        }
        else
        {
            _sub(result, tmpT);
        }
    }

    pDeviceData[uiSiteIndex] = result;
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2_CR(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bEven,
    UBOOL bDDagger,
    Real fCoeff,
    Real fXiF)
{
    intokernalEOHalf;

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    deviceVector result = _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
        _oppo(result);
    }

    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _adddagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _subdagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }

#pragma unroll
    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        deviceVector tmpT = _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _sub(result, tmpT);
        }
        else
        {
            _add(result, tmpT);
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        tmpT = _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _add(result, tmpT);
        }
        else
        {
            _sub(result, tmpT);
        }
    }

    _mul(result, fCoeff);
    pDeviceData[uiSiteIndex] = result;
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2_CC(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bEven,
    UBOOL bDDagger,
    CLGComplex cCoeff,
    Real fXiF)
{
    intokernalEOHalf;

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    deviceVector result = _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
        _oppo(result);
    }

    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _adddagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _subdagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }

#pragma unroll
    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        deviceVector tmpT = _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _sub(result, tmpT);
        }
        else
        {
            _add(result, tmpT);
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        tmpT = _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _add(result, tmpT);
        }
        else
        {
            _sub(result, tmpT);
        }
    }

    _mul(result, cCoeff);
    pDeviceData[uiSiteIndex] = result;
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST_DOnEvenOrOdd2_NOMC(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bEven,
    UBOOL bDDagger,
    Real fXiF)
{
    intokernalEOHalf;

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    deviceVector result = _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
        _oppo(result);
    }

    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), eta & 1U, bDDagger))
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _adddagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _add(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }
    else
    {
#if _CLG_SAVE_REGISTOR_STRATEGY
        _subdagmulP(&result, pGauge + (x_move_Fermion.m_uiSiteIndex << 2U), pDeviceData + x_move_Fermion.m_uiSiteIndex);
#else
        _sub(result, _dagmulVec(pGauge[x_move_Fermion.m_uiSiteIndex << 2U], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
#endif
    }

#pragma unroll
    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        deviceVector tmpT = _mulVec(pGauge[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _sub(result, tmpT);
        }
        else
        {
            _add(result, tmpT);
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        tmpT = _dagmulVec(pGauge[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR3(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U, bDDagger))
        {
            _add(result, tmpT);
        }
        else
        {
            _sub(result, tmpT);
        }
    }

    pDeviceData[uiSiteIndex] = result;
}

#pragma endregion

/**
 * For some strange boundary condition
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSPlusEtaT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    Real f2am,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    Real fXiF)
{
    intokernaldir;

    deviceVector result = _makeZero<deviceVector>();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //This is in fact, -1 * eta(n + mu)
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        const Real eta_mu2 = (1 == ((pEtaTable[x_m_mu_Fermion.m_uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        const deviceGauge& x_Gauge_element = pGauge[linkIndex];
        deviceGauge x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            _dagger(x_m_mu_Gauge_element);
        }

        //U(x,mu) phi(x+ mu)
        deviceVector u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            _mul(u_phi_x_p_m, F(-1.0) * eta_mu);
        }
        else
        {
            _mul(u_phi_x_p_m, eta_mu);
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceVector u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        _mul(u_dagger_phi_x_m_m, eta_mu2);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            _add(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        else
        {
            _sub(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        //_mul(u_phi_x_p_m, eta_mu);
        if (3 == idir)
        {
            _mul(u_phi_x_p_m, fXiF);
        }
        _add(result, u_phi_x_p_m);
    }

    _mul(pResultData[uiSiteIndex], f2am);
    if (bDDagger)
    {
        _sub(pResultData[uiSiteIndex], result);
    }
    else
    {
        _add(pResultData[uiSiteIndex], result);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(pResultData[uiSiteIndex], fCoeff);
        break;
    case EOCT_Complex:
        _mul(pResultData[uiSiteIndex], cCoeff);
        break;
    default:
        break;
    }
}

/**
 * Calculate Force
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForceT(
    //const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    Real fXiF)
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
            const deviceVector* phi_i = pFermionPointers[uiR];
            const deviceVector* phi_id = pFermionPointers[uiR + uiRational];

            //deviceVector toContract = _mulVec(pGauge[linkIndex], phi_i[x_p_mu_Fermion.m_uiSiteIndex]);
            deviceGauge thisTerm = _makeContract<deviceGauge, deviceVector>(phi_i[x_p_mu_Fermion.m_uiSiteIndex], phi_id[uiSiteIndex]);

            //toContract = _mulVec(pGauge[linkIndex], phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
            _sub(thisTerm, _makeContract<deviceGauge, deviceVector>(phi_id[x_p_mu_Fermion.m_uiSiteIndex], phi_i[uiSiteIndex]));

            if (x_p_mu_Fermion.NeedToOpposite())
            {
                _mul(thisTerm, eta_mu * pNumerators[uiR] * (3 == idir ? fXiF : F(1.0)) * F(-1.0));
            }
            else
            {
                _mul(thisTerm, eta_mu * pNumerators[uiR] * (3 == idir ? fXiF : F(1.0)));
            }

            _sub(pForce[linkIndex], thisTerm);
        }
    }
}

#pragma endregion

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSHISQAnisotropicKernel<deviceVector, deviceGauge, vectorN>::DOperatorKS(UBOOL bEachSiteEta, deviceVector* pTarget, const deviceVector* pSource,
    const deviceGauge* pGauge, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff, Real fXiF)
{
    preparethread;
    if (bEachSiteEta)
    {
        _LAUNCH_KERNEL(_kernelDFermionKSPlusEtaT TMPARG(deviceVector, deviceGauge), block, threads,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            pTarget,
            f2am,
            byFieldId,
            bDagger,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            fXiF);
    }
    else
    {
        _LAUNCH_KERNEL(_kernelDFermionKST TMPARG(deviceVector, deviceGauge), block, threads,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            pTarget,
            f2am,
            byFieldId,
            bDagger,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            fXiF);
    }
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSHISQAnisotropicKernel<deviceVector, deviceGauge, vectorN>::DOperatorKSOnEvenOrOdd(deviceVector* pTargetBuffer,
    const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEvenOdd, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff, Real fXiF)
{
    if (bEvenOdd && abs(f2am) > _CLG_FLT_MIN_)
    {
        appGeneral(_T("Usually D operator on even usually does not have mass, but call with mass=%f\n"), f2am);
    }

    //preparethreadDir;
    //_LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOddDir, block, threads, 
    //    pTargetBuffer,
    //    pGaugeBuffer,
    //    appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pEtaMu,
    //    byFieldId,
    //    f2am,
    //    bEvenOdd,
    //    bDagger,
    //    eOCT,
    //    fRealCoeff,
    //    cCmpCoeff);

    //very slow
    //preparethreadE(8);
    //_LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOddDirPM, block, threads, 
    //    pTargetBuffer,
    //    pGaugeBuffer,
    //    appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pEtaMu,
    //    byFieldId,
    //    f2am,
    //    bEvenOdd ? 1 : 0,
    //    bDagger ? 1 : 0,
    //    eOCT,
    //    fRealCoeff,
    //    cCmpCoeff);


    preparethreadHalf;
    //_LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2, block, threads, 
    //        pTargetBuffer,
    //        pGaugeBuffer,
    //        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //        appGetLattice()->m_pIndexCache->m_pEtaMu,
    //        byFieldId,
    //        f2am,
    //        bEvenOdd ? 1 : 0,
    //        bDagger ? 1 : 0,
    //        eOCT,
    //        fRealCoeff,
    //        cCmpCoeff);

    if (abs(f2am) > _CLG_FLT_MIN_)
    {
        switch (eOCT)
        {
        case EOCT_None:
            _LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2_M TMPARG(deviceVector, deviceGauge), block, threads,
                pTargetBuffer,
                pGaugeBuffer,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                f2am,
                bEvenOdd ? 1 : 0,
                bDagger ? 1 : 0,
                fXiF);
            break;
        case EOCT_Real:
            _LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2_MCR TMPARG(deviceVector, deviceGauge), block, threads,
                pTargetBuffer,
                pGaugeBuffer,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                f2am,
                bEvenOdd ? 1 : 0,
                bDagger ? 1 : 0,
                fRealCoeff,
                fXiF);
            break;
        case EOCT_Complex:
            _LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2_MCC TMPARG(deviceVector, deviceGauge), block, threads,
                pTargetBuffer,
                pGaugeBuffer,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                f2am,
                bEvenOdd ? 1 : 0,
                bDagger ? 1 : 0,
                cCmpCoeff,
                fXiF);
            break;
        default:
            appCrucial(_T("should not call this with eoc:%d\n"), eOCT);
            break;
        }
    }
    else
    {
        switch (eOCT)
        {
        case EOCT_None:
            _LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2_NOMC TMPARG(deviceVector, deviceGauge), block, threads,
                pTargetBuffer,
                pGaugeBuffer,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bEvenOdd ? 1 : 0,
                bDagger ? 1 : 0,
                fXiF);
            break;
        case EOCT_Real:
            _LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2_CR TMPARG(deviceVector, deviceGauge), block, threads,
                pTargetBuffer,
                pGaugeBuffer,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bEvenOdd ? 1 : 0,
                bDagger ? 1 : 0,
                fRealCoeff,
                fXiF);
            break;
        case EOCT_Complex:
            _LAUNCH_KERNEL(_kernelDFermionKST_DOnEvenOrOdd2_CC TMPARG(deviceVector, deviceGauge), block, threads,
                pTargetBuffer,
                pGaugeBuffer,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bEvenOdd ? 1 : 0,
                bDagger ? 1 : 0,
                cCmpCoeff,
                fXiF);
            break;
        default:
            appCrucial(_T("should not call this with eoc:%d\n"), eOCT);
            break;
        }
    }

    //UINT block, threads; 
    //appBlockThreadsE(_HC_VolumeHalf, 8, block, threads);
    //_kernelDFermionKST_DOnEvenOrOddDirPM2 << <block, threads, sizeof(deviceVector)*(threads/8) >> > (
    //    pTargetBuffer,
    //    pGaugeBuffer,
    //    appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pEtaMu,
    //    byFieldId,
    //    f2am,
    //    bEvenOdd ? 1 : 0,
    //    bDagger ? 1 : 0,
    //    eOCT,
    //    fRealCoeff,
    //    cCmpCoeff);

    //if (abs(f2am) > _CLG_FLT_MIN_)
    //{
    //    if (eOCT == EOCT_None)
    //    {
    //        _kernelDFermionKST_DOnEvenOrOddDirPM2_M << <block, threads, sizeof(deviceVector)* (threads / 8) >> > (
    //            pTargetBuffer,
    //            pGaugeBuffer,
    //            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //            appGetLattice()->m_pIndexCache->m_pEtaMu,
    //            f2am,
    //            bEvenOdd ? 1 : 0,
    //            bDagger ? 1 : 0);
    //    }
    //    else
    //    {
    //        _kernelDFermionKST_DOnEvenOrOddDirPM2_MC << <block, threads, sizeof(deviceVector)* (threads / 8) >> > (
    //            pTargetBuffer,
    //            pGaugeBuffer,
    //            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //            appGetLattice()->m_pIndexCache->m_pEtaMu,
    //            f2am,
    //            bEvenOdd ? 1 : 0,
    //            bDagger ? 1 : 0,
    //            eOCT,
    //            fRealCoeff,
    //            cCmpCoeff);
    //    }
    //}
    //else
    //{
    //    if (eOCT == EOCT_None)
    //    {
    //        _kernelDFermionKST_DOnEvenOrOddDirPM2_NOMC << <block, threads, sizeof(deviceVector)* (threads / 8) >> > (
    //            pTargetBuffer,
    //            pGaugeBuffer,
    //            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //            appGetLattice()->m_pIndexCache->m_pEtaMu,
    //            bEvenOdd ? 1 : 0,
    //            bDagger ? 1 : 0);
    //    }
    //    else
    //    {
    //        _kernelDFermionKST_DOnEvenOrOddDirPM2_C << <block, threads, sizeof(deviceVector)* (threads / 8) >> > (
    //            pTargetBuffer,
    //            pGaugeBuffer,
    //            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
    //            appGetLattice()->m_pIndexCache->m_pEtaMu,
    //            bEvenOdd ? 1 : 0,
    //            bDagger ? 1 : 0,
    //            eOCT,
    //            fRealCoeff,
    //            cCmpCoeff);
    //    }
    //}

}

/**
 * partial D_{st0} / partial omega
 * Make sure m_pMDNumerator and m_pRationalFieldPointers are filled
 */
template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSHISQAnisotropicKernel<deviceVector, deviceGauge, vectorN>::DerivateD0(
    const deviceVector* pFermion,
    BYTE byFieldId,
    deviceGauge* pForce,
    const deviceGauge* pGaugeBuffer,
    BYTE byGaugeFieldId,
    const deviceVector* const* pRationalFields,
    const Real* pNumerator,
    UINT uiRationApproxOrder,
    Real fXiF)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionKSForceT TMPARG(deviceVector, deviceGauge), block, threads,
        //pGaugeBuffer,
        pForce,
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pRationalFields,
        pNumerator,
        uiRationApproxOrder,
        byFieldId,
        fXiF);
}

#pragma endregion

#pragma region HISQ

#pragma region kernels

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSNaikCached(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pNaikLink,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bDDagger,
    Real fCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    Real fXiF)
{
    intokernalDir;

    const SIndex& x_p_mu_Fermion = pFermionMove[2 * uiLinkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * uiLinkIndex + 1];

    deviceVector result = _mulVec(pNaikLink[uiLinkIndex], pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
    const BYTE etamu = ((pEtaTable[uiSiteIndex] >> dir) & 1);
    if (_UBOOLXOR(x_p_mu_Fermion.NeedToOpposite(), etamu))
    {
        _oppo(result);
    }
    if (_UBOOLXOR(x_m_mu_Fermion.NeedToOpposite(), etamu))
    {
        _add(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_m_mu_Fermion.m_uiSiteIndex, x_m_mu_Fermion.m_byDir)], pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _sub(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_m_mu_Fermion.m_uiSiteIndex, x_m_mu_Fermion.m_byDir)], pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]));
    }

    if (bDDagger)
    {
        fCoefficient = fCoefficient * F(-1.0);
    }
    _mul(result, fCoefficient);
    if (3 == dir)
    {
        _mul(result, fXiF);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    if (0 == dir)
    {
        _add(pResultData[uiSiteIndex], result);
    }
    __syncthreads();
    if (1 == dir)
    {
        _add(pResultData[uiSiteIndex], result);
    }
    __syncthreads();
    if (2 == dir)
    {
        _add(pResultData[uiSiteIndex], result);
    }
    __syncthreads();
    if (3 == dir)
    {
        _add(pResultData[uiSiteIndex], result);
    }
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSNaikCachedEvenOdd2(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pNaikLink,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bEven,
    UBOOL bDDagger,
    Real fCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    Real fXiF)
{
    intokernalEOHalf;

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    deviceVector result = _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), eta & 1))
    {
        _oppo(result);
    }
    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), eta & 1))
    {
        _add(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _sub(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }

    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        deviceVector tmpT = _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U))
        {
            _sub(result, tmpT);
        }
        else
        {
            _add(result, tmpT);
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U))
        {
            tmpT = _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
            if (3 == idir)
            {
                _mul(tmpT, fXiF);
            }
            _add(result, tmpT);
        }
        else
        {
            tmpT = _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, idir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
            if (3 == idir)
            {
                _mul(tmpT, fXiF);
            }
            _sub(result, tmpT);
        }
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(fCoeff, fCoefficient * (1 - (bDDagger << 1)));
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(cCoeff, fCoefficient * (1 - (bDDagger << 1)));
        _mul(result, cCoeff);
        break;
    default:
        _mul(result, fCoefficient * (1 - (bDDagger << 1)));
        break;
    }

    _add(pDeviceData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSNaikCachedEvenOdd2_R(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pNaikLink,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bEven,
    Real fCoeff,
    Real fXiF)
{
    intokernalEOHalf;

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    deviceVector result = _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), eta & 1))
    {
        _oppo(result);
    }
    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), eta & 1))
    {
        _add(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _sub(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }

    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        deviceVector tmpT = _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U))
        {
            _sub(result, tmpT);
        }
        else
        {
            _add(result, tmpT);
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        tmpT = _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U))
        {
            _add(result, tmpT);
        }
        else
        {
            _sub(result, tmpT);
        }
    }

    _mul(result, fCoeff);
    _add(pDeviceData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSNaikCachedEvenOdd2_C(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pNaikLink,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bEven,
    CLGComplex cCoeff,
    Real fXiF)
{
    intokernalEOHalf;

    UINT linkIndex = (uiSiteIndex << 2U);
    UINT dblinkIndex = (linkIndex << 1U);
    SIndex x_move_Fermion = pFermionMove[dblinkIndex];
    deviceVector result = _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
    if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), eta & 1))
    {
        _oppo(result);
    }
    x_move_Fermion = pFermionMove[dblinkIndex | 1];
    if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), eta & 1))
    {
        _add(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _sub(result, _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]));
    }

    for (BYTE idir = 1U; idir < 4U; ++idir)
    {
        ++linkIndex;
        dblinkIndex = (linkIndex << 1U);
        x_move_Fermion = pFermionMove[dblinkIndex];
        deviceVector tmpT = _mulVec(pNaikLink[linkIndex], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U))
        {
            _sub(result, tmpT);
        }
        else
        {
            _add(result, tmpT);
        }
        x_move_Fermion = pFermionMove[dblinkIndex | 1];
        tmpT = _dagmulVec(pNaikLink[_deviceGetLinkIndex(x_move_Fermion.m_uiSiteIndex, x_move_Fermion.m_byDir)], pDeviceData[x_move_Fermion.m_uiSiteIndex]);
        if (3 == idir)
        {
            _mul(tmpT, fXiF);
        }
        if (_UBOOLXOR(x_move_Fermion.NeedToOpposite(), (eta >> idir) & 1U))
        {
            _add(result, tmpT);
        }
        else
        {
            _sub(result, tmpT);
        }
    }

    _mul(result, cCoeff);

    _add(pDeviceData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddConnectionNaik(
    Real fCoeff,
    const deviceVector* __restrict__ n,
    const BYTE* __restrict__ pEtaTable,
    const SIndex* __restrict__ pFermionMove,
    deviceGauge* res,
    BYTE byFieldId,
    Real fXiF)
{
    intokernalDir;

    SIndex sn1 = pFermionMove[uiLinkIndex << 1];
    deviceGauge resv = _makeContract<deviceGauge, deviceVector>(n[sn1.m_uiSiteIndex], n[uiSiteIndex]);
    _mul(resv, fCoeff);
    if (3 == dir)
    {
        _mul(resv, fXiF);
    }
    BYTE eta = pEtaTable[uiSiteIndex] >> dir;

    //This is not eta, but even-odd
    if (pEtaTable[uiSiteIndex] >> 4)
    {
        eta = eta + 1;
    }
    if (sn1.NeedToOpposite())
    {
        eta = eta + 1;
    }

    if (eta & 1)
    {
        _sub(res[uiLinkIndex], resv);
    }
    else
    {
        _add(res[uiLinkIndex], resv);
    }
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddConnectionNaik2(
    Real fCoeff,
    const deviceVector* __restrict__ n,
    const deviceVector* __restrict__ n_p_mu,
    const BYTE* __restrict__ pEtaTable,
    const SIndex* __restrict__ pFermionMove,
    deviceGauge* res,
    BYTE byFieldId,
    Real fXiF)
{
    intokernalDir;

    SIndex sn1 = pFermionMove[uiLinkIndex << 1];
    deviceGauge resv = _makeContract<deviceGauge, deviceVector>(n_p_mu[sn1.m_uiSiteIndex], n[uiSiteIndex]);
    _sub(resv, _makeContract<deviceGauge, deviceVector>(n[sn1.m_uiSiteIndex], n_p_mu[uiSiteIndex]));
    _mul(resv, fCoeff);
    if (3 == dir)
    {
        _mul(resv, fXiF);
    }
    BYTE eta = pEtaTable[uiSiteIndex] >> dir;

    if (sn1.NeedToOpposite())
    {
        eta = eta + 1;
    }

    if (eta & 1)
    {
        _sub(res[uiLinkIndex], resv);
    }
    else
    {
        _add(res[uiLinkIndex], resv);
    }
}

#pragma endregion

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSHISQAnisotropicKernel<deviceVector, deviceGauge, vectorN>::DOperatorNaik(UBOOL bEachSiteEta, deviceVector* pTargetBuffer, const deviceVector* pBuffer,
    BYTE byFieldId, BYTE byGaugeFieldId, Real fNaik, Real fEpsilon,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff, Real fXiF)
{
    const deviceGauge* pNaikLink = (const deviceGauge*)(appGetGaugeSmearing(byGaugeFieldId)->GetNaikLink()->GetData());
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelDFermionKSNaikCached TMPARG(deviceVector, deviceGauge), block, threads,
        (const deviceVector*)pBuffer,
        pNaikLink,
        appGetLattice()->m_pIndexCache->m_pNaikCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        (deviceVector*)pTargetBuffer,
        byFieldId,
        byGaugeFieldId,
        bDagger,
        fNaik,
        eOCT,
        fRealCoeff,
        cCmpCoeff,
        fXiF
    );
    if (abs(fEpsilon) > _CLG_FLT_EPSILON)
    {
        const deviceGauge* pEffectiveGaugeLevel1 = (const deviceGauge*)(appGetGaugeSmearing(byGaugeFieldId)->GetEffectiveGaugeLevel1()->GetData());
        _LAUNCH_KERNEL(_kernelDFermionKSNaikCached TMPARG(deviceVector, deviceGauge), block, threads,
            (const deviceVector*)pBuffer,
            pEffectiveGaugeLevel1,
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            (deviceVector*)pTargetBuffer,
            byFieldId,
            byGaugeFieldId,
            bDagger,
            fEpsilon,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            fXiF
        );
    }
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSHISQAnisotropicKernel<deviceVector, deviceGauge, vectorN>::DOperatorNaikOnEvenOrOdd(deviceVector* pTargetBuffer,
    BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEven, Real fNaik, Real fEpsilon,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff, Real fXiF)
{
    const deviceGauge* pNaikLink = (const deviceGauge*)(appGetGaugeSmearing(byGaugeFieldId)->GetNaikLink()->GetData());
    //preparethreadDir;
    //_LAUNCH_KERNEL(_kernelDFermionKSNaikCachedEvenOdd, block, threads, 
    //    (deviceSU3Vector*)pTargetBuffer,
    //    pNaikLink,
    //    appGetLattice()->m_pIndexCache->m_pNaikCache[m_byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pEtaMu,
    //    m_byFieldId,
    //    byGaugeFieldId,
    //    bEven,
    //    bDagger,
    //    m_fNaik,
    //    eOCT,
    //    fRealCoeff,
    //    cCmpCoeff
    //    );

    preparethreadHalf;
    if (EOCT_None == eOCT || EOCT_Real == eOCT)
    {

        Real fCoef = bDagger ? (-fNaik) : (fNaik);
        if (EOCT_Real == eOCT)
        {
            fCoef *= fRealCoeff;
        }
        _LAUNCH_KERNEL(_kernelDFermionKSNaikCachedEvenOdd2_R TMPARG(deviceVector, deviceGauge), block, threads,
            (deviceVector*)pTargetBuffer,
            pNaikLink,
            appGetLattice()->m_pIndexCache->m_pNaikCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            byFieldId,
            byGaugeFieldId,
            bEven,
            fCoef,
            fXiF
        );
    }
    else if (EOCT_Complex == eOCT)
    {
        CLGComplex cCoeff = cuCmulf_cr(cCmpCoeff, bDagger ? -fNaik : fNaik);
        _LAUNCH_KERNEL(_kernelDFermionKSNaikCachedEvenOdd2_C TMPARG(deviceVector, deviceGauge), block, threads,
            (deviceVector*)pTargetBuffer,
            pNaikLink,
            appGetLattice()->m_pIndexCache->m_pNaikCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            byFieldId,
            byGaugeFieldId,
            bEven,
            cCoeff,
            fXiF
        );
    }
    if (abs(fEpsilon) > _CLG_FLT_EPSILON)
    {
        const deviceGauge* pEffectiveGaugeLevel1 = (const deviceGauge*)(appGetGaugeSmearing(byGaugeFieldId)->GetEffectiveGaugeLevel1()->GetData());
        if (EOCT_None == eOCT || EOCT_Real == eOCT)
        {
            Real fCoef = bDagger ? (-fEpsilon) : (fEpsilon);
            if (EOCT_Real == eOCT)
            {
                fCoef *= fRealCoeff;
            }
            _LAUNCH_KERNEL(_kernelDFermionKSNaikCachedEvenOdd2_R TMPARG(deviceVector, deviceGauge), block, threads,
                (deviceVector*)pTargetBuffer,
                pEffectiveGaugeLevel1,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                byFieldId,
                byGaugeFieldId,
                bEven,
                fCoef,
                fXiF
            );
        }
        else if (EOCT_Complex == eOCT)
        {
            CLGComplex cCoeff = cuCmulf_cr(cCmpCoeff, bDagger ? -fEpsilon : fEpsilon);
            _LAUNCH_KERNEL(_kernelDFermionKSNaikCachedEvenOdd2_C TMPARG(deviceVector, deviceGauge), block, threads,
                (deviceVector*)pTargetBuffer,
                pEffectiveGaugeLevel1,
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                byFieldId,
                byGaugeFieldId,
                bEven,
                cCoeff,
                fXiF
            );
        }
    }
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSHISQAnisotropicKernel<deviceVector, deviceGauge, vectorN>::NaikConnection(Real fNumerator, const deviceVector* rfield, deviceGauge* naikforce, BYTE byFieldId, Real fXiF)
{
    _RECORD(CFieldFermionKSHISQAnisotropicKernel::NaikConnection);
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelAddConnectionNaik TMPARG(deviceVector, deviceGauge), block, threads,
        fNumerator,
        rfield,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        appGetLattice()->m_pIndexCache->m_pNaikCache[byFieldId],
        naikforce,
        byFieldId,
        fXiF
    );
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSHISQAnisotropicKernel<deviceVector, deviceGauge, vectorN>::NaikConnection2(Real fNumerator, const deviceVector* rfield, const deviceVector* rfield_p_mu, deviceGauge* naikforce, BYTE byFieldId, Real fXiF)
{
    _RECORD(CFieldFermionKSHISQAnisotropicKernel::NaikConnection);
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelAddConnectionNaik2 TMPARG(deviceVector, deviceGauge), block, threads,
        fNumerator,
        rfield,
        rfield_p_mu,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        appGetLattice()->m_pIndexCache->m_pNaikCache[byFieldId],
        naikforce,
        byFieldId,
        fXiF
    );
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSHISQAnisotropicKernel<deviceVector, deviceGauge, vectorN>::AddConnection2(Real fNumerator, const deviceVector* rfield, const deviceVector* rfield_p_mu, deviceGauge* naikforce, BYTE byFieldId, Real fXiF)
{
    _RECORD(CFieldFermionKSHISQAnisotropicKernel::NaikConnection);
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelAddConnectionNaik2 TMPARG(deviceVector, deviceGauge), block, threads,
        fNumerator,
        rfield,
        rfield_p_mu,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        naikforce,
        byFieldId,
        fXiF
    );
}

#pragma endregion

#pragma region Connection

#pragma region kernel

template<typename T, typename gaugetype>
__global__ void _CLG_LAUNCH_BOUND
_kernelConnectionOneFieldStaggered(
    Real fCoeff,
    const T* __restrict__ n, 
    const BYTE* __restrict__ pEtaTable,
    gaugetype* res, 
    const SIndex* __restrict__ move,
    Real fXiF)
{
    intokernalDirInt4;
    
    const SIndex& x_p_mu_Fermion = move[2 * uiLinkIndex];
    res[uiLinkIndex] = _makeContract<gaugetype, T>(n[x_p_mu_Fermion.m_uiSiteIndex], n[uiSiteIndex]);
    _mul(res[uiLinkIndex], fCoeff);
    if (3 == dir)
    {
        _mul(res[uiLinkIndex], fXiF);
    }
    BYTE eta = pEtaTable[uiSiteIndex] >> dir;
    if (sSite4.IsOdd())
    {
        eta = eta + 1;
    }
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        eta = eta + 1;
    }

    if (eta & 1)
    {
        _oppo(res[uiLinkIndex]);
    }
}

template<typename T, typename gaugetype>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddConnectionOneFieldStaggered(
    Real fCoeff,
    const T* __restrict__ n,
    const BYTE* __restrict__ pEtaTable,
    gaugetype* res,
    const SIndex* __restrict__ move,
    Real fXiF)
{
    intokernalDirInt4;

    const SIndex& x_p_mu_Fermion = move[2 * uiLinkIndex];
    BYTE eta = pEtaTable[uiSiteIndex] >> dir;
    if (sSite4.IsOdd())
    {
        eta = eta + 1;
    }
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        eta = eta + 1;
    }
    gaugetype toadd = _makeContract<gaugetype, T>(n[x_p_mu_Fermion.m_uiSiteIndex], n[uiSiteIndex]);
    _mul(toadd, fCoeff);
    if (3 == dir)
    {
        _mul(toadd, fXiF);
    }
    if (eta & 1)
    {
        _sub(res[uiLinkIndex], toadd);
    }
    else
    {
        _add(res[uiLinkIndex], toadd);
    }
}

#pragma endregion

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSHISQAnisotropicKernel<deviceVector, deviceGauge, vectorN>::ConnectionOneFieldStaggered(const deviceVector* v, deviceGauge* res, BYTE byFieldId, Real fCoeff, Real fXiF)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelConnectionOneFieldStaggered TMPARG(deviceVector, deviceGauge), block, threads,
        fCoeff,
        v, 
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        res, 
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        fXiF);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSHISQAnisotropicKernel<deviceVector, deviceGauge, vectorN>::AddConnectionOneFieldStaggered(const deviceVector* v, deviceGauge* res, BYTE byFieldId, Real fCoeff, Real fXiF)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelAddConnectionOneFieldStaggered TMPARG(deviceVector, deviceGauge), block, threads,
        fCoeff,
        v,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        res,
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        fXiF);
}

#pragma endregion

template class CFieldFermionKSHISQAnisotropicKernel<CLGComplex, CLGComplex, 1>;
template class CFieldFermionKSHISQAnisotropicKernel<deviceSU2Vector, deviceSU2, 2>;
template class CFieldFermionKSHISQAnisotropicKernel<deviceSU3Vector, deviceSU3, 3>;
template class CFieldFermionKSHISQAnisotropicKernel<deviceSU4Vector, deviceSU4, 4>;
//template class CFieldFermionKSHISQAnisotropicKernel<deviceSU5Vector, deviceSU5, 5>;
//template class CFieldFermionKSHISQAnisotropicKernel<deviceSU6Vector, deviceSU6, 6>;
//template class CFieldFermionKSHISQAnisotropicKernel<deviceSU7Vector, deviceSU7, 7>;
//template class CFieldFermionKSHISQAnisotropicKernel<deviceSU8Vector, deviceSU8, 8>;

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
