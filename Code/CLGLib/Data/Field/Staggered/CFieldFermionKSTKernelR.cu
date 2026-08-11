//=============================================================================
// FILENAME : CFieldFermionKSTKernel.cu
// 
// DESCRIPTION:
// This is the device implementations of Wilson fermion
//
// This implementation assumes SU3 and square lattice
//
// REVISION:
//  [mm/dd/yy]
//  [07/21/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "Data/Field/Gauge/CFieldGaugeLink.h"
#include "CFieldFermionKST.h"

__BEGIN_NAMESPACE

#pragma region DOperator Rotation

#pragma region kernel

/**
* When link n and n+mu, the coordinate is stick with n
* When link n and n-mu, the coordinate is stick with n-mu
* Irrelavent with tau
* Optimization: bXorY removed, block.x *= 2
*/
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTermT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    SSmallInt4 sCenter,
    Real shift,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;
    //if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sSite4)].IsDirichlet())
    //{
    //    return;
    //}
    deviceVector result = _makeZero<deviceVector>();
    //const INT eta_tau = ((pEtaTable[uiSiteIndex] >> 3) & 1);
    const INT eta_tau = pEtaTable[uiSiteIndex] >> 3;

    #pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusMu = idx & 2;
        const UBOOL bPlusTau = idx & 4;
        //x or y, and y or x is the derivate, not coefficient
        const UINT bXorY = idx & 1;
        const UINT bYorX = 1 - bXorY;
        SSmallInt4 sTargetSite = sSite4;
        SSmallInt4 sMidSite = sSite4;
        sTargetSite.m_byData4[bYorX] = sTargetSite.m_byData4[bYorX] + (bPlusMu ? 2 : -2);
        sMidSite.m_byData4[bYorX] = sMidSite.m_byData4[bYorX] + (bPlusMu ? 1 : -1);
        sTargetSite.w = sTargetSite.w + (bPlusTau ? 1 : -1);
        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sTargetSite)];
        if (sTargetBigIndex.IsDirichlet())
        {
            continue;
        }
        const SIndex& sMiddleBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sMidSite)];
        sMidSite = __deviceSiteIndexToInt4(sMiddleBigIndex.m_uiSiteIndex);

        //note that bYorX = 1, it is x partial_y term, therefore is '-'
        //INT this_eta_tau = (bPlusTau ? eta_tau : ((pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3) & 1))
        INT this_eta_tau = (bPlusTau ? eta_tau : (pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3))
            + bYorX;

        if (sTargetBigIndex.NeedToOpposite())
        {
            this_eta_tau = this_eta_tau + 1;
        }

        deviceVector right = _mulVec(_deviceVXXTauOptimizedT(pGauge, sSite4, byGaugeFieldId, bXorY, bPlusMu, bPlusTau),
            pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        //when bXorY = 1, it is y partial _x, so is [1]
        //when bXorY = 0, it is x partial _y, so is [0]
        _mul(right, sMidSite.m_byData4[bXorY] - sCenter.m_byData4[bXorY] + shift);

        if (!bPlusMu)
        {
            //for -2x, -2y terms, there is another minus sign
            this_eta_tau = this_eta_tau + 1;
        }

        if (this_eta_tau & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, F(-0.25) * fOmega);
    }
    else
    {
        _mul(result, F(0.25) * fOmega);
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

    _add(pResultData[uiSiteIndex], result);
}

/**
 * Even-odd version of _kernelDFermionKS_PR_XYTermT
 * in-place: read opposite parity sites, add to the sites of one parity
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTermT_EO(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    SSmallInt4 sCenter,
    Real shift,
    UBOOL bEven,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalEOHalf;
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    deviceVector result = _makeZero<deviceVector>();
    //const INT eta_tau = ((pEtaTable[uiSiteIndex] >> 3) & 1);
    const INT eta_tau = eta >> 3;

    #pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusMu = idx & 2;
        const UBOOL bPlusTau = idx & 4;
        //x or y, and y or x is the derivate, not coefficient
        const UINT bXorY = idx & 1;
        const UINT bYorX = 1 - bXorY;
        SSmallInt4 sTargetSite = sSite4;
        SSmallInt4 sMidSite = sSite4;
        sTargetSite.m_byData4[bYorX] = sTargetSite.m_byData4[bYorX] + (bPlusMu ? 2 : -2);
        sMidSite.m_byData4[bYorX] = sMidSite.m_byData4[bYorX] + (bPlusMu ? 1 : -1);
        sTargetSite.w = sTargetSite.w + (bPlusTau ? 1 : -1);
        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sTargetSite)];
        if (sTargetBigIndex.IsDirichlet())
        {
            continue;
        }
        const SIndex& sMiddleBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sMidSite)];
        sMidSite = __deviceSiteIndexToInt4(sMiddleBigIndex.m_uiSiteIndex);

        //note that bYorX = 1, it is x partial_y term, therefore is '-'
        //INT this_eta_tau = (bPlusTau ? eta_tau : ((pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3) & 1))
        INT this_eta_tau = (bPlusTau ? eta_tau : (pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3))
            + bYorX;

        if (sTargetBigIndex.NeedToOpposite())
        {
            this_eta_tau = this_eta_tau + 1;
        }

        deviceVector right = _mulVec(_deviceVXXTauOptimizedT(pGauge, sSite4, byGaugeFieldId, bXorY, bPlusMu, bPlusTau),
            pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        //when bXorY = 1, it is y partial _x, so is [1]
        //when bXorY = 0, it is x partial _y, so is [0]
        _mul(right, sMidSite.m_byData4[bXorY] - sCenter.m_byData4[bXorY] + shift);

        if (!bPlusMu)
        {
            //for -2x, -2y terms, there is another minus sign
            this_eta_tau = this_eta_tau + 1;
        }

        if (this_eta_tau & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, F(-0.25) * fOmega);
    }
    else
    {
        _mul(result, F(0.25) * fOmega);
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

    _add(pDeviceData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_Cached_XYTermT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pCachedGauge,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    SSmallInt4 sCenter,
    Real shift,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;
    //if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sSite4)].IsDirichlet())
    //{
    //    return;
    //}
    deviceVector result = _makeZero<deviceVector>();
    //const INT eta_tau = ((pEtaTable[uiSiteIndex] >> 3) & 1);
    const INT eta_tau = pEtaTable[uiSiteIndex] >> 3;

#pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusMu = idx & 2;
        const UBOOL bPlusTau = idx & 4;
        //x or y, and y or x is the derivate, not coefficient
        const UINT bXorY = idx & 1;
        const UINT bYorX = 1 - bXorY;
        SSmallInt4 sTargetSite = sSite4;
        SSmallInt4 sMidSite = sSite4;
        sTargetSite.m_byData4[bYorX] = sTargetSite.m_byData4[bYorX] + (bPlusMu ? 2 : -2);
        sMidSite.m_byData4[bYorX] = sMidSite.m_byData4[bYorX] + (bPlusMu ? 1 : -1);
        sTargetSite.w = sTargetSite.w + (bPlusTau ? 1 : -1);
        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sTargetSite)];
        if (sTargetBigIndex.IsDirichlet())
        {
            continue;
        }
        const SIndex& sMiddleBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sMidSite)];
        sMidSite = __deviceSiteIndexToInt4(sMiddleBigIndex.m_uiSiteIndex);

        //note that bYorX = 1, it is x partial_y term, therefore is '-'
        //INT this_eta_tau = (bPlusTau ? eta_tau : ((pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3) & 1))
        INT this_eta_tau = (bPlusTau ? eta_tau : (pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3))
            + bYorX;

        if (sTargetBigIndex.NeedToOpposite())
        {
            this_eta_tau = this_eta_tau + 1;
        }

        deviceVector right = _mulVec(pCachedGauge[uiSiteIndex * 8 + idx], pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        //when bXorY = 1, it is y partial _x, so is [1]
        //when bXorY = 0, it is x partial _y, so is [0]
        _mul(right, sMidSite.m_byData4[bXorY] - sCenter.m_byData4[bXorY] + shift);

        if (!bPlusMu)
        {
            //for -2x, -2y terms, there is another minus sign
            this_eta_tau = this_eta_tau + 1;
        }

        if (this_eta_tau & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, F(-0.25) * fOmega);
    }
    else
    {
        _mul(result, F(0.25) * fOmega);
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

    _add(pResultData[uiSiteIndex], result);
}

/**
 * Even-odd version of _kernelDFermionKS_PR_Cached_XYTermT
 * in-place: read opposite parity sites, add to the sites of one parity
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_Cached_XYTermT_EO(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pCachedGauge,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    SSmallInt4 sCenter,
    Real shift,
    UBOOL bEven,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalEOHalf;
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    deviceVector result = _makeZero<deviceVector>();
    //const INT eta_tau = ((pEtaTable[uiSiteIndex] >> 3) & 1);
    const INT eta_tau = eta >> 3;

#pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusMu = idx & 2;
        const UBOOL bPlusTau = idx & 4;
        //x or y, and y or x is the derivate, not coefficient
        const UINT bXorY = idx & 1;
        const UINT bYorX = 1 - bXorY;
        SSmallInt4 sTargetSite = sSite4;
        SSmallInt4 sMidSite = sSite4;
        sTargetSite.m_byData4[bYorX] = sTargetSite.m_byData4[bYorX] + (bPlusMu ? 2 : -2);
        sMidSite.m_byData4[bYorX] = sMidSite.m_byData4[bYorX] + (bPlusMu ? 1 : -1);
        sTargetSite.w = sTargetSite.w + (bPlusTau ? 1 : -1);
        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sTargetSite)];
        if (sTargetBigIndex.IsDirichlet())
        {
            continue;
        }
        const SIndex& sMiddleBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sMidSite)];
        sMidSite = __deviceSiteIndexToInt4(sMiddleBigIndex.m_uiSiteIndex);

        //note that bYorX = 1, it is x partial_y term, therefore is '-'
        //INT this_eta_tau = (bPlusTau ? eta_tau : ((pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3) & 1))
        INT this_eta_tau = (bPlusTau ? eta_tau : (pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3))
            + bYorX;

        if (sTargetBigIndex.NeedToOpposite())
        {
            this_eta_tau = this_eta_tau + 1;
        }

        deviceVector right = _mulVec(pCachedGauge[uiSiteIndex * 8 + idx], pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        //when bXorY = 1, it is y partial _x, so is [1]
        //when bXorY = 0, it is x partial _y, so is [0]
        _mul(right, sMidSite.m_byData4[bXorY] - sCenter.m_byData4[bXorY] + shift);

        if (!bPlusMu)
        {
            //for -2x, -2y terms, there is another minus sign
            this_eta_tau = this_eta_tau + 1;
        }

        if (this_eta_tau & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, F(-0.25) * fOmega);
    }
    else
    {
        _mul(result, F(0.25) * fOmega);
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

    _add(pDeviceData[uiSiteIndex], result);
}

/**
* almost copy of _kernelDFermionKS_PR_XYTermT
*/
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTermEMT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pPhase,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    Real fCharge,
    SSmallInt4 sCenter,
    Real shift,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();
    //const INT eta_tau = ((pEtaTable[uiSiteIndex] >> 3) & 1);
    const INT eta_tau = pEtaTable[uiSiteIndex] >> 3;

    #pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusMu = idx & 2;
        const UBOOL bPlusTau = idx & 4;
        //x or y, and y or x is the derivate, not coefficient
        const UINT bXorY = idx & 1;
        const UINT bYorX = 1 - bXorY;
        SSmallInt4 sTargetSite = sSite4;
        SSmallInt4 sMidSite = sSite4;
        sTargetSite.m_byData4[bYorX] = sTargetSite.m_byData4[bYorX] + (bPlusMu ? 2 : -2);
        sMidSite.m_byData4[bYorX] = sMidSite.m_byData4[bYorX] + (bPlusMu ? 1 : -1);
        sTargetSite.w = sTargetSite.w + (bPlusTau ? 1 : -1);
        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sTargetSite)];
        const SIndex& sMiddleBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sMidSite)];
        sMidSite = __deviceSiteIndexToInt4(sMiddleBigIndex.m_uiSiteIndex);

        //note that bYorX = 1, it is x partial_y term, therefore is '-'
        //INT this_eta_tau = (bPlusTau ? eta_tau : ((pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3) & 1))
        INT this_eta_tau = (bPlusTau ? eta_tau : (pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3))
            + bYorX;

        if (sTargetBigIndex.NeedToOpposite())
        {
            this_eta_tau = this_eta_tau + 1;
        }

        deviceVector right = _mulVec(_deviceVXXTauOptimizedEMT(pGauge, pPhase, sSite4, fCharge, byGaugeFieldId, bXorY, bPlusMu, bPlusTau),
            pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        //when bXorY = 1, it is y partial _x, so is [1]
        //when bXorY = 0, it is x partial _y, so is [0]
        _mul(right, sMidSite.m_byData4[bXorY] - sCenter.m_byData4[bXorY] + shift);

        if (!bPlusMu)
        {
            //for -2x, -2y terms, there is another minus sign
            this_eta_tau = this_eta_tau + 1;
        }

        if (this_eta_tau & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, F(-0.25) * fOmega);
    }
    else
    {
        _mul(result, F(0.25) * fOmega);
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

    _add(pResultData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTermRealT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    SSmallInt4 sCenter,
    Real shift,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();
    //const INT eta_tau = ((pEtaTable[uiSiteIndex] >> 3) & 1);
    const INT eta_tau = pEtaTable[uiSiteIndex] >> 3;

#pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusMu = idx & 2;
        const UBOOL bPlusTau = idx & 4;
        //x or y, and y or x is the derivate, not coefficient
        const UINT bXorY = idx & 1;
        const UINT bYorX = 1 - bXorY;
        SSmallInt4 sTargetSite = sSite4;
        SSmallInt4 sMidSite = sSite4;
        sTargetSite.m_byData4[bYorX] = sTargetSite.m_byData4[bYorX] + (bPlusMu ? 2 : -2);
        sMidSite.m_byData4[bYorX] = sMidSite.m_byData4[bYorX] + (bPlusMu ? 1 : -1);
        sTargetSite.w = sTargetSite.w + (bPlusTau ? 1 : -1);
        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sTargetSite)];
        const SIndex& sMiddleBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sMidSite)];
        sMidSite = __deviceSiteIndexToInt4(sMiddleBigIndex.m_uiSiteIndex);

        //note that bYorX = 1, it is x partial_y term, therefore is '-'
        //INT this_eta_tau = (bPlusTau ? eta_tau : ((pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3) & 1))
        INT this_eta_tau = (bPlusTau ? eta_tau : (pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3))
            + bYorX;

        if (sTargetBigIndex.NeedToOpposite())
        {
            this_eta_tau = this_eta_tau + 1;
        }

        deviceVector right = _mulVec(_deviceVXXTauOptimizedT(pGauge, sSite4, byGaugeFieldId, bXorY, bPlusMu, bPlusTau),
            pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        //when bXorY = 1, it is y partial _x, so is [1]
        //when bXorY = 0, it is x partial _y, so is [0]
        _mul(right, sMidSite.m_byData4[bXorY] - sCenter.m_byData4[bXorY] + shift);

        if (!bPlusMu)
        {
            //for -2x, -2y terms, there is another minus sign
            this_eta_tau = this_eta_tau + 1;
        }

        if (this_eta_tau & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    //if (bDDagger)
    //{
    //    _mul(result, F(-0.25) * fOmega);
    //}
    //else
    //{
        //_mul(result, F(0.25) * fOmega);
    //}
    _mul(result, _make_cuComplex(F(0.0), F(-0.25) * fOmega));

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

    _add(pResultData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTau_TermT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

#ifndef _CLG_DTK
    #pragma unroll
#endif
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusX = (0 != (idx & 1));
        const UBOOL bPlusY = (0 != (idx & 2));
        const UBOOL bPlusT = (0 != (idx & 4));

        SSmallInt4 sOffset = sSite4;
        sOffset.x = sOffset.x + (bPlusX ? 1 : -1);
        sOffset.y = sOffset.y + (bPlusY ? 1 : -1);
        sOffset.w = sOffset.w + (bPlusT ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];
        if (sTargetBigIndex.IsDirichlet())
        {
            continue;
        }

        const deviceVector right = _mulVec(_deviceVXYTOptimizedT(pGauge, sSite4, byGaugeFieldId, bPlusX, bPlusY, bPlusT), pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta124 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta124 = bPlusT ? _deviceEta3(sSite4, 2) : (_deviceEta3(site_target, 2) + 1);

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta124 = eta124 + 1;
        }

        if (eta124 & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, -F(0.125) * fOmega);
    }
    else
    {
        _mul(result, F(0.125) * fOmega);
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

    _add(pResultData[uiSiteIndex], result);
}

/**
 * Even-odd version of _kernelDFermionKS_PR_XYTau_TermT
 * in-place: read opposite parity sites, add to the sites of one parity
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTau_TermT_EO(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    UBOOL bEven,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalEOHalf;
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    deviceVector result = _makeZero<deviceVector>();

#ifndef _CLG_DTK
    #pragma unroll
#endif
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusX = (0 != (idx & 1));
        const UBOOL bPlusY = (0 != (idx & 2));
        const UBOOL bPlusT = (0 != (idx & 4));

        SSmallInt4 sOffset = sSite4;
        sOffset.x = sOffset.x + (bPlusX ? 1 : -1);
        sOffset.y = sOffset.y + (bPlusY ? 1 : -1);
        sOffset.w = sOffset.w + (bPlusT ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];
        if (sTargetBigIndex.IsDirichlet())
        {
            continue;
        }

        const deviceVector right = _mulVec(_deviceVXYTOptimizedT(pGauge, sSite4, byGaugeFieldId, bPlusX, bPlusY, bPlusT), pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta124 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta124 = bPlusT ? _deviceEta3(sSite4, 2) : (_deviceEta3(site_target, 2) + 1);

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta124 = eta124 + 1;
        }

        if (eta124 & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, -F(0.125) * fOmega);
    }
    else
    {
        _mul(result, F(0.125) * fOmega);
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

    _add(pDeviceData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_Cached_XYTau_TermT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pCachedGauge,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

#ifndef _CLG_DTK
#pragma unroll
#endif
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusX = (0 != (idx & 1));
        const UBOOL bPlusY = (0 != (idx & 2));
        const UBOOL bPlusT = (0 != (idx & 4));

        SSmallInt4 sOffset = sSite4;
        sOffset.x = sOffset.x + (bPlusX ? 1 : -1);
        sOffset.y = sOffset.y + (bPlusY ? 1 : -1);
        sOffset.w = sOffset.w + (bPlusT ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];
        if (sTargetBigIndex.IsDirichlet())
        {
            continue;
        }

        const deviceVector right = _mulVec(pCachedGauge[uiSiteIndex * 8 + idx], pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta124 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta124 = bPlusT ? _deviceEta3(sSite4, 2) : (_deviceEta3(site_target, 2) + 1);

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta124 = eta124 + 1;
        }

        if (eta124 & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, -F(0.125) * fOmega);
    }
    else
    {
        _mul(result, F(0.125) * fOmega);
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

    _add(pResultData[uiSiteIndex], result);
}

/**
 * Even-odd version of _kernelDFermionKS_PR_Cached_XYTau_TermT
 * in-place: read opposite parity sites, add to the sites of one parity
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_Cached_XYTau_TermT_EO(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pCachedGauge,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    UBOOL bEven,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalEOHalf;
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    deviceVector result = _makeZero<deviceVector>();

#ifndef _CLG_DTK
#pragma unroll
#endif
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusX = (0 != (idx & 1));
        const UBOOL bPlusY = (0 != (idx & 2));
        const UBOOL bPlusT = (0 != (idx & 4));

        SSmallInt4 sOffset = sSite4;
        sOffset.x = sOffset.x + (bPlusX ? 1 : -1);
        sOffset.y = sOffset.y + (bPlusY ? 1 : -1);
        sOffset.w = sOffset.w + (bPlusT ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];
        if (sTargetBigIndex.IsDirichlet())
        {
            continue;
        }

        const deviceVector right = _mulVec(pCachedGauge[uiSiteIndex * 8 + idx], pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta124 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta124 = bPlusT ? _deviceEta3(sSite4, 2) : (_deviceEta3(site_target, 2) + 1);

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta124 = eta124 + 1;
        }

        if (eta124 & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, -F(0.125) * fOmega);
    }
    else
    {
        _mul(result, F(0.125) * fOmega);
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

    _add(pDeviceData[uiSiteIndex], result);
}


template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTau_TermEMT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pPhase,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    Real fCharge,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

#ifndef _CLG_DTK
#pragma unroll
#endif
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusX = (0 != (idx & 1));
        const UBOOL bPlusY = (0 != (idx & 2));
        const UBOOL bPlusT = (0 != (idx & 4));

        SSmallInt4 sOffset = sSite4;
        sOffset.x = sOffset.x + (bPlusX ? 1 : -1);
        sOffset.y = sOffset.y + (bPlusY ? 1 : -1);
        sOffset.w = sOffset.w + (bPlusT ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceVector right = _mulVec(_deviceVXYTOptimizedEMT(pGauge, pPhase, sSite4, fCharge, byGaugeFieldId, bPlusX, bPlusY, bPlusT), pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta124 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta124 = bPlusT ? _deviceEta3(sSite4, 2) : (_deviceEta3(site_target, 2) + 1);

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta124 = eta124 + 1;
        }

        if (eta124 & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, -F(0.125) * fOmega);
    }
    else
    {
        _mul(result, F(0.125) * fOmega);
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

    _add(pResultData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTau_TermRealT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

#ifndef _CLG_DTK
    #pragma unroll
#endif
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusX = (0 != (idx & 1));
        const UBOOL bPlusY = (0 != (idx & 2));
        const UBOOL bPlusT = (0 != (idx & 4));

        SSmallInt4 sOffset = sSite4;
        sOffset.x = sOffset.x + (bPlusX ? 1 : -1);
        sOffset.y = sOffset.y + (bPlusY ? 1 : -1);
        sOffset.w = sOffset.w + (bPlusT ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceVector right = _mulVec(_deviceVXYTOptimizedT(pGauge, sSite4, byGaugeFieldId, bPlusX, bPlusY, bPlusT), pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta124 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta124 = bPlusT ? _deviceEta3(sSite4, 2) : (_deviceEta3(site_target, 2) + 1);

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta124 = eta124 + 1;
        }

        if (eta124 & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    //if (bDDagger)
    //{
    //    _mul(result, -F(0.125) * fOmega);
    //}
    //else
    //{
    //    _mul(result, F(0.125) * fOmega);
    //}
    _mul(result, _make_cuComplex(F(0.0), F(-0.125) * fOmega));

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

    _add(pResultData[uiSiteIndex], result);
}

#pragma endregion

#pragma region Derivate

/**
 * Have n, n->n1, n->n2,
 * 1. we need to obtain V_(n, n1) , V_(n, n2)
 * 2. we need phi(n1), phi(n2), phid(n1), phid(n2)
 *
 * byContribution: 0 for mu, 1 for tau, 2 for both mu and tau
 * 
 * 0 for right-mu
 * 1 for right-tau
 * 2 for right-mu&left-t
 * 3 for only left-t
 * 4 for nothing
 *
 * iTau = 1 for +t, -1 for -t
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PR_XYTermT(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const BYTE* __restrict__ pEtaTable,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    Real shift,
    BYTE byMu, INT iTau,
    SCHAR pathLdir1, SCHAR pathLdir2, SCHAR pathLdir3, BYTE Llength,
    SCHAR pathRdir1, SCHAR pathRdir2, SCHAR pathRdir3, BYTE Rlength,
    BYTE byContribution)
{
    intokernalInt4;
    //const UINT uiBigIdx = __bi(sSite4);

    //=================================
    // 1. Find n1, n2
    SCHAR Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    SCHAR Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    SSmallInt4 site_n1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(_deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength))];
    if (sn1.IsDirichlet() || sn2.IsDirichlet())
    {
        return;
    }
    //const SSmallInt4 middleSite = _deviceSmallInt4OffsetC(site_n1, byMu + 1);
    //From now on, site_n1 is smiddle
    site_n1 = _deviceSmallInt4OffsetC(site_n1, byMu + 1);
    const SIndex& smiddle = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];

    site_n1 = __deviceSiteIndexToInt4(smiddle.m_uiSiteIndex);
    //y Dx and -x Dy
    const Real fNv = (0 == byMu)
        ? static_cast<Real>(site_n1.y - _DC_Centery + shift)
        : static_cast<Real>(_DC_Centerx - site_n1.x - shift);

    const Real eta_tau = (iTau > 0 ?
        ((pEtaTable[sn1.m_uiSiteIndex] >> 3) & 1)
        : ((pEtaTable[sn2.m_uiSiteIndex] >> 3) & 1))
        ? F(-1.0) : F(1.0);

    //=================================
    // 2. Find V(n,n1), V(n,n2)
    //left-T
    if (2 == byContribution || 3 == byContribution)
    {
        const deviceGauge vnn1 = _deviceLinkTSkipOne(pGauge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            //=================================
            // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
            deviceVector phi1 = (Llength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
            deviceVector phi2 = (Rlength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
            _mul(phi1, sn1.NeedToOppositeCoeff());
            _mul(phi2, sn2.NeedToOppositeCoeff());

            deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
            //This Add is required by partial(D^+D)
            phi1 = (Rlength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
            phi2 = (Llength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
            _mul(phi1, sn2.NeedToOppositeCoeff());
            _mul(phi2, sn1.NeedToOppositeCoeff());

            _sub(res, _makeContract<deviceGauge, deviceVector>(phi2, phi1));
            _mul(res, OneOver12 * fOmega * fNv * pNumerators[rfieldId] * eta_tau);
            _sub(pForce[linkIndex], res);
        }
    }

    if (0 == byContribution || 1 == byContribution || 2 == byContribution)
    {
        const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkTSkipOne(pGauge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, (1 == byContribution) ? 3 : byMu);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            //=================================
            // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
            deviceVector phi1 = (Llength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
            deviceVector phi2 = (Rlength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
            _mul(phi1, sn1.NeedToOppositeCoeff());
            _mul(phi2, sn2.NeedToOppositeCoeff());

            deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
            //This Add is required by partial(D^+D)
            phi1 = (Rlength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
            phi2 = (Llength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
            _mul(phi1, sn2.NeedToOppositeCoeff());
            _mul(phi2, sn1.NeedToOppositeCoeff());

            _sub(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
            _mul(res, OneOver12 * fOmega * fNv * pNumerators[rfieldId] * eta_tau);
            _sub(pForce[linkIndex], res);
        }
    }
#if 0
    const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, Llength, 1, Ldirs);
    const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, Rlength, 1, Rdirs);

    for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
    {
        const deviceVector* phi_i = pFermionPointers[rfieldId];
        const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];
        //=================================
        // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
        deviceVector phi1 = _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]);
        deviceVector phi2 = _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]);
        deviceVector phi3 = _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]);
        deviceVector phi4 = _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]);
        if (sn1.NeedToOpposite())
        {
            _mul(phi1, F(-1.0));
            _mul(phi3, F(-1.0));
        }
        if (sn2.NeedToOpposite())
        {
            _mul(phi2, F(-1.0));
            _mul(phi4, F(-1.0));
        }
        deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
        _add(res, _makeContract<deviceGauge, deviceVector>(phi4, phi3));
        _ta(res);
        const Real eta_tau = (iTau > 0 ?
            ((pEtaTable[sn1.m_uiSiteIndex] >> 3) & 1)
            : ((pEtaTable[sn2.m_uiSiteIndex] >> 3) & 1))
            ? F(-1.0) : F(1.0);
        _mul(res, OneOver12 * fOmega * fNv * pNumerators[rfieldId] * eta_tau);

        //For mu
        if (0 == byContribution || 2 == byContribution)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byMu);
            _sub(pForce[linkIndex], res);
        }

        //For tau
        if (1 == byContribution || 2 == byContribution)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);
            if (iTau > 0)
            {
                _sub(pForce[linkIndex], res);
            }
            else
            {
                _add(pForce[linkIndex], res);
            }
        }
    }
#endif
}

/**
 * Even-odd split version of _kernelDFermionKSForce_PR_XYTermT.
 * phi_i is non-zero only on even sites and phi_id only on odd sites, so for any
 * (n1, n2) pair exactly one of the two contracts is zero. Only the surviving
 * contract is computed, saving about half of the FLOPS.
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PR_XYTermT_EO(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const BYTE* __restrict__ pEtaTable,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    Real shift,
    BYTE byMu, INT iTau,
    SCHAR pathLdir1, SCHAR pathLdir2, SCHAR pathLdir3, BYTE Llength,
    SCHAR pathRdir1, SCHAR pathRdir2, SCHAR pathRdir3, BYTE Rlength,
    BYTE byContribution)
{
    intokernalInt4;

    //=================================
    // 1. Find n1, n2
    SCHAR Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    SCHAR Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    SSmallInt4 site_n1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(_deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength))];
    if (sn1.IsDirichlet() || sn2.IsDirichlet())
    {
        return;
    }
    site_n1 = _deviceSmallInt4OffsetC(site_n1, byMu + 1);
    const SIndex& smiddle = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];

    site_n1 = __deviceSiteIndexToInt4(smiddle.m_uiSiteIndex);
    //y Dx and -x Dy
    const Real fNv = (0 == byMu)
        ? static_cast<Real>(site_n1.y - _DC_Centery + shift)
        : static_cast<Real>(_DC_Centerx - site_n1.x - shift);

    const Real eta_tau = (iTau > 0 ?
        ((pEtaTable[sn1.m_uiSiteIndex] >> 3) & 1)
        : ((pEtaTable[sn2.m_uiSiteIndex] >> 3) & 1))
        ? F(-1.0) : F(1.0);

    //sn1/sn2 have opposite parities (path length is odd), so exactly one contract survives
    const UBOOL bSn1Odd = ((pEtaTable[sn1.m_uiSiteIndex] >> 4U) & 1U);

    //=================================
    // 2. Find V(n,n1), V(n,n2)
    //left-T
    if (2 == byContribution || 3 == byContribution)
    {
        const deviceGauge vnn1 = _deviceLinkTSkipOne(pGauge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            deviceGauge res;
            if (bSn1Odd)
            {
                deviceVector phi1 = (Llength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                deviceVector phi2 = (Rlength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
                _mul(phi1, sn1.NeedToOppositeCoeff());
                _mul(phi2, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
            }
            else
            {
                deviceVector phi2 = (Llength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                deviceVector phi1 = (Rlength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                _mul(phi2, sn1.NeedToOppositeCoeff());
                _mul(phi1, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
                //The other contract is zero, so res = 0 - res
                _mul(res, F(-1.0));
            }
            _mul(res, OneOver12 * fOmega * fNv * pNumerators[rfieldId] * eta_tau);
            _sub(pForce[linkIndex], res);
        }
    }

    if (0 == byContribution || 1 == byContribution || 2 == byContribution)
    {
        const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkTSkipOne(pGauge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, (1 == byContribution) ? 3 : byMu);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            deviceGauge res;
            if (bSn1Odd)
            {
                deviceVector phi1 = (Llength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                deviceVector phi2 = (Rlength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
                _mul(phi1, sn1.NeedToOppositeCoeff());
                _mul(phi2, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
            }
            else
            {
                deviceVector phi2 = (Llength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                deviceVector phi1 = (Rlength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                _mul(phi2, sn1.NeedToOppositeCoeff());
                _mul(phi1, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                //The other contract is zero, so res = 0 - res
                _mul(res, F(-1.0));
            }
            _mul(res, OneOver12 * fOmega * fNv * pNumerators[rfieldId] * eta_tau);
            _sub(pForce[linkIndex], res);
        }
    }
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PR_XYTermT_EM(
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pPhase,
    deviceGauge* pForce,
    const BYTE* __restrict__ pEtaTable,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    Real shift,
    Real fCharge,
    BYTE byMu, INT iTau,
    SCHAR pathLdir1, SCHAR pathLdir2, SCHAR pathLdir3, BYTE Llength,
    SCHAR pathRdir1, SCHAR pathRdir2, SCHAR pathRdir3, BYTE Rlength,
    BYTE byContribution)
{
    intokernalInt4;
    //const UINT uiBigIdx = __bi(sSite4);

    //=================================
    // 1. Find n1, n2
    SCHAR Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    SCHAR Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    SSmallInt4 site_n1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(_deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength))];
    //const SSmallInt4 middleSite = _deviceSmallInt4OffsetC(site_n1, byMu + 1);
    //From now on, site_n1 is smiddle
    site_n1 = _deviceSmallInt4OffsetC(site_n1, byMu + 1);
    const SIndex& smiddle = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];

    site_n1 = __deviceSiteIndexToInt4(smiddle.m_uiSiteIndex);
    //y Dx and -x Dy
    const Real fNv = (0 == byMu)
        ? static_cast<Real>(site_n1.y - _DC_Centery + shift)
        : static_cast<Real>(_DC_Centerx - site_n1.x - shift);

    const Real eta_tau = (iTau > 0 ?
        ((pEtaTable[sn1.m_uiSiteIndex] >> 3) & 1)
        : ((pEtaTable[sn2.m_uiSiteIndex] >> 3) & 1))
        ? F(-1.0) : F(1.0);

    //=================================
    // 2. Find V(n,n1), V(n,n2)
    //left-T
    if (2 == byContribution || 3 == byContribution)
    {
        const deviceGauge vnn1 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            //=================================
            // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
            deviceVector phi1 = (Llength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
            deviceVector phi2 = (Rlength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
            _mul(phi1, sn1.NeedToOppositeCoeff());
            _mul(phi2, sn2.NeedToOppositeCoeff());

            deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
            //This Add is required by partial(D^+D)
            phi1 = (Rlength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
            phi2 = (Llength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
            _mul(phi1, sn2.NeedToOppositeCoeff());
            _mul(phi2, sn1.NeedToOppositeCoeff());

            _sub(res, _makeContract<deviceGauge, deviceVector>(phi2, phi1));
            _mul(res, OneOver12 * fOmega * fNv * pNumerators[rfieldId] * eta_tau);
            _sub(pForce[linkIndex], res);
        }
    }

    if (0 == byContribution || 1 == byContribution || 2 == byContribution)
    {
        const deviceGauge vnn1 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, (1 == byContribution) ? 3 : byMu);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            //=================================
            // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
            deviceVector phi1 = (Llength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
            deviceVector phi2 = (Rlength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
            _mul(phi1, sn1.NeedToOppositeCoeff());
            _mul(phi2, sn2.NeedToOppositeCoeff());

            deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
            //This Add is required by partial(D^+D)
            phi1 = (Rlength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
            phi2 = (Llength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
            _mul(phi1, sn2.NeedToOppositeCoeff());
            _mul(phi2, sn1.NeedToOppositeCoeff());

            _sub(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
            _mul(res, OneOver12 * fOmega * fNv * pNumerators[rfieldId] * eta_tau);
            _sub(pForce[linkIndex], res);
        }
    }
}

/**
 *
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PR_XYTau_TermT(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    SCHAR pathLdir1, SCHAR pathLdir2, SCHAR pathLdir3, BYTE Llength,
    SCHAR pathRdir1, SCHAR pathRdir2, SCHAR pathRdir3, BYTE Rlength)
{
    intokernalInt4;
    //const UINT uiBigIdx = __bi(sSite4);

    //=================================
    // 1. Find n1, n2
    SCHAR Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    SCHAR Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];
    if (sn1.IsDirichlet() || sn2.IsDirichlet())
    { 
        return;
    }
    //Why use sn2? shouldn't it be sn1?
    const Real eta124 = _deviceEta124(__deviceSiteIndexToInt4(sn1.m_uiSiteIndex));
    //=================================
    // 2. Find V(n,n1), V(n,n2)
    if (pathLdir1 > 0)
    {
        const deviceGauge vnn1 = _deviceLinkTSkipOne(pGauge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLdir1 - 1);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            //=================================
            // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
            deviceVector phi1 = (Llength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
            deviceVector phi2 = (Rlength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
            _mul(phi1, sn1.NeedToOppositeCoeff());
            _mul(phi2, sn2.NeedToOppositeCoeff());

            deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
            //This Add is required by partial(D^+D)
            phi1 = (Rlength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
            phi2 = (Llength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
            _mul(phi1, sn2.NeedToOppositeCoeff());
            _mul(phi2, sn1.NeedToOppositeCoeff());

            _sub(res, _makeContract<deviceGauge, deviceVector>(phi2, phi1));
            _mul(res, OneOver48 * static_cast<Real>(fOmega) * pNumerators[rfieldId] * eta124);
            _sub(pForce[linkIndex], res);
        }
    }

    if (pathRdir1 > 0)
    {
        const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkTSkipOne(pGauge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRdir1 - 1);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            //=================================
            // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
            deviceVector phi1 = (Llength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
            deviceVector phi2 = (Rlength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
            _mul(phi1, sn1.NeedToOppositeCoeff());
            _mul(phi2, sn2.NeedToOppositeCoeff());

            deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
            //This Add is required by partial(D^+D)
            phi1 = (Rlength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
            phi2 = (Llength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
            _mul(phi1, sn2.NeedToOppositeCoeff());
            _mul(phi2, sn1.NeedToOppositeCoeff());

            _sub(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
            _mul(res, OneOver48 * static_cast<Real>(fOmega) * pNumerators[rfieldId] * eta124);
            _sub(pForce[linkIndex], res);
        }
    }
#if 0
    const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, Llength, 1, Ldirs);
    const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, Rlength, 1, Rdirs);

    for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
    {
        const deviceVector* phi_i = pFermionPointers[rfieldId];
        const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

        //=================================
        // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
        deviceVector phi1 = _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]);
        deviceVector phi2 = _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]);
        deviceVector phi3 = _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]);
        deviceVector phi4 = _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]);
        if (sn1.NeedToOpposite())
        {
            _mul(phi1, F(-1.0));
            _mul(phi3, F(-1.0));
        }
        if (sn2.NeedToOpposite())
        {
            _mul(phi2, F(-1.0));
            _mul(phi4, F(-1.0));
        }
        deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
        //This was phi2 phi1+ * eta124(n1) - phi3 phi4+ * eta124(n2)
        //The sign of the second term is because of 'dagger'
        //However, eta124(n1) = -eta124(n2), so use Add directly.
        _add(res, _makeContract<deviceGauge, deviceVector>(phi4, phi3));
        _ta(res);
        _mul(res, OneOver48 * static_cast<Real>(fOmega) * pNumerators[rfieldId] * eta124);

        //Use eta124 of n2 so Add left Sub right
        //Change to use eta124 of n1, Sub left and Add right
        if (pathLdir1 > 0)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLdir1 - 1);
            _add(pForce[linkIndex], res);
        }

        if (pathRdir1 > 0)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRdir1 - 1);
            _sub(pForce[linkIndex], res);
        }
    }
#endif
}

/**
 * Even-odd split version of _kernelDFermionKSForce_PR_XYTau_TermT.
 * Same as _kernelDFermionKSForce_PR_XYTermT_EO: only the surviving contract is
 * computed (phi_i lives on even sites, phi_id lives on odd sites).
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PR_XYTau_TermT_EO(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    SCHAR pathLdir1, SCHAR pathLdir2, SCHAR pathLdir3, BYTE Llength,
    SCHAR pathRdir1, SCHAR pathRdir2, SCHAR pathRdir3, BYTE Rlength)
{
    intokernalInt4;

    //=================================
    // 1. Find n1, n2
    SCHAR Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    SCHAR Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];
    if (sn1.IsDirichlet() || sn2.IsDirichlet())
    { 
        return;
    }
    const Real eta124 = _deviceEta124(__deviceSiteIndexToInt4(sn1.m_uiSiteIndex));
    //sn1/sn2 have opposite parities (path length is odd), so exactly one contract survives
    const UBOOL bSn1Odd = ((__idx->m_pEtaMu[sn1.m_uiSiteIndex] >> 4U) & 1U);

    //=================================
    // 2. Find V(n,n1), V(n,n2)
    if (pathLdir1 > 0)
    {
        const deviceGauge vnn1 = _deviceLinkTSkipOne(pGauge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLdir1 - 1);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            deviceGauge res;
            if (bSn1Odd)
            {
                deviceVector phi1 = (Llength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                deviceVector phi2 = (Rlength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
                _mul(phi1, sn1.NeedToOppositeCoeff());
                _mul(phi2, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
            }
            else
            {
                deviceVector phi2 = (Llength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                deviceVector phi1 = (Rlength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                _mul(phi2, sn1.NeedToOppositeCoeff());
                _mul(phi1, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
                //The other contract is zero, so res = 0 - res
                _mul(res, F(-1.0));
            }
            _mul(res, OneOver48 * static_cast<Real>(fOmega) * pNumerators[rfieldId] * eta124);
            _sub(pForce[linkIndex], res);
        }
    }

    if (pathRdir1 > 0)
    {
        const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkTSkipOne(pGauge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRdir1 - 1);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            deviceGauge res;
            if (bSn1Odd)
            {
                deviceVector phi1 = (Llength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                deviceVector phi2 = (Rlength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
                _mul(phi1, sn1.NeedToOppositeCoeff());
                _mul(phi2, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
            }
            else
            {
                deviceVector phi2 = (Llength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                deviceVector phi1 = (Rlength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                _mul(phi2, sn1.NeedToOppositeCoeff());
                _mul(phi1, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                //The other contract is zero, so res = 0 - res
                _mul(res, F(-1.0));
            }
            _mul(res, OneOver48 * static_cast<Real>(fOmega) * pNumerators[rfieldId] * eta124);
            _sub(pForce[linkIndex], res);
        }
    }
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTermEMT_EO(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pPhase,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    Real fCharge,
    SSmallInt4 sCenter,
    Real shift,
    UBOOL bEven,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalEOHalf;
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    deviceVector result = _makeZero<deviceVector>();
    //const INT eta_tau = ((pEtaTable[uiSiteIndex] >> 3) & 1);
    const INT eta_tau = eta >> 3;

    #pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusMu = idx & 2;
        const UBOOL bPlusTau = idx & 4;
        //x or y, and y or x is the derivate, not coefficient
        const UINT bXorY = idx & 1;
        const UINT bYorX = 1 - bXorY;
        SSmallInt4 sTargetSite = sSite4;
        SSmallInt4 sMidSite = sSite4;
        sTargetSite.m_byData4[bYorX] = sTargetSite.m_byData4[bYorX] + (bPlusMu ? 2 : -2);
        sMidSite.m_byData4[bYorX] = sMidSite.m_byData4[bYorX] + (bPlusMu ? 1 : -1);
        sTargetSite.w = sTargetSite.w + (bPlusTau ? 1 : -1);
        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sTargetSite)];
        if (sTargetBigIndex.IsDirichlet())
        {
            continue;
        }
        const SIndex& sMiddleBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sMidSite)];
        sMidSite = __deviceSiteIndexToInt4(sMiddleBigIndex.m_uiSiteIndex);

        //note that bYorX = 1, it is x partial_y term, therefore is '-'
        //INT this_eta_tau = (bPlusTau ? eta_tau : ((pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3) & 1))
        INT this_eta_tau = (bPlusTau ? eta_tau : (pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3))
            + bYorX;

        if (sTargetBigIndex.NeedToOpposite())
        {
            this_eta_tau = this_eta_tau + 1;
        }

        deviceVector right = _mulVec(_deviceVXXTauOptimizedEMT(pGauge, pPhase, sSite4, fCharge, byGaugeFieldId, bXorY, bPlusMu, bPlusTau),
            pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        //when bXorY = 1, it is y partial _x, so is [1]
        //when bXorY = 0, it is x partial _y, so is [0]
        _mul(right, sMidSite.m_byData4[bXorY] - sCenter.m_byData4[bXorY] + shift);

        if (!bPlusMu)
        {
            //for -2x, -2y terms, there is another minus sign
            this_eta_tau = this_eta_tau + 1;
        }

        if (this_eta_tau & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, F(-0.25) * fOmega);
    }
    else
    {
        _mul(result, F(0.25) * fOmega);
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

    _add(pDeviceData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTau_TermEMT_EO(
    deviceVector* pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pPhase,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    Real fCharge,
    UBOOL bEven,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalEOHalf;
    const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

    deviceVector result = _makeZero<deviceVector>();

#ifndef _CLG_DTK
    #pragma unroll
#endif
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusX = (0 != (idx & 1));
        const UBOOL bPlusY = (0 != (idx & 2));
        const UBOOL bPlusT = (0 != (idx & 4));

        SSmallInt4 sOffset = sSite4;
        sOffset.x = sOffset.x + (bPlusX ? 1 : -1);
        sOffset.y = sOffset.y + (bPlusY ? 1 : -1);
        sOffset.w = sOffset.w + (bPlusT ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];
        if (sTargetBigIndex.IsDirichlet())
        {
            continue;
        }

        const deviceVector right = _mulVec(_deviceVXYTOptimizedEMT(pGauge, pPhase, sSite4, fCharge, byGaugeFieldId, bPlusX, bPlusY, bPlusT), pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta124 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta124 = bPlusT ? _deviceEta3(sSite4, 2) : (_deviceEta3(site_target, 2) + 1);

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta124 = eta124 + 1;
        }

        if (eta124 & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, -F(0.125) * fOmega);
    }
    else
    {
        _mul(result, F(0.125) * fOmega);
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

    _add(pDeviceData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PR_XYTermT_EM_EO(
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pPhase,
    deviceGauge* pForce,
    const BYTE* __restrict__ pEtaTable,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    Real shift,
    Real fCharge,
    BYTE byMu, INT iTau,
    SCHAR pathLdir1, SCHAR pathLdir2, SCHAR pathLdir3, BYTE Llength,
    SCHAR pathRdir1, SCHAR pathRdir2, SCHAR pathRdir3, BYTE Rlength,
    BYTE byContribution)
{
    intokernalInt4;

    //=================================
    // 1. Find n1, n2
    SCHAR Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    SCHAR Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    SSmallInt4 site_n1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(_deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength))];
    if (sn1.IsDirichlet() || sn2.IsDirichlet())
    {
        return;
    }
    site_n1 = _deviceSmallInt4OffsetC(site_n1, byMu + 1);
    const SIndex& smiddle = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];

    site_n1 = __deviceSiteIndexToInt4(smiddle.m_uiSiteIndex);
    //y Dx and -x Dy
    const Real fNv = (0 == byMu)
        ? static_cast<Real>(site_n1.y - _DC_Centery + shift)
        : static_cast<Real>(_DC_Centerx - site_n1.x - shift);

    const Real eta_tau = (iTau > 0 ?
        ((pEtaTable[sn1.m_uiSiteIndex] >> 3) & 1)
        : ((pEtaTable[sn2.m_uiSiteIndex] >> 3) & 1))
        ? F(-1.0) : F(1.0);

    //sn1/sn2 have opposite parities (path length is odd), so exactly one contract survives
    const UBOOL bSn1Odd = ((pEtaTable[sn1.m_uiSiteIndex] >> 4U) & 1U);

    //=================================
    // 2. Find V(n,n1), V(n,n2)
    //left-T
    if (2 == byContribution || 3 == byContribution)
    {
        const deviceGauge vnn1 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            deviceGauge res;
            if (bSn1Odd)
            {
                deviceVector phi1 = (Llength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                deviceVector phi2 = (Rlength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
                _mul(phi1, sn1.NeedToOppositeCoeff());
                _mul(phi2, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
            }
            else
            {
                deviceVector phi2 = (Llength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                deviceVector phi1 = (Rlength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                _mul(phi2, sn1.NeedToOppositeCoeff());
                _mul(phi1, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
                //The other contract is zero, so res = 0 - res
                _mul(res, F(-1.0));
            }
            _mul(res, OneOver12 * fOmega * fNv * pNumerators[rfieldId] * eta_tau);
            _sub(pForce[linkIndex], res);
        }
    }

    if (0 == byContribution || 1 == byContribution || 2 == byContribution)
    {
        const deviceGauge vnn1 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, (1 == byContribution) ? 3 : byMu);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            deviceGauge res;
            if (bSn1Odd)
            {
                deviceVector phi1 = (Llength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                deviceVector phi2 = (Rlength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
                _mul(phi1, sn1.NeedToOppositeCoeff());
                _mul(phi2, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
            }
            else
            {
                deviceVector phi2 = (Llength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                deviceVector phi1 = (Rlength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                _mul(phi2, sn1.NeedToOppositeCoeff());
                _mul(phi1, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                //The other contract is zero, so res = 0 - res
                _mul(res, F(-1.0));
            }
            _mul(res, OneOver12 * fOmega * fNv * pNumerators[rfieldId] * eta_tau);
            _sub(pForce[linkIndex], res);
        }
    }
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PR_XYTau_TermT_EM_EO(
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pPhase,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    Real fCharge,
    SCHAR pathLdir1, SCHAR pathLdir2, SCHAR pathLdir3, BYTE Llength,
    SCHAR pathRdir1, SCHAR pathRdir2, SCHAR pathRdir3, BYTE Rlength)
{
    intokernalInt4;

    //=================================
    // 1. Find n1, n2
    SCHAR Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    SCHAR Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];
    if (sn1.IsDirichlet() || sn2.IsDirichlet())
    { 
        return;
    }
    const Real eta124 = _deviceEta124(__deviceSiteIndexToInt4(sn1.m_uiSiteIndex));
    //sn1/sn2 have opposite parities (path length is odd), so exactly one contract survives
    const UBOOL bSn1Odd = ((__idx->m_pEtaMu[sn1.m_uiSiteIndex] >> 4U) & 1U);

    //=================================
    // 2. Find V(n,n1), V(n,n2)
    if (pathLdir1 > 0)
    {
        const deviceGauge vnn1 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLdir1 - 1);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            deviceGauge res;
            if (bSn1Odd)
            {
                deviceVector phi1 = (Llength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                deviceVector phi2 = (Rlength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
                _mul(phi1, sn1.NeedToOppositeCoeff());
                _mul(phi2, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
            }
            else
            {
                deviceVector phi2 = (Llength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                deviceVector phi1 = (Rlength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                _mul(phi2, sn1.NeedToOppositeCoeff());
                _mul(phi1, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
                //The other contract is zero, so res = 0 - res
                _mul(res, F(-1.0));
            }
            _mul(res, OneOver48 * static_cast<Real>(fOmega) * pNumerators[rfieldId] * eta124);
            _sub(pForce[linkIndex], res);
        }
    }

    if (pathRdir1 > 0)
    {
        const deviceGauge vnn1 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRdir1 - 1);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            deviceGauge res;
            if (bSn1Odd)
            {
                deviceVector phi1 = (Llength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
                deviceVector phi2 = (Rlength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
                _mul(phi1, sn1.NeedToOppositeCoeff());
                _mul(phi2, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
            }
            else
            {
                deviceVector phi2 = (Llength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
                deviceVector phi1 = (Rlength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
                _mul(phi2, sn1.NeedToOppositeCoeff());
                _mul(phi1, sn2.NeedToOppositeCoeff());
                res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                //The other contract is zero, so res = 0 - res
                _mul(res, F(-1.0));
            }
            _mul(res, OneOver48 * static_cast<Real>(fOmega) * pNumerators[rfieldId] * eta124);
            _sub(pForce[linkIndex], res);
        }
    }
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PR_XYTau_TermT_EM(
    const deviceGauge* __restrict__ pGauge,
    const Real* __restrict__ pPhase,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    Real fCharge,
    SCHAR pathLdir1, SCHAR pathLdir2, SCHAR pathLdir3, BYTE Llength,
    SCHAR pathRdir1, SCHAR pathRdir2, SCHAR pathRdir3, BYTE Rlength)
{
    intokernalInt4;
    //const UINT uiBigIdx = __bi(sSite4);

    //=================================
    // 1. Find n1, n2
    SCHAR Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    SCHAR Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];

    //Why use sn2? shouldn't it be sn1?
    const Real eta124 = _deviceEta124(__deviceSiteIndexToInt4(sn1.m_uiSiteIndex));
    //=================================
    // 2. Find V(n,n1), V(n,n2)
    if (pathLdir1 > 0)
    {
        const deviceGauge vnn1 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLdir1 - 1);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            //=================================
            // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
            deviceVector phi1 = (Llength > 1) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
            deviceVector phi2 = (Rlength > 0) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
            _mul(phi1, sn1.NeedToOppositeCoeff());
            _mul(phi2, sn2.NeedToOppositeCoeff());

            deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
            //This Add is required by partial(D^+D)
            phi1 = (Rlength > 0) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
            phi2 = (Llength > 1) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
            _mul(phi1, sn2.NeedToOppositeCoeff());
            _mul(phi2, sn1.NeedToOppositeCoeff());

            _sub(res, _makeContract<deviceGauge, deviceVector>(phi2, phi1));
            _mul(res, OneOver48 * static_cast<Real>(fOmega) * pNumerators[rfieldId] * eta124);
            _sub(pForce[linkIndex], res);
        }
    }

    if (pathRdir1 > 0)
    {
        const deviceGauge vnn1 = _deviceLinkEMT(pGauge, pPhase, fCharge, sSite4, Llength, byGaugeFieldId, Ldirs);
        const deviceGauge vnn2 = _deviceLinkEMTSkipOne(pGauge, pPhase, fCharge, sSite4, Rlength, byGaugeFieldId, Rdirs);

        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRdir1 - 1);

        for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
        {
            const deviceVector* phi_i = pFermionPointers[rfieldId];
            const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

            //=================================
            // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
            deviceVector phi1 = (Llength > 0) ? _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]) : phi_id[sn1.m_uiSiteIndex];
            deviceVector phi2 = (Rlength > 1) ? _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]) : phi_i[sn2.m_uiSiteIndex];
            _mul(phi1, sn1.NeedToOppositeCoeff());
            _mul(phi2, sn2.NeedToOppositeCoeff());

            deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi2, phi1);
            //This Add is required by partial(D^+D)
            phi1 = (Rlength > 1) ? _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]) : phi_id[sn2.m_uiSiteIndex];
            phi2 = (Llength > 0) ? _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]) : phi_i[sn1.m_uiSiteIndex];
            _mul(phi1, sn2.NeedToOppositeCoeff());
            _mul(phi2, sn1.NeedToOppositeCoeff());

            _sub(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
            _mul(res, OneOver48 * static_cast<Real>(fOmega) * pNumerators[rfieldId] * eta124);
            _sub(pForce[linkIndex], res);
        }
    }
}

#pragma endregion

#pragma region D and derivate

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelR<deviceVector, deviceGauge, vectorN>::DOperatorKS_R_RealRotation(
    DOUBLE fOmega, 
    UBOOL bShiftHalfCoord, deviceVector* pTarget, const deviceVector* pSource,
    const deviceGauge* pGauge, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    //CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, byGaugeFieldId, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);

    preparethread;
    if (bDagger)
    {
        appCrucial(_T("D dagger is not supported for real rotation!\n"));
    }

    _LAUNCH_KERNEL(_kernelDFermionKS_PR_XYTermRealT TMPARG(deviceVector, deviceGauge), block, threads,
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        _HC_Center,
        bShiftHalfCoord ? F(0.5) : F(0.0),
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    _LAUNCH_KERNEL(_kernelDFermionKS_PR_XYTau_TermRealT TMPARG(deviceVector, deviceGauge), block, threads,
        pSource,
        pGauge,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelR<deviceVector, deviceGauge, vectorN>::DOperatorKS_R_ImaginaryRotation(
    DOUBLE fOmega,
    UBOOL bShiftHalfCoord, deviceVector* pTarget, const deviceVector* pSource,
    const deviceGauge* pGauge, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionKS_PR_XYTermT TMPARG(deviceVector, deviceGauge), block, threads,
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        _HC_Center,
        bShiftHalfCoord ? F(0.5) : F(0.0),
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    _LAUNCH_KERNEL(_kernelDFermionKS_PR_XYTau_TermT TMPARG(deviceVector, deviceGauge), block, threads,
        pSource,
        pGauge,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelR<deviceVector, deviceGauge, vectorN>::DOperatorKSOnEvenOrOdd_R_ImaginaryRotation(
    DOUBLE fOmega,
    UBOOL bShiftHalfCoord, deviceVector* pData,
    const deviceGauge* pGauge, BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEven,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    //same convention as CFieldFermionKSTKernel::DOperatorKSOnEvenOrOdd, bEven is passed to intokernalEOHalf directly
    preparethreadHalf;
    _LAUNCH_KERNEL(_kernelDFermionKS_PR_XYTermT_EO TMPARG(deviceVector, deviceGauge), block, threads,
        pData,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        _HC_Center,
        bShiftHalfCoord ? F(0.5) : F(0.0),
        bEven,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    _LAUNCH_KERNEL(_kernelDFermionKS_PR_XYTau_TermT_EO TMPARG(deviceVector, deviceGauge), block, threads,
        pData,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        bEven,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelR<deviceVector, deviceGauge, vectorN>::DOperatorKSOnEvenOrOdd_R_ImaginaryRotation_Cached(
    DOUBLE fOmega,
    UBOOL bShiftHalfCoord, deviceVector* pData,
    const deviceGauge* pCachedGauge, BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEven,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    //same convention as CFieldFermionKSTKernel::DOperatorKSOnEvenOrOdd, bEven is passed to intokernalEOHalf directly
    preparethreadHalf;
    _LAUNCH_KERNEL(_kernelDFermionKS_PR_Cached_XYTermT_EO TMPARG(deviceVector, deviceGauge), block, threads,
        pData,
        pCachedGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        _HC_Center,
        bShiftHalfCoord ? F(0.5) : F(0.0),
        bEven,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    _LAUNCH_KERNEL(_kernelDFermionKS_PR_Cached_XYTau_TermT_EO TMPARG(deviceVector, deviceGauge), block, threads,
        pData,
        pCachedGauge + _HC_Volume * 8,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        bEven,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelR<deviceVector, deviceGauge, vectorN>::DOperatorKSOnEvenOdd_R_ImaginaryRotation_EM(
    DOUBLE fOmega,
    UBOOL bShiftHalfCoord, deviceVector* pData,
    const deviceGauge* pGauge, const Real* pPhase, Real fCharge,
    BYTE byFieldId, BYTE byGaugeFieldId, UBOOL bEven,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    //same convention as CFieldFermionKSTKernel::DOperatorKSOnEvenOrOdd, bEven is passed to intokernalEOHalf directly
    preparethreadHalf;
    _LAUNCH_KERNEL(_kernelDFermionKS_PR_XYTermEMT_EO TMPARG(deviceVector, deviceGauge), block, threads,
        pData,
        pGauge,
        pPhase,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        fCharge,
        _HC_Center,
        bShiftHalfCoord ? F(0.5) : F(0.0),
        bEven,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    _LAUNCH_KERNEL(_kernelDFermionKS_PR_XYTau_TermEMT_EO TMPARG(deviceVector, deviceGauge), block, threads,
        pData,
        pGauge,
        pPhase,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        fCharge,
        bEven,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelR<deviceVector, deviceGauge, vectorN>::DOperatorKS_R_ImaginaryRotation_Cached(
    DOUBLE fOmega,
    UBOOL bShiftHalfCoord, deviceVector* pTarget, const deviceVector* pSource,
    const deviceGauge* pGaugeBuffer, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionKS_PR_Cached_XYTermT TMPARG(deviceVector, deviceGauge), block, threads,
        pSource,
        pGaugeBuffer,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        _HC_Center,
        bShiftHalfCoord ? F(0.5) : F(0.0),
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    _LAUNCH_KERNEL(_kernelDFermionKS_PR_Cached_XYTau_TermT TMPARG(deviceVector, deviceGauge), block, threads,
        pSource,
        pGaugeBuffer + _HC_Volume * 8,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelR<deviceVector, deviceGauge, vectorN>::DOperatorKS_R_ImaginaryRotation_EM(
    DOUBLE fOmega, Real fCharge,
    UBOOL bShiftHalfCoord, deviceVector* pTarget, const deviceVector* pSource,
    const deviceGauge* pGaugeBuffer, const Real* pPhase, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelDFermionKS_PR_XYTermEMT TMPARG(deviceVector, deviceGauge), block, threads,
        pSource,
        pGaugeBuffer,
        pPhase,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        fCharge,
        _HC_Center,
        bShiftHalfCoord ? F(0.5) : F(0.0),
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    _LAUNCH_KERNEL(_kernelDFermionKS_PR_XYTau_TermEMT TMPARG(deviceVector, deviceGauge), block, threads,
        pSource,
        pGaugeBuffer,
        pPhase,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        fCharge,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelR<deviceVector, deviceGauge, vectorN>::DerivateD0_R(
    DOUBLE fOmega,
    UBOOL bShiftCenter,
    const deviceVector* pFermion,
    BYTE byFieldId,
    deviceGauge* pForce,
    const deviceGauge* pGaugeBuffer,
    BYTE byGaugeFieldId,
    const deviceVector* const* pRationalFields,
    const Real* pNumerator,
    UINT uiRationApproxOrder)
{
    preparethread;

    #pragma region X Y Term
    SCHAR mu[2] = { 0, 1 };
    for (INT imu = 0; imu < 2; ++imu)
    {
        SCHAR dirs[6][3] =
        {
            {static_cast<SCHAR>(4), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(4), static_cast<SCHAR>(mu[imu] + 1)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(4)},
            //{4, -mu[imu] - 1, -mu[imu] - 1},
            //{-mu[imu] - 1, 4, -mu[imu] - 1},
            //{-mu[imu] - 1, -mu[imu] - 1, 4},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(-4)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(-4), static_cast<SCHAR>(mu[imu] + 1)},
            {static_cast<SCHAR>(-4), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1)},
        };

        INT iTau[6] = { 1, 1, 1, -1, -1, -1 };

        //0 for right-mu
        //1 for right-tau
        //2 for right-mu&left-t
        //3 for only left-t
        //4 for nothing
        BYTE contributionOf[6][4] =
        {
            {1, 0, 0, 4},
            {0, 1, 0, 4},
            {0, 0, 1, 4},
            //{1, 4, 0, 0},
            //{4, 2, 4, 0},
            //{4, 0, 2, 4},
            {0, 0, 4, 3},
            {0, 4, 2, 4},
            {4, 2, 0, 4},
        };

        for (INT pathidx = 0; pathidx < 6; ++pathidx)
        {
            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                if (4 == contributionOf[pathidx][iSeperation])
                {
                    continue;
                }

                SCHAR L[3] = { 0, 0, 0 };
                SCHAR R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(dirs[pathidx], iSeperation, L, R, LLength, RLength);

                _LAUNCH_KERNEL(_kernelDFermionKSForce_PR_XYTermT TMPARG(deviceVector, deviceGauge), block, threads,
                    pGaugeBuffer,
                    pForce,
                    appGetLattice()->m_pIndexCache->m_pEtaMu,
                    pRationalFields,
                    pNumerator,
                    uiRationApproxOrder,
                    byFieldId,
                    byGaugeFieldId,
                    fOmega,
                    bShiftCenter ? F(0.5) : F(0.0),
                    static_cast<BYTE>(imu), iTau[pathidx],
                    L[0], L[1], L[2], LLength,
                    R[0], R[1], R[2], RLength,
                    contributionOf[pathidx][iSeperation]
                    );
            }
        }
    }

#pragma endregion

#pragma region Polarization term

    //===========================
    //polarization terms
    //ilinkType is +-x +-y +t,
    //INT linkTypes[4][3] =
    //{
    //    {1, 2, 4},
    //    {1, 2, -4},
    //    {-1, 2, 4},
    //    {-1, 2, -4}
    //};
    SCHAR linkTypes[4][3] =
    {
        {1, 2, 4},
        {1, -2, 4},
        {-1, 2, 4},
        {-1, -2, 4}
    };

    for (INT ilinkType = 0; ilinkType < 4; ++ilinkType)
    {
        SCHAR sixlinks[6][3] =
        {
            {linkTypes[ilinkType][0], linkTypes[ilinkType][1], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][0], linkTypes[ilinkType][2], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][0], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][2], linkTypes[ilinkType][0]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][0], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][1], linkTypes[ilinkType][0]}
        };

        for (INT isixtype = 0; isixtype < 6; ++isixtype)
        {
            //bearly no change of time, because force calculation is not frequent
            /*
            _LAUNCH_KERNEL(_giveupkernelDFermionKSForce_PR_XYTau_Term2, block, threads, 
                (const deviceGauge*)pGaugeBuffer,
                (deviceGauge*)pForce,
                m_pRationalFieldPointers,
                m_pMDNumerator,
                m_rMD.m_uiDegree,
                m_byFieldId,
                CCommonData::m_fOmega,
                sixlinks[isixtype][0], sixlinks[isixtype][1], sixlinks[isixtype][2]
                );
            */

            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                SCHAR L[3] = { 0, 0, 0 };
                SCHAR R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(sixlinks[isixtype], iSeperation, L, R, LLength, RLength);

                const UBOOL bHasLeft = (LLength > 0) && (L[0] > 0);
                const UBOOL bHasRight = (RLength > 0) && (R[0] > 0);

                if (bHasLeft || bHasRight)
                {
                    _LAUNCH_KERNEL(_kernelDFermionKSForce_PR_XYTau_TermT TMPARG(deviceVector, deviceGauge), block, threads,
                        pGaugeBuffer,
                        pForce,
                        pRationalFields,
                        pNumerator,
                        uiRationApproxOrder,
                        byFieldId,
                        byGaugeFieldId,
                        fOmega,
                        L[0], L[1], L[2], LLength,
                        R[0], R[1], R[2], RLength
                        );
                }
            }
        }
    }

#pragma endregion
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelR<deviceVector, deviceGauge, vectorN>::DerivateD0_ROnEvenOdd(
    DOUBLE fOmega,
    UBOOL bShiftCenter,
    const deviceVector* pFermion,
    BYTE byFieldId,
    deviceGauge* pForce,
    const deviceGauge* pGaugeBuffer,
    BYTE byGaugeFieldId,
    const deviceVector* const* pRationalFields,
    const Real* pNumerator,
    UINT uiRationApproxOrder)
{
    preparethread;

    #pragma region X Y Term
    SCHAR mu[2] = { 0, 1 };
    for (INT imu = 0; imu < 2; ++imu)
    {
        SCHAR dirs[6][3] =
        {
            {static_cast<SCHAR>(4), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(4), static_cast<SCHAR>(mu[imu] + 1)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(4)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(-4)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(-4), static_cast<SCHAR>(mu[imu] + 1)},
            {static_cast<SCHAR>(-4), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1)},
        };

        INT iTau[6] = { 1, 1, 1, -1, -1, -1 };

        BYTE contributionOf[6][4] =
        {
            {1, 0, 0, 4},
            {0, 1, 0, 4},
            {0, 0, 1, 4},
            {0, 0, 4, 3},
            {0, 4, 2, 4},
            {4, 2, 0, 4},
        };

        for (INT pathidx = 0; pathidx < 6; ++pathidx)
        {
            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                if (4 == contributionOf[pathidx][iSeperation])
                {
                    continue;
                }

                SCHAR L[3] = { 0, 0, 0 };
                SCHAR R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(dirs[pathidx], iSeperation, L, R, LLength, RLength);

                _LAUNCH_KERNEL(_kernelDFermionKSForce_PR_XYTermT_EO TMPARG(deviceVector, deviceGauge), block, threads,
                    pGaugeBuffer,
                    pForce,
                    appGetLattice()->m_pIndexCache->m_pEtaMu,
                    pRationalFields,
                    pNumerator,
                    uiRationApproxOrder,
                    byFieldId,
                    byGaugeFieldId,
                    fOmega,
                    bShiftCenter ? F(0.5) : F(0.0),
                    static_cast<BYTE>(imu), iTau[pathidx],
                    L[0], L[1], L[2], LLength,
                    R[0], R[1], R[2], RLength,
                    contributionOf[pathidx][iSeperation]
                    );
            }
        }
    }

#pragma endregion

#pragma region Polarization term

    SCHAR linkTypes[4][3] =
    {
        {1, 2, 4},
        {1, -2, 4},
        {-1, 2, 4},
        {-1, -2, 4}
    };

    for (INT ilinkType = 0; ilinkType < 4; ++ilinkType)
    {
        SCHAR sixlinks[6][3] =
        {
            {linkTypes[ilinkType][0], linkTypes[ilinkType][1], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][0], linkTypes[ilinkType][2], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][0], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][2], linkTypes[ilinkType][0]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][0], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][1], linkTypes[ilinkType][0]}
        };

        for (INT isixtype = 0; isixtype < 6; ++isixtype)
        {
            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                SCHAR L[3] = { 0, 0, 0 };
                SCHAR R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(sixlinks[isixtype], iSeperation, L, R, LLength, RLength);

                const UBOOL bHasLeft = (LLength > 0) && (L[0] > 0);
                const UBOOL bHasRight = (RLength > 0) && (R[0] > 0);

                if (bHasLeft || bHasRight)
                {
                    _LAUNCH_KERNEL(_kernelDFermionKSForce_PR_XYTau_TermT_EO TMPARG(deviceVector, deviceGauge), block, threads,
                        pGaugeBuffer,
                        pForce,
                        pRationalFields,
                        pNumerator,
                        uiRationApproxOrder,
                        byFieldId,
                        byGaugeFieldId,
                        fOmega,
                        L[0], L[1], L[2], LLength,
                        R[0], R[1], R[2], RLength
                        );
                }
            }
        }
    }

#pragma endregion
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelR<deviceVector, deviceGauge, vectorN>::DerivateD0_REM(
    DOUBLE fOmega,
    UBOOL bShiftCenter,
    Real fCharge,
    const deviceVector* pFermion,
    BYTE byFieldId,
    deviceGauge* pForce,
    const deviceGauge* pGaugeBuffer,
    const Real* pPhaseBuffer,
    BYTE byGaugeFieldId,
    const deviceVector* const* pRationalFields,
    const Real* pNumerator,
    UINT uiRationApproxOrder)
{
    preparethread;

#pragma region X Y Term
    SCHAR mu[2] = { 0, 1 };
    for (INT imu = 0; imu < 2; ++imu)
    {
        SCHAR dirs[6][3] =
        {
            {static_cast<SCHAR>(4), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(4), static_cast<SCHAR>(mu[imu] + 1)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(4)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(-4)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(-4), static_cast<SCHAR>(mu[imu] + 1)},
            {static_cast<SCHAR>(-4), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1)},
        };

        INT iTau[6] = { 1, 1, 1, -1, -1, -1 };

        BYTE contributionOf[6][4] =
        {
            {1, 0, 0, 4},
            {0, 1, 0, 4},
            {0, 0, 1, 4},
            {0, 0, 4, 3},
            {0, 4, 2, 4},
            {4, 2, 0, 4},
        };

        for (INT pathidx = 0; pathidx < 6; ++pathidx)
        {
            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                if (4 == contributionOf[pathidx][iSeperation])
                {
                    continue;
                }

                SCHAR L[3] = { 0, 0, 0 };
                SCHAR R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(dirs[pathidx], iSeperation, L, R, LLength, RLength);

                _LAUNCH_KERNEL(_kernelDFermionKSForce_PR_XYTermT_EM TMPARG(deviceVector, deviceGauge), block, threads,
                    pGaugeBuffer,
                    pPhaseBuffer,
                    pForce,
                    appGetLattice()->m_pIndexCache->m_pEtaMu,
                    pRationalFields,
                    pNumerator,
                    uiRationApproxOrder,
                    byFieldId,
                    byGaugeFieldId,
                    fOmega,
                    bShiftCenter ? F(0.5) : F(0.0),
                    fCharge,
                    static_cast<BYTE>(imu), iTau[pathidx],
                    L[0], L[1], L[2], LLength,
                    R[0], R[1], R[2], RLength,
                    contributionOf[pathidx][iSeperation]
                    );
            }
        }
    }

#pragma endregion

#pragma region Polarization term

    //===========================
    //polarization terms
    //ilinkType is +-x +-y +t,
    //INT linkTypes[4][3] =
    //{
    //    {1, 2, 4},
    //    {1, 2, -4},
    //    {-1, 2, 4},
    //    {-1, 2, -4}
    //};
    SCHAR linkTypes[4][3] =
    {
        {1, 2, 4},
        {1, -2, 4},
        {-1, 2, 4},
        {-1, -2, 4}
    };

    for (INT ilinkType = 0; ilinkType < 4; ++ilinkType)
    {
        SCHAR sixlinks[6][3] =
        {
            {linkTypes[ilinkType][0], linkTypes[ilinkType][1], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][0], linkTypes[ilinkType][2], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][0], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][2], linkTypes[ilinkType][0]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][0], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][1], linkTypes[ilinkType][0]}
        };

        for (INT isixtype = 0; isixtype < 6; ++isixtype)
        {
            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                SCHAR L[3] = { 0, 0, 0 };
                SCHAR R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(sixlinks[isixtype], iSeperation, L, R, LLength, RLength);

                const UBOOL bHasLeft = (LLength > 0) && (L[0] > 0);
                const UBOOL bHasRight = (RLength > 0) && (R[0] > 0);

                if (bHasLeft || bHasRight)
                {
                    _LAUNCH_KERNEL(_kernelDFermionKSForce_PR_XYTau_TermT_EM TMPARG(deviceVector, deviceGauge), block, threads,
                        pGaugeBuffer,
                        pPhaseBuffer,
                        pForce,
                        pRationalFields,
                        pNumerator,
                        uiRationApproxOrder,
                        byFieldId,
                        byGaugeFieldId,
                        fOmega,
                        fCharge,
                        L[0], L[1], L[2], LLength,
                        R[0], R[1], R[2], RLength
                        );
                }
            }
        }
    }

#pragma endregion
}

#pragma endregion

#pragma endregion


template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernelR<deviceVector, deviceGauge, vectorN>::DerivateD0_REMOnEvenOdd(
    DOUBLE fOmega,
    UBOOL bShiftCenter,
    Real fCharge,
    const deviceVector* pFermion,
    BYTE byFieldId,
    deviceGauge* pForce,
    const deviceGauge* pGaugeBuffer,
    const Real* pPhase,
    BYTE byGaugeFieldId,
    const deviceVector* const* pRationalFields,
    const Real* pNumerator,
    UINT uiRationApproxOrder)
{
    preparethread;

#pragma region X Y Term
    SCHAR mu[2] = { 0, 1 };
    for (INT imu = 0; imu < 2; ++imu)
    {
        SCHAR dirs[6][3] =
        {
            {static_cast<SCHAR>(4), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(4), static_cast<SCHAR>(mu[imu] + 1)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(4)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(-4)},
            {static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(-4), static_cast<SCHAR>(mu[imu] + 1)},
            {static_cast<SCHAR>(-4), static_cast<SCHAR>(mu[imu] + 1), static_cast<SCHAR>(mu[imu] + 1)},
        };

        INT iTau[6] = { 1, 1, 1, -1, -1, -1 };

        BYTE contributionOf[6][4] =
        {
            {1, 0, 0, 4},
            {0, 1, 0, 4},
            {0, 0, 1, 4},
            {0, 0, 4, 3},
            {0, 4, 2, 4},
            {4, 2, 0, 4},
        };

        for (INT pathidx = 0; pathidx < 6; ++pathidx)
        {
            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                if (4 == contributionOf[pathidx][iSeperation])
                {
                    continue;
                }

                SCHAR L[3] = { 0, 0, 0 };
                SCHAR R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(dirs[pathidx], iSeperation, L, R, LLength, RLength);

                _LAUNCH_KERNEL(_kernelDFermionKSForce_PR_XYTermT_EM_EO TMPARG(deviceVector, deviceGauge), block, threads,
                    pGaugeBuffer,
                    pPhase,
                    pForce,
                    appGetLattice()->m_pIndexCache->m_pEtaMu,
                    pRationalFields,
                    pNumerator,
                    uiRationApproxOrder,
                    byFieldId,
                    byGaugeFieldId,
                    fOmega,
                    bShiftCenter ? F(0.5) : F(0.0),
                    fCharge,
                    static_cast<BYTE>(imu), iTau[pathidx],
                    L[0], L[1], L[2], LLength,
                    R[0], R[1], R[2], RLength,
                    contributionOf[pathidx][iSeperation]
                    );
            }
        }
    }

#pragma endregion

#pragma region Polarization term

    //===========================
    //polarization terms
    //ilinkType is +-x +-y +t,
    SCHAR linkTypes[4][3] =
    {
        {1, 2, 4},
        {1, -2, 4},
        {-1, 2, 4},
        {-1, -2, 4}
    };

    for (INT ilinkType = 0; ilinkType < 4; ++ilinkType)
    {
        SCHAR sixlinks[6][3] =
        {
            {linkTypes[ilinkType][0], linkTypes[ilinkType][1], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][0], linkTypes[ilinkType][2], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][0], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][2], linkTypes[ilinkType][0]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][0], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][1], linkTypes[ilinkType][0]}
        };

        for (INT isixtype = 0; isixtype < 6; ++isixtype)
        {
            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                SCHAR L[3] = { 0, 0, 0 };
                SCHAR R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(sixlinks[isixtype], iSeperation, L, R, LLength, RLength);

                const UBOOL bHasLeft = (LLength > 0) && (L[0] > 0);
                const UBOOL bHasRight = (RLength > 0) && (R[0] > 0);

                if (bHasLeft || bHasRight)
                {
                    _LAUNCH_KERNEL(_kernelDFermionKSForce_PR_XYTau_TermT_EM_EO TMPARG(deviceVector, deviceGauge), block, threads,
                        pGaugeBuffer,
                        pPhase,
                        pForce,
                        pRationalFields,
                        pNumerator,
                        uiRationApproxOrder,
                        byFieldId,
                        byGaugeFieldId,
                        fOmega,
                        fCharge,
                        L[0], L[1], L[2], LLength,
                        R[0], R[1], R[2], RLength
                        );
                }
            }
        }
    }

#pragma endregion
}

template class CFieldFermionKSTKernelR<CLGComplex, CLGComplex, 1>;
template class CFieldFermionKSTKernelR<deviceSU2Vector, deviceSU2, 2>;
template class CFieldFermionKSTKernelR<deviceSU3Vector, deviceSU3, 3>;
template class CFieldFermionKSTKernelR<deviceSU4Vector, deviceSU4, 4>;
//template class CFieldFermionKSTKernelR<deviceSU5Vector, deviceSU5, 5>;
//template class CFieldFermionKSTKernelR<deviceSU6Vector, deviceSU6, 6>;
//template class CFieldFermionKSTKernelR<deviceSU7Vector, deviceSU7, 7>;
//template class CFieldFermionKSTKernelR<deviceSU8Vector, deviceSU8, 8>;

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
