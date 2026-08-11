//=============================================================================
// FILENAME : CActionGaugePlaquetteRotatingT3D.cu
// 
// DESCRIPTION:
// This is the class for rotating in the case of 3D
// 
// see:
// https://cboard.cprogramming.com/cplusplus-programming/113400-gcc-template-class-child-cant-directly-access-parent-fields.html
//
// REVISION:
//  [mm/dd/yy]
//  [27/10/2022 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CActionGaugePlaquetteRotatingT3D.h"

__BEGIN_NAMESPACE

#pragma region kernels

#pragma region Clover terms

/**
* Using plaqutte and (f(n)+f(n+mu)+f(n+nu)+f(n+mu+nu))/4 
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermT3D(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
    DOUBLE betaOverN, DOUBLE fOmegaSq,
    DOUBLE* results
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);
#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLength = __idx->m_pSmallData[CIndexData::kPlaqLengthIdx];
    const UINT plaqCountAllSite = __idx->m_pSmallData[CIndexData::kPlaqPerSiteIdx] * plaqLength;
#endif

    DOUBLE res = 0.0;
    #pragma unroll
    for (BYTE idx0 = 0; idx0 < 2; ++idx0)
    {
        //i=0: 12
        //  1: 13
        //  2: 23
        //0->1, 1->2
        //0-> x^2, 1->y^2
        const BYTE idx = idx0 + 1;

        //========================================
        //find plaqutte 1-3, or 2-3
        SIndex first = pCachedPlaqutte[idx * plaqLength + uiSiteIndex * plaqCountAllSite];
        deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }
        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedPlaqutte[idx * plaqLength + j + uiSiteIndex * plaqCountAllSite];
            deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
            if (first.NeedToDagger())
            {
                _muldag(toAdd, toMul);
            }
            else
            {
                _mul(toAdd, toMul);
            }
        }

        //0 -> xz, 1 -> yz
        const BYTE mushift = idx0;
        const BYTE nushift = 2;
        res += static_cast<DOUBLE>(betaOverN * fOmegaSq * (3.0 - _retr(toAdd)) * _deviceFi(byFieldId, sSite4, uiN, 2 - idx0, mushift, nushift));
    }

    results[uiSiteIndex] = res;
}


/**
*
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermT_XY3D(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pForceData,
    DOUBLE betaOverN,
    DOUBLE fOmegaSq
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(-0.5);

    //idx is the coefficient, 0 -> r^2, 1 -> y^2, 2 -> x^2
    BYTE idx[4] = {2, 1, 2, 1};
    //byOtherDir is the other direction of the staple
    BYTE byOtherDir[4] = {2, 2, 0, 1};
    //deviceGauge plaqSum = deviceGauge::makeSU3Zero();
    #pragma unroll
    for (BYTE idir = 0; idir < 3; ++idir)
    {

        if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, idir))
        {
            continue;
        }
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        deviceGauge stap(_deviceStapleTermGfactorT(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
            idir, 
            byOtherDir[idir],
            idx[idir]));
        if (2 == idir)
        {
            _add(stap, _deviceStapleTermGfactorT(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
                idir,
                byOtherDir[idir + 1],
                idx[idir + 1]));
        }

        _mul(stap, betaOverN);
        _add(pForceData[linkIndex], stap);
    }
}

#pragma endregion

#pragma region Chair Energy

/**
* Split into 3 functions to avoid max-register problem
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermT_Term1234_3D(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pDeviceData,
    DOUBLE betaOverN, DOUBLE fOmega,
    DOUBLE* results
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        return;
    }

    betaOverN = F(0.125) * betaOverN;
    //Multi-GPU: the coordinates entering the Omega factors must be GLOBAL
    //(P4-1.2/1.4-R1); identity on single-GPU. sSite4 itself stays local for
    //the index-table lookups in _deviceChairTermT.
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));
    const Real fXOmega = (sSite4G.x - _DC_Centerx) * fOmega;
    const Real fYOmega = (sSite4G.y - _DC_Centery) * fOmega;

    //===============
    //+x Omega V312
    const Real fV312 = -fXOmega * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 2, 0, 1, uiN);

    //===============
    //+y Omega V321
    const Real fV321 = fYOmega * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 2, 1, 0, uiN);

    results[uiSiteIndex] = (fV312  + fV321) * betaOverN;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermT_Term53D(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pDeviceData,
    DOUBLE betaOverN, DOUBLE fOmegaSq,
    DOUBLE* results
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        return;
    }

    betaOverN = F(0.125) * betaOverN;
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));
    const Real fXYOmega2 = (sSite4G.x - _DC_Centerx) * (sSite4G.y - _DC_Centery) * fOmegaSq;

    //===============
    //+Omega^2 xy V132
    const Real fV132 = fXYOmega2 * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);

    results[uiSiteIndex] = fV132 * betaOverN;
}


#pragma endregion

#pragma region Chair force

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term13D(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pDeviceData,
    deviceGauge *pForceData,
    DOUBLE betaOverN, DOUBLE fOmega
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * fOmega * F(0.0625);

    //===============
    //+x Omega V312
    //add force for dir=3
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 2);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 2))
    {
        deviceGauge staple_term1_4 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            2, 0, 1, 0);

        _mul(staple_term1_4, betaOverN);
        _add(pForceData[uiLink4], staple_term1_4);
    }

    //===============
    //+x Omega V312
    //add force for dir=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 1))
    {
        deviceGauge staple_term1_2 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            1, 0, 2, 0);

        _mul(staple_term1_2, betaOverN);
        _add(pForceData[uiLink2], staple_term1_2);
    }

    //===============
    //+x Omega V312
    //add force for dir=x
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 0))
    {
        deviceGauge staple_term1_1 = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            2, 0, 1, 0);

        _mul(staple_term1_1, betaOverN);
        _add(pForceData[uiLink1], staple_term1_1);
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term33D(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pDeviceData,
    deviceGauge *pForceData,
    DOUBLE betaOverN, DOUBLE fOmega
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * fOmega * F(0.0625);

    //===============
    //+y Omega V421
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 2);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 2))
    {
        deviceGauge staple_term3_4 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            2, 1, 0, _deviceHi1);

        _mul(staple_term3_4, betaOverN);
        _add(pForceData[uiLink4], staple_term3_4);
    }

    //===============
    // y Omega V421
    //add force for mu=4
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 0))
    {
        deviceGauge staple_term3_1 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            0, 1, 2, _deviceHi1);

        _mul(staple_term3_1, betaOverN);
        _add(pForceData[uiLink1], staple_term3_1);
    }

    //===============
    // y Omega V421
    //add force for mu=4
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 1))
    {
        deviceGauge staple_term3_2 = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            2, 1, 0, _deviceHi1);

        _mul(staple_term3_2, betaOverN);
        _add(pForceData[uiLink2], staple_term3_2);
    }

}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term53D(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pDeviceData,
    deviceGauge *pForceData,
    DOUBLE betaOverN, DOUBLE fOmegaSq
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmegaSq * F(0.125);

    //===============
    //+Omega^2 xy V132
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 0))
    {
        deviceGauge staple_term5_1 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            0, 2, 1, _deviceHi2);

        _mul(staple_term5_1, betaOverN);
        _add(pForceData[uiLink1], staple_term5_1);
    }

    //===============
    //+Omega^2 xy V132
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 1))
    {
        deviceGauge staple_term5_2 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            1, 2, 0, _deviceHi2);

        _mul(staple_term5_2, betaOverN);
        _add(pForceData[uiLink2], staple_term5_2);
    }

    //===============
    //+Omega^2 xy V132
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 2))
    {
        deviceGauge staple_term5_3 = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            0, 2, 1, _deviceHi2);

        _mul(staple_term5_3, betaOverN);
        _add(pForceData[uiLink3], staple_term5_3);
    }

}

#pragma endregion

#pragma region Projective plane

#pragma region Clover

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermT_Shifted3D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE betaOverN, DOUBLE fOmegaSq,
    DOUBLE* results
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //Multi-GPU: global coordinates for the Omega factors (P4-1.2/1.4-R1);
    //identity on single-GPU.
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));
    DOUBLE fXSq = (sSite4G.x - _DC_Centerx + 0.5);
    fXSq = fXSq * fXSq;
    DOUBLE fYSq = (sSite4G.y - _DC_Centery + 0.5);
    fYSq = fYSq * fYSq;

    //======================================================
    //4-plaqutte terms
    //Omega^2 x^2 Retr[1 - U_2,3]
    const DOUBLE fU13 = fXSq * _device4PlaqutteTermT(pDeviceData, 0, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 y^2 Retr[1 - U_1,3]
    const DOUBLE fU23 = fYSq * _device4PlaqutteTermT(pDeviceData, 1, 2, uiBigIdx, sSite4, byFieldId);

    results[uiSiteIndex] = (fU23 + fU13) * betaOverN * fOmegaSq;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermT_XYZ_Shifted3D(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pForceData,
    DOUBLE betaOverN, DOUBLE fOmegaSq
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(-0.5);
    //deviceGauge plaqSum = deviceGauge::makeSU3Zero();
    BYTE idx[4] = { 2, 1, 2, 1 };
    BYTE byOtherDir[4] = { 2, 2, 0, 1 };

    #pragma unroll
    for (UINT idir = 0; idir < 3; ++idir)
    {
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        deviceGauge stap(_deviceStapleTermGfactorT(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
            idir,
            byOtherDir[idir],
            idx[idir],
            TRUE));
        if (2 == idir)
        {
            _add(stap, _deviceStapleTermGfactorT(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
                idir,
                byOtherDir[idir + 1],
                idx[idir + 1],
                TRUE));
        }
        
        _mul(stap, betaOverN);
        _add(pForceData[linkIndex], stap);
    }
}

#pragma endregion

#pragma region Chair energy

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermT_Term1234_Shifted3D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE betaOverN, DOUBLE fOmega,
    DOUBLE* results
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = 0.125 * betaOverN;
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));
    const DOUBLE fXOmega = (sSite4G.x - _DC_Centerx + 0.5) * fOmega;
    const DOUBLE fYOmega = (sSite4G.y - _DC_Centery + 0.5) * fOmega;

    //===============
    // F01F12 term  x Omega V312
    const DOUBLE fV312 = -fXOmega * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 2, 0, 1, uiN);
    // F02F12 term  y Omega V321
    const DOUBLE fV321 = fYOmega * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 2, 1, 0, uiN);

    results[uiSiteIndex] = (fV312 + fV321) * betaOverN;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermT_Term5_Shifted3D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE betaOverN, DOUBLE fOmegaSq,
    DOUBLE* results
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = 0.125 * betaOverN;
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));
    const DOUBLE fXYOmega2 = (sSite4G.x - _DC_Centerx + 0.5) * (sSite4G.y - _DC_Centery + 0.5) * fOmegaSq;

    //===============
    //+Omega^2 xy V142
    const DOUBLE fV132 = fXYOmega2 * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);

    results[uiSiteIndex] = fV132 * betaOverN;
}

#pragma endregion

#pragma region Chair force

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term1_Shifted3D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pForceData,
    DOUBLE betaOverN, DOUBLE fOmega
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+x Omega V312
    //add force for dir=3
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 2);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    //{
    deviceGauge staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 0, 1, _deviceHiShifted0);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink4], staple_term);
    //}


    //===============
    //+x Omega V312
    //add force for dir=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    //{
    staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        1, 0, 2, _deviceHiShifted0);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink2], staple_term);
   // }

    //===============
    //+x Omega V312
    //add force for dir=x
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    //{
    staple_term = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 0, 1, _deviceHiShifted0);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink1], staple_term);
    //}
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term3_Shifted3D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pForceData,
    DOUBLE betaOverN, DOUBLE fOmega
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+ y Omega V321
    //add force for mu=3
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 2);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    //{
    deviceGauge staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 1, 0, _deviceHiShifted1);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink4], staple_term);
    //}

    //===============
    //+ y Omega V321
    //add force for mu=1
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    //{
    staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 1, 2, _deviceHiShifted1);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink1], staple_term);
    //}


    //===============
    //+ y Omega V321
    //add force for mu=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    //{
    staple_term = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 1, 0, _deviceHiShifted1);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink2], staple_term);
    //}

}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term5_Shifted3D(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pForceData,
    DOUBLE betaOverN, DOUBLE fOmegaSq
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmegaSq * F(0.125);

    //===============
    //+ Omega^2 xy V132
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    //{
    deviceGauge staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 2, 1, _deviceHiShifted2);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink1], staple_term);
    //}

    //===============
    //+ Omega^2 xy V132
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    //{
    staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        1, 2, 0, _deviceHiShifted2);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink2], staple_term);
    //}

    //===============
    //+ Omega^2 xy V132
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    //{
    staple_term = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 2, 1, _deviceHiShifted2);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink3], staple_term);
    //}
}

#pragma endregion

#pragma endregion

#pragma endregion

template<typename deviceGauge, INT matrixN>
UBOOL CActionGaugePlaquetteRotatingT3D<deviceGauge, matrixN>::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    pGauge->CalculateForceAndStaple(pForce, pStaple, static_cast<Real>(this->GetBetaOverN()));

    const CFieldGaugeLink<deviceGauge, matrixN>* pGaugeSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(pGauge);
    CFieldGaugeLink<deviceGauge, matrixN>* pForceSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pForce);
    if (NULL == pGaugeSU3 || NULL == pForceSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3 now.\n"));
        return TRUE;
    }

    preparethread;
    const DOUBLE beta = this->GetBetaOverN();
    const DOUBLE omega = this->GetOmega();
    const DOUBLE omegasq = omega * omega;

    if (!this->m_bShiftHalfCoord)
    {
        _LAUNCH_KERNEL(_kernelAddForce4PlaqutteTermT_XY3D<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, beta, omegasq);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term13D<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, beta, omega);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term33D<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, beta, omega);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term53D<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, beta, omegasq);
    }
    else
    {

        _LAUNCH_KERNEL(_kernelAddForce4PlaqutteTermT_XYZ_Shifted3D<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, beta, omegasq);
        
        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term1_Shifted3D<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, beta, omega);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term3_Shifted3D<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, beta, omega);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term5_Shifted3D<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, beta, omegasq);
    }

    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
template<typename deviceGauge, INT matrixN>
DOUBLE CActionGaugePlaquetteRotatingT3D<deviceGauge, matrixN>::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple)
{
    //see:
    //https://cboard.cprogramming.com/cplusplus-programming/113400-gcc-template-class-child-cant-directly-access-parent-fields.html
    if (bBeforeEvolution)
    {
        return this->m_fLastEnergy;
    }

    const CFieldGaugeLink<deviceGauge, matrixN>* pGaugeSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(pGauge);
    if (NULL == pGaugeSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3 now.\n"));
        return this->m_fNewEnergy;
    }

    if (this->IsCloverEnergy())
    {
        this->m_fNewEnergy = pGauge->CalculatePlaqutteEnergyUseClover(this->GetBetaOverN());
    }
    else
    {
        this->m_fNewEnergy = pGauge->CalculatePlaqutteEnergy(this->GetBetaOverN());
    }

    preparethread;

    appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

    if (this->m_bShiftHalfCoord)
    {

        _LAUNCH_KERNEL(_kernelAdd4PlaqutteTermT_Shifted3D<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            this->GetBetaOverN(),
            this->GetOmega() * this->GetOmega(),
            _D_RealThreadBuffer);

        this->m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _LAUNCH_KERNEL(_kernelAddChairTermT_Term1234_Shifted3D<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId, 
            pGaugeSU3->m_pDeviceData, 
            this->GetBetaOverN(),
            this->GetOmega(), 
            _D_RealThreadBuffer);

        this->m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _LAUNCH_KERNEL(_kernelAddChairTermT_Term5_Shifted3D<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            this->GetBetaOverN(),
            this->GetOmega() * this->GetOmega(),
            _D_RealThreadBuffer);

        this->m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    }
    else
    {
        _LAUNCH_KERNEL(_kernelAdd4PlaqutteTermT3D<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGaugeSU3->m_byFieldId],
            this->GetBetaOverN(),
            this->GetOmega() * this->GetOmega(),
            _D_RealThreadBuffer);

        this->m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _LAUNCH_KERNEL(_kernelAddChairTermT_Term1234_3D<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            this->GetBetaOverN(),
            this->GetOmega(),
            _D_RealThreadBuffer);

        this->m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _LAUNCH_KERNEL(_kernelAddChairTermT_Term53D<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            this->GetBetaOverN(),
            this->GetOmega() * this->GetOmega(),
            _D_RealThreadBuffer);

        this->m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    }

    //The rotating-term kernels write per-site partial sums into the thread
    //buffer and ThreadBufferSum reduces only this rank's local sub-lattice;
    //the plaquette part (CalculatePlaqutteEnergy*) is likewise local. Sum the
    //total across ranks for the global action (P4-1.2). No-op on a lone rank.
    appGlobalSum(this->m_fNewEnergy);

    return this->m_fNewEnergy;
}

template class CActionGaugePlaquetteRotatingT3D<CLGComplex, 1>;
template class CActionGaugePlaquetteRotatingT3D<deviceSU3, 3>;

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotatingU1_3D)
__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotating3D)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================