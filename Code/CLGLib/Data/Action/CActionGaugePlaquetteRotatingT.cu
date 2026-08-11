//=============================================================================
// FILENAME : CActionGaugePlaquetteRotatingT.cu
// 
// DESCRIPTION:
// This is the class for rotating su3
//
// REVISION:
//  [mm/dd/yy]
//  [07/08/2024 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CActionGaugePlaquetteRotatingT.h"

__BEGIN_NAMESPACE

#pragma region kernels

#pragma region Clover terms

/**
* Using plaqutte and (f(n)+f(n+mu)+f(n+nu)+f(n+mu+nu))/4 
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermT(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
    DOUBLE betaOverN, DOUBLE fOmegaSq, INT matrixN,
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
    for (BYTE idx0 = 0; idx0 < 3; ++idx0)
    {
        //i=0: 12
        //  1: 13
        //  2: 14
        //  3: 23
        //  4: 24
        //  5: 34
        //0->0, 1->1, 2->3
        //0-> r^2, 1->y^2, 2(or 3)-> x^2
        const BYTE idx = (2 == idx0) ? (idx0 + 1) : idx0;

        //Real resThisThread = F(0.0);

        //========================================
        //find plaqutte 1-4, or 2-4, or 3-4
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

        //0 -> xy, 1 -> xz, 2 -> yz
        //x x y
        const BYTE mushift = (idx0 / 2);
        //y z z
        const BYTE nushift = ((idx0 + 1) / 2) + 1;
        res += static_cast<DOUBLE>(betaOverN * fOmegaSq * (matrixN - _retr(toAdd)) * _deviceFi(byFieldId, sSite4, uiN, idx0, mushift, nushift));
    }

    results[uiSiteIndex] = res;
}


/**
*
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermT_XY(
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
    BYTE idx[6] = { 1, 0, 2, 0, 1, 2};
    BYTE byOtherDir[6] = {2, 1, 2, 0, 0, 1};

    #pragma unroll
    for (BYTE idir = 0; idir < 3; ++idir)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, idir))
        {
            continue;
        }
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        //for xz, yz, i=1,2
        //for xy, i = 0

        //idir = 0: xz, xy with coefficient 1 0
        //idir = 1: yz, yx with coefficient 2 0
        //idir = 2: zx, zy with coefficient 1 2

        deviceGauge stap(_deviceStapleTermGfactorT(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
            idir, 
            byOtherDir[2 * idir],
            idx[2 * idir]));
        _add(stap, _deviceStapleTermGfactorT(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
            idir,
            byOtherDir[2 * idir + 1],
            idx[2 * idir + 1]));

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
_kernelAddChairTermT_Term12(
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
    //Multi-GPU: coordinates entering the Omega factors must be GLOBAL
    //(P4-1.3/1.4-R1); identity on single-GPU. sSite4 itself stays local for
    //the index-table lookups in _deviceChairTermT.
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));
    const Real fXOmega = -(sSite4G.x - _DC_Centerx) * fOmega;

    //===============
    //-x Omega V412
    const Real fV412 = fXOmega * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 3, 0, 1, uiN);

    //===============
    //-x Omega V432
    const Real fV432 = fXOmega * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 3, 2, 1, uiN);

    results[uiSiteIndex] = (fV412 + fV432) * betaOverN;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermT_Term34(
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
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));
    const Real fYOmega = (sSite4G.y - _DC_Centery) * fOmega;

    //===============
    //+y Omega V421
    const Real fV421 = fYOmega * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 3, 1, 0, uiN);

    //===============
    //+y Omega V431
    const Real fV431 = fYOmega * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 3, 2, 0, uiN);

    results[uiSiteIndex] = (fV421 + fV431) * betaOverN;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermT_Term5(
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
    const Real fXYOmega2 = -(sSite4G.x - _DC_Centerx) * (sSite4G.y - _DC_Centery) * fOmegaSq;

    //===============
    //-Omega^2 xy V132
    const Real fV132 = fXYOmega2 * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);

    results[uiSiteIndex] = fV132 * betaOverN;
}


#pragma endregion

#pragma region Chair force

/**
* Split to 15 functions to avoid max-regcount
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term1(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceGauge * __restrict__ pDeviceData,
    deviceGauge *pForceData,
    DOUBLE betaOverN, DOUBLE fOmega
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //-x Omega V412
    //add force for dir=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 3))
    {
        deviceGauge staple_term1_4 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 0, 1, bTorus ? _deviceHi0T : _deviceHi0);

        _mul(staple_term1_4, betaOverN);
        _add(pForceData[uiLink4], staple_term1_4);
    }

    //===============
    //+x Omega V412
    //add force for dir=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 1))
    {
        deviceGauge staple_term1_2 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            1, 0, 3, bTorus ? _deviceHi0T : _deviceHi0);

        _mul(staple_term1_2, betaOverN);
        _add(pForceData[uiLink2], staple_term1_2);
    }

    //===============
    //+x Omega V412
    //add force for dir=x
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 0))
    {
        deviceGauge staple_term1_1 = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 0, 1, bTorus ? _deviceHi0T : _deviceHi0);

        _mul(staple_term1_1, betaOverN);
        _add(pForceData[uiLink1], staple_term1_1);
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term2(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceGauge * __restrict__ pDeviceData,
    deviceGauge *pForceData,
    DOUBLE betaOverN, DOUBLE fOmega
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //-x Omega V432
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 3))
    {
        deviceGauge staple_term2_4 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 2, 1, bTorus ? _deviceHi0T : _deviceHi0);

        _mul(staple_term2_4, betaOverN);
        _add(pForceData[uiLink4], staple_term2_4);
    }

    //===============
    //+x Omega V432
    //add force for mu=4
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 1))
    {
        deviceGauge staple_term2_2 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            1, 2, 3, bTorus ? _deviceHi0T : _deviceHi0);

        _mul(staple_term2_2, betaOverN);
        _add(pForceData[uiLink2], staple_term2_2);
    }

    //===============
    //+x Omega V432
    //add force for mu=4
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 2))
    {
        deviceGauge staple_term2_3 = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 2, 1, bTorus ? _deviceHi0T : _deviceHi0);

        _mul(staple_term2_3, betaOverN);
        _add(pForceData[uiLink3], staple_term2_3);
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term3(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceGauge * __restrict__ pDeviceData,
    deviceGauge *pForceData,
    DOUBLE betaOverN, DOUBLE fOmega
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+y Omega V421
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 3))
    {
        deviceGauge staple_term3_4 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 1, 0, bTorus ? _deviceHi1T : _deviceHi1);

        _mul(staple_term3_4, betaOverN);
        _add(pForceData[uiLink4], staple_term3_4);
    }

    //===============
    //-y Omega V421
    //add force for mu=4
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 0))
    {
        deviceGauge staple_term3_1 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            0, 1, 3, bTorus ? _deviceHi1T : _deviceHi1);

        _mul(staple_term3_1, betaOverN);
        _add(pForceData[uiLink1], staple_term3_1);
    }

    //===============
    //-y Omega V421
    //add force for mu=4
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 1))
    {
        deviceGauge staple_term3_2 = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 1, 0, bTorus ? _deviceHi1T : _deviceHi1);

        _mul(staple_term3_2, betaOverN);
        _add(pForceData[uiLink2], staple_term3_2);
    }

}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term4(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceGauge * __restrict__ pDeviceData,
    deviceGauge *pForceData,
    DOUBLE betaOverN, DOUBLE fOmega
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+y Omega V431
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 3))
    {
        deviceGauge staple_term4_4 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 2, 0, bTorus ? _deviceHi1T : _deviceHi1);

        _mul(staple_term4_4, betaOverN);
        _add(pForceData[uiLink4], staple_term4_4);
    }

    //===============
    //-y Omega V431
    //add force for mu=4
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 0))
    {
        deviceGauge staple_term4_1 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            0, 2, 3, bTorus ? _deviceHi1T : _deviceHi1);

        _mul(staple_term4_1, betaOverN);
        _add(pForceData[uiLink1], staple_term4_1);
    }

    //===============
    //-y Omega V431
    //add force for mu=4
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 2))
    {
        deviceGauge staple_term4_3 = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 2, 0, bTorus ? _deviceHi1T : _deviceHi1);

        _mul(staple_term4_3, betaOverN);
        _add(pForceData[uiLink3], staple_term4_3);
    }

}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term5(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceGauge * __restrict__ pDeviceData,
    deviceGauge *pForceData,
    DOUBLE betaOverN, DOUBLE fOmegaSq
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmegaSq * F(0.125);

    //===============
    //-Omega^2 xy V132
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 0))
    {
        deviceGauge staple_term5_1 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            0, 2, 1, bTorus ? _deviceHi2T : _deviceHi2);

        _mul(staple_term5_1, betaOverN);
        _add(pForceData[uiLink1], staple_term5_1);
    }

    //===============
    //+Omega^2 xy V132
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 1))
    {
        deviceGauge staple_term5_2 = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            1, 2, 0, bTorus ? _deviceHi2T : _deviceHi2);

        _mul(staple_term5_2, betaOverN);
        _add(pForceData[uiLink2], staple_term5_2);
    }

    //===============
    //+Omega^2 xy V132
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 2))
    {
        deviceGauge staple_term5_3 = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            0, 2, 1, bTorus ? _deviceHi2T : _deviceHi2);

        _mul(staple_term5_3, betaOverN);
        _add(pForceData[uiLink3], staple_term5_3);
    }

}

#pragma endregion

#pragma region Projective plane

#pragma region Clover

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermT_Shifted(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE betaOverN, DOUBLE fOmegaSq,
    DOUBLE* results
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //Multi-GPU: global coordinates for the Omega factors (P4-1.3/1.4-R1);
    //identity on single-GPU.
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));
    DOUBLE fXSq = (sSite4G.x - _DC_Centerx + 0.5);
    fXSq = fXSq * fXSq;
    DOUBLE fYSq = (sSite4G.y - _DC_Centery + 0.5);
    fYSq = fYSq * fYSq;

    //======================================================
    //4-plaqutte terms
    //Omega^2 x^2 Retr[1 - U_2,3]
    const DOUBLE fU23 = fXSq * _device4PlaqutteTermT(pDeviceData, 1, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 y^2 Retr[1 - U_1,3]
    const DOUBLE fU13 = fYSq * _device4PlaqutteTermT(pDeviceData, 0, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 (x^2 + y^2) Retr[1 - U_1,2]
    const DOUBLE fU12 = (fXSq + fYSq) * _device4PlaqutteTermT(pDeviceData, 0, 1, uiBigIdx, sSite4, byFieldId);

    results[uiSiteIndex] = (fU23 + fU13 + fU12) * betaOverN * fOmegaSq;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermT_XYZ_Shifted(
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
    BYTE idx[6] = { 1, 0, 2, 0, 1, 2 };
    BYTE byOtherDir[6] = { 2, 1, 2, 0, 0, 1 };

    #pragma unroll
    for (UINT idir = 0; idir < 3; ++idir)
    {
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        //mu = idir, nu = 4, i = mu
        deviceGauge stap(_deviceStapleTermGfactorT(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
            idir,
            byOtherDir[2 * idir],
            idx[2 * idir],
            TRUE));
        _add(stap, _deviceStapleTermGfactorT(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
            idir,
            byOtherDir[2 * idir + 1],
            idx[2 * idir + 1],
            TRUE));
        
        _mul(stap, betaOverN);
        _add(pForceData[linkIndex], stap);
    }
}

#pragma endregion

#pragma region Chair energy

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermT_Term12_Shifted(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE betaOverN, DOUBLE fOmega,
    DOUBLE* results
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -0.125 * betaOverN;
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));
    const DOUBLE fXOmega = (sSite4G.x - _DC_Centerx + 0.5) * fOmega;

    //===============
    //- x Omega V412
    const DOUBLE fV412 = fXOmega * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 3, 0, 1, uiN);

    //===============
    //- x Omega V432
    const DOUBLE fV432 = fXOmega * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 3, 2, 1, uiN);

    results[uiSiteIndex] = (fV412 + fV432) * betaOverN;
}


template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermT_Term34_Shifted(
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
    const DOUBLE fYOmega = (sSite4G.y - _DC_Centery + 0.5) * fOmega;

    //===============
    //+ y Omega V421
    const DOUBLE fV421 = fYOmega * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 3, 1, 0, uiN);

    //===============
    //+ y Omega V431
    const DOUBLE fV431 = fYOmega * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 3, 2, 0, uiN);

    results[uiSiteIndex] = (fV421 + fV431) * betaOverN;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermT_Term5_Shifted(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    DOUBLE betaOverN, DOUBLE fOmegaSq,
    DOUBLE* results
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -0.125 * betaOverN;
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));
    const DOUBLE fXYOmega2 = (sSite4G.x - _DC_Centerx + 0.5) * (sSite4G.y - _DC_Centery + 0.5) * fOmegaSq;

    //===============
    //-Omega^2 xy V142
    const DOUBLE fV132 = fXYOmega2 * _deviceChairTermT(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);

    results[uiSiteIndex] = fV132 * betaOverN;
}

#pragma endregion

#pragma region Chair force

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term1_Shifted(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pForceData,
    DOUBLE betaOverN, DOUBLE fOmega
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //-x Omega V412
    //add force for dir=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    deviceGauge staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 0, 1, bTorus ? _deviceHiShiftedT0 : _deviceHiShifted0);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink4], staple_term);

    //===============
    //+x Omega V412
    //add force for dir=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        1, 0, 3, bTorus ? _deviceHiShiftedT0 : _deviceHiShifted0);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink2], staple_term);

    //===============
    //+x Omega V412
    //add force for dir=x
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    staple_term = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 0, 1, bTorus ? _deviceHiShiftedT0 : _deviceHiShifted0);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink1], staple_term);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term2_Shifted(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pForceData,
    DOUBLE betaOverN, DOUBLE fOmega
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //-x Omega V432
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    deviceGauge staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 2, 1, bTorus ? _deviceHiShiftedT0 : _deviceHiShifted0);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink4], staple_term);

    //===============
    //+x Omega V432
    //add force for mu=4
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        1, 2, 3, bTorus ? _deviceHiShiftedT0 : _deviceHiShifted0);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink2], staple_term);

    //===============
    //+x Omega V432
    //add force for mu=4
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    staple_term = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 2, 1, bTorus ? _deviceHiShiftedT0 : _deviceHiShifted0);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink3], staple_term);
}


template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term3_Shifted(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pForceData,
    DOUBLE betaOverN, DOUBLE fOmega
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+ y Omega V421
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    deviceGauge staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 1, 0, bTorus ? _deviceHiShiftedT1 : _deviceHiShifted1);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink4], staple_term);

    //===============
    //+ y Omega V421
    //add force for mu=1
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 1, 3, bTorus ? _deviceHiShiftedT1 : _deviceHiShifted1);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink1], staple_term);

    //===============
    //+ y Omega V421
    //add force for mu=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    staple_term = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 1, 0, bTorus ? _deviceHiShiftedT1 : _deviceHiShifted1);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink2], staple_term);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term4_Shifted(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pForceData,
    DOUBLE betaOverN, DOUBLE fOmega
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+ y Omega V431
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    deviceGauge staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 2, 0, bTorus ? _deviceHiShiftedT1 : _deviceHiShifted1);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink4], staple_term);

    //===============
    //+ y Omega V431
    //add force for mu=4
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 2, 3, bTorus ? _deviceHiShiftedT1 : _deviceHiShifted1);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink1], staple_term);

    //===============
    //+ y Omega V431
    //add force for mu=3
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    staple_term = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 2, 0, bTorus ? _deviceHiShiftedT1 : _deviceHiShifted1);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink3], staple_term);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermT_Term5_Shifted(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pForceData,
    DOUBLE betaOverN, DOUBLE fOmegaSq
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * F(0.5) * fOmegaSq * F(0.125);

    //===============
    //- Omega^2 xy V132
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    deviceGauge staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 2, 1, bTorus ? _deviceHiShiftedT2 : _deviceHiShifted2);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink1], staple_term);

    //===============
    //- Omega^2 xy V132
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        1, 2, 0, bTorus ? _deviceHiShiftedT2 : _deviceHiShifted2);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink2], staple_term);

    //===============
    //- Omega^2 xy V132
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    staple_term = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 2, 1, bTorus ? _deviceHiShiftedT2 : _deviceHiShifted2);

    _mul(staple_term, betaOverN);
    _add(pForceData[uiLink3], staple_term);
}

#pragma endregion

#pragma endregion

#pragma endregion

template<typename deviceGauge, INT matrixN>
void CActionGaugePlaquetteRotatingT<deviceGauge, matrixN>::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        m_fLastEnergy = EnergySingleField(FALSE, pGauge, NULL);
    }
}

template<typename deviceGauge, INT matrixN>
void CActionGaugePlaquetteRotatingT<deviceGauge, matrixN>::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    DOUBLE fOmega = 0.1;
    param.FetchValueDOUBLE(_T("Omega"), fOmega);
    m_fOmega = fOmega;

    DOUBLE fXi = 1.0;
    param.FetchValueDOUBLE(_T("Xi"), fXi);
    m_fXi = fXi;

    INT iUsing4Plaq = 0;
    if (param.FetchValueINT(_T("CloverEnergy"), iUsing4Plaq))
    {
        if (1 == iUsing4Plaq)
        {
            m_bCloverEnergy = TRUE;
        }
    }

    INT iShiftCoord = 0;
    param.FetchValueINT(_T("ShiftCoord"), iShiftCoord);
    m_bShiftHalfCoord = (0 != iShiftCoord);

    INT iTorus = 0;
    param.FetchValueINT(_T("Torus"), iTorus);
    m_bTorus = (0 != iTorus);

    // The rotating (S1/S2) terms assume the Wilson beta = 6/g^2 convention, while
    // improved gauge fields feed beta = 10/g^2. Rescale beta for the rotating terms
    // only; the S0 plaquette/rectangle term keeps the raw beta (handled by the field).
    CCString sBetaScale = _T("ERBSF_Naive");
    param.FetchStringValue(_T("RotationBetaScale"), sBetaScale);
    m_eBetaScaleType = __STRING_TO_ENUM(ERotationBetaScaleFactor, sBetaScale);

    switch (m_eBetaScaleType)
    {
    case ERBSF_TreeImprove:
    case ERBSF_OneLoopImprove:
        // beta = 10/g^2 -> 6/g^2, factor 6/10
        m_fRotationScaleFactor = 0.6;
        break;
    case ERBSF_Custom:
        {
            DOUBLE fScaleFactor = 1.0;
            param.FetchValueDOUBLE(_T("ScaleFactor"), fScaleFactor);
            m_fRotationScaleFactor = fScaleFactor;
        }
        break;
    case ERBSF_Naive:
    default:
        m_fRotationScaleFactor = 1.0;
        break;
    }
}

template<typename deviceGauge, INT matrixN>
UBOOL CActionGaugePlaquetteRotatingT<deviceGauge, matrixN>::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    if (abs(m_fXi - 1.0) > _CLG_FLT_EPSILON)
    {
        if (NULL != pStaple)
        {
            appCrucial(_T("CActionGaugePlaquetteRotatingT<deviceGauge, matrixN> do not support staple now.\n"));
        }
        pGauge->CalculateForceAnisotropy(pForce, m_fBetaOverN, m_fXi);
    }
    else
    {
        pGauge->CalculateForceAndStaple(pForce, pStaple, static_cast<Real>(m_fBetaOverN));
    }

    const CFieldGaugeLink<deviceGauge, matrixN>* pGaugeSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(pGauge);
    CFieldGaugeLink<deviceGauge, matrixN>* pForceSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pForce);
    if (NULL == pGaugeSU3 || NULL == pForceSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotatingT<deviceGauge, matrixN> only work with SU3 now.\n"));
        return TRUE;
    }

    const CBoundaryCondition* pBC = appGetLattice()->m_pIndex->GetBoudanryCondition();

    if (m_bTorus)
    {
        CalculateForceOnGaugeTorus(pGaugeSU3, pForceSU3);
    }
    else if (NULL != dynamic_cast<const CBoundaryConditionProjectivePlaneSquare*>(pBC))
    {
        CalculateForceOnGaugeProjectivePlane(pGaugeSU3, pForceSU3);
    }
    else
    {
        CalculateForceOnGaugeDirichlet(pGaugeSU3, pForceSU3);
    }

    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
template<typename deviceGauge, INT matrixN>
DOUBLE CActionGaugePlaquetteRotatingT<deviceGauge, matrixN>::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

    const CFieldGaugeLink<deviceGauge, matrixN>* pGaugeSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(pGauge);
    if (NULL == pGaugeSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotatingT<deviceGauge, matrixN> only work with SU3 now.\n"));
        return m_fNewEnergy;
    }

    const CBoundaryCondition* pBC = appGetLattice()->m_pIndex->GetBoudanryCondition();

    if (m_bTorus)
    {
        EnergyTorus(pGaugeSU3);
    }
    else if (NULL != dynamic_cast<const CBoundaryConditionProjectivePlaneSquare*>(pBC))
    {
        EnergyProjectivePlane(pGaugeSU3);
    }
    else
    {
        EnergyDirichlet(pGaugeSU3);
    }

    //m_fNewEnergy is assembled from local partial ThreadBufferSums
    //(plaquette part included); sum the total across ranks for the global
    //action (P4-1.3). No-op on a lone rank.
    appGlobalSum(m_fNewEnergy);
    //P4-1.6: S0/S1/S2 are also local partial sums (S0 via the base
    //CalculatePlaqutteEnergy path, S1/S2 via ThreadBufferSum); the
    //CMeasureRotatingAction measurement reports them, so they must be
    //global too. Identity on single-GPU / lone rank.
    appGlobalSum(m_fS0Energy);
    appGlobalSum(m_fS1Energy);
    appGlobalSum(m_fS2Energy);

    return m_fNewEnergy;
}

/**
* In the case of Dirichlet, 
* m_bShiftHalfCoord = FALSE
* m_bCloverEnergy = FALSE
* 
*/
template<typename deviceGauge, INT matrixN>
void CActionGaugePlaquetteRotatingT<deviceGauge, matrixN>::EnergyDirichlet(const class CFieldGaugeLink<deviceGauge, matrixN>* pGaugeSU3)
{
    appAssert(!m_bShiftHalfCoord && !m_bCloverEnergy);

    if (abs(m_fXi - 1.0) > _CLG_FLT_EPSILON)
    {
        m_fS0Energy = pGaugeSU3->CalculatePlaqutteEnergyAnisotropy(m_fBetaOverN, m_fXi);
    }
    else
    {
        m_fS0Energy = pGaugeSU3->CalculatePlaqutteEnergy(m_fBetaOverN);
    }

    preparethread;
    appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

    //xy, xz, yz plaquettes (raw, omega applied later)
    _LAUNCH_KERNEL(_kernelAdd4PlaqutteTermT<deviceGauge>, block, threads,
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGaugeSU3->m_byFieldId],
        GetRotationBetaOverN() * m_fXi,
        1.0,
        matrixN,
        _D_RealThreadBuffer);

    m_fS2Energy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    //412, 432 chair terms (raw, omega applied later)
    _LAUNCH_KERNEL(_kernelAddChairTermT_Term12<deviceGauge>, block, threads,
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        GetRotationBetaOverN(),
        1.0,
        _D_RealThreadBuffer);

    m_fS1Energy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    //421, 431 chair terms (raw, omega applied later)
    _LAUNCH_KERNEL(_kernelAddChairTermT_Term34<deviceGauge>, block, threads,
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        GetRotationBetaOverN(),
        1.0,
        _D_RealThreadBuffer);

    m_fS1Energy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    //132 chair term (raw, omega applied later)
    _LAUNCH_KERNEL(_kernelAddChairTermT_Term5<deviceGauge>, block, threads,
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        GetRotationBetaOverN() * m_fXi,
        1.0,
        _D_RealThreadBuffer);

    m_fS2Energy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    m_fNewEnergy = m_fS0Energy + m_fOmega * m_fS1Energy + m_fOmega * m_fOmega * m_fS2Energy;
}

/**
* In the case of ProjectivePlane,
* m_bShiftHalfCoord = TRUE
* m_bCloverEnergy = TRUE
*
*/
template<typename deviceGauge, INT matrixN>
void CActionGaugePlaquetteRotatingT<deviceGauge, matrixN>::EnergyProjectivePlane(const class CFieldGaugeLink<deviceGauge, matrixN>* pGaugeSU3)
{
    appAssert(m_bShiftHalfCoord && m_bCloverEnergy);
    if (abs(m_fXi - 1.0) > _CLG_FLT_EPSILON)
    {
        m_fS0Energy = pGaugeSU3->CalculatePlaqutteEnergyUseCloverAnisotropy(m_fBetaOverN, m_fXi);
    }
    else
    {
        m_fS0Energy = pGaugeSU3->CalculatePlaqutteEnergyUseClover(m_fBetaOverN);
    }

    preparethread;
    appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

    //xy, xz, yz plaquettes (Omega^2 terms)
    _LAUNCH_KERNEL(_kernelAdd4PlaqutteTermT_Shifted<deviceGauge>, block, threads,
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        GetRotationBetaOverN() * m_fXi,
        1.0,
        _D_RealThreadBuffer);

    m_fS2Energy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    //412, 432 chair terms (Omega terms)
    _LAUNCH_KERNEL(_kernelAddChairTermT_Term12_Shifted<deviceGauge>, block, threads,
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        GetRotationBetaOverN(),
        1.0,
        _D_RealThreadBuffer);
    m_fS1Energy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    //421, 431 chair terms (Omega terms)
    _LAUNCH_KERNEL(_kernelAddChairTermT_Term34_Shifted<deviceGauge>, block, threads,
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        GetRotationBetaOverN(),
        1.0,
        _D_RealThreadBuffer);
    m_fS1Energy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    //132 chair term (Omega^2 term)
    _LAUNCH_KERNEL(_kernelAddChairTermT_Term5_Shifted<deviceGauge>, block, threads,
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        GetRotationBetaOverN() * m_fXi,
        1.0,
        _D_RealThreadBuffer);
    m_fS2Energy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    m_fNewEnergy = m_fS0Energy + m_fOmega * m_fS1Energy + m_fOmega * m_fOmega * m_fS2Energy;
}

template<typename deviceGauge, INT matrixN>
void CActionGaugePlaquetteRotatingT<deviceGauge, matrixN>::EnergyTorus(const class CFieldGaugeLink<deviceGauge, matrixN>* pGaugeSU3)
{
    if (abs(m_fXi - 1.0) > _CLG_FLT_EPSILON)
    {
        m_fS0Energy = pGaugeSU3->CalculatePlaqutteEnergyAnisotropy(m_fBetaOverN, m_fXi);
    }
    else
    {
        m_fS0Energy = pGaugeSU3->CalculatePlaqutteEnergy(m_fBetaOverN);
    }

    preparethread;
    appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

    if (m_bShiftHalfCoord)
    {
        _LAUNCH_KERNEL(_kernelAdd4PlaqutteTermT_Shifted<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            GetRotationBetaOverN() * m_fXi,
            1.0,
            _D_RealThreadBuffer);

        m_fS2Energy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _LAUNCH_KERNEL(_kernelAddChairTermT_Term12_Shifted<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            GetRotationBetaOverN(),
            1.0,
            _D_RealThreadBuffer);
        m_fS1Energy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _LAUNCH_KERNEL(_kernelAddChairTermT_Term34_Shifted<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            GetRotationBetaOverN(),
            1.0,
            _D_RealThreadBuffer);
        m_fS1Energy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _LAUNCH_KERNEL(_kernelAddChairTermT_Term5_Shifted<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            GetRotationBetaOverN() * m_fXi,
            1.0,
            _D_RealThreadBuffer);
        m_fS2Energy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
    }
    else
    {
        _LAUNCH_KERNEL(_kernelAdd4PlaqutteTermT<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGaugeSU3->m_byFieldId],
            GetRotationBetaOverN() * m_fXi,
            1.0,
            matrixN,
            _D_RealThreadBuffer);

        m_fS2Energy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _LAUNCH_KERNEL(_kernelAddChairTermT_Term12<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            GetRotationBetaOverN(),
            1.0,
            _D_RealThreadBuffer);

        m_fS1Energy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _LAUNCH_KERNEL(_kernelAddChairTermT_Term34<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            GetRotationBetaOverN(),
            1.0,
            _D_RealThreadBuffer);

        m_fS1Energy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _LAUNCH_KERNEL(_kernelAddChairTermT_Term5<deviceGauge>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            GetRotationBetaOverN() * m_fXi,
            1.0,
            _D_RealThreadBuffer);

        m_fS2Energy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
    }

    m_fNewEnergy = m_fS0Energy + m_fOmega * m_fS1Energy + m_fOmega * m_fOmega * m_fS2Energy;
}

template<typename deviceGauge, INT matrixN>
void CActionGaugePlaquetteRotatingT<deviceGauge, matrixN>::CalculateForceOnGaugeDirichlet(const class CFieldGaugeLink<deviceGauge, matrixN>* pGaugeSU3, class CFieldGaugeLink<deviceGauge, matrixN>* pForceSU3) const
{
    DOUBLE fOmegaSq = m_fOmega * m_fOmega;
    preparethread;

    //xz, yx, xy plaqutte
    _LAUNCH_KERNEL(_kernelAddForce4PlaqutteTermT_XY<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, GetRotationBetaOverN() * m_fXi, fOmegaSq);

    _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term1<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

    _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term2<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

    _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term3<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

    _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term4<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

    _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term5<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, GetRotationBetaOverN() * m_fXi, fOmegaSq);
}

template<typename deviceGauge, INT matrixN>
void CActionGaugePlaquetteRotatingT<deviceGauge, matrixN>::CalculateForceOnGaugeProjectivePlane(const class CFieldGaugeLink<deviceGauge, matrixN>* pGaugeSU3, class CFieldGaugeLink<deviceGauge, matrixN>* pForceSU3) const
{
    DOUBLE fOmegaSq = m_fOmega * m_fOmega;
    preparethread;
    _LAUNCH_KERNEL(_kernelAddForce4PlaqutteTermT_XYZ_Shifted<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, GetRotationBetaOverN() * m_fXi, fOmegaSq);

    _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term1_Shifted<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

    _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term2_Shifted<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

    _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term3_Shifted<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

    _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term4_Shifted<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

    _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term5_Shifted<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, GetRotationBetaOverN() * m_fXi, fOmegaSq);
}

template<typename deviceGauge, INT matrixN>
void CActionGaugePlaquetteRotatingT<deviceGauge, matrixN>::CalculateForceOnGaugeTorus(const class CFieldGaugeLink<deviceGauge, matrixN>* pGaugeSU3, class CFieldGaugeLink<deviceGauge, matrixN>* pForceSU3) const
{
    preparethread;
    DOUBLE fOmegaSq = m_fOmega * m_fOmega;
    if (m_bShiftHalfCoord)
    {
        _LAUNCH_KERNEL(_kernelAddForce4PlaqutteTermT_XYZ_Shifted<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, GetRotationBetaOverN() * m_fXi, fOmegaSq);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term1_Shifted<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term2_Shifted<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term3_Shifted<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term4_Shifted<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term5_Shifted<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, GetRotationBetaOverN() * m_fXi, fOmegaSq);
    }
    else
    {
        _LAUNCH_KERNEL(_kernelAddForce4PlaqutteTermT_XY<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, GetRotationBetaOverN() * m_fXi, fOmegaSq);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term1<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term2<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term3<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term4<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, GetRotationBetaOverN(), m_fOmega);

        _LAUNCH_KERNEL(_kernelAddForceChairTermT_Term5<deviceGauge>, block, threads, pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, GetRotationBetaOverN() * m_fXi, fOmegaSq);
    }
}

template<typename deviceGauge, INT matrixN>
CCString CActionGaugePlaquetteRotatingT<deviceGauge, matrixN>::GetInfos(const CCString &tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("Omega : ") + appToString(m_fOmega) + _T("\n");
    sRet = sRet + tab + _T("Xi : ") + appToString(m_fXi) + _T("\n");

    sRet = sRet + tab + _T("ShiftCenter : ") + (m_bShiftHalfCoord ? _T("1") : _T("0")) + _T("\n");
    sRet = sRet + tab + _T("Clover : ") + (m_bCloverEnergy ? _T("1") : _T("0")) + _T("\n");
    sRet = sRet + tab + _T("Torus : ") + (m_bTorus ? _T("1") : _T("0")) + _T("\n");

    sRet = sRet + tab + _T("RotationBetaScale : ") + __ENUM_TO_STRING(ERotationBetaScaleFactor, m_eBetaScaleType) + _T("\n");
    sRet = sRet + tab + _T("RotationScaleFactor : ") + appToString(m_fRotationScaleFactor) + _T("\n");

    return sRet;
}

template class CActionGaugePlaquetteRotatingT<CLGComplex, 1>;
template class CActionGaugePlaquetteRotatingT<deviceSU2, 2>;
template class CActionGaugePlaquetteRotatingT<deviceSU3, 3>;
template class CActionGaugePlaquetteRotatingT<deviceSU4, 4>;

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotatingU1)
__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotatingSU2)
__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotating)
__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotatingSU4)
//__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotatingSU5)
//__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotatingSU6)
//__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotatingSU7)
//__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotatingSU8)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================