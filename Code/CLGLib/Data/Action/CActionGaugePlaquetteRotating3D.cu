//=============================================================================
// FILENAME : CActionGaugePlaquetteRotating3D.cu
// 
// DESCRIPTION:
// This is the class for rotating su3
//
// REVISION:
//  [27/10/2022 nbale]
//=============================================================================
#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotating3D)

#pragma region kernels

#pragma region Clover terms

/**
* Using plaqutte and (f(n)+f(n+mu)+f(n+nu)+f(n+mu+nu))/4 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermSU33D(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmegaSq,
    DOUBLE* results
#else
    Real betaOverN, Real fOmegaSq,
    Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    const UINT plaqLength = __idx->m_pSmallData[CIndexData::kPlaqLengthIdx];
    const UINT plaqCountAll = __idx->m_pSmallData[CIndexData::kPlaqPerSiteIdx] * plaqLength;

#if !_CLG_DOUBLEFLOAT
    DOUBLE res = 0.0;
#else
    Real res = F(0.0);
#endif
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
        SIndex first = pCachedPlaqutte[idx * plaqLength + uiSiteIndex * plaqCountAll];
        deviceSU3 toAdd(_deviceGetGaugeBCSU3(pDeviceData, first));
        if (first.NeedToDagger())
        {
            toAdd.Dagger();
        }
        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedPlaqutte[idx * plaqLength + j + uiSiteIndex * plaqCountAll];
            deviceSU3 toMul(_deviceGetGaugeBCSU3(pDeviceData, first));
            if (first.NeedToDagger())
            {
                toAdd.MulDagger(toMul);
            }
            else
            {
                toAdd.Mul(toMul);
            }
        }

        //0 -> xz, 1 -> yz
        const BYTE mushift = idx0;
        const BYTE nushift = 2;
#if !_CLG_DOUBLEFLOAT
        res += static_cast<DOUBLE>(betaOverN * fOmegaSq * (3.0 - toAdd.ReTr()) * _deviceFi(byFieldId, sSite4, uiN, 2 - idx0, mushift, nushift));
#else
        res += betaOverN * fOmegaSq * (F(3.0) - toAdd.ReTr()) * _deviceFi(byFieldId, sSite4, uiN, 2 - idx0, mushift, nushift);
#endif
    }

    results[uiSiteIndex] = res;
}


/**
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermSU3_XY3D(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN,
    DOUBLE fOmegaSq
#else
    Real betaOverN,
    Real fOmegaSq
#endif
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(-0.5);

    //idx is the coefficient, 0 -> r^2, 1 -> y^2, 2 -> x^2
    BYTE idx[4] = {2, 1, 2, 1};
    //byOtherDir is the other direction of the staple
    BYTE byOtherDir[4] = {2, 2, 0, 1};
    //deviceSU3 plaqSum = deviceSU3::makeSU3Zero();
    #pragma unroll
    for (BYTE idir = 0; idir < 3; ++idir)
    {

        if (__idx->_deviceIsBondOnSurface(uiBigIdx, idir))
        {
            continue;
        }
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        deviceSU3 stap(_deviceStapleTermGfactor(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
            idir, 
            byOtherDir[idir],
            idx[idir]));
        if (2 == idir)
        {
            stap.Add(_deviceStapleTermGfactor(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
                idir,
                byOtherDir[idir + 1],
                idx[idir + 1]));
        }

        deviceSU3 force(pDeviceData[linkIndex]);

        force.MulDagger(stap);
        force.Ta();
        force.MulReal(betaOverN);
        pForceData[linkIndex].Add(force);
    }
}

#pragma endregion

#pragma region Chair Energy

/**
* Split into 3 functions to avoid max-register problem
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term1234_3D(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmega,
    DOUBLE* results
#else
    Real betaOverN, Real fOmega,
    Real* results
#endif
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
    const Real fXOmega = (sSite4.x - _DC_Centerx) * fOmega;
    const Real fYOmega = (sSite4.y - _DC_Centery) * fOmega;

    //===============
    //+x Omega V312
    const Real fV312 = -fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 2, 0, 1, uiN);

    //===============
    //+y Omega V321
    const Real fV321 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 2, 1, 0, uiN);

    results[uiSiteIndex] = (fV312  + fV321) * betaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term53D(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmegaSq,
    DOUBLE* results
#else
    Real betaOverN, Real fOmegaSq,
    Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

    if (__idx->m_pDeviceIndexPositionToSIndex[1][uiN].IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        return;
    }

    betaOverN = F(0.125) * betaOverN;
    const Real fXYOmega2 = (sSite4.x - _DC_Centerx) * (sSite4.y - _DC_Centery) * fOmegaSq;

    //===============
    //+Omega^2 xy V132
    const Real fV132 = fXYOmega2 * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);

    results[uiSiteIndex] = fV132 * betaOverN;
}


#pragma endregion

#pragma region Chair force

/**
* Split to 15 functions to avoid max-regcount
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term13D(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
    deviceSU3 *pForceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmega
#else
    Real betaOverN, Real fOmega
#endif
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * fOmega * F(0.0625);

    //===============
    //+x Omega V312
    //add force for dir=3
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 2);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 2))
    {
        const deviceSU3 staple_term1_4 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            2, 0, 1, 0);
        deviceSU3 force4(pDeviceData[uiLink4]);
        force4.MulDagger(staple_term1_4);
        force4.Ta();
        force4.MulReal(betaOverN);
        pForceData[uiLink4].Add(force4);
    }

    //===============
    //+x Omega V312
    //add force for dir=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    {
        const deviceSU3 staple_term1_2 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            1, 0, 2, 0);
        deviceSU3 force2(pDeviceData[uiLink2]);
        force2.MulDagger(staple_term1_2);
        force2.Ta();
        force2.MulReal(betaOverN);
        pForceData[uiLink2].Add(force2);
    }

    //===============
    //+x Omega V312
    //add force for dir=x
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    {
        const deviceSU3 staple_term1_1 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            2, 0, 1, 0);
        deviceSU3 force1(pDeviceData[uiLink1]);
        force1.MulDagger(staple_term1_1);
        force1.Ta();
        force1.MulReal(betaOverN);
        pForceData[uiLink1].Add(force1);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term33D(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
    deviceSU3 *pForceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmega
#else
    Real betaOverN, Real fOmega
#endif
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = -betaOverN * fOmega * F(0.0625);

    //===============
    //+y Omega V421
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 2);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 2))
    {
        const deviceSU3 staple_term3_4 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            2, 1, 0, _deviceHi1);
        deviceSU3 force4(pDeviceData[uiLink4]);
        force4.MulDagger(staple_term3_4);
        force4.Ta();
        force4.MulReal(betaOverN);
        pForceData[uiLink4].Add(force4);
    }

    //===============
    // y Omega V421
    //add force for mu=4
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    {
        const deviceSU3 staple_term3_1 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            0, 1, 2, _deviceHi1);
        deviceSU3 force1(pDeviceData[uiLink1]);
        force1.MulDagger(staple_term3_1);
        force1.Ta();
        force1.MulReal(betaOverN);
        pForceData[uiLink1].Add(force1);
    }

    //===============
    // y Omega V421
    //add force for mu=4
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    {
        const deviceSU3 staple_term3_2 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            2, 1, 0, _deviceHi1);
        deviceSU3 force2(pDeviceData[uiLink2]);
        force2.MulDagger(staple_term3_2);
        force2.Ta();
        force2.MulReal(betaOverN);
        pForceData[uiLink2].Add(force2);
    }

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term53D(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
    deviceSU3 *pForceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmegaSq
#else
    Real betaOverN, Real fOmegaSq
#endif
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmegaSq * F(0.125);

    //===============
    //+Omega^2 xy V132
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    {
        const deviceSU3 staple_term5_1 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            0, 2, 1, _deviceHi2);
        deviceSU3 force1(pDeviceData[uiLink1]);
        force1.MulDagger(staple_term5_1);
        force1.Ta();
        force1.MulReal(betaOverN);
        pForceData[uiLink1].Add(force1);
    }

    //===============
    //+Omega^2 xy V132
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    {
        const deviceSU3 staple_term5_2 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            1, 2, 0, _deviceHi2);
        deviceSU3 force2(pDeviceData[uiLink2]);
        force2.MulDagger(staple_term5_2);
        force2.Ta();
        force2.MulReal(betaOverN);
        pForceData[uiLink2].Add(force2);
    }

    //===============
    //+Omega^2 xy V132
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 2))
    {
        const deviceSU3 staple_term5_3 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            0, 2, 1, _deviceHi2);
        deviceSU3 force3(pDeviceData[uiLink3]);
        force3.MulDagger(staple_term5_3);
        force3.Ta();
        force3.MulReal(betaOverN);
        pForceData[uiLink3].Add(force3);
    }

}

#pragma endregion

#pragma region Projective plane

#pragma region Clover

__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermSU3_Shifted3D(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmegaSq,
    DOUBLE* results
#else
    Real betaOverN, Real fOmegaSq,
    Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

#if !_CLG_DOUBLEFLOAT
    DOUBLE fXSq = (sSite4.x - _DC_Centerx + 0.5);
    fXSq = fXSq * fXSq;
    DOUBLE fYSq = (sSite4.y - _DC_Centery + 0.5);
    fYSq = fYSq * fYSq;

    //======================================================
    //4-plaqutte terms
    //Omega^2 x^2 Retr[1 - U_2,3]
    const DOUBLE fU13 = fXSq * _device4PlaqutteTerm(pDeviceData, 0, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 y^2 Retr[1 - U_1,3]
    const DOUBLE fU23 = fYSq * _device4PlaqutteTerm(pDeviceData, 1, 2, uiBigIdx, sSite4, byFieldId);

#else
    Real fXSq = (sSite4.x - _DC_Centerx + F(0.5));
    fXSq = fXSq * fXSq;
    Real fYSq = (sSite4.y - _DC_Centery + F(0.5));
    fYSq = fYSq * fYSq;

    //======================================================
    //4-plaqutte terms
    //Omega^2 x^2 Retr[1 - U_2,3]
    const Real fU13 = fXSq * _device4PlaqutteTerm(pDeviceData, 0, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 y^2 Retr[1 - U_1,3]
    const Real fU23 = fYSq * _device4PlaqutteTerm(pDeviceData, 1, 2, uiBigIdx, sSite4, byFieldId);

#endif

    results[uiSiteIndex] = (fU23 + fU13) * betaOverN * fOmegaSq;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermSU3_XYZ_Shifted3D(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmegaSq
#else
    Real betaOverN, Real fOmegaSq
#endif
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(-0.5);
    //deviceSU3 plaqSum = deviceSU3::makeSU3Zero();
    BYTE idx[4] = { 2, 1, 2, 1 };
    BYTE byOtherDir[4] = { 2, 2, 0, 1 };

    #pragma unroll
    for (UINT idir = 0; idir < 3; ++idir)
    {
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        deviceSU3 stap(_deviceStapleTermGfactor(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
            idir,
            byOtherDir[idir],
            idx[idir],
            TRUE));
        if (2 == idir)
        {
            stap.Add(_deviceStapleTermGfactor(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
                idir,
                byOtherDir[idir + 1],
                idx[idir + 1],
                TRUE));
        }
        
        deviceSU3 force(pDeviceData[linkIndex]);
        force.MulDagger(stap);
        force.Ta();
        force.MulReal(betaOverN);
        pForceData[linkIndex].Add(force);
    }
}

#pragma endregion

#pragma region Chair energy

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term1234_Shifted3D(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmega,
    DOUBLE* results
#else
    Real betaOverN, Real fOmega,
    Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

#if !_CLG_DOUBLEFLOAT
    betaOverN = 0.125 * betaOverN;
    const DOUBLE fXOmega = (sSite4.x - _DC_Centerx + 0.5) * fOmega;
    const DOUBLE fYOmega = (sSite4.y - _DC_Centery + 0.5) * fOmega;

    //===============
    // F01F12 term  x Omega V312
    const DOUBLE fV312 = -fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 2, 0, 1, uiN);
    // F02F12 term  y Omega V321
    const DOUBLE fV321 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 2, 1, 0, uiN);

#else

    betaOverN = F(0.125) * betaOverN;
    const Real fXOmega = (sSite4.x - _DC_Centerx + F(0.5)) * fOmega;
    const Real fYOmega = (sSite4.y - _DC_Centery + F(0.5)) * fOmega;

    //===============
    // F01F12 term  x Omega V312
    const Real fV312 = -fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 2, 0, 1, uiN);
    // F02F12 term  y Omega V321
    const Real fV321 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 2, 1, 0, uiN);

#endif

    results[uiSiteIndex] = (fV312 + fV321) * betaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term5_Shifted3D(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmegaSq,
    DOUBLE* results
#else
    Real betaOverN, Real fOmegaSq,
    Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

#if !_CLG_DOUBLEFLOAT
    betaOverN = 0.125 * betaOverN;
    const DOUBLE fXYOmega2 = (sSite4.x - _DC_Centerx + 0.5) * (sSite4.y - _DC_Centery + 0.5) * fOmegaSq;

    //===============
    //+Omega^2 xy V142
    const DOUBLE fV132 = fXYOmega2 * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);
#else
    betaOverN = F(0.125) * betaOverN;
    const Real fXYOmega2 = (sSite4.x - _DC_Centerx + F(0.5)) * (sSite4.y - _DC_Centery + F(0.5)) * fOmegaSq;

    //===============
    //+Omega^2 xy V142
    const Real fV132 = fXYOmega2 * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);
#endif

    results[uiSiteIndex] = fV132 * betaOverN;
}

#pragma endregion

#pragma region Chair force

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term1_Shifted3D(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmega
#else
    Real betaOverN, Real fOmega
#endif
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
    const deviceSU3 staple_term1_4 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 0, 1, _deviceHiShifted0);
    deviceSU3 force4(pDeviceData[uiLink4]);
    force4.MulDagger(staple_term1_4);
    force4.Ta();
    force4.MulReal(betaOverN);
    pForceData[uiLink4].Add(force4);
    //}


    //===============
    //+x Omega V312
    //add force for dir=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    //{
    const deviceSU3 staple_term1_2 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        1, 0, 2, _deviceHiShifted0);
    deviceSU3 force2(pDeviceData[uiLink2]);
    force2.MulDagger(staple_term1_2);
    force2.Ta();
    force2.MulReal(betaOverN);
    pForceData[uiLink2].Add(force2);
   // }

    //===============
    //+x Omega V312
    //add force for dir=x
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    //{
    const deviceSU3 staple_term1_1 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 0, 1, _deviceHiShifted0);
    deviceSU3 force1(pDeviceData[uiLink1]);
    force1.MulDagger(staple_term1_1);
    force1.Ta();
    force1.MulReal(betaOverN);
    pForceData[uiLink1].Add(force1);
    //}
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term3_Shifted3D(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmega
#else
    Real betaOverN, Real fOmega
#endif
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
    const deviceSU3 staple_term3_4 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 1, 0, _deviceHiShifted1);
    deviceSU3 force4(pDeviceData[uiLink4]);
    force4.MulDagger(staple_term3_4);
    force4.Ta();
    force4.MulReal(betaOverN);
    pForceData[uiLink4].Add(force4);
    //}

    //===============
    //+ y Omega V321
    //add force for mu=1
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    //{
    const deviceSU3 staple_term3_1 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 1, 2, _deviceHiShifted1);
    deviceSU3 force1(pDeviceData[uiLink1]);
    force1.MulDagger(staple_term3_1);
    force1.Ta();
    force1.MulReal(betaOverN);
    pForceData[uiLink1].Add(force1);
    //}


    //===============
    //+ y Omega V321
    //add force for mu=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    //{
    const deviceSU3 staple_term3_2 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 1, 0, _deviceHiShifted1);
    deviceSU3 force2(pDeviceData[uiLink2]);
    force2.MulDagger(staple_term3_2);
    force2.Ta();
    force2.MulReal(betaOverN);
    pForceData[uiLink2].Add(force2);
    //}

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term5_Shifted3D(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE fOmegaSq
#else
    Real betaOverN, Real fOmegaSq
#endif
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
    const deviceSU3 staple_term5_1 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 2, 1, _deviceHiShifted2);
    deviceSU3 force1(pDeviceData[uiLink1]);
    force1.MulDagger(staple_term5_1);
    force1.Ta();
    force1.MulReal(betaOverN);
    pForceData[uiLink1].Add(force1);
    //}

    //===============
    //+ Omega^2 xy V132
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    //{
    const deviceSU3 staple_term5_2 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        1, 2, 0, _deviceHiShifted2);
    deviceSU3 force2(pDeviceData[uiLink2]);
    force2.MulDagger(staple_term5_2);
    force2.Ta();
    force2.MulReal(betaOverN);
    pForceData[uiLink2].Add(force2);
    //}

    //===============
    //+ Omega^2 xy V132
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    //{
    const deviceSU3 staple_term5_3 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 2, 1, _deviceHiShifted2);
    deviceSU3 force3(pDeviceData[uiLink3]);
    force3.MulDagger(staple_term5_3);
    force3.Ta();
    force3.MulReal(betaOverN);
    pForceData[uiLink3].Add(force3);
    //}
}

#pragma endregion

#pragma endregion

#pragma endregion

CActionGaugePlaquetteRotating3D::CActionGaugePlaquetteRotating3D()
    : CActionGaugePlaquetteRotating()
{
}

UBOOL CActionGaugePlaquetteRotating3D::CalculateForceOnGauge(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
#if !_CLG_DOUBLEFLOAT
    pGauge->CalculateForceAndStaple(pForce, pStaple, static_cast<Real>(m_fBetaOverN));
#else
    pGauge->CalculateForceAndStaple(pForce, pStaple, m_fBetaOverN);
#endif

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    if (NULL == pGaugeSU3 || NULL == pForceSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3 now.\n"));
        return TRUE;
    }

    preparethread;


    if (!m_bShiftHalfCoord)
    {
        _kernelAddForce4PlaqutteTermSU3_XY3D << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

        _kernelAddForceChairTermSU3_Term13D << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term33D << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term53D << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);
    }
    else
    {

        _kernelAddForce4PlaqutteTermSU3_XYZ_Shifted3D << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);
        
        _kernelAddForceChairTermSU3_Term1_Shifted3D << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term3_Shifted3D << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term5_Shifted3D << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);
    }

    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
#if !_CLG_DOUBLEFLOAT
DOUBLE CActionGaugePlaquetteRotating3D::Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
#else
Real CActionGaugePlaquetteRotating3D::Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
#endif
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

    if (m_bCloverEnergy)
    {
        m_fNewEnergy = pGauge->CalculatePlaqutteEnergyUseClover(m_fBetaOverN);
    }
    else
    {
        m_fNewEnergy = pGauge->CalculatePlaqutteEnergy(m_fBetaOverN);
    }
    
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    if (NULL == pGaugeSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3 now.\n"));
        return m_fNewEnergy;
    }

    preparethread;

    appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

    if (m_bShiftHalfCoord)
    {

        _kernelAdd4PlaqutteTermSU3_Shifted3D << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _kernelAddChairTermSU3_Term1234_Shifted3D << <block, threads >> > (
            pGaugeSU3->m_byFieldId, 
            pGaugeSU3->m_pDeviceData, 
            m_fBetaOverN, 
            m_fOmega, 
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _kernelAddChairTermSU3_Term5_Shifted3D << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    }
    else
    {
        _kernelAdd4PlaqutteTermSU33D << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _kernelAddChairTermSU3_Term1234_3D << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            m_fBetaOverN,
            m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _kernelAddChairTermSU3_Term53D << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    }


    return m_fNewEnergy;
}

CCString CActionGaugePlaquetteRotating3D::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CActionGaugePlaquetteRotating3D\n");
    sRet = sRet + tab + _T("Beta : ") + appFloatToString(CCommonData::m_fBeta) + _T("\n");
    sRet = sRet + tab + _T("Omega : ") + appFloatToString(m_fOmega) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================