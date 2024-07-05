//=============================================================================
// FILENAME : CActionGaugePlaquetteRotating.cu
// 
// DESCRIPTION:
// This is the class for rotating su3
//
// REVISION:
//  [05/07/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionGaugePlaquetteRotating.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotating)



#pragma region kernels

#pragma region Clover terms

/**
* This is slower, just for testing
* directly calculate Retr[1 - \hat{U}]
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermSU3_Test(
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

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx].IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        return;
    }

    Real fXSq = (sSite4.x - _DC_Centerx);
    fXSq = fXSq * fXSq;
    Real fYSq = (sSite4.y - _DC_Centery);
    fYSq = fYSq * fYSq;

    //======================================================
    //4-plaqutte terms
    //Omega^2 x^2 Retr[1 - U_2,3]
    const Real fU23 = fOmegaSq * fXSq * _device4PlaqutteTerm(pDeviceData, 1, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 y^2 Retr[1 - U_1,3]
    const Real fU13 = fOmegaSq * fYSq * _device4PlaqutteTerm(pDeviceData, 0, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 (x^2 + y^2) Retr[1 - U_1,2]
    const Real fU12 = fOmegaSq * (fXSq + fYSq) * _device4PlaqutteTerm(pDeviceData, 0, 1, uiBigIdx, sSite4, byFieldId);

    results[uiSiteIndex] = (fU23 + fU13 + fU12) * betaOverN;
}

/**
* Using plaqutte and (f(n)+f(n+mu)+f(n+nu)+f(n+mu+nu))/4 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermSU3(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
    DOUBLE betaOverN, DOUBLE fOmegaSq,
    DOUBLE* results
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
        SIndex first = pCachedPlaqutte[idx * plaqLength + uiSiteIndex * plaqCountAll];
        deviceSU3 toAdd(_deviceGetGaugeBCSU3(byFieldId, pDeviceData, first));
        if (first.NeedToDagger())
        {
            toAdd.Dagger();
        }
        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedPlaqutte[idx * plaqLength + j + uiSiteIndex * plaqCountAll];
            deviceSU3 toMul(_deviceGetGaugeBCSU3(byFieldId, pDeviceData, first));
            if (first.NeedToDagger())
            {
                toAdd.MulDagger(toMul);
            }
            else
            {
                toAdd.Mul(toMul);
            }
        }

        //0 -> xy, 1 -> xz, 2 -> yz
        //x x y
        const BYTE mushift = (idx0 / 2);
        //y z z
        const BYTE nushift = ((idx0 + 1) / 2) + 1;
        res += static_cast<DOUBLE>(betaOverN * fOmegaSq * (3.0 - toAdd.ReTr()) * _deviceFi(byFieldId, sSite4, uiN, idx0, mushift, nushift));
    }

    results[uiSiteIndex] = res;
}


/**
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermSU3_XY(
    BYTE byFieldId,
    UBOOL bTorus,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
    DOUBLE betaOverN,
    DOUBLE fOmegaSq
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(-0.5);
    BYTE idx[6] = { 1, 0, 2, 0, 1, 2};
    BYTE byOtherDir[6] = {2, 1, 2, 0, 0, 1};
    //deviceSU3 plaqSum = deviceSU3::makeSU3Zero();
    #pragma unroll
    for (BYTE idir = 0; idir < 3; ++idir)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, idir))
        {
            continue;
        }
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        //for xz, yz, i=1,2
        //for xy, i = 0

        //idir = 0: xz, xy with coefficient 1 0
        //idir = 1: yz, yx with coefficient 2 0
        //idir = 2: zx, zy with coefficient 1 2

        deviceSU3 stap(_deviceStapleTermGfactor(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
            idir, 
            byOtherDir[2 * idir],
            idx[2 * idir]));
        stap.Add(_deviceStapleTermGfactor(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
            idir,
            byOtherDir[2 * idir + 1],
            idx[2 * idir + 1]));
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
_kernelAddChairTermSU3_Term12(
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
    const Real fXOmega = -(sSite4.x - _DC_Centerx) * fOmega;

    //===============
    //-x Omega V412
    const Real fV412 = fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 0, 1, uiN);

    //===============
    //-x Omega V432
    const Real fV432 = fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 1, uiN);

    results[uiSiteIndex] = (fV412 + fV432) * betaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term34(
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

    if (__idx->m_pDeviceIndexPositionToSIndex[1][uiN].IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        return;
    }

    betaOverN = F(0.125) * betaOverN;
    const Real fYOmega = (sSite4.y - _DC_Centery) * fOmega;

    //===============
    //+y Omega V421
    const Real fV421 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 1, 0, uiN);

    //===============
    //+y Omega V431
    const Real fV431 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 0, uiN);

    results[uiSiteIndex] = (fV421 + fV431) * betaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term5(
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
    const Real fXYOmega2 = -(sSite4.x - _DC_Centerx) * (sSite4.y - _DC_Centery) * fOmegaSq;

    //===============
    //-Omega^2 xy V132
    const Real fV132 = fXYOmega2 * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);

    results[uiSiteIndex] = fV132 * betaOverN;
}


#pragma endregion

#pragma region Chair force

/**
* Split to 15 functions to avoid max-regcount
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term1(
    BYTE byFieldId,
    UBOOL bTorus,
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

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //-x Omega V412
    //add force for dir=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        const deviceSU3 staple_term1_4 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 0, 1, bTorus ? _deviceHi0T : _deviceHi0);
        deviceSU3 force4(pDeviceData[uiLink4]);
        force4.MulDagger(staple_term1_4);
        force4.Ta();
        force4.MulReal(betaOverN);
        pForceData[uiLink4].Add(force4);
    }

    //===============
    //+x Omega V412
    //add force for dir=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    {
        const deviceSU3 staple_term1_2 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            1, 0, 3, bTorus ? _deviceHi0T : _deviceHi0);
        deviceSU3 force2(pDeviceData[uiLink2]);
        force2.MulDagger(staple_term1_2);
        force2.Ta();
        force2.MulReal(betaOverN);
        pForceData[uiLink2].Add(force2);
    }

    //===============
    //+x Omega V412
    //add force for dir=x
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    {
        const deviceSU3 staple_term1_1 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 0, 1, bTorus ? _deviceHi0T : _deviceHi0);
        deviceSU3 force1(pDeviceData[uiLink1]);
        force1.MulDagger(staple_term1_1);
        force1.Ta();
        force1.MulReal(betaOverN);
        pForceData[uiLink1].Add(force1);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term2(
    BYTE byFieldId,
    UBOOL bTorus,
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

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //-x Omega V432
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        const deviceSU3 staple_term2_4 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 2, 1, bTorus ? _deviceHi0T : _deviceHi0);
        deviceSU3 force4(pDeviceData[uiLink4]);
        force4.MulDagger(staple_term2_4);
        force4.Ta();
        force4.MulReal(betaOverN);
        pForceData[uiLink4].Add(force4);
    }

    //===============
    //+x Omega V432
    //add force for mu=4
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    {
        const deviceSU3 staple_term2_2 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            1, 2, 3, bTorus ? _deviceHi0T : _deviceHi0);
        deviceSU3 force2(pDeviceData[uiLink2]);
        force2.MulDagger(staple_term2_2);
        force2.Ta();
        force2.MulReal(betaOverN);
        pForceData[uiLink2].Add(force2);
    }

    //===============
    //+x Omega V432
    //add force for mu=4
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 2))
    {
        const deviceSU3 staple_term2_3 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 2, 1, bTorus ? _deviceHi0T : _deviceHi0);
        deviceSU3 force3(pDeviceData[uiLink3]);
        force3.MulDagger(staple_term2_3);
        force3.Ta();
        force3.MulReal(betaOverN);
        pForceData[uiLink3].Add(force3);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term3(
    BYTE byFieldId,
    UBOOL bTorus,
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

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+y Omega V421
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        const deviceSU3 staple_term3_4 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 1, 0, bTorus ? _deviceHi1T : _deviceHi1);
        deviceSU3 force4(pDeviceData[uiLink4]);
        force4.MulDagger(staple_term3_4);
        force4.Ta();
        force4.MulReal(betaOverN);
        pForceData[uiLink4].Add(force4);
    }

    //===============
    //-y Omega V421
    //add force for mu=4
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    {
        const deviceSU3 staple_term3_1 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            0, 1, 3, bTorus ? _deviceHi1T : _deviceHi1);
        deviceSU3 force1(pDeviceData[uiLink1]);
        force1.MulDagger(staple_term3_1);
        force1.Ta();
        force1.MulReal(betaOverN);
        pForceData[uiLink1].Add(force1);
    }

    //===============
    //-y Omega V421
    //add force for mu=4
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    {
        const deviceSU3 staple_term3_2 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 1, 0, bTorus ? _deviceHi1T : _deviceHi1);
        deviceSU3 force2(pDeviceData[uiLink2]);
        force2.MulDagger(staple_term3_2);
        force2.Ta();
        force2.MulReal(betaOverN);
        pForceData[uiLink2].Add(force2);
    }

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term4(
    BYTE byFieldId,
    UBOOL bTorus,
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

    betaOverN = -betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+y Omega V431
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        const deviceSU3 staple_term4_4 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 2, 0, bTorus ? _deviceHi1T : _deviceHi1);
        deviceSU3 force4(pDeviceData[uiLink4]);
        force4.MulDagger(staple_term4_4);
        force4.Ta();
        force4.MulReal(betaOverN);
        pForceData[uiLink4].Add(force4);
    }

    //===============
    //-y Omega V431
    //add force for mu=4
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    {
        const deviceSU3 staple_term4_1 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            0, 2, 3, bTorus ? _deviceHi1T : _deviceHi1);
        deviceSU3 force1(pDeviceData[uiLink1]);
        force1.MulDagger(staple_term4_1);
        force1.Ta();
        force1.MulReal(betaOverN);
        pForceData[uiLink1].Add(force1);
    }

    //===============
    //-y Omega V431
    //add force for mu=4
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 2))
    {
        const deviceSU3 staple_term4_3 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            3, 2, 0, bTorus ? _deviceHi1T : _deviceHi1);
        deviceSU3 force3(pDeviceData[uiLink3]);
        force3.MulDagger(staple_term4_3);
        force3.Ta();
        force3.MulReal(betaOverN);
        pForceData[uiLink3].Add(force3);
    }

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term5(
    BYTE byFieldId,
    UBOOL bTorus,
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

    betaOverN = -betaOverN * F(0.5) * fOmegaSq * F(0.125);

    //===============
    //-Omega^2 xy V132
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    {
        const deviceSU3 staple_term5_1 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
            0, 2, 1, bTorus ? _deviceHi2T : _deviceHi2);
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
            1, 2, 0, bTorus ? _deviceHi2T : _deviceHi2);
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
            0, 2, 1, bTorus ? _deviceHi2T : _deviceHi2);
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
_kernelAdd4PlaqutteTermSU3_Shifted(
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
    const DOUBLE fU23 = fXSq * _device4PlaqutteTerm(pDeviceData, 1, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 y^2 Retr[1 - U_1,3]
    const DOUBLE fU13 = fYSq * _device4PlaqutteTerm(pDeviceData, 0, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 (x^2 + y^2) Retr[1 - U_1,2]
    const DOUBLE fU12 = (fXSq + fYSq) * _device4PlaqutteTerm(pDeviceData, 0, 1, uiBigIdx, sSite4, byFieldId);
#else
    Real fXSq = (sSite4.x - _DC_Centerx + F(0.5));
    fXSq = fXSq * fXSq;
    Real fYSq = (sSite4.y - _DC_Centery + F(0.5));
    fYSq = fYSq * fYSq;

    //======================================================
    //4-plaqutte terms
    //Omega^2 x^2 Retr[1 - U_2,3]
    const Real fU23 = fXSq * _device4PlaqutteTerm(pDeviceData, 1, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 y^2 Retr[1 - U_1,3]
    const Real fU13 = fYSq * _device4PlaqutteTerm(pDeviceData, 0, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 (x^2 + y^2) Retr[1 - U_1,2]
    const Real fU12 = (fXSq + fYSq) * _device4PlaqutteTerm(pDeviceData, 0, 1, uiBigIdx, sSite4, byFieldId);
#endif

    results[uiSiteIndex] = (fU23 + fU13 + fU12) * betaOverN * fOmegaSq;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermSU3_XYZ_Shifted(
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
    BYTE idx[6] = { 1, 0, 2, 0, 1, 2 };
    BYTE byOtherDir[6] = { 2, 1, 2, 0, 0, 1 };

    #pragma unroll
    for (UINT idir = 0; idir < 3; ++idir)
    {
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        //mu = idir, nu = 4, i = mu
        deviceSU3 stap(_deviceStapleTermGfactor(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
            idir,
            byOtherDir[2 * idir],
            idx[2 * idir],
            TRUE));
        stap.Add(_deviceStapleTermGfactor(byFieldId, bTorus, pDeviceData, sSite4, fOmegaSq, uiBigIdx,
            idir,
            byOtherDir[2 * idir + 1],
            idx[2 * idir + 1],
            TRUE));
        
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
_kernelAddChairTermSU3_Term12_Shifted(
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
    betaOverN = -0.125 * betaOverN;
    const DOUBLE fXOmega = (sSite4.x - _DC_Centerx + 0.5) * fOmega;

    //===============
    //- x Omega V412
    const DOUBLE fV412 = fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 0, 1, uiN);

    //===============
    //- x Omega V432
    const DOUBLE fV432 = fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 1, uiN);

#else
    betaOverN = -F(0.125) * betaOverN;
    const Real fXOmega = (sSite4.x - _DC_Centerx + F(0.5)) * fOmega;

    //===============
    //+x Omega V412
    const Real fV412 = fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 0, 1, uiN);

    //===============
    //+x Omega V432
    const Real fV432 = fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 1, uiN);
#endif

    results[uiSiteIndex] = (fV412 + fV432) * betaOverN;
}


__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term34_Shifted(
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
    const DOUBLE fYOmega = (sSite4.y - _DC_Centery + 0.5) * fOmega;

    //===============
    //+ y Omega V421
    const DOUBLE fV421 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 1, 0, uiN);

    //===============
    //+ y Omega V431
    const DOUBLE fV431 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 0, uiN);
#else
    betaOverN = F(0.125) * betaOverN;
    const Real fYOmega = (sSite4.y - _DC_Centery + F(0.5)) * fOmega;

    //===============
    //-y Omega V421
    const Real fV421 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 1, 0, uiN);

    //===============
    //-y Omega V431
    const Real fV431 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 0, uiN);
#endif

    results[uiSiteIndex] = (fV421 + fV431) * betaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term5_Shifted(
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
    betaOverN = -0.125 * betaOverN;
    const DOUBLE fXYOmega2 = (sSite4.x - _DC_Centerx + 0.5) * (sSite4.y - _DC_Centery + 0.5) * fOmegaSq;

    //===============
    //-Omega^2 xy V142
    const DOUBLE fV132 = fXYOmega2 * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);
#else
    betaOverN = -F(0.125) * betaOverN;
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
_kernelAddForceChairTermSU3_Term1_Shifted(
    BYTE byFieldId,
    UBOOL bTorus,
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
    //-x Omega V412
    //add force for dir=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    //{
    const deviceSU3 staple_term1_4 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 0, 1, bTorus ? _deviceHiShiftedT0 : _deviceHiShifted0);
    deviceSU3 force4(pDeviceData[uiLink4]);
    force4.MulDagger(staple_term1_4);
    force4.Ta();
    force4.MulReal(betaOverN);
    pForceData[uiLink4].Add(force4);
    //}


    //===============
    //+x Omega V412
    //add force for dir=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    //{
    const deviceSU3 staple_term1_2 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        1, 0, 3, bTorus ? _deviceHiShiftedT0 : _deviceHiShifted0);
    deviceSU3 force2(pDeviceData[uiLink2]);
    force2.MulDagger(staple_term1_2);
    force2.Ta();
    force2.MulReal(betaOverN);
    pForceData[uiLink2].Add(force2);
   // }

    //===============
    //+x Omega V412
    //add force for dir=x
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    //{
    const deviceSU3 staple_term1_1 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 0, 1, bTorus ? _deviceHiShiftedT0 : _deviceHiShifted0);
    deviceSU3 force1(pDeviceData[uiLink1]);
    force1.MulDagger(staple_term1_1);
    force1.Ta();
    force1.MulReal(betaOverN);
    pForceData[uiLink1].Add(force1);
    //}
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term2_Shifted(
    BYTE byFieldId,
    UBOOL bTorus,
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
    //-x Omega V432
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    //{
    const deviceSU3 staple_term2_4 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 2, 1, bTorus ? _deviceHiShiftedT0 : _deviceHiShifted0);
    deviceSU3 force4(pDeviceData[uiLink4]);
    force4.MulDagger(staple_term2_4);
    force4.Ta();
    force4.MulReal(betaOverN);
    pForceData[uiLink4].Add(force4);
    //}

    //===============
    //+x Omega V432
    //add force for mu=4
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    //{
    const deviceSU3 staple_term2_2 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        1, 2, 3, bTorus ? _deviceHiShiftedT0 : _deviceHiShifted0);
    deviceSU3 force2(pDeviceData[uiLink2]);
    force2.MulDagger(staple_term2_2);
    force2.Ta();
    force2.MulReal(betaOverN);
    pForceData[uiLink2].Add(force2);
    //}

    //===============
    //+x Omega V432
    //add force for mu=4
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 2))
    //{
    const deviceSU3 staple_term2_3 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 2, 1, bTorus ? _deviceHiShiftedT0 : _deviceHiShifted0);
    deviceSU3 force3(pDeviceData[uiLink3]);
    force3.MulDagger(staple_term2_3);
    force3.Ta();
    force3.MulReal(betaOverN);
    pForceData[uiLink3].Add(force3);
    //}
}


__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term3_Shifted(
    BYTE byFieldId,
    UBOOL bTorus,
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
    //+ y Omega V421
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    //{
    const deviceSU3 staple_term3_4 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 1, 0, bTorus ? _deviceHiShiftedT1 : _deviceHiShifted1);
    deviceSU3 force4(pDeviceData[uiLink4]);
    force4.MulDagger(staple_term3_4);
    force4.Ta();
    force4.MulReal(betaOverN);
    pForceData[uiLink4].Add(force4);
    //}

    //===============
    //+ y Omega V421
    //add force for mu=1
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    //{
    const deviceSU3 staple_term3_1 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 1, 3, bTorus ? _deviceHiShiftedT1 : _deviceHiShifted1);
    deviceSU3 force1(pDeviceData[uiLink1]);
    force1.MulDagger(staple_term3_1);
    force1.Ta();
    force1.MulReal(betaOverN);
    pForceData[uiLink1].Add(force1);
    //}


    //===============
    //+ y Omega V421
    //add force for mu=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    //{
    const deviceSU3 staple_term3_2 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 1, 0, bTorus ? _deviceHiShiftedT1 : _deviceHiShifted1);
    deviceSU3 force2(pDeviceData[uiLink2]);
    force2.MulDagger(staple_term3_2);
    force2.Ta();
    force2.MulReal(betaOverN);
    pForceData[uiLink2].Add(force2);
    //}

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term4_Shifted(
    BYTE byFieldId,
    UBOOL bTorus,
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
    //+ y Omega V431
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    //{
    const deviceSU3 staple_term4_4 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 2, 0, bTorus ? _deviceHiShiftedT1 : _deviceHiShifted1);
    deviceSU3 force4(pDeviceData[uiLink4]);
    force4.MulDagger(staple_term4_4);
    force4.Ta();
    force4.MulReal(betaOverN);
    pForceData[uiLink4].Add(force4);
    //}


    //===============
    //+ y Omega V431
    //add force for mu=4
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    //{
    const deviceSU3 staple_term4_1 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 2, 3, bTorus ? _deviceHiShiftedT1 : _deviceHiShifted1);
    deviceSU3 force1(pDeviceData[uiLink1]);
    force1.MulDagger(staple_term4_1);
    force1.Ta();
    force1.MulReal(betaOverN);
    pForceData[uiLink1].Add(force1);
    //}

    //===============
    //+ y Omega V431
    //add force for mu=3
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 2))
    //{
    const deviceSU3 staple_term4_3 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 2, 0, bTorus ? _deviceHiShiftedT1 : _deviceHiShifted1);
    deviceSU3 force3(pDeviceData[uiLink3]);
    force3.MulDagger(staple_term4_3);
    force3.Ta();
    force3.MulReal(betaOverN);
    pForceData[uiLink3].Add(force3);
    //}

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term5_Shifted(
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

    betaOverN = -betaOverN * F(0.5) * fOmegaSq * F(0.125);

    //===============
    //- Omega^2 xy V132
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    //{
    const deviceSU3 staple_term5_1 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 2, 1, bTorus ? _deviceHiShiftedT2 : _deviceHiShifted2);
    deviceSU3 force1(pDeviceData[uiLink1]);
    force1.MulDagger(staple_term5_1);
    force1.Ta();
    force1.MulReal(betaOverN);
    pForceData[uiLink1].Add(force1);
    //}

    //===============
    //- Omega^2 xy V132
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    //{
    const deviceSU3 staple_term5_2 = _deviceStapleChairTerm1(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        1, 2, 0, bTorus ? _deviceHiShiftedT2 : _deviceHiShifted2);
    deviceSU3 force2(pDeviceData[uiLink2]);
    force2.MulDagger(staple_term5_2);
    force2.Ta();
    force2.MulReal(betaOverN);
    pForceData[uiLink2].Add(force2);
    //}

    //===============
    //- Omega^2 xy V132
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    //if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    //{
    const deviceSU3 staple_term5_3 = _deviceStapleChairTerm2(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        0, 2, 1, bTorus ? _deviceHiShiftedT2 : _deviceHiShifted2);
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

#if 0

#pragma region Detailed about chair terms

static __device__ __inline__ Real _deviceOneChairLoop(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 start,
    INT mu, INT nu, INT rho
    )
{
    const INT path[6] = { mu, nu, -mu, rho, -nu, -rho};
    return _deviceLink(pDeviceData, start, 6, byFieldId, path).ReTr();
}

static __device__ __inline__ deviceSU3 _deviceChairStapleMu_412(
    BYTE byFieldId, BYTE byCoeffType,
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 start, const SSmallInt4& sCenter,
    INT mu, INT nu, INT rho)
{
    const SSmallInt4 siteA = start;
    const SSmallInt4 siteB = _deviceSmallInt4OffsetC(start, -nu);

    const INT path1[5] = { rho, nu, -rho, mu, -nu };
    deviceSU3 linkV = _deviceLink(pDeviceData, start, 5, byFieldId, path1);
    linkV.MulReal(_deviceSiteCoeff(siteA, sCenter, byFieldId, byCoeffType));

    const INT path2[5] = { rho, -nu, -rho, mu, nu };
    linkV.Add(_deviceLink(pDeviceData, start, 5, byFieldId, path2)
        .MulRealC(_deviceSiteCoeff(siteB, sCenter, byFieldId, byCoeffType)));
    return linkV;
}

static __device__ __inline__ deviceSU3 _deviceChairStapleMu2_412(
    BYTE byFieldId, BYTE byCoeffType,
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 start, const SSmallInt4& sCenter,
    INT mu, INT nu, INT rho)
{
    const SSmallInt4 siteA = _deviceSmallInt4OffsetC(start, -mu);
    const SSmallInt4 siteB = _deviceSmallInt4OffsetC(siteA, -nu);
    mu = -mu;

    const INT path1[5] = { nu, mu, rho, -nu, -rho };
    deviceSU3 linkV = _deviceLink(pDeviceData, start, 5, byFieldId, path1);
    linkV.MulReal(_deviceSiteCoeff(siteA, sCenter, byFieldId, byCoeffType));

    const INT path2[5] = { -nu, mu, rho, nu, -rho };
    linkV.Add(_deviceLink(pDeviceData, start, 5, byFieldId, path2)
        .MulRealC(_deviceSiteCoeff(siteB, sCenter, byFieldId, byCoeffType)));
    return linkV;
}

static __device__ __inline__ deviceSU3 _deviceChairStapleNu_412(
    BYTE byFieldId, BYTE byCoeffType,
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 start, const SSmallInt4& sCenter,
    INT mu, INT nu, INT rho)
{
    SSmallInt4 sStart2 = start;
    if (nu < 0)
    {
        _deviceSmallInt4Offset(sStart2, -nu);
        nu = -nu;
    }
    SSmallInt4 siteA = _deviceSmallInt4OffsetC(sStart2, -mu);
    SSmallInt4 siteB = _deviceSmallInt4OffsetC(sStart2, -rho);
    const INT path1[5] = { -mu, rho, nu, -rho, mu };
    deviceSU3 linkV = _deviceLink(pDeviceData, start, 5, byFieldId, path1);
    linkV.MulReal(_deviceSiteCoeff(siteA, sCenter, byFieldId, byCoeffType));

    const INT path2[5] = { -rho, mu, nu, -mu, rho };
    linkV.Add(_deviceLink(pDeviceData, start, 5, byFieldId, path2)
        .MulRealC(_deviceSiteCoeff(siteB, sCenter, byFieldId, byCoeffType)));

    return linkV;
}

static __device__ __inline__ deviceSU3 _deviceChairStapleRho_412(
    BYTE byFieldId, BYTE byCoeffType,
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 start, const SSmallInt4& sCenter,
    INT mu, INT nu, INT rho)
{
    const SSmallInt4 siteA = start;
    const SSmallInt4 siteB = _deviceSmallInt4OffsetC(start, -nu);

    const INT path1[5] = { mu, nu, -mu, rho, -nu }; 
    deviceSU3 linkV = _deviceLink(pDeviceData, start, 5, byFieldId, path1);
    linkV.MulReal(_deviceSiteCoeff(siteA, sCenter, byFieldId, byCoeffType));

    const INT path2[5] = { mu, -nu, -mu, rho, nu };
    linkV.Add(_deviceLink(pDeviceData, start, 5, byFieldId, path2)
        .MulRealC(_deviceSiteCoeff(siteB, sCenter, byFieldId, byCoeffType)));
    return linkV;
}

static __device__ __inline__ deviceSU3 _deviceChairStapleRho2_412(
    BYTE byFieldId, BYTE byCoeffType,
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 start, const SSmallInt4& sCenter,
    INT mu, INT nu, INT rho)
{
    //rho < 0,
    const SSmallInt4 siteA = _deviceSmallInt4OffsetC(start, -rho);
    const SSmallInt4 siteB = _deviceSmallInt4OffsetC(siteA, -nu);

    rho = -rho;
    const INT path1[5] = { nu, rho, mu, -nu, -mu };
    deviceSU3 linkV = _deviceLink(pDeviceData, start, 5, byFieldId, path1);
    linkV.MulReal(_deviceSiteCoeff(siteA, sCenter, byFieldId, byCoeffType));

    const INT path2[5] = { -nu, rho, mu, nu, -mu };
    linkV.Add(_deviceLink(pDeviceData, start, 5, byFieldId, path2)
        .MulRealC(_deviceSiteCoeff(siteB, sCenter, byFieldId, byCoeffType)));
    return linkV;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddLoopsForEachSite_412(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    Real betaOverN, Real fOmega,
    Real* results)
{
    intokernalInt4;

    const Real fXOmega = _deviceSiteCoeff(sSite4, sCenterSite, byFieldId, 0) * fOmega;
    Real fLoop1 = //F(0.0);
    _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, 4, 1, 2);
    fLoop1 = fLoop1 + _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, 4, -1, 2);
    fLoop1 = fLoop1 - _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, 4, -1, -2);
    fLoop1 = fLoop1 - _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, 4, 1, -2);

    fLoop1 = fLoop1 + _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, -4, 1, -2);
    fLoop1 = fLoop1 + _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, -4, -1, -2);
    fLoop1 = fLoop1 - _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, -4, -1, 2);
    fLoop1 = fLoop1 - _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, -4, 1, 2);

    results[uiSiteIndex] = fLoop1 * betaOverN * fXOmega;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddLoopsForceForEachSite_412(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3* pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    betaOverN = betaOverN * F(0.5);

    //Mu
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);
    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        deviceSU3 force4 = // deviceSU3::makeSU3Zero(); 
        _deviceChairStapleMu_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, 4, 1, 2);
        force4.Add(_deviceChairStapleMu_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, 4, -1, 2));
        force4.Sub(_deviceChairStapleMu_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, 4, -1, -2));
        force4.Sub(_deviceChairStapleMu_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, 4, 1, -2));

        force4.Add(_deviceChairStapleMu2_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, -4, 1, -2));
        force4.Add(_deviceChairStapleMu2_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, -4, -1, -2));
        force4.Sub(_deviceChairStapleMu2_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, -4, -1, 2));
        force4.Sub(_deviceChairStapleMu2_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, -4, 1, 2));

        force4.MulDagger(pDeviceData[uiLink4]);
        force4.Ta();
        force4.MulReal(betaOverN * fOmega);
        pForceData[uiLink4].Sub(force4);
    }

    //Nu
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);
    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    {
        deviceSU3 force1 = // deviceSU3::makeSU3Zero();
        _deviceChairStapleNu_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, 4, 1, 2);
        force1.Add(_deviceChairStapleNu_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, 4, -1, 2));
        force1.Sub(_deviceChairStapleNu_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, 4, -1, -2));
        force1.Sub(_deviceChairStapleNu_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, 4, 1, -2));

        force1.Add(_deviceChairStapleNu_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, -4, 1, -2));
        force1.Add(_deviceChairStapleNu_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, -4, -1, -2));
        force1.Sub(_deviceChairStapleNu_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, -4, -1, 2));
        force1.Sub(_deviceChairStapleNu_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, -4, 1, 2));

        force1.MulDagger(pDeviceData[uiLink1]);
        force1.Ta();
        force1.MulReal(betaOverN * fOmega);
        pForceData[uiLink1].Sub(force1);
    }

    //Rho
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);
    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    {
        deviceSU3 force2 = // deviceSU3::makeSU3Zero();
        _deviceChairStapleRho_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, 4, 1, 2);
        force2.Add(_deviceChairStapleRho_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, 4, -1, 2));
        force2.Sub(_deviceChairStapleRho2_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, 4, -1, -2));
        force2.Sub(_deviceChairStapleRho2_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, 4, 1, -2));

        force2.Add(_deviceChairStapleRho2_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, -4, 1, -2));
        force2.Add(_deviceChairStapleRho2_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, -4, -1, -2));
        force2.Sub(_deviceChairStapleRho_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, -4, -1, 2));
        force2.Sub(_deviceChairStapleRho_412(byFieldId, 0, pDeviceData, sSite4, sCenterSite, -4, 1, 2));

        force2.MulDagger(pDeviceData[uiLink2]);
        force2.Ta();
        force2.MulReal(betaOverN * fOmega);
        pForceData[uiLink2].Sub(force2);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddLoopsForEachSite_421(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    Real betaOverN, Real fOmega,
    Real* results)
{
    intokernalInt4;

    const Real fXOmega = _deviceSiteCoeff(sSite4, sCenterSite, byFieldId, 1) * fOmega;
    Real fLoop1 = //F(0.0);
        _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, 4, 2, 1);
    fLoop1 = fLoop1 + _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, 4, -2, 1);
    fLoop1 = fLoop1 - _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, 4, -2, -1);
    fLoop1 = fLoop1 - _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, 4, 2, -1);

    fLoop1 = fLoop1 + _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, -4, 2, -1);
    fLoop1 = fLoop1 + _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, -4, -2, -1);
    fLoop1 = fLoop1 - _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, -4, -2, 1);
    fLoop1 = fLoop1 - _deviceOneChairLoop(byFieldId, pDeviceData, sSite4, -4, 2, 1);

    results[uiSiteIndex] = fLoop1 * betaOverN * fXOmega;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddLoopsForceForEachSite_421(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3* pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    betaOverN = betaOverN * F(0.5);

    //Mu
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);
    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        deviceSU3 force4 = // deviceSU3::makeSU3Zero(); 
            _deviceChairStapleMu_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, 4, 2, 1);
        force4.Add(_deviceChairStapleMu_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, 4, -2, 1));
        force4.Sub(_deviceChairStapleMu_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, 4, -2, -1));
        force4.Sub(_deviceChairStapleMu_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, 4, 2, -1));

        force4.Add(_deviceChairStapleMu2_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, -4, 2, -1));
        force4.Add(_deviceChairStapleMu2_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, -4, -2, -1));
        force4.Sub(_deviceChairStapleMu2_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, -4, -2, 1));
        force4.Sub(_deviceChairStapleMu2_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, -4, 2, 1));

        force4.MulDagger(pDeviceData[uiLink4]);
        force4.Ta();
        force4.MulReal(betaOverN * fOmega);
        pForceData[uiLink4].Sub(force4);
    }

    //Nu
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);
    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    {
        deviceSU3 force2 = // deviceSU3::makeSU3Zero();
            _deviceChairStapleNu_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, 4, 2, 1);
        force2.Add(_deviceChairStapleNu_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, 4, -2, 1));
        force2.Sub(_deviceChairStapleNu_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, 4, -2, -1));
        force2.Sub(_deviceChairStapleNu_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, 4, 2, -1));

        force2.Add(_deviceChairStapleNu_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, -4, 2, -1));
        force2.Add(_deviceChairStapleNu_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, -4, -2, -1));
        force2.Sub(_deviceChairStapleNu_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, -4, -2, 1));
        force2.Sub(_deviceChairStapleNu_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, -4, 2, 1));

        force2.MulDagger(pDeviceData[uiLink2]);
        force2.Ta();
        force2.MulReal(betaOverN * fOmega);
        pForceData[uiLink2].Sub(force2);
    }

    //Rho
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);
    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    {
        deviceSU3 force1 = // deviceSU3::makeSU3Zero();
            _deviceChairStapleRho_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, 4, 2, 1);
        force1.Add(_deviceChairStapleRho_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, 4, -2, 1));
        force1.Sub(_deviceChairStapleRho2_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, 4, -2, -1));
        force1.Sub(_deviceChairStapleRho2_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, 4, 2, -1));

        force1.Add(_deviceChairStapleRho2_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, -4, 2, -1));
        force1.Add(_deviceChairStapleRho2_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, -4, -2, -1));
        force1.Sub(_deviceChairStapleRho_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, -4, -2, 1));
        force1.Sub(_deviceChairStapleRho_412(byFieldId, 1, pDeviceData, sSite4, sCenterSite, -4, 2, 1));

        force1.MulDagger(pDeviceData[uiLink1]);
        force1.Ta();
        force1.MulReal(betaOverN * fOmega);
        pForceData[uiLink1].Sub(force1);
    }
}

#pragma endregion

#endif

CActionGaugePlaquetteRotating::CActionGaugePlaquetteRotating()
    : CAction()
    , m_fOmega(F(0.0))
    , m_bCloverEnergy(FALSE)
    , m_bShiftHalfCoord(FALSE)
    , m_bTorus(FALSE)
    , m_uiPlaqutteCount(0)
{
}

void CActionGaugePlaquetteRotating::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        m_fLastEnergy = EnergySingleField(FALSE, pGauge, NULL);
    }
}

void CActionGaugePlaquetteRotating::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    m_fBetaOverN = CCommonData::m_fBeta / static_cast<DOUBLE>(GetDefaultMatrixN());

#if !_CLG_DOUBLEFLOAT
    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

    DOUBLE fOmega = 0.1;
    param.FetchValueDOUBLE(_T("Omega"), fOmega);
    m_fOmega = fOmega;
    CCommonData::m_fOmega = fOmega;
#else
    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

    Real fOmega = F(0.1);
    param.FetchValueReal(_T("Omega"), fOmega);
    m_fOmega = fOmega;
    CCommonData::m_fOmega = fOmega;
#endif


    //TArray<INT> centerArray;
    //param.FetchValueArrayINT(_T("Center"), centerArray);
    //if (centerArray.Num() > 3)
    //{
    //    SSmallInt4 sCenter;
    //    sCenter.x = static_cast<SBYTE>(centerArray[0]);
    //    sCenter.y = static_cast<SBYTE>(centerArray[1]);
    //    sCenter.z = static_cast<SBYTE>(centerArray[2]);
    //    sCenter.w = static_cast<SBYTE>(centerArray[3]);
    //    CCommonData::m_sCenter = sCenter;
    //}

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
}

void CActionGaugePlaquetteRotating::SetBeta(DOUBLE fBeta)
{
    CCommonData::m_fBeta = fBeta;
    m_fBetaOverN = fBeta / static_cast<DOUBLE>(GetDefaultMatrixN());
}

UBOOL CActionGaugePlaquetteRotating::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
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

#if old_cold

    preparethread;


    if (!m_bShiftHalfCoord)
    {
        _kernelAddForce4PlaqutteTermSU3_XY << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData, CCommonData::m_sCenter,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

        _kernelAddForceChairTermSU3_Term1 << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData, CCommonData::m_sCenter,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term2 << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData, CCommonData::m_sCenter,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term3 << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData, CCommonData::m_sCenter,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term4 << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData, CCommonData::m_sCenter,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term5 << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData, CCommonData::m_sCenter,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);
    }
    else
    {

        _kernelAddForce4PlaqutteTermSU3_XYZ_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData, CCommonData::m_sCenter,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

        
        _kernelAddForceChairTermSU3_Term1_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData, CCommonData::m_sCenter,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term2_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData, CCommonData::m_sCenter,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);
        
        _kernelAddForceChairTermSU3_Term3_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData, CCommonData::m_sCenter,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term4_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData, CCommonData::m_sCenter,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term5_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData, CCommonData::m_sCenter,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);
    }

#endif

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
        CalculateForceOnGaugeTorus(pGaugeSU3, pForceSU3);
    }

    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
DOUBLE CActionGaugePlaquetteRotating::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

#if old_code

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

        _kernelAdd4PlaqutteTermSU3_Shifted << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            CCommonData::m_sCenter,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        
        _kernelAddChairTermSU3_Term12_Shifted << <block, threads >> > (
            pGaugeSU3->m_byFieldId, 
            pGaugeSU3->m_pDeviceData, 
            CCommonData::m_sCenter,
            m_fBetaOverN, 
            m_fOmega, 
            _D_RealThreadBuffer);
        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        
        _kernelAddChairTermSU3_Term34_Shifted << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            CCommonData::m_sCenter,
            m_fBetaOverN,
            m_fOmega,
            _D_RealThreadBuffer);
        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _kernelAddChairTermSU3_Term5_Shifted << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            CCommonData::m_sCenter,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);
        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    }
    else
    {
        //======== this is only for test ================
        //_kernelAdd4PlaqutteTermSU3_Test << <block, threads >> > (
        //    pGaugeSU3->m_byFieldId,
        //    pGaugeSU3->m_pDeviceData,
        //    CCommonData::m_sCenter,
        //    m_fBetaOverN,
        //    m_fOmega * m_fOmega,
        //    _D_RealThreadBuffer);

        _kernelAdd4PlaqutteTermSU3 << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
            CCommonData::m_sCenter,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _kernelAddChairTermSU3_Term12 << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            CCommonData::m_sCenter,
            m_fBetaOverN,
            m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _kernelAddChairTermSU3_Term34 << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            CCommonData::m_sCenter,
            m_fBetaOverN,
            m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _kernelAddChairTermSU3_Term5 << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            CCommonData::m_sCenter,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    }

#endif

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    if (NULL == pGaugeSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3 now.\n"));
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

    return m_fNewEnergy;
}

DOUBLE CActionGaugePlaquetteRotating::XYTerm1(const class CFieldGauge* pGauge)
{
    preparethread;
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    _kernelAdd4PlaqutteTermSU3_Test << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        m_fBetaOverN,
        m_fOmega * m_fOmega,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

DOUBLE CActionGaugePlaquetteRotating::XYTerm2(const class CFieldGauge* pGauge)
{
    preparethread;
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    _kernelAdd4PlaqutteTermSU3 << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        m_fBetaOverN,
        m_fOmega * m_fOmega,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}


/**
* In the case of Dirichlet, 
* m_bShiftHalfCoord = FALSE
* m_bCloverEnergy = FALSE
* 
*/
void CActionGaugePlaquetteRotating::EnergyDirichlet(const class CFieldGaugeSU3* pGaugeSU3)
{
    assert(!m_bShiftHalfCoord && !m_bCloverEnergy);
    m_fNewEnergy = pGaugeSU3->CalculatePlaqutteEnergy(m_fBetaOverN);

    preparethread;
    appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

    //======== this is only for test ================
    //_kernelAdd4PlaqutteTermSU3_Test << <block, threads >> > (
    //    pGaugeSU3->m_byFieldId,
    //    pGaugeSU3->m_pDeviceData,
    //    m_fBetaOverN,
    //    m_fOmega * m_fOmega,
    //    _D_RealThreadBuffer);

    _kernelAdd4PlaqutteTermSU3 << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        m_fBetaOverN,
        m_fOmega * m_fOmega,
        _D_RealThreadBuffer);

    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    _kernelAddChairTermSU3_Term12 << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        m_fBetaOverN,
        m_fOmega,
        _D_RealThreadBuffer);

    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    _kernelAddChairTermSU3_Term34 << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        m_fBetaOverN,
        m_fOmega,
        _D_RealThreadBuffer);

    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    _kernelAddChairTermSU3_Term5 << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        m_fBetaOverN,
        m_fOmega * m_fOmega,
        _D_RealThreadBuffer);

    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

/**
* In the case of ProjectivePlane,
* m_bShiftHalfCoord = TRUE
* m_bCloverEnergy = TRUE
*
*/
void CActionGaugePlaquetteRotating::EnergyProjectivePlane(const class CFieldGaugeSU3* pGaugeSU3)
{
    assert(m_bShiftHalfCoord && m_bCloverEnergy);

    m_fNewEnergy = pGaugeSU3->CalculatePlaqutteEnergyUseClover(m_fBetaOverN);
    preparethread;
    appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

    _kernelAdd4PlaqutteTermSU3_Shifted << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        m_fBetaOverN,
        m_fOmega * m_fOmega,
        _D_RealThreadBuffer);

    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);


    _kernelAddChairTermSU3_Term12_Shifted << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        m_fBetaOverN,
        m_fOmega,
        _D_RealThreadBuffer);
    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);


    _kernelAddChairTermSU3_Term34_Shifted << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        m_fBetaOverN,
        m_fOmega,
        _D_RealThreadBuffer);
    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    _kernelAddChairTermSU3_Term5_Shifted << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        m_fBetaOverN,
        m_fOmega * m_fOmega,
        _D_RealThreadBuffer);
    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

void CActionGaugePlaquetteRotating::EnergyTorus(const class CFieldGaugeSU3* pGaugeSU3)
{
    //assert(!m_bShiftHalfCoord);

    m_fNewEnergy = pGaugeSU3->CalculatePlaqutteEnergy(m_fBetaOverN);

    preparethread;
    appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

    if (m_bShiftHalfCoord)
    {
        _kernelAdd4PlaqutteTermSU3_Shifted << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);


        _kernelAddChairTermSU3_Term12_Shifted << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            m_fBetaOverN,
            m_fOmega,
            _D_RealThreadBuffer);
        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);


        _kernelAddChairTermSU3_Term34_Shifted << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            m_fBetaOverN,
            m_fOmega,
            _D_RealThreadBuffer);
        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _kernelAddChairTermSU3_Term5_Shifted << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);
        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
    }
    else
    {
        _kernelAdd4PlaqutteTermSU3 << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _kernelAddChairTermSU3_Term12 << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            m_fBetaOverN,
            m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _kernelAddChairTermSU3_Term34 << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            m_fBetaOverN,
            m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

        _kernelAddChairTermSU3_Term5 << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);

        m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
    }
}

void CActionGaugePlaquetteRotating::CalculateForceOnGaugeDirichlet(const class CFieldGaugeSU3* pGaugeSU3, class CFieldGaugeSU3* pForceSU3) const
{
    preparethread;

    _kernelAddForce4PlaqutteTermSU3_XY << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

    _kernelAddForceChairTermSU3_Term1 << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term2 << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term3 << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term4 << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term5 << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);
}

void CActionGaugePlaquetteRotating::CalculateForceOnGaugeProjectivePlane(const class CFieldGaugeSU3* pGaugeSU3, class CFieldGaugeSU3* pForceSU3) const
{
    preparethread;
    _kernelAddForce4PlaqutteTermSU3_XYZ_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

    _kernelAddForceChairTermSU3_Term1_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term2_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term3_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term4_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term5_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, FALSE, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);
}

void CActionGaugePlaquetteRotating::CalculateForceOnGaugeTorus(const class CFieldGaugeSU3* pGaugeSU3, class CFieldGaugeSU3* pForceSU3) const
{
    preparethread;

    if (m_bShiftHalfCoord)
    {
        _kernelAddForce4PlaqutteTermSU3_XYZ_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

        _kernelAddForceChairTermSU3_Term1_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term2_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term3_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term4_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term5_Shifted << <block, threads >> > (pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);
    }
    else
    {
        _kernelAddForce4PlaqutteTermSU3_XY << <block, threads >> > (pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

        _kernelAddForceChairTermSU3_Term1 << <block, threads >> > (pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term2 << <block, threads >> > (pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term3 << <block, threads >> > (pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term4 << <block, threads >> > (pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

        _kernelAddForceChairTermSU3_Term5 << <block, threads >> > (pGaugeSU3->m_byFieldId, TRUE, pGaugeSU3->m_pDeviceData,
            pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);
    }
}


//Real CActionGaugePlaquetteRotating::GetEnergyPerPlaqutte() const
//{
//    return m_pOwner->m_pGaugeField->CalculatePlaqutteEnergy(m_fBetaOverN) / m_uiPlaqutteCount;
//}

void CActionGaugePlaquetteRotating::SetOmega(DOUBLE fOmega)
{ 
    m_fOmega = fOmega; 
    CCommonData::m_fOmega = fOmega;
}

//void CActionGaugePlaquetteRotating::SetCenter(const SSmallInt4 &newCenter) 
//{
//    CCommonData::m_sCenter = newCenter;
//}

CCString CActionGaugePlaquetteRotating::GetInfos(const CCString &tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("Beta : ") + appToString(CCommonData::m_fBeta) + _T("\n");
    sRet = sRet + tab + _T("Omega : ") + appToString(m_fOmega) + _T("\n");

    sRet = sRet + tab + _T("ShiftCenter : ") + (m_bShiftHalfCoord ? _T("1") : _T("0")) + _T("\n");
    sRet = sRet + tab + _T("Clover : ") + (m_bCloverEnergy ? _T("1") : _T("0")) + _T("\n");
    sRet = sRet + tab + _T("Torus : ") + (m_bTorus ? _T("1") : _T("0")) + _T("\n");

    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================