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


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteRotating)


#pragma region kernels

/**
* This is slower, just for testing
* directly calculate Retr[1 - \hat{U}]
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermSU3_Test(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    Real betaOverN, Real fOmegaSq, 
    Real* results)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    if (__idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx].IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        return;
    }

    Real fXSq = (sSite4.x - sCenterSite.x);
    fXSq = fXSq * fXSq;
    Real fYSq = (sSite4.y - sCenterSite.y);
    fYSq = fYSq * fYSq;

    //======================================================
    //4-plaqutte terms
    //Omega^2 x^2 Retr[1 - U_1,4]
    Real fU14 = fOmegaSq * fXSq * _device4PlaqutteTerm(pDeviceData, 0, 3, uiBigIdx);

    //Omega^2 y^2 Retr[1 - U_2,4]
    Real fU24 = fOmegaSq * fYSq * _device4PlaqutteTerm(pDeviceData, 1, 3, uiBigIdx);

    //Omega^2 (x^2 + y^2) Retr[1 - U_3,4]
    Real fU34 = fOmegaSq * (fXSq + fYSq) * _device4PlaqutteTerm(pDeviceData, 2, 3, uiBigIdx);

    results[uiSiteIndex] = (fU14 + fU24 + fU34) * betaOverN;
}

/**
* Using plaqutte and (f(n)+f(n+mu)+f(n+nu)+f(n+mu+nu))/4 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermSU3(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
    SSmallInt4 sCenterSite,
    Real betaOverN, Real fOmegaSq,
    Real* results)
{
    //intokernalInt4;
    SSmallInt4 sSite4;
    UINT _ixy = (threadIdx.x + blockIdx.x * blockDim.x);
    UINT _iz_idx = (threadIdx.y + blockIdx.y * blockDim.y);

    sSite4.x = static_cast<SBYTE> (_ixy / _DC_Lx);
    sSite4.y = static_cast<SBYTE> (_ixy % _DC_Lx);
    sSite4.z = static_cast<SBYTE>(_iz_idx / 3);
    sSite4.w = static_cast<SBYTE>(threadIdx.z + blockIdx.z * blockDim.z);
    UINT uiSiteIndex = _ixy * _DC_GridDimZT + sSite4.z * _DC_Lt + sSite4.w;
    BYTE idx0 = _iz_idx % 3;

    UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    UINT plaqLength = __idx->m_pSmallData[CIndexData::kPlaqLengthIdx];
    UINT plaqCountAll = __idx->m_pSmallData[CIndexData::kPlaqPerSiteIdx] * plaqLength;
    
    //i=0: 12
    //  1: 13
    //  2: 14
    //  3: 23
    //  4: 24
    //  5: 34
    //0->2, 1->4, 2->5
    BYTE idx = (2 == idx0) ? 5 : ((idx0 + 1) * 2);

    //Real resThisThread = F(0.0);

    //========================================
    //find plaqutte 1-4, or 2-4, or 3-4
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

    atomicAdd(&results[uiSiteIndex], betaOverN * fOmegaSq * (F(3.0) - toAdd.ReTr()) * _deviceFi(sCenterSite, uiN, idx0, idx0, 3));

}


/**
* Split into 3 functions to avoid max-register problem
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term12(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    Real betaOverN, Real fOmega,
    Real* results)
{
    intokernalInt4;

    UINT uiN = __idx->_deviceGetBigIndex(sSite4);

    if (__idx->m_pDeviceIndexPositionToSIndex[1][uiN].IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        return;
    }

    betaOverN = F(0.125) * betaOverN;
    Real fXOmega = (sSite4.x - sCenterSite.x) * fOmega;

    //===============
    //+x Omega V412
    Real fV412 = fXOmega * _deviceChairTerm(pDeviceData, 3, 0, 1, uiN);

    //===============
    //+x Omega V432
    Real fV432 = fXOmega * _deviceChairTerm(pDeviceData, 3, 2, 1, uiN);

    results[uiSiteIndex] = (fV412 + fV432) * betaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term34(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    Real betaOverN, Real fOmega,
    Real* results)
{
    intokernalInt4;

    UINT uiN = __idx->_deviceGetBigIndex(sSite4);

    if (__idx->m_pDeviceIndexPositionToSIndex[1][uiN].IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        return;
    }

    betaOverN = F(0.125) * betaOverN;
    Real fYOmega = -(sSite4.y - sCenterSite.y) * fOmega;

    //===============
    //-y Omega V421
    Real fV421 = fYOmega * _deviceChairTerm(pDeviceData, 3, 1, 0, uiN);

    //===============
    //-y Omega V431
    Real fV431 = fYOmega * _deviceChairTerm(pDeviceData, 3, 2, 0, uiN);

    results[uiSiteIndex] = (fV421 + fV431) * betaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term5(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    Real betaOverN, Real fOmegaSq,
    Real* results)
{
    intokernalInt4;

    UINT uiN = __idx->_deviceGetBigIndex(sSite4);

    if (__idx->m_pDeviceIndexPositionToSIndex[1][uiN].IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        return;
    }

    betaOverN = F(0.125) * betaOverN;
    Real fXYOmega2 = (sSite4.x - sCenterSite.x) * (sSite4.y - sCenterSite.y) * fOmegaSq;

    //===============
    //+Omega^2 xy V142
    Real fV142 = fXYOmega2 * _deviceChairTerm(pDeviceData, 0, 3, 1, uiN);

    results[uiSiteIndex] = fV142 * betaOverN;
}

/**
* 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermSU3_XYZ(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmegaSq)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(-0.5);
    deviceSU3 plaqSum = deviceSU3::makeSU3Zero();

    for (UINT idir = 0; idir < 3; ++idir)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, idir))
        {
            continue;
        }
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        //mu = idir, nu = 4, i = mu
        deviceSU3 stap(_deviceStapleTerm123(pDeviceData, sCenterSite, sSite4, fOmegaSq, uiBigIdx, idir, 3, idir));
        deviceSU3 force(pDeviceData[linkIndex]);
        force.MulDagger(stap);
        force.Ta();
        force.MulReal(betaOverN);
        pForceData[linkIndex].Add(force);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermSU3_T(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmegaSq)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    //UINT uiDir = _DC_Dir;

    betaOverN = betaOverN * F(-0.5);
    deviceSU3 plaqSum = deviceSU3::makeSU3Zero();

    BYTE idir = 3;
    if (__idx->_deviceIsBondOnSurface(uiBigIdx, idir))
    {
        return;
    }
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

    //mu = idir, nu = i = sum _1-3
    deviceSU3 stap(_deviceStapleTerm4(pDeviceData, sCenterSite, sSite4, fOmegaSq, uiBigIdx, idir, 0));
    stap.Add(_deviceStapleTerm4(pDeviceData, sCenterSite, sSite4, fOmegaSq, uiBigIdx, idir, 1));
    stap.Add(_deviceStapleTerm123(pDeviceData, sCenterSite, sSite4, fOmegaSq, uiBigIdx, idir, 2, 2));
    deviceSU3 force(pDeviceData[linkIndex]);
    force.MulDagger(stap);
    force.Ta();
    force.MulReal(betaOverN);
    pForceData[linkIndex].Add(force);
}

/**
* Split to 15 functions to avoid max-regcount
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term1_1(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+x Omega V412
    //add force for mu=4
    UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        deviceSU3 staple_term1_4 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            3, 0, 1, 0);
        deviceSU3 force4(pDeviceData[uiLink4]);
        force4.MulDagger(staple_term1_4);
        force4.Ta();
        force4.MulReal(betaOverN);
        pForceData[uiLink4].Add(force4);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term1_2(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+x Omega V412
    //add force for mu=4
    UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    {
        deviceSU3 staple_term1_2 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            1, 0, 3, 0);
        deviceSU3 force2(pDeviceData[uiLink2]);
        force2.MulDagger(staple_term1_2);
        force2.Ta();
        force2.MulReal(betaOverN);
        pForceData[uiLink2].Add(force2);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term1_3(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+x Omega V412
    //add force for mu=4
    UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    {
        deviceSU3 staple_term1_1 = _deviceStapleChairTerm2(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            3, 0, 1, 0);
        deviceSU3 force1(pDeviceData[uiLink1]);
        force1.MulDagger(staple_term1_1);
        force1.Ta();
        force1.MulReal(betaOverN);
        pForceData[uiLink1].Add(force1);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term2_1(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+x Omega V432
    //add force for mu=4
    UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        deviceSU3 staple_term2_4 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            3, 2, 1, 0);
        deviceSU3 force4(pDeviceData[uiLink4]);
        force4.MulDagger(staple_term2_4);
        force4.Ta();
        force4.MulReal(betaOverN);
        pForceData[uiLink4].Add(force4);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term2_2(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+x Omega V432
    //add force for mu=4
    UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    {
        deviceSU3 staple_term2_2 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            1, 2, 3, 0);
        deviceSU3 force2(pDeviceData[uiLink2]);
        force2.MulDagger(staple_term2_2);
        force2.Ta();
        force2.MulReal(betaOverN);
        pForceData[uiLink2].Add(force2);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term2_3(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //+x Omega V432
    //add force for mu=4
    UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 2))
    {
        deviceSU3 staple_term2_3 = _deviceStapleChairTerm2(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            3, 2, 1, 0);
        deviceSU3 force3(pDeviceData[uiLink3]);
        force3.MulDagger(staple_term2_3);
        force3.Ta();
        force3.MulReal(betaOverN);
        pForceData[uiLink3].Add(force3);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term3_1(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //-y Omega V421
    //add force for mu=4
    UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        deviceSU3 staple_term3_4 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            3, 1, 0, 1);
        deviceSU3 force4(pDeviceData[uiLink4]);
        force4.MulDagger(staple_term3_4);
        force4.Ta();
        force4.MulReal(betaOverN);
        pForceData[uiLink4].Add(force4);
    }

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term3_2(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //-y Omega V421
    //add force for mu=4
    UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    {
        deviceSU3 staple_term3_1 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            0, 1, 3, 1);
        deviceSU3 force1(pDeviceData[uiLink1]);
        force1.MulDagger(staple_term3_1);
        force1.Ta();
        force1.MulReal(betaOverN);
        pForceData[uiLink1].Add(force1);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term3_3(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //-y Omega V421
    //add force for mu=4
    UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    {
        deviceSU3 staple_term3_2 = _deviceStapleChairTerm2(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            3, 1, 0, 1);
        deviceSU3 force2(pDeviceData[uiLink2]);
        force2.MulDagger(staple_term3_2);
        force2.Ta();
        force2.MulReal(betaOverN);
        pForceData[uiLink2].Add(force2);
    }

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term4_1(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //-y Omega V431
    //add force for mu=4
    UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        deviceSU3 staple_term4_4 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            3, 2, 0, 1);
        deviceSU3 force4(pDeviceData[uiLink4]);
        force4.MulDagger(staple_term4_4);
        force4.Ta();
        force4.MulReal(betaOverN);
        pForceData[uiLink4].Add(force4);
    }

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term4_2(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //-y Omega V431
    //add force for mu=4
    UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    {
        deviceSU3 staple_term4_1 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            0, 2, 3, 1);
        deviceSU3 force1(pDeviceData[uiLink1]);
        force1.MulDagger(staple_term4_1);
        force1.Ta();
        force1.MulReal(betaOverN);
        pForceData[uiLink1].Add(force1);
    }

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term4_3(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmega)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmega * F(0.125);

    //===============
    //-y Omega V431
    //add force for mu=4
    UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 2))
    {
        deviceSU3 staple_term4_3 = _deviceStapleChairTerm2(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            3, 2, 0, 1);
        deviceSU3 force3(pDeviceData[uiLink3]);
        force3.MulDagger(staple_term4_3);
        force3.Ta();
        force3.MulReal(betaOverN);
        pForceData[uiLink3].Add(force3);
    }

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term5_1(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmegaSq)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmegaSq * F(0.125);

    //===============
    //+Omega^2 xy V142
    UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 0))
    {
        deviceSU3 staple_term5_1 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            0, 3, 1, 2);
        deviceSU3 force1(pDeviceData[uiLink1]);
        force1.MulDagger(staple_term5_1);
        force1.Ta();
        force1.MulReal(betaOverN);
        pForceData[uiLink1].Add(force1);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term5_2(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmegaSq)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmegaSq * F(0.125);

    //===============
    //+Omega^2 xy V142
    UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 1))
    {
        deviceSU3 staple_term5_2 = _deviceStapleChairTerm1(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            1, 3, 0, 2);
        deviceSU3 force2(pDeviceData[uiLink2]);
        force2.MulDagger(staple_term5_2);
        force2.Ta();
        force2.MulReal(betaOverN);
        pForceData[uiLink2].Add(force2);
    }

}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term5_3(
    const deviceSU3 * __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
    deviceSU3 *pForceData,
    Real betaOverN, Real fOmegaSq)
{
    intokernalInt4;

    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fOmegaSq * F(0.125);

    //===============
    //+Omega^2 xy V142
    UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);

    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
    {
        deviceSU3 staple_term5_4 = _deviceStapleChairTerm2(pDeviceData, sCenterSite, sSite4, uiSiteIndex, uiBigIdx,
            0, 3, 1, 2);
        deviceSU3 force4(pDeviceData[uiLink4]);
        force4.MulDagger(staple_term5_4);
        force4.Ta();
        force4.MulReal(betaOverN);
        pForceData[uiLink4].Add(force4);
    }

}

#pragma endregion


CActionGaugePlaquetteRotating::CActionGaugePlaquetteRotating()
    : CAction()
    , m_uiPlaqutteCount(0)
    , m_fLastEnergy(F(0.0))
    , m_fNewEnergy(F(0.0))
    , m_fBetaOverN(F(0.1))
    , m_fOmega(F(0.0))
{
}

void CActionGaugePlaquetteRotating::PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        m_fLastEnergy = Energy(FALSE, pGauge, NULL);
    }
}

void CActionGaugePlaquetteRotating::OnFinishTrajectory(UBOOL bAccepted)
{
    if (bAccepted)
    {
        m_fLastEnergy = m_fNewEnergy;
    }
}

void CActionGaugePlaquetteRotating::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    m_pOwner = pOwner;
    m_byActionId = byId;
    Real fBeta = 0.1f;
    param.FetchValueReal(_T("Beta"), fBeta);
    CCommonData::m_fBeta = fBeta;
    if (NULL != pOwner->m_pGaugeField && EFT_GaugeSU3 == pOwner->m_pGaugeField->GetFieldType())
    {
        fBeta = fBeta / F(3.0);
    }
    m_fBetaOverN = fBeta;
    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

    Real fOmega = 0.1f;
    param.FetchValueReal(_T("Omega"), fOmega);
    m_fOmega = fOmega;
    CCommonData::m_fOmega = fOmega;

    TArray<INT> centerArray;
    param.FetchValueArrayINT(_T("Center"), centerArray);
    if (centerArray.Num() > 3)
    {
        m_sCenter.x = static_cast<SBYTE>(centerArray[0]);
        m_sCenter.y = static_cast<SBYTE>(centerArray[1]);
        m_sCenter.z = static_cast<SBYTE>(centerArray[2]);
        m_sCenter.w = static_cast<SBYTE>(centerArray[3]);
    }
    CCommonData::m_sCenter = m_sCenter;
}

void CActionGaugePlaquetteRotating::SetBeta(Real fBeta)
{
    CCommonData::m_fBeta = fBeta;
    if (NULL != m_pOwner->m_pGaugeField && EFT_GaugeSU3 == m_pOwner->m_pGaugeField->GetFieldType())
    {
        fBeta = fBeta / F(3.0);
    }
    m_fBetaOverN = fBeta;
}

UBOOL CActionGaugePlaquetteRotating::CalculateForceOnGauge(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    pGauge->CalculateForceAndStaple(pForce, pStaple, m_fBetaOverN);

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    if (NULL == pGaugeSU3 || NULL == pForceSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3 now.\n"));
        return TRUE;
    }

    preparethread;

    _kernelAddForce4PlaqutteTermSU3_XYZ << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

    _kernelAddForce4PlaqutteTermSU3_T << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

    _kernelAddForceChairTermSU3_Term1_1 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term1_2 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term1_3 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term2_1 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term2_2 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term2_3 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term3_1 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term3_2 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term3_3 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term4_1 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term4_2 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term4_3 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega);

    _kernelAddForceChairTermSU3_Term5_1 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

    _kernelAddForceChairTermSU3_Term5_2 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

    _kernelAddForceChairTermSU3_Term5_3 << <block, threads >> >(pGaugeSU3->m_pDeviceData, m_sCenter,
        pForceSU3->m_pDeviceData, m_fBetaOverN, m_fOmega * m_fOmega);

    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
Real CActionGaugePlaquetteRotating::Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }
    m_fNewEnergy = pGauge->CalculatePlaqutteEnergy(m_fBetaOverN);

    const CFieldGaugeSU3D* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3D*>(pGauge);
    if (NULL == pGaugeSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3-Dirichlet now.\n"));
        return m_fNewEnergy;
    }

    preparethread;

    //======== this is only for test ================
    //_kernelAdd4PlaqutteTermSU3_Test << <block, threads >> > (
    //    pGaugeSU3->m_pDeviceData,
    //    m_sCenter,
    //    m_fBetaOverN,
    //    m_fOmega * m_fOmega,
    //    _D_RealThreadBuffer);

    appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

    dim3 block2 = block;
    block2.y = block.y * 3;
    _kernelAdd4PlaqutteTermSU3 << <block2, threads >> > (
            pGaugeSU3->m_pDeviceData, 
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
            m_sCenter,
            m_fBetaOverN,
            m_fOmega * m_fOmega,
            _D_RealThreadBuffer);

    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);


    _kernelAddChairTermSU3_Term12 << <block, threads >> > (
        pGaugeSU3->m_pDeviceData,
        m_sCenter,
        m_fBetaOverN,
        m_fOmega,
        _D_RealThreadBuffer);

    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    _kernelAddChairTermSU3_Term34 << <block, threads >> > (
        pGaugeSU3->m_pDeviceData,
        m_sCenter,
        m_fBetaOverN,
        m_fOmega,
        _D_RealThreadBuffer);

    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    _kernelAddChairTermSU3_Term5 << <block, threads >> > (
        pGaugeSU3->m_pDeviceData,
        m_sCenter,
        m_fBetaOverN,
        m_fOmega * m_fOmega,
        _D_RealThreadBuffer);

    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    return m_fNewEnergy;
}

//Real CActionGaugePlaquetteRotating::GetEnergyPerPlaqutte() const
//{
//    return m_pOwner->m_pGaugeField->CalculatePlaqutteEnergy(m_fBetaOverN) / m_uiPlaqutteCount;
//}

void CActionGaugePlaquetteRotating::SetOmega(Real fOmega) 
{ 
    m_fOmega = fOmega; 
    CCommonData::m_fOmega = fOmega;
}

void CActionGaugePlaquetteRotating::SetCenter(const SSmallInt4 &newCenter) 
{
    m_sCenter = newCenter;
    CCommonData::m_sCenter = newCenter;
}

CCString CActionGaugePlaquetteRotating::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CActionGaugePlaquetteRotating\n");
    sRet = sRet + tab + _T("Beta : ") + appFloatToString(CCommonData::m_fBeta) + _T("\n");
    sRet = sRet + tab + _T("Omega : ") + appFloatToString(m_fOmega) + _T("\n");
    CCString sCenter;
    sCenter.Format(_T("Center: [%d, %d, %d, %d]\n")
        , static_cast<INT>(m_sCenter.x)
        , static_cast<INT>(m_sCenter.y)
        , static_cast<INT>(m_sCenter.z)
        , static_cast<INT>(m_sCenter.w));
    sRet = sRet + tab + sCenter;
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================