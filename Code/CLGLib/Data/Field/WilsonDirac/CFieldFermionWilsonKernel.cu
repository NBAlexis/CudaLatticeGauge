//=============================================================================
// FILENAME : CFieldFermionWilsonKernel.cu
// 
// DESCRIPTION:
// This is the device implementations of Wilson fermion
//
// This implementation assumes SU3 and square lattice
//
// REVISION:
//  [mm/dd/yy]
//  [08/04/2025 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "Data/Field/Gauge/CFieldGaugeLink.h"
#include "CFieldFermionWilsonKernel.h"


__BEGIN_NAMESPACE

#pragma region Clover

#pragma region Clover kernel

/**
* kaisw = kai * csw
*
*
* Sigma:
* SIGMA12E, 20
* SIGMA23E, 21
* SIGMA31E, 22
* SIGMA41 (actually it was SIGMA14), 17
* SIGMA42, 18
* SIGMA43, 19
*
* 14, 16, 17, 15, 18, 19
* +   -   -   +   -   -
*
*/
__device__ __constant__ constexpr SCHAR _gamma_munu[6][2] = {
    {20, 1},
    {22, -1},
    {17, 1},
    {21, 1},
    {18, 1},
    {19, 1},
};

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareCloverSU3(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pFmunu,
    deviceWilsonVectorSU3* pResultData,
    DOUBLE kaisw,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalE(6);

    const gammaMatrix& gamma5 = __chiralGamma[GAMMA5];
    //-kappa csw sum_{mu>nu} F, with F = (P...)/8i
    //8i was not divided in calculation of F
    //+(0, +i kaisw/8)
    kaisw = kaisw * _gamma_munu[elementIdx][1] * (0.125);

    deviceWilsonVectorSU3 res = pDeviceData[uiSiteIndex];
    if (bDDagger)
    {
        res = gamma5.MulWilsonC(res);
    }
    res = pFmunu[elementIdx * _DC_Volume + uiSiteIndex].MulWilsonVector(__chiralGamma[_gamma_munu[elementIdx][0]].MulWilsonC(res));
    res.MulComp(_make_cuComplex(F(0.0), static_cast<Real>(kaisw)));

    if (bDDagger)
    {
        res = gamma5.MulWilsonC(res);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        res.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        res.MulComp(cCoeff);
        break;
    default:
        break;
    }

    #pragma unroll
    for (BYTE i = 0U; i < 6U; ++i)
    {
        if (i == elementIdx)
        {
            pResultData[uiSiteIndex].Add(res);
        }
        __syncthreads();
    }
}

/**
* calculate (sigma munu phi) phi^dagger
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelForceFermionWilsonSquareCloverSU3_PrepareSigmaMunu(
    const deviceWilsonVectorSU3* __restrict__ phi,
    const deviceWilsonVectorSU3* __restrict__ phid,
    deviceSU3* pResPhidPhi)
{
    intokernalE(6);
    //Improve-1 (I6): site-major storage -- the sigma index is the fast index
    //inside one site slot so the buffer can carry a per-site halo tail managed
    //by the generic halo machinery (was: elementIdx * _DC_Volume + uiSiteIndex).
    const UINT residx = uiSiteIndex * 6U + elementIdx;
    deviceWilsonVectorSU3 right = __chiralGamma[_gamma_munu[elementIdx][0]].MulWilsonC(phi[uiSiteIndex]);
    pResPhidPhi[residx] = _makeContract<deviceSU3, deviceWilsonVectorSU3>(phid[uiSiteIndex], right);

    right = __chiralGamma[_gamma_munu[elementIdx][0]].MulWilsonC(phid[uiSiteIndex]);
    _add(pResPhidPhi[residx], _makeContract<deviceSU3, deviceWilsonVectorSU3>(phi[uiSiteIndex], right));

    if (_gamma_munu[elementIdx][1] < 0)
    {
        _oppo(pResPhidPhi[residx]);
    }
}

/**
*   link   p3_1    p3_2    p3_3
*    0      yx-0    zx-1    tx-2
*    1      xy-0    zy-3    ty-4
*    2      xz-1    yz-3    tz-5
*    3      xt-2    yt-4    zt-5
*/
__device__ __constant__ constexpr SCHAR _staple_to_sigma[4][3] = {
    {0, 1, 2},
    {0, 3, 4},
    {1, 3, 5},
    {2, 4, 5}
};
/**
*
*/
__device__ __constant__ constexpr SCHAR _sigma_opposite[4][3] = {
    { 1,  1,  1},
    {-1,  1,  1},
    {-1, -1,  1},
    {-1, -1, -1}
};
__global__ void _CLG_LAUNCH_BOUND
_kernelCloverForceCacheIndex(
    const deviceSU3* __restrict__ pResPhidPhi,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pCachedIndex,
    deviceSU3* pForce,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE kaisw)
{
    // 4 dirs, 6 plaqs, 
    // forward-backward go together
    // 4 types comes together
    intokernalEDir(6U);
    const BYTE stapleidx = elementIdx >> 1U; //0,1,2
    const BYTE forwardbackward = (elementIdx & 1U); //forward-0 backward-1
    //Improve-1 (I6): site-major -- slot = site * 6 + sigma (was sigma * V + site).
    const SCHAR bySigmaIdx = _staple_to_sigma[dir][stapleidx];
    //const UINT uiCacheIndex = stapleidx * 6U + uiLinkIndex * 18U + forwardbackward * 3U; //18 for each link, for this link there are 6 staples, and each has 3 index, 6 staples are forward&backward, therefore 6 index
    const UINT uiCacheIndex = elementIdx * 3U + uiLinkIndex * 18U; //18 for each link, for this link there are 6 staples, and each has 3 index, 6 staples are forward&backward, therefore 6 index

    //0.125 is added because we need Fmunu = (Gmunu-Gmunud)/8i, but we calculate it as (Gmunu-Gmunud)
    kaisw = kaisw * (0.125) * (1 - 2 * forwardbackward) * _sigma_opposite[dir][stapleidx];
    const CLGComplex csw = _make_cuComplex(F(0.0), static_cast<Real>(kaisw));

    //add1 = o - - -
    //add2 = - o - -
    //add3 = - - o -
    //add4 = - - - o
    SIndex link = pCachedIndex[uiCacheIndex];
    SSmallInt4 startsite = __deviceSiteIndexToInt4(uiSiteIndex);
    deviceSU3 toAdd1 = pResPhidPhi[uiSiteIndex * 6U + bySigmaIdx];
    deviceSU3 toAdd2 = _deviceGetGaugeBCT(byGaugeFieldId, pGauge, link);
    if (link.NeedToDagger())
    {
        _dagger(toAdd2);
        startsite.m_byData4[link.m_byDir] = startsite.m_byData4[link.m_byDir] - 1;
    }
    else
    {
        startsite.m_byData4[link.m_byDir] = startsite.m_byData4[link.m_byDir] + 1;
    }
    _mul(toAdd1, toAdd2);
    _add(toAdd1, _mulC(toAdd2, pResPhidPhi[__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(startsite)].m_uiSiteIndex * 6U + bySigmaIdx]));

    link = pCachedIndex[uiCacheIndex + 1];
    const deviceSU3& tomul1 = _deviceGetGaugeBCT(byGaugeFieldId, pGauge, link);
    if (link.NeedToDagger())
    {
        _muldag(toAdd1, tomul1);
        _muldag(toAdd2, tomul1);
        startsite.m_byData4[link.m_byDir] = startsite.m_byData4[link.m_byDir] - 1;
    }
    else
    {
        _mul(toAdd1, tomul1);
        _mul(toAdd2, tomul1);
        startsite.m_byData4[link.m_byDir] = startsite.m_byData4[link.m_byDir] + 1;
    }
    _add(toAdd1, _mulC(toAdd2, pResPhidPhi[__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(startsite)].m_uiSiteIndex * 6U + bySigmaIdx])); //add3 is free
    link = pCachedIndex[uiCacheIndex + 2];
    const deviceSU3& tomul2 = _deviceGetGaugeBCT(byGaugeFieldId, pGauge, link);
    if (link.NeedToDagger())
    {
        _muldag(toAdd1, tomul2);
        _muldag(toAdd2, tomul2);
        startsite.m_byData4[link.m_byDir] = startsite.m_byData4[link.m_byDir] - 1;
    }
    else
    {
        _mul(toAdd1, tomul2);
        _mul(toAdd2, tomul2);
        startsite.m_byData4[link.m_byDir] = startsite.m_byData4[link.m_byDir] + 1;
    }
    _mul(toAdd2, pResPhidPhi[__idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(startsite)].m_uiSiteIndex * 6U + bySigmaIdx]);

    //===================================
    //add
    _add(toAdd1, toAdd2);
    _mul(toAdd1, csw);
    //_mul(toAdd1, _make_cuComplex(F(0.0), static_cast<Real>(kaisw * (-0.03125) * (1 - 2 * forwardbackward) * _sigma_opposite[dir][stapleidx])));

    if (0 == elementIdx)
    {
        _add(pForce[uiLinkIndex], toAdd1);
    }
    __syncthreads();
    if (1 == elementIdx)
    {
        _add(pForce[uiLinkIndex], toAdd1);
    }
    __syncthreads();
    if (2 == elementIdx)
    {
        _add(pForce[uiLinkIndex], toAdd1);
    }
    __syncthreads();
    if (3 == elementIdx)
    {
        _add(pForce[uiLinkIndex], toAdd1);
    }
    __syncthreads();
    if (4 == elementIdx)
    {
        _add(pForce[uiLinkIndex], toAdd1);
    }
    __syncthreads();
    if (5 == elementIdx)
    {
        _add(pForce[uiLinkIndex], toAdd1);
    }
    //__syncthreads();
}

#pragma endregion

void CFieldFermionWilsonKernel::DOperatorClover(deviceWilsonVectorSU3* pTarget, const deviceWilsonVectorSU3* pSource,
    const deviceSU3* pGauge, const deviceSU3* pFmunu, BYTE byFieldId, BYTE byGaugeFieldId, DOUBLE fCoef,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethreadE(6);
    _LAUNCH_KERNEL(_kernelDFermionWilsonSquareCloverSU3, block, threads,
        pSource,
        pFmunu,
        pTarget,
        fCoef,
        byFieldId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}


void CFieldFermionWilsonKernel::PrepareSigmaMunu(const deviceWilsonVectorSU3* phi, const deviceWilsonVectorSU3* phid, deviceSU3* res)
{
    preparethreadE(6);
    _LAUNCH_KERNEL(_kernelForceFermionWilsonSquareCloverSU3_PrepareSigmaMunu, block, threads,
        phi,
        phid,
        res
    );
}

void CFieldFermionWilsonKernel::CloverForce(const deviceSU3* phidphi, const deviceSU3* pGauge, deviceSU3* force, BYTE byFermionId, BYTE byGaugeId, DOUBLE fCoef)
{
    preparethreadEDir(6);
    _LAUNCH_KERNEL(_kernelCloverForceCacheIndex, block, threads,
        phidphi,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pStappleCache[byGaugeId],
        force,
        byFermionId,
        byGaugeId,
        fCoef
    );
}

#pragma endregion


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================