//=============================================================================
// FILENAME : CFieldGaugeSU3.cu
// 
// DESCRIPTION:
// This is the device implementations of gauge SU3
//
// The SU3 Matrix is
// 0 1 2
// 3 4 5
// 6 7 8
//
// Number of threads: < 1024
// Number of blocks: V / 1024
//
// threadIdx.xyz = xyz, and we loop for t and dir
//
// REVISION:
//  [12/4/2018 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldGaugeSU3)

#pragma region Kernels

/**
* Initial SU3 Field with a value
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialSU3Feield(deviceSU3 *pDevicePtr, EFieldInitialType eInitialType)
{
    deviceSU3 id = deviceSU3::makeSU3Id();
    deviceSU3 zero = deviceSU3::makeSU3Zero();

    intokernaldir;
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        switch (eInitialType)
        {
        case EFIT_Zero:
        {
            pDevicePtr[uiLinkIndex] = zero;
        }
        break;
        case EFIT_Identity:
        {
            pDevicePtr[uiLinkIndex] = id;
        }
        break;
        case EFIT_Random:
        {
            pDevicePtr[uiLinkIndex] = deviceSU3::makeSU3Random(_deviceGetFatIndex(uiSiteIndex, idir + 1));
            //Real fArg = __cuCargf(pDevicePtr[uiLinkIndex].Tr());
            //pDevicePtr[uiLinkIndex].MulComp(_make_cuComplex(_cos(fArg), -_sin(fArg)));
            //pDevicePtr[uiLinkIndex].Norm();
            //printf("arg=%f\n", __cuCargf(pDevicePtr[uiLinkIndex].Tr()));
        }
        break;
        case EFIT_RandomGenerator:
        {
            pDevicePtr[uiLinkIndex] = deviceSU3::makeSU3RandomGenerator(_deviceGetFatIndex(uiSiteIndex, idir + 1));
        }
        break;
        case EFIT_SumGenerator:
        {
            pDevicePtr[uiLinkIndex] = deviceSU3::makeSU3SumGenerator(F(1.0));
        }
        break;
        default:
        {
            printf("SU3 Field cannot be initialized with this type!");
        }
        break;
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDaggerSU3(deviceSU3* pDevicePtr)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Dagger();

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpySU3A(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x, CLGComplex a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].MulCompC(a));

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMulSU3(deviceSU3* pDevicePtr, const deviceSU3* __restrict__ x, UBOOL bDagger)
{
    gaugeSU3KernelFuncionStart

    if (bDagger)
    {
        pDevicePtr[uiLinkIndex].DaggerMul(x[uiLinkIndex]);
    }
    else
    {
        pDevicePtr[uiLinkIndex].Mul(x[uiLinkIndex]);
    }

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpySU3Real(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x, Real a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].MulRealC(a));

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusSU3(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusSU3(deviceSU3 *pDevicePtr, const deviceSU3* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Sub(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplySU3Complex(deviceSU3 *pDevicePtr, CLGComplex a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].MulComp(a);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplySU3Real(deviceSU3 *pDevicePtr, Real a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].MulReal(a);

    gaugeSU3KernelFuncionEnd
}

/**
* debug kernel
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelPrintSU3(const deviceSU3 * __restrict__ pDeviceData)
{
    gaugeSU3KernelFuncionStart;

    printf("link at %d: %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i, %f+%f i\n",
        uiLinkIndex,
        pDeviceData[uiLinkIndex].m_me[0].x, pDeviceData[uiLinkIndex].m_me[0].y,
        pDeviceData[uiLinkIndex].m_me[1].x, pDeviceData[uiLinkIndex].m_me[1].y,
        pDeviceData[uiLinkIndex].m_me[2].x, pDeviceData[uiLinkIndex].m_me[2].y,
        pDeviceData[uiLinkIndex].m_me[3].x, pDeviceData[uiLinkIndex].m_me[3].y,
        pDeviceData[uiLinkIndex].m_me[4].x, pDeviceData[uiLinkIndex].m_me[4].y,
        pDeviceData[uiLinkIndex].m_me[5].x, pDeviceData[uiLinkIndex].m_me[5].y,
        pDeviceData[uiLinkIndex].m_me[6].x, pDeviceData[uiLinkIndex].m_me[6].y,
        pDeviceData[uiLinkIndex].m_me[7].x, pDeviceData[uiLinkIndex].m_me[7].y,
        pDeviceData[uiLinkIndex].m_me[8].x, pDeviceData[uiLinkIndex].m_me[8].y
    );

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteSU3CacheIndex(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSU3 *pStapleData, //can be NULL
    deviceSU3 *pForceData,
    Real betaOverN)
{
    intokernaldir;

    //Real test_force = F(0.0);
    betaOverN = betaOverN * F(-0.5);
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU3 res = deviceSU3::makeSU3Zero();

        //there are 6 staples, each is sum of two plaquttes
        for (INT i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);

                if (nextlink.NeedToDagger())
                {
                    toAdd.MulDagger(toMul);
                }
                else
                {
                    toAdd.Mul(toMul);
                }
            }
            res.Add(toAdd);
        }
        if (NULL != pStapleData)
        {
            pStapleData[linkIndex] = res;
        }

        //staple calculated
        deviceSU3 force(pDeviceData[linkIndex]);
        force.MulDagger(res);
        //test_force += F(-2.0) * betaOverN * __SU3Generators[8].MulC(force).ImTr();
        force.Ta();
        force.MulReal(betaOverN);

        //force is additive
        pForceData[linkIndex].Add(force);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateOnlyStaple(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSU3 *pStapleData)
{
    intokernaldir;

    //Real test_force = F(0.0);
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU3 res = deviceSU3::makeSU3Zero();

        //there are 6 staples, each is sum of two plaquttes
        for (INT i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);

                if (nextlink.NeedToDagger())
                {
                    toAdd.MulDagger(toMul);
                }
                else
                {
                    toAdd.Mul(toMul);
                }
            }
            res.Add(toAdd);
        }
        pStapleData[linkIndex] = res;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergySU3CacheIndex(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN,
    DOUBLE* results
#else
    Real betaOverN,
    Real* results
#endif
)
{
    intokernal;

#if !_CLG_DOUBLEFLOAT
    DOUBLE resThisThread = 0.0;
#else
    Real resThisThread = F(0.0);
#endif
    const UINT indexSkip = plaqCount * plaqLength * uiSiteIndex;
    for (BYTE i = 0; i < plaqCount; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLength + indexSkip];
        deviceSU3 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

        if (first.NeedToDagger())
        {
            toAdd.Dagger();
        }

        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedIndex[i * plaqLength + j + indexSkip];
            deviceSU3 toMul(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
            if (first.NeedToDagger())
            {
                toAdd.MulDagger(toMul);
            }
            else
            {
                toAdd.Mul(toMul);
            }
        }

#if _CLG_DEBUG
        Real reTr = toAdd.ReTr();
        assert(reTr > -F(1.50001));
        assert(reTr < F(3.00001));
#endif

#if !_CLG_DOUBLEFLOAT
        resThisThread += (3.0 - toAdd.ReTr());
#else
        resThisThread += (F(3.0) - toAdd.ReTr());
#endif
    }

    results[uiSiteIndex] = resThisThread * betaOverN;

    //printf("  ---- energy: thread=%d, res=%f\n", __thread_id, results[__thread_id]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergySU3_UseClover(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE stapleConstant,
    DOUBLE fBetaOverN,
    DOUBLE* results
#else
    Real stapleConstant,
    Real fBetaOverN,
    Real* results
#endif
)
{
    intokernalInt4;

#if !_CLG_DOUBLEFLOAT
    DOUBLE fRes = 0.0;
#else
    Real fRes = F(0.0);
#endif
    for (BYTE byDir1 = 0; byDir1 < _DC_Dir; ++byDir1)
    {
        for (BYTE byDir2 = byDir1 + 1; byDir2 < _DC_Dir; ++byDir2)
        {
            //For any strange boundary condition, U_{mu,nu}=U_{nu,mu}^+ should be TRUE!
            //if (byDir1 != byDir2)
            //{
                fRes += _deviceCloverRetr(pDeviceData, sSite4, __bi(sSite4), byDir1, byDir2, byFieldId);
            //}
        }
    }
#if !_CLG_DOUBLEFLOAT
    fRes = stapleConstant - 0.25 * fRes;
#else
    fRes = stapleConstant - F(0.25) * fRes;
#endif
    results[uiSiteIndex] = fRes * fBetaOverN;// *F(0.5);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyUsingStableSU3(
    const deviceSU3 * __restrict__ pDeviceData,
    const deviceSU3 * __restrict__ pStableData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE stapleConstant,
    DOUBLE betaOverN,
    DOUBLE* results
#else
    Real stapleConstant,
    Real betaOverN,
    Real* results
#endif
)
{
    intokernaldir;

#if !_CLG_DOUBLEFLOAT
    DOUBLE resThisThread = 0.0;
#else
    Real resThisThread = F(0.0);
#endif
    
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        //For each link, there are 6 staples
        resThisThread += (stapleConstant - pDeviceData[linkIndex].MulDaggerC(pStableData[linkIndex]).ReTr());
    }

#if !_CLG_DOUBLEFLOAT
    results[uiSiteIndex] = resThisThread * betaOverN * 0.25;
#else
    results[uiSiteIndex] = resThisThread * betaOverN * F(0.25);
#endif

    //printf("  ---- energy: thread=%d, res=%f\n", __thread_id, results[__thread_id]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelExpMultSU3RealQ(
    const deviceSU3 * __restrict__ pMyDeviceData,
    Real a,
    deviceSU3 *pU)
{
    intokernaldir;

    for (BYTE idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        //deviceSU3 expP = pMyDeviceData[linkIndex].ExpReal(a, _DC_ExpPrecision);
        deviceSU3 expP = pMyDeviceData[linkIndex].QuickExp(a);
        expP.Mul(pU[linkIndex]);
        pU[linkIndex] = expP;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelExpMultSU3Real(
    const deviceSU3 * __restrict__ pMyDeviceData,
    Real a,
    deviceSU3 *pU,
    BYTE prec)
{
    intokernaldir;

    for (BYTE idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        //deviceSU3 expP = pMyDeviceData[linkIndex].ExpReal(a, _DC_ExpPrecision);
        deviceSU3 expP = pMyDeviceData[linkIndex].ExpReal(a, prec);
        expP.Mul(pU[linkIndex]);
        pU[linkIndex] = expP;
    }
}

/**
* Trace (P^2)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateKinematicEnergySU3(const deviceSU3 * __restrict__ pDeviceData, 
#if !_CLG_DOUBLEFLOAT
    DOUBLE* results
#else
    Real* results
#endif
)
{
    intokernaldir;

#if !_CLG_DOUBLEFLOAT
    DOUBLE resThisThread = 0.0;
#else
    Real resThisThread = F(0.0);
#endif
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        resThisThread += pDeviceData[linkIndex].DaggerMulC(pDeviceData[linkIndex]).ReTr();
    }

    results[uiSiteIndex] = resThisThread;
}


__global__ void _CLG_LAUNCH_BOUND
_kernelNormalizeSU3(deviceSU3 * pMyDeviceData)
{
    intokernaldir;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        pMyDeviceData[linkIndex].Norm();
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDotSU3(
    const deviceSU3 * __restrict__ pMyDeviceData, 
    const deviceSU3 * __restrict__ pOtherDeviceData,
#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex* result
#else
    CLGComplex* result
#endif
)
{
    intokernaldir;

#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex resThisThread = make_cuDoubleComplex(0, 0);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        resThisThread = cuCadd(resThisThread, 
            _cToDouble(pMyDeviceData[linkIndex].DaggerMulC(pOtherDeviceData[linkIndex]).Tr())
        );
    }

    result[uiSiteIndex] = resThisThread;
#else
    CLGComplex resThisThread = _make_cuComplex(0,0);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        resThisThread = _cuCaddf(resThisThread, pMyDeviceData[linkIndex].DaggerMulC(pOtherDeviceData[linkIndex]).Tr());
    }

    result[uiSiteIndex] = resThisThread;
#endif
    //printf("res = %f %f\n", pOtherDeviceData[uiSiteIndex * 4].m_me[0].x, pMyDeviceData[uiSiteIndex * 4].m_me[0].x);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSetConfigurationSU3(
    deviceSU3* pDeviceData,
    const Real* __restrict__ pRealData)
{
    gaugeSU3KernelFuncionStart

        //In Bridge, it is t,z,y,x
        //x + y * nx + z * nx * ny + t * nx * ny * nz
    SSmallInt4 xyzt = __deviceSiteIndexToInt4(uiSiteIndex);
    UINT uiBridgeSiteIndex = xyzt.w * _DC_Lx * _DC_Ly * _DC_Lz + xyzt.z * _DC_Lx * _DC_Ly + xyzt.y * _DC_Lx + xyzt.x;
    UINT uiBridgeLinkIndex = uiBridgeSiteIndex * _DC_Dir + idir;

    //0 1 2
    //3 4 5
    //6 7 8
    pDeviceData[uiLinkIndex].m_me[0] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex +  0], pRealData[18 * uiBridgeLinkIndex +  1]);
    pDeviceData[uiLinkIndex].m_me[1] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex +  2], pRealData[18 * uiBridgeLinkIndex +  3]);
    pDeviceData[uiLinkIndex].m_me[2] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex +  4], pRealData[18 * uiBridgeLinkIndex +  5]);
    pDeviceData[uiLinkIndex].m_me[3] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex +  6], pRealData[18 * uiBridgeLinkIndex +  7]);
    pDeviceData[uiLinkIndex].m_me[4] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex +  8], pRealData[18 * uiBridgeLinkIndex +  9]);
    pDeviceData[uiLinkIndex].m_me[5] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex + 10], pRealData[18 * uiBridgeLinkIndex + 11]);
    pDeviceData[uiLinkIndex].m_me[6] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex + 12], pRealData[18 * uiBridgeLinkIndex + 13]);
    pDeviceData[uiLinkIndex].m_me[7] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex + 14], pRealData[18 * uiBridgeLinkIndex + 15]);
    pDeviceData[uiLinkIndex].m_me[8] = _make_cuComplex(pRealData[18 * uiBridgeLinkIndex + 16], pRealData[18 * uiBridgeLinkIndex + 17]);

    //pDeviceData[uiLinkIndex].DebugPrint();
    gaugeSU3KernelFuncionEnd
}

/**
 * iA = U.TA() 
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToIA(
    deviceSU3* pDeviceData)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        pDeviceData[uiLinkIndex].Ta();
    }
}

#if _CLG_DEBUG
__global__ void _CLG_LAUNCH_BOUND_HALF
#else
__global__ void _CLG_LAUNCH_BOUND
#endif
_kernelTransformToIALog(
    deviceSU3* pDeviceData)
{
    intokernaldir;

    //printf("I should have be in here\n");
    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        //if (52 == uiLinkIndex)
        //{
        //    pDeviceData[uiLinkIndex].DebugPrint("a");
        //}
        pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].Log();
        //if (52 == uiLinkIndex)
        //{
        //    pDeviceData[uiLinkIndex].DebugPrint("b");
        //}
    }
}

/**
 * E_mu = F_{0 mu}
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToE(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pRes)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    
    //we only need uiDir = 0,1,2
    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        deviceSU3 res = deviceSU3::makeSU3Zero();
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        if (dir < uiDir - 1)
        {
            //find clover F
            //res = _deviceClover(pDeviceData, sSite4, uiBigIdx, 3, dir, byFieldId);
            //test not using clover
            res = _device1PlaqutteTermPP(pDeviceData, 3, dir, uiBigIdx, sSite4, byFieldId);
            //res.iIm2();
            //not using clover not multiply 0.25
            //res.MulReal(F(0.125));

            //test clover
            //INT dirs[4];
            //dirs[0] = 4;
            //dirs[1] = dir + 1;
            //dirs[2] = -4;
            //dirs[3] = -static_cast<INT>(dir) - 1;
            //SSmallInt4 sSite4_m_mu = sSite4;
            //sSite4_m_mu.m_byData4[dir] = sSite4_m_mu.m_byData4[dir] - 1;
            //SSmallInt4 sSite4_m_mu_m_t = sSite4_m_mu;
            //sSite4_m_mu_m_t.m_byData4[3] = sSite4_m_mu_m_t.m_byData4[3] - 1;
            //SSmallInt4 sSite4_m_t = sSite4;
            //sSite4_m_t.m_byData4[3] = sSite4_m_t.m_byData4[3] - 1;

            //res = _deviceLink(pDeviceData, sSite4, 4, byFieldId, dirs);
            //res.Add(_deviceLink(pDeviceData, sSite4_m_mu, 4, byFieldId, dirs));
            //res.Add(_deviceLink(pDeviceData, sSite4_m_mu_m_t, 4, byFieldId, dirs));
            //res.Add(_deviceLink(pDeviceData, sSite4_m_t, 4, byFieldId, dirs));
            //res.MulReal(0.25);

            res.Ta();
            res.MulReal(F(-1.0));
            //res.MulReal(F(0.25));
        }

        pRes[uiLinkIndex] = res;
    }
}

/**
 * This is wrong! the order of the plaqutte must be considered
 * This is to make sure gauge transform is g(x) nabla E g^+(n)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateNablaE(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byFieldId, deviceSU3* pRes)
{
    intokernalInt4;
    //const BYTE uiDir2 = static_cast<BYTE>(_DC_Dir) * 2;
    //const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);


    //i=0: 12
    //  1: 13
    //  2: 14
    //  3: 23
    //  4: 24
    //  5: 34  

    UINT uiResLinkIdx = _deviceGetLinkIndex(uiSiteIndex, 3);

    pRes[uiResLinkIdx] = deviceSU3::makeSU3Zero();
    #pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {
        //we need 2, 4 and 5
        //BYTE byPlaqIdx = (dir + 1) << 1;
        //if (byPlaqIdx > 5) byPlaqIdx = 5;

        INT dirs[4];
        //deviceSU3 toMul(_devicePlaqutte(pDeviceData, pCachedPlaqutte, uiSiteIndex, byPlaqIdx, plaqLength, plaqCount));
        dirs[0] = 4;
        dirs[1] = dir + 1;
        dirs[2] = -4;
        dirs[3] = -static_cast<INT>(dir) - 1;
        deviceSU3 toMul(
            //_device1PlaqutteTermPP(pDeviceData, 3, dir, uiBigIdx)
            _deviceLink(pDeviceData, sSite4, 4, byFieldId, dirs)
        );

        //
        dirs[0] = 4;
        dirs[1] = -static_cast<INT>(dir) - 1;
        dirs[2] = -4;
        dirs[3] = dir + 1;

        toMul.Mul(
            _deviceLink(pDeviceData, sSite4, 4, byFieldId, dirs)
        );
        //toMul.Ta();
        pRes[uiResLinkIdx].Sub(toMul);
    }
    pRes[uiResLinkIdx].Ta();
    //pRes[uiResLinkIdx].SubReal(F(3.0));
}

/**
 * Larger than the above
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateNablaENaive(
    const deviceSU3* __restrict__ pDeviceData,
    BYTE byFieldId, deviceSU3* pRes)
{
    intokernalInt4;
    //const BYTE uiDir2 = static_cast<BYTE>(_DC_Dir) * 2;
    //const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);


    //i=0: 12
    //  1: 13
    //  2: 14
    //  3: 23
    //  4: 24
    //  5: 34  

    UINT uiResLinkIdx = _deviceGetLinkIndex(uiSiteIndex, 3);

    pRes[uiResLinkIdx] = deviceSU3::makeSU3Zero();
    #pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {
        //we need 2, 4 and 5
        //BYTE byPlaqIdx = (dir + 1) << 1;
        //if (byPlaqIdx > 5) byPlaqIdx = 5;

        INT dirs[4];
        //deviceSU3 toMul(_devicePlaqutte(pDeviceData, pCachedPlaqutte, uiSiteIndex, byPlaqIdx, plaqLength, plaqCount));
        dirs[0] = 4;
        dirs[1] = dir + 1;
        dirs[2] = -4;
        dirs[3] = -static_cast<INT>(dir) - 1;
        //SSmallInt4 sSite4_m_mu = sSite4;
        //sSite4_m_mu.m_byData4[dir] = sSite4_m_mu.m_byData4[dir] - 1;
        //SSmallInt4 sSite4_m_mu_m_t = sSite4_m_mu;
        //sSite4_m_mu_m_t.m_byData4[3] = sSite4_m_mu_m_t.m_byData4[3] - 1;
        //SSmallInt4 sSite4_m_t = sSite4;
        //sSite4_m_t.m_byData4[3] = sSite4_m_t.m_byData4[3] - 1;

        deviceSU3 a(
            //_device1PlaqutteTermPP(pDeviceData, 3, dir, uiBigIdx)
            _deviceLink(pDeviceData, sSite4, 4, byFieldId, dirs)
            //_deviceClover(pDeviceData, sSite4, __bi(sSite4), 3, dir, byFieldId)
        );
        //a.Add(_deviceLink(pDeviceData, sSite4_m_mu, 4, byFieldId, dirs));
        //a.Add(_deviceLink(pDeviceData, sSite4_m_mu_m_t, 4, byFieldId, dirs));
        //a.Add(_deviceLink(pDeviceData, sSite4_m_t, 4, byFieldId, dirs));
        //a.Ta();
        
        //SSmallInt4 sSite4_m_2mu = sSite4_m_mu;
        //sSite4_m_2mu.m_byData4[dir] = sSite4_m_2mu.m_byData4[dir] - 1;
        //SSmallInt4 sSite4_m_2mu_m_t = sSite4_m_2mu;
        //sSite4_m_2mu_m_t.m_byData4[3] = sSite4_m_2mu_m_t.m_byData4[3] - 1;
        dirs[0] = -static_cast<INT>(dir) - 1;
        dirs[1] = 4;
        dirs[2] = dir + 1;
        dirs[3] = -4;

        deviceSU3 b(
            //_device1PlaqutteTermPP(pDeviceData, 3, dir, uiBigIdx)
            _deviceLink(pDeviceData, sSite4, 4, byFieldId, dirs)
            //_deviceClover(pDeviceData, sSite4_m_mu, __bi(sSite4_m_mu), 3, dir, byFieldId)
        );
        //b.Add(_deviceLink(pDeviceData, sSite4_m_mu_m_t, 4, byFieldId, dirs));
        //b.Add(_deviceLink(pDeviceData, sSite4_m_2mu, 4, byFieldId, dirs));
        //b.Add(_deviceLink(pDeviceData, sSite4_m_2mu_m_t, 4, byFieldId, dirs));

        //b.Ta();
        a.Sub(b);
        //a.Ta();
        //a.MulReal(F(0.25));
        pRes[uiResLinkIdx].Sub(a);
    }
}

/**
 * U = exp(A)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToU(
    deviceSU3* pDeviceData)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        pDeviceData[uiLinkIndex] = (0 == _DC_ExpPrecision)
            ? pDeviceData[uiLinkIndex].QuickExp(F(1.0))
            : pDeviceData[uiLinkIndex].ExpReal(F(1.0), _DC_ExpPrecision);
    }
}

#if _CLG_DEBUG
__global__ void _CLG_LAUNCH_BOUND_HALF
#else
__global__ void _CLG_LAUNCH_BOUND
#endif
_kernelTransformToULog(
    deviceSU3* pDeviceData)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        //if (1023 == uiLinkIndex)
        //{
        //    pDeviceData[uiLinkIndex].DebugPrint("a");
        //}
        //pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].StrictExp(1023 == uiLinkIndex);
        //if (1023 == uiLinkIndex)
        //{
        //    pDeviceData[uiLinkIndex].DebugPrint("b");
        //}
        pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].StrictExp();
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSetOneDirUnity(deviceSU3* pDeviceData, BYTE byDir)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = deviceSU3::makeSU3Id();
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSetOneDirZero(deviceSU3* pDeviceData, BYTE byDir)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = deviceSU3::makeSU3Zero();
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelSetOneDirUnityPoint(deviceSU3* pDeviceData, UINT uiSiteIndex, BYTE byDir)
{
    for (BYTE dir = 0; dir < 4; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = deviceSU3::makeSU3Id();
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelSetOneDirZeroPoint(deviceSU3* pDeviceData, UINT uiSiteIndex, BYTE byDir)
{
    for (BYTE dir = 0; dir < 4; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = deviceSU3::makeSU3Zero();
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovLoopOfSiteSU3(
    const deviceSU3* __restrict__ pDeviceBuffer,
    cuDoubleComplex* res,
    BYTE byFieldId)
{
    UINT uiXYZ = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    const UINT uiSiteIndex = uiXYZ * _DC_Lt;
    UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, _DC_Dir - 1);
    SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

    deviceSU3 tmp = deviceSU3::makeSU3Zero();
    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, _DC_Dir - 1))
    {
        tmp = pDeviceBuffer[uiLinkIdx];
    }

    for (UINT uiT = 1; uiT < _DC_Lt; ++uiT)
    {
        UINT newSiteIndex = uiSiteIndex + uiT;
        uiLinkIdx = _deviceGetLinkIndex(newSiteIndex, _DC_Dir - 1);
        site4 = __deviceSiteIndexToInt4(newSiteIndex);
        uiBigIdx = __idx->_deviceGetBigIndex(site4);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, _DC_Dir - 1))
        {
            tmp.Mul(pDeviceBuffer[uiLinkIdx]);
        }
    }

    res[uiXYZ] = _cToDouble(tmp.Tr());
}

#pragma endregion

void CFieldGaugeSU3::AxpyPlus(const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    preparethread;
    _kernelAxpyPlusSU3 << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData);
}

void CFieldGaugeSU3::AxpyMinus(const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    preparethread;
    _kernelAxpyMinusSU3 << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData);

}

void CFieldGaugeSU3::ScalarMultply(const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplySU3Complex << <block, threads >> > (m_pDeviceData, a);
}

void CFieldGaugeSU3::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplySU3Real << <block, threads >> > (m_pDeviceData, a);
}

void CFieldGaugeSU3::Axpy(Real a, const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    preparethread;
    _kernelAxpySU3Real << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData, a);
}

void CFieldGaugeSU3::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || EFT_GaugeSU3 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(x);
    preparethread;
    _kernelAxpySU3A << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData, a);
}

void CFieldGaugeSU3::Mul(const CField* other, UBOOL bDagger)
{
    if (NULL == other || EFT_GaugeSU3 != other->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeSU3* pSU3x = dynamic_cast<const CFieldGaugeSU3*>(other);
    preparethread;
    _kernelMulSU3 << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData, bDagger);
}


void CFieldGaugeSU3::Zero()
{
    preparethread;
    _kernelInitialSU3Feield << <block, threads >> > (m_pDeviceData, EFIT_Zero);
}

void CFieldGaugeSU3::Identity()
{
    preparethread;
    _kernelInitialSU3Feield << <block, threads >> > (m_pDeviceData, EFIT_Identity);
}

void CFieldGaugeSU3::Dagger()
{
    preparethread;
    _kernelDaggerSU3 << <block, threads >> > (m_pDeviceData);
}

void CFieldGaugeSU3::MakeRandomGenerator()
{
    preparethread;
    _kernelInitialSU3Feield << <block, threads >> > (m_pDeviceData, EFIT_RandomGenerator);
}

/**
*
*/
void CFieldGaugeSU3::InitialField(EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialSU3Feield << <block, threads >> > (m_pDeviceData, eInitialType);
}

void CFieldGaugeSU3::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eType)
{
    if (!CFileSystem::IsFileExist(sFileName))
    {
        appCrucial(_T("File not exist!!! %s \n"), sFileName.c_str());
        _FAIL_EXIT;
    }

    switch (eType)
    {
    case EFFT_BridgePPTXT:
    {
        const CCString sContent = appGetFileSystem()->ReadAllText(sFileName);
        TArray<INT> seps;
        seps.AddItem(_T('\n'));
        seps.AddItem(_T('\r'));
        TArray<CCString> sStringlist = appGetStringList(sContent, seps, 0x7fffffff);
        assert(static_cast<UINT>(sStringlist.Num()) == _HC_LinkCount * 18);

        Real* pData = (Real*)malloc(sizeof(Real) * sStringlist.Num());
        for (INT i = 0; i < sStringlist.Num(); ++i)
        {
            pData[i] = static_cast<Real>(appStrToDOUBLE(sStringlist[i]));
        }

        SetByArray(pData);
    }
    break;
    case EFFT_BridgePPBin:
    {
        UINT uiSize = 0;
        BYTE* allBytes = appGetFileSystem()->ReadAllBytes(sFileName, uiSize);
        assert(uiSize == 8 * 18 * _HC_LinkCount);
        Real* pData = (Real*)malloc(sizeof(Real) * 18 * _HC_LinkCount);
        for (UINT i = 0; i < 18 * _HC_LinkCount; ++i)
        {
            BYTE data[8];
            for (INT j = 0; j < 8; ++j)
            {
                data[j] = allBytes[i * 8 + (7 - j)];
            }
            DOUBLE * dbData = (DOUBLE*)data;
            pData[i] = static_cast<Real>(*dbData);
        }
        free(allBytes);

        SetByArray(pData);
    }
    break;
    case EFFT_CLGBin:
#if _CLG_DOUBLEFLOAT
    case EFFT_CLGBinDouble:
#else
    case EFFT_CLGBinFloat:
#endif
    {
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 18 * m_uiLinkeCount);
        BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (uiSize != sizeof(Real) * 18 * _HC_LinkCount)
        {
            appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(Real) * 18 * _HC_LinkCount));
        }
        InitialWithByte(data);
        free(data);
    }
    break;
#if _CLG_DOUBLEFLOAT
    case EFFT_CLGBinFloat:
    {
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 18 * m_uiLinkeCount);
        BYTE* data = (BYTE*)malloc(uiSize);
        Real* rdata = (Real*)data;
        FLOAT* fdata = (FLOAT*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (uiSize != sizeof(FLOAT) * 18 * _HC_LinkCount)
        {
            appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(FLOAT) * 18 * _HC_LinkCount));
        }
        for (UINT i = 0; i < 18 * m_uiLinkeCount; ++i)
        {
            rdata[i] = static_cast<Real>(fdata[i]);
        }
        InitialWithByte(data);
        free(fdata);
        free(data);
    }
    break;
#else
    case EFFT_CLGBinDouble:
    {
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 18 * m_uiLinkeCount);
        BYTE* data = (BYTE*)malloc(uiSize);
        Real* rdata = (Real*)data;
        DOUBLE* ddata = (DOUBLE*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (uiSize != sizeof(DOUBLE) * 18 * _HC_LinkCount)
        {
            appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(DOUBLE) * 18 * _HC_LinkCount));
        }
        for (UINT i = 0; i < 18 * m_uiLinkeCount; ++i)
        {
            rdata[i] = static_cast<Real>(ddata[i]);
        }
        InitialWithByte(data);
        free(ddata);
        free(data);
    }
    break;
#endif
    case EFFT_CLGBinCompressed:
    {
        InitialWithByteCompressed(sFileName);
    }
    break;
    default:
        appCrucial(_T("Not supported input file type %s\n"), __ENUM_TO_STRING(EFieldFileType, eType).c_str());
        break;

    }
}

void CFieldGaugeSU3::InitialWithByte(BYTE* byData)
{
    deviceSU3* readData = (deviceSU3*)malloc(sizeof(deviceSU3) * m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[18];
        memcpy(oneLink, byData + sizeof(Real) * 18 * i, sizeof(Real) * 18);
        for (UINT j = 0; j < 16; ++j)
        {
            if (j < 9)
            {
                readData[i].m_me[j] =
                    _make_cuComplex(
                        oneLink[2 * j],
                        oneLink[2 * j + 1]);
            }
            else
            {
                readData[i].m_me[j] = _make_cuComplex(F(0.0), F(0.0));
            }
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceSU3) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    free(readData);
}

void CFieldGaugeSU3::InitialWithByteCompressed(const CCString& sFileName)
{
    UINT uiSize = static_cast<UINT>(sizeof(Real) * 9 * m_uiLinkeCount);
    BYTE* byData = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);

    deviceSU3* readData = (deviceSU3*)malloc(sizeof(deviceSU3) * m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[9];
        memcpy(oneLink, byData + sizeof(Real) * 9 * i, sizeof(Real) * 9);

        readData[i].m_me[1] = _make_cuComplex(oneLink[0], oneLink[1]);
        readData[i].m_me[3] = _make_cuComplex(-oneLink[0], oneLink[1]);

        readData[i].m_me[2] = _make_cuComplex(oneLink[2], oneLink[3]);
        readData[i].m_me[6] = _make_cuComplex(-oneLink[2], oneLink[3]);

        readData[i].m_me[5] = _make_cuComplex(oneLink[4], oneLink[5]);
        readData[i].m_me[7] = _make_cuComplex(-oneLink[4], oneLink[5]);

        readData[i].m_me[0] = _make_cuComplex(F(0.0), oneLink[6]);
        readData[i].m_me[4] = _make_cuComplex(F(0.0), oneLink[7]);
        readData[i].m_me[8] = _make_cuComplex(F(0.0), oneLink[8]);

        for (UINT j = 9; j < 16; ++j)
        {
            readData[i].m_me[j] = _zeroc;
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceSU3) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    free(readData);

    //DebugPrintMe();

    preparethread;
    _kernelTransformToULog << <block, threads >> > (m_pDeviceData);
    checkCudaErrors(cudaDeviceSynchronize());

    //DebugPrintMe();
    free(byData);
}

void CFieldGaugeSU3::SetByArray(Real* array)
{
    assert(NULL != array);
    //we algin the su3 now
    //assert(sizeof(deviceSU3) == 32 * sizeof(Real));

    //checkCudaErrors(cudaMemcpy(m_pDeviceData, array, sizeof(Real) * _HC_LinkCount * 18, cudaMemcpyHostToDevice));
    
    Real* pDeviceArray;
    checkCudaErrors(__cudaMalloc((void**)&pDeviceArray, sizeof(Real) * _HC_LinkCount * 18));
    checkCudaErrors(cudaMemcpy(pDeviceArray, array, sizeof(Real) * _HC_LinkCount * 18, cudaMemcpyHostToDevice));
    preparethread;
    _kernelSetConfigurationSU3 << <block, threads >> > (m_pDeviceData, pDeviceArray);
    checkCudaErrors(__cudaFree(pDeviceArray));

    free(array);

    ElementNormalize();
}

/**
* (1) calculate staples
* (2) calculate force(additive)
* (3) calculate energy
*/
void CFieldGaugeSU3::CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, Real betaOverN) const
{
    if (NULL == pForce || EFT_GaugeSU3 != pForce->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: force field is not SU3");
        return;
    }
    if (NULL != pStable && EFT_GaugeSU3 != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stape field is not SU3");
        return;
    }

    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    CFieldGaugeSU3* pStableSU3 = NULL == pStable ? NULL : dynamic_cast<CFieldGaugeSU3*>(pStable);

    preparethread;

    assert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache);

    _kernelStapleAtSiteSU3CacheIndex << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData,
        betaOverN);
}

void CFieldGaugeSU3::CalculateOnlyStaple(CFieldGauge* pStable) const
{
    if (NULL == pStable || EFT_GaugeSU3 != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stable field is not SU3");
        return;
    }
    CFieldGaugeSU3* pStableSU3 = dynamic_cast<CFieldGaugeSU3*>(pStable);

    preparethread;
    _kernelCalculateOnlyStaple << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStableSU3->m_pDeviceData);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeSU3::CalculatePlaqutteEnergy(DOUBLE betaOverN) const
#else
Real CFieldGaugeSU3::CalculatePlaqutteEnergy(Real betaOverN) const
#endif
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);

    preparethread;
    _kernelPlaqutteEnergySU3CacheIndex << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer
        );

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeSU3::CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const
#else
Real CFieldGaugeSU3::CalculatePlaqutteEnergyUseClover(Real betaOverN) const
#endif
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);
    //appGeneral(_T("const %f\n"), 3.0 * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite);
    preparethread;
    _kernelPlaqutteEnergySU3_UseClover << <block, threads >> > (
        m_byFieldId,
        m_pDeviceData,
#if !_CLG_DOUBLEFLOAT
        3.0 * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
#else
        F(3.0) * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
#endif
        betaOverN,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeSU3::CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge *pStable) const
#else
Real CFieldGaugeSU3::CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge* pStable) const
#endif
{
    if (NULL == pStable || EFT_GaugeSU3 != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stape field is not SU3");
        return F(0.0);
    }
    const CFieldGaugeSU3* pStableSU3 = dynamic_cast<const CFieldGaugeSU3*>(pStable);

    //appGeneral(_T("const = %f\n"), 3.0 * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink);

    preparethread;
    _kernelPlaqutteEnergyUsingStableSU3 << <block, threads >> > (
        m_pDeviceData, 
        pStableSU3->m_pDeviceData,
#if !_CLG_DOUBLEFLOAT
        3.0 * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
#else
        F(3.0) * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
#endif
        betaOverN, 
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeSU3::CalculateKinematicEnergy() const
#else
Real CFieldGaugeSU3::CalculateKinematicEnergy() const
#endif
{
    preparethread;
    _kernelCalculateKinematicEnergySU3 << <block, threads >> > (m_pDeviceData, _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

void CFieldGaugeSU3::SetOneDirectionUnity(BYTE byDir)
{
    if (0 == (byDir & 15))
    {
        return;
    }
    preparethread;
    _kernelSetOneDirUnity << <block, threads >> >(m_pDeviceData, byDir);

    //for (SBYTE byz = 0; byz < _HC_Lz; ++byz)
    //{
    //    for (SBYTE byw = 0; byw < _HC_Lt; ++byw)
    //    {
    //        _kernelSetOneDirUnityPoint <<<1,1>>>(m_pDeviceData, _hostGetSiteIndex(SSmallInt4(0, 0, byz, byw)), 15);
    //        _kernelSetOneDirUnityPoint << <1, 1 >> > (m_pDeviceData, _hostGetSiteIndex(SSmallInt4(0, 3, byz, byw)), 15);
    //        _kernelSetOneDirUnityPoint << <1, 1 >> > (m_pDeviceData, _hostGetSiteIndex(SSmallInt4(3, 0, byz, byw)), 15);
    //        _kernelSetOneDirUnityPoint << <1, 1 >> > (m_pDeviceData, _hostGetSiteIndex(SSmallInt4(3, 3, byz, byw)), 15);
    //    }
    //}
}

void CFieldGaugeSU3::SetOneDirectionZero(BYTE byDir)
{
    if (0 == (byDir & 15))
    {
        return;
    }
    preparethread;
    _kernelSetOneDirZero << <block, threads >> > (m_pDeviceData, byDir);
    //for (SBYTE byz = 0; byz < _HC_Lz; ++byz)
    //{
    //    for (SBYTE byw = 0; byw < _HC_Lt; ++byw)
    //    {
    //        _kernelSetOneDirZeroPoint << <1, 1 >> > (m_pDeviceData, _hostGetSiteIndex(SSmallInt4(0, 0, byz, byw)), 15);
    //        _kernelSetOneDirZeroPoint << <1, 1 >> > (m_pDeviceData, _hostGetSiteIndex(SSmallInt4(0, 3, byz, byw)), 15);
    //        _kernelSetOneDirZeroPoint << <1, 1 >> > (m_pDeviceData, _hostGetSiteIndex(SSmallInt4(3, 0, byz, byw)), 15);
    //        _kernelSetOneDirZeroPoint << <1, 1 >> > (m_pDeviceData, _hostGetSiteIndex(SSmallInt4(3, 3, byz, byw)), 15);
    //    }
    //}
}

CFieldGaugeSU3::CFieldGaugeSU3() : CFieldGauge()
{
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount));
}

CFieldGaugeSU3::~CFieldGaugeSU3()
{
    checkCudaErrors(__cudaFree(m_pDeviceData));
}

void CFieldGaugeSU3::ExpMult(Real a, CField* U) const
{
    if (NULL == U || EFT_GaugeSU3 != U->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: U field is not SU3");
        return;
    }

    CFieldGaugeSU3* pUField = dynamic_cast<CFieldGaugeSU3*>(U);

    preparethread;
    if (0 == _HC_ExpPrecision)
    {
        _kernelExpMultSU3RealQ << < block, threads >> > (m_pDeviceData, a, pUField->m_pDeviceData);
    }
    else
    {
        _kernelExpMultSU3Real << < block, threads >> > (m_pDeviceData, a, pUField->m_pDeviceData, static_cast<BYTE>(_HC_ExpPrecision));
    }
    
}

void CFieldGaugeSU3::ElementNormalize()
{
    preparethread;
    _kernelNormalizeSU3 << < block, threads >> > (m_pDeviceData);
}

#if !_CLG_DOUBLEFLOAT
cuDoubleComplex CFieldGaugeSU3::Dot(const CField* other) const
#else
CLGComplex CFieldGaugeSU3::Dot(const CField* other) const
#endif
{
    if (NULL == other || EFT_GaugeSU3 != other->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: U field is not SU3");
        return make_cuDoubleComplex(0,0);
    }

    const CFieldGaugeSU3* pUField = dynamic_cast<const CFieldGaugeSU3*>(other);

    preparethread;
    _kernelDotSU3 << < block, threads >> > (m_pDeviceData, pUField->m_pDeviceData, _D_ComplexThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

void CFieldGaugeSU3::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || EFT_GaugeSU3 != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: target field is not SU3");
        return;
    }

    CFieldGauge::CopyTo(pTarget);

    CFieldGaugeSU3* pTargetField = dynamic_cast<CFieldGaugeSU3*>(pTarget);
    checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceData, m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount, cudaMemcpyDeviceToDevice));
}

void CFieldGaugeSU3::TransformToIA()
{
    preparethread;
    if (0 == _HC_ALog)
    {
        _kernelTransformToIA << <block, threads >> > (m_pDeviceData);
    }
    else
    {
        _kernelTransformToIALog << <block, threads >> > (m_pDeviceData);
    }
}


void CFieldGaugeSU3::TransformToU()
{
    preparethread;
    if (0 == _HC_ALog)
    {
        _kernelTransformToU << <block, threads >> > (m_pDeviceData);
    }
    else
    {
        _kernelTransformToULog << <block, threads >> > (m_pDeviceData);
    }
    
}

void CFieldGaugeSU3::CalculateE_Using_U(CFieldGauge* pResoult) const
{
    if (NULL == pResoult || EFT_GaugeSU3 != pResoult->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: U field is not SU3");
        return;
    }

    CFieldGaugeSU3* pUField = dynamic_cast<CFieldGaugeSU3*>(pResoult);

    preparethread;
    _kernelTransformToE << <block, threads >> > (pUField->m_byFieldId, m_pDeviceData, pUField->m_pDeviceData);
}

void CFieldGaugeSU3::CalculateNablaE_Using_U(CFieldGauge* pResoult, UBOOL bNaive) const
{
    if (NULL == pResoult || EFT_GaugeSU3 != pResoult->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: U field is not SU3");
        return;
    }

    CFieldGaugeSU3* pUField = dynamic_cast<CFieldGaugeSU3*>(pResoult);

    preparethread;
    if (bNaive)
    {
        _kernelCalculateNablaENaive << <block, threads >> > (
            m_pDeviceData,
            m_byFieldId,
            pUField->m_pDeviceData);
    }
    else
    {
        _kernelCalculateNablaE << <block, threads >> > (
            m_pDeviceData,
            m_byFieldId,
            pUField->m_pDeviceData);
    }
}

void CFieldGaugeSU3::DebugPrintMe() const
{
    //preparethread;
    //_kernelPrintSU3 << < block, threads >> > (m_pDeviceData);

    //===================================================
    //Since Debug Print Me is only used to debug, we do it slow but convinient
    deviceSU3* pToPrint = (deviceSU3*)malloc(sizeof(deviceSU3) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(pToPrint, m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount, cudaMemcpyDeviceToHost));

    for (UINT uiLink = 0; uiLink < m_uiLinkeCount; ++uiLink)
    {
        UINT uiSite = uiLink / _HC_Dir;
        UINT uiDir = uiLink % _HC_Dir;
        SSmallInt4 site = __hostSiteIndexToInt4(uiSite);
        appGeneral(_T(" --- %d(%d, %d, %d, %d)_%d ---\n {{%f %s %f I, %f %s %f I, %f %s %f I},\n {%f %s %f I, %f %s %f I, %f %s %f I},\n {%f %s %f I, %f %s %f I, %f %s %f I}}\n"),
            uiLink,
            static_cast<INT>(site.x),
            static_cast<INT>(site.y),
            static_cast<INT>(site.z),
            static_cast<INT>(site.w),
            uiDir,

            pToPrint[uiLink].m_me[0].x,
            pToPrint[uiLink].m_me[0].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(pToPrint[uiLink].m_me[0].y),

            pToPrint[uiLink].m_me[1].x,
            pToPrint[uiLink].m_me[1].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(pToPrint[uiLink].m_me[1].y),

            pToPrint[uiLink].m_me[2].x,
            pToPrint[uiLink].m_me[2].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(pToPrint[uiLink].m_me[2].y),

            pToPrint[uiLink].m_me[3].x,
            pToPrint[uiLink].m_me[3].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(pToPrint[uiLink].m_me[3].y),

            pToPrint[uiLink].m_me[4].x,
            pToPrint[uiLink].m_me[4].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(pToPrint[uiLink].m_me[4].y),

            pToPrint[uiLink].m_me[5].x,
            pToPrint[uiLink].m_me[5].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(pToPrint[uiLink].m_me[5].y),

            pToPrint[uiLink].m_me[6].x,
            pToPrint[uiLink].m_me[6].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(pToPrint[uiLink].m_me[6].y),

            pToPrint[uiLink].m_me[7].x,
            pToPrint[uiLink].m_me[7].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(pToPrint[uiLink].m_me[7].y),

            pToPrint[uiLink].m_me[8].x,
            pToPrint[uiLink].m_me[8].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(pToPrint[uiLink].m_me[8].y)
        );
    }

    free(pToPrint);
}

CCString CFieldGaugeSU3::SaveToCompressedFile(const CCString& fileName) const
{
    CFieldGaugeSU3* pPooledGauge = dynamic_cast<CFieldGaugeSU3*>(GetCopy());
    //pPooledGauge->DebugPrintMe();

    preparethread;
    _kernelTransformToIALog << <block, threads >> > (pPooledGauge->m_pDeviceData);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    deviceSU3* toSave = (deviceSU3*)malloc(sizeof(deviceSU3) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, pPooledGauge->m_pDeviceData, sizeof(deviceSU3)* m_uiLinkeCount, cudaMemcpyDeviceToHost));

    //This is a traceless anti-Hermitian now, so we only save part of them
    const UINT uiSize = static_cast<UINT>(sizeof(Real)* m_uiLinkeCount * 9);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[9];
        oneLink[0] = static_cast<Real>(toSave[i].m_me[1].x);
        oneLink[1] = static_cast<Real>(toSave[i].m_me[1].y);
        oneLink[2] = static_cast<Real>(toSave[i].m_me[2].x);
        oneLink[3] = static_cast<Real>(toSave[i].m_me[2].y);
        oneLink[4] = static_cast<Real>(toSave[i].m_me[5].x);
        oneLink[5] = static_cast<Real>(toSave[i].m_me[5].y);
        oneLink[6] = static_cast<Real>(toSave[i].m_me[0].y);
        oneLink[7] = static_cast<Real>(toSave[i].m_me[4].y);
        //The element is in fact can be NOT traceless!!!!, the trace can be 2 Pi or -2 Pi !!!
        oneLink[8] = static_cast<Real>(toSave[i].m_me[8].y);

        memcpy(byToSave + i * sizeof(Real) * 9, oneLink, sizeof(Real) * 9);
    }

    appGetFileSystem()->WriteAllBytes(fileName.c_str(), byToSave, uiSize);
    //pPooledGauge->DebugPrintMe();
    free(toSave);
    CCString MD5 = CLGMD5Hash(byToSave, uiSize);
    free(byToSave);
    appSafeDelete(pPooledGauge);
    return MD5;
}

BYTE* CFieldGaugeSU3::CopyDataOut(UINT &uiSize) const
{
    deviceSU3* toSave = (deviceSU3*)malloc(sizeof(deviceSU3) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount * 18);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[18];
        for (UINT j = 0; j < 9; ++j)
        {
            oneLink[2 * j] = static_cast<Real>(toSave[i].m_me[j].x);
            oneLink[2 * j + 1] = static_cast<Real>(toSave[i].m_me[j].y);
        }
        memcpy(byToSave + i * sizeof(Real) * 18, oneLink, sizeof(Real) * 18);
    }
    free(toSave);

    return byToSave;
}

BYTE* CFieldGaugeSU3::CopyDataOutFloat(UINT& uiSize) const
{
    deviceSU3* toSave = (deviceSU3*)malloc(sizeof(deviceSU3) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(FLOAT) * m_uiLinkeCount * 18);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        FLOAT oneLink[18];
        for (UINT j = 0; j < 9; ++j)
        {
            oneLink[2 * j] = static_cast<FLOAT>(toSave[i].m_me[j].x);
            oneLink[2 * j + 1] = static_cast<FLOAT>(toSave[i].m_me[j].y);
        }
        memcpy(byToSave + i * sizeof(FLOAT) * 18, oneLink, sizeof(FLOAT) * 18);
    }
    free(toSave);

    return byToSave;
}

BYTE* CFieldGaugeSU3::CopyDataOutDouble(UINT& uiSize) const
{
    deviceSU3* toSave = (deviceSU3*)malloc(sizeof(deviceSU3) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(DOUBLE) * m_uiLinkeCount * 18);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        DOUBLE oneLink[18];
        for (UINT j = 0; j < 9; ++j)
        {
            oneLink[2 * j] = static_cast<DOUBLE>(toSave[i].m_me[j].x);
            oneLink[2 * j + 1] = static_cast<DOUBLE>(toSave[i].m_me[j].y);
        }
        memcpy(byToSave + i * sizeof(DOUBLE) * 18, oneLink, sizeof(DOUBLE) * 18);
    }
    free(toSave);

    return byToSave;
}

void CFieldGaugeSU3::PolyakovOnSpatialSite(cuDoubleComplex* buffer) const
{
    dim3 block(_HC_DecompX, _HC_DecompY, 1);
    dim3 threads(_HC_DecompLx, _HC_DecompLy, 1);
    _kernelPolyakovLoopOfSiteSU3 << <block, threads >> > (m_pDeviceData, buffer, m_byFieldId);
}



__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================