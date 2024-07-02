//=============================================================================
// FILENAME : CFieldGaugeSU2.cu
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

__CLGIMPLEMENT_CLASS(CFieldGaugeSU2)

#pragma region Kernels

/**
* Initial SU3 Field with a value
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialSU2Feield(deviceSU2* pDevicePtr, EFieldInitialType eInitialType)
{
    deviceSU2 id = deviceSU2::makeSU2Id();
    deviceSU2 zero = deviceSU2::makeSU2Zero();

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
            pDevicePtr[uiLinkIndex] = deviceSU2::makeSU2Random(_deviceGetFatIndex(uiSiteIndex, idir + 1));
            //Real fArg = __cuCargf(pDevicePtr[uiLinkIndex].Tr());
            //pDevicePtr[uiLinkIndex].MulComp(_make_cuComplex(_cos(fArg), -_sin(fArg)));
            //pDevicePtr[uiLinkIndex].Norm();
            //printf("arg=%f\n", __cuCargf(pDevicePtr[uiLinkIndex].Tr()));
        }
        break;
        case EFIT_RandomGenerator:
        {
            pDevicePtr[uiLinkIndex] = deviceSU2::makeSU2RandomGenerator(_deviceGetFatIndex(uiSiteIndex, idir + 1));
        }
        break;
        case EFIT_SumGenerator:
        {
            pDevicePtr[uiLinkIndex] = deviceSU2::makeSU2SumGenerator(F(1.0));
        }
        break;
        default:
        {
            printf("SU2 Field cannot be initialized with this type!");
        }
        break;
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDaggerSU2(deviceSU2* pDevicePtr)
{
    gaugeSU3KernelFuncionStart

        pDevicePtr[uiLinkIndex].Dagger();

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpySU2A(deviceSU2* pDevicePtr, const deviceSU2* __restrict__ x, CLGComplex a)
{
    gaugeSU3KernelFuncionStart

        pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].MulCompC(a));

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMulSU2(deviceSU2* pDevicePtr, const deviceSU2* __restrict__ x, UBOOL bDagger)
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
_kernelAxpySU2Real(deviceSU2* pDevicePtr, const deviceSU2* __restrict__ x, Real a)
{
    gaugeSU3KernelFuncionStart

        pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].MulRealC(a));

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusSU2(deviceSU2* pDevicePtr, const deviceSU2* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

        pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusSU2(deviceSU2* pDevicePtr, const deviceSU2* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

        pDevicePtr[uiLinkIndex].Sub(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplySU2Complex(deviceSU2* pDevicePtr, CLGComplex a)
{
    gaugeSU3KernelFuncionStart

        pDevicePtr[uiLinkIndex].MulComp(a);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplySU2Real(deviceSU2* pDevicePtr, Real a)
{
    gaugeSU3KernelFuncionStart

        pDevicePtr[uiLinkIndex].MulReal(a);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteSU2CacheIndex(
    const deviceSU2* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSU2* pStapleData, //can be NULL
    deviceSU2* pForceData,
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
        deviceSU2 res = deviceSU2::makeSU2Zero();

        //there are 6 staples, each is sum of two plaquttes
        for (INT i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceSU2 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                deviceSU2 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);

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
        deviceSU2 force(pDeviceData[linkIndex]);
        force.MulDagger(res);
        //test_force += F(-2.0) * betaOverN * __SU3Generators[8].MulC(force).ImTr();
        force.Ta();
        force.MulReal(betaOverN);

        //force is additive
        pForceData[linkIndex].Add(force);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateOnlyStapleSU2(
    const deviceSU2* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSU2* pStapleData)
{
    intokernaldir;

    //Real test_force = F(0.0);
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU2 res = deviceSU2::makeSU2Zero();

        //there are 6 staples, each is sum of two plaquttes
        for (INT i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceSU2 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                deviceSU2 toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);

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
_kernelPlaqutteEnergySU2CacheIndex(
    const deviceSU2* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    DOUBLE betaOverN,
    DOUBLE* results
)
{
    intokernal;

    DOUBLE resThisThread = 0.0;
    const UINT indexSkip = plaqCount * plaqLength * uiSiteIndex;
    for (BYTE i = 0; i < plaqCount; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLength + indexSkip];
        deviceSU2 toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

        if (first.NeedToDagger())
        {
            toAdd.Dagger();
        }

        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedIndex[i * plaqLength + j + indexSkip];
            deviceSU2 toMul(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
            if (first.NeedToDagger())
            {
                toAdd.MulDagger(toMul);
            }
            else
            {
                toAdd.Mul(toMul);
            }
        }

        //printf("retr = %f\n", toAdd.ReTr());
        resThisThread += (2.0 - toAdd.ReTr());
    }

    results[uiSiteIndex] = resThisThread * betaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergySU2_UseClover(
    BYTE byFieldId,
    const deviceSU2* __restrict__ pDeviceData,
    DOUBLE stapleConstant,
    DOUBLE fBetaOverN,
    DOUBLE* results
)
{
    intokernalInt4;

    DOUBLE fRes = 0.0;
    for (BYTE byDir1 = 0; byDir1 < _DC_Dir; ++byDir1)
    {
        for (BYTE byDir2 = byDir1 + 1; byDir2 < _DC_Dir; ++byDir2)
        {
            fRes += _deviceCloverRetrSU2(pDeviceData, sSite4, __bi(sSite4), byDir1, byDir2, byFieldId);
        }
    }
    fRes = stapleConstant - 0.25 * fRes;
    results[uiSiteIndex] = fRes * fBetaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyUsingStableSU2(
    const deviceSU2* __restrict__ pDeviceData,
    const deviceSU2* __restrict__ pStableData,
    DOUBLE stapleConstant,
    DOUBLE betaOverN,
    DOUBLE* results
)
{
    intokernaldir;

    DOUBLE resThisThread = 0.0;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        //For each link, there are 6 staples
        resThisThread += (stapleConstant - pDeviceData[linkIndex].MulDaggerC(pStableData[linkIndex]).ReTr());
    }

    results[uiSiteIndex] = resThisThread * betaOverN * 0.25;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelExpMultSU2Real(
    const deviceSU2* __restrict__ pMyDeviceData,
    Real a,
    deviceSU2* pU,
    BYTE prec)
{
    intokernaldir;

    for (BYTE idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU2 expP = pMyDeviceData[linkIndex].QuickExp(a);
        expP.Mul(pU[linkIndex]);
        pU[linkIndex] = expP;
    }
}

/**
* Trace (P^2)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateKinematicEnergySU2(const deviceSU2* __restrict__ pDeviceData,
    DOUBLE* results
)
{
    intokernaldir;

    DOUBLE resThisThread = 0.0;
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        resThisThread += pDeviceData[linkIndex].DaggerMulC(pDeviceData[linkIndex]).ReTr();
    }
    results[uiSiteIndex] = resThisThread;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelNormalizeSU2(deviceSU2* pMyDeviceData)
{
    intokernaldir;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        pMyDeviceData[linkIndex].Norm();
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDotSU2(
    const deviceSU2* __restrict__ pMyDeviceData,
    const deviceSU2* __restrict__ pOtherDeviceData,
    cuDoubleComplex* result
)
{
    intokernaldir;

    cuDoubleComplex resThisThread = make_cuDoubleComplex(0, 0);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        resThisThread = cuCadd(resThisThread,
            _cToDouble(pMyDeviceData[linkIndex].DaggerMulC(pOtherDeviceData[linkIndex]).Tr())
        );
    }

    result[uiSiteIndex] = resThisThread;
}

/**
 * iA = U.TA()
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToIASU2(
    deviceSU2* pDeviceData)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        pDeviceData[uiLinkIndex].Ta();
    }
}

//template<INT N, INT NofE>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelTransformToIALog(
//    deviceSU2* pDeviceData)
//{
//    intokernaldir;
//
//    //printf("I should have be in here\n");
//    for (BYTE dir = 0; dir < uiDir; ++dir)
//    {
//        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
//        pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].Log();
//    }
//}

/**
 * E_mu = F_{0 mu}
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToESU2(
    BYTE byFieldId,
    const deviceSU2* __restrict__ pDeviceData,
    deviceSU2* pRes)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //we only need uiDir = 0,1,2
    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        deviceSU2 res = deviceSU2::makeSU2Zero();
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        if (dir < uiDir - 1)
        {
            res = _device1PlaqutteTermPPSU2(pDeviceData, 3, dir, uiBigIdx, sSite4, byFieldId);
            res.Ta();
            res.MulReal(F(-1.0));
        }

        pRes[uiLinkIndex] = res;
    }
}

/**
 * This is wrong! the order of the plaqutte must be considered
 * This is to make sure gauge transform is g(x) nabla E g^+(n)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateNablaESU2(
    const deviceSU2* __restrict__ pDeviceData,
    BYTE byFieldId, deviceSU2* pRes)
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

    pRes[uiResLinkIdx] = deviceSU2::makeSU2Zero();
#pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {
        //we need 2, 4 and 5
        //BYTE byPlaqIdx = (dir + 1) << 1;
        //if (byPlaqIdx > 5) byPlaqIdx = 5;

        INT dirs[4];
        //deviceSU2 toMul(_devicePlaqutte(pDeviceData, pCachedPlaqutte, uiSiteIndex, byPlaqIdx, plaqLength, plaqCount));
        dirs[0] = 4;
        dirs[1] = dir + 1;
        dirs[2] = -4;
        dirs[3] = -static_cast<INT>(dir) - 1;
        deviceSU2 toMul(
            //_device1PlaqutteTermPP(pDeviceData, 3, dir, uiBigIdx)
            _deviceLinkSU2(pDeviceData, sSite4, 4, byFieldId, dirs)
        );

        //
        dirs[0] = 4;
        dirs[1] = -static_cast<INT>(dir) - 1;
        dirs[2] = -4;
        dirs[3] = dir + 1;

        toMul.Mul(
            _deviceLinkSU2(pDeviceData, sSite4, 4, byFieldId, dirs)
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
_kernelCalculateNablaENaiveSU2(
    const deviceSU2* __restrict__ pDeviceData,
    BYTE byFieldId, deviceSU2* pRes)
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

    pRes[uiResLinkIdx] = deviceSU2::makeSU2Zero();
    #pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {

        INT dirs[4];
        dirs[0] = 4;
        dirs[1] = dir + 1;
        dirs[2] = -4;
        dirs[3] = -static_cast<INT>(dir) - 1;


        deviceSU2 a(
            _deviceLinkSU2(pDeviceData, sSite4, 4, byFieldId, dirs)
        );

        dirs[0] = -static_cast<INT>(dir) - 1;
        dirs[1] = 4;
        dirs[2] = dir + 1;
        dirs[3] = -4;

        deviceSU2 b(
            _deviceLinkSU2(pDeviceData, sSite4, 4, byFieldId, dirs)
        );

        a.Sub(b);
        pRes[uiResLinkIdx].Sub(a);
    }
}

/**
 * U = exp(A)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToUSU2(deviceSU2* pDeviceData)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].QuickExp(F(1.0));
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToULogSU2(deviceSU2* pDeviceData)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].QuickExp(F(1.0));
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSetOneDirUnitySU2(deviceSU2* pDeviceData, BYTE byDir)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = deviceSU2::makeSU2Id();
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSetOneDirZeroSU2(deviceSU2* pDeviceData, BYTE byDir)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = deviceSU2::makeSU2Zero();
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelSetOneDirUnityPointSU2(deviceSU2* pDeviceData, UINT uiSiteIndex, BYTE byDir)
{
    for (BYTE dir = 0; dir < 4; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = deviceSU2::makeSU2Id();
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelSetOneDirZeroPointSU2(deviceSU2* pDeviceData, UINT uiSiteIndex, BYTE byDir)
{
    for (BYTE dir = 0; dir < 4; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = deviceSU2::makeSU2Zero();
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovLoopOfSiteSU2(
    const deviceSU2* __restrict__ pDeviceBuffer,
    cuDoubleComplex* res)
{
    UINT uiXYZ = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    const UINT uiSiteIndex = uiXYZ * _DC_Lt;
    UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, _DC_Dir - 1);
    SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

    deviceSU2 tmp = deviceSU2::makeSU2Zero();
    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, _DC_Dir - 1))
    {
        tmp = pDeviceBuffer[uiLinkIdx];
    }

    for (UINT uiT = 1; uiT < _DC_Lt; ++uiT)
    {
        UINT newSiteIndex = uiSiteIndex + uiT;
        uiLinkIdx = _deviceGetLinkIndex(newSiteIndex, _DC_Dir - 1);
        site4 = __deviceSiteIndexToInt4(newSiteIndex);
        uiBigIdx = __idx->_deviceGetBigIndex(site4);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, _DC_Dir - 1))
        {
            tmp.Mul(pDeviceBuffer[uiLinkIdx]);
        }
    }

    res[uiXYZ] = _cToDouble(tmp.Tr());
}

#pragma endregion

void CFieldGaugeSU2::AxpyPlus(const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU2: axpy failed because the otherfield is not SU2 with same N\n");
        return;
    }
    const CFieldGaugeSU2* pSU2x = dynamic_cast<const CFieldGaugeSU2*>(x);
    if (NULL == pSU2x)
    {
        appCrucial("CFieldGaugeSU2: axpy failed because the otherfield is not SU2 with same N\n");
        return;
    }

    preparethread;
    _kernelAxpyPlusSU2 << <block, threads >> > (m_pDeviceData, pSU2x->m_pDeviceData);
}

void CFieldGaugeSU2::AxpyMinus(const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU2: axpy failed because the otherfield is not SU2 with same N\n");
        return;
    }
    const CFieldGaugeSU2* pSU2x = dynamic_cast<const CFieldGaugeSU2*>(x);
    if (NULL == pSU2x)
    {
        appCrucial("CFieldGaugeSU2: axpy failed because the otherfield is not SU2 with same N\n");
        return;
    }

    preparethread;
    _kernelAxpyMinusSU2 << <block, threads >> > (m_pDeviceData, pSU2x->m_pDeviceData);
}

void CFieldGaugeSU2::ScalarMultply(const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplySU2Complex << <block, threads >> > (m_pDeviceData, a);
}

void CFieldGaugeSU2::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplySU2Real << <block, threads >> > (m_pDeviceData, a);
}

void CFieldGaugeSU2::Axpy(Real a, const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU2: axpy failed because the otherfield is not SU2 with same N\n");
        return;
    }
    const CFieldGaugeSU2* pSU2x = dynamic_cast<const CFieldGaugeSU2*>(x);
    if (NULL == pSU2x)
    {
        appCrucial("CFieldGaugeSU2: axpy failed because the otherfield is not SU2 with same N\n");
        return;
    }

    preparethread;
    _kernelAxpySU2Real << <block, threads >> > (m_pDeviceData, pSU2x->m_pDeviceData, a);

}

void CFieldGaugeSU2::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU2: axpy failed because the otherfield is not SU2 with same N\n");
        return;
    }
    const CFieldGaugeSU2* pSU2x = dynamic_cast<const CFieldGaugeSU2*>(x);
    if (NULL == pSU2x)
    {
        appCrucial("CFieldGaugeSU2: axpy failed because the otherfield is not SU2 with same N\n");
        return;
    }

    preparethread;
    _kernelAxpySU2A << <block, threads >> > (m_pDeviceData, pSU2x->m_pDeviceData, a);
}

void CFieldGaugeSU2::Mul(const CField* x, UBOOL bDagger)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU2: axpy failed because the otherfield is not SU2 with same N\n");
        return;
    }
    const CFieldGaugeSU2* pSU2x = dynamic_cast<const CFieldGaugeSU2*>(x);
    if (NULL == pSU2x)
    {
        appCrucial("CFieldGaugeSU2: axpy failed because the otherfield is not SU2 with same N\n");
        return;
    }

    preparethread;
    _kernelMulSU2 << <block, threads >> > (m_pDeviceData, pSU2x->m_pDeviceData, bDagger);
}

void CFieldGaugeSU2::Zero()
{
    preparethread;
    _kernelInitialSU2Feield << <block, threads >> > (m_pDeviceData, EFIT_Zero);
}

void CFieldGaugeSU2::Identity()
{
    preparethread;
    _kernelInitialSU2Feield << <block, threads >> > (m_pDeviceData, EFIT_Identity);
}

void CFieldGaugeSU2::Dagger()
{
    preparethread;
    _kernelDaggerSU2 << <block, threads >> > (m_pDeviceData);
}

void CFieldGaugeSU2::MakeRandomGenerator()
{
    preparethread;
    _kernelInitialSU2Feield << <block, threads >> > (m_pDeviceData, EFIT_RandomGenerator);
}

void CFieldGaugeSU2::InitialField(EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialSU2Feield << <block, threads >> > (m_pDeviceData, eInitialType);
}

DOUBLE CFieldGaugeSU2::CalculateKinematicEnergy() const
{
    preparethread;
    _kernelCalculateKinematicEnergySU2 << <block, threads >> > (m_pDeviceData, _D_RealThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

void CFieldGaugeSU2::SetOneDirectionUnity(BYTE byDir)
{
    if (0 == (byDir & 15))
    {
        return;
    }
    preparethread;
    _kernelSetOneDirUnitySU2 << <block, threads >> > (m_pDeviceData, byDir);
}

void CFieldGaugeSU2::SetOneDirectionZero(BYTE byDir)
{
    if (0 == (byDir & 15))
    {
        return;
    }
    preparethread;
    _kernelSetOneDirZeroSU2 << <block, threads >> > (m_pDeviceData, byDir);
}

CFieldGaugeSU2::CFieldGaugeSU2() : CFieldGauge()
{
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData, sizeof(deviceSU2) * m_uiLinkeCount));
}

CFieldGaugeSU2::~CFieldGaugeSU2()
{
    checkCudaErrors(__cudaFree(m_pDeviceData));
}

void CFieldGaugeSU2::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eType)
{
    if (!CFileSystem::IsFileExist(sFileName))
    {
        appCrucial(_T("File not exist!!! %s \n"), sFileName.c_str());
        _FAIL_EXIT;
    }

    switch (eType)
    {
    case EFFT_CLGBin:
#if _CLG_DOUBLEFLOAT
    case EFFT_CLGBinDouble:
#else
    case EFFT_CLGBinFloat:
#endif
    {
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 8 * m_uiLinkeCount);
        BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (uiSize != sizeof(Real) * 8 * _HC_LinkCount)
        {
            appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(Real) * 8 * _HC_LinkCount));
        }
        InitialWithByte(data);
        free(data);
    }
    break;
#if _CLG_DOUBLEFLOAT
    case EFFT_CLGBinFloat:
    {
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 8 * m_uiLinkeCount);
        BYTE* data = (BYTE*)malloc(uiSize);
        Real* rdata = (Real*)data;
        FLOAT* fdata = (FLOAT*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (uiSize != sizeof(FLOAT) * 8 * _HC_LinkCount)
        {
            appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, sizeof(FLOAT) * 8 * _HC_LinkCount);
        }
        for (UINT i = 0; i < 8 * m_uiLinkeCount; ++i)
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
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 8 * m_uiLinkeCount);
        BYTE* data = (BYTE*)malloc(uiSize);
        Real* rdata = (Real*)data;
        DOUBLE* ddata = (DOUBLE*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (uiSize != sizeof(DOUBLE) * 8 * _HC_LinkCount)
        {
            appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(DOUBLE) * 8 * _HC_LinkCount));
        }
        for (UINT i = 0; i < 8 * m_uiLinkeCount; ++i)
        {
            rdata[i] = static_cast<Real>(ddata[i]);
        }
        InitialWithByte(data);
        free(ddata);
        free(data);
    }
    break;
#endif
    default:
        appCrucial(_T("Not supported input file type %s\n"), __ENUM_TO_STRING(EFieldFileType, eType).c_str());
        break;

    }
}

void CFieldGaugeSU2::InitialWithByte(BYTE* byData)
{
    deviceSU2* readData = (deviceSU2*)malloc(sizeof(deviceSU2) * m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[8];
        memcpy(oneLink, byData + sizeof(Real) * 8 * i, sizeof(Real) * 8);
        for (UINT j = 0; j < 4; ++j)
        {
            readData[i].m_me[j] = _make_cuComplex(oneLink[2 * j], oneLink[2 * j + 1]);
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceSU2) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    free(readData);
}

/**
* (1) calculate staples
* (2) calculate force(additive)
* (3) calculate energy
*/
void CFieldGaugeSU2::CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, Real betaOverN) const
{
    if (NULL == pForce || GetFieldType() != pForce->GetFieldType())
    {
        appCrucial("CFieldGaugeSU2: force field is not SU3");
        return;
    }
    if (NULL != pStable && GetFieldType() != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU2: stape field is not SU2");
        return;
    }

    CFieldGaugeSU2* pForceSU3 = dynamic_cast<CFieldGaugeSU2*>(pForce);
    CFieldGaugeSU2* pStableSU3 = NULL == pStable ? NULL : dynamic_cast<CFieldGaugeSU2*>(pStable);

    preparethread;

    assert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache);

    _kernelStapleAtSiteSU2CacheIndex << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData,
        betaOverN);
}

void CFieldGaugeSU2::CalculateOnlyStaple(CFieldGauge* pStable) const
{
    if (NULL == pStable || GetFieldType() != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU2: stable field is not SU2");
        return;
    }
    CFieldGaugeSU2* pStableSU3 = dynamic_cast<CFieldGaugeSU2*>(pStable);

    preparethread;
    _kernelCalculateOnlyStapleSU2 << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStableSU3->m_pDeviceData);
}

DOUBLE CFieldGaugeSU2::CalculatePlaqutteEnergy(DOUBLE betaOverN) const
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);

    preparethread;
    _kernelPlaqutteEnergySU2CacheIndex << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer
        );

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

DOUBLE CFieldGaugeSU2::CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);
    preparethread;
    _kernelPlaqutteEnergySU2_UseClover << <block, threads >> > (
        m_byFieldId,
        m_pDeviceData,
        2.0 * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

DOUBLE CFieldGaugeSU2::CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge* pStable) const
{
    if (NULL == pStable || GetFieldType() != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stape field is not SU3");
        return F(0.0);
    }
    const CFieldGaugeSU2* pStableSU2 = dynamic_cast<const CFieldGaugeSU2*>(pStable);
    preparethread;
    _kernelPlaqutteEnergyUsingStableSU2 << <block, threads >> > (
        m_pDeviceData,
        pStableSU2->m_pDeviceData,
        2.0 * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        betaOverN,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

void CFieldGaugeSU2::ExpMult(Real a, CField* U) const
{
    if (NULL == U || GetFieldType() != U->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: U field is not SU3");
        return;
    }

    CFieldGaugeSU2* pUField = dynamic_cast<CFieldGaugeSU2*>(U);

    preparethread;
    _kernelExpMultSU2Real << < block, threads >> > (m_pDeviceData, a, pUField->m_pDeviceData, static_cast<BYTE>(_HC_ExpPrecision));
}

void CFieldGaugeSU2::ElementNormalize()
{
    preparethread;
    _kernelNormalizeSU2 << < block, threads >> > (m_pDeviceData);
}

cuDoubleComplex CFieldGaugeSU2::Dot(const CField* other) const
{
    if (NULL == other || GetFieldType() != other->GetFieldType())
    {
        appCrucial("CFieldGaugeSU2: U field is not SU3");
        return make_cuDoubleComplex(0, 0);
    }

    const CFieldGaugeSU2* pUField = dynamic_cast<const CFieldGaugeSU2*>(other);

    preparethread;
    _kernelDotSU2 << < block, threads >> > (m_pDeviceData, pUField->m_pDeviceData, _D_ComplexThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

void CFieldGaugeSU2::DebugPrintMe() const
{
    //preparethread;
    //_kernelPrintSU3 << < block, threads >> > (m_pDeviceData);

    //===================================================
    //Since Debug Print Me is only used to debug, we do it slow but convinient
    deviceSU2* pToPrint = (deviceSU2*)malloc(sizeof(deviceSU2) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(pToPrint, m_pDeviceData, sizeof(deviceSU2) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    appPushLogDate(FALSE);

    for (UINT uiLink = 0; uiLink < m_uiLinkeCount; ++uiLink)
    {
        UINT uiSite = uiLink / _HC_Dir;
        UINT uiDir = uiLink % _HC_Dir;
        SSmallInt4 site = __hostSiteIndexToInt4(uiSite);
        appGeneral(_T(" --- %d(%d, %d, %d, %d)_%d ---\n{{%f %s %f I, %f %s %f I},{%f %s %f I, %f %s %f I}};\n"), 
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
            appAbs(pToPrint[uiLink].m_me[3].y)
            );
    }
    appPopLogDate();
    free(pToPrint);
}

BYTE* CFieldGaugeSU2::CopyDataOut(UINT& uiSize) const
{
    deviceSU2* toSave = (deviceSU2*)malloc(sizeof(deviceSU2) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceSU2) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount * 8);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[8];
        for (UINT j = 0; j < 4; ++j)
        {
            oneLink[2 * j] = static_cast<Real>(toSave[i].m_me[j].x);
            oneLink[2 * j + 1] = static_cast<Real>(toSave[i].m_me[j].y);
        }
        memcpy(byToSave + i * sizeof(Real) * 8, oneLink, sizeof(Real) * 8);
    }
    free(toSave);

    return byToSave;
}

BYTE* CFieldGaugeSU2::CopyDataOutFloat(UINT& uiSize) const
{
    deviceSU2* toSave = (deviceSU2*)malloc(sizeof(deviceSU2) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceSU2) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(FLOAT) * m_uiLinkeCount * 8);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        FLOAT oneLink[8];
        for (UINT j = 0; j < 4; ++j)
        {
            oneLink[2 * j] = static_cast<FLOAT>(toSave[i].m_me[j].x);
            oneLink[2 * j + 1] = static_cast<FLOAT>(toSave[i].m_me[j].y);
        }
        memcpy(byToSave + i * sizeof(FLOAT) * 8, oneLink, sizeof(FLOAT) * 8);
    }
    free(toSave);

    return byToSave;
}

BYTE* CFieldGaugeSU2::CopyDataOutDouble(UINT& uiSize) const
{
    deviceSU2* toSave = (deviceSU2*)malloc(sizeof(deviceSU2) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceSU2) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(DOUBLE) * m_uiLinkeCount * 8);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        DOUBLE oneLink[8];
        for (UINT j = 0; j < 4; ++j)
        {
            oneLink[2 * j] = static_cast<DOUBLE>(toSave[i].m_me[j].x);
            oneLink[2 * j + 1] = static_cast<DOUBLE>(toSave[i].m_me[j].y);
        }
        memcpy(byToSave + i * sizeof(DOUBLE) * 8, oneLink, sizeof(DOUBLE) * 8);
    }
    free(toSave);

    return byToSave;
}

void CFieldGaugeSU2::PolyakovOnSpatialSite(cuDoubleComplex* buffer) const
{
    dim3 block(_HC_DecompX, _HC_DecompY, 1);
    dim3 threads(_HC_DecompLx, _HC_DecompLy, 1);
    _kernelPolyakovLoopOfSiteSU2 << <block, threads >> > (m_pDeviceData, buffer);
}

void CFieldGaugeSU2::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || EFT_GaugeSU2 != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: target field is not SU3");
        return;
    }

    CFieldGauge::CopyTo(pTarget);

    CFieldGaugeSU2* pTargetField = dynamic_cast<CFieldGaugeSU2*>(pTarget);
    checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceData, m_pDeviceData, sizeof(deviceSU2) * m_uiLinkeCount, cudaMemcpyDeviceToDevice));
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================