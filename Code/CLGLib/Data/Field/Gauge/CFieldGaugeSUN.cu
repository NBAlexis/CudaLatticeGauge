//=============================================================================
// FILENAME : CFieldGaugeSUN.cu
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
#include "CFieldGaugeSUN.h"

__BEGIN_NAMESPACE

//__CLGIMPLEMENT_CLASS(CFieldGaugeSU3)

#pragma region Kernels

/**
* Initial SU3 Field with a value
*/
template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialSUNFeield(deviceSUN<N, NofE> *pDevicePtr, EFieldInitialType eInitialType)
{
    deviceSUN<N, NofE> id = deviceSUN<N, NofE>::makeSUNId();
    deviceSUN<N, NofE> zero = deviceSUN<N, NofE>::makeSUNZero();

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
            pDevicePtr[uiLinkIndex] = deviceSUN<N, NofE>::makeSUNRandom(_deviceGetFatIndex(uiSiteIndex, idir + 1));
            //Real fArg = __cuCargf(pDevicePtr[uiLinkIndex].Tr());
            //pDevicePtr[uiLinkIndex].MulComp(_make_cuComplex(_cos(fArg), -_sin(fArg)));
            //pDevicePtr[uiLinkIndex].Norm();
            //printf("arg=%f\n", __cuCargf(pDevicePtr[uiLinkIndex].Tr()));
        }
        break;
        case EFIT_RandomGenerator:
        {
            pDevicePtr[uiLinkIndex] = deviceSUN<N, NofE>::makeSUNRandomGenerator(_deviceGetFatIndex(uiSiteIndex, idir + 1));
        }
        break;
        case EFIT_SumGenerator:
        {
            pDevicePtr[uiLinkIndex] = deviceSUN<N, NofE>::makeSUNSumGenerator(F(1.0));
        }
        break;
        default:
        {
            printf("SUN Field cannot be initialized with this type!");
        }
        break;
        }
    }
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelDaggerSUN(deviceSUN<N, NofE>* pDevicePtr)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Dagger();

    gaugeSU3KernelFuncionEnd
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpySUNA(deviceSUN<N, NofE> *pDevicePtr, const deviceSUN<N, NofE>* __restrict__ x, CLGComplex a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].MulCompC(a));

    gaugeSU3KernelFuncionEnd
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelMulSUN(deviceSUN<N, NofE>* pDevicePtr, const deviceSUN<N, NofE>* __restrict__ x, UBOOL bDagger)
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

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpySUNReal(deviceSUN<N, NofE> *pDevicePtr, const deviceSUN<N, NofE>* __restrict__ x, Real a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex].MulRealC(a));

    gaugeSU3KernelFuncionEnd
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusSUN(deviceSUN<N, NofE> *pDevicePtr, const deviceSUN<N, NofE>* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Add(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusSUN(deviceSUN<N, NofE> *pDevicePtr, const deviceSUN<N, NofE>* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].Sub(x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplySUNComplex(deviceSUN<N, NofE> *pDevicePtr, CLGComplex a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].MulComp(a);

    gaugeSU3KernelFuncionEnd
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplySUNReal(deviceSUN<N, NofE> *pDevicePtr, Real a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex].MulReal(a);

    gaugeSU3KernelFuncionEnd
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteSUNCacheIndex(
    const deviceSUN<N, NofE> * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSUN<N, NofE> *pStapleData, //can be NULL
    deviceSUN<N, NofE> *pForceData,
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
        deviceSUN<N, NofE> res = deviceSUN<N, NofE>::makeSUNZero();

        //there are 6 staples, each is sum of two plaquttes
        for (INT i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceSUN<N, NofE> toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                deviceSUN<N, NofE> toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);

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
        deviceSUN<N, NofE> force(pDeviceData[linkIndex]);
        force.MulDagger(res);
        //test_force += F(-2.0) * betaOverN * __SU3Generators[8].MulC(force).ImTr();
        force.Ta();
        force.MulReal(betaOverN);

        //force is additive
        pForceData[linkIndex].Add(force);
    }
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateOnlyStapleSUN(
    const deviceSUN<N, NofE> * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceSUN<N, NofE> *pStapleData)
{
    intokernaldir;

    //Real test_force = F(0.0);
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSUN<N, NofE> res = deviceSUN<N, NofE>::makeSUNZero();

        //there are 6 staples, each is sum of two plaquttes
        for (INT i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceSUN<N, NofE> toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                deviceSUN<N, NofE> toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);

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

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergySUNCacheIndex(
    const deviceSUN<N, NofE> * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
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
        deviceSUN<N, NofE> toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

        if (first.NeedToDagger())
        {
            toAdd.Dagger();
        }

        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedIndex[i * plaqLength + j + indexSkip];
            deviceSUN<N, NofE> toMul(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
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
        resThisThread += (static_cast<DOUBLE>(N) - toAdd.ReTr());
    }

    results[uiSiteIndex] = resThisThread * betaOverN;
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergySUN_UseClover(
    BYTE byFieldId,
    const deviceSUN<N, NofE>* __restrict__ pDeviceData,
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
            fRes += _deviceCloverRetrSUN(pDeviceData, sSite4, __bi(sSite4), byDir1, byDir2, byFieldId);
        }
    }
    fRes = stapleConstant - 0.25 * fRes;
    results[uiSiteIndex] = fRes * fBetaOverN;
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyUsingStableSUN(
    const deviceSUN<N, NofE> * __restrict__ pDeviceData,
    const deviceSUN<N, NofE> * __restrict__ pStableData,
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

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelExpMultSUNReal(
    const deviceSUN<N, NofE> * __restrict__ pMyDeviceData,
    Real a,
    deviceSUN<N, NofE> *pU,
    BYTE prec)
{
    intokernaldir;

    for (BYTE idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSUN<N, NofE> expP = pMyDeviceData[linkIndex].ExpReal(a, (prec <= N) ? (N + 1) : prec);
        expP.Mul(pU[linkIndex]);
        pU[linkIndex] = expP;
    }
}

/**
* Trace (P^2)
*/
template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateKinematicEnergySUN(const deviceSUN<N, NofE> * __restrict__ pDeviceData, 
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

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelNormalizeSUN(deviceSUN<N, NofE> * pMyDeviceData)
{
    intokernaldir;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        pMyDeviceData[linkIndex].Norm();
    }
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelDotSUN(
    const deviceSUN<N, NofE> * __restrict__ pMyDeviceData, 
    const deviceSUN<N, NofE> * __restrict__ pOtherDeviceData,
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
template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToIASUN(
    deviceSUN<N, NofE>* pDeviceData)
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
//    deviceSUN<N, NofE>* pDeviceData)
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
template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToESUN(
    BYTE byFieldId,
    const deviceSUN<N, NofE>* __restrict__ pDeviceData,
    deviceSUN<N, NofE>* pRes)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    
    //we only need uiDir = 0,1,2
    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        deviceSUN<N, NofE> res = deviceSUN<N, NofE>::makeSUNZero();
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        if (dir < uiDir - 1)
        {
            res = _device1PlaqutteTermPPSUN(pDeviceData, 3, dir, uiBigIdx, sSite4, byFieldId);
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
template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateNablaESUN(
    const deviceSUN<N, NofE>* __restrict__ pDeviceData,
    BYTE byFieldId, deviceSUN<N, NofE>* pRes)
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

    pRes[uiResLinkIdx] = deviceSUN<N, NofE>::makeSUNZero();
    #pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {
        //we need 2, 4 and 5
        //BYTE byPlaqIdx = (dir + 1) << 1;
        //if (byPlaqIdx > 5) byPlaqIdx = 5;

        INT dirs[4];
        //deviceSUN<N, NofE> toMul(_devicePlaqutte(pDeviceData, pCachedPlaqutte, uiSiteIndex, byPlaqIdx, plaqLength, plaqCount));
        dirs[0] = 4;
        dirs[1] = dir + 1;
        dirs[2] = -4;
        dirs[3] = -static_cast<INT>(dir) - 1;
        deviceSUN<N, NofE> toMul(
            //_device1PlaqutteTermPP(pDeviceData, 3, dir, uiBigIdx)
            _deviceLinkSUN(pDeviceData, sSite4, 4, byFieldId, dirs)
        );

        //
        dirs[0] = 4;
        dirs[1] = -static_cast<INT>(dir) - 1;
        dirs[2] = -4;
        dirs[3] = dir + 1;

        toMul.Mul(
            _deviceLinkSUN(pDeviceData, sSite4, 4, byFieldId, dirs)
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
template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateNablaENaiveSUN(
    const deviceSUN<N, NofE>* __restrict__ pDeviceData,
    BYTE byFieldId, deviceSUN<N, NofE>* pRes)
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

    pRes[uiResLinkIdx] = deviceSUN<N, NofE>::makeSUNZero();
    #pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {

        INT dirs[4];
        dirs[0] = 4;
        dirs[1] = dir + 1;
        dirs[2] = -4;
        dirs[3] = -static_cast<INT>(dir) - 1;


        deviceSUN<N, NofE> a(
            _deviceLinkSUN(pDeviceData, sSite4, 4, byFieldId, dirs)
        );
        
        dirs[0] = -static_cast<INT>(dir) - 1;
        dirs[1] = 4;
        dirs[2] = dir + 1;
        dirs[3] = -4;

        deviceSUN<N, NofE> b(
            _deviceLinkSUN(pDeviceData, sSite4, 4, byFieldId, dirs)
        );

        a.Sub(b);
        pRes[uiResLinkIdx].Sub(a);
    }
}

/**
 * U = exp(A)
 */
template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToUSUN(deviceSUN<N, NofE>* pDeviceData)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].ExpReal(F(1.0), (_DC_ExpPrecision <= N) ? (N + 1) : _DC_ExpPrecision);
    }
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToULogSUN(deviceSUN<N, NofE>* pDeviceData)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].ExpReal(F(1.0));
    }
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelSetOneDirUnitySUN(deviceSUN<N, NofE>* pDeviceData, BYTE byDir)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = deviceSUN<N, NofE>::makeSUNId();
        }
    }
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND
_kernelSetOneDirZeroSUN(deviceSUN<N, NofE>* pDeviceData, BYTE byDir)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = deviceSUN<N, NofE>::makeSUNZero();
        }
    }
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelSetOneDirUnityPointSUN(deviceSUN<N, NofE>* pDeviceData, UINT uiSiteIndex, BYTE byDir)
{
    for (BYTE dir = 0; dir < 4; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = deviceSUN<N, NofE>::makeSUNId();
        }
    }
}

template<INT N, INT NofE>
__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelSetOneDirZeroPointSUN(deviceSUN<N, NofE>* pDeviceData, UINT uiSiteIndex, BYTE byDir)
{
    for (BYTE dir = 0; dir < 4; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = deviceSUN<N, NofE>::makeSUNZero();
        }
    }
}

template<INT N, INT NoE>
__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovLoopOfSiteSUN(
    const deviceSUN<N, NoE>* __restrict__ pDeviceBuffer,
    cuDoubleComplex* res)
{
    UINT uiXYZ = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    const UINT uiSiteIndex = uiXYZ * _DC_Lt;
    UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, _DC_Dir - 1);
    SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

    deviceSUN<N, NoE> tmp = deviceSUN<N, NoE>::makeSUNZero();
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

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::AxpyPlus(const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSUN: axpy failed because the otherfield is not SUN with same N\n"); 
        return;
    }
    const CFieldGaugeSUN<N, NoE>* pSUNx = dynamic_cast<const CFieldGaugeSUN<N, NoE>*>(x);
    if (NULL == pSUNx)
    {
        appCrucial("CFieldGaugeSUN: axpy failed because the otherfield is not SUN with same N\n"); 
        return;
    }

    preparethread;
    _kernelAxpyPlusSUN << <block, threads >> > (m_pDeviceData, pSUNx->m_pDeviceData);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::AxpyMinus(const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSUN: axpy failed because the otherfield is not SUN with same N\n");
        return;
    }
    const CFieldGaugeSUN<N, NoE>* pSUNx = dynamic_cast<const CFieldGaugeSUN<N, NoE>*>(x);
    if (NULL == pSUNx)
    {
        appCrucial("CFieldGaugeSUN: axpy failed because the otherfield is not SUN with same N\n");
        return;
    }

    preparethread;
    _kernelAxpyMinusSUN << <block, threads >> > (m_pDeviceData, pSUNx->m_pDeviceData);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::ScalarMultply(const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplySUNComplex << <block, threads >> > (m_pDeviceData, a);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplySUNReal << <block, threads >> > (m_pDeviceData, a);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::Axpy(Real a, const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSUN: axpy failed because the otherfield is not SUN with same N\n");
        return;
    }
    const CFieldGaugeSUN<N, NoE>* pSUNx = dynamic_cast<const CFieldGaugeSUN<N, NoE>*>(x);
    if (NULL == pSUNx)
    {
        appCrucial("CFieldGaugeSUN: axpy failed because the otherfield is not SUN with same N\n");
        return;
    }

    preparethread;
    _kernelAxpySUNReal << <block, threads >> > (m_pDeviceData, pSUNx->m_pDeviceData, a);

}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSUN: axpy failed because the otherfield is not SUN with same N\n");
        return;
    }
    const CFieldGaugeSUN<N, NoE>* pSUNx = dynamic_cast<const CFieldGaugeSUN<N, NoE>*>(x);
    if (NULL == pSUNx)
    {
        appCrucial("CFieldGaugeSUN: axpy failed because the otherfield is not SUN with same N\n");
        return;
    }

    preparethread;
    _kernelAxpySUNA << <block, threads >> > (m_pDeviceData, pSUNx->m_pDeviceData, a);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::Mul(const CField* x, UBOOL bDagger)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSUN: axpy failed because the otherfield is not SUN with same N\n");
        return;
    }
    const CFieldGaugeSUN<N, NoE>* pSUNx = dynamic_cast<const CFieldGaugeSUN<N, NoE>*>(x);
    if (NULL == pSUNx)
    {
        appCrucial("CFieldGaugeSUN: axpy failed because the otherfield is not SUN with same N\n");
        return;
    }

    preparethread;
    _kernelMulSUN << <block, threads >> > (m_pDeviceData, pSUNx->m_pDeviceData, bDagger);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::Zero()
{
    preparethread;
    _kernelInitialSUNFeield << <block, threads >> > (m_pDeviceData, EFIT_Zero);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::Identity()
{
    preparethread;
    _kernelInitialSUNFeield << <block, threads >> > (m_pDeviceData, EFIT_Identity);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::Dagger()
{
    preparethread;
    _kernelDaggerSUN << <block, threads >> > (m_pDeviceData);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::MakeRandomGenerator()
{
    preparethread;
    _kernelInitialSUNFeield << <block, threads >> > (m_pDeviceData, EFIT_RandomGenerator);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::InitialField(EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialSUNFeield << <block, threads >> > (m_pDeviceData, eInitialType);
}

template<INT N, INT NoE>
DOUBLE CFieldGaugeSUN<N, NoE>::CalculateKinematicEnergy() const
{
    preparethread;
    _kernelCalculateKinematicEnergySUN << <block, threads >> > (m_pDeviceData, _D_RealThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::SetOneDirectionUnity(BYTE byDir)
{
    if (0 == (byDir & 15))
    {
        return;
    }
    preparethread;
    _kernelSetOneDirUnitySUN << <block, threads >> > (m_pDeviceData, byDir);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::SetOneDirectionZero(BYTE byDir)
{
    if (0 == (byDir & 15))
    {
        return;
    }
    preparethread;
    _kernelSetOneDirZeroSUN << <block, threads >> > (m_pDeviceData, byDir);
}

template<INT N, INT NoE>
CFieldGaugeSUN<N, NoE>::CFieldGaugeSUN() : CFieldGauge()
{
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData, sizeof(deviceSUN<N, NoE>) * m_uiLinkeCount));
}

template<INT N, INT NoE>
CFieldGaugeSUN<N, NoE>::~CFieldGaugeSUN()
{
    checkCudaErrors(__cudaFree(m_pDeviceData));
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eType)
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
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * N * N * m_uiLinkeCount);
        BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (uiSize != sizeof(Real) * 2 * N * N * _HC_LinkCount)
        {
            appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(Real) * 2 * N * N * _HC_LinkCount));
        }
        InitialWithByte(data);
        free(data);
    }
    break;
#if _CLG_DOUBLEFLOAT
    case EFFT_CLGBinFloat:
    {
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * N * N * m_uiLinkeCount);
        BYTE* data = (BYTE*)malloc(uiSize);
        Real* rdata = (Real*)data;
        FLOAT* fdata = (FLOAT*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (uiSize != sizeof(FLOAT) * 2 * N * N * _HC_LinkCount)
        {
            appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, sizeof(FLOAT) * 2 * N * N * _HC_LinkCount);
        }
        for (UINT i = 0; i < 2 * N * N * m_uiLinkeCount; ++i)
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
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * N * N * m_uiLinkeCount);
        BYTE* data = (BYTE*)malloc(uiSize);
        Real* rdata = (Real*)data;
        DOUBLE* ddata = (DOUBLE*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (uiSize != sizeof(DOUBLE) * 2 * N * N * _HC_LinkCount)
        {
            appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(DOUBLE) * 2 * N * N * _HC_LinkCount));
        }
        for (UINT i = 0; i < 2 * N * N * m_uiLinkeCount; ++i)
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

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::InitialWithByte(BYTE* byData)
{
    deviceSUN<N, NoE>* readData = (deviceSUN<N, NoE>*)malloc(sizeof(deviceSUN<N, NoE>) * m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[2 * N * N];
        memcpy(oneLink, byData + sizeof(Real) * 2 * N * N * i, sizeof(Real) * 2 * N * N);
        for (UINT j = 0; j < NoE; ++j)
        {
            if (j < N * N)
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
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceSUN<N, NoE>) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    free(readData);
}

/**
* (1) calculate staples
* (2) calculate force(additive)
* (3) calculate energy
*/
template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, Real betaOverN) const
{
    if (NULL == pForce || GetFieldType() != pForce->GetFieldType())
    {
        appCrucial("CFieldGaugeSUN: force field is not SU3");
        return;
    }
    if (NULL != pStable && GetFieldType() != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSUN: stape field is not SUN");
        return;
    }

    CFieldGaugeSUN<N, NoE>* pForceSU3 = dynamic_cast<CFieldGaugeSUN<N, NoE>*>(pForce);
    CFieldGaugeSUN<N, NoE>* pStableSU3 = NULL == pStable ? NULL : dynamic_cast<CFieldGaugeSUN<N, NoE>*>(pStable);

    preparethread;

    assert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache);

    _kernelStapleAtSiteSUNCacheIndex << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData,
        betaOverN);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::CalculateOnlyStaple(CFieldGauge* pStable) const
{
    if (NULL == pStable || GetFieldType() != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSUN: stable field is not SUN");
        return;
    }
    CFieldGaugeSUN<N, NoE>* pStableSU3 = dynamic_cast<CFieldGaugeSUN<N, NoE>*>(pStable);

    preparethread;
    _kernelCalculateOnlyStapleSUN << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStableSU3->m_pDeviceData);
}

template<INT N, INT NoE>
DOUBLE CFieldGaugeSUN<N, NoE>::CalculatePlaqutteEnergy(DOUBLE betaOverN) const
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);

    preparethread;
    _kernelPlaqutteEnergySUNCacheIndex << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer
        );

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<INT N, INT NoE>
DOUBLE CFieldGaugeSUN<N, NoE>::CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);
    preparethread;
    _kernelPlaqutteEnergySUN_UseClover << <block, threads >> > (
        m_byFieldId,
        m_pDeviceData,
        static_cast<DOUBLE>(N) * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<INT N, INT NoE>
DOUBLE CFieldGaugeSUN<N, NoE>::CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge *pStable) const
{
    if (NULL == pStable || GetFieldType() != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stape field is not SU3");
        return F(0.0);
    }
    const CFieldGaugeSUN<N, NoE>* pStableSUN = dynamic_cast<const CFieldGaugeSUN<N, NoE>*>(pStable);
    preparethread;
    _kernelPlaqutteEnergyUsingStableSUN << <block, threads >> > (
        m_pDeviceData, 
        pStableSUN->m_pDeviceData,
        static_cast<DOUBLE>(N) * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        betaOverN, 
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::ExpMult(Real a, CField* U) const
{
    if (NULL == U || GetFieldType() != U->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: U field is not SU3");
        return;
    }

    CFieldGaugeSUN<N, NoE>* pUField = dynamic_cast<CFieldGaugeSUN<N, NoE>*>(U);

    preparethread;
    _kernelExpMultSUNReal << < block, threads >> > (m_pDeviceData, a, pUField->m_pDeviceData, static_cast<BYTE>(_HC_ExpPrecision));
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::ElementNormalize()
{
    preparethread;
    _kernelNormalizeSUN << < block, threads >> > (m_pDeviceData);
    //DebugPrintMe();
}

template<INT N, INT NoE>
cuDoubleComplex CFieldGaugeSUN<N, NoE>::Dot(const CField* other) const
{
    if (NULL == other || GetFieldType() != other->GetFieldType())
    {
        appCrucial("CFieldGaugeSUN: U field is not SU3");
        return make_cuDoubleComplex(0,0);
    }

    const CFieldGaugeSUN<N, NoE>* pUField = dynamic_cast<const CFieldGaugeSUN<N, NoE>*>(other);

    preparethread;
    _kernelDotSUN << < block, threads >> > (m_pDeviceData, pUField->m_pDeviceData, _D_ComplexThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::DebugPrintMe() const
{
    //preparethread;
    //_kernelPrintSU3 << < block, threads >> > (m_pDeviceData);

    //===================================================
    //Since Debug Print Me is only used to debug, we do it slow but convinient
    deviceSUN<N, NoE>* pToPrint = (deviceSUN<N, NoE>*)malloc(sizeof(deviceSUN<N, NoE>) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(pToPrint, m_pDeviceData, sizeof(deviceSUN<N, NoE>) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    appPushLogDate(FALSE);

    for (UINT uiLink = 0; uiLink < m_uiLinkeCount; ++uiLink)
    {
        UINT uiSite = uiLink / _HC_Dir;
        UINT uiDir = uiLink % _HC_Dir;
        SSmallInt4 site = __hostSiteIndexToInt4(uiSite);
        appGeneral(_T(" --- %d(%d, %d, %d, %d)_%d ---\n {{"), uiLink, static_cast<INT>(site.x), static_cast<INT>(site.y), static_cast<INT>(site.z), static_cast<INT>(site.w), uiDir);
        for (INT y = 0; y < N; ++y)
        {
            for (INT x = 0; x < N; ++x)
            {
                appGeneral(_T("%f %s %f I"), pToPrint[uiLink].m_me[y * N + x].x, pToPrint[uiLink].m_me[y * N + x].y > F(0.0) ? _T("+") : _T("-"), appAbs(pToPrint[uiLink].m_me[y * N + x].y));
                if (x != (N - 1))
                {
                    appGeneral(_T(", "));
                }
                else
                {
                    if (y != (N - 1))
                    {
                        appGeneral(_T("},\n {"));
                    }
                }
            }
        }
        appGeneral(_T("}}\n"));
    }
    appPopLogDate();
    free(pToPrint);
}

template<INT N, INT NoE>
BYTE* CFieldGaugeSUN<N, NoE>::CopyDataOut(UINT &uiSize) const
{
    deviceSUN<N, NoE>* toSave = (deviceSUN<N, NoE>*)malloc(sizeof(deviceSUN<N, NoE>) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceSUN<N, NoE>) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount * 2 * N * N);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[2 * N * N];
        for (UINT j = 0; j < N * N; ++j)
        {
            oneLink[2 * j] = static_cast<Real>(toSave[i].m_me[j].x);
            oneLink[2 * j + 1] = static_cast<Real>(toSave[i].m_me[j].y);
        }
        memcpy(byToSave + i * sizeof(Real) * 2 * N * N, oneLink, sizeof(Real) * 2 * N * N);
    }
    free(toSave);

    return byToSave;
}

template<INT N, INT NoE>
BYTE* CFieldGaugeSUN<N, NoE>::CopyDataOutFloat(UINT& uiSize) const
{
    deviceSUN<N, NoE>* toSave = (deviceSUN<N, NoE>*)malloc(sizeof(deviceSUN<N, NoE>) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceSUN<N, NoE>) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(FLOAT) * m_uiLinkeCount * 2 * N * N);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        FLOAT oneLink[2 * N * N];
        for (UINT j = 0; j < N * N; ++j)
        {
            oneLink[2 * j] = static_cast<FLOAT>(toSave[i].m_me[j].x);
            oneLink[2 * j + 1] = static_cast<FLOAT>(toSave[i].m_me[j].y);
        }
        memcpy(byToSave + i * sizeof(FLOAT) * 2 * N * N, oneLink, sizeof(FLOAT) * 2 * N * N);
    }
    free(toSave);

    return byToSave;
}

template<INT N, INT NoE>
BYTE* CFieldGaugeSUN<N, NoE>::CopyDataOutDouble(UINT& uiSize) const
{
    deviceSUN<N, NoE>* toSave = (deviceSUN<N, NoE>*)malloc(sizeof(deviceSUN<N, NoE>) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceSUN<N, NoE>) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(DOUBLE) * m_uiLinkeCount * 2 * N * N);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        DOUBLE oneLink[2 * N * N];
        for (UINT j = 0; j < N * N; ++j)
        {
            oneLink[2 * j] = static_cast<DOUBLE>(toSave[i].m_me[j].x);
            oneLink[2 * j + 1] = static_cast<DOUBLE>(toSave[i].m_me[j].y);
        }
        memcpy(byToSave + i * sizeof(DOUBLE) * 2 * N * N, oneLink, sizeof(DOUBLE) * 2 * N * N);
    }
    free(toSave);

    return byToSave;
}

template<INT N, INT NoE>
void CFieldGaugeSUN<N, NoE>::PolyakovOnSpatialSite(cuDoubleComplex* buffer) const
{
    dim3 block(_HC_DecompX, _HC_DecompY, 1);
    dim3 threads(_HC_DecompLx, _HC_DecompLy, 1);
    _kernelPolyakovLoopOfSiteSUN << <block, threads >> > (m_pDeviceData, buffer);
}

__CLGIMPLEMENT_CLASS(CFieldGaugeSU4)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU5)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU6)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU7)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU8)

void CFieldGaugeSU4::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || EFT_GaugeSU4 != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: target field is not SU3");
        return;
    }

    CFieldGauge::CopyTo(pTarget);

    CFieldGaugeSU4* pTargetField = dynamic_cast<CFieldGaugeSU4*>(pTarget);
    checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceData, m_pDeviceData, sizeof(deviceSU4) * m_uiLinkeCount, cudaMemcpyDeviceToDevice));
}

void CFieldGaugeSU5::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || EFT_GaugeSU5 != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: target field is not SU3");
        return;
    }

    CFieldGauge::CopyTo(pTarget);

    CFieldGaugeSU5* pTargetField = dynamic_cast<CFieldGaugeSU5*>(pTarget);
    checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceData, m_pDeviceData, sizeof(deviceSU5) * m_uiLinkeCount, cudaMemcpyDeviceToDevice));
}

void CFieldGaugeSU6::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || EFT_GaugeSU6 != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: target field is not SU3");
        return;
    }

    CFieldGauge::CopyTo(pTarget);

    CFieldGaugeSU6* pTargetField = dynamic_cast<CFieldGaugeSU6*>(pTarget);
    checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceData, m_pDeviceData, sizeof(deviceSU6) * m_uiLinkeCount, cudaMemcpyDeviceToDevice));
}

void CFieldGaugeSU7::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || EFT_GaugeSU7 != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: target field is not SU3");
        return;
    }

    CFieldGauge::CopyTo(pTarget);

    CFieldGaugeSU7* pTargetField = dynamic_cast<CFieldGaugeSU7*>(pTarget);
    checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceData, m_pDeviceData, sizeof(deviceSU7) * m_uiLinkeCount, cudaMemcpyDeviceToDevice));
}

void CFieldGaugeSU8::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || EFT_GaugeSU8 != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: target field is not SU3");
        return;
    }

    CFieldGauge::CopyTo(pTarget);

    CFieldGaugeSU8* pTargetField = dynamic_cast<CFieldGaugeSU8*>(pTarget);
    checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceData, m_pDeviceData, sizeof(deviceSU8) * m_uiLinkeCount, cudaMemcpyDeviceToDevice));
}


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================