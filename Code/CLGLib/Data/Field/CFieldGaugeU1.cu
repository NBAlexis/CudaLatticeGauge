//=============================================================================
// FILENAME : CFieldGaugeU1.cu
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
//  [10/13/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldGaugeU1)

#pragma region Kernels

/**
* Initial SU3 Field with a value
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialU1Field(CLGComplex *pDevicePtr, EFieldInitialType eInitialType)
{
    CLGComplex id = _onec;
    CLGComplex zero = _zeroc;

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
            const Real fArg = _deviceRandomF(_deviceGetFatIndex(uiSiteIndex, idir + 1)) * PI2;
            pDevicePtr[uiLinkIndex] = _make_cuComplex(_cos(fArg), -_sin(fArg));
        }
        break;
        case EFIT_RandomGenerator:
        {
            const Real r1 = _deviceRandomGaussFSqrt2(_deviceGetFatIndex(uiSiteIndex, idir + 1)) * PI2;
            pDevicePtr[uiLinkIndex] = _make_cuComplex(_cos(r1), -_sin(r1));
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
_kernelDaggerU1(CLGComplex* pDevicePtr)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex] = _cuConjf(pDevicePtr[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyU1A(CLGComplex*pDevicePtr, const CLGComplex* __restrict__ x, CLGComplex a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex] = _cuCaddf(pDevicePtr[uiLinkIndex], cuCmulf(x[uiLinkIndex], a));

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyU1Real(CLGComplex*pDevicePtr, const CLGComplex* __restrict__ x, Real a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex] = _cuCaddf(pDevicePtr[uiLinkIndex], cuCmulf_cr(x[uiLinkIndex], a));

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusU1(CLGComplex*pDevicePtr, const CLGComplex* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex] = _cuCaddf(pDevicePtr[uiLinkIndex], x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusU1(CLGComplex*pDevicePtr, const CLGComplex* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex] = _cuCsubf(pDevicePtr[uiLinkIndex], x[uiLinkIndex]);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyU1Complex(CLGComplex*pDevicePtr, CLGComplex a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex] = _cuCmulf(pDevicePtr[uiLinkIndex], a);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyU1Real(CLGComplex*pDevicePtr, Real a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex] = cuCmulf_cr(pDevicePtr[uiLinkIndex], a);

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteU1CacheIndex(
    const CLGComplex * __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    CLGComplex*pStapleData, //can be NULL
    CLGComplex*pForceData,
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
        CLGComplex res = _zeroc;

        //there are 6 staples, each is sum of two plaquttes
        for (INT i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            CLGComplex toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

            if (first.NeedToDagger())
            {
                toAdd.y = -toAdd.y;
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                CLGComplex toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);

                if (nextlink.NeedToDagger())
                {
                    toAdd = _cuCmulf(toAdd, _make_cuComplex(toMul.x, -toMul.y));
                }
                else
                {
                    toAdd = _cuCmulf(toAdd, toMul);
                }
            }
            res = _cuCaddf(res, toAdd);
        }
        if (NULL != pStapleData)
        {
            pStapleData[linkIndex] = res;
        }

        //staple calculated
        CLGComplex force(pDeviceData[linkIndex]);
        //force = _cuCmulf(_make_cuComplex(res.x, -res.y));
        //force.Ta();
        //force.MulReal(betaOverN);
        force = _make_cuComplex(F(0.0), betaOverN * (pDeviceData[linkIndex].y * res.x - pDeviceData[linkIndex].x * res.y));
        

        //force is additive
        pForceData[linkIndex] = _cuCaddf(pForceData[linkIndex], force);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateOnlyStapleU1(
    const CLGComplex* __restrict__ pDeviceData,
    const SIndex * __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    CLGComplex*pStapleData)
{
    intokernaldir;

    //Real test_force = F(0.0);
    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        CLGComplex res = _zeroc;

        //there are 6 staples, each is sum of two plaquttes
        for (INT i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            CLGComplex toAdd(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);

            if (first.NeedToDagger())
            {
                toAdd.y = -toAdd.y;
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                CLGComplex toMul(pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)]);

                if (nextlink.NeedToDagger())
                {
                    toAdd = _cuCmulf(toAdd, _make_cuComplex(toMul.x, -toMul.y));
                }
                else
                {
                    toAdd = _cuCmulf(toAdd, toMul);
                }
            }
            res = _cuCaddf(res, toAdd);
        }
        if (NULL != pStapleData)
        {
            pStapleData[linkIndex] = res;
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyU1CacheIndex(
    const CLGComplex* __restrict__ pDeviceData,
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
    UINT plaqCountAll = plaqCount * plaqLength;
    for (BYTE i = 0; i < plaqCount; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLength + uiSiteIndex * plaqCountAll];
#if !_CLG_DOUBLEFLOAT
        cuDoubleComplex toAdd = _cToDouble(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
#else
        CLGComplex toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];
#endif

        if (first.NeedToDagger())
        {
            toAdd.y = -toAdd.y;
        }

        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedIndex[i * plaqLength + j + uiSiteIndex * plaqCountAll];

#if !_CLG_DOUBLEFLOAT
            cuDoubleComplex toMul = _cToDouble(pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)]);
            if (first.NeedToDagger())
            {
                toAdd = cuCmul(toAdd, make_cuDoubleComplex(toMul.x, -toMul.y));
            }
            else
            {
                toAdd = cuCmul(toAdd, toMul);
            }
#else
            CLGComplex toMul = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];
            if (first.NeedToDagger())
            {
                toAdd = _cuCmulf(toAdd, _make_cuComplex(toMul.x, -toMul.y));
            }
            else
            {
                toAdd = _cuCmulf(toAdd, toMul);
            }
#endif
        }

#if !_CLG_DOUBLEFLOAT
        resThisThread += (1.0 - toAdd.x);
#else
        resThisThread += (F(1.0) - toAdd.x);
#endif
    }

    results[uiSiteIndex] = resThisThread * betaOverN;

    //printf("  ---- energy: thread=%d, res=%f\n", __thread_id, results[__thread_id]);
}

/*
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyU1_UseClover(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE fBetaOverN,
    DOUBLE* results
#else
    Real fBetaOverN,
    Real* results
#endif
)
{
    intokernalInt4;

    Real fRes = F(0.0);
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
    fRes = F(18.0) - F(0.25) * fRes;
    results[uiSiteIndex] = fRes * fBetaOverN;// *F(0.5);
}
*/

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyUsingStableU1(
    const CLGComplex* __restrict__ pDeviceData,
    const CLGComplex* __restrict__ pStableData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN,
    DOUBLE* results
#else
    Real betaOverN,
    Real* results
#endif
)
{
    intokernaldir;

#if !_CLG_DOUBLEFLOAT
    DOUBLE resThisThread = 0.0;
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        //For each link, there are 6 staples
        resThisThread += (6.0 -(pDeviceData[linkIndex].x * pStableData[linkIndex].x + pDeviceData[linkIndex].y * pStableData[linkIndex].y));
    }

    results[uiSiteIndex] = resThisThread * betaOverN * 0.25;
#else
    Real resThisThread = F(0.0);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        //For each link, there are 6 staples
        resThisThread += (F(6.0) - (pDeviceData[linkIndex].x * pStableData[linkIndex].x + pDeviceData[linkIndex].y * pStableData[linkIndex].y));
    }

    results[uiSiteIndex] = resThisThread * betaOverN * F(0.25);
#endif
    //printf("  ---- energy: thread=%d, res=%f\n", __thread_id, results[__thread_id]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelExpMultU1Real(
    const CLGComplex* __restrict__ pMyDeviceData,
    Real a,
    CLGComplex*pU)
{
    intokernaldir;

    for (BYTE idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        pU[linkIndex] = _cuCmulf(__cuCexpf(cuCmulf_cr(pMyDeviceData[linkIndex], a)), pU[linkIndex]);
    }
}

/**
* Trace (P^2)
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateKinematicEnergyU1(const CLGComplex* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE* results
#else
    Real* results
#endif
)
{
    intokernaldir;

#if !_CLG_DOUBLEFLOAT
    DOUBLE resThisThread = F(0.0);
#else
    Real resThisThread = F(0.0);
#endif
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        resThisThread += __cuCabsSqf(pDeviceData[linkIndex]);
    }

    results[uiSiteIndex] = resThisThread;
}


__global__ void _CLG_LAUNCH_BOUND
_kernelNormalizeU1(CLGComplex* pMyDeviceData)
{
    intokernaldir;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        Real fArg = __cuCargf(pMyDeviceData[linkIndex]);
        pMyDeviceData[linkIndex] = _make_cuComplex(_cos(fArg), _sin(fArg));
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDotU1(
    const CLGComplex* __restrict__ pMyDeviceData,
    const CLGComplex* __restrict__ pOtherDeviceData,
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
            make_cuDoubleComplex(pMyDeviceData[linkIndex].x * pOtherDeviceData[linkIndex].x + pMyDeviceData[linkIndex].y * pOtherDeviceData[linkIndex].y, 
                pMyDeviceData[linkIndex].x * pOtherDeviceData[linkIndex].y - pMyDeviceData[linkIndex].y * pOtherDeviceData[linkIndex].x)
        );
    }

    result[uiSiteIndex] = resThisThread;
#else
    CLGComplex resThisThread = _make_cuComplex(0,0);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        resThisThread = _cuCaddf(resThisThread, 
            _make_cuComplex(pMyDeviceData[linkIndex].x * pOtherDeviceData[linkIndex].x + pMyDeviceData[linkIndex].y * pOtherDeviceData[linkIndex].y,
                pMyDeviceData[linkIndex].x * pOtherDeviceData[linkIndex].y - pMyDeviceData[linkIndex].y * pOtherDeviceData[linkIndex].x)
        );
    }

    result[uiSiteIndex] = resThisThread;
#endif
    //printf("res = %f %f\n", pOtherDeviceData[uiSiteIndex * 4].m_me[0].x, pMyDeviceData[uiSiteIndex * 4].m_me[0].x);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSetConfigurationU1(
    CLGComplex* pDeviceData,
    const Real* __restrict__ pRealData)
{
    gaugeSU3KernelFuncionStart

    pDeviceData[uiLinkIndex] = _make_cuComplex(pRealData[2 * uiLinkIndex +  0], pRealData[2 * uiLinkIndex +  1]);

    gaugeSU3KernelFuncionEnd
}

/**
 * iA = U.TA() 
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToIAU1(
    CLGComplex* pDeviceData)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        pDeviceData[uiLinkIndex] = __cuClogf(pDeviceData[uiLinkIndex]);
    }
}

/**
 * E_mu = F_{0 mu}
 */
//__global__ void _CLG_LAUNCH_BOUND
//_kernelTransformToEU1(
//    BYTE byFieldId,
//    const CLGComplex* __restrict__ pDeviceData,
//    CLGComplex* pRes)
//{
//    intokernalInt4;
//    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
//    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
//    
//    //we only need uiDir = 0,1,2
//    for (BYTE dir = 0; dir < uiDir; ++dir)
//    {
//        CLGComplex res = _zeroc;
//        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
//        if (dir < uiDir - 1)
//        {
//            //find clover F
//            //res = _deviceClover(pDeviceData, uiBigIdx, 3, dir);
//            //test not using clover
//            res = _device1PlaqutteTermPP(pDeviceData, 3, dir, uiBigIdx, sSite4, byFieldId);
//            //res.iIm2();
//            //not using clover not multiply 0.25
//            //res.MulReal(F(0.125));
//
//            res.Ta();
//            //res.MulReal(F(0.25));
//        }
//
//        pRes[uiLinkIndex] = res;
//    }
//}

/**
 * This is wrong! the order of the plaqutte must be considered
 * This is to make sure gauge transform is g(x) nabla E g^+(n)
 */
//__global__ void _CLG_LAUNCH_BOUND
//_kernelCalculateNablaE(
//    const deviceSU3* __restrict__ pDeviceData,
//    BYTE byFieldId, deviceSU3* pRes)
//{
//    intokernalInt4;
//    //const BYTE uiDir2 = static_cast<BYTE>(_DC_Dir) * 2;
//    //const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
//
//
//    //i=0: 12
//    //  1: 13
//    //  2: 14
//    //  3: 23
//    //  4: 24
//    //  5: 34  
//
//    UINT uiResLinkIdx = _deviceGetLinkIndex(uiSiteIndex, 3);
//
//    pRes[uiResLinkIdx] = deviceSU3::makeSU3Zero();
//    #pragma unroll
//    for (BYTE dir = 0; dir < 3; ++dir)
//    {
//        //we need 2, 4 and 5
//        //BYTE byPlaqIdx = (dir + 1) << 1;
//        //if (byPlaqIdx > 5) byPlaqIdx = 5;
//
//        INT dirs[4];
//        //deviceSU3 toMul(_devicePlaqutte(pDeviceData, pCachedPlaqutte, uiSiteIndex, byPlaqIdx, plaqLength, plaqCount));
//        dirs[0] = dir + 1;
//        dirs[1] = 4;
//        dirs[2] = dir;
//        dirs[2] = -dirs[2] - 1;
//        dirs[3] = -4;
//        deviceSU3 toMul(
//            //_device1PlaqutteTermPP(pDeviceData, 3, dir, uiBigIdx)
//            _deviceLink(pDeviceData, sSite4, 4, byFieldId, dirs)
//        );
//
//        //
//        dirs[0] = 4;
//        dirs[1] = dir;
//        dirs[1] = -dirs[1] - 1;
//        dirs[2] = -4;
//        dirs[3] = dir + 1;
//        toMul.MulDagger(
//            _deviceLink(pDeviceData, sSite4, 4, byFieldId, dirs)
//        );
//        pRes[uiResLinkIdx].Add(toMul);
//    }
//
//    pRes[uiResLinkIdx].SubReal(F(3.0));
//}

/**
 * U = exp(A)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToU_U1(CLGComplex* pDeviceData)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        pDeviceData[uiLinkIndex] = __cuCexpf(pDeviceData[uiLinkIndex]);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSetOneDirUnity_U1(CLGComplex* pDeviceData, BYTE byDir)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = _onec;
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSetOneDirZero_U1(CLGComplex* pDeviceData, BYTE byDir)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = _zeroc;
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelSetOneDirUnityPoint_U1(CLGComplex* pDeviceData, UINT uiSiteIndex, BYTE byDir)
{
    for (BYTE dir = 0; dir < 4; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = _onec;
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelSetOneDirZeroPoint_U1(CLGComplex* pDeviceData, UINT uiSiteIndex, BYTE byDir)
{
    for (BYTE dir = 0; dir < 4; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = _zeroc;
        }
    }
}

#pragma endregion

void CFieldGaugeU1::AxpyPlus(const CField* x)
{
    if (NULL == x || EFT_GaugeU1 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeU1* pSU3x = dynamic_cast<const CFieldGaugeU1*>(x);
    preparethread;
    _kernelAxpyPlusU1 << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData);
}

void CFieldGaugeU1::AxpyMinus(const CField* x)
{
    if (NULL == x || EFT_GaugeU1 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeU1* pSU3x = dynamic_cast<const CFieldGaugeU1*>(x);
    preparethread;
    _kernelAxpyMinusU1 << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData);

}

void CFieldGaugeU1::ScalarMultply(const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplyU1Complex << <block, threads >> > (m_pDeviceData, a);
}

void CFieldGaugeU1::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplyU1Real << <block, threads >> > (m_pDeviceData, a);
}

void CFieldGaugeU1::Axpy(Real a, const CField* x)
{
    if (NULL == x || EFT_GaugeU1 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeU1* pSU3x = dynamic_cast<const CFieldGaugeU1*>(x);
    preparethread;
    _kernelAxpyU1Real << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData, a);
}

void CFieldGaugeU1::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || EFT_GaugeU1 != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeU1* pSU3x = dynamic_cast<const CFieldGaugeU1*>(x);
    preparethread;
    _kernelAxpyU1A << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData, a);
}


void CFieldGaugeU1::Zero()
{
    preparethread;
    _kernelInitialU1Field << <block, threads >> > (m_pDeviceData, EFIT_Zero);
}

void CFieldGaugeU1::Identity()
{
    preparethread;
    _kernelInitialU1Field << <block, threads >> > (m_pDeviceData, EFIT_Identity);
}

void CFieldGaugeU1::Dagger()
{
    preparethread;
    _kernelDaggerU1 << <block, threads >> > (m_pDeviceData);
}

void CFieldGaugeU1::MakeRandomGenerator()
{
    preparethread;
    _kernelInitialU1Field << <block, threads >> > (m_pDeviceData, EFIT_RandomGenerator);
}

/**
*
*/
void CFieldGaugeU1::InitialField(EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialU1Field << <block, threads >> > (m_pDeviceData, eInitialType);
}

void CFieldGaugeU1::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eType)
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
        appCrucial(_T("U1 EFFT_BridgePPTXT Not supported!\n"));
        //const CCString sContent = appGetFileSystem()->ReadAllText(sFileName);
        //TArray<INT> seps;
        //seps.AddItem(_T('\n'));
        //seps.AddItem(_T('\r'));
        //TArray<CCString> sStringlist = appGetStringList(sContent, seps, 0x7fffffff);
        //assert(static_cast<UINT>(sStringlist.Num()) == _HC_LinkCount * 18);

        //Real* pData = (Real*)malloc(sizeof(Real) * sStringlist.Num());
        //for (INT i = 0; i < sStringlist.Num(); ++i)
        //{
        //    pData[i] = static_cast<Real>(appStrToDOUBLE(sStringlist[i]));
        //}

        //SetByArray(pData);
    }
    break;
    case EFFT_BridgePPBin:
    {
        appCrucial(_T("U1 EFFT_BridgePPBin Not supported!\n"));
        //UINT uiSize = 0;
        //BYTE* allBytes = appGetFileSystem()->ReadAllBytes(sFileName, uiSize);
        //assert(uiSize == 8 * 18 * _HC_LinkCount);
        //Real* pData = (Real*)malloc(sizeof(Real) * 18 * _HC_LinkCount);
        //for (UINT i = 0; i < 18 * _HC_LinkCount; ++i)
        //{
        //    BYTE data[8];
        //    for (INT j = 0; j < 8; ++j)
        //    {
        //        data[j] = allBytes[i * 8 + (7 - j)];
        //    }
        //    DOUBLE * dbData = (DOUBLE*)data;
        //    pData[i] = static_cast<Real>(*dbData);
        //}
        //free(allBytes);

        //SetByArray(pData);
    }
    break;
    case EFFT_CLGBin:
#if _CLG_DOUBLEFLOAT
    case EFFT_CLGBinDouble:
#else
    case EFFT_CLGBinFloat:
#endif
    {
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * m_uiLinkeCount);
        BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
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
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * m_uiLinkeCount);
        BYTE* data = (BYTE*)malloc(uiSize);
        Real* rdata = (Real*)data;
        DOUBLE* ddata = (DOUBLE*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        for (UINT i = 0; i < 2 * m_uiLinkeCount; ++i)
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
        appCrucial(_T("U1 EFFT_CLGBinCompressed Not supported!\n"));
        //UINT uiSize = static_cast<UINT>(sizeof(Real) * 9 * m_uiLinkeCount);
        //BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        //InitialWithByteCompressed(data);
        //free(data);
    }
    break;
    default:
        appCrucial(_T("Not supported input file type %s\n"), __ENUM_TO_STRING(EFieldFileType, eType).c_str());
        break;

    }
}

void CFieldGaugeU1::InitialWithByte(BYTE* byData)
{
    CLGComplex* readData = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[2];
        memcpy(oneLink, byData + sizeof(Real) * 2 * i, sizeof(Real) * 2);
        readData[i] = _make_cuComplex(
            oneLink[2 * i],
            oneLink[2 * i + 1]);
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(CLGComplex) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    free(readData);
}

void CFieldGaugeU1::InitialWithByteCompressed(BYTE* byData)
{
    appCrucial(_T("U1 InitialWithByteCompressed Not supported!\n"));
    _FAIL_EXIT;

    //deviceSU3* readData = (deviceSU3*)malloc(sizeof(deviceSU3) * m_uiLinkeCount);
    //for (UINT i = 0; i < m_uiLinkeCount; ++i)
    //{
    //    Real oneLink[9];
    //    memcpy(oneLink, byData + sizeof(Real) * 9 * i, sizeof(Real) * 9);

    //    readData[i].m_me[1] = _make_cuComplex(oneLink[0], oneLink[1]);
    //    readData[i].m_me[3] = _make_cuComplex(-oneLink[0], oneLink[1]);

    //    readData[i].m_me[2] = _make_cuComplex(oneLink[2], oneLink[3]);
    //    readData[i].m_me[6] = _make_cuComplex(-oneLink[2], oneLink[3]);

    //    readData[i].m_me[5] = _make_cuComplex(oneLink[4], oneLink[5]);
    //    readData[i].m_me[7] = _make_cuComplex(-oneLink[4], oneLink[5]);

    //    readData[i].m_me[0] = _make_cuComplex(F(0.0), oneLink[6]);
    //    readData[i].m_me[4] = _make_cuComplex(F(0.0), oneLink[7]);
    //    readData[i].m_me[8] = _make_cuComplex(F(0.0), oneLink[8]);

    //    for (UINT j = 9; j < 16; ++j)
    //    {
    //        readData[i].m_me[j] = _zeroc;
    //    }
    //}
    //checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(CLGComplex) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    //checkCudaErrors(cudaDeviceSynchronize());
    //free(readData);

    //preparethread;
    //_kernelTransformToU_U1 << <block, threads >> > (m_pDeviceData);
    //checkCudaErrors(cudaDeviceSynchronize());
}

void CFieldGaugeU1::SetByArray(Real* array)
{
    assert(NULL != array);
    //we algin the su3 now
    //assert(sizeof(deviceSU3) == 32 * sizeof(Real));

    //checkCudaErrors(cudaMemcpy(m_pDeviceData, array, sizeof(Real) * _HC_LinkCount * 18, cudaMemcpyHostToDevice));
    
    Real* pDeviceArray;
    checkCudaErrors(__cudaMalloc((void**)&pDeviceArray, sizeof(Real) * _HC_LinkCount * 2));
    checkCudaErrors(cudaMemcpy(pDeviceArray, array, sizeof(Real) * _HC_LinkCount * 2, cudaMemcpyHostToDevice));
    preparethread;
    _kernelSetConfigurationU1 << <block, threads >> > (m_pDeviceData, pDeviceArray);
    checkCudaErrors(__cudaFree(pDeviceArray));

    free(array);

    ElementNormalize();
}

/**
* (1) calculate staples
* (2) calculate force(additive)
* (3) calculate energy
*/
void CFieldGaugeU1::CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, Real betaOverN) const
{
    if (NULL == pForce || EFT_GaugeU1 != pForce->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: force field is not SU3");
        return;
    }
    if (NULL != pStable && EFT_GaugeU1 != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stape field is not SU3");
        return;
    }

    CFieldGaugeU1* pForceSU3 = dynamic_cast<CFieldGaugeU1*>(pForce);
    CFieldGaugeU1* pStableSU3 = NULL == pStable ? NULL : dynamic_cast<CFieldGaugeU1*>(pStable);

    preparethread;

    assert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache);

    _kernelStapleAtSiteU1CacheIndex << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData,
        betaOverN);
}

void CFieldGaugeU1::CalculateOnlyStaple(CFieldGauge* pStable) const
{
    if (NULL == pStable || EFT_GaugeU1 != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stable field is not SU3");
        return;
    }
    CFieldGaugeU1* pStableSU3 = dynamic_cast<CFieldGaugeU1*>(pStable);

    preparethread;
    _kernelCalculateOnlyStapleU1 << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStableSU3->m_pDeviceData);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeU1::CalculatePlaqutteEnergy(DOUBLE betaOverN) const
#else
Real CFieldGaugeU1::CalculatePlaqutteEnergy(Real betaOverN) const
#endif
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);

    preparethread;
    _kernelPlaqutteEnergyU1CacheIndex << <block, threads >> > (
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
DOUBLE CFieldGaugeU1::CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const
#else
Real CFieldGaugeU1::CalculatePlaqutteEnergyUseClover(Real betaOverN) const
#endif
{
    appCrucial(_T("U1 CalculatePlaqutteEnergyUseClover Not supported!\n"));
    _FAIL_EXIT;
    //assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);

    //preparethread;
    //_kernelPlaqutteEnergySU3_UseClover << <block, threads >> > (
    //    m_byFieldId,
    //    m_pDeviceData,
    //    betaOverN,
    //    _D_RealThreadBuffer);

    //return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeU1::CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge *pStable) const
#else
Real CFieldGaugeU1::CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge* pStable) const
#endif
{
    if (NULL == pStable || EFT_GaugeU1 != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stape field is not SU3");
        return F(0.0);
    }
    const CFieldGaugeU1* pStableSU3 = dynamic_cast<const CFieldGaugeU1*>(pStable);

    preparethread;
    _kernelPlaqutteEnergyUsingStableU1 << <block, threads >> > (
        m_pDeviceData, 
        pStableSU3->m_pDeviceData, 
        betaOverN, 
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeU1::CalculateKinematicEnergy() const
#else
Real CFieldGaugeU1::CalculateKinematicEnergy() const
#endif
{
    preparethread;
    _kernelCalculateKinematicEnergyU1 << <block, threads >> > (m_pDeviceData, _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

void CFieldGaugeU1::SetOneDirectionUnity(BYTE byDir)
{
    if (0 == (byDir & 15))
    {
        return;
    }
    preparethread;
    _kernelSetOneDirUnity_U1 << <block, threads >> >(m_pDeviceData, byDir);

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

void CFieldGaugeU1::SetOneDirectionZero(BYTE byDir)
{
    if (0 == (byDir & 15))
    {
        return;
    }
    preparethread;
    _kernelSetOneDirZero_U1 << <block, threads >> > (m_pDeviceData, byDir);
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

CFieldGaugeU1::CFieldGaugeU1() : CFieldGauge()
{
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceData, sizeof(deviceSU3) * m_uiLinkeCount));
}

CFieldGaugeU1::~CFieldGaugeU1()
{
    checkCudaErrors(__cudaFree(m_pDeviceData));
}

void CFieldGaugeU1::ExpMult(Real a, CField* U) const
{
    if (NULL == U || EFT_GaugeU1 != U->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: U field is not SU3");
        return;
    }

    CFieldGaugeU1* pUField = dynamic_cast<CFieldGaugeU1*>(U);

    preparethread;
    _kernelExpMultU1Real << < block, threads >> > (m_pDeviceData, a, pUField->m_pDeviceData);
    
}

void CFieldGaugeU1::ElementNormalize()
{
    preparethread;
    _kernelNormalizeU1 << < block, threads >> > (m_pDeviceData);
}

#if !_CLG_DOUBLEFLOAT
cuDoubleComplex CFieldGaugeU1::Dot(const CField* other) const
#else
CLGComplex CFieldGaugeU1::Dot(const CField* other) const
#endif
{
    if (NULL == other || EFT_GaugeU1 != other->GetFieldType())
    {
        appCrucial("CFieldGaugeU1: U field is not SU3");
        return make_cuDoubleComplex(0,0);
    }

    const CFieldGaugeU1* pUField = dynamic_cast<const CFieldGaugeU1*>(other);

    preparethread;
    _kernelDotU1 << < block, threads >> > (m_pDeviceData, pUField->m_pDeviceData, _D_ComplexThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

void CFieldGaugeU1::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || EFT_GaugeU1 != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: target field is not SU3");
        return;
    }

    CFieldGauge::CopyTo(pTarget);

    CFieldGaugeU1* pTargetField = dynamic_cast<CFieldGaugeU1*>(pTarget);
    checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceData, m_pDeviceData, sizeof(CLGComplex) * m_uiLinkeCount, cudaMemcpyDeviceToDevice));
}

void CFieldGaugeU1::TransformToIA()
{
    preparethread;
    _kernelTransformToIAU1 << <block, threads >> > (m_pDeviceData);
}


void CFieldGaugeU1::TransformToU()
{
    preparethread;
    _kernelTransformToU_U1 << <block, threads >> > (m_pDeviceData);
}

void CFieldGaugeU1::CalculateE_Using_U(CFieldGauge* pResoult) const
{
    appCrucial(_T("U1 CalculateE_Using_U Not supported!\n"));
    _FAIL_EXIT;

    //if (NULL == pResoult || EFT_GaugeSU3 != pResoult->GetFieldType())
    //{
    //    appCrucial("CFieldGaugeSU3: U field is not SU3");
    //    return;
    //}

    //CFieldGaugeU1* pUField = dynamic_cast<CFieldGaugeU1*>(pResoult);

    //preparethread;
    //_kernelTransformToE << <block, threads >> > (pUField->m_byFieldId, m_pDeviceData, pUField->m_pDeviceData);
}

void CFieldGaugeU1::CalculateNablaE_Using_U(CFieldGauge* pResoult) const
{
    appCrucial(_T("U1 CalculateNablaE_Using_U Not supported!\n"));
    _FAIL_EXIT;

    //if (NULL == pResoult || EFT_GaugeSU3 != pResoult->GetFieldType())
    //{
    //    appCrucial("CFieldGaugeSU3: U field is not SU3");
    //    return;
    //}

    //CFieldGaugeU1* pUField = dynamic_cast<CFieldGaugeU1*>(pResoult);

    //preparethread;
    //_kernelCalculateNablaE << <block, threads >> > (
    //    m_pDeviceData,
    //    m_byFieldId,
    //    pUField->m_pDeviceData);
}

void CFieldGaugeU1::DebugPrintMe() const
{
    //preparethread;
    //_kernelPrintSU3 << < block, threads >> > (m_pDeviceData);

    //===================================================
    //Since Debug Print Me is only used to debug, we do it slow but convinient
    CLGComplex* pToPrint = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(pToPrint, m_pDeviceData, sizeof(CLGComplex) * m_uiLinkeCount, cudaMemcpyDeviceToHost));

    for (UINT uiSite = 0; uiSite < m_uiLinkeCount / _HC_Dir; ++uiSite)
    {
        appGeneral(_T(" --- site: %d --- "), uiSite);
        for (UINT uiDir = 0; uiDir < _HC_Dir; ++uiDir)
        {
            UINT uiLink = uiSite * _HC_Dir + uiDir;
            appGeneral(_T(" %f %s %f I, "),
                pToPrint[uiLink].x,
                pToPrint[uiLink].y > F(0.0) ? _T("+") : _T("-"),
                appAbs(pToPrint[uiLink].y)
                );
        }
    }

    free(pToPrint);
}

void CFieldGaugeU1::SaveToFile(const CCString &fileName) const
{
    UINT uiSize = 0;
    BYTE* byToSave = CopyDataOut(uiSize);
    appGetFileSystem()->WriteAllBytes(fileName.c_str(), byToSave, uiSize);
    free(byToSave);
}

void CFieldGaugeU1::SaveToCompressedFile(const CCString& fileName) const
{
    appCrucial(_T("U1 SaveToCompressedFile Not supported!\n"));
    _FAIL_EXIT;
}

BYTE* CFieldGaugeU1::CopyDataOut(UINT &uiSize) const
{
    CLGComplex* toSave = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(CLGComplex) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount * 2);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[2];
        for (UINT j = 0; j < 9; ++j)
        {
            oneLink[2 * j] = static_cast<Real>(toSave[i].x);
            oneLink[2 * j + 1] = static_cast<Real>(toSave[i].y);
        }
        memcpy(byToSave + i * sizeof(Real) * 2, oneLink, sizeof(Real) * 2);
    }
    free(toSave);

    return byToSave;
}

CCString CFieldGaugeU1::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CFieldGaugeU1\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================