//=============================================================================
// FILENAME : CFieldGaugeU1Real.cu
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
#include "CFieldGaugeU1.h"
#include "CFieldGaugeU1Real.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldGaugeU1Real)

#pragma region Kernels

/**
* Initial SU3 Field with a value
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialU1RealField(Real *pDevicePtr, EFieldInitialType eInitialType)
{
    intokernaldir;
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        switch (eInitialType)
        {
        case EFIT_Zero:
        {
            pDevicePtr[uiLinkIndex] = F(0.0);
        }
        break;
        case EFIT_Identity:
        {
            pDevicePtr[uiLinkIndex] = F(1.0);
        }
        break;
        case EFIT_Random:
        {
            pDevicePtr[uiLinkIndex] = _deviceRandomF(_deviceGetFatIndex(uiSiteIndex, idir + 1)) * PI2 - PI;
        }
        break;
        case EFIT_RandomGenerator:
        {
            pDevicePtr[uiLinkIndex] = _deviceRandomGaussFSqrt2(_deviceGetFatIndex(uiSiteIndex, idir + 1)) * PI2 - PI;
        }
        break;
        default:
        {
            printf("U1 Field cannot be initialized with this type!");
        }
        break;
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDaggerU1Real(Real* pDevicePtr)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex] = -pDevicePtr[uiLinkIndex];

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyU1Real_R(Real*pDevicePtr, const Real* __restrict__ x, Real a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex] = pDevicePtr[uiLinkIndex] + x[uiLinkIndex] * a;

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusU1Real(Real*pDevicePtr, const Real* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex] = pDevicePtr[uiLinkIndex] + x[uiLinkIndex];

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusU1Real(Real*pDevicePtr, const Real* __restrict__ x)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex] = pDevicePtr[uiLinkIndex] - x[uiLinkIndex];

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyU1Real_R(Real*pDevicePtr, Real a)
{
    gaugeSU3KernelFuncionStart

    pDevicePtr[uiLinkIndex] = pDevicePtr[uiLinkIndex] * a;

    gaugeSU3KernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteU1RealCacheIndex(
    const Real * __restrict__ pDeviceData,
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
            Real toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

            if (first.NeedToDagger())
            {
                toAdd = -toAdd;
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                Real toMul = pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)];


                if (nextlink.NeedToDagger())
                {
                    toAdd = toAdd - toMul;
                }
                else
                {
                    toAdd = toAdd + toMul;
                }
            }
            res = _cuCaddf(res, _make_cuComplex(_cos(toAdd), _sin(toAdd)));
        }
        if (NULL != pStapleData)
        {
            pStapleData[linkIndex] = res;
        }

        //staple calculated
        //CLGComplex force(pDeviceData[linkIndex]);
        //force = _cuCmulf(_make_cuComplex(res.x, -res.y));
        //force.Ta();
        //force.MulReal(betaOverN);
        const Real thisLinkArg = pDeviceData[linkIndex];
        const CLGComplex thisLink = _make_cuComplex(_cos(thisLinkArg), _sin(thisLinkArg));
        CLGComplex force = _make_cuComplex(F(0.0), betaOverN * 
            (thisLink.y * res.x - thisLink.x * res.y));
        

        //force is additive
        pForceData[linkIndex] = _cuCaddf(pForceData[linkIndex], force);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateOnlyStapleU1Real(
    const Real* __restrict__ pDeviceData,
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
            Real toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

            if (first.NeedToDagger())
            {
                toAdd = -toAdd;
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                Real toMul = pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)];

                if (nextlink.NeedToDagger())
                {
                    toAdd = toAdd - toMul;
                }
                else
                {
                    toAdd = toAdd + toMul;
                }
            }
            res = _cuCaddf(res, _make_cuComplex(_cos(toAdd), _sin(toAdd)));
        }
        if (NULL != pStapleData)
        {
            pStapleData[linkIndex] = res;
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyU1RealCacheIndex(
    const Real* __restrict__ pDeviceData,
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
        Real toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];
        if (first.NeedToDagger())
        {
            toAdd = -toAdd;
        }

        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedIndex[i * plaqLength + j + uiSiteIndex * plaqCountAll];

            Real toMul = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];
            if (first.NeedToDagger())
            {
                toAdd = toAdd - toMul;
            }
            else
            {
                toAdd = toAdd + toMul;
            }
        }

#if !_CLG_DOUBLEFLOAT
        resThisThread += (1.0 - cos(static_cast<DOUBLE>(toAdd)));
#else
        resThisThread += (F(1.0) - _cos(toAdd));
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
_kernelPlaqutteEnergyUsingStableU1Real(
    const Real* __restrict__ pDeviceData,
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
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        const DOUBLE thisLinkArg = static_cast<DOUBLE>(pDeviceData[linkIndex]);
        //For each link, there are 6 staples
        resThisThread += (6.0 -(cos(thisLinkArg) * pStableData[linkIndex].x + sin(thisLinkArg) * pStableData[linkIndex].y));
    }

    results[uiSiteIndex] = resThisThread * betaOverN * 0.25;
#else
    Real resThisThread = F(0.0);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        const Real thisLinkArg = pDeviceData[linkIndex];
        //For each link, there are 6 staples
        resThisThread += (F(6.0) - (_cos(thisLinkArg) * pStableData[linkIndex].x + _sin(thisLinkArg) * pStableData[linkIndex].y));
    }

    results[uiSiteIndex] = resThisThread * betaOverN * F(0.25);
#endif
    //printf("  ---- energy: thread=%d, res=%f\n", __thread_id, results[__thread_id]);
}

/**
 * What is this for?
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelExpMultU1Real_R(
    const Real* __restrict__ pMyDeviceData,
    Real a,
    CLGComplex*pU)
{
    intokernaldir;

    for (BYTE idir = 0; idir < uiDir; ++idir)
    {
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        const Real fArg = pMyDeviceData[linkIndex];
        pU[linkIndex] = _cuCmulf(__cuCexpf(_make_cuComplex(_cos(fArg) * a, _sin(fArg) * a)), pU[linkIndex]);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDotU1RealComplex(
    const Real* __restrict__ pMyDeviceData,
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
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        const DOUBLE fArg = static_cast<DOUBLE>(pMyDeviceData[linkIndex]);
        const DOUBLE cs = cos(fArg);
        const DOUBLE sn = sin(fArg);
        resThisThread = cuCadd(resThisThread,
            make_cuDoubleComplex(cs * pOtherDeviceData[linkIndex].x + sn * pOtherDeviceData[linkIndex].y,
                cs * pOtherDeviceData[linkIndex].y - sn * pOtherDeviceData[linkIndex].x)
        );
    }

    result[uiSiteIndex] = resThisThread;
#else
    CLGComplex resThisThread = _make_cuComplex(0,0);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        const Real fArg = pMyDeviceData[linkIndex];
        const Real cs = _cos(fArg);
        const Real sn = _sin(fArg);
        resThisThread = _cuCaddf(resThisThread, 
            _make_cuComplex(cs * pOtherDeviceData[linkIndex].x + sn * pOtherDeviceData[linkIndex].y,
                cs * pOtherDeviceData[linkIndex].y - sn * pOtherDeviceData[linkIndex].x)
        );
    }

    result[uiSiteIndex] = resThisThread;
#endif
    //printf("res = %f %f\n", pOtherDeviceData[uiSiteIndex * 4].m_me[0].x, pMyDeviceData[uiSiteIndex * 4].m_me[0].x);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDotU1RealReal(
    const Real* __restrict__ pMyDeviceData,
    const Real* __restrict__ pOtherDeviceData,
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
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        const DOUBLE fArg1 = static_cast<DOUBLE>(pMyDeviceData[linkIndex]);
        const DOUBLE cs1 = cos(fArg1);
        const DOUBLE sn1 = sin(fArg1);
        const DOUBLE fArg2 = static_cast<DOUBLE>(pOtherDeviceData[linkIndex]);
        const DOUBLE cs2 = cos(fArg2);
        const DOUBLE sn2 = sin(fArg2);
        resThisThread = cuCadd(resThisThread,
            make_cuDoubleComplex(cs1 * cs2 + sn1 * sn2,
                cs1 * sn2 - sn1 * cs2)
        );
    }

    result[uiSiteIndex] = resThisThread;
#else
    CLGComplex resThisThread = _make_cuComplex(0, 0);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        const Real fArg1 = pMyDeviceData[linkIndex];
        const Real cs1 = _cos(fArg1);
        const Real sn1 = _sin(fArg1);
        const Real fArg2 = pOtherDeviceData[linkIndex];
        const Real cs2 = _cos(fArg2);
        const Real sn2 = _sin(fArg2);
        resThisThread = _cuCaddf(resThisThread,
            _make_cuComplex(cs1 * cs2 + sn1 * sn2,
                cs1 * sn2 - sn1 * cs2)
        );
    }

    result[uiSiteIndex] = resThisThread;
#endif
    //printf("res = %f %f\n", pOtherDeviceData[uiSiteIndex * 4].m_me[0].x, pMyDeviceData[uiSiteIndex * 4].m_me[0].x);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelSetConfigurationU1Real(
    Real* pDeviceData,
    const Real* __restrict__ pRealData)
{
    gaugeSU3KernelFuncionStart

    pDeviceData[uiLinkIndex] = pRealData[uiLinkIndex];

    gaugeSU3KernelFuncionEnd
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


__global__ void _CLG_LAUNCH_BOUND
_kernelSetOneDirUnity_U1Real(Real* pDeviceData, BYTE byDir)
{
    intokernaldir;

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = F(0.0);
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelSetOneDirUnityPoint_U1Real(Real* pDeviceData, UINT uiSiteIndex, BYTE byDir)
{
    for (BYTE dir = 0; dir < 4; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = F(0.0);
        }
    }
}

/**
 * At(n) = constant
 * Az(Lz) = Lz Ez t
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialAsImagineChemical(Real* pDeviceData, Real fValue)
{
    intokernal;
    pDeviceData[_deviceGetLinkIndex(uiSiteIndex, 3)] = fValue;
}

/**
 * At(n) = - Ez z
 * Az(Lz) = Lz Ez t
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialAsEz_Type0(Real* pDeviceData, Real fEz)
{
    intokernalInt4;
    const UINT uiLinkT = _deviceGetLinkIndex(uiSiteIndex, 3);
    const Real fZ = sSite4.z - _DC_Centerz;
    pDeviceData[uiLinkT] = -fEz * fZ;

    if (sSite4.z == _DC_Lz - 1)
    {
        const UINT uiLinkZ = _deviceGetLinkIndex(uiSiteIndex, 2);
        const Real fT = sSite4.w - _DC_Centert;
        pDeviceData[uiLinkZ] = _DC_Lz * fEz * fT;
    }
}

/**
 * Az(n) = Ez t
 * At(Lt) = -Lt Ez z
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialAsEz_Type1(Real* pDeviceData, Real fEz)
{
    intokernalInt4;
    const UINT uiLinkZ = _deviceGetLinkIndex(uiSiteIndex, 2);
    const Real fT = sSite4.w - _DC_Centert;
    pDeviceData[uiLinkZ] = fEz * fT;

    if (sSite4.w == _DC_Lt - 1)
    {
        const UINT uiLinkT = _deviceGetLinkIndex(uiSiteIndex, 3);
        const Real fZ = sSite4.z - _DC_Centerz;
        pDeviceData[uiLinkT] = -_DC_Lti * fEz * fZ;
    }
}

/**
 * Projective plane, magnetic field
 *
 * Ay(n) = Bz nx
 * if twisted Ax(Lx) = - q B ny
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialAsBz_Type0(Real* pDeviceData, Real fBz, UBOOL bTwisted, UBOOL bProjectivePlane, UBOOL bShiftCenter)
{
    intokernalInt4;
    const UINT uiLinkY = _deviceGetLinkIndex(uiSiteIndex, 1);
    const Real fX = sSite4.x - _DC_Centerx + (bShiftCenter ? F(0.5) : F(0.0));
    pDeviceData[uiLinkY] = fBz * fX;

    if (bProjectivePlane)
    {
        const UINT uiLinkX = _deviceGetLinkIndex(uiSiteIndex, 0);
        const Real fY = sSite4.y - _DC_Centery + (bShiftCenter ? F(0.5) : F(0.0));
        if (bTwisted && sSite4.x == _DC_Lx - 1)
        {
            pDeviceData[uiLinkX] = -fBz * fY;
        }
        if (bTwisted && sSite4.y == _DC_Ly - 1)
        {
            pDeviceData[uiLinkY] = fBz * fX;
        }
    }
    else
    {
        if (bTwisted && sSite4.x == _DC_Lx - 1)
        {
            const UINT uiLinkX = _deviceGetLinkIndex(uiSiteIndex, 0);
            const Real fY = sSite4.y - _DC_Centery + (bShiftCenter ? F(0.5) : F(0.0));
            pDeviceData[uiLinkX] = -fBz * fY * _DC_Lxi;
        }
    }
}

/**
 * Projective plane, magnetic field
 *
 * Ax(n) = - Bz ny
 * if twisted Ay(Ly) = q B nx
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialAsBz_Type1(Real* pDeviceData, Real fBz, UBOOL bTwisted, UBOOL bProjectivePlane, UBOOL bShiftCenter)
{
    intokernalInt4;
    const UINT uiLinkX = _deviceGetLinkIndex(uiSiteIndex, 0);
    const Real fY = sSite4.y - _DC_Centery + (bShiftCenter ? F(0.5) : F(0.0));
    pDeviceData[uiLinkX] = -fBz * fY;

    if (bProjectivePlane)
    {
        const UINT uiLinkY = _deviceGetLinkIndex(uiSiteIndex, 1);
        const Real fX = sSite4.x - _DC_Centerx + (bShiftCenter ? F(0.5) : F(0.0));
        if (bTwisted && sSite4.x == _DC_Lx - 1)
        {
            pDeviceData[uiLinkX] = -fBz * fY;
        }
        if (bTwisted && sSite4.y == _DC_Ly - 1)
        {
            pDeviceData[uiLinkY] = fBz * fX;
        }
    }
    else
    {
        if (bTwisted && sSite4.y == _DC_Ly - 1)
        {
            const UINT uiLinkY = _deviceGetLinkIndex(uiSiteIndex, 1);
            const Real fX = sSite4.x - _DC_Centerx + (bShiftCenter ? F(0.5) : F(0.0));
            pDeviceData[uiLinkY] = fBz * fX * _DC_Lyi;
        }
    }
}

/**
 * Projective plane, magnetic field
 *
 * Ax(n) = - q B ny/2
 * Ay(n) = qB nx / 2
 *
 * if twisted
 * Ax(Lx) = - q B ny
 * Ay(Ly) = q B nx
 * 
 * Torus:
 * 
 * Link x: -y*Bz/2
 * Link y:  x*Bz/2
 * Edge x: -(Lx+1)*y*Bz/2
 * Edge y:  (Ly+1)*x*Bz/2
 * Corner: -(Lx*Ly-1)*Bz
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialAsBz_Type2(Real* pDeviceData, Real fBz, UBOOL bTwisted, UBOOL bProjectivePlane, UBOOL bShiftCenter)
{
    intokernalInt4;

    const UINT uiLinkX = _deviceGetLinkIndex(uiSiteIndex, 0);
    const UINT uiLinkY = _deviceGetLinkIndex(uiSiteIndex, 1);
    const Real fX = sSite4.x - _DC_Centerx + (bShiftCenter ? F(0.5) : F(0.0));
    const Real fY = sSite4.y - _DC_Centery + (bShiftCenter ? F(0.5) : F(0.0));
    pDeviceData[uiLinkX] = -fBz * fY * F(0.5);
    pDeviceData[uiLinkY] = fBz * fX * F(0.5);

    if (bProjectivePlane)
    {
        if (bTwisted && sSite4.x == _DC_Lx - 1)
        {
            pDeviceData[uiLinkX] = -fBz * fY;
        }
        if (bTwisted && sSite4.y == _DC_Ly - 1)
        {
            pDeviceData[uiLinkY] = fBz * fX;
        }
    }
    else
    {
        if (bTwisted && sSite4.x == _DC_Lx - 1)
        {
            pDeviceData[uiLinkX] = -fBz * fY * (_DC_Lx + 1) * F(0.5);
        }
        if (bTwisted && sSite4.y == _DC_Ly - 1)
        {
            pDeviceData[uiLinkY] = fBz * fX * (_DC_Ly + 1) * F(0.5);
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovLoopOfSiteU1Real(
    const Real *__restrict__ pDeviceBuffer,
    cuDoubleComplex* res)
{
    UINT uiXYZ = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    const UINT uiSiteIndex = uiXYZ * _DC_Lt;
    UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, _DC_Dir - 1);
    SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

    Real tmp = F(0.0);
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
            tmp = tmp + pDeviceBuffer[uiLinkIdx];
        }
    }

    res[uiXYZ] = make_cuDoubleComplex(_cos(tmp), _sin(tmp));
}

#pragma endregion

void CFieldGaugeU1Real::AxpyPlus(const CField* x)
{
    if (NULL == x || EFT_GaugeReal != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeU1Real* pSU3x = dynamic_cast<const CFieldGaugeU1Real*>(x);
    preparethread;
    _kernelAxpyPlusU1Real << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData);
}

void CFieldGaugeU1Real::AxpyMinus(const CField* x)
{
    if (NULL == x || EFT_GaugeReal != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeU1Real* pSU3x = dynamic_cast<const CFieldGaugeU1Real*>(x);
    preparethread;
    _kernelAxpyMinusU1Real << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData);

}

void CFieldGaugeU1Real::ScalarMultply(const CLGComplex& a)
{
    appCrucial(_T("CFieldGaugeU1Real::ScalarMultply with complex Not supported yet!\n"));
    _FAIL_EXIT;
}

void CFieldGaugeU1Real::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplyU1Real_R << <block, threads >> > (m_pDeviceData, a);
}

void CFieldGaugeU1Real::Axpy(Real a, const CField* x)
{
    if (NULL == x || EFT_GaugeReal != x->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeU1Real* pSU3x = dynamic_cast<const CFieldGaugeU1Real*>(x);
    preparethread;
    _kernelAxpyU1Real_R << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData, a);
}

void CFieldGaugeU1Real::Axpy(const CLGComplex& a, const CField* x)
{
    appCrucial(_T("CFieldGaugeU1Real::Axpy with complex Not supported yet!\n"));
    _FAIL_EXIT;
}

void CFieldGaugeU1Real::Mul(const CField* other, UBOOL bDagger)
{
    appCrucial(_T("CFieldGaugeU1Real::Mul with complex Not supported yet!\n"));
    _FAIL_EXIT;
}

void CFieldGaugeU1Real::Zero()
{
    appCrucial(_T("CFieldGaugeU1Real::Zero Not supported yet!\n"));
    _FAIL_EXIT;
}

void CFieldGaugeU1Real::Identity()
{
    preparethread;
    _kernelInitialU1RealField << <block, threads >> > (m_pDeviceData, EFIT_Identity);
}

void CFieldGaugeU1Real::Dagger()
{
    preparethread;
    _kernelDaggerU1Real << <block, threads >> > (m_pDeviceData);
}

void CFieldGaugeU1Real::MakeRandomGenerator()
{
    preparethread;
    _kernelInitialU1RealField << <block, threads >> > (m_pDeviceData, EFIT_RandomGenerator);
}

void CFieldGaugeU1Real::InitialOtherParameters(CParameters& param)
{
    CFieldGauge::InitialOtherParameters(param);

    if (EFIT_U1Real == m_eInitialType)
    {
        CCString sType;
        Real fChemical = F(0.0);
        Real fEz = F(0.0);
        Real fBz = F(0.0);

        EU1RealType eChemical = EURT_None;
        EU1RealType eEz = EURT_None;
        EU1RealType eBz = EURT_None;

        sType = _T("EURT_None");

        if (param.FetchStringValue(_T("ChemicalType"), sType))
        {
            eChemical = __STRING_TO_ENUM(EU1RealType, sType);
            if (EURT_None != eChemical)
            {
                param.FetchValueReal(_T("ChemicalValue"), fChemical);
            }
        }

        sType = _T("EURT_None");
        if (param.FetchStringValue(_T("EzType"), sType))
        {
            eEz = __STRING_TO_ENUM(EU1RealType, sType);
            if (EURT_None != eEz)
            {
                param.FetchValueReal(_T("EzValue"), fEz);
            }
        }

        sType = _T("EURT_None");
        if (param.FetchStringValue(_T("BzType"), sType))
        {
            eBz = __STRING_TO_ENUM(EU1RealType, sType);
            if (EURT_None != eBz)
            {
                param.FetchValueReal(_T("BzValue"), fBz);
            }
        }

        INT iValue = 1;
        param.FetchValueINT(_T("XYShiftCenter"), iValue);
        UBOOL bXYShiftCenter = (0 != iValue);

        InitialU1Real(eChemical, eEz, eBz, fChemical, fEz, fBz, bXYShiftCenter);
    }
}

/**
*
*/
void CFieldGaugeU1Real::InitialField(EFieldInitialType eInitialType)
{
    if (EFIT_U1Real == eInitialType)
    {
        m_eInitialType = EFIT_U1Real;
        return;
    }

    preparethread;
    _kernelInitialU1RealField << <block, threads >> > (m_pDeviceData, eInitialType);
}

void CFieldGaugeU1Real::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eType)
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
        UINT uiSize = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount);
        BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        InitialWithByte(data);
        free(data);
    }
    break;
#if _CLG_DOUBLEFLOAT
    case EFFT_CLGBinFloat:
    {
        UINT uiSize = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount);
        BYTE* data = (BYTE*)malloc(uiSize);
        Real* rdata = (Real*)data;
        FLOAT* fdata = (FLOAT*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        for (UINT i = 0; i < m_uiLinkeCount; ++i)
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
        UINT uiSize = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount);
        BYTE* data = (BYTE*)malloc(uiSize);
        Real* rdata = (Real*)data;
        DOUBLE* ddata = (DOUBLE*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        for (UINT i = 0; i < m_uiLinkeCount; ++i)
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

void CFieldGaugeU1Real::InitialWithByte(BYTE* byData)
{
    Real* readData = (Real*)malloc(sizeof(Real) * m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[1];
        memcpy(oneLink, byData + sizeof(Real) * i, sizeof(Real));
        readData[i] = oneLink[i];
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(Real) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    free(readData);
}

void CFieldGaugeU1Real::InitialU1Real(EU1RealType eChemicalType, EU1RealType eEType, EU1RealType eBType, Real fChemical, Real feEz, Real feBz, UBOOL bXYShiftCenter)
{
    m_eChemical = eChemicalType;
    m_eE = eEType;
    m_eB = eBType;
    m_fChemical = fChemical;
    m_feEz = feEz;
    m_feBz = feBz;
    m_bXYShiftCenter = bXYShiftCenter;

    UBOOL bProjective = (NULL != dynamic_cast<const CBoundaryConditionProjectivePlaneSquare*>(appGetLattice()->m_pIndex->GetBoudanryCondition()));

    preparethread;
    _kernelInitialU1RealField << <block, threads >> > (m_pDeviceData, EFIT_Zero);

    switch (eChemicalType)
    {
    case EURT_ImagineChemical:
        _kernelInitialAsImagineChemical << <block, threads >> > (m_pDeviceData, fChemical);
        break;
    default:
        appParanoiac(_T("No chemical set\n"));
        break;
    }

    switch (eEType)
    {
    case EURT_E_t:
        _kernelInitialAsEz_Type0 << <block, threads >> > (m_pDeviceData, feEz);
        break;
    case EURT_E_z:
        _kernelInitialAsEz_Type1 << <block, threads >> > (m_pDeviceData, feEz);
        break;
    default:
        appParanoiac(_T("No electric set\n"));
        break;
    }

    switch (eBType)
    {
    case EURT_Bp_x:
        _kernelInitialAsBz_Type0 << <block, threads >> > (m_pDeviceData, feBz, TRUE, bProjective, bXYShiftCenter);
        break;
    case EURT_Bp_y:
        _kernelInitialAsBz_Type1 << <block, threads >> > (m_pDeviceData, feBz, TRUE, bProjective, bXYShiftCenter);
        break;
    case EURT_Bp_xy:
        _kernelInitialAsBz_Type2 << <block, threads >> > (m_pDeviceData, feBz, TRUE, bProjective, bXYShiftCenter);
        break;
    case EURT_Bp_x_notwist:
        _kernelInitialAsBz_Type0 << <block, threads >> > (m_pDeviceData, feBz, FALSE, bProjective, bXYShiftCenter);
        break;
    case EURT_Bp_y_notwist:
        _kernelInitialAsBz_Type1 << <block, threads >> > (m_pDeviceData, feBz, FALSE, bProjective, bXYShiftCenter);
        break;
    case EURT_Bp_xy_notwist:
        _kernelInitialAsBz_Type2 << <block, threads >> > (m_pDeviceData, feBz, FALSE, bProjective, bXYShiftCenter);
        break;
    default:
        appParanoiac(_T("No magnetic set\n"));
        break;
    }

}

void CFieldGaugeU1Real::InitialWithByteCompressed(const CCString& sFileName)
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

void CFieldGaugeU1Real::SetByArray(Real* array)
{
    assert(NULL != array);
    //we algin the su3 now
    //assert(sizeof(deviceSU3) == 32 * sizeof(Real));

    //checkCudaErrors(cudaMemcpy(m_pDeviceData, array, sizeof(Real) * _HC_LinkCount * 18, cudaMemcpyHostToDevice));
    
    Real* pDeviceArray;
    checkCudaErrors(__cudaMalloc((void**)&pDeviceArray, sizeof(Real) * _HC_LinkCount));
    checkCudaErrors(cudaMemcpy(pDeviceArray, array, sizeof(Real) * _HC_LinkCount, cudaMemcpyHostToDevice));
    preparethread;
    _kernelSetConfigurationU1Real << <block, threads >> > (m_pDeviceData, pDeviceArray);
    checkCudaErrors(__cudaFree(pDeviceArray));

    free(array);
}

/**
* (1) calculate staples
* (2) calculate force(additive)
* (3) calculate energy
*/
void CFieldGaugeU1Real::CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, Real betaOverN) const
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

    _kernelStapleAtSiteU1RealCacheIndex << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData,
        betaOverN);
}

void CFieldGaugeU1Real::CalculateOnlyStaple(CFieldGauge* pStable) const
{
    if (NULL == pStable || EFT_GaugeU1 != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stable field is not SU3");
        return;
    }
    CFieldGaugeU1* pStableSU3 = dynamic_cast<CFieldGaugeU1*>(pStable);

    preparethread;
    _kernelCalculateOnlyStapleU1Real << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStableSU3->m_pDeviceData);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeU1Real::CalculatePlaqutteEnergy(DOUBLE betaOverN) const
#else
Real CFieldGaugeU1Real::CalculatePlaqutteEnergy(Real betaOverN) const
#endif
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);

    preparethread;
    _kernelPlaqutteEnergyU1RealCacheIndex << <block, threads >> > (
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
DOUBLE CFieldGaugeU1Real::CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const
#else
Real CFieldGaugeU1Real::CalculatePlaqutteEnergyUseClover(Real betaOverN) const
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
DOUBLE CFieldGaugeU1Real::CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge *pStable) const
#else
Real CFieldGaugeU1Real::CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge* pStable) const
#endif
{
    if (NULL == pStable || EFT_GaugeU1 != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: stape field is not SU3");
        return F(0.0);
    }
    const CFieldGaugeU1* pStableSU3 = dynamic_cast<const CFieldGaugeU1*>(pStable);

    preparethread;
    _kernelPlaqutteEnergyUsingStableU1Real << <block, threads >> > (
        m_pDeviceData, 
        pStableSU3->m_pDeviceData, 
        betaOverN, 
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CFieldGaugeU1Real::CalculateKinematicEnergy() const
#else
Real CFieldGaugeU1Real::CalculateKinematicEnergy() const
#endif
{
    appCrucial(_T("U1Real CalculateKinematicEnergy not supported\n"));
#if !_CLG_DOUBLEFLOAT
    return 0.0;
#else
    return F(0.0);
#endif
}

void CFieldGaugeU1Real::SetOneDirectionUnity(BYTE byDir)
{
    if (0 == (byDir & 15))
    {
        return;
    }
    preparethread;
    _kernelSetOneDirUnity_U1Real << <block, threads >> >(m_pDeviceData, byDir);

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

void CFieldGaugeU1Real::SetOneDirectionZero(BYTE byDir)
{
    if (0 == (byDir & 15))
    {
        return;
    }
    appCrucial(_T("U1Real Zero not supported\n"));
}

CFieldGaugeU1Real::CFieldGaugeU1Real()
    : CFieldGauge()
    , m_eInitialType(EFIT_Random)
    , m_eChemical(EURT_None)
    , m_eE(EURT_None)
    , m_eB(EURT_None)
    , m_fChemical(F(0.0))
    , m_feEz(F(0.0))
    , m_feBz(F(0.0))
    , m_bXYShiftCenter(TRUE)
{
    checkCudaErrors(__cudaMalloc((void **)&m_pDeviceData, sizeof(Real) * m_uiLinkeCount));
}

CFieldGaugeU1Real::~CFieldGaugeU1Real()
{
    checkCudaErrors(__cudaFree(m_pDeviceData));
}

void CFieldGaugeU1Real::ExpMult(Real a, CField* U) const
{
    if (NULL == U || EFT_GaugeReal != U->GetFieldType())
    {
        appCrucial("ExpMult: U field is not EFT_GaugeReal");
        return;
    }

    CFieldGaugeU1* pUField = dynamic_cast<CFieldGaugeU1*>(U);

    preparethread;
    _kernelExpMultU1Real_R << < block, threads >> > (m_pDeviceData, a, pUField->m_pDeviceData);
    
}

#if !_CLG_DOUBLEFLOAT
cuDoubleComplex CFieldGaugeU1Real::Dot(const CField* other) const
#else
CLGComplex CFieldGaugeU1Real::Dot(const CField* other) const
#endif
{
    if (NULL == other || (EFT_GaugeU1 != other->GetFieldType() && EFT_GaugeReal != other->GetFieldType()))
    {
        appCrucial("CFieldGaugeU1: U field is not SU3");
        return make_cuDoubleComplex(0,0);
    }

    if (EFT_GaugeU1 == other->GetFieldType())
    {
        const CFieldGaugeU1* pUField = dynamic_cast<const CFieldGaugeU1*>(other);

        preparethread;
        _kernelDotU1RealComplex << < block, threads >> > (m_pDeviceData, pUField->m_pDeviceData, _D_ComplexThreadBuffer);
        return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
    }

    const CFieldGaugeU1Real* pUField = dynamic_cast<const CFieldGaugeU1Real*>(other);

    preparethread;
    _kernelDotU1RealReal << < block, threads >> > (m_pDeviceData, pUField->m_pDeviceData, _D_ComplexThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

void CFieldGaugeU1Real::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || EFT_GaugeReal != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeSU3: target field is not EFT_GaugeReal");
        return;
    }

    CFieldGauge::CopyTo(pTarget);

    CFieldGaugeU1Real* pTargetField = dynamic_cast<CFieldGaugeU1Real*>(pTarget);
    checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceData, m_pDeviceData, sizeof(Real) * m_uiLinkeCount, cudaMemcpyDeviceToDevice));
}

void CFieldGaugeU1Real::TransformToIA()
{
    appCrucial(_T("U1Real itself is A, so no transform can be done\n"));
}


void CFieldGaugeU1Real::TransformToU()
{
    appCrucial(_T("U1Real itself is A, so no transform can be done\n"));
}

void CFieldGaugeU1Real::CalculateE_Using_U(CFieldGauge* pResoult) const
{
    appCrucial(_T("U1Real CalculateE_Using_U Not supported!\n"));
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

void CFieldGaugeU1Real::CalculateNablaE_Using_U(CFieldGauge* pResoult, UBOOL bNaive) const
{
    appCrucial(_T("U1Real CalculateNablaE_Using_U Not supported!\n"));
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

void CFieldGaugeU1Real::DebugPrintMe() const
{
    //preparethread;
    //_kernelPrintSU3 << < block, threads >> > (m_pDeviceData);

    //===================================================
    //Since Debug Print Me is only used to debug, we do it slow but convinient
    Real* pToPrint = (Real*)malloc(sizeof(Real) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(pToPrint, m_pDeviceData, sizeof(Real) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    appPushLogDate(FALSE);
    for (UINT uiSite = 0; uiSite < m_uiLinkeCount / _HC_Dir; ++uiSite)
    {
        SSmallInt4 site = __hostSiteIndexToInt4(uiSite);
        appGeneral(_T(" --- site: %d %d %d %d --- "), site.x, site.y, site.z, site.w);
        for (UINT uiDir = 0; uiDir < _HC_Dir; ++uiDir)
        {
            const UINT uiLink = uiSite * _HC_Dir + uiDir;
            appGeneral(_T(" %f, "), pToPrint[uiLink]);
        }
        appGeneral(_T("\n"));
    }
    appPopLogDate();
    free(pToPrint);
}

BYTE* CFieldGaugeU1Real::CopyDataOut(UINT &uiSize) const
{
    uiSize = sizeof(Real) * m_uiLinkeCount;
    BYTE* byToSave = (BYTE*)malloc(uiSize);
    checkCudaErrors(cudaMemcpy(byToSave, m_pDeviceData, uiSize, cudaMemcpyDeviceToHost));
    return byToSave;
}

BYTE* CFieldGaugeU1Real::CopyDataOutFloat(UINT& uiSize) const
{
    Real* toSave = (Real*)malloc(sizeof(Real) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(Real) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(FLOAT) * m_uiLinkeCount);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        FLOAT oneLink = static_cast<FLOAT>(toSave[i]);
        memcpy(byToSave + i * sizeof(FLOAT), &oneLink, sizeof(FLOAT));
    }
    free(toSave);

    return byToSave;
}

BYTE* CFieldGaugeU1Real::CopyDataOutDouble(UINT& uiSize) const
{
    Real* toSave = (Real*)malloc(sizeof(Real) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(Real) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(DOUBLE) * m_uiLinkeCount);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        DOUBLE oneLink = static_cast<DOUBLE>(toSave[i]);
        memcpy(byToSave + i * sizeof(DOUBLE), &oneLink, sizeof(DOUBLE));
    }
    free(toSave);

    return byToSave;
}

CCString CFieldGaugeU1Real::GetInfos(const CCString &tab) const
{
    CCString sRet = CFieldGauge::GetInfos(tab);
    sRet = sRet + tab + _T("Initialed : ") + __ENUM_TO_STRING(EFieldInitialType, m_eInitialType) + _T("\n");

    sRet = sRet + tab + _T("Chemical : ") + __ENUM_TO_STRING(EU1RealType, m_eChemical) + _T(" , v = ") + appToString(m_fChemical) + _T("\n");
    sRet = sRet + tab + _T("Electric : ") + __ENUM_TO_STRING(EU1RealType, m_eE) + _T(" , v = ") + appToString(m_feEz) + _T("\n");
    sRet = sRet + tab + _T("Magnetic : ") + __ENUM_TO_STRING(EU1RealType, m_eB) + _T(" , v = ") + appToString(m_feBz) + _T("\n");
    sRet = sRet + tab + _T("XYShiftCenter : ") + appToString(m_bXYShiftCenter) + _T("\n");

    return sRet;
}

Real CFieldGaugeU1Real::CheckSliceSame(BYTE dir1, BYTE dir2) const
{
    Real* toCheck = (Real*)malloc(sizeof(Real) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toCheck, m_pDeviceData, sizeof(Real) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    TArray<INT> ext;
    ext.AddItem(_HC_Lxi);
    ext.AddItem(_HC_Lyi);
    ext.AddItem(_HC_Lzi);
    ext.AddItem(_HC_Lti);
    TArray<BYTE> otherdir;
    for (BYTE d = 0; d < 4; ++d)
    {
        if (d != dir1 && d != dir2)
        {
            otherdir.AddItem(d);
        }
    }

    Real fDelta = F(0.0);
    for (INT iInSliceX = 0; iInSliceX < ext[dir1]; ++iInSliceX)
    {
        for (INT iInSliceY = 0; iInSliceY < ext[dir2]; ++iInSliceY)
        {
            SSmallInt4 site1;
            for (SBYTE dr = 0; dr < 4; ++dr)
            {
                if (dr == dir1)
                {
                    site1.m_byData4[dr] = static_cast<SBYTE>(iInSliceX);
                }
                else if (dr == dir2)
                {
                    site1.m_byData4[dr] = static_cast<SBYTE>(iInSliceY);
                }
                else
                {
                    site1.m_byData4[dr] = 0;
                }
            }
            UINT uiSiteIndex = site1._hostToSiteIndex();
            Real lx1 = toCheck[uiSiteIndex * 4];
            Real ly1 = toCheck[uiSiteIndex * 4 + 1];
            Real lz1 = toCheck[uiSiteIndex * 4 + 2];
            Real lw1 = toCheck[uiSiteIndex * 4 + 3];
            for (INT iOtherX = 1; iOtherX < otherdir[0]; ++iOtherX)
            {
                for (INT iOtherY = 1; iOtherY < otherdir[1]; ++iOtherY)
                {
                    SSmallInt4 site2;
                    for (SBYTE dr2 = 0; dr2 < 4; ++dr2)
                    {
                        if (dr2 == dir1)
                        {
                            site2.m_byData4[dr2] = static_cast<SBYTE>(iInSliceX);
                        }
                        else if (dr2 == dir2)
                        {
                            site2.m_byData4[dr2] = static_cast<SBYTE>(iInSliceY);
                        }
                        else if (dr2 == otherdir[0])
                        {
                            site2.m_byData4[dr2] = static_cast<SBYTE>(iOtherX);
                        }
                        else if (dr2 == otherdir[1])
                        {
                            site2.m_byData4[dr2] = static_cast<SBYTE>(iOtherY);
                        }
                    }

                    UINT uiSiteIndex2 = site2._hostToSiteIndex();
                    Real lx2 = toCheck[uiSiteIndex2 * 4];
                    Real ly2 = toCheck[uiSiteIndex2 * 4 + 1];
                    Real lz2 = toCheck[uiSiteIndex2 * 4 + 2];
                    Real lw2 = toCheck[uiSiteIndex2 * 4 + 3];

                    fDelta = fDelta + (lx1 - lx2) * (lx1 - lx2) + (ly1 - ly2) * (ly1 - ly2) + (lz1 - lz2) * (lz1 - lz2) + (lw1 - lw2) * (lw1 - lw2);
                }
            }
        }
    }
    appSafeFree(toCheck);
    return fDelta;
}

Real CFieldGaugeU1Real::CheckZero(BYTE dir1, BYTE dir2, const TArray<BYTE>& linkdirs) const
{
    Real* toCheck = (Real*)malloc(sizeof(Real) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toCheck, m_pDeviceData, sizeof(Real) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    TArray<INT> ext;
    ext.AddItem(_HC_Lxi);
    ext.AddItem(_HC_Lyi);
    ext.AddItem(_HC_Lzi);
    ext.AddItem(_HC_Lti);

    Real fDelta = F(0.0);
    for (INT iInSliceX = 0; iInSliceX < ext[dir1]; ++iInSliceX)
    {
        for (INT iInSliceY = 0; iInSliceY < ext[dir2]; ++iInSliceY)
        {
            SSmallInt4 site1;
            for (SBYTE dr = 0; dr < 4; ++dr)
            {
                if (dr == dir1)
                {
                    site1.m_byData4[dr] = static_cast<SBYTE>(iInSliceX);
                }
                else if (dr == dir2)
                {
                    site1.m_byData4[dr] = static_cast<SBYTE>(iInSliceY);
                }
                else
                {
                    site1.m_byData4[dr] = 0;
                }
            }
            UINT uiSiteIndex = site1._hostToSiteIndex();
            for (INT ld = 0; ld < linkdirs.Num(); ++ld)
            {
                Real fV = toCheck[uiSiteIndex * 4 + linkdirs[ld]];
                fDelta = fDelta + fV * fV;
            }
        }
    }
    appSafeFree(toCheck);
    return fDelta;
}

void CFieldGaugeU1Real::DebugPrintSlice(BYTE dir1, BYTE dir2, const TArray<BYTE>& linkdirs) const
{
    appPushLogDate(FALSE);
    Real* toCheck = (Real*)malloc(sizeof(Real) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toCheck, m_pDeviceData, sizeof(Real) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    TArray<INT> ext;
    ext.AddItem(_HC_Lxi);
    ext.AddItem(_HC_Lyi);
    ext.AddItem(_HC_Lzi);
    ext.AddItem(_HC_Lti);

    for (INT iInSliceX = 0; iInSliceX < ext[dir1]; ++iInSliceX)
    {
        for (INT iInSliceY = 0; iInSliceY < ext[dir2]; ++iInSliceY)
        {
            SSmallInt4 site1;
            for (SBYTE dr = 0; dr < 4; ++dr)
            {
                if (dr == dir1)
                {
                    site1.m_byData4[dr] = static_cast<SBYTE>(iInSliceX);
                }
                else if (dr == dir2)
                {
                    site1.m_byData4[dr] = static_cast<SBYTE>(iInSliceY);
                }
                else
                {
                    site1.m_byData4[dr] = 0;
                }
            }
            UINT uiSiteIndex = site1._hostToSiteIndex();

            appGeneral(_T("- %s (%s): "), appToString(site1).c_str(), appToString(linkdirs).c_str());

            for (INT ld = 0; ld < linkdirs.Num(); ++ld)
            {
                appGeneral(_T("%.8f "), toCheck[uiSiteIndex * 4 + linkdirs[ld]]);
            }

            appGeneral(_T("\n"));
        }
    }

    appPopLogDate();
    appSafeFree(toCheck);
}

void CFieldGaugeU1Real::PolyakovOnSpatialSite(cuDoubleComplex* buffer) const
{
    dim3 block(_HC_DecompX, _HC_DecompY, 1);
    dim3 threads(_HC_DecompLx, _HC_DecompLy, 1);
    _kernelPolyakovLoopOfSiteU1Real << <block, threads >> > (m_pDeviceData, buffer);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================