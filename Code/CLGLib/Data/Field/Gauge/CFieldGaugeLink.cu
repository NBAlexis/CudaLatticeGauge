//=============================================================================
// FILENAME : CFieldGaugeLink<deviceGauge, matrixN>.cu
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
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldGaugeLink.h"

__BEGIN_NAMESPACE


#pragma region Kernels

/**
* Initial SU3 Field with a value
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialGaugeField(deviceGauge* pDevicePtr, EFieldInitialType eInitialType)
{
    deviceGauge id = _makeId<deviceGauge>();
    deviceGauge zero = _makeZero<deviceGauge>();

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
            pDevicePtr[uiLinkIndex] = _makeRandom<deviceGauge>(_deviceGetFatIndex(uiSiteIndex, idir + 1));
        }
        break;
        case EFIT_RandomGenerator:
        {
            pDevicePtr[uiLinkIndex] = _makeGaussian<deviceGauge>(_deviceGetFatIndex(uiSiteIndex, idir + 1));
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

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDaggerGauge(deviceGauge* pDevicePtr)
{
    gaugeLinkKernelFuncionStart

        _dagger(pDevicePtr[uiLinkIndex]);

    gaugeLinkKernelFuncionEnd
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyGaugeA(deviceGauge* pDevicePtr, const deviceGauge* __restrict__ x, CLGComplex a)
{
    gaugeLinkKernelFuncionStart

        _add(pDevicePtr[uiLinkIndex], _mulC(x[uiLinkIndex], a));

    gaugeLinkKernelFuncionEnd
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelMulGauge(deviceGauge* pDevicePtr, const deviceGauge* __restrict__ x, UBOOL bDagger)
{
    gaugeLinkKernelFuncionStart

        if (bDagger)
        {
            _dagmul(pDevicePtr[uiLinkIndex], x[uiLinkIndex]);
        }
        else
        {
            _mul(pDevicePtr[uiLinkIndex], x[uiLinkIndex]);
        }

    gaugeLinkKernelFuncionEnd
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyGaugeReal(deviceGauge* pDevicePtr, const deviceGauge* __restrict__ x, Real a)
{
    gaugeLinkKernelFuncionStart

        _add(pDevicePtr[uiLinkIndex], _mulC(x[uiLinkIndex], a));

    gaugeLinkKernelFuncionEnd
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusGauge(deviceGauge* pDevicePtr, const deviceGauge* __restrict__ x)
{
    gaugeLinkKernelFuncionStart

        _add(pDevicePtr[uiLinkIndex], x[uiLinkIndex]);

    gaugeLinkKernelFuncionEnd
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusGauge(deviceGauge* pDevicePtr, const deviceGauge* __restrict__ x)
{
    gaugeLinkKernelFuncionStart

        _sub(pDevicePtr[uiLinkIndex], x[uiLinkIndex]);

    gaugeLinkKernelFuncionEnd
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyGaugeComplex(deviceGauge* pDevicePtr, CLGComplex a)
{
    gaugeLinkKernelFuncionStart

        _mul(pDevicePtr[uiLinkIndex], a);

    gaugeLinkKernelFuncionEnd
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyGaugeReal(deviceGauge* pDevicePtr, Real a)
{
    gaugeLinkKernelFuncionStart

        _mul(pDevicePtr[uiLinkIndex], a);

    gaugeLinkKernelFuncionEnd
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStapleAtSiteGaugeCacheIndex(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceGauge* pStapleData, //can be NULL
    deviceGauge* pForceData,
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
        deviceGauge res = _makeZero<deviceGauge>();

        //there are 6 staples, each is sum of two plaquttes
        for (INT i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceGauge toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

            if (first.NeedToDagger())
            {
                _dagger(toAdd);
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                const deviceGauge& toMul = pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)];

                if (nextlink.NeedToDagger())
                {
                    _muldag(toAdd, toMul);
                }
                else
                {
                    _mul(toAdd, toMul);
                }
            }
            _add(res, toAdd);
        }
        if (NULL != pStapleData)
        {
            pStapleData[linkIndex] = res;
        }

        //staple calculated
        deviceGauge force = pDeviceData[linkIndex];
        _muldag(force, res);
        _ta(force);
        _mul(force, betaOverN);

        //force is additive
        _add(pForceData[linkIndex], force);
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateOnlyStapleGauge(
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
    UINT plaqLength, UINT plaqCount,
    deviceGauge* pStapleData)
{
    intokernaldir;

    const UINT plaqLengthm1 = plaqLength - 1;
    UINT plaqCountAll = plaqCount * plaqLengthm1;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceGauge res = _makeZero<deviceGauge>();

        //there are 6 staples, each is sum of two plaquttes
        for (INT i = 0; i < plaqCount; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAll];
            deviceGauge toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

            if (first.NeedToDagger())
            {
                _dagger(toAdd);
            }

            for (INT j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAll];
                const deviceGauge& toMul = pDeviceData[_deviceGetLinkIndex(nextlink.m_uiSiteIndex, nextlink.m_byDir)];

                if (nextlink.NeedToDagger())
                {
                    _muldag(toAdd, toMul);
                }
                else
                {
                    _mul(toAdd, toMul);
                }
            }
            _add(res, toAdd);
        }
        pStapleData[linkIndex] = res;
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyGaugeCacheIndex(
    const deviceGauge* __restrict__ pDeviceData,
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
        deviceGauge toAdd = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];

        if (first.NeedToDagger())
        {
            _dagger(toAdd);
        }

        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedIndex[i * plaqLength + j + indexSkip];
            const deviceGauge& toMul = pDeviceData[_deviceGetLinkIndex(first.m_uiSiteIndex, first.m_byDir)];
            if (first.NeedToDagger())
            {
                _muldag(toAdd, toMul);
            }
            else
            {
                _mul(toAdd, toMul);
            }
        }

        resThisThread += (_dim<deviceGauge>() - _retr(toAdd));
    }

    results[uiSiteIndex] = resThisThread * betaOverN;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyGauge_UseClover(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
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
            fRes += _deviceCloverRetrT(pDeviceData, sSite4, __bi(sSite4), byDir1, byDir2, byFieldId);
        }
    }
    fRes = stapleConstant - 0.25 * fRes;
    results[uiSiteIndex] = fRes * fBetaOverN;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteEnergyUsingStableGauge(
    const deviceGauge* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pStableData,
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
        resThisThread += (stapleConstant - _retr(_muldagC(pDeviceData[linkIndex], pStableData[linkIndex])));
    }

    results[uiSiteIndex] = resThisThread * betaOverN * 0.25;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelExpMultGaugeReal(
    const deviceGauge* __restrict__ pMyDeviceData,
    Real a,
    deviceGauge* pU)
{
    gaugeLinkKernelFuncionStart

        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceGauge expP = _expreal(pMyDeviceData[linkIndex], a);
        _mul(expP, pU[linkIndex]);
        pU[linkIndex] = expP;

    gaugeLinkKernelFuncionEnd
}

/**
* Trace (P^2)
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateKinematicEnergyGauge(const deviceGauge* __restrict__ pDeviceData, DOUBLE* results)
{
    intokernaldir;

    DOUBLE resThisThread = 0.0;
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        resThisThread += _retr(_dagmulC(pDeviceData[linkIndex], pDeviceData[linkIndex]));
    }

    results[uiSiteIndex] = resThisThread;
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelNormalizeGauge(deviceGauge* pMyDeviceData)
{
    gaugeLinkKernelFuncionStart

        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        _norm(pMyDeviceData[linkIndex]);

    gaugeLinkKernelFuncionEnd
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDotGauge(
    const deviceGauge* __restrict__ pMyDeviceData,
    const deviceGauge* __restrict__ pOtherDeviceData,
    cuDoubleComplex* result
)
{
    intokernaldir;
    cuDoubleComplex resThisThread = make_cuDoubleComplex(0, 0);
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        resThisThread = cuCadd(resThisThread,
            _cToDouble(_tr(_dagmulC(pMyDeviceData[linkIndex], pOtherDeviceData[linkIndex])))
        );
    }

    result[uiSiteIndex] = resThisThread;
}

/**
 * iA = U.TA()
 */
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToIAGauge(
    deviceGauge* pDeviceData)
{
    gaugeLinkKernelFuncionStart

        _ta(pDeviceData[uiLinkIndex]);

    gaugeLinkKernelFuncionEnd
}

/**
 * E_mu = F_{0 mu}
 */
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToEGauge(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pRes)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    //we only need uiDir = 0,1,2
    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        deviceGauge res = _makeZero<deviceGauge>();
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
        if (dir < uiDir - 1)
        {
            //find clover F
            res = _device1PlaqutteTermPPT(pDeviceData, 3, dir, uiBigIdx, sSite4, byFieldId);
            _ta(res);
            _mul(res, F(-1.0));
        }

        pRes[uiLinkIndex] = res;
    }
}

/**
 * This is wrong! the order of the plaqutte must be considered
 * This is to make sure gauge transform is g(x) nabla E g^+(n)
 */
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateNablaEGauge(
    const deviceGauge* __restrict__ pDeviceData,
    BYTE byFieldId, deviceGauge* pRes)
{
    intokernalInt4;

    //i=0: 12
    //  1: 13
    //  2: 14
    //  3: 23
    //  4: 24
    //  5: 34  

    UINT uiResLinkIdx = _deviceGetLinkIndex(uiSiteIndex, 3);

    pRes[uiResLinkIdx] = _makeZero<deviceGauge>();
    #pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {
        //we need 2, 4 and 5
        //BYTE byPlaqIdx = (dir + 1) << 1;
        //if (byPlaqIdx > 5) byPlaqIdx = 5;

        INT dirs[4];
        //deviceGauge toMul(_devicePlaqutte(pDeviceData, pCachedPlaqutte, uiSiteIndex, byPlaqIdx, plaqLength, plaqCount));
        dirs[0] = 4;
        dirs[1] = dir + 1;
        dirs[2] = -4;
        dirs[3] = -static_cast<INT>(dir) - 1;
        deviceGauge toMul(
            //_device1PlaqutteTermPP(pDeviceData, 3, dir, uiBigIdx)
            _deviceLinkT(pDeviceData, sSite4, 4, byFieldId, dirs)
        );

        //
        dirs[0] = 4;
        dirs[1] = -static_cast<INT>(dir) - 1;
        dirs[2] = -4;
        dirs[3] = dir + 1;

        _mul(toMul,
            _deviceLinkT(pDeviceData, sSite4, 4, byFieldId, dirs)
        );
        //toMul.Ta();
        _sub(pRes[uiResLinkIdx], toMul);
    }
    _ta(pRes[uiResLinkIdx]);
    //pRes[uiResLinkIdx].SubReal(F(3.0));
}

/**
 * Larger than the above
 */
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateNablaENaiveGauge(
    const deviceGauge* __restrict__ pDeviceData,
    BYTE byFieldId, deviceGauge* pRes)
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

    pRes[uiResLinkIdx] = _makeZero<deviceGauge>();
    #pragma unroll
    for (BYTE dir = 0; dir < 3; ++dir)
    {
        //we need 2, 4 and 5
        //BYTE byPlaqIdx = (dir + 1) << 1;
        //if (byPlaqIdx > 5) byPlaqIdx = 5;

        INT dirs[4];
        //deviceGauge toMul(_devicePlaqutte(pDeviceData, pCachedPlaqutte, uiSiteIndex, byPlaqIdx, plaqLength, plaqCount));
        dirs[0] = 4;
        dirs[1] = dir + 1;
        dirs[2] = -4;
        dirs[3] = -static_cast<INT>(dir) - 1;

        deviceGauge a(
            //_device1PlaqutteTermPP(pDeviceData, 3, dir, uiBigIdx)
            _deviceLinkT(pDeviceData, sSite4, 4, byFieldId, dirs)
            //_deviceClover(pDeviceData, sSite4, __bi(sSite4), 3, dir, byFieldId)
        );

        dirs[0] = -static_cast<INT>(dir) - 1;
        dirs[1] = 4;
        dirs[2] = dir + 1;
        dirs[3] = -4;

        deviceGauge b(
            //_device1PlaqutteTermPP(pDeviceData, 3, dir, uiBigIdx)
            _deviceLinkT(pDeviceData, sSite4, 4, byFieldId, dirs)
            //_deviceClover(pDeviceData, sSite4_m_mu, __bi(sSite4_m_mu), 3, dir, byFieldId)
        );
        //b.Add(_deviceLink(pDeviceData, sSite4_m_mu_m_t, 4, byFieldId, dirs));
        //b.Add(_deviceLink(pDeviceData, sSite4_m_2mu, 4, byFieldId, dirs));
        //b.Add(_deviceLink(pDeviceData, sSite4_m_2mu_m_t, 4, byFieldId, dirs));

        //b.Ta();
        _sub(a, b);
        _sub(pRes[uiResLinkIdx], a);
    }
}

/**
 * U = exp(A)
 */
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToUGauge(
    deviceGauge* pDeviceData)
{
    gaugeLinkKernelFuncionStart

        pDeviceData[uiLinkIndex] = _expreal(pDeviceData[uiLinkIndex], F(1.0));

    gaugeLinkKernelFuncionEnd
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelSetOneDirUnityGauge(deviceGauge* pDeviceData, BYTE byDir)
{
    gaugeLinkKernelFuncionStart

        if (0 != ((1 << idir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
            pDeviceData[uiLinkIndex] = _makeId<deviceGauge>();
        }

    gaugeLinkKernelFuncionEnd
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelSetOneDirZeroGauge(deviceGauge* pDeviceData, BYTE byDir)
{
    gaugeLinkKernelFuncionStart

        if (0 != ((1 << idir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
            pDeviceData[uiLinkIndex] = _makeZero<deviceGauge>();
        }

    gaugeLinkKernelFuncionEnd
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelSetOneDirUnityPointGauge(deviceGauge* pDeviceData, UINT uiSiteIndex, BYTE byDir)
{
    for (BYTE dir = 0; dir < 4; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = _makeId<deviceGauge>();
        }
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelSetOneDirZeroPointGauge(deviceGauge* pDeviceData, UINT uiSiteIndex, BYTE byDir)
{
    for (BYTE dir = 0; dir < 4; ++dir)
    {
        if (0 != ((1 << dir) & byDir))
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pDeviceData[uiLinkIndex] = _makeZero<deviceGauge>();
        }
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovLoopOfSiteGauge(
    const deviceGauge* __restrict__ pDeviceBuffer,
    cuDoubleComplex* res)
{
    UINT uiXYZ = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    const UINT uiSiteIndex = uiXYZ * _DC_Lt;
    UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, _DC_Dir - 1);
    SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

    deviceGauge tmp = _makeZero<deviceGauge>();
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
            _mul(tmp, pDeviceBuffer[uiLinkIdx]);
        }
    }

    res[uiXYZ] = _cToDouble(_tr(tmp));
}

#pragma endregion

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::AxpyPlus(const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
    preparethread;
    _kernelAxpyPlusGauge << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::AxpyMinus(const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
    preparethread;
    _kernelAxpyMinusGauge << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData);

}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::ScalarMultply(const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplyGaugeComplex << <block, threads >> > (m_pDeviceData, a);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplyGaugeReal << <block, threads >> > (m_pDeviceData, a);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::Axpy(Real a, const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
    preparethread;
    _kernelAxpyGaugeReal << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData, a);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(x);
    preparethread;
    _kernelAxpyGaugeA << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData, a);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::Mul(const CField* other, UBOOL bDagger)
{
    if (NULL == other || GetFieldType() != other->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: axpy failed because the otherfield is not SU3");
        return;
    }

    const CFieldGaugeLink<deviceGauge, matrixN>* pSU3x = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(other);
    preparethread;
    _kernelMulGauge << <block, threads >> > (m_pDeviceData, pSU3x->m_pDeviceData, bDagger);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::Zero()
{
    preparethread;
    _kernelInitialGaugeField << <block, threads >> > (m_pDeviceData, EFIT_Zero);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::Identity()
{
    preparethread;
    _kernelInitialGaugeField << <block, threads >> > (m_pDeviceData, EFIT_Identity);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::Dagger()
{
    preparethread;
    _kernelDaggerGauge << <block, threads >> > (m_pDeviceData);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::MakeRandomGenerator()
{
    preparethread;
    _kernelInitialGaugeField << <block, threads >> > (m_pDeviceData, EFIT_RandomGenerator);
}

/**
*
*/
template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::InitialField(EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialGaugeField << <block, threads >> > (m_pDeviceData, eInitialType);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eType)
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
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * MatrixN() * MatrixN() * m_uiLinkeCount);
        BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (uiSize != sizeof(Real) * 2 * MatrixN() * MatrixN() * _HC_LinkCount)
        {
            appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(Real) * 2 * MatrixN() * MatrixN() * _HC_LinkCount));
        }
        InitialWithByte(data);
        free(data);
    }
    break;
#if _CLG_DOUBLEFLOAT
    case EFFT_CLGBinFloat:
    {
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * MatrixN() * MatrixN() * m_uiLinkeCount);
        BYTE* data = (BYTE*)malloc(uiSize);
        Real* rdata = (Real*)data;
        FLOAT* fdata = (FLOAT*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (uiSize != sizeof(FLOAT) * 2 * MatrixN() * MatrixN() * _HC_LinkCount)
        {
            appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, sizeof(FLOAT) * 2 * MatrixN() * MatrixN() * _HC_LinkCount);
        }
        for (UINT i = 0; i < 2 * MatrixN() * MatrixN() * m_uiLinkeCount; ++i)
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
        UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * MatrixN() * MatrixN() * m_uiLinkeCount);
        BYTE* data = (BYTE*)malloc(uiSize);
        Real* rdata = (Real*)data;
        DOUBLE* ddata = (DOUBLE*)appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
        if (uiSize != sizeof(DOUBLE) * 2 * MatrixN() * MatrixN() * _HC_LinkCount)
        {
            appCrucial(_T("Loading file size not match: %s, %d, expecting %d"), sFileName.c_str(), uiSize, static_cast<INT>(sizeof(DOUBLE) * 2 * MatrixN() * MatrixN() * _HC_LinkCount));
        }
        for (UINT i = 0; i < 2 * MatrixN() * MatrixN() * m_uiLinkeCount; ++i)
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

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::InitialWithByte(BYTE* byData)
{
    deviceGauge* readData = (deviceGauge*)malloc(sizeof(deviceGauge) * m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[2 * matrixN * matrixN];
        memcpy(oneLink, byData + sizeof(Real) * 2 * MatrixN() * MatrixN() * i, sizeof(Real) * 2 * MatrixN() * MatrixN());
        for (UINT j = 0; j < 2 * MatrixN() * MatrixN(); ++j)
        {
            _setelement(readData[i], j, oneLink[j]);
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceGauge) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    free(readData);
}

/**
* (1) calculate staples
* (2) calculate force(additive)
* (3) calculate energy
*/
template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, Real betaOverN) const
{
    if (NULL == pForce || GetFieldType() != pForce->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: force field is not SU3");
        return;
    }
    if (NULL != pStable && GetFieldType() != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: stape field is not SU3");
        return;
    }

    CFieldGaugeLink<deviceGauge, matrixN>* pForceSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pForce);
    CFieldGaugeLink<deviceGauge, matrixN>* pStableSU3 = NULL == pStable ? NULL : dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pStable);

    preparethread;

    assert(NULL != appGetLattice()->m_pIndexCache->m_pStappleCache);

    _kernelStapleAtSiteGaugeCacheIndex << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        NULL == pStableSU3 ? NULL : pStableSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData,
        betaOverN);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::CalculateOnlyStaple(CFieldGauge* pStable) const
{
    if (NULL == pStable || GetFieldType() != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: stable field is not SU3");
        return;
    }
    CFieldGaugeLink<deviceGauge, matrixN>* pStableSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pStable);

    preparethread;
    _kernelCalculateOnlyStapleGauge << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pStableSU3->m_pDeviceData);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeLink<deviceGauge, matrixN>::CalculatePlaqutteEnergy(DOUBLE betaOverN) const
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);

    preparethread;
    _kernelPlaqutteEnergyGaugeCacheIndex << <block, threads >> > (
        m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer
        );

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeLink<deviceGauge, matrixN>::CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const
{
    assert(NULL != appGetLattice()->m_pIndexCache->m_pPlaqutteCache);
    //appGeneral(_T("const %f\n"), 3.0 * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite);
    preparethread;
    _kernelPlaqutteEnergyGauge_UseClover << <block, threads >> > (
        m_byFieldId,
        m_pDeviceData,
        MatrixN() * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
        betaOverN,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeLink<deviceGauge, matrixN>::CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge* pStable) const
{
    if (NULL == pStable || GetFieldType() != pStable->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: stape field is not SU3");
        return F(0.0);
    }
    const CFieldGaugeLink<deviceGauge, matrixN>* pStableSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(pStable);

    //appGeneral(_T("const = %f\n"), 3.0 * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink);

    preparethread;
    _kernelPlaqutteEnergyUsingStableGauge << <block, threads >> > (
        m_pDeviceData,
        pStableSU3->m_pDeviceData,
        MatrixN() * appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        betaOverN,
        _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
DOUBLE CFieldGaugeLink<deviceGauge, matrixN>::CalculateKinematicEnergy() const
{
    preparethread;
    _kernelCalculateKinematicEnergyGauge << <block, threads >> > (m_pDeviceData, _D_RealThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::SetOneDirectionUnity(BYTE byDir)
{
    if (0 == (byDir & 15))
    {
        return;
    }
    preparethread;
    _kernelSetOneDirUnityGauge << <block, threads >> > (m_pDeviceData, byDir);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::SetOneDirectionZero(BYTE byDir)
{
    if (0 == (byDir & 15))
    {
        return;
    }
    preparethread;
    _kernelSetOneDirZeroGauge << <block, threads >> > (m_pDeviceData, byDir);
}

template<typename deviceGauge, INT matrixN>
CFieldGaugeLink<deviceGauge, matrixN>::CFieldGaugeLink<deviceGauge, matrixN>() : CFieldGauge()
{
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData, sizeof(deviceGauge) * m_uiLinkeCount));
}

template<typename deviceGauge, INT matrixN>
CFieldGaugeLink<deviceGauge, matrixN>::~CFieldGaugeLink<deviceGauge, matrixN>()
{
    checkCudaErrors(__cudaFree(m_pDeviceData));
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::ExpMult(Real a, CField* U) const
{
    if (NULL == U || GetFieldType() != U->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: U field is not SU3");
        return;
    }

    CFieldGaugeLink<deviceGauge, matrixN>* pUField = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(U);

    preparethread;
    _kernelExpMultGaugeReal << < block, threads >> > (m_pDeviceData, a, pUField->m_pDeviceData);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::ElementNormalize()
{
    preparethread;
    _kernelNormalizeGauge << < block, threads >> > (m_pDeviceData);
}

template<typename deviceGauge, INT matrixN>
cuDoubleComplex CFieldGaugeLink<deviceGauge, matrixN>::Dot(const CField* other) const
{
    if (NULL == other || GetFieldType() != other->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: U field is not SUN");
        return make_cuDoubleComplex(0, 0);
    }

    const CFieldGaugeLink<deviceGauge, matrixN>* pUField = dynamic_cast<const CFieldGaugeLink<deviceGauge, matrixN>*>(other);

    preparethread;
    _kernelDotGauge << < block, threads >> > (m_pDeviceData, pUField->m_pDeviceData, _D_ComplexThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::CopyTo(CField* pTarget) const
{
    if (NULL == pTarget || GetFieldType() != pTarget->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: target field is not SUN");
        return;
    }

    CFieldGauge::CopyTo(pTarget);

    CFieldGaugeLink<deviceGauge, matrixN>* pTargetField = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pTarget);
    checkCudaErrors(cudaMemcpy(pTargetField->m_pDeviceData, m_pDeviceData, sizeof(deviceGauge) * m_uiLinkeCount, cudaMemcpyDeviceToDevice));
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::TransformToIA()
{
    preparethread;
    _kernelTransformToIAGauge << <block, threads >> > (m_pDeviceData);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::TransformToU()
{
    preparethread;
    _kernelTransformToUGauge << <block, threads >> > (m_pDeviceData);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::CalculateE_Using_U(CFieldGauge* pResoult) const
{
    if (NULL == pResoult || GetFieldType() != pResoult->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: U field is not SU3");
        return;
    }

    CFieldGaugeLink<deviceGauge, matrixN>* pUField = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pResoult);

    preparethread;
    _kernelTransformToEGauge << <block, threads >> > (pUField->m_byFieldId, m_pDeviceData, pUField->m_pDeviceData);
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::CalculateNablaE_Using_U(CFieldGauge* pResoult, UBOOL bNaive) const
{
    if (NULL == pResoult || GetFieldType() != pResoult->GetFieldType())
    {
        appCrucial("CFieldGaugeLink<deviceGauge, matrixN>: U field is not SU3");
        return;
    }

    CFieldGaugeLink<deviceGauge, matrixN>* pUField = dynamic_cast<CFieldGaugeLink<deviceGauge, matrixN>*>(pResoult);

    preparethread;
    if (bNaive)
    {
        _kernelCalculateNablaENaiveGauge << <block, threads >> > (
            m_pDeviceData,
            m_byFieldId,
            pUField->m_pDeviceData);
    }
    else
    {
        _kernelCalculateNablaEGauge << <block, threads >> > (
            m_pDeviceData,
            m_byFieldId,
            pUField->m_pDeviceData);
    }
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::DebugPrintMe() const
{
    //preparethread;
    //_kernelPrintSU3 << < block, threads >> > (m_pDeviceData);

    //===================================================
    //Since Debug Print Me is only used to debug, we do it slow but convinient
    deviceGauge* pToPrint = (deviceGauge*)malloc(sizeof(deviceGauge) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(pToPrint, m_pDeviceData, sizeof(deviceGauge) * m_uiLinkeCount, cudaMemcpyDeviceToHost));

    for (UINT uiLink = 0; uiLink < m_uiLinkeCount; ++uiLink)
    {
        UINT uiSite = uiLink / _HC_Dir;
        UINT uiDir = uiLink % _HC_Dir;
        SSmallInt4 site = __hostSiteIndexToInt4(uiSite);
        appGeneral(_T(" --- %d(%d, %d, %d, %d)_%d ---\n %s\n"),
            uiLink,
            static_cast<INT>(site.x),
            static_cast<INT>(site.y),
            static_cast<INT>(site.z),
            static_cast<INT>(site.w),
            uiDir,
            appToString(pToPrint[uiLink]).c_str());
    }

    free(pToPrint);
}

template<typename deviceGauge, INT matrixN>
BYTE* CFieldGaugeLink<deviceGauge, matrixN>::CopyDataOut(UINT& uiSize) const
{
    deviceGauge* toSave = (deviceGauge*)malloc(sizeof(deviceGauge) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceGauge) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount * 2 * MatrixN() * MatrixN());
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[2 * matrixN * matrixN];
        for (UINT j = 0; j < 2 * MatrixN() * MatrixN(); ++j)
        {
            oneLink[j] = _element(toSave[i], j);
        }
        memcpy(byToSave + i * sizeof(Real) * 2 * MatrixN() * MatrixN(), oneLink, sizeof(Real) * 2 * MatrixN() * MatrixN());
    }
    free(toSave);

    return byToSave;
}

template<typename deviceGauge, INT matrixN>
BYTE* CFieldGaugeLink<deviceGauge, matrixN>::CopyDataOutFloat(UINT& uiSize) const
{
    deviceGauge* toSave = (deviceGauge*)malloc(sizeof(deviceGauge) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceGauge) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(FLOAT) * m_uiLinkeCount * 2 * MatrixN() * MatrixN());
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[2 * matrixN * matrixN];
        for (UINT j = 0; j < 2 * MatrixN() * MatrixN(); ++j)
        {
            oneLink[j] = static_cast<FLOAT>(_element(toSave[i], j));
        }
        memcpy(byToSave + i * sizeof(FLOAT) * 2 * MatrixN() * MatrixN(), oneLink, sizeof(FLOAT) * 2 * MatrixN() * MatrixN());
    }
    free(toSave);

    return byToSave;
}

template<typename deviceGauge, INT matrixN>
BYTE* CFieldGaugeLink<deviceGauge, matrixN>::CopyDataOutDouble(UINT& uiSize) const
{
    deviceGauge* toSave = (deviceGauge*)malloc(sizeof(deviceGauge) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceGauge) * m_uiLinkeCount, cudaMemcpyDeviceToHost));
    //fuck ofstream
    uiSize = static_cast<UINT>(sizeof(DOUBLE) * m_uiLinkeCount * 2 * MatrixN() * MatrixN());
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        DOUBLE oneLink[2 * matrixN * matrixN];
        for (UINT j = 0; j < 2 * MatrixN() * MatrixN(); ++j)
        {
            oneLink[j] = static_cast<DOUBLE>(_element(toSave[i], j));
        }
        memcpy(byToSave + i * sizeof(DOUBLE) * 2 * MatrixN() * MatrixN(), oneLink, sizeof(DOUBLE) * 2 * MatrixN() * MatrixN());
    }
    free(toSave);

    return byToSave;
}

template<typename deviceGauge, INT matrixN>
void CFieldGaugeLink<deviceGauge, matrixN>::PolyakovOnSpatialSite(cuDoubleComplex* buffer) const
{
    dim3 block(_HC_DecompX, _HC_DecompY, 1);
    dim3 threads(_HC_DecompLx, _HC_DecompLy, 1);
    _kernelPolyakovLoopOfSiteGauge << <block, threads >> > (m_pDeviceData, buffer);
}

#pragma region SU2 functions

__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToULogSU2(deviceSU2* pDeviceData)
{
    gaugeLinkKernelFuncionStart

        pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].StrictExp();

    gaugeLinkKernelFuncionEnd
}

__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToIALogSU2(deviceSU2* pDeviceData)
{
    gaugeLinkKernelFuncionStart

        pDeviceData[uiLinkIndex] = pDeviceData[uiLinkIndex].Log();

    gaugeLinkKernelFuncionEnd
}

#pragma endregion

void CFieldGaugeSU2::InitialWithByteCompressed(const CCString& sFileName)
{
    UINT uiSize = static_cast<UINT>(sizeof(Real) * 3 * m_uiLinkeCount);
    BYTE* byData = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);

    deviceSU2* readData = (deviceSU2*)malloc(sizeof(deviceSU2) * m_uiLinkeCount);
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[3];
        memcpy(oneLink, byData + sizeof(Real) * 3 * i, sizeof(Real) * 3);

        readData[i].m_me[1] = _make_cuComplex(oneLink[0], oneLink[1]);
        readData[i].m_me[2] = _make_cuComplex(-oneLink[0], oneLink[1]);

        readData[i].m_me[0] = _make_cuComplex(F(0.0), oneLink[2]);
        readData[i].m_me[3] = _make_cuComplex(F(0.0), -oneLink[2]);
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceSU2) * m_uiLinkeCount, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaDeviceSynchronize());
    free(readData);

    //DebugPrintMe();

    preparethread;
    _kernelTransformToULogSU2 << <block, threads >> > (m_pDeviceData);
    checkCudaErrors(cudaDeviceSynchronize());

    //DebugPrintMe();
    free(byData);
}

CCString CFieldGaugeSU2::SaveToCompressedFile(const CCString& fileName) const
{
    CFieldGaugeSU2* pPooledGauge = dynamic_cast<CFieldGaugeSU2*>(GetCopy());

    preparethread;
    _kernelTransformToIALogSU2 << <block, threads >> > (pPooledGauge->m_pDeviceData);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    deviceSU2* toSave = (deviceSU2*)malloc(sizeof(deviceSU2) * m_uiLinkeCount);
    checkCudaErrors(cudaMemcpy(toSave, pPooledGauge->m_pDeviceData, sizeof(deviceSU2) * m_uiLinkeCount, cudaMemcpyDeviceToHost));

    //This is a traceless anti-Hermitian now, so we only save part of them
    const UINT uiSize = static_cast<UINT>(sizeof(Real) * m_uiLinkeCount * 3);
    BYTE* byToSave = (BYTE*)malloc(static_cast<size_t>(uiSize));
    for (UINT i = 0; i < m_uiLinkeCount; ++i)
    {
        Real oneLink[3];
        oneLink[0] = static_cast<Real>(toSave[i].m_me[1].x);
        oneLink[1] = static_cast<Real>(toSave[i].m_me[1].y);
        oneLink[2] = static_cast<Real>(toSave[i].m_me[0].x);

        memcpy(byToSave + i * sizeof(Real) * 3, oneLink, sizeof(Real) * 3);
    }

    appGetFileSystem()->WriteAllBytes(fileName.c_str(), byToSave, uiSize);
    //pPooledGauge->DebugPrintMe();
    free(toSave);
    CCString MD5 = CLGMD5Hash(byToSave, uiSize);
    free(byToSave);
    appSafeDelete(pPooledGauge);
    return MD5;
}

__CLGIMPLEMENT_CLASS(CFieldGaugeSU2)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================