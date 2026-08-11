//=============================================================================
// FILENAME : CFieldBosonvn.h
// 
// DESCRIPTION:
// This is the class for the spin fields
//
// REVISION:
//  [mm/dd/yy]
//  [07/04/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldCommonKernel.h"

__BEGIN_NAMESPACE

void CLGAPI appSimpleCopyDD(void* target, const void* source, size_t size)
{
    checkCudaErrors(cudaMemcpy(target, source, size, cudaMemcpyDeviceToDevice));
}

void CLGAPI appSimpleCopyHD(void* target, const void* source, size_t size)
{
    checkCudaErrors(cudaMemcpy(target, source, size, cudaMemcpyHostToDevice));
}

void CLGAPI appSimpleCopyDH(void* target, const void* source, size_t size)
{
    checkCudaErrors(cudaMemcpy(target, source, size, cudaMemcpyDeviceToHost));
}

void CLGAPI appSetGaugeLink(CFieldGauge* pGauge, UINT uiLinkIndex, const CLGComplex* pMatrix)
{
    appAssert(NULL != pGauge && NULL != pMatrix);
    const UINT matrixN = pGauge->MatrixN();
    const size_t uiOffset = static_cast<size_t>(uiLinkIndex) * matrixN * matrixN * sizeof(CLGComplex);
    appSimpleCopyHD(static_cast<BYTE*>(pGauge->GetData()) + uiOffset,
                    pMatrix, matrixN * matrixN * sizeof(CLGComplex));
    _CHECKCUDA;
}

#if _CLG_MULTI_GPU

#pragma region Halo pack / gather

/**
 * Multi-GPU (Phase 1, Improve-1 I8/3.7): gather local boundary sites into a
 * contiguous buffer, one thread per (halo slot, byte-word).
 * pHaloGatherIndex[slot].m_uiSiteIndex is this rank's LOCAL interior source
 * site; uiBytesPerSite bytes are copied from pField[src] to pDst[slot]. Working
 * in 4-byte words keeps every field element (Real / complex / SU3) aligned since
 * every field element is a multiple of 4 bytes.
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelHaloGatherBytes(
    const SIndex* __restrict__ pHaloGatherIndex,
    const BYTE* __restrict__ pFieldBytes,
    BYTE* pDstBytes,
    UINT uiHaloSiteCount,
    UINT uiWordsPerSite)
{
    const UINT uiThread = threadIdx.x + blockDim.x * blockIdx.x;
    const UINT uiTotalWords = uiHaloSiteCount * uiWordsPerSite;
    if (uiThread >= uiTotalWords)
    {
        return;
    }
    const UINT uiSlot = uiThread / uiWordsPerSite;
    const UINT uiWord = uiThread % uiWordsPerSite;
    const SIndex sSource = pHaloGatherIndex[uiSlot];
    //Improve-1 (3.7): debug builds verify the gather source is a legal interior
    //site -- never invalid, never a halo slot (gather reads interior only).
    assert(!sSource.IsInvalid() && sSource.m_uiSiteIndex < _DC_Volume);
    const UINT uiSrcSite = sSource.m_uiSiteIndex;

    const UINT* pSrc = reinterpret_cast<const UINT*>(pFieldBytes) + uiSrcSite * uiWordsPerSite + uiWord;
    UINT* pDst = reinterpret_cast<UINT*>(pDstBytes) + uiSlot * uiWordsPerSite + uiWord;
    *pDst = *pSrc;
}

/**
 * Gather this rank's boundary sites (the ones a neighbour would need) into
 * pDstContiguous, laid out halo-slot-major. bytesPerSite must be a multiple of 4.
 * Returns the number of halo sites gathered.
 */
UINT CLGAPI appHaloGatherLocal(const BYTE* pFieldDevice, BYTE* pDstContiguousDevice, UINT uiBytesPerSite)
{
    const CIndexData* pIdx = appGetLattice()->m_pIndexCache;
    if (NULL == pIdx || 0 == pIdx->m_uiHaloSiteCount || NULL == pIdx->m_pHaloGatherIndex)
    {
        return 0;
    }
    const UINT uiWordsPerSite = uiBytesPerSite / static_cast<UINT>(sizeof(UINT));
    const UINT uiTotalWords = pIdx->m_uiHaloSiteCount * uiWordsPerSite;
    const UINT uiThreadPerBlock = 256;
    const UINT uiBlocks = (uiTotalWords + uiThreadPerBlock - 1) / uiThreadPerBlock;
    //Improve-1 (3.3): infrastructure kernel invoked from CHaloManager refills
    //-- must use the raw backend, never the guarded public macro, or the
    //guard would Ensure -> refill -> guard ... recursively (the manager
    //fail-fasts on that re-entry).
    _CLG_LAUNCH_KERNEL_RAW(_kernelHaloGatherBytes, uiBlocks, uiThreadPerBlock,
        pIdx->m_pHaloGatherIndex, pFieldDevice, pDstContiguousDevice,
        pIdx->m_uiHaloSiteCount, uiWordsPerSite);
    _CHECKCUDA;
    return pIdx->m_uiHaloSiteCount;
}

#pragma endregion

#endif //#if _CLG_MULTI_GPU

#pragma region common

template<typename T>
void CCommonKernel<T>::CopyBuffer(T* dest, const T* source, UINT count)
{
    checkCudaErrors(cudaMemcpy(dest, source, sizeof(T) * count, cudaMemcpyDeviceToDevice));
}

template class CCommonKernel<INT>;
#if _CLG_DOUBLEFLOAT
template class CCommonKernel<FLOAT>;
#else
template class CCommonKernel<DOUBLE>;
#endif
template class CCommonKernel<Real>;
#if _CLG_DOUBLEFLOAT
template class CCommonKernel<cuComplex>;
#else
template class CCommonKernel<cuDoubleComplex>;
#endif
template class CCommonKernel<CLGComplex>;
template class CCommonKernel<deviceSU2>;
template class CCommonKernel<deviceSU3>;
template class CCommonKernel<deviceSU2Vector>;
template class CCommonKernel<deviceSU3Vector>;
template class CCommonKernel<deviceWilsonVectorSU3>;

#if _CLG_SU4_GAUGE
template class CCommonKernel<deviceSU4>;
#endif
#if _CLG_SU5_GAUGE
template class CCommonKernel<deviceSU5>;
#endif
#if _CLG_SU6_GAUGE
template class CCommonKernel<deviceSU6>;
#endif
#if _CLG_SU7_GAUGE
template class CCommonKernel<deviceSU7>;
#endif
#if _CLG_SU8_GAUGE
template class CCommonKernel<deviceSU8>;
#endif
#if _CLG_Z2_GAUGE
template class CCommonKernel<deviceZN<2>>;
#endif
#if _CLG_Z3_GAUGE
template class CCommonKernel<deviceZN<3>>;
#endif
#if _CLG_Z4_GAUGE
template class CCommonKernel<deviceZN<4>>;
#endif
#if _CLG_Z5_GAUGE
template class CCommonKernel<deviceZN<5>>;
#endif
#if _CLG_Z6_GAUGE
template class CCommonKernel<deviceZN<6>>;
#endif
#if _CLG_D3_GAUGE
template class CCommonKernel<deviceDN<3>>;
#endif
#if _CLG_D4_GAUGE
template class CCommonKernel<deviceDN<4>>;
#endif
#if _CLG_D8_GAUGE
template class CCommonKernel<deviceDN<8>>;
#endif

#if _CLG_SU4_KS || _CLG_SU4_BOSON
template class CCommonKernel<deviceSU4Vector>;
#endif
#if _CLG_SU5_KS || _CLG_SU5_BOSON
template class CCommonKernel<deviceSU5Vector>;
#endif
#if _CLG_SU6_KS || _CLG_SU6_BOSON
template class CCommonKernel<deviceSU6Vector>;
#endif
#if _CLG_SU7_KS || _CLG_SU7_BOSON
template class CCommonKernel<deviceSU7Vector>;
#endif
#if _CLG_SU8_KS || _CLG_SU8_BOSON
template class CCommonKernel<deviceSU8Vector>;
#endif

#if _CLG_SL3C_GAUGE
template class CCommonKernel<deviceSL3C>;
#endif
#if _CLG_U3_GAUGE
template class CCommonKernel<deviceU3>;
#endif
#if _CLG_O3_GAUGE
template class CCommonKernel<deviceO3>;
#endif
#if _CLG_SO3_GAUGE
template class CCommonKernel<deviceSO3>;
#endif

#pragma endregion

#pragma region field

#pragma region kernel

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialCommon(T* pDevicePtr, UINT count, EFieldInitialType eInitialType)
{
    __simplekernel(
        switch (eInitialType)
        {
        case EFIT_Zero:
            {
                _Zero(pDevicePtr[idx]);
            }
            break;
        case EFIT_Identity:
            {
                _Id(pDevicePtr[idx]);
            }
            break;
        case EFIT_RandomGaussian:
        case EFIT_RandomGenerator:
            {
                pDevicePtr[idx] = _makeGaussian<T>(_deviceGetLinkIndex(idx, 0));
            }
            break;
        case EFIT_Random:
            {
                pDevicePtr[idx] = _makeRandom<T>(_deviceGetLinkIndex(idx, 0));
            }
            break;
        default:
            {
                printf("Field cannot be initialized with this type! %d\n", eInitialType);
            }
            break;
        }
    )
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialCommonEO(T* pDevicePtr, UINT count, EFieldInitialType eInitialType, 
    UBOOL bSite, UBOOL bEven, UBOOL bZeroOtherSites,
    const BYTE* __restrict__ etatable)
{
    __simplekernel(
        const BYTE eta5 = bSite ? (etatable[idx] >> 4) : (etatable[idx / _DC_Dir] >> 4);
        if (eta5)
        {
            if (bEven)
            {
                if (bZeroOtherSites)
                {
                    _Zero(pDevicePtr[idx]);
                }
                return;
            }
        }
        else
        {
            if (!bEven)
            {
                if (bZeroOtherSites)
                {
                    _Zero(pDevicePtr[idx]);
                }
                return;
            }
        }
        switch (eInitialType)
        {
        case EFIT_Zero:
            {
                _Zero(pDevicePtr[idx]);
            }
            break;
        case EFIT_Identity:
            {
                _Id(pDevicePtr[idx]);
            }
            break;
        case EFIT_RandomGaussian:
        case EFIT_RandomGenerator:
            {
                pDevicePtr[idx] = _makeGaussian<T>(_deviceGetLinkIndex(idx, 0));
            }
            break;
        case EFIT_Random:
            {
                pDevicePtr[idx] = _makeRandom<T>(_deviceGetLinkIndex(idx, 0));
            }
            break;
        default:
            {
                printf("Field cannot be initialized with this type! %d\n", eInitialType);
            }
            break;
        }
    )
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelZeroEO(T* pDevicePtr, UINT count, UBOOL bEven,
    const BYTE* __restrict__ pEtaTable)
{
    //__simplekernel(
    //    if (etatable[idx] >> 4)
    //    {
    //        if (bEven)
    //        {
    //            return;
    //        }
    //    }
    //    else
    //    {
    //        if (!bEven)
    //        {
    //            return;
    //        }
    //    }
    //    _Zero(pDevicePtr[idx]);
    //    )

    intokernalEOHalf;
    _Zero(pDevicePtr[uiSiteIndex]);
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd(T* pMe, const T* __restrict__ pOther, UINT count)
{
    __simplekernel(_add(pMe[idx], pOther[idx]))
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelSub(T* pMe, const T* __restrict__ pOther, UINT count)
{
    __simplekernel(_sub(pMe[idx], pOther[idx]))
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelMul(T* pMe, const T* __restrict__ pOther, UINT count, UBOOL bDaggerLeft, UBOOL bDaggerRight)
{
    __simplekernel(
        if (bDaggerLeft)
        {
            _dagger(pMe[idx]);
        }
        if (bDaggerRight)
        {
            _muldag(pMe[idx], pOther[idx]);
        }
        else
        {
            _mul(pMe[idx], pOther[idx]);
        }
    )
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelApplyPhaseC(T* pMe, const CLGComplex* __restrict__ pOther, UINT count)
{
    __simplekernel(
        _mul(pMe[idx], pOther[idx]);
        )
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelApplyPhaseR(T* pMe, const Real* __restrict__ pOther, Real fCharge, UINT count)
{
    __simplekernel(
        const Real fPhase = pOther[idx] * fCharge;
        _mul(pMe[idx], _make_cuComplex(_cos(fPhase), _sin(fPhase)));
        )
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelLeftMul(T* pMe, const T* __restrict__ pOther, UINT count, UBOOL bDaggerLeft, UBOOL bDaggerRight)
{
    __simplekernel(
        if (bDaggerLeft)
        {
            if (bDaggerRight)
            {
                pMe[idx] = _muldagC(_daggerC(pOther[idx]), pMe[idx]);
            }
            else
            {
                pMe[idx] = _mulC(_daggerC(pOther[idx]), pMe[idx]);
            }
            
        }
        else
        {
            if (bDaggerRight)
            {
                pMe[idx] = _muldagC(pOther[idx], pMe[idx]);
            }
            else
            {
                pMe[idx] = _mulC(pOther[idx], pMe[idx]);
            }
        }
    )
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpy(T* pMe, const T* __restrict__ pOther, CLGComplex a, UINT count)
{
    __simplekernel(_add(pMe[idx], _mulC(pOther[idx], a)))
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyReal(T* pMe, const T* __restrict__ pOther, Real a, UINT count)
{
    __simplekernel(_add(pMe[idx], _mulC(pOther[idx], a)))
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddEO(T* pMe, const T* __restrict__ pOther, UBOOL bEven,
    const BYTE* __restrict__ pEtaTable)
{
    intokernalEOHalf;
    _add(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelSubEO(T* pMe, const T* __restrict__ pOther, UBOOL bEven,
    const BYTE* __restrict__ pEtaTable)
{
    intokernalEOHalf;
    _sub(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyEO(T* pMe, const T* __restrict__ pOther, CLGComplex a, UBOOL bEven,
    const BYTE* __restrict__ pEtaTable)
{
    intokernalEOHalf;
    _add(pMe[uiSiteIndex], _mulC(pOther[uiSiteIndex], a));
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyRealEO(T* pMe, const T* __restrict__ pOther, Real a, UBOOL bEven,
    const BYTE* __restrict__ pEtaTable)
{
    intokernalEOHalf;
    _add(pMe[uiSiteIndex], _mulC(pOther[uiSiteIndex], a));
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelDot(const T* __restrict__ pMe, const T* __restrict__ pOther, cuDoubleComplex* result, UINT count)
{
    __simplekernel(result[idx] = _cToDouble(_dot(pMe[idx], pOther[idx])))
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelDotDir(const T* __restrict__ pMe, const T* __restrict__ pOther, cuDoubleComplex* result, UINT count)
{
    __simplekernel(
        const UINT siteidx = idx / _DC_Dir;
        const UINT dir = idx % _DC_Dir;
        if (0 == dir)
        {
            result[siteidx] = _cToDouble(_dot(pMe[idx], pOther[idx]));
        }
        __syncthreads();
        if (1 == dir)
        {
            result[siteidx] = cuCadd(result[siteidx], _cToDouble(_dot(pMe[idx], pOther[idx])));
        }
        __syncthreads();
        if (2 == dir)
        {
            result[siteidx] = cuCadd(result[siteidx], _cToDouble(_dot(pMe[idx], pOther[idx])));
        }
        __syncthreads();
        if (3 == dir)
        {
            result[siteidx] = cuCadd(result[siteidx], _cToDouble(_dot(pMe[idx], pOther[idx])));
        }
    )
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelLengthSq(const T* __restrict__ pMe, DOUBLE* result, UINT count)
{
    __simplekernel(result[idx] = static_cast<DOUBLE>(_lensq(pMe[idx])))
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelLengthSqDir(const T* __restrict__ pMe, DOUBLE* result, UINT count)
{
    __simplekernel(
        const UINT siteidx = idx / _DC_Dir;
        const UINT dir = idx % _DC_Dir;

        //see quda, minus constant 4 is for numerical stability
        const DOUBLE generatorNum = static_cast<DOUBLE>(_generator_count<T>());
        if (0 == dir)
        {
            result[siteidx] = static_cast<DOUBLE>(_lensq(pMe[idx])) - generatorNum;
        }
        __syncthreads();
        if (1 == dir)
        {
            result[siteidx] += static_cast<DOUBLE>(_lensq(pMe[idx])) - generatorNum;
        }
        __syncthreads();
        if (2 == dir)
        {
            result[siteidx] += static_cast<DOUBLE>(_lensq(pMe[idx])) - generatorNum;
        }
        __syncthreads();
        if (3 == dir)
        {
            result[siteidx] += static_cast<DOUBLE>(_lensq(pMe[idx])) - generatorNum;
        }
    )
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelElement(const T* __restrict__ pMe, UINT elementidx, DOUBLE* result, UINT count)
{
    __simplekernel(result[idx] = static_cast<DOUBLE>(_element(pMe[idx], elementidx)))
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelElementDir(const T* __restrict__ pMe, UINT elementidx, DOUBLE* result, UINT count)
{
    __simplekernel(
        const UINT siteidx = idx / _DC_Dir;
        const UINT dir = idx % _DC_Dir;
        if (0 == dir)
        {
            result[siteidx] = static_cast<DOUBLE>(_element(pMe[idx], elementidx));
        }
        __syncthreads();
        if (1 == dir)
        {
            result[siteidx] += static_cast<DOUBLE>(_element(pMe[idx], elementidx));
        }
        __syncthreads();
        if (2 == dir)
        {
            result[siteidx] += static_cast<DOUBLE>(_element(pMe[idx], elementidx));
        }
        __syncthreads();
        if (3 == dir)
        {
            result[siteidx] += static_cast<DOUBLE>(_element(pMe[idx], elementidx));
        }
    )
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiply(T* pMe, CLGComplex a, UINT count)
{
    __simplekernel(_mul(pMe[idx], a))
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyReal(T* pMe, Real a, UINT count)
{
    __simplekernel(_mul(pMe[idx], a))
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelConjugate(T* pDeviceData, UINT count)
{
    __simplekernel(_dagger(pDeviceData[idx]))
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelNorm(T* pDeviceData, UINT count)
{
    __simplekernel(_norm(pDeviceData[idx]))
}

#pragma endregion

template<typename T>
void CCommonKernelField<T>::Initial(T* pointer, UINT count, EFieldInitialType eInitialType)
{
    if (count > _HC_Volume && (eInitialType == EFIT_Random || eInitialType == EFIT_RandomGaussian || eInitialType == EFIT_RandomZ4))
    {
        appCrucial(_T("calling initial with random but count is more than site number which is not supported!\n"))
    }
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelInitialCommon<T>, block, thread, pointer, count, eInitialType);
    //checkCudaErrors(cudaDeviceSynchronize());
    //checkCudaErrors(cudaGetLastError());
}

template<typename T>
void CCommonKernelField<T>::InitialEvenOdd(T* pointer, UINT count, EFieldInitialType eInitialType, UBOOL bSite, UBOOL bEven, UBOOL bZeroOtherSites)
{
    if (count > _HC_Volume && (eInitialType == EFIT_Random || eInitialType == EFIT_RandomGaussian || eInitialType == EFIT_RandomZ4))
    {
        appCrucial(_T("calling initial with random but count is more than site number which is not supported!\n"))
    }
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelInitialCommonEO<T>, block, thread, pointer, count, eInitialType,
        bSite, bEven, bZeroOtherSites,
        appGetLattice()->m_pIndexCache->m_pEtaMu);
    //checkCudaErrors(cudaDeviceSynchronize());
    //checkCudaErrors(cudaGetLastError());
}

template<typename T>
void CCommonKernelField<T>::ZeroEvenOddSite(T* pointer, UINT count, UBOOL bEven)
{
    //__SIMPLEDECOMPOSE(count);
    //_LAUNCH_KERNEL(_kernelZeroEO, block, thread, pointer, count, bEven,
    //    appGetLattice()->m_pIndexCache->m_pEtaMu);

    //NOTE! when use intokernalEOHalf, bEven is to skip even, so here we use !bEven
    preparethreadHalf;
    const UBOOL bOdd = !bEven;
    _LAUNCH_KERNEL(_kernelZeroEO<T>, block, threads, pointer, count, bOdd,
        appGetLattice()->m_pIndexCache->m_pEtaMu);
}

template<typename T>
BYTE* CCommonKernelField<T>::CopyDataOut(T* pointer, UINT count, UINT& uiSize)
{
    T* toSave = (T*)malloc(sizeof(T) * count);
    uiSize = static_cast<UINT>(sizeof(Real) * count * _elementdim<T>());
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    Real* fsaveData = (Real*)saveData;
    checkCudaErrors(cudaMemcpy(toSave, pointer, sizeof(T) * count, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < count; ++i)
    {
        for (UINT j = 0; j < _elementdim<T>(); ++j)
        {
            fsaveData[_elementdim<T>() * i + j] = _element(toSave[i], j);
        }
    }
    free(toSave);
    return saveData;
}

template<typename T>
BYTE* CCommonKernelField<T>::CopyDataOutFloat(T* pointer, UINT count, UINT& uiSize)
{
    T* toSave = (T*)malloc(sizeof(T) * count);
    uiSize = static_cast<UINT>(sizeof(FLOAT) * count * _elementdim<T>());
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    FLOAT* fsaveData = (FLOAT*)saveData;
    checkCudaErrors(cudaMemcpy(toSave, pointer, sizeof(T) * count, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < count; ++i)
    {
        for (UINT j = 0; j < _elementdim<T>(); ++j)
        {
            fsaveData[_elementdim<T>() * i + j] = static_cast<FLOAT>(_element(toSave[i], j));
        }
    }
    free(toSave);
    return saveData;
}

template<typename T>
BYTE* CCommonKernelField<T>::CopyDataOutDouble(T* pointer, UINT count, UINT& uiSize)
{
    T* toSave = (T*)malloc(sizeof(T) * count);
    uiSize = static_cast<UINT>(sizeof(DOUBLE) * count * _elementdim<T>());
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    DOUBLE* fsaveData = (DOUBLE*)saveData;
    checkCudaErrors(cudaMemcpy(toSave, pointer, sizeof(T) * count, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < count; ++i)
    {
        for (UINT j = 0; j < _elementdim<T>(); ++j)
        {
            fsaveData[_elementdim<T>() * i + j] = static_cast<DOUBLE>(_element(toSave[i], j));
        }
    }
    free(toSave);
    return saveData;
}

template<typename T>
void CCommonKernelField<T>::InitialWithByte(T* pointer, UINT count, const BYTE* byData)
{
    T* readData = (T*)malloc(sizeof(T) * count);
    for (UINT i = 0; i < count; ++i)
    {
        appAssert(_elementdim<T>() <= CCString::_CLG_MAX_PATH);
        Real thisSite[CCString::_CLG_MAX_PATH];
        memcpy(thisSite, byData + i * sizeof(Real) * _elementdim<T>(), sizeof(Real) * _elementdim<T>());
        for (UINT k = 0; k < _elementdim<T>(); ++k)
        {
            _setelement(readData[i], k, thisSite[k]);
        }
    }
    checkCudaErrors(cudaMemcpy(pointer, readData, sizeof(T) * count, cudaMemcpyHostToDevice));
    free(readData);
}

template<typename T>
void CCommonKernelField<T>::Dagger(T* dest, UINT count)
{
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelConjugate<T>, block, thread, dest, count);
}

template<typename T>
void CCommonKernelField<T>::AxpyPlus(T* dest, UINT count, const T* x)
{
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelAdd<T>, block, thread, dest, x, count);
}

template<typename T>
void CCommonKernelField<T>::AxpyMinus(T* dest, UINT count, const T* x)
{
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelSub<T>, block, thread, dest, x, count);
}

template<typename T>
void CCommonKernelField<T>::Axpy(T* dest, UINT count, Real a, const T* x)
{
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelAxpyReal<T>, block, thread, dest, x, a, count);
}

template<typename T>
void CCommonKernelField<T>::Axpy(T* dest, UINT count, const CLGComplex& a, const T* x)
{
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelAxpy<T>, block, thread, dest, x, a, count);
}

template<typename T>
void CCommonKernelField<T>::AxpyPlusEvenOdd(T* dest, UINT count, const T* x, UBOOL bEven)
{
    //NOTE! when use intokernalEOHalf, bEven is to skip even, so here we use !bEven
    preparethreadHalf;
    const UBOOL bOdd = !bEven;
    _LAUNCH_KERNEL(_kernelAddEO<T>, block, threads, dest, x, bOdd,
        appGetLattice()->m_pIndexCache->m_pEtaMu);
}

template<typename T>
void CCommonKernelField<T>::AxpyMinusEvenOdd(T* dest, UINT count, const T* x, UBOOL bEven)
{
    //NOTE! when use intokernalEOHalf, bEven is to skip even, so here we use !bEven
    preparethreadHalf;
    const UBOOL bOdd = !bEven;
    _LAUNCH_KERNEL(_kernelSubEO<T>, block, threads, dest, x, bOdd,
        appGetLattice()->m_pIndexCache->m_pEtaMu);
}

template<typename T>
void CCommonKernelField<T>::AxpyEvenOdd(T* dest, UINT count, Real a, const T* x, UBOOL bEven)
{
    //NOTE! when use intokernalEOHalf, bEven is to skip even, so here we use !bEven
    preparethreadHalf;
    const UBOOL bOdd = !bEven;
    _LAUNCH_KERNEL(_kernelAxpyRealEO<T>, block, threads, dest, x, a, bOdd,
        appGetLattice()->m_pIndexCache->m_pEtaMu);
}

template<typename T>
void CCommonKernelField<T>::AxpyEvenOdd(T* dest, UINT count, const CLGComplex& a, const T* x, UBOOL bEven)
{
    //NOTE! when use intokernalEOHalf, bEven is to skip even, so here we use !bEven
    preparethreadHalf;
    const UBOOL bOdd = !bEven;
    _LAUNCH_KERNEL(_kernelAxpyEO<T>, block, threads, dest, x, a, bOdd,
        appGetLattice()->m_pIndexCache->m_pEtaMu);
}

template<typename T>
void CCommonKernelField<T>::Mul(T* dest, UINT count, const T* x, UBOOL bDaggerLeft, UBOOL bDaggerRight)
{
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelMul<T>, block, thread, dest, x, count, bDaggerLeft, bDaggerRight);
}

template<typename T>
void CCommonKernelField<T>::ApplyPhaseC(T* dest, UINT count, const CLGComplex* other)
{
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelApplyPhaseC<T>, block, thread, dest, other, count);
}

template<typename T>
void CCommonKernelField<T>::ApplyPhaseR(T* dest, UINT count, const Real* other, Real fCharge)
{
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelApplyPhaseR<T>, block, thread, dest, other, fCharge, count);
}

template<typename T>
void CCommonKernelField<T>::LeftMul(T* dest, UINT count, const T* x, UBOOL bDaggerLeft, UBOOL bDaggerRight)
{
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelLeftMul<T>, block, thread, dest, x, count, bDaggerLeft, bDaggerRight);
}

template<typename T>
void CCommonKernelField<T>::ScalarMultply(T* dest, UINT count, const CLGComplex& a)
{
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelScalarMultiply<T>, block, thread, dest, a, count);
}

template<typename T>
void CCommonKernelField<T>::ScalarMultply(T* dest, UINT count, Real a)
{
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelScalarMultiplyReal<T>, block, thread, dest, a, count);
}

template<typename T>
cuDoubleComplex CCommonKernelField<T>::Dot(const T* me, UINT count, const T* other)
{
    if (count == _HC_Volume)
    {
        __SIMPLEDECOMPOSE(count);
        _LAUNCH_KERNEL(_kernelDot<T>, block, thread, me, other, _D_ComplexThreadBuffer, count);
        _RECORD2(CCommonKernelField::Dot::ThreadBufferSum, a);
        return _clgGlobalThreadBufferSum(_D_ComplexThreadBuffer);
    }
    else if (count == _HC_Volume * _HC_Dir)
    {
        //__SIMPLEDECOMPOSE(count);
        preparethreadDir;
        _LAUNCH_KERNEL(_kernelDotDir<T>, block, threads, me, other, _D_ComplexThreadBuffer, count);
        return _clgGlobalThreadBufferSum(_D_ComplexThreadBuffer);
    }
    appCrucial(_T("not supported!\n"));
    return make_cuDoubleComplex(0.0, 0.0);
}

template<typename T>
DOUBLE CCommonKernelField<T>::LengthSq(const T* me, UINT count)
{
    if (count == _HC_Volume)
    {
        __SIMPLEDECOMPOSE(count);
        _LAUNCH_KERNEL(_kernelLengthSq<T>, block, thread, me, _D_RealThreadBuffer, count);
        return _clgGlobalThreadBufferSum(_D_RealThreadBuffer);
    }
    else if (count == _HC_Volume * _HC_Dir)
    {
        //__SIMPLEDECOMPOSE(count);
        preparethreadDir;
        _LAUNCH_KERNEL(_kernelLengthSqDir<T>, block, threads, me, _D_RealThreadBuffer, count);
        return _clgGlobalThreadBufferSum(_D_RealThreadBuffer);
    }
    appCrucial(_T("not supported!\n"));
    return 0.0;
}

template<typename T>
TArray<DOUBLE> CCommonKernelField<T>::Sum(const T* me, UINT count)
{
    if (count == _HC_Volume)
    {
        __SIMPLEDECOMPOSE(count);
        TArray<DOUBLE> ret;
        for (UINT i = 0; i < _elementdim<T>(); ++i)
        {
            _LAUNCH_KERNEL(_kernelElement<T>, block, thread, me, i, _D_RealThreadBuffer, count);
            ret.AddItem(_clgGlobalThreadBufferSum(_D_RealThreadBuffer));
        }
        return ret;
    }
    else if (count == _HC_Volume * _HC_Dir)
    {
        //__SIMPLEDECOMPOSE(count);
        preparethreadDir;
        TArray<DOUBLE> ret;
        for (UINT i = 0; i < _elementdim<T>(); ++i)
        {
            _LAUNCH_KERNEL(_kernelElementDir<T>, block, threads, me, i, _D_RealThreadBuffer, count);
            ret.AddItem(_clgGlobalThreadBufferSum(_D_RealThreadBuffer));
        }
        return ret;
    }
    appCrucial(_T("not supported!\n"));
    return TArray<DOUBLE>();
}

template<typename T>
void CCommonKernelField<T>::Norm(T* dest, UINT count)
{
    __SIMPLEDECOMPOSE(count);
    _LAUNCH_KERNEL(_kernelNorm<T>, block, thread, dest, count);
}

template class CCommonKernelField<Real>;
template class CCommonKernelField<CLGComplex>;
template class CCommonKernelField<deviceSU2>;
template class CCommonKernelField<deviceSU3>;
template class CCommonKernelField<deviceSU2Vector>;
template class CCommonKernelField<deviceSU3Vector>;
template class CCommonKernelField<deviceWilsonVectorSU3>;

#if _CLG_SU4_GAUGE
template class CCommonKernelField<deviceSU4>;
#endif
#if _CLG_SU5_GAUGE
template class CCommonKernelField<deviceSU5>;
#endif
#if _CLG_SU6_GAUGE
template class CCommonKernelField<deviceSU6>;
#endif
#if _CLG_SU7_GAUGE
template class CCommonKernelField<deviceSU7>;
#endif
#if _CLG_SU8_GAUGE
template class CCommonKernelField<deviceSU8>;
#endif

#if _CLG_Z2_GAUGE
template class CCommonKernelField<deviceZN<2>>;
#endif
#if _CLG_Z3_GAUGE
template class CCommonKernelField<deviceZN<3>>;
#endif
#if _CLG_Z4_GAUGE
template class CCommonKernelField<deviceZN<4>>;
#endif
#if _CLG_Z5_GAUGE
template class CCommonKernelField<deviceZN<5>>;
#endif
#if _CLG_Z6_GAUGE
template class CCommonKernelField<deviceZN<6>>;
#endif
#if _CLG_D3_GAUGE
template class CCommonKernelField<deviceDN<3>>;
#endif
#if _CLG_D4_GAUGE
template class CCommonKernelField<deviceDN<4>>;
#endif
#if _CLG_D8_GAUGE
template class CCommonKernelField<deviceDN<8>>;
#endif

#if _CLG_SU4_KS || _CLG_SU4_BOSON
template class CCommonKernelField<deviceSU4Vector>;
#endif
#if _CLG_SU5_KS || _CLG_SU5_BOSON
template class CCommonKernelField<deviceSU5Vector>;
#endif
#if _CLG_SU6_KS || _CLG_SU6_BOSON
template class CCommonKernelField<deviceSU6Vector>;
#endif
#if _CLG_SU7_KS || _CLG_SU7_BOSON
template class CCommonKernelField<deviceSU7Vector>;
#endif
#if _CLG_SU8_KS || _CLG_SU8_BOSON
template class CCommonKernelField<deviceSU8Vector>;
#endif

#if _CLG_SL3C_GAUGE
template class CCommonKernelField<deviceSL3C>;
#endif
#if _CLG_U3_GAUGE
template class CCommonKernelField<deviceU3>;
#endif
#if _CLG_O3_GAUGE
template class CCommonKernelField<deviceO3>;
#endif
#if _CLG_SO3_GAUGE
template class CCommonKernelField<deviceSO3>;
#endif


#pragma endregion

#pragma region site

#pragma region kernels

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialSite(T* pDevicePtr, BYTE byFieldId, EFieldInitialType eInitialType)
{
    intokernalInt4;
    const UINT uiBigIndex = __bi(sSite4);

    switch (eInitialType)
    {
    case EFIT_Zero:
    {
        pDevicePtr[uiSiteIndex] = _makeZero<T>();
    }
    break;
    case EFIT_Identity:
    {
        pDevicePtr[uiSiteIndex] = _makeId<T>();
    }
    break;
    case EFIT_RandomGaussian:
    case EFIT_RandomGenerator:
    {
        if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
        {
            if (NULL != __boundaryFieldPointers[byFieldId])
            {
                UINT uiRegion = __idx->_devcieExchangeBoundaryFieldSiteIndexBI(byFieldId, uiBigIndex);
                const T* buffer = ((CFieldBoundary<T>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData;
                pDevicePtr[uiSiteIndex] = buffer[uiRegion];
            }
            pDevicePtr[uiSiteIndex] = _makeZero<T>();
            return;
        }
        pDevicePtr[uiSiteIndex] = _makeGaussian<T>(_deviceGetLinkIndex(uiSiteIndex, 0));
    }
    break;
    case EFIT_Random:
    {
        if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
        {
            if (NULL != __boundaryFieldPointers[byFieldId])
            {
                UINT uiRegion = __idx->_devcieExchangeBoundaryFieldSiteIndexBI(byFieldId, uiBigIndex);
                const T* buffer = ((CFieldBoundary<T>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData;
                pDevicePtr[uiSiteIndex] = buffer[uiRegion];
            }
            pDevicePtr[uiSiteIndex] = _makeZero<T>();
            return;
        }
        pDevicePtr[uiSiteIndex] = _makeRandom<T>(_deviceGetLinkIndex(uiSiteIndex, 0));
    }
    break;
    case EFIT_RandomZ4:
    {
        if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
        {
            if (NULL != __boundaryFieldPointers[byFieldId])
            {
                UINT uiRegion = __idx->_devcieExchangeBoundaryFieldSiteIndexBI(byFieldId, uiBigIndex);
                const T* buffer = ((CFieldBoundary<T>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData;
                pDevicePtr[uiSiteIndex] = buffer[uiRegion];
            }
            pDevicePtr[uiSiteIndex] = _makeZero<T>();
            return;
        }
        pDevicePtr[uiSiteIndex] = _makeZ4<T>(_deviceGetLinkIndex(uiSiteIndex, 0));
    }
    break;
    default:
    {
        printf("_kernelInitialSite cannot be initialized with this type! %d\n", eInitialType);
    }
    break;
    }
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialSiteEvenOdd(T* pDevicePtr, BYTE byFieldId, EFieldInitialType eInitialType, UBOOL bEven, UBOOL bZeroOtherSites)
{
    intokernalInt4;
    if (sSite4.IsOdd())
    {
        if (bEven)
        {
            if (bZeroOtherSites)
            {
                pDevicePtr[uiSiteIndex] = _makeZero<T>();
            }
            return;
        }
    }
    else
    {
        if (!bEven)
        {
            if (bZeroOtherSites)
            {
                pDevicePtr[uiSiteIndex] = _makeZero<T>();
            }
            return;
        }
    }

    const UINT uiBigIndex = __bi(sSite4);

    switch (eInitialType)
    {
    case EFIT_Zero:
    {
        pDevicePtr[uiSiteIndex] = _makeZero<T>();
    }
    break;
    case EFIT_Identity:
    {
        pDevicePtr[uiSiteIndex] = _makeId<T>();
    }
    break;
    case EFIT_RandomGaussian:
    case EFIT_RandomGenerator:
    {
        if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
        {
            if (NULL != __boundaryFieldPointers[byFieldId])
            {
                UINT uiRegion = __idx->_devcieExchangeBoundaryFieldSiteIndexBI(byFieldId, uiBigIndex);
                const T* buffer = ((CFieldBoundary<T>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData;
                pDevicePtr[uiSiteIndex] = buffer[uiRegion];
            }
            pDevicePtr[uiSiteIndex] = _makeZero<T>();
            return;
        }
        pDevicePtr[uiSiteIndex] = _makeGaussian<T>(_deviceGetLinkIndex(uiSiteIndex, 0));
    }
    break;
    case EFIT_Random:
    {
        if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
        {
            if (NULL != __boundaryFieldPointers[byFieldId])
            {
                UINT uiRegion = __idx->_devcieExchangeBoundaryFieldSiteIndexBI(byFieldId, uiBigIndex);
                const T* buffer = ((CFieldBoundary<T>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData;
                pDevicePtr[uiSiteIndex] = buffer[uiRegion];
            }
            pDevicePtr[uiSiteIndex] = _makeZero<T>();
            return;
        }
        pDevicePtr[uiSiteIndex] = _makeRandom<T>(_deviceGetLinkIndex(uiSiteIndex, 0));
    }
    break;
    case EFIT_RandomZ4:
    {
        if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
        {
            if (NULL != __boundaryFieldPointers[byFieldId])
            {
                UINT uiRegion = __idx->_devcieExchangeBoundaryFieldSiteIndexBI(byFieldId, uiBigIndex);
                const T* buffer = ((CFieldBoundary<T>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData;
                pDevicePtr[uiSiteIndex] = buffer[uiRegion];
            }
            pDevicePtr[uiSiteIndex] = _makeZero<T>();
            return;
        }
        pDevicePtr[uiSiteIndex] = _makeZ4<T>(_deviceGetLinkIndex(uiSiteIndex, 0));
    }
    break;
    default:
    {
        printf("_kernelInitialSite cannot be initialized with this type! %d\n", eInitialType);
    }
    break;
    }
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelFixBoundarySite(T* pDeviceData, BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIndex = __idx->_deviceGetBigIndex(sSite4);

    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
    {
        if (NULL != __boundaryFieldPointers[byFieldId])
        {
            const UINT uiRegion = __idx->_devcieExchangeBoundaryFieldSiteIndexBI(byFieldId, uiBigIndex);
            pDeviceData[uiSiteIndex] = ((CFieldBoundary<T>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData[uiRegion];
        }
        else
        {
            pDeviceData[uiSiteIndex] = _makeZero<T>();
        }
        return;
    }
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelMakePointSourceSite(T* pDeviceData, UINT uiDesiredSite, BYTE spin, BYTE byColor)
{
    intokernal;
    if (uiSiteIndex == uiDesiredSite)
    {
        pDeviceData[uiSiteIndex] = _makeColorVector<T>(spin, byColor);
    }
    else
    {
        pDeviceData[uiSiteIndex] = _makeZero<T>();
    }
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelMakeFullWallSourceSite(T* pDeviceData, INT uiDesiredT, BYTE spin, BYTE color, BYTE byFieldID)
{
    intokernalInt4NoConstant;

    //Since there is a 'uiShift', which set values for neighbours
    //We should not set every site to zero here!

    if (uiDesiredT == sSite4.w)
    {
        //sSite4 is no longer used
        const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldID][__bi(sSite4)];
        if (!sIdx.IsDirichlet())
        {
            pDeviceData[sIdx.m_uiSiteIndex] = _makeColorVector<T>(spin, color);
        }
    }

}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelMakeWallSourceSite(T* pDeviceData,
    INT uiDesiredT, UINT uiShift, BYTE color, BYTE byFieldID)
{
    intokernalInt4NoConstant;

    //Since there is a 'uiShift', which set values for neighbours
    //We should not set every site to zero here!

    if ((0 == (sSite4.x & 1))
     && (0 == (sSite4.y & 1))
     && (0 == (sSite4.z & 1))
     && (uiDesiredT < 0 || uiDesiredT == sSite4.w))
    {
        //sSite4 is no longer used
        sSite4.x = sSite4.x + static_cast<SCHAR>(uiShift & 1);
        sSite4.y = sSite4.y + static_cast<SCHAR>((uiShift >> 1) & 1);
        sSite4.z = sSite4.z + static_cast<SCHAR>((uiShift >> 2) & 1);
        const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldID][__bi(sSite4)];
        if (!sIdx.IsDirichlet())
        {
            pDeviceData[sIdx.m_uiSiteIndex] = _makeColorVector<T>(0, color);
        }
    }
    
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelMakeZWallSourceSite(T* pDeviceData,
    INT uiDesiredZ, UINT uiShift, BYTE color, BYTE byFieldID)
{
    intokernalInt4NoConstant;

    //Since there is a 'uiShift', which set values for neighbours
    //We should not set every site to zero here!

    if ((0 == (sSite4.x & 1))
     && (0 == (sSite4.y & 1))
     && (0 == (sSite4.w & 1))
     && (uiDesiredZ < 0 || uiDesiredZ == sSite4.z))
    {
        //sSite4 is no longer used
        sSite4.x = sSite4.x + static_cast<SCHAR>(uiShift & 1);
        sSite4.y = sSite4.y + static_cast<SCHAR>((uiShift >> 1) & 1);
        sSite4.w = sSite4.w + static_cast<SCHAR>((uiShift >> 2) & 1);
        const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldID][__bi(sSite4)];
        if (!sIdx.IsDirichlet())
        {
            pDeviceData[sIdx.m_uiSiteIndex] = _makeColorVector<T>(0, color);
        }
    }
}

/**
* f(x) phi^*(x)phi(x)
*/
template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelDiagnalTerm(
    const T* __restrict__ pDeviceData,
    T* pResultData,
    BYTE byFieldId,
    DOUBLE fCoefficient,
    _deviceCoeffFunctionPointer pfcoeff,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = _makeZero<T>();
        return;
    }

    const Real fCoefficient1 = (*pfcoeff)(byFieldId, sSite4, sIdx) * fCoefficient;
    T result = _mulC(pDeviceData[uiSiteIndex], fCoefficient1);

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    default:
        break;
    }

    _add(pResultData[uiSiteIndex], result);
}

#pragma endregion

template<typename T> 
void CCommonKernelSite<T>::InitialBuffer(T* dest, BYTE byFieldId, EFieldInitialType eInitialType)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelInitialSite<T>, block, threads, dest, byFieldId, eInitialType);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

template<typename T>
void CCommonKernelSite<T>::InitialBufferEvenOdd(T* dest, BYTE byFieldId, UBOOL bEven, UBOOL bZeroOtherSites, EFieldInitialType eInitialType)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelInitialSiteEvenOdd<T>, block, threads, dest, byFieldId, eInitialType, bEven, bZeroOtherSites);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

template<typename T>
void CCommonKernelSite<T>::FixBoundary(T* dest, BYTE byFieldId)
{
    SSmallInt4 bc = appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(byFieldId);
    if (bc.x != 0 && bc.y != 0 && bc.z != 0 && bc.w != 0)
    {
        return;
    }

    if (!appGetLattice()->HasBoundaryField(byFieldId))
    {
        appWarning(_T("Call fix boundary but without set boundary field!!\n"));
        return;
    }
    preparethread;
    _LAUNCH_KERNEL(_kernelFixBoundarySite<T>, block, threads, dest, byFieldId);
}

template<typename T>
void CCommonKernelSite<T>::InitialSource(T* data, BYTE byFieldId, const SFermionBosonSource& sourceData)
{
    
    switch (sourceData.m_eSourceType)
    {
    case EFS_Point:
    {
        preparethread;
        const UINT uiSiteIndex = _hostGetSiteIndex(sourceData.m_sSourcePoint);
        _LAUNCH_KERNEL(_kernelMakePointSourceSite<T>, block, threads, data, uiSiteIndex, sourceData.m_bySpinIndex, sourceData.m_byColorIndex);
    }
    break;
    case EFS_Wall:
    {
        preparethread;
        _LAUNCH_KERNEL(_kernelInitialSite<T>, block, threads, data, byFieldId, EFIT_Zero);
        _LAUNCH_KERNEL(_kernelMakeFullWallSourceSite<T>, block, threads,
            data,
            static_cast<INT>(sourceData.m_sSourcePoint.w),
            sourceData.m_bySpinIndex,
            sourceData.m_byColorIndex,
            byFieldId);
    }
    break;
    case EFS_StaggeredWall:
    {
        preparethread;
        _LAUNCH_KERNEL(_kernelInitialSite<T>, block, threads, data, byFieldId, EFIT_Zero);
        _LAUNCH_KERNEL(_kernelMakeWallSourceSite<T>, block, threads,
            data,
            static_cast<INT>(sourceData.m_sSourcePoint.w),
            static_cast<UINT>(sourceData.m_bySpinIndex),
            sourceData.m_byColorIndex,
            byFieldId);
    }
    break;
    case EFS_StaggeredZWall:
    {
        preparethread;
        _LAUNCH_KERNEL(_kernelInitialSite<T>, block, threads, data, byFieldId, EFIT_Zero);
        _LAUNCH_KERNEL(_kernelMakeZWallSourceSite<T>, block, threads,
            data,
            static_cast<INT>(sourceData.m_sSourcePoint.z),
            static_cast<UINT>(sourceData.m_bySpinIndex),
            sourceData.m_byColorIndex,
            byFieldId);
    }
    break;
    default:
        appCrucial(_T("The source type %s not implemented yet!\n"), __ENUM_TO_STRING(EFermionBosonSource, sourceData.m_eSourceType).c_str());
        break;
    }
}

template<typename T>
void CCommonKernelSite<T>::DebugPrint(const T* data, UINT sitecount)
{
    T* toprint = (T*)malloc(sizeof(T) * sitecount);
    checkCudaErrors(cudaMemcpy(toprint, data, sizeof(T) * sitecount, cudaMemcpyDeviceToHost));

    appPushLogDate(FALSE);
    for (UINT uiSite = 0; uiSite < sitecount; ++uiSite)
    {
        if (0 == (uiSite % _HC_Lt))
        {
            appGeneral(_T("\n"));
        }
        const SSmallInt4 site4 = __hostSiteIndexToInt4(uiSite);
        appGeneral(_T(" (%d,%d,%d,%d) = %s, "),
            site4.x, site4.y, site4.z, site4.w,
            appToString(toprint[uiSite]).c_str());
    }
    appPopLogDate();

    appSafeFree(toprint);
}

template<typename T>
void CCommonKernelSite<T>::DiagnalTerm(
    T* pTarget,
    BYTE byFieldId,
    const T* pSource, DOUBLE fCoeffiecient, _deviceCoeffFunctionPointer fpCoeff,
    EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelDiagnalTerm<T>, block, threads,
        pSource,
        pTarget,
        byFieldId,
        fCoeffiecient,
        fpCoeff,
        eOCT,
        fRealCoeff,
        cCmpCoeff
        );
}

template class CCommonKernelSite<Real>;
template class CCommonKernelSite<CLGComplex>;
template class CCommonKernelSite<deviceSU2Vector>;
template class CCommonKernelSite<deviceSU3Vector>;

#if _CLG_SU4_KS || _CLG_SU4_BOSON
template class CCommonKernelSite<deviceSU4Vector>;
#endif
#if _CLG_SU5_KS || _CLG_SU5_BOSON
template class CCommonKernelSite<deviceSU5Vector>;
#endif
#if _CLG_SU6_KS || _CLG_SU6_BOSON
template class CCommonKernelSite<deviceSU6Vector>;
#endif
#if _CLG_SU7_KS || _CLG_SU7_BOSON
template class CCommonKernelSite<deviceSU7Vector>;
#endif
#if _CLG_SU8_KS || _CLG_SU8_BOSON
template class CCommonKernelSite<deviceSU8Vector>;
#endif
//template class CCommonKernelSite<deviceWilsonVectorSU3>;

#pragma endregion

#pragma region link

#pragma region kernels

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialLink(T* pDevicePtr, BYTE byFieldId, EFieldInitialType eInitialType)
{
    const T id = _makeId<T>();
    const T zero = _makeZero<T>();

    intokernaldir;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

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
                pDevicePtr[uiLinkIndex] = _makeRandom<T>(_deviceGetLinkIndex(uiSiteIndex, idir));
            }
            break;
        case EFIT_RandomGenerator:
            {
                pDevicePtr[uiLinkIndex] = _makeGaussian<T>(_deviceGetLinkIndex(uiSiteIndex, idir));
            }
            break;
        case EFIT_SumGenerator:
            {
                pDevicePtr[uiLinkIndex] = _makeSumGenerator<T>(F(1.0));
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

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialLinkD(T* pDevicePtr, BYTE byFieldId, EFieldInitialType eInitialType)
{
    const T id = _makeId<T>();
    const T zero = _makeZero<T>();

    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const UINT uiDir = _DC_Dir;

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
            if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, idir))
            {
                if (NULL != __boundaryFieldPointers[byFieldId])
                {
                    UINT uiRegion = __idx->_devcieExchangeBoundaryFieldSiteIndexBI(byFieldId, uiBigIdx);
                    const T* buffer = ((CFieldBoundary<T>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData;
                    pDevicePtr[uiLinkIndex] = buffer[uiRegion * uiDir + idir];
                    continue;
                }
                pDevicePtr[uiLinkIndex] = id;
                continue;
            }
            pDevicePtr[uiLinkIndex] = _makeRandom<T>(_deviceGetLinkIndex(uiSiteIndex, idir));
        }
        break;
        case EFIT_RandomGenerator:
        {
            if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, idir))
            {
                if (NULL != __boundaryFieldPointers[byFieldId])
                {
                    UINT uiRegion = __idx->_devcieExchangeBoundaryFieldSiteIndexBI(byFieldId, uiBigIdx);
                    const T* buffer = ((CFieldBoundary<T>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData;
                    pDevicePtr[uiLinkIndex] = buffer[uiRegion * uiDir + idir];
                    continue;
                }
                pDevicePtr[uiLinkIndex] = id;
                continue;
            }
            pDevicePtr[uiLinkIndex] = _makeGaussian<T>(_deviceGetLinkIndex(uiSiteIndex, idir));
        }
        break;
        case EFIT_SumGenerator:
        {
            pDevicePtr[uiLinkIndex] = _makeSumGenerator<T>(F(1.0));
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

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelFixBoundaryLink(T* pDeviceData, BYTE byFieldId, UBOOL bId)
{
    intokernalDirInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    const T id = _makeId<T>();
    const T zero = _makeZero<T>();

    if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
    {
        if (NULL != __boundaryFieldPointers[byFieldId])
        {
            UINT uiRegion = __idx->_devcieExchangeBoundaryFieldSiteIndexBI(byFieldId, uiBigIdx);
            const T* buffer = ((CFieldBoundary<T>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData;
            pDeviceData[uiSiteIndex] = buffer[uiRegion * _DC_Dir + dir];
            return;
        }
        pDeviceData[uiSiteIndex] = bId ? id : zero;
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStrictExp(deviceGauge* pDeviceData, BYTE byFieldId)
{
    __gaugeKernel(pDeviceData[uiLinkIndex] = _strictexp(pDeviceData[uiLinkIndex]))
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelStrictLog(deviceGauge* pDeviceData, BYTE byFieldId)
{
    __gaugeKernel(pDeviceData[uiLinkIndex] = _strictlog(pDeviceData[uiLinkIndex]))
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelQuickExp(deviceGauge* pDeviceData, BYTE byFieldId)
{
    __gaugeKernel(_expreal(pDeviceData[uiLinkIndex], F(1.0)))
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelQuickLog(deviceGauge* pDeviceData, BYTE byFieldId)
{
    __gaugeKernel(_ta(pDeviceData[uiLinkIndex]))
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelOneDirId(deviceGauge* pDeviceData, BYTE byFieldId, BYTE byDir)
{
    __gaugeKernel(
        if (0 != ((1 << dir) & byDir))
        {
            pDeviceData[uiLinkIndex] = _makeId<deviceGauge>();
        }
    )
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelOneDirZero(deviceGauge* pDeviceData, BYTE byFieldId, BYTE byDir)
{
    __gaugeKernel(
        if (0 != ((1 << dir) & byDir))
        {
            pDeviceData[uiLinkIndex] = _makeZero<deviceGauge>();
        }
    )
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelExp(deviceGauge* pTarget, const deviceGauge* __restrict__ pSource, BYTE byFieldId, Real a)
{
    __gaugeKernel(
        deviceGauge expP = _expreal(pSource[uiLinkIndex], a);
        _mul(expP, pTarget[uiLinkIndex]);
        pTarget[uiLinkIndex] = expP;
    )
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovLoopOfSite(
    const deviceGauge* __restrict__ pDeviceBuffer,
    cuDoubleComplex* res,
    BYTE byFieldId)
{
    intokernalInt4_S(0);
    UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, _DC_Dir - 1);
    SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

    deviceGauge tmp = _makeZero<deviceGauge>();
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
            _mul(tmp, pDeviceBuffer[uiLinkIdx]);
        }
    }

    res[uiSiteIndex3D] = _cToDouble(_tr(tmp));
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovLoopOfSiteX(
    const deviceGauge* __restrict__ pDeviceBuffer,
    cuDoubleComplex* res,
    BYTE byFieldId)
{
    intokernalInt4_Syzt(0);
    UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, 0);
    SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

    deviceGauge tmp = _makeZero<deviceGauge>();
    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 0))
    {
        tmp = pDeviceBuffer[uiLinkIdx];
    }

    for (UINT uiX = 1; uiX < _DC_Lx; ++uiX)
    {
        UINT newSiteIndex = uiSiteIndex + uiX * _DC_MultX;
        uiLinkIdx = _deviceGetLinkIndex(newSiteIndex, 0);
        site4 = __deviceSiteIndexToInt4(newSiteIndex);
        uiBigIdx = __idx->_deviceGetBigIndex(site4);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 0))
        {
            _mul(tmp, pDeviceBuffer[uiLinkIdx]);
        }
    }

    res[uiSiteIndex3DYZT] = _cToDouble(_tr(tmp));
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovLoopOfSiteY(
    const deviceGauge* __restrict__ pDeviceBuffer,
    cuDoubleComplex* res,
    BYTE byFieldId)
{
    intokernalInt4_Sxzt(0);
    UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, 1);
    SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

    deviceGauge tmp = _makeZero<deviceGauge>();
    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 1))
    {
        tmp = pDeviceBuffer[uiLinkIdx];
    }

    for (UINT uiY = 1; uiY < _DC_Ly; ++uiY)
    {
        UINT newSiteIndex = uiSiteIndex + uiY * _DC_MultY;
        uiLinkIdx = _deviceGetLinkIndex(newSiteIndex, 1);
        site4 = __deviceSiteIndexToInt4(newSiteIndex);
        uiBigIdx = __idx->_deviceGetBigIndex(site4);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 1))
        {
            _mul(tmp, pDeviceBuffer[uiLinkIdx]);
        }
    }

    res[uiSiteIndex3DXZT] = _cToDouble(_tr(tmp));
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovLoopOfSiteZ(
    const deviceGauge* __restrict__ pDeviceBuffer,
    cuDoubleComplex* res,
    BYTE byFieldId)
{
    intokernalInt4_Sxyt(0);
    UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, 2);
    SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

    deviceGauge tmp = _makeZero<deviceGauge>();
    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 2))
    {
        tmp = pDeviceBuffer[uiLinkIdx];
    }

    for (UINT uiZ = 1; uiZ < _DC_Lz; ++uiZ)
    {
        UINT newSiteIndex = uiSiteIndex + uiZ * _DC_MultZ;
        uiLinkIdx = _deviceGetLinkIndex(newSiteIndex, 2);
        site4 = __deviceSiteIndexToInt4(newSiteIndex);
        uiBigIdx = __idx->_deviceGetBigIndex(site4);

        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, 2))
        {
            _mul(tmp, pDeviceBuffer[uiLinkIdx]);
        }
    }

    res[uiSiteIndex3DXYT] = _cToDouble(_tr(tmp));
}

/**
 * E_mu = F_{0 mu}
 */
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelTransformToE(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    deviceGauge* pRes)
{
    intokernalDirInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    deviceGauge res = _makeZero<deviceGauge>();
    if (dir < uiDir - 1)
    {
        //find clover F
        res = _device1PlaqutteTermPPT(pDeviceData, 3, dir, uiBigIdx, sSite4, byFieldId);
        _ta(res);
        _mul(res, F(-1.0));
    }

    pRes[uiLinkIndex] = res;
}

/**
 * This is wrong! the order of the plaqutte must be considered
 * This is to make sure gauge transform is g(x) nabla E g^+(n)
 */
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateNablaE(
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

        SCHAR dirs[4];
        //deviceGauge toMul(_devicePlaqutte(pDeviceData, pCachedPlaqutte, uiSiteIndex, byPlaqIdx, plaqLength, plaqCount));
        dirs[0] = 4;
        dirs[1] = dir + 1;
        dirs[2] = -4;
        dirs[3] = -static_cast<SCHAR>(dir) - 1;
        deviceGauge toMul(
            //_device1PlaqutteTermPP(pDeviceData, 3, dir, uiBigIdx)
            _deviceLinkT(pDeviceData, sSite4, 4, byFieldId, dirs)
        );

        //
        dirs[0] = 4;
        dirs[1] = -static_cast<SCHAR>(dir) - 1;
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
_kernelCalculateNablaENaive(
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

        SCHAR dirs[4];
        //deviceGauge toMul(_devicePlaqutte(pDeviceData, pCachedPlaqutte, uiSiteIndex, byPlaqIdx, plaqLength, plaqCount));
        dirs[0] = 4;
        dirs[1] = dir + 1;
        dirs[2] = -4;
        dirs[3] = -static_cast<SCHAR>(dir) - 1;

        deviceGauge a(
            //_device1PlaqutteTermPP(pDeviceData, 3, dir, uiBigIdx)
            _deviceLinkT(pDeviceData, sSite4, 4, byFieldId, dirs)
            //_deviceClover(pDeviceData, sSite4, __bi(sSite4), 3, dir, byFieldId)
        );

        dirs[0] = -static_cast<SCHAR>(dir) - 1;
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


template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateKinematicEnergyT_D(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pDeviceData,
    DOUBLE* results
)
{
    intokernalDirInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    DOUBLE resThisThread = __idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir) ? 0.0 
        : static_cast<DOUBLE>(_retr(_dagmulC(pDeviceData[uiLinkIndex], pDeviceData[uiLinkIndex])));

    if (0 == dir)
    {
        results[uiSiteIndex] = resThisThread;
    }
    __syncthreads();

    if (1 == dir)
    {
        results[uiSiteIndex] += resThisThread;
    }
    __syncthreads();

    if (2 == dir)
    {
        results[uiSiteIndex] += resThisThread;
    }
    __syncthreads();

    if (3 == dir)
    {
        results[uiSiteIndex] += resThisThread;
    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeAddLink(
    const deviceGauge* __restrict__ pSource,
    deviceGauge* pTarget,
    const SCHAR* __restrict__ devicePath,
    BYTE byPathLen,
    BYTE byMu,
    BYTE byGaugeFieldId,
    Real fCoeff
)
{
    intokernalInt4;
    deviceGauge vn = _deviceLinkT(pSource, sSite4, byPathLen, byGaugeFieldId, devicePath);
    _mul(vn, fCoeff);
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byMu);
    _add(pTarget[linkIndex], vn);
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelApplyStaggeredPhase(deviceGauge* pDeviceData, const BYTE* __restrict__ etatable, BYTE byFieldId)
{
    __gaugeKernel(
        if ((etatable[uiSiteIndex] >> dir) & 1)
        {
            _oppo(pDeviceData[uiLinkIndex]);
        }
    )
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelApplyStaggeredPhaseMILC(deviceGauge* pDeviceData, BYTE byFieldId)
{
    intokernalDirInt4;
    SCHAR coordlst[4] = { sSite4.w, sSite4.x, sSite4.y, sSite4.z };
    BYTE shiftedDir = (dir + 1U) & 3U;
    BYTE parity = 0U;
    for (BYTE i = 0; i < shiftedDir; ++i)
    {
        parity += coordlst[i];
    }
    if (parity & 1U)
    {
        _oppo(pDeviceData[uiLinkIndex]);
    }
}

#pragma endregion

template<typename T>
void CCommonKernelLink<T>::InitialBuffer(T* dest, BYTE byFieldId, EFieldInitialType eInitialType)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelInitialLink<T>, block, threads, dest, byFieldId, eInitialType);
    _CHECKCUDA;
}

template<typename T>
void CCommonKernelLink<T>::InitialBufferD(T* dest, BYTE byFieldId, EFieldInitialType eInitialType)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelInitialLinkD<T>, block, threads, dest, byFieldId, eInitialType);
    _CHECKCUDA;
}

template<typename T>
void CCommonKernelLink<T>::FixBoundary(T* dest, BYTE byFieldId)
{
    SSmallInt4 bc = appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(byFieldId);
    if (bc.x != 0 && bc.y != 0 && bc.z != 0 && bc.w != 0)
    {
        return;
    }

    if (!appGetLattice()->HasBoundaryField(byFieldId))
    {
        appWarning(_T("Call fix boundary but without set boundary field!!\n"));
        return;
    }
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelFixBoundaryLink<T>, block, threads, dest, byFieldId, TRUE);
}

template<typename T>
void CCommonKernelLink<T>::FixBoundaryZero(T* dest, BYTE byFieldId)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelFixBoundaryLink<T>, block, threads, dest, byFieldId, FALSE);
}

template<typename T>
void CCommonKernelLink<T>::DebugPrint(const T* data, UINT uiLinkCount)
{
    //preparethread;
    //_LAUNCH_KERNEL(_kernelPrintSU3, block, threads, m_pDeviceData);

    //===================================================
    //Since Debug Print Me is only used to debug, we do it slow but convinient
    T* pToPrint = (T*)malloc(sizeof(T) * uiLinkCount);
    checkCudaErrors(cudaMemcpy(pToPrint, data, sizeof(T) * uiLinkCount, cudaMemcpyDeviceToHost));

    for (UINT uiLink = 0; uiLink < uiLinkCount; ++uiLink)
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

template<typename T>
void CCommonKernelLink<T>::ExpMul(T* other, BYTE byFieldId, const T* me, Real a)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelExp<T>, block, threads, other, me, byFieldId, a);
}

template<typename T>
void CCommonKernelLink<T>::QuickLog(T* data, BYTE byFieldId)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelQuickLog<T>, block, threads, data, byFieldId);
}

template<typename T>
void CCommonKernelLink<T>::QuickExp(T* data, BYTE byFieldId)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelQuickExp<T>, block, threads, data, byFieldId);
}

template<typename T>
void CCommonKernelLink<T>::StrictLog(T* data, BYTE byFieldId)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelStrictLog<T>, block, threads, data, byFieldId);
}

template<typename T>
void CCommonKernelLink<T>::StrictExp(T* data, BYTE byFieldId)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelStrictExp<T>, block, threads, data, byFieldId);
}

template<typename T>
void CCommonKernelLink<T>::SetOneDirectionUnity(T* data, BYTE byFieldId, BYTE byDir)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelOneDirId<T>, block, threads, data, byFieldId, byDir);
}

template<typename T>
void CCommonKernelLink<T>::SetOneDirectionZero(T* data, BYTE byFieldId, BYTE byDir)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelOneDirZero<T>, block, threads, data, byFieldId, byDir);
}

template<typename T>
void CCommonKernelLink<T>::PolyakovOnSpatialSite(const T* data, BYTE byFieldId, cuDoubleComplex* buffer, BYTE byDir)
{
    if (0 == byDir)
    {
        preparethread_Syzt;
        _LAUNCH_KERNEL(_kernelPolyakovLoopOfSiteX<T>, block3dyzt, threads3dyzt, data, buffer, byFieldId);
    }
    else if (1 == byDir)
    {
        preparethread_Sxzt;
        _LAUNCH_KERNEL(_kernelPolyakovLoopOfSiteY<T>, block3dxzt, threads3dxzt, data, buffer, byFieldId);
    }
    else if (2 == byDir)
    {
        preparethread_Sxyt;
        _LAUNCH_KERNEL(_kernelPolyakovLoopOfSiteZ<T>, block3dxyt, threads3dxyt, data, buffer, byFieldId);
    }
    else
    {
        preparethread_S;
        _LAUNCH_KERNEL(_kernelPolyakovLoopOfSite<T>, block3d, threads3d, data, buffer, byFieldId);
    }
}

template<typename T>
void CCommonKernelLink<T>::CalculateE_Using_U(const T* deviceData, BYTE byFieldId, T* pResoult)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelTransformToE<T>, block, threads, byFieldId, deviceData, pResoult);
}

template<typename T>
void CCommonKernelLink<T>::CalculateNablaE_Using_U(const T* deviceData, BYTE byFieldId, T* pResoult, UBOOL bNaive)
{
    preparethread;
    if (bNaive)
    {
        _LAUNCH_KERNEL(_kernelCalculateNablaENaive<T>, block, threads,
            deviceData,
            byFieldId,
            pResoult);
    }
    else
    {
        _LAUNCH_KERNEL(_kernelCalculateNablaE<T>, block, threads,
            deviceData,
            byFieldId,
            pResoult);
    }
}

template<typename T>
DOUBLE CCommonKernelLink<T>::CalcKineticEnery(const T* me, BYTE byFieldId)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelCalculateKinematicEnergyT_D<T>, block, threads, byFieldId, me, _D_RealThreadBuffer);
    return _clgGlobalThreadBufferSum(_D_RealThreadBuffer);
}

template<typename T>
void CCommonKernelLink<T>::AddLink(const T* source, T* target, Real fCoeff, const SCHAR* devicePath, BYTE byPathLen, BYTE byMu, BYTE byFieldId)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelGaugeAddLink<T>, block, threads,
        source,
        target,
        devicePath,
        byPathLen,
        byMu,
        byFieldId,
        fCoeff
        );
}

template<typename T>
void CCommonKernelLink<T>::ApplyStaggeredPhase(T* deviceData, BYTE byFieldId)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelApplyStaggeredPhase<T>, block, threads, deviceData, appGetLattice()->m_pIndexCache->m_pEtaMu, byFieldId);
}

template<typename T>
void CCommonKernelLink<T>::ApplyStaggeredPhaseMILC(T* deviceData, BYTE byFieldId)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelApplyStaggeredPhaseMILC<T>, block, threads, deviceData, byFieldId);
}

template class CCommonKernelLink<CLGComplex>;
template class CCommonKernelLink<deviceSU2>;
template class CCommonKernelLink<deviceSU3>;

#if _CLG_SU4_GAUGE
template class CCommonKernelLink<deviceSU4>;
#endif
#if _CLG_SU5_GAUGE
template class CCommonKernelLink<deviceSU5>;
#endif
#if _CLG_SU6_GAUGE
template class CCommonKernelLink<deviceSU6>;
#endif
#if _CLG_SU7_GAUGE
template class CCommonKernelLink<deviceSU7>;
#endif
#if _CLG_SU8_GAUGE
template class CCommonKernelLink<deviceSU8>;
#endif

#if _CLG_Z2_GAUGE
template class CCommonKernelLink<deviceZN<2>>;
#endif
#if _CLG_Z3_GAUGE
template class CCommonKernelLink<deviceZN<3>>;
#endif
#if _CLG_Z4_GAUGE
template class CCommonKernelLink<deviceZN<4>>;
#endif
#if _CLG_Z5_GAUGE
template class CCommonKernelLink<deviceZN<5>>;
#endif
#if _CLG_Z6_GAUGE
template class CCommonKernelLink<deviceZN<6>>;
#endif
#if _CLG_D3_GAUGE
template class CCommonKernelLink<deviceDN<3>>;
#endif
#if _CLG_D4_GAUGE
template class CCommonKernelLink<deviceDN<4>>;
#endif
#if _CLG_D8_GAUGE
template class CCommonKernelLink<deviceDN<8>>;
#endif

#if _CLG_SL3C_GAUGE
template class CCommonKernelLink<deviceSL3C>;
#endif
#if _CLG_U3_GAUGE
template class CCommonKernelLink<deviceU3>;
#endif
#if _CLG_O3_GAUGE
template class CCommonKernelLink<deviceO3>;
#endif
#if _CLG_SO3_GAUGE
template class CCommonKernelLink<deviceSO3>;
#endif

#pragma endregion

#pragma region MV

#pragma region kernels

template<typename T, typename gaugetype>
__global__ void _CLG_LAUNCH_BOUND
_kernelConnectionOneField(const T* __restrict__ n, gaugetype* res, const SIndex* __restrict__ move)
{
    intokernalDir_NoDir;
    res[uiLinkIndex] = _makeContract<gaugetype, T>(n[uiSiteIndex], n[move[2 * uiLinkIndex].m_uiSiteIndex]);
}

/**
* n_p_m.n^+
*/
template<typename T, typename gaugetype>
__global__ void _CLG_LAUNCH_BOUND
_kernelConnection(const T* __restrict__ n, const T* __restrict__ n_p_m, gaugetype* res, const SIndex* __restrict__ move)
{
    intokernalDir_NoDir;
    res[uiLinkIndex] = _makeContract<gaugetype, T>(n[uiSiteIndex], n_p_m[move[2 * uiLinkIndex].m_uiSiteIndex]);
}

template<typename T, typename gaugetype>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddConnection(const T* __restrict__ n, const T* __restrict__ n_p_m, gaugetype* res, const SIndex* __restrict__ move, Real fCoeff)
{
    intokernalDir_NoDir;
    gaugetype res2 = _makeContract<gaugetype, T>(n[uiSiteIndex], n_p_m[move[2 * uiLinkIndex].m_uiSiteIndex]);
    _mul(res2, fCoeff);
    _add(res[uiLinkIndex], res2);
}

template<typename T, typename gaugetype>
__global__ void _CLG_LAUNCH_BOUND
_kernelConnectionOneFieldStaggered(
    Real fCoeff,
    const T* __restrict__ n, 
    const BYTE* __restrict__ pEtaTable,
    gaugetype* res, 
    const SIndex* __restrict__ move)
{
    intokernalDirInt4;
    
    const SIndex& x_p_mu_Fermion = move[2 * uiLinkIndex];
    res[uiLinkIndex] = _makeContract<gaugetype, T>(n[x_p_mu_Fermion.m_uiSiteIndex], n[uiSiteIndex]);
    _mul(res[uiLinkIndex], fCoeff);
    BYTE eta = pEtaTable[uiSiteIndex] >> dir;
    if (sSite4.IsOdd())
    {
        eta = eta + 1;
    }
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        eta = eta + 1;
    }

    if (eta & 1)
    {
        _oppo(res[uiLinkIndex]);
    }
}

template<typename T, typename gaugetype>
__global__ void _CLG_LAUNCH_BOUND
_kernelAddConnectionOneFieldStaggered(
    Real fCoeff,
    const T* __restrict__ n,
    const BYTE* __restrict__ pEtaTable,
    gaugetype* res,
    const SIndex* __restrict__ move)
{
    intokernalDirInt4;

    const SIndex& x_p_mu_Fermion = move[2 * uiLinkIndex];
    BYTE eta = pEtaTable[uiSiteIndex] >> dir;
    if (sSite4.IsOdd())
    {
        eta = eta + 1;
    }
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        eta = eta + 1;
    }
    gaugetype toadd = _makeContract<gaugetype, T>(n[x_p_mu_Fermion.m_uiSiteIndex], n[uiSiteIndex]);
    _mul(toadd, fCoeff);
    if (eta & 1)
    {
        _sub(res[uiLinkIndex], toadd);
    }
    else
    {
        _add(res[uiLinkIndex], toadd);
    }
}

//template<typename T, typename gaugetype>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelConnectionOneFieldStaggeredTest(
//    Real fCoeff,
//    const T* __restrict__ n,
//    const BYTE* __restrict__ pEtaTable,
//    gaugetype* res,
//    const SIndex* __restrict__ move)
//{
//    intokernalDirInt4;
//
//    const SIndex& x_p_mu_Fermion = move[2 * uiLinkIndex];
//    res[uiLinkIndex] = _makeContract<gaugetype, T>(n[x_p_mu_Fermion.m_uiSiteIndex], n[uiSiteIndex]);
//    _mul(res[uiLinkIndex], fCoeff);
//    BYTE eta = pEtaTable[uiSiteIndex] >> dir;
//    if (sSite4.IsOdd())
//    {
//        eta = eta + 1;
//    }
//    if (x_p_mu_Fermion.NeedToOpposite())
//    {
//        eta = eta + 1;
//    }
//
//    if (eta & 1)
//    {
//        _oppo(res[uiLinkIndex]);
//    }
//}
//
//template<typename T, typename gaugetype>
//__global__ void _CLG_LAUNCH_BOUND
//_kernelAddConnectionOneFieldStaggeredTest(
//    Real fCoeff,
//    const T* __restrict__ n,
//    const BYTE* __restrict__ pEtaTable,
//    gaugetype* res,
//    const SIndex* __restrict__ move)
//{
//    intokernalDirInt4;
//
//    const SIndex& x_p_mu_Fermion = move[2 * uiLinkIndex];
//    BYTE eta = pEtaTable[uiSiteIndex] >> dir;
//    if (sSite4.IsOdd())
//    {
//        eta = eta + 1;
//    }
//    if (x_p_mu_Fermion.NeedToOpposite())
//    {
//        eta = eta + 1;
//    }
//    gaugetype toadd = _makeContract<gaugetype, T>(n[x_p_mu_Fermion.m_uiSiteIndex], n[uiSiteIndex]);
//    _mul(toadd, fCoeff);
//    if (eta & 1)
//    {
//        _sub(res[uiLinkIndex], toadd);
//    }
//    else
//    {
//        _add(res[uiLinkIndex], toadd);
//    }
//}

#pragma endregion

template<typename vector, typename matrix, INT vectorN>
void CCommonKernelMV<vector, matrix, vectorN>::ConnectionOneField(const vector* v, matrix* res, BYTE byFieldId)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelConnectionOneField TMPARG(vector, matrix), block, threads, v, res, appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId]);
}

template<typename vector, typename matrix, INT vectorN>
void CCommonKernelMV<vector, matrix, vectorN>::ConnectionOneFieldStaggered(const vector* v, matrix* res, BYTE byFieldId, Real fCoeff)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelConnectionOneFieldStaggered TMPARG(vector, matrix), block, threads,
        fCoeff,
        v, 
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        res, 
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId]);
}

template<typename vector, typename matrix, INT vectorN>
void CCommonKernelMV<vector, matrix, vectorN>::AddConnectionOneFieldStaggered(const vector* v, matrix* res, BYTE byFieldId, Real fCoeff)
{
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelAddConnectionOneFieldStaggered TMPARG(vector, matrix), block, threads,
        fCoeff,
        v,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        res,
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId]);
}

//template<typename vector, typename matrix, INT vectorN>
//void CCommonKernelMV<vector, matrix, vectorN>::ConnectionOneFieldStaggeredTest(const vector* v, matrix* res, BYTE byFieldId, Real fCoeff)
//{
//    preparethreadDir;
//    _LAUNCH_KERNEL(_kernelConnectionOneFieldStaggeredTest, block, threads, 
//        fCoeff,
//        v,
//        appGetLattice()->m_pIndexCache->m_pEtaMu,
//        res,
//        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId]);
//}
//
//template<typename vector, typename matrix, INT vectorN>
//void CCommonKernelMV<vector, matrix, vectorN>::AddConnectionOneFieldStaggeredTest(const vector* v, matrix* res, BYTE byFieldId, Real fCoeff)
//{
//    preparethreadDir;
//    _LAUNCH_KERNEL(_kernelAddConnectionOneFieldStaggeredTest, block, threads, 
//        fCoeff,
//        v,
//        appGetLattice()->m_pIndexCache->m_pEtaMu,
//        res,
//        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId]);
//}

/**
* n_p_m.n^+
* Note: This function is wrong and should never be called
* In CLGLib, we always calculate f0^+ for simplicity, which should be n.n_p_m^+
*/
template<typename vector, typename matrix, INT vectorN>
void CCommonKernelMV<vector, matrix, vectorN>::Connection(const vector* n, const vector* n_p_m, matrix* res, BYTE byFieldId)
{
    appGeneral(_T("CCommonKernelMV<vector, matrix, vectorN>::Connection should not be here.\n"));
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelConnection TMPARG(vector, matrix), block, threads, n, n_p_m, res, appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId]);
}

/**
* n_p_m.n^+
* Note: This function is wrong and should never be called
* In CLGLib, we always calculate f0^+ for simplicity, which should be n.n_p_m^+
*/
template<typename vector, typename matrix, INT vectorN>
void CCommonKernelMV<vector, matrix, vectorN>::AddConnection(const vector* n, const vector* n_p_m, matrix* res, Real fCoeff, BYTE byFieldId)
{
    appGeneral(_T("CCommonKernelMV<vector, matrix, vectorN>::AddConnection should not be here.\n"));
    preparethreadDir;
    _LAUNCH_KERNEL(_kernelAddConnection TMPARG(vector, matrix), block, threads, n, n_p_m, res, appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId], fCoeff);
}

template class CCommonKernelMV<CLGComplex, CLGComplex, 1>;
template class CCommonKernelMV<deviceSU2Vector, deviceSU2, 2>;
template class CCommonKernelMV<deviceSU3Vector, deviceSU3, 3>;

#if _CLG_SU4_KS || _CLG_SU4_BOSON
template class CCommonKernelMV<deviceSU4Vector, deviceSU4, 4>;
#endif
#if _CLG_SU5_KS || _CLG_SU5_BOSON
template class CCommonKernelMV<deviceSU5Vector, deviceSU5, 5>;
#endif
#if _CLG_SU6_KS || _CLG_SU6_BOSON
template class CCommonKernelMV<deviceSU6Vector, deviceSU6, 6>;
#endif
#if _CLG_SU7_KS || _CLG_SU7_BOSON
template class CCommonKernelMV<deviceSU7Vector, deviceSU7, 7>;
#endif
#if _CLG_SU8_KS || _CLG_SU8_BOSON
template class CCommonKernelMV<deviceSU8Vector, deviceSU8, 8>;
#endif

#pragma endregion


__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================