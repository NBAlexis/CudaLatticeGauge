//=============================================================================
// FILENAME : CFieldBosonvn.h
// 
// DESCRIPTION:
// This is the class for the spin fields
//
// REVISION:
//  [07/04/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldCommonKernel.h"

__BEGIN_NAMESPACE

void CLGAPI appSimpleMalloc(void** ptr, size_t size)
{
    checkCudaErrors(cudaMalloc(ptr, size));
}

void CLGAPI appSimpleFree(void* ptr)
{
    checkCudaErrors(cudaFree(ptr));
}

void CLGAPI appSimpleCopyDD(void* target, void* source, size_t size)
{
    checkCudaErrors(cudaMemcpy(target, source, size, cudaMemcpyDeviceToDevice));
}

void CLGAPI appSimpleCopyHD(void* target, void* source, size_t size)
{
    checkCudaErrors(cudaMemcpy(target, source, size, cudaMemcpyHostToDevice));
}

#pragma region common

template<typename T>
void CCommonKernel<T>::AllocateBuffer(T** pointer, UINT count)
{
    checkCudaErrors(__cudaMalloc((void**)pointer, sizeof(T) * count));
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

template<typename T>
void CCommonKernel<T>::FreeBuffer(T** pointer)
{
    checkCudaErrors(__cudaFree(*pointer));
    *pointer = NULL;
}

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
template class CCommonKernel<deviceSU4>;
template class CCommonKernel<deviceSU5>;
template class CCommonKernel<deviceSU6>;
template class CCommonKernel<deviceSU7>;
template class CCommonKernel<deviceSU8>;
template class CCommonKernel<deviceSU2Vector>;
template class CCommonKernel<deviceSU3Vector>;
template class CCommonKernel<deviceSU4Vector>;
template class CCommonKernel<deviceSU5Vector>;
template class CCommonKernel<deviceSU6Vector>;
template class CCommonKernel<deviceSU7Vector>;
template class CCommonKernel<deviceSU8Vector>;
template class CCommonKernel<deviceWilsonVectorSU3>;

#pragma endregion

#pragma region field

#pragma region kernel

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialCommon(T* pDevicePtr, UINT count, EFieldInitialType eInitialType)
{
    __simplekernel(
        const UINT idx = (threadIdx.x + blockIdx.x * blockDim.x);
        if (idx < count)
        {
            switch (eInitialType)
            {
            case EFIT_Zero:
            {
                pDevicePtr[idx] = _makeZero<T>();
            }
            break;
            case EFIT_Identity:
            {
                pDevicePtr[idx] = _makeId<T>();
            }
            break;
            default:
            {
                printf("Field cannot be initialized with this type! %d\n", eInitialType);
            }
            break;
            }
        }
    )
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
_kernelMul(T* pMe, const T* __restrict__ pOther, UBOOL bConj, UINT count)
{
    __simplekernel(
        if (bConj)
        {
            _dagger(pMe[idx]);
        }
        _mul(pMe[idx], pOther[idx]);
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
_kernelDot(const T* __restrict__ pMe, const T* __restrict__ pOther, cuDoubleComplex* result, UINT count)
{
    __simplekernel(result[idx] = _cToDouble(_dot(pMe[idx], pOther[idx])))
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelDotDir(const T* __restrict__ pMe, const T* __restrict__ pOther, cuDoubleComplex* result, UINT count)
{
    __simplekernel(
        const UINT uiDir = _DC_Dir;
        result[idx] = _cToDouble(_dot(pMe[idx * uiDir], pOther[idx * uiDir]));
        for (UINT idir = 1; idir < uiDir; ++idir)
        {
            result[idx] = cuCadd(result[idx], _cToDouble(_dot(pMe[idx * uiDir + idir], pOther[idx * uiDir + idir])));
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
        const UINT uiDir = _DC_Dir;
        result[idx] = static_cast<DOUBLE>(_lensq(pMe[idx * uiDir]));
        for (UINT idir = 1; idir < uiDir; ++idir)
        {
            result[idx] += static_cast<DOUBLE>(_lensq(pMe[idx * uiDir + idir]));
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
        const UINT uiDir = _DC_Dir;
        result[idx] = static_cast<DOUBLE>(_element(pMe[idx * uiDir], elementidx));
        for (UINT idir = 1; idir < uiDir; ++idir)
        {
            result[idx] += static_cast<DOUBLE>(_element(pMe[idx * uiDir + idir], elementidx));
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
    __SIMPLEDECOMPOSE(count);
    _kernelInitialCommon << <block, thread >> > (pointer, count, eInitialType);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
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
        assert(_elementdim<T>() <= CCString::_CLG_MAX_PATH);
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
    _kernelConjugate << <block, thread >> > (dest, count);
}

template<typename T>
void CCommonKernelField<T>::AxpyPlus(T* dest, UINT count, const T* x)
{
    __SIMPLEDECOMPOSE(count);
    _kernelAdd << <block, thread >> > (dest, x, count);
}

template<typename T>
void CCommonKernelField<T>::AxpyMinus(T* dest, UINT count, const T* x)
{
    __SIMPLEDECOMPOSE(count);
    _kernelSub << <block, thread >> > (dest, x, count);
}

template<typename T>
void CCommonKernelField<T>::Axpy(T* dest, UINT count, Real a, const T* x)
{
    __SIMPLEDECOMPOSE(count);
    _kernelAxpyReal << <block, thread >> > (dest, x, a, count);
}

template<typename T>
void CCommonKernelField<T>::Axpy(T* dest, UINT count, const CLGComplex& a, const T* x)
{
    __SIMPLEDECOMPOSE(count);
    _kernelAxpy << <block, thread >> > (dest, x, a, count);
}

template<typename T>
void CCommonKernelField<T>::Mul(T* dest, UINT count, const T* x, UBOOL bDagger)
{
    __SIMPLEDECOMPOSE(count);
    _kernelMul << <block, thread >> > (dest, x, bDagger, count);
}

template<typename T>
void CCommonKernelField<T>::ScalarMultply(T* dest, UINT count, const CLGComplex& a)
{
    __SIMPLEDECOMPOSE(count);
    _kernelScalarMultiply << <block, thread >> > (dest, a, count);
}

template<typename T>
void CCommonKernelField<T>::ScalarMultply(T* dest, UINT count, Real a)
{
    __SIMPLEDECOMPOSE(count);
    _kernelScalarMultiplyReal << <block, thread >> > (dest, a, count);
}

template<typename T>
cuDoubleComplex CCommonKernelField<T>::Dot(const T* me, UINT count, const T* other)
{
    if (count == _HC_Volume)
    {
        __SIMPLEDECOMPOSE(count);
        _kernelDot << <block, thread >> > (me, other, _D_ComplexThreadBuffer, count);
        return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
    }
    else if (count == _HC_Volume * _HC_Dir)
    {
        count = count / _HC_Dir;
        __SIMPLEDECOMPOSE(count);
        _kernelDotDir << <block, thread >> > (me, other, _D_ComplexThreadBuffer, count);
        return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
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
        _kernelLengthSq << <block, thread >> > (me, _D_RealThreadBuffer, count);
        return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
    }
    else if (count == _HC_Volume * _HC_Dir)
    {
        count = count / _HC_Dir;
        __SIMPLEDECOMPOSE(count);
        _kernelLengthSqDir << <block, thread >> > (me, _D_RealThreadBuffer, count);
        return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
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
            _kernelElement << <block, thread >> > (me, i, _D_RealThreadBuffer, count);
            ret.AddItem(appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer));
        }
        return ret;
    }
    else if (count == _HC_Volume * _HC_Dir)
    {
        count = count / _HC_Dir;
        __SIMPLEDECOMPOSE(count);
        TArray<DOUBLE> ret;
        for (UINT i = 0; i < _elementdim<T>(); ++i)
        {
            _kernelElementDir << <block, thread >> > (me, i, _D_RealThreadBuffer, count);
            ret.AddItem(appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer));
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
    _kernelNorm << <block, thread >> > (dest, count);
}

template class CCommonKernelField<CLGComplex>;
template class CCommonKernelField<deviceSU2>;
template class CCommonKernelField<deviceSU3>;
template class CCommonKernelField<deviceSU4>;
template class CCommonKernelField<deviceSU5>;
template class CCommonKernelField<deviceSU6>;
template class CCommonKernelField<deviceSU7>;
template class CCommonKernelField<deviceSU8>;
template class CCommonKernelField<deviceSU2Vector>;
template class CCommonKernelField<deviceSU3Vector>;
template class CCommonKernelField<deviceSU4Vector>;
template class CCommonKernelField<deviceSU5Vector>;
template class CCommonKernelField<deviceSU6Vector>;
template class CCommonKernelField<deviceSU7Vector>;
template class CCommonKernelField<deviceSU8Vector>;
template class CCommonKernelField<deviceWilsonVectorSU3>;

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
        pDevicePtr[uiSiteIndex] = _makeGaussian<T>(_deviceGetFatIndex(uiSiteIndex, 0));
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
        pDevicePtr[uiSiteIndex] = _makeRandom<T>(_deviceGetFatIndex(uiSiteIndex, 0));
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
        pDevicePtr[uiSiteIndex] = _makeZ4<T>(_deviceGetFatIndex(uiSiteIndex, 0));
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
            UINT uiRegion = __idx->_devcieExchangeBoundaryFieldSiteIndexBI(byFieldId, uiBigIndex);
            const T* buffer = ((CFieldBoundary<T>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData;
            pDeviceData[uiSiteIndex] = buffer[uiRegion];
        }
        pDeviceData[uiSiteIndex] = _makeZero<T>();
        return;
    }
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelMakePointSourceSite(T* pDeviceData, UINT uiDesiredSite, BYTE byColor)
{
    intokernal;
    if (uiSiteIndex == uiDesiredSite)
    {
        pDeviceData[uiSiteIndex] = _makeColorVector<T>(byColor);
    }
    else
    {
        pDeviceData[uiSiteIndex] = _makeZero<T>();
    }
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelMakeWallSourceSite(T* pDeviceData,
    INT uiDesiredT, UINT uiShift, BYTE color, BYTE byFieldID)
{
    intokernalOnlyInt4;

    //Since there is a 'uiShift', which set values for neighbours
    //We should not set every site to zero here!

    if ((0 == (sSite4.x & 1))
     && (0 == (sSite4.y & 1))
     && (0 == (sSite4.z & 1))
     && (uiDesiredT < 0 || uiDesiredT == sSite4.w))
    {
        //sSite4 is no longer used
        sSite4.x = sSite4.x + static_cast<SBYTE>(uiShift & 1);
        sSite4.y = sSite4.y + static_cast<SBYTE>((uiShift >> 1) & 1);
        sSite4.z = sSite4.z + static_cast<SBYTE>((uiShift >> 2) & 1);
        const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldID][__bi(sSite4)];
        if (!sIdx.IsDirichlet())
        {
            pDeviceData[sIdx.m_uiSiteIndex] = _makeColorVector<T>(color);
        }
    }
    
}

#pragma endregion

template<typename T> 
void CCommonKernelSite<T>::InitialBuffer(T* dest, BYTE byFieldId, EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialSite << <block, threads >> > (dest, byFieldId, eInitialType);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

template<typename T>
void CCommonKernelSite<T>::FixBoundary(T* dest, BYTE byFieldId)
{
    preparethread;
    _kernelFixBoundarySite << <block, threads >> > (dest, byFieldId);
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
        _kernelMakePointSourceSite << <block, threads >> > (data, uiSiteIndex, sourceData.m_byColorIndex);
    }
    break;
    case EFS_Wall:
    {
        preparethread;
        _kernelInitialSite << <block, threads >> > (data, byFieldId, EFIT_Zero);
        _kernelMakeWallSourceSite << <block, threads >> > (
            data,
            static_cast<INT>(sourceData.m_sSourcePoint.w),
            sourceData.m_bySpinIndex,
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

template class CCommonKernelSite<CLGComplex>;
template class CCommonKernelSite<deviceSU2Vector>;
template class CCommonKernelSite<deviceSU3Vector>;
template class CCommonKernelSite<deviceSU4Vector>;
template class CCommonKernelSite<deviceSU5Vector>;
template class CCommonKernelSite<deviceSU6Vector>;
template class CCommonKernelSite<deviceSU7Vector>;
template class CCommonKernelSite<deviceSU8Vector>;
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
                    pDevicePtr[uiSiteIndex] = buffer[uiRegion * uiDir + idir];
                }
                pDevicePtr[uiSiteIndex] = id;
                continue;
            }
            pDevicePtr[uiLinkIndex] = _makeRandom<T>(_deviceGetFatIndex(uiSiteIndex, idir + 1));
        }
        break;
        case EFIT_RandomGenerator:
        {
            if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, idir))
            {
                pDevicePtr[uiSiteIndex] = zero;
                continue;
            }
            pDevicePtr[uiLinkIndex] = _makeGaussian<T>(_deviceGetFatIndex(uiSiteIndex, idir + 1));
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
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const UINT uiDir = _DC_Dir;

    const T id = _makeId<T>();
    const T zero = _makeZero<T>();

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, idir))
        {
            if (NULL != __boundaryFieldPointers[byFieldId])
            {
                UINT uiRegion = __idx->_devcieExchangeBoundaryFieldSiteIndexBI(byFieldId, uiBigIdx);
                const T* buffer = ((CFieldBoundary<T>*)__boundaryFieldPointers[byFieldId])->m_pDeviceData;
                pDeviceData[uiSiteIndex] = buffer[uiRegion * uiDir + idir];
            }
            pDeviceData[uiSiteIndex] = bId ? id : zero;
            return;
        }

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
        if (0 != ((1 << idir) & byDir))
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
        if (0 != ((1 << idir) & byDir))
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
    UINT uiXYZ = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    const UINT uiSiteIndex = uiXYZ * _DC_Lt;
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

    res[uiXYZ] = _cToDouble(_tr(tmp));
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


template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateKinematicEnergyT_D(
    BYTE byFieldId,
    const deviceGauge * __restrict__ pDeviceData,
    DOUBLE* results
)
{
    intokernalInt4;
    const UINT uiDir = _DC_Dir;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    DOUBLE resThisThread = 0.0;
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, idir))
        {
            UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
            resThisThread += static_cast<DOUBLE>(_retr(_dagmulC(pDeviceData[linkIndex], pDeviceData[linkIndex])));
        }
    }
    results[uiSiteIndex] = resThisThread;
}

#pragma endregion

template<typename T>
void CCommonKernelLink<T>::InitialBuffer(T* dest, BYTE byFieldId, EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialLink << <block, threads >> > (dest, byFieldId, eInitialType);
}

template<typename T>
void CCommonKernelLink<T>::FixBoundary(T* dest, BYTE byFieldId)
{
    preparethread;
    _kernelFixBoundaryLink << <block, threads >> > (dest, byFieldId, TRUE);
}

template<typename T>
void CCommonKernelLink<T>::FixBoundaryZero(T* dest, BYTE byFieldId)
{
    preparethread;
    _kernelFixBoundaryLink << <block, threads >> > (dest, byFieldId, FALSE);
}

template<typename T>
void CCommonKernelLink<T>::DebugPrint(const T* data, UINT uiLinkCount)
{
    //preparethread;
    //_kernelPrintSU3 << < block, threads >> > (m_pDeviceData);

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
    preparethread;
    _kernelExp << <block, threads >> > (other, me, byFieldId, a);
}

template<typename T>
void CCommonKernelLink<T>::QuickLog(T* data, BYTE byFieldId)
{
    preparethread;
    _kernelQuickLog << <block, threads >> > (data, byFieldId);
}

template<typename T>
void CCommonKernelLink<T>::QuickExp(T* data, BYTE byFieldId)
{
    preparethread;
    _kernelQuickLog << <block, threads >> > (data, byFieldId);
}

template<typename T>
void CCommonKernelLink<T>::StrictLog(T* data, BYTE byFieldId)
{
    preparethread;
    _kernelStrictLog << <block, threads >> > (data, byFieldId);
}

template<typename T>
void CCommonKernelLink<T>::StrictExp(T* data, BYTE byFieldId)
{
    preparethread;
    _kernelStrictExp << <block, threads >> > (data, byFieldId);
}

template<typename T>
void CCommonKernelLink<T>::SetOneDirectionUnity(T* data, BYTE byFieldId, BYTE byDir)
{
    preparethread;
    _kernelOneDirId << <block, threads >> > (data, byFieldId, byDir);
}

template<typename T>
void CCommonKernelLink<T>::SetOneDirectionZero(T* data, BYTE byFieldId, BYTE byDir)
{
    preparethread;
    _kernelOneDirZero << <block, threads >> > (data, byFieldId, byDir);
}

template<typename T>
void CCommonKernelLink<T>::PolyakovOnSpatialSite(const T* data, BYTE byFieldId, cuDoubleComplex* buffer)
{
    dim3 block(_HC_DecompX, _HC_DecompY, 1);
    dim3 threads(_HC_DecompLx, _HC_DecompLy, 1);
    _kernelPolyakovLoopOfSite << <block, threads >> > (data, buffer, byFieldId);
}

template<typename T>
void CCommonKernelLink<T>::CalculateE_Using_U(const T* deviceData, BYTE byFieldId, T* pResoult)
{
    preparethread;
    _kernelTransformToE << <block, threads >> > (byFieldId, deviceData, pResoult);
}

template<typename T>
void CCommonKernelLink<T>::CalculateNablaE_Using_U(const T* deviceData, BYTE byFieldId, T* pResoult, UBOOL bNaive)
{
    preparethread;
    if (bNaive)
    {
        _kernelCalculateNablaENaive << <block, threads >> > (
            deviceData,
            byFieldId,
            pResoult);
    }
    else
    {
        _kernelCalculateNablaE << <block, threads >> > (
            deviceData,
            byFieldId,
            pResoult);
    }
}

template<typename T>
DOUBLE CCommonKernelLink<T>::CalcKineticEnery(const T* me, BYTE byFieldId)
{
    preparethread;
    _kernelCalculateKinematicEnergyT_D << <block, threads >> > (byFieldId, me, _D_RealThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
}

template class CCommonKernelLink<CLGComplex>;
template class CCommonKernelLink<deviceSU2>;
template class CCommonKernelLink<deviceSU3>;
template class CCommonKernelLink<deviceSU4>;
template class CCommonKernelLink<deviceSU5>;
template class CCommonKernelLink<deviceSU6>;
template class CCommonKernelLink<deviceSU7>;
template class CCommonKernelLink<deviceSU8>;

#pragma endregion



__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================