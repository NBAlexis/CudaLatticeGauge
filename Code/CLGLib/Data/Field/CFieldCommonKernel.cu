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
_kernelAddSite(T* pMe, const T* __restrict__ pOther)
{
    intokernal;
    _add(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelSubSite(T* pMe, const T* __restrict__ pOther)
{
    intokernal;
    _sub(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelMulSite(T* pMe, const T* __restrict__ pOther)
{
    intokernal;
    _mul(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpySite(T* pMe, const T* __restrict__ pOther, CLGComplex a)
{
    intokernal;
    _add(pMe[uiSiteIndex], _mulC(pOther[uiSiteIndex], a));
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelMulComplexSite(T* pMe, const T* __restrict__ pOther, UBOOL bConj)
{
    intokernal;
    if (bConj)
    {
        _dagger(pMe[uiSiteIndex]);
    }
    _mul(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyRealSite(T* pMe, const T* __restrict__ pOther, Real a)
{
    intokernal;
    _add(pMe[uiSiteIndex], _mulC(pOther[uiSiteIndex], a));
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelDotSite(const T* __restrict__ pMe, const T* __restrict__ pOther, cuDoubleComplex* result)
{
    intokernal;
    result[uiSiteIndex] = _cToDouble(_dot(pMe[uiSiteIndex], pOther[uiSiteIndex]));
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelElementSite(const T* __restrict__ pMe, UINT idx, DOUBLE* result)
{
    intokernal;
    result[uiSiteIndex] = static_cast<DOUBLE>(_element(pMe[uiSiteIndex], idx));
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplySite(T* pMe, CLGComplex a)
{
    intokernal;
    _mul(pMe[uiSiteIndex], a);
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyRealSite(T* pMe, Real a)
{
    intokernal;
    _mul(pMe[uiSiteIndex], a);
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelConjugateSite(T* pDeviceData)
{
    intokernal;
    _dagger(pDeviceData[uiSiteIndex]);
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
void CCommonKernelSite<T>::Dagger(T* dest, BYTE byFieldId)
{
    preparethread;
    _kernelConjugateSite << <block, threads >> > (dest);
}

template<typename T>
void CCommonKernelSite<T>::FixBoundary(T* dest, BYTE byFieldId)
{
    preparethread;
    _kernelFixBoundarySite << <block, threads >> > (dest, byFieldId);
}

template<typename T>
void CCommonKernelSite<T>::AxpyPlus(T* dest, BYTE byFieldId, const T* x)
{
    preparethread;
    _kernelAddSite << <block, threads >> > (dest, x);
}

template<typename T>
void CCommonKernelSite<T>::AxpyMinus(T* dest, BYTE byFieldId, const T* x)
{
    preparethread;
    _kernelSubSite << <block, threads >> > (dest, x);
}

template<typename T>
void CCommonKernelSite<T>::Axpy(T* dest, BYTE byFieldId, Real a, const T* x)
{
    preparethread;
    _kernelAxpyRealSite << <block, threads >> > (dest, x, a);
}

template<typename T>
void CCommonKernelSite<T>::Axpy(T* dest, BYTE byFieldId, const CLGComplex& a, const T* x)
{
    preparethread;
    _kernelAxpySite << <block, threads >> > (dest, x, a);
}

template<typename T>
void CCommonKernelSite<T>::Mul(T* dest, BYTE byFieldId, const T* x, UBOOL bDagger)
{
    preparethread;
    _kernelMulComplexSite << <block, threads >> > (dest, x, bDagger);
}

template<typename T>
void CCommonKernelSite<T>::ScalarMultply(T* dest, BYTE byFieldId, const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplySite << <block, threads >> > (dest, a);
}

template<typename T>
void CCommonKernelSite<T>::ScalarMultply(T* dest, BYTE byFieldId, Real a)
{
    preparethread;
    _kernelScalarMultiplyRealSite << <block, threads >> > (dest, a);
}

template<typename T>
cuDoubleComplex CCommonKernelSite<T>::Dot(const T* me, BYTE byFieldId, const T* other)
{
    preparethread;
    _kernelDotSite << <block, threads >> > (me, other, _D_ComplexThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

template<typename T>
TArray<DOUBLE> CCommonKernelSite<T>::Sum(const T* me, BYTE byFieldId)
{
    preparethread;
    TArray<DOUBLE> ret;
    for (UINT i = 0; i < _elementdim<T>(); ++i)
    {
        _kernelElementSite << <block, threads >> > (me, i, _D_RealThreadBuffer);
        ret.AddItem(appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer));
    }
    return ret;
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
    T id = _makeId<T>();
    T zero = _makeZero<T>();

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
                pDevicePtr[uiSiteIndex] = _makeZero<T>();
                continue;
            }
            pDevicePtr[uiLinkIndex] = _makeRandom<T>(_deviceGetFatIndex(uiSiteIndex, idir + 1));
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
                    pDevicePtr[uiSiteIndex] = buffer[uiRegion * uiDir + idir];
                }
                pDevicePtr[uiSiteIndex] = _makeId<T>();
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
_kernelDaggerLink(T* pDevicePtr)
{
    gaugeLinkKernelFuncionStart

        _dagger(pDevicePtr[uiLinkIndex]);

    gaugeLinkKernelFuncionEnd
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyLink(T* pDevicePtr, const T* __restrict__ x, CLGComplex a)
{
    gaugeLinkKernelFuncionStart

        _add(pDevicePtr[uiLinkIndex], _mulC(x[uiLinkIndex], a));

    gaugeLinkKernelFuncionEnd
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelMulLink(T* pDevicePtr, const T* __restrict__ x, UBOOL bDagger)
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

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyLinkReal(T* pDevicePtr, const T* __restrict__ x, Real a)
{
    gaugeLinkKernelFuncionStart

        _add(pDevicePtr[uiLinkIndex], _mulC(x[uiLinkIndex], a));

    gaugeLinkKernelFuncionEnd
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusLink(T* pDevicePtr, const T* __restrict__ x)
{
    gaugeLinkKernelFuncionStart

        _add(pDevicePtr[uiLinkIndex], x[uiLinkIndex]);

    gaugeLinkKernelFuncionEnd
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusLink(T* pDevicePtr, const T* __restrict__ x)
{
    gaugeLinkKernelFuncionStart

        _sub(pDevicePtr[uiLinkIndex], x[uiLinkIndex]);

    gaugeLinkKernelFuncionEnd
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyLinkComplex(T* pDevicePtr, CLGComplex a)
{
    gaugeLinkKernelFuncionStart

        _mul(pDevicePtr[uiLinkIndex], a);

    gaugeLinkKernelFuncionEnd
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyLinkReal(T* pDevicePtr, Real a)
{
    gaugeLinkKernelFuncionStart

        _mul(pDevicePtr[uiLinkIndex], a);

    gaugeLinkKernelFuncionEnd
}

template<typename T>
__global__ void _CLG_LAUNCH_BOUND
_kernelFixBoundaryLink(T* pDeviceData, BYTE byFieldId)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const UINT uiDir = _DC_Dir;

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
            pDeviceData[uiSiteIndex] = _makeZero<T>();
            return;
        }

    }
}

template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDotLink(
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

#pragma endregion

template<typename T>
void CCommonKernelLink<T>::InitialBuffer(T* dest, BYTE byFieldId, EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialLink << <block, threads >> > (dest, byFieldId, eInitialType);
}

template<typename T>
void CCommonKernelLink<T>::Dagger(T* dest, BYTE byFieldId)
{
    preparethread;
    _kernelDaggerLink << <block, threads >> > (dest);
}

template<typename T>
void CCommonKernelLink<T>::FixBoundary(T* dest, BYTE byFieldId)
{
    preparethread;
    _kernelFixBoundaryLink << <block, threads >> > (dest, byFieldId);
}

template<typename T>
void CCommonKernelLink<T>::AxpyPlus(T* dest, BYTE byFieldId, const T* x)
{
    preparethread;
    _kernelAxpyPlusLink << <block, threads >> > (dest, x);
}

template<typename T>
void CCommonKernelLink<T>::AxpyMinus(T* dest, BYTE byFieldId, const T* x)
{
    preparethread;
    _kernelAxpyMinusLink << <block, threads >> > (dest, x);
}

template<typename T>
void CCommonKernelLink<T>::Axpy(T* dest, BYTE byFieldId, Real a, const T* x)
{
    preparethread;
    _kernelAxpyLinkReal << <block, threads >> > (dest, x, a);
}

template<typename T>
void CCommonKernelLink<T>::Axpy(T* dest, BYTE byFieldId, const CLGComplex& a, const T* x)
{
    preparethread;
    _kernelAxpyLink << <block, threads >> > (dest, x, a);
}

template<typename T>
void CCommonKernelLink<T>::Mul(T* dest, BYTE byFieldId, const T* x, UBOOL bDagger)
{
    preparethread;
    _kernelMulLink << <block, threads >> > (dest, x, bDagger);
}

template<typename T>
void CCommonKernelLink<T>::ScalarMultply(T* dest, BYTE byFieldId, const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplyLinkComplex << <block, threads >> > (dest, a);
}

template<typename T>
void CCommonKernelLink<T>::ScalarMultply(T* dest, BYTE byFieldId, Real a)
{
    preparethread;
    _kernelScalarMultiplyLinkReal << <block, threads >> > (dest, a);
}

template<typename T>
cuDoubleComplex CCommonKernelLink<T>::Dot(const T* me, BYTE byFieldId, const T* other)
{
    preparethread;
    _kernelDotLink << < block, threads >> > (me, other, _D_ComplexThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
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