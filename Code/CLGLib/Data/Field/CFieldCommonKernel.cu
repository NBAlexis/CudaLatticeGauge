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

#pragma region common

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

template<typename T>
void CCommonKernel<T>::Initial(T* pointer, UINT count, EFieldInitialType eInitialType)
{
    __SIMPLEDECOMPOSE(count);
    _kernelInitialCommon << <block, thread >> > (pointer, count, eInitialType);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
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
        printf("Wilson Fermion Field cannot be initialized with this type! %d\n", eInitialType);
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
cuDoubleComplex CCommonKernelSite<T>::Dot(T* dest, BYTE byFieldId, const T* other)
{
    preparethread;
    _kernelDotSite << <block, threads >> > (dest, other, _D_ComplexThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

template<typename T>
TArray<DOUBLE> CCommonKernelSite<T>::Sum(T* dest, BYTE byFieldId)
{
    preparethread;
    TArray<DOUBLE> ret;
    for (UINT i = 0; i < _elementdim<T>(); ++i)
    {
        _kernelElementSite << <block, threads >> > (dest, i, _D_RealThreadBuffer);
        ret.AddItem(appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer));
    }
    return ret;
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
cuDoubleComplex CCommonKernelLink<T>::Dot(T* dest, BYTE byFieldId, const T* other)
{
    preparethread;
    _kernelDotLink << < block, threads >> > (dest, other, _D_ComplexThreadBuffer);
    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
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