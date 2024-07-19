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

template<typename T>
void CCommonKernel<T>::AllocateBuffer(T*& pointer, UINT count)
{
    checkCudaErrors(__cudaMalloc((void**)&pointer, sizeof(T) * count));
}
template<typename T>
void CCommonKernel<T>::FreeBuffer(T*& pointer)
{
    checkCudaErrors(__cudaFree(pointer));
    pointer = NULL;
}

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

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================