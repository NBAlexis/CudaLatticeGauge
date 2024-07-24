//=============================================================================
// FILENAME : CFieldCommonKernel.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [07/20/2024 nbale]
//=============================================================================
#include "Tools/Math/DeviceTemplates/DeviceInlineUseNoTemplateFunction.h"

#ifndef _CFIELDCOMMON_KERNEL_H_
#define _CFIELDCOMMON_KERNEL_H_

#define gaugeLinkKernelFuncionStart \
    intokernaldir; \
    for (UINT idir = 0; idir < uiDir; ++idir) \
    { \
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir); 


#define gaugeLinkKernelFuncionEnd \
    } 


#define __simplekernel(...) \
const UINT idx = (threadIdx.x + blockIdx.x * blockDim.x); \
if (idx < count) \
{ \
    __VA_ARGS__; \
}

#define __gaugeKernel(...) \
intokernalInt4; \
const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4); \
const UINT uiDir = _DC_Dir; \
for (UINT idir = 0; idir < uiDir; ++idir) \
{ \
    if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, idir)) \
    { \
        const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, idir); \
        __VA_ARGS__;  \
    } \
}


#define __SIMPLEDECOMPOSE(length) \
UINT block = (length) > CCommonData::m_uiMaxThreadPerBlock ? Ceil(length, CCommonData::m_uiMaxThreadPerBlock) : 1; \
UINT thread = (length) > CCommonData::m_uiMaxThreadPerBlock ? Ceil(length, block) : (length);

__BEGIN_NAMESPACE

static inline INT Ceil(LONGLONG a, LONGLONG b)
{
    const LONGLONG c = a / b;
    return static_cast<INT>((a == b * c) ? c : (c + 1));
}

extern void CLGAPI appSimpleMalloc(void** ptr, size_t size);
extern void CLGAPI appSimpleFree(void* ptr);
extern void CLGAPI appSimpleCopyDD(void* target, void* source, size_t size);
extern void CLGAPI appSimpleCopyHD(void* target, void* source, size_t size);


//to avoid explict instantiation of every fucntion, we just explict instantiation a class
template<typename T>
class __DLL_EXPORT CCommonKernel
{
public:
    static void AllocateBuffer(T** pointer, UINT count);
    static void FreeBuffer(T** pointer);
    static void CopyBuffer(T* dest, const T* source, UINT count);
};

template<typename T>
class __DLL_EXPORT CCommonKernelField
{
public:
    static void Initial(T* pointer, UINT count, EFieldInitialType eInitialType);
    static BYTE* CopyDataOut(T* pointer, UINT count, UINT& uiSize);
    static BYTE* CopyDataOutFloat(T* pointer, UINT count, UINT& uiSize);
    static BYTE* CopyDataOutDouble(T* pointer, UINT count, UINT& uiSize);
    static void InitialWithByte(T* pointer, UINT count, const BYTE* data);

    static void Dagger(T* dest, UINT count);
    static void AxpyPlus(T* dest, UINT count, const T* x);
    static void AxpyMinus(T* dest, UINT count, const T* x);
    static void Axpy(T* dest, UINT count, Real a, const T* x);
    static void Axpy(T* dest, UINT count, const CLGComplex& a, const T* x);
    static void Mul(T* dest, UINT count, const T* other, UBOOL bDagger = TRUE);
    static void ScalarMultply(T* dest, UINT count, const CLGComplex& a);
    static void ScalarMultply(T* dest, UINT count, Real a);
    static cuDoubleComplex Dot(const T* me, UINT count, const T* other);
    static DOUBLE LengthSq(const T* me, UINT count);
    static TArray<DOUBLE> Sum(const T* me, UINT count);
    static void Norm(T* dest, UINT count);
};

//Common for boson and fermion
template<typename T>
class __DLL_EXPORT CCommonKernelSite
{
public:
    static void InitialBuffer(T* dest, BYTE byFieldId, EFieldInitialType eInitialType);
    static void FixBoundary(T* dest, BYTE byFieldId);
    static void InitialSource(T* data, BYTE byFieldId, const SFermionBosonSource& sourceData);
    static void DebugPrint(const T* data, UINT sitecount);
};

template<typename T>
class __DLL_EXPORT CCommonKernelLink
{
public:
    static void InitialBuffer(T* dest, BYTE byFieldId, EFieldInitialType eInitialType);
    static void FixBoundary(T* dest, BYTE byFieldId);
    static void FixBoundaryZero(T* dest, BYTE byFieldId);
    static void DebugPrint(const T* data, UINT uiLinkCount);

    static void ExpMul(T* other, BYTE byFieldId, const T* me, Real a);
    static DOUBLE CalcKineticEnery(const T* me, BYTE byFieldId);

    //TA
    static void QuickLog(T* data, BYTE byFieldId);
    //
    static void QuickExp(T* data, BYTE byFieldId);

    static void StrictLog(T* data, BYTE byFieldId);
    static void StrictExp(T* data, BYTE byFieldId);
    static void SetOneDirectionUnity(T* data, BYTE byFieldId, BYTE byDir);
    static void SetOneDirectionZero(T* data, BYTE byFieldId, BYTE byDir);

    static void PolyakovOnSpatialSite(const T* data, BYTE byFieldId, cuDoubleComplex* buffer);

    static void CalculateE_Using_U(const T* deviceData, BYTE byFieldId, T* pResoult);
    static void CalculateNablaE_Using_U(const T* deviceData, BYTE byFieldId, T* pResoult, UBOOL bNaive = FALSE);
};

__END_NAMESPACE

#endif //#ifndef _CFIELDCOMMON_KERNEL_H_

//=============================================================================
// END OF FILE
//=============================================================================