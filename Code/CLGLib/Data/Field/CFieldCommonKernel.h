//=============================================================================
// FILENAME : CFieldCommonKernel.h
// 
// DESCRIPTION:
// This is the class for all boson fields
//
// REVISION:
//  [mm/dd/yy]
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
intokernalDirInt4; \
const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4); \
if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir)) \
{ \
    __VA_ARGS__;  \
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

extern void CLGAPI appSimpleCopyDD(void* target, const void* source, size_t size);
extern void CLGAPI appSimpleCopyHD(void* target, const void* source, size_t size);
extern void CLGAPI appSimpleCopyDH(void* target, const void* source, size_t size);

//=============================================================================
// Global (cross-rank) reductions for the multi-GPU build.
//
// The device-tree reduction (ThreadBufferSum) sums only over THIS rank's local
// sub-lattice (_HC_Volume; halo slots are appended after the volume and are
// never touched by the reduction kernels). A field-wide Dot / LengthSq / Sum /
// kinetic energy therefore needs one MPI_Allreduce on the resulting host scalar
// so every rank sees the same global sum. We deliberately wrap ABOVE
// ThreadBufferSum -- injecting inside it would also globalise measurement and
// gauge-fixing paths that reduce local-only quantities (some divide by the
// local _HC_Volume). AllreduceSum is a no-op when _CLG_MULTI_GPU is 0 or the
// world size is 1, so single-GPU behaviour is unchanged. Reduce in DOUBLE /
// cuDoubleComplex here, before any narrowing to Real, to keep the global sum
// full precision. See Docs/MultiGPU-Plan.md Phase 3 and section 2.6.
//
// I9: promoted from CFieldCommonKernel.cu (was file-static) so the field and
// action layers share this single implementation instead of mirroring it.
//=============================================================================
inline DOUBLE _clgGlobalThreadBufferSum(DOUBLE* pDeviceBuffer)
{
    DOUBLE dResult = appGetCudaHelper()->ThreadBufferSum(pDeviceBuffer);
#if _CLG_MULTI_GPU
    if (NULL != appGetComm())
    {
        appGetComm()->AllreduceSum(dResult);
    }
#endif
    return dResult;
}

inline cuDoubleComplex _clgGlobalThreadBufferSum(cuDoubleComplex* pDeviceBuffer)
{
    cuDoubleComplex cResult = appGetCudaHelper()->ThreadBufferSum(pDeviceBuffer);
#if _CLG_MULTI_GPU
    if (NULL != appGetComm())
    {
        appGetComm()->AllreduceSum(cResult);
    }
#endif
    return cResult;
}

#if _CLG_MULTI_GPU
//Multi-GPU (Phase 1): gather this rank's boundary sites (halo-slot order) from a
//field's device buffer into a contiguous device buffer. bytesPerSite must be a
//multiple of 4. Returns the halo site count (0 when nothing is split).
extern UINT CLGAPI appHaloGatherLocal(const BYTE* pFieldDevice, BYTE* pDstContiguousDevice, UINT uiBytesPerSite);
#endif


//to avoid explict instantiation of every fucntion, we just explict instantiation a class
template<typename T>
#if _CLG_WIN
class __DLL_EXPORT CCommonKernel
#else
class CCommonKernel
#endif
{
public:
    static void CopyBuffer(T* dest, const T* source, UINT count);
};

#if !_CLG_WIN
extern template class CCommonKernel<INT>;
#if _CLG_DOUBLEFLOAT
extern template class CCommonKernel<FLOAT>;
#else
extern template class CCommonKernel<DOUBLE>;
#endif
extern template class CCommonKernel<Real>;
#if _CLG_DOUBLEFLOAT
extern template class CCommonKernel<cuComplex>;
#else
extern template class CCommonKernel<cuDoubleComplex>;
#endif
extern template class CCommonKernel<CLGComplex>;
extern template class CCommonKernel<deviceSU2>;
extern template class CCommonKernel<deviceSU3>;
extern template class CCommonKernel<deviceSU2Vector>;
extern template class CCommonKernel<deviceSU3Vector>;
extern template class CCommonKernel<deviceWilsonVectorSU3>;

#if _CLG_SU4_GAUGE
extern template class CCommonKernel<deviceSU4>;
#endif
#if _CLG_SU5_GAUGE
extern template class CCommonKernel<deviceSU5>;
#endif
#if _CLG_SU6_GAUGE
extern template class CCommonKernel<deviceSU6>;
#endif
#if _CLG_SU7_GAUGE
extern template class CCommonKernel<deviceSU7>;
#endif
#if _CLG_SU8_GAUGE
extern template class CCommonKernel<deviceSU8>;
#endif

#if _CLG_Z2_GAUGE
extern template class CCommonKernel<deviceZN<2>>;
#endif
#if _CLG_Z3_GAUGE
extern template class CCommonKernel<deviceZN<3>>;
#endif
#if _CLG_Z4_GAUGE
extern template class CCommonKernel<deviceZN<4>>;
#endif
#if _CLG_Z5_GAUGE
extern template class CCommonKernel<deviceZN<5>>;
#endif
#if _CLG_Z6_GAUGE
extern template class CCommonKernel<deviceZN<6>>;
#endif

#if _CLG_D3_GAUGE
extern template class CCommonKernel<deviceDN<3>>;
#endif
#if _CLG_D4_GAUGE
extern template class CCommonKernel<deviceDN<4>>;
#endif
#if _CLG_D8_GAUGE
extern template class CCommonKernel<deviceDN<8>>;
#endif

#if _CLG_SU4_KS || _CLG_SU4_BOSON
extern template class CCommonKernel<deviceSU4Vector>;
#endif
#if _CLG_SU5_KS || _CLG_SU5_BOSON
extern template class CCommonKernel<deviceSU5Vector>;
#endif
#if _CLG_SU6_KS || _CLG_SU6_BOSON
extern template class CCommonKernel<deviceSU6Vector>;
#endif
#if _CLG_SU7_KS || _CLG_SU7_BOSON
extern template class CCommonKernel<deviceSU7Vector>;
#endif
#if _CLG_SU8_KS || _CLG_SU8_BOSON
extern template class CCommonKernel<deviceSU8Vector>;
#endif

#if _CLG_SL3C_GAUGE
extern template class CCommonKernel<deviceSL3C>;
#endif
#if _CLG_U3_GAUGE
extern template class CCommonKernel<deviceU3>;
#endif
#if _CLG_O3_GAUGE
extern template class CCommonKernel<deviceO3>;
#endif
#if _CLG_SO3_GAUGE
extern template class CCommonKernel<deviceSO3>;
#endif
#endif

template<typename T>
#if _CLG_WIN
class __DLL_EXPORT CCommonKernelField
#else
class CCommonKernelField
#endif
{
public:
    static void Initial(T* pointer, UINT count, EFieldInitialType eInitialType);
    static void InitialEvenOdd(T* pointer, UINT count, EFieldInitialType eInitialType, UBOOL bSite, UBOOL bEven, UBOOL bZeroOtherSites);
    static void ZeroEvenOddSite(T* pointer, UINT count, UBOOL bEven);
    static BYTE* CopyDataOut(T* pointer, UINT count, UINT& uiSize);
    static BYTE* CopyDataOutFloat(T* pointer, UINT count, UINT& uiSize);
    static BYTE* CopyDataOutDouble(T* pointer, UINT count, UINT& uiSize);
    static void InitialWithByte(T* pointer, UINT count, const BYTE* data);

    static void Dagger(T* dest, UINT count);
    static void AxpyPlus(T* dest, UINT count, const T* x);
    static void AxpyMinus(T* dest, UINT count, const T* x);
    static void Axpy(T* dest, UINT count, Real a, const T* x);
    static void Axpy(T* dest, UINT count, const CLGComplex& a, const T* x);
    //Axpy on even or odd sites only, used by even pseudo-fermion, whose other sites are always zero
    static void AxpyPlusEvenOdd(T* dest, UINT count, const T* x, UBOOL bEven);
    static void AxpyMinusEvenOdd(T* dest, UINT count, const T* x, UBOOL bEven);
    static void AxpyEvenOdd(T* dest, UINT count, Real a, const T* x, UBOOL bEven);
    static void AxpyEvenOdd(T* dest, UINT count, const CLGComplex& a, const T* x, UBOOL bEven);
    static void Mul(T* dest, UINT count, const T* other, UBOOL bDaggerLeft = TRUE, UBOOL bDaggerRight = FALSE);
    static void LeftMul(T* dest, UINT count, const T* other, UBOOL bDaggerLeft = TRUE, UBOOL bDaggerRight = FALSE);
    static void ScalarMultply(T* dest, UINT count, const CLGComplex& a);
    static void ScalarMultply(T* dest, UINT count, Real a);
    static void ApplyPhaseC(T* dest, UINT count, const CLGComplex* other);
    static void ApplyPhaseR(T* dest, UINT count, const Real* other, Real fCharge);
    static cuDoubleComplex Dot(const T* me, UINT count, const T* other);
    static DOUBLE LengthSq(const T* me, UINT count);
    static TArray<DOUBLE> Sum(const T* me, UINT count);
    static void Norm(T* dest, UINT count);
};

#if !_CLG_WIN
extern template class CCommonKernelField<Real>;
extern template class CCommonKernelField<CLGComplex>;
extern template class CCommonKernelField<deviceSU2>;
extern template class CCommonKernelField<deviceSU3>;
extern template class CCommonKernelField<deviceSU2Vector>;
extern template class CCommonKernelField<deviceSU3Vector>;
extern template class CCommonKernelField<deviceWilsonVectorSU3>;

#if _CLG_SU4_GAUGE
extern template class CCommonKernelField<deviceSU4>;
#endif
#if _CLG_SU5_GAUGE
extern template class CCommonKernelField<deviceSU5>;
#endif
#if _CLG_SU6_GAUGE
extern template class CCommonKernelField<deviceSU6>;
#endif
#if _CLG_SU7_GAUGE
extern template class CCommonKernelField<deviceSU7>;
#endif
#if _CLG_SU8_GAUGE
extern template class CCommonKernelField<deviceSU8>;
#endif

#if _CLG_Z2_GAUGE
extern template class CCommonKernelField<deviceZN<2>>;
#endif
#if _CLG_Z3_GAUGE
extern template class CCommonKernelField<deviceZN<3>>;
#endif
#if _CLG_Z4_GAUGE
extern template class CCommonKernelField<deviceZN<4>>;
#endif
#if _CLG_Z5_GAUGE
extern template class CCommonKernelField<deviceZN<5>>;
#endif
#if _CLG_Z6_GAUGE
extern template class CCommonKernelField<deviceZN<6>>;
#endif

#if _CLG_D3_GAUGE
extern template class CCommonKernelField<deviceDN<3>>;
#endif
#if _CLG_D4_GAUGE
extern template class CCommonKernelField<deviceDN<4>>;
#endif
#if _CLG_D8_GAUGE
extern template class CCommonKernelField<deviceDN<8>>;
#endif

#if _CLG_SL3C_GAUGE
extern template class CCommonKernelField<deviceSL3C>;
#endif
#if _CLG_U3_GAUGE
extern template class CCommonKernelField<deviceU3>;
#endif
#if _CLG_O3_GAUGE
extern template class CCommonKernelField<deviceO3>;
#endif
#if _CLG_SO3_GAUGE
extern template class CCommonKernelField<deviceSO3>;
#endif
#endif

//Common for boson and fermion
template<typename T>
#if _CLG_WIN
class __DLL_EXPORT CCommonKernelSite
#else
class CCommonKernelSite
#endif
{
public:
    static void InitialBuffer(T* dest, BYTE byFieldId, EFieldInitialType eInitialType);
    static void InitialBufferEvenOdd(T* dest, BYTE byFieldId, UBOOL bEven, UBOOL bZeroOtherSites, EFieldInitialType eInitialType);
    static void FixBoundary(T* dest, BYTE byFieldId);
    static void InitialSource(T* data, BYTE byFieldId, const SFermionBosonSource& sourceData);
    static void DebugPrint(const T* data, UINT sitecount);
    static void DiagnalTerm(
        T* pTarget,
        BYTE byFieldId,
        const T* pSource, DOUBLE fCoeffiecient, _deviceCoeffFunctionPointer fpCoeff,
        EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff);
};

#if !_CLG_WIN
extern template class CCommonKernelSite<Real>;
extern template class CCommonKernelSite<CLGComplex>;
extern template class CCommonKernelSite<deviceSU2Vector>;
extern template class CCommonKernelSite<deviceSU3Vector>;

#if _CLG_SU4_KS || _CLG_SU4_BOSON
extern template class CCommonKernelSite<deviceSU4Vector>;
#endif
#if _CLG_SU5_KS || _CLG_SU5_BOSON
extern template class CCommonKernelSite<deviceSU5Vector>;
#endif
#if _CLG_SU6_KS || _CLG_SU6_BOSON
extern template class CCommonKernelSite<deviceSU6Vector>;
#endif
#if _CLG_SU7_KS || _CLG_SU7_BOSON
extern template class CCommonKernelSite<deviceSU7Vector>;
#endif
#if _CLG_SU8_KS || _CLG_SU8_BOSON
extern template class CCommonKernelSite<deviceSU8Vector>;
#endif
#endif

template<typename T>
#if _CLG_WIN
class __DLL_EXPORT CCommonKernelLink
#else
class CCommonKernelLink
#endif
{
public:
    static void InitialBuffer(T* dest, BYTE byFieldId, EFieldInitialType eInitialType);
    static void InitialBufferD(T* dest, BYTE byFieldId, EFieldInitialType eInitialType);
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

    static void PolyakovOnSpatialSite(const T* data, BYTE byFieldId, cuDoubleComplex* buffer, BYTE byDir = 3);

    static void CalculateE_Using_U(const T* deviceData, BYTE byFieldId, T* pResoult);
    static void CalculateNablaE_Using_U(const T* deviceData, BYTE byFieldId, T* pResoult, UBOOL bNaive = FALSE);

    static void AddLink(const T* source, T* target, Real fCoeff, const SCHAR* devicePath, BYTE byPathLen, BYTE byMu, BYTE byFieldId);

    /**
    * add for test
    */
    static void ApplyStaggeredPhase(T* deviceData, BYTE byFieldId);
    static void ApplyStaggeredPhaseMILC(T* deviceData, BYTE byFieldId);
};

#if !_CLG_WIN
extern template class CCommonKernelLink<CLGComplex>;
extern template class CCommonKernelLink<deviceSU2>;
extern template class CCommonKernelLink<deviceSU3>;

#if _CLG_SU4_GAUGE
extern template class CCommonKernelLink<deviceSU4>;
#endif
#if _CLG_SU5_GAUGE
extern template class CCommonKernelLink<deviceSU5>;
#endif
#if _CLG_SU6_GAUGE
extern template class CCommonKernelLink<deviceSU6>;
#endif
#if _CLG_SU7_GAUGE
extern template class CCommonKernelLink<deviceSU7>;
#endif
#if _CLG_SU8_GAUGE
extern template class CCommonKernelLink<deviceSU8>;
#endif

#if _CLG_Z2_GAUGE
extern template class CCommonKernelLink<deviceZN<2>>;
#endif
#if _CLG_Z3_GAUGE
extern template class CCommonKernelLink<deviceZN<3>>;
#endif
#if _CLG_Z4_GAUGE
extern template class CCommonKernelLink<deviceZN<4>>;
#endif
#if _CLG_Z5_GAUGE
extern template class CCommonKernelLink<deviceZN<5>>;
#endif
#if _CLG_Z6_GAUGE
extern template class CCommonKernelLink<deviceZN<6>>;
#endif

#if _CLG_D3_GAUGE
extern template class CCommonKernelLink<deviceDN<3>>;
#endif
#if _CLG_D4_GAUGE
extern template class CCommonKernelLink<deviceDN<4>>;
#endif
#if _CLG_D8_GAUGE
extern template class CCommonKernelLink<deviceDN<8>>;
#endif

#if _CLG_SL3C_GAUGE
extern template class CCommonKernelLink<deviceSL3C>;
#endif
#if _CLG_U3_GAUGE
extern template class CCommonKernelLink<deviceU3>;
#endif
#if _CLG_O3_GAUGE
extern template class CCommonKernelLink<deviceO3>;
#endif
#if _CLG_SO3_GAUGE
extern template class CCommonKernelLink<deviceSO3>;
#endif
#endif

/**
* simple calculations involving matrix and vector (common functions for fermion and boson with gauge)
*/
template<typename vector, typename matrix, INT vectorN>
#if _CLG_WIN
class __DLL_EXPORT CCommonKernelMV
#else
class CCommonKernelMV
#endif
{
public:

    /**
    * v[even] = X
    * v[odd] = Y
    * when n is even, res[n_mu] = Y[n+mu]X[n]^+ = v[n+mu]v[n]^+
    * when n is odd,  res[n_mu] = X[n+mu]Y[n]^+ = v[n+mu]v[n]^+
    */
    static void ConnectionOneField(const vector* v, matrix* res, BYTE byFieldId);

    /**
    * similar as ConnectionOneField, but for calculation of force, it is reversed for the convinience of calculation of force:
    * when n is even, res[n_mu] = X[n]Y[n+mu]^+ = v[n]v[n+mu]^+
    * when n is odd,  res[n_mu] = Y[n]X[n+mu]^+ = v[n]v[n+mu]^+
    */
    static void ConnectionOneFieldStaggered(const vector* v, matrix* res, BYTE byFieldId, Real fCoeff);
    static void AddConnectionOneFieldStaggered(const vector* v, matrix* res, BYTE byFieldId, Real fCoeff);

    //static void ConnectionOneFieldStaggeredTest(const vector* v, matrix* res, BYTE byFieldId, Real fCoeff);
    //static void AddConnectionOneFieldStaggeredTest(const vector* v, matrix* res, BYTE byFieldId, Real fCoeff);

    /**
    * calculate g_mu(n) = n_p_m(x+mu). n(x)^+
    */
    static void Connection(const vector* n, const vector* n_p_m, matrix* res, BYTE byFieldId);

    /**
    * calculate g_mu(n) = g_mu(n) + coef * (n_p_m(x+mu). n(x)^+)
    */
    static void AddConnection(const vector* n, const vector* n_p_m, matrix* res, Real fCoeef, BYTE byFieldId);

};

#if !_CLG_WIN
//Deepseek said put them also in .h will speed up linker
extern template class CCommonKernelMV<CLGComplex, CLGComplex, 1>;
extern template class CCommonKernelMV<deviceSU2Vector, deviceSU2, 2>;
extern template class CCommonKernelMV<deviceSU3Vector, deviceSU3, 3>;

#if _CLG_SU4_KS || _CLG_SU4_BOSON
extern template class CCommonKernelMV<deviceSU4Vector, deviceSU4, 4>;
#endif
#if _CLG_SU5_KS || _CLG_SU5_BOSON
extern template class CCommonKernelMV<deviceSU5Vector, deviceSU5, 5>;
#endif
#if _CLG_SU6_KS || _CLG_SU6_BOSON
extern template class CCommonKernelMV<deviceSU6Vector, deviceSU6, 6>;
#endif
#if _CLG_SU7_KS || _CLG_SU7_BOSON
extern template class CCommonKernelMV<deviceSU7Vector, deviceSU7, 7>;
#endif
#if _CLG_SU8_KS || _CLG_SU8_BOSON
extern template class CCommonKernelMV<deviceSU8Vector, deviceSU8, 8>;
#endif
#endif

__END_NAMESPACE

#endif //#ifndef _CFIELDCOMMON_KERNEL_H_

//=============================================================================
// END OF FILE
//=============================================================================