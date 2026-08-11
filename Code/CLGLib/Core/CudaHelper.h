//=============================================================================
// FILENAME : CudaHelper.h
// 
// DESCRIPTION:
// This is the file for some common CUDA usage
//
// REVISION:
//  [mm/dd/yy]
//  [12/3/2018 nbale]
//=============================================================================
#pragma once

#ifndef _CUDAHELPER_H_
#define _CUDAHELPER_H_

__BEGIN_NAMESPACE

inline class CCudaHelper* appGetCudaHelper();

static inline UINT appCeil(UINT a, UINT b)
{
    return (a + b - 1) / b;
}

static inline void appBlockThreads(UINT threads, UINT& blockvar, UINT& threadvar);
static inline void appBlockThreadsE(UINT threads, UINT elements, UINT& blockvar, UINT& threadvar);

extern __device__ __constant__ UINT _constIntegers[kContentLength];
extern __device__ __constant__ INT _constSignedIntegers[kContentLength];
extern __device__ __constant__ Real _constFloats[kContentLength];
__device__ __constant__ constexpr SCHAR _plaq_idx[6][2] = {
    {1, 2},
    {1, 3},
    {1, 4},
    {2, 3},
    {2, 4},
    {3, 4},
};

/**
* Note that, the pointers are copyied here. So the virtual functions should not be used!
*/
extern __device__ __constant__ class CField* __fieldPointers[kMaxFieldCount];
extern __device__ __constant__ class CFieldBoundaryParent* __boundaryFieldPointers[kMaxFieldCount];

extern __device__ __constant__ class CRandom* __r;
extern __device__ __constant__ class CIndexData* __idx;

//NOTE!!!!
//SIGMA41, SIGMA42, SIGMA43 ARE ACCUATELY SIGMA14, SIGMA24, SIGMA34!!
__DEFINE_ENUM(EGammaMatrix,
    UNITY,
    GAMMA1,
    GAMMA2,
    GAMMA3,
    GAMMA4,
    GAMMA5,
    GAMMA51,
    GAMMA52,
    GAMMA53,
    GAMMA54,
    GAMMA15,
    GAMMA25,
    GAMMA35,
    GAMMA45,
    SIGMA12,
    SIGMA23,
    SIGMA31,
    SIGMA41,
    SIGMA42,
    SIGMA43,
    SIGMA12E,
    SIGMA23E,
    SIGMA31E,
    CHARGECONJG,
    EGM_MAX,
    )

//extern __constant__ struct gammaMatrix __diracGamma[EGM_MAX];
extern __device__ __constant__ struct gammaMatrix __chiralGamma[EGM_MAX];

extern __device__ __constant__ struct deviceSU3 __SU3Generators[9];

#define _ARG_PACK_IMPL(N, ...) EXPAND(_ARG_PACK_##N(__VA_ARGS__))

#define _ARG_PACK_1(a) {_vs(a)}
#define _ARG_PACK_2(a, b) {_vs(a), _vs(b)}
#define _ARG_PACK_3(a, b, c) {_vs(a), _vs(b), _vs(c)}
#define _ARG_PACK_4(a1, a2, a3, a4) {_vs(a1), _vs(a2), _vs(a3), _vs(a4)}
#define _ARG_PACK_5(a1, a2, a3, a4, a5) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5)}
#define _ARG_PACK_6(a1, a2, a3, a4, a5, a6) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6)}
#define _ARG_PACK_7(a1, a2, a3, a4, a5, a6, a7) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7)}
#define _ARG_PACK_8(a1, a2, a3, a4, a5, a6, a7, a8) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8)}
#define _ARG_PACK_9(a1, a2, a3, a4, a5, a6, a7, a8, a9) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9)}
#define _ARG_PACK_10(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10)}
#define _ARG_PACK_11(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11)}
#define _ARG_PACK_12(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12)}
#define _ARG_PACK_13(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12), _vs(a13)}
#define _ARG_PACK_14(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12), _vs(a13), _vs(a14)}
#define _ARG_PACK_15(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12), _vs(a13), _vs(a14), _vs(a15)}
#define _ARG_PACK_16(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12), _vs(a13), _vs(a14), _vs(a15), _vs(a16)}
#define _ARG_PACK_17(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12), _vs(a13), _vs(a14), _vs(a15), _vs(a16), _vs(a17)}
#define _ARG_PACK_18(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12), _vs(a13), _vs(a14), _vs(a15), _vs(a16), _vs(a17), _vs(a18)}
#define _ARG_PACK_19(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12), _vs(a13), _vs(a14), _vs(a15), _vs(a16), _vs(a17), _vs(a18), _vs(a19)}
#define _ARG_PACK_20(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12), _vs(a13), _vs(a14), _vs(a15), _vs(a16), _vs(a17), _vs(a18), _vs(a19), _vs(a20)}
#define _ARG_PACK_21(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12), _vs(a13), _vs(a14), _vs(a15), _vs(a16), _vs(a17), _vs(a18), _vs(a19), _vs(a20), _vs(a21)}
#define _ARG_PACK_22(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12), _vs(a13), _vs(a14), _vs(a15), _vs(a16), _vs(a17), _vs(a18), _vs(a19), _vs(a20), _vs(a21), _vs(a22)}
#define _ARG_PACK_23(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12), _vs(a13), _vs(a14), _vs(a15), _vs(a16), _vs(a17), _vs(a18), _vs(a19), _vs(a20), _vs(a21), _vs(a22), _vs(a23)}
#define _ARG_PACK_24(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12), _vs(a13), _vs(a14), _vs(a15), _vs(a16), _vs(a17), _vs(a18), _vs(a19), _vs(a20), _vs(a21), _vs(a22), _vs(a23), _vs(a24)}
#define _ARG_PACK_25(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25) {_vs(a1), _vs(a2), _vs(a3), _vs(a4), _vs(a5), _vs(a6), _vs(a7), _vs(a8), _vs(a9), _vs(a10), _vs(a11), _vs(a12), _vs(a13), _vs(a14), _vs(a15), _vs(a16), _vs(a17), _vs(a18), _vs(a19), _vs(a20), _vs(a21), _vs(a22), _vs(a23), _vs(a24), _vs(a25)}


#define _LAUN_KERNEL_IMPL(func, block, thread, count, ...) \
void* _arg_##func[] = _ARG_PACK_IMPL(count, __VA_ARGS__); \
launchKernel(func, block, thread, _arg_##func); 

#define _LAUN_KERNEL_SHAREMEM_IMPL(func, block, thread, sharemem, count, ...) \
void* _arg_##func[] = _ARG_PACK_IMPL(count, __VA_ARGS__); \
launchKernel(func, block, thread, sharemem, _arg_##func);

#define _LAUN_KERNEL(func, block, thread, ...) \
{ \
    const void* funcvar = (const void*) func; \
    _LAUN_KERNEL_IMPL(funcvar, block, thread, GET_ARG_COUNT(__VA_ARGS__), __VA_ARGS__) \
    appGetCudaHelper()->FreeStack(); \
}

#define _LAUN_KERNEL_SHAREMEM(func, block, thread, sharemem, ...) \
{ \
    const void* funcvar = (const void*) func; \
    _LAUN_KERNEL_SHAREMEM_IMPL(funcvar, block, thread, sharemem, GET_ARG_COUNT(__VA_ARGS__), __VA_ARGS__) \
    appGetCudaHelper()->FreeStack(); \
}

#define _LAUN_KERNEL0(func, block, thread) \
{ \
    launchKernel((const void*)func, block, thread, (void**)NULL);  \
}


#pragma region Improve-1 launch backends

/**
* Improve-1 (multi-GPU-improve1.md 3.3): the two historical launch flavours
* with NO halo guard, kept as private raw backends. Reserved for the
* halo-exchange infrastructure itself (e.g. the gather kernel invoked from a
* CHaloManager refill, where re-entering the guard would recurse); every
* other launch must use the public macros below so the guard can
* Ensure/invalidate. These expand byte-for-byte like the pre-Improve-1
* _LAUNCH_KERNEL family on both backends.
*/
#if _CLG_LAUNCH_KERNEL
#define _CLG_LAUNCH_KERNEL_RAW _LAUN_KERNEL
#define _CLG_LAUNCH_KERNEL0_RAW _LAUN_KERNEL0
#define _CLG_LAUNCH_KERNELS_RAW _LAUN_KERNEL_SHAREMEM
#else
#define _CLG_LAUNCH_KERNEL_RAW(func, bl, th, ...) func<<<bl, th>>>(__VA_ARGS__);
#define _CLG_LAUNCH_KERNEL0_RAW(func, bl, th, ...) func<<<bl, th>>>();
#define _CLG_LAUNCH_KERNELS_RAW(func, bl, th, sh, ...) func<<<bl, th, sh>>>(__VA_ARGS__);
#endif

/**
* Public launch macros.
*
* Single GPU: plain aliases of the raw backends -- byte-for-byte the
* historical behaviour, no halo bookkeeping at all.
*
* Multi-GPU (3.3): every launch argument is inspected exactly once BEFORE the
* backend launch (CHaloManager::BeginFromArguments); a pointer argument
* hitting a registered HaloCapable buffer Ensures its halo to the full
* configured width, and every managed hit is conservatively invalidated
* (NotifyWritten) AFTER the successful launch. _LAUNCH_KERNEL0 has no
* arguments, so it stays raw. _LAUNCH_KERNEL_RANK_LOCAL (3.4.4) is for
* launches whose grid/work is intentionally rank-asymmetric; it NEVER Ensures
* (no MPI may happen there) and a HaloCapable argument fails fast.
*/
#if !_CLG_MULTI_GPU

#define _LAUNCH_KERNEL _CLG_LAUNCH_KERNEL_RAW
#define _LAUNCH_KERNEL0 _CLG_LAUNCH_KERNEL0_RAW
#define _LAUNCH_KERNELS _CLG_LAUNCH_KERNELS_RAW
#define _LAUNCH_KERNEL_RANK_LOCAL _CLG_LAUNCH_KERNEL_RAW

#elif _CLG_LAUNCH_KERNEL

//launchKernel() backend (required to compile on DCU/DTK): the guard runs in
//the free helpers defined next to _vs() below; each macro argument
//expression is still evaluated exactly once.
#define _LAUNCH_KERNEL(func, block, thread, ...) \
{ \
    const void* _clgFuncVar = (const void*)func; \
    LaunchGuarded(_clgFuncVar, block, thread, __VA_ARGS__); \
}
#define _LAUNCH_KERNELS(func, block, thread, sharemem, ...) \
{ \
    const void* _clgFuncVar = (const void*)func; \
    LaunchGuardedShareMem(_clgFuncVar, block, thread, sharemem, __VA_ARGS__); \
}
#define _LAUNCH_KERNEL_RANK_LOCAL(func, block, thread, ...) \
{ \
    const void* _clgFuncVar = (const void*)func; \
    LaunchGuardedRankLocal(_clgFuncVar, block, thread, __VA_ARGS__); \
}
#define _LAUNCH_KERNEL0(func, block, thread) _CLG_LAUNCH_KERNEL0_RAW(func, block, thread)

#else

//<<<>>> backend: keeps compile-time type checking of the kernel arguments.
//The guard runs before the launch, Commit() after it; each macro argument
//expression is bound once through the immediately-invoked generic lambda.
#define _LAUNCH_KERNEL(func, bl, th, ...) \
[&](auto&&... _clgArgs) { \
    CHaloLaunchGuard _clgGuard = appGetHaloManager()->BeginFromArguments(_clgArgs...); \
    func<<<bl, th>>>(_clgArgs...); \
    _clgGuard.Commit(); \
}(__VA_ARGS__);
#define _LAUNCH_KERNELS(func, bl, th, sh, ...) \
[&](auto&&... _clgArgs) { \
    CHaloLaunchGuard _clgGuard = appGetHaloManager()->BeginFromArguments(_clgArgs...); \
    func<<<bl, th, sh>>>(_clgArgs...); \
    _clgGuard.Commit(); \
}(__VA_ARGS__);
#define _LAUNCH_KERNEL_RANK_LOCAL(func, bl, th, ...) \
[&](auto&&... _clgArgs) { \
    CHaloLaunchGuard _clgGuard = appGetHaloManager()->BeginFromArgumentsRankLocal(_clgArgs...); \
    func<<<bl, th>>>(_clgArgs...); \
    _clgGuard.Commit(); \
}(__VA_ARGS__);
#define _LAUNCH_KERNEL0(func, bl, th, ...) _CLG_LAUNCH_KERNEL0_RAW(func, bl, th)

#endif

#pragma endregion



extern CLGAPI void launchKernel(const void* function, dim3 block, dim3 thread, void** args);
extern CLGAPI void launchKernel(const void* function, UINT block, UINT thread, void** args);
extern CLGAPI void launchKernel(const void* function, dim3 block, dim3 thread, size_t shareMem, void** args);
extern CLGAPI void launchKernel(const void* function, UINT block, UINT thread, size_t shareMem, void** args);

enum EConstIntId
{
    ECI_Dim,
    ECI_Dir,
    ECI_Lx,
    ECI_Ly,
    ECI_Lz,
    ECI_Lt,
    ECI_Volume,
    ECI_VolumeHalf,
    ECI_Volume_xyz,
    ECI_Volume_xyt,
    ECI_Volume_xzt,
    ECI_Volume_yzt,
    ECI_PlaqutteCount,
    ECI_LinkCount,
    ECI_MultX,
    ECI_MultY,
    ECI_MultZ,
    ECI_DecompX, //number of blocks
    ECI_DecompY,
    ECI_DecompZ,
    ECI_DecompLx, //threads per block (Also known as blockDim.x)
    ECI_DecompLy, //threads per block (Also known as blockDim.y)
    ECI_DecompLz, //threads per block (Also known as blockDim.z)
    ECI_GridDimZT, // ECI_Lz*ECI_Lt
    ECI_DecompAllBlock, // use one dimension decompose
    ECI_DecompAllThread,
    ECI_DecompAllBlockDir,
    ECI_DecompAllThreadDir,
    ECI_DecompAllBlockHalf, // use one dimension decompose and even-odd
    ECI_DecompAllThreadHalf,
    ECI_DecompAllBlockDirHalf,
    ECI_DecompAllThreadDirHalf,
    ECI_ThreadCountPerBlock, //thread per block
    ECI_DecompX3D, //number of blocks
    ECI_DecompY3D,
    ECI_DecompZ3D,
    ECI_DecompLx3D, //threads per block (Also known as blockDim.x)
    ECI_DecompLy3D, //threads per block (Also known as blockDim.y)
    ECI_DecompLz3D, //threads per block (Also known as blockDim.z)
    ECI_DecompX3DXYT, //number of blocks
    ECI_DecompY3DXYT,
    ECI_DecompZ3DXYT,
    ECI_DecompLx3DXYT, //threads per block (Also known as blockDim.x)
    ECI_DecompLy3DXYT, //threads per block (Also known as blockDim.y)
    ECI_DecompLz3DXYT, //threads per block (Also known as blockDim.z)
    ECI_DecompX3DXZT, //number of blocks
    ECI_DecompY3DXZT,
    ECI_DecompZ3DXZT,
    ECI_DecompLx3DXZT, //threads per block (Also known as blockDim.x)
    ECI_DecompLy3DXZT, //threads per block (Also known as blockDim.y)
    ECI_DecompLz3DXZT, //threads per block (Also known as blockDim.z)
    ECI_DecompX3DYZT, //number of blocks
    ECI_DecompY3DYZT,
    ECI_DecompZ3DYZT,
    ECI_DecompLx3DYZT, //threads per block (Also known as blockDim.x)
    ECI_DecompLy3DYZT, //threads per block (Also known as blockDim.y)
    ECI_DecompLz3DYZT, //threads per block (Also known as blockDim.z)
    ECI_RandomSeed,
    ECI_ExponentPrecision,
    ECI_ActionListLength,
    ECI_FermionFieldLength,
    ECI_MeasureListLength,
    ECI_ThreadConstaint,
    ECI_ThreadConstaintX,
    ECI_ThreadConstaintY,
    ECI_ThreadConstaintZ,
    ECI_SummationDecompose,
    ECI_UseLogADefinition, // A = U.TA() ? or A = Log(U)
    ECI_OtherGaugeField,
    ECI_GaugeFieldCount,
    ECI_BosonFieldCount,
    ECI_Tensor2FieldCount,

    ECI_Center,

    ECI_Profiler,

    ECI_MILC_StaggeredPhase, //Use t,x,y,z convention

    //Multi-GPU (Phase 1): global lattice, process grid, per-rank offset, halo.
    //These are always set; on single-GPU builds they equal the full lattice /
    //grid=[1,1,1,1] / offset=0 so existing code that reads them is safe.
    ECI_GlobalLx,  // global lattice x -- must be consecutive, used as ECI_GlobalLx+i
    ECI_GlobalLy,
    ECI_GlobalLz,
    ECI_GlobalLt,
    ECI_GpuGridX,  // process grid x -- must be consecutive, used as ECI_GpuGridX+i
    ECI_GpuGridY,
    ECI_GpuGridZ,
    ECI_GpuGridT,
    ECI_GlobalOffsetX, // this rank's sub-lattice origin in global coords, consecutive
    ECI_GlobalOffsetY,
    ECI_GlobalOffsetZ,
    ECI_GlobalOffsetT,
    ECI_HaloWidth,  // widest stencil reach (default 2 for HISQ Naik)

    ECI_ForceDWORD = 0x7fffffff,
};

enum EConstSignedIntId
{
    ECSI_CenterX,
    ECSI_CenterY,
    ECSI_CenterZ,
    ECSI_CenterT,

    ECSI_ForceDWORD = 0x7fffffff,
};

enum EConstFloatId
{
    ECF_GaugeMomentumFactor, //This is not using... remember to remove it
};

class CLGAPI CCudaHelper
{
public:
    CCudaHelper()
        : m_pDevicePtrIndexData(NULL)
        , m_pFunctionStackSpace(NULL)
        , m_uiFunctionStackSpacePointer(0)
    {
        memset(m_ConstIntegers, 0, sizeof(UINT) * kContentLength);
        memset(m_ConstSignedIntegers, 0, sizeof(INT) * kContentLength);
        memset(m_ConstFloats, 0, sizeof(Real) * kContentLength);

        //4k
        m_pFunctionStackSpace = (BYTE*)malloc(1 << 12);
        assert(NULL != m_pFunctionStackSpace);
    }
    ~CCudaHelper();

    static void DeviceQuery();
    static void MemoryQuery();

    static void DebugFunction();

    static inline UINT GetReduceDim(UINT uiLength)
    {
        UINT iRet = 0;
        while ((1U << iRet) < uiLength)
        {
            ++iRet;
        }
        return iRet;
    }

    static DOUBLE ReduceReal(DOUBLE* deviceBuffer, UINT uiLength);
    DOUBLE ReduceRealWithThreadCount(DOUBLE* deviceBuffer);
    static cuDoubleComplex ReduceComplex(cuDoubleComplex* deviceBuffer, UINT uiLength);
    cuDoubleComplex ReduceComplexWithThreadCount(cuDoubleComplex* deviceBuffer);

    void CopyConstants() const;
    void CopyRandomPointer(const class CRandom* r) const;
    void SetDeviceIndex(class CIndexData* ppIdx) const;

    class CIndexData* m_pDevicePtrIndexData;

    //we never need gamma matrix on host, so this is purely hiden in device
    void CreateGammaMatrix() const;

    void SetFieldPointers();

    /**ret[0] = max thread count, ret[1,2,3] = max thread for x,y,z per block*/
    static TArray<UINT> GetMaxThreadCountAndThreadPerblock(INT deviceId);

    UINT m_ConstIntegers[kContentLength];
    INT m_ConstSignedIntegers[kContentLength];
    Real m_ConstFloats[kContentLength];

    #pragma region global temperary buffers

    /**
    * The buffer size is NOT thread count of a block, but thread count of a grid
    * make sure this is called after thread is partitioned
    */
    void AllocateTemeraryBuffers(UINT uiThreadCount);

    void ReleaseTemeraryBuffers()
    {
        if (NULL != m_pDevicePtrIndexData)
        {
            checkCudaErrors(__cudaFree(m_pDevicePtrIndexData));
        }

        checkCudaErrors(__cudaFree(m_pRealBufferThreadCount));
        checkCudaErrors(__cudaFree(m_pComplexBufferThreadCount));

        //checkCudaErrors(cudaFree(m_pIndexBuffer));

        for (UINT i = 0; i < kMaxFieldCount; ++i)
        {
            if (NULL != m_deviceFieldPointers[i])
            {
                checkCudaErrors(__cudaFree(m_deviceFieldPointers[i]));
                m_deviceFieldPointers[i] = NULL;
            }
            if (NULL != m_deviceBoundaryFieldPointers[i])
            {
                checkCudaErrors(__cudaFree(m_deviceBoundaryFieldPointers[i]));
                m_deviceBoundaryFieldPointers[i] = NULL;
            }
        }
    }

    void ThreadBufferZero(cuDoubleComplex* pDeviceBuffer, cuDoubleComplex cInitial = make_cuDoubleComplex(0.0, 0.0)) const;
    void ThreadBufferZero(DOUBLE* pDeviceBuffer, DOUBLE fInitial = 0.0) const;

    //m_uiThreadCount = Volumn, so this is in fact volumn sum
    cuDoubleComplex ThreadBufferSum(cuDoubleComplex* pDeviceBuffer);
    DOUBLE ThreadBufferSum(DOUBLE* pDeviceBuffer);

    //struct SIndex* m_pIndexBuffer;
    cuDoubleComplex* m_pComplexBufferThreadCount;
    DOUBLE* m_pRealBufferThreadCount;

    class CField * m_deviceFieldPointers[kMaxFieldCount];
    class CFieldBoundaryParent* m_deviceBoundaryFieldPointers[kMaxFieldCount];

    //thread per grid ( = volumn)
    UINT m_uiThreadCount;
    UINT m_uiReducePower;

    BYTE* m_pFunctionStackSpace;
    UINT m_uiFunctionStackSpacePointer = 0;

    BYTE* AllocateStack(UINT size)
    {
        UINT currentPointer = m_uiFunctionStackSpacePointer;
        m_uiFunctionStackSpacePointer += size;
        return (m_pFunctionStackSpace + currentPointer);
    }

    void FreeStack()
    {
        m_uiFunctionStackSpacePointer = 0;
    }

    #pragma endregion
};

//change everything to void star
template<typename T>
void* _vs(const T& a)
{
    BYTE* ptr = appGetCudaHelper()->AllocateStack(sizeof(T));
    memcpy(ptr, &a, sizeof(T));
    return (void*)ptr;
}

#pragma region Improve-1 guarded launch helpers

#if _CLG_MULTI_GPU

/**
* Improve-1 (multi-GPU-improve1.md 3.3): launchKernel()-backend bodies of the
* public launch macros. The halo guard runs BEFORE launchKernel (Ensuring the
* halo of every registered HaloCapable buffer the arguments hit), and
* Commit() runs AFTER it (conservatively NotifyWritten-ing every recorded
* handle). FreeStack() keeps its historical position right after the launch.
* Every argument expression has already been evaluated exactly once by the
* calling macro.
*/
template<typename TBlock, typename TThread, typename... TArgs>
inline void LaunchGuarded(const void* funcvar, TBlock block, TThread thread, const TArgs&... args)
{
    CHaloLaunchGuard guard = appGetHaloManager()->BeginFromArguments(args...);
    void* pArgArray[] = { _vs(args)... };
    launchKernel(funcvar, block, thread, pArgArray);
    appGetCudaHelper()->FreeStack();
    guard.Commit();
}

template<typename TBlock, typename TThread, typename... TArgs>
inline void LaunchGuardedShareMem(const void* funcvar, TBlock block, TThread thread, size_t shareMem, const TArgs&... args)
{
    CHaloLaunchGuard guard = appGetHaloManager()->BeginFromArguments(args...);
    void* pArgArray[] = { _vs(args)... };
    launchKernel(funcvar, block, thread, shareMem, pArgArray);
    appGetCudaHelper()->FreeStack();
    guard.Commit();
}

/** 3.4.4 rank-local variant: never Ensures; a HaloCapable hit fails fast. */
template<typename TBlock, typename TThread, typename... TArgs>
inline void LaunchGuardedRankLocal(const void* funcvar, TBlock block, TThread thread, const TArgs&... args)
{
    CHaloLaunchGuard guard = appGetHaloManager()->BeginFromArgumentsRankLocal(args...);
    void* pArgArray[] = { _vs(args)... };
    launchKernel(funcvar, block, thread, pArgArray);
    appGetCudaHelper()->FreeStack();
    guard.Commit();
}

#endif

#pragma endregion

extern void CLGAPI appExistCuda();

#if _CLG_MULTI_GPU

/**
* Improve-1 I5 gate helper (single evaluation, multi-GPU-improve1.md 3.3):
* launches a trivial probe kernel through the PUBLIC _LAUNCH_KERNEL macro
* with side-effecting block / thread / argument expressions. On return,
* puiCounts[0..3] (host array) must each be exactly 1, and the kernel has
* written 7 + 8 to puiDeviceMarker[0] (device array). Exported so the MG test
* suite can drive it on either backend (launchKernel or raw <<<>>>).
*/
extern void CLGAPI appLaunchGuardSingleEvalProbe(UINT* puiCounts, UINT* puiDeviceMarker);

#endif

#define cudaSafeFree(ptr) if (NULL != ptr) { checkCudaErrors(__cudaFree(ptr)); ptr = NULL; }

__END_NAMESPACE


#endif //#ifndef _CUDAHELPER_H_

//=============================================================================
// END OF FILE
//=============================================================================