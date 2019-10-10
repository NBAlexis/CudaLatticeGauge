//=============================================================================
// FILENAME : CudaHelper.h
// 
// DESCRIPTION:
// This is the file for CUDA Testing usage
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__constant__ UINT _constIntegers[kContentLength];
__constant__ Real _constFloats[kContentLength];
__constant__ CRandom* __r;
__constant__ CIndexData* __idx;
//why we create Dirac gamma matrix?
//__constant__ gammaMatrix __diracGamma[EGM_MAX];
__constant__ gammaMatrix __chiralGamma[EGM_MAX];
__constant__ deviceSU3 __SU3Generators[9];

__constant__ CField* __fieldPointers[_kMaxFieldCount];
__constant__ CFieldBoundary* __boundaryFieldPointers[_kMaxFieldCount];

#pragma region Kernels

/**
* The construction is on device
*/
__global__ void 
_CLG_LAUNCH_BOUND_SINGLE
_kernelCreateMatrix(/*gammaMatrix* pDirac,*/ gammaMatrix* pChiral, deviceSU3* pGenerator)
{
    gammaMatrixSet::CreateGammaMatrix(EGMS_Chiral, pChiral);
    for (int i = 0; i < 9; ++i)
    {
        pGenerator[i] = deviceSU3::makeSU3Generator(i);
    }
}

__global__ void 
_CLG_LAUNCH_BOUND_SINGLE
_kernelDebugFunction()
{
    //ULONGLONG test[4];
    //test[0] = 0;
    //test[1] = 0;
    //test[2] = 0;
    //test[3] = 0;
    //((SIndex*)&test[2])->m_uiSiteIndex = 1;
    //((SIndex*)&test[2])->m_byDir = 1;

    //printf("testres %lld\n", test[0]);
    //printf("testres %lld\n", test[1]);
    //printf("testres %lld\n", test[2]);
    //printf("testres %lld\n", test[3]);

    //printf("testres %d\n", ((SIndex*)&test[2])->m_uiSiteIndex);
    //printf("testres %lld\n", ((SIndex*)&test[2])->m_ullData);
    //__idx->DebugPrintWalkingTable();
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelDebugFunctionForSites()
{
    //intokernal;
    //printf("%d\n", uiSiteIndex);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelThreadBufferZeroReal(Real * arr, Real initial)
{
    intokernal;
    arr[uiSiteIndex] = initial;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelThreadBufferZeroComplex(CLGComplex * arr, CLGComplex initial)
{
    intokernal;
    arr[uiSiteIndex] = initial;
}

__global__ void 
_CLG_LAUNCH_BOUND
_kernelReduceReal(Real* arr, UINT uiJump, UINT uiMax)
{
    //for length 16 array
    //for jump = 1, this is 1->0, 3->2, 5->4, 7->6, 9->10, 11->10, 13->12, 15->14 
    //for jump = 2, this is 2->0, 6->4, 10->8, 14->12 
    //for jump = 4, this is 4->0, 12->8 
    //for jump = 8, this is 8->0, and is finished.

    //id target = idx * (jump << 1)
    //id from = target + jump
    UINT uiIdFrom = (threadIdx.x + blockIdx.x * blockDim.x) * (uiJump << 1) + uiJump;
    if (uiIdFrom < uiMax)
    {
        arr[uiIdFrom - uiJump] += arr[uiIdFrom];
    }
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelReduceComp(CLGComplex* arr, UINT uiJump, UINT uiMax)
{
    UINT uiIdFrom = (threadIdx.x + blockIdx.x * blockDim.x) * (uiJump << 1) + uiJump;
    if (uiIdFrom < uiMax)
    {
        arr[uiIdFrom - uiJump] = _cuCaddf(arr[uiIdFrom - uiJump], arr[uiIdFrom]);
    }
}

#pragma endregion

CCudaHelper::~CCudaHelper()
{
    ReleaseTemeraryBuffers();
}

void CCudaHelper::DeviceQuery()
{
    appGeneral(" CUDA Device Query (Runtime API) version (CUDART static linking)\n\n");

    INT deviceCount = 0;
    cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

    if (error_id != cudaSuccess) 
    {
        appGeneral("cudaGetDeviceCount returned %d\n-> %s\n",
            static_cast<INT>(error_id), cudaGetErrorString(error_id));
        appCrucial("Result = FAIL\n");
        _FAIL_EXIT;
    }

    // This function call returns 0 if there are no CUDA capable devices.
    if (deviceCount == 0) 
    {
        appGeneral("There are no available device(s) that support CUDA\n");
    }
    else 
    {
        appGeneral("Detected %d CUDA Capable device(s)\n", deviceCount);
    }

    INT dev, driverVersion = 0, runtimeVersion = 0;

    for (dev = 0; dev < deviceCount; ++dev) 
    {
        cudaSetDevice(dev);
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);

        appGeneral("\nDevice %d: \"%s\"\n", dev, deviceProp.name);

        // Console log
        cudaDriverGetVersion(&driverVersion);
        cudaRuntimeGetVersion(&runtimeVersion);
        appGeneral("  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n",
            driverVersion / 1000, (driverVersion % 100) / 10,
            runtimeVersion / 1000, (runtimeVersion % 100) / 10);
        appGeneral("  CUDA Capability Major/Minor version number:    %d.%d\n",
            deviceProp.major, deviceProp.minor);

        char msg[256];
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        appSprintf(msg, sizeof(msg),
            "  Total amount of global memory:                 %.0f MBytes "
            "(%llu bytes)\n",
            static_cast<float>(deviceProp.totalGlobalMem / 1048576.0f),
            (unsigned long long)deviceProp.totalGlobalMem);
#else
        appSprintf(msg, sizeof(msg),
            "  Total amount of global memory:                 %.0f MBytes "
            "(%llu bytes)\n",
            static_cast<float>(deviceProp.totalGlobalMem / 1048576.0f),
            (unsigned long long)deviceProp.totalGlobalMem);
#endif
        appGeneral("%s", msg);

        appGeneral("  (%2d) Multiprocessors, (%3d) CUDA Cores/MP:     %d CUDA Cores\n",
            deviceProp.multiProcessorCount,
            _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor),
            _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) *
            deviceProp.multiProcessorCount);
        appGeneral(
            "  GPU Max Clock rate:                            %.0f MHz (%0.2f "
            "GHz)\n",
            deviceProp.clockRate * 1e-3f, deviceProp.clockRate * 1e-6f);

#if CUDART_VERSION >= 5000
        // This is supported in CUDA 5.0 (runtime API device properties)
        appGeneral("  Memory Clock rate:                             %.0f Mhz\n",
            deviceProp.memoryClockRate * 1e-3f);
        appGeneral("  Memory Bus Width:                              %d-bit\n",
            deviceProp.memoryBusWidth);

        if (deviceProp.l2CacheSize) 
        {
            appGeneral("  L2 Cache Size:                                 %d bytes\n",
                deviceProp.l2CacheSize);
        }

#else
        // This only available in CUDA 4.0-4.2 (but these were only exposed in the
        // CUDA Driver API)
        int memoryClock;
        getCudaAttribute<int>(&memoryClock, CU_DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE,
            dev);
        appGeneral("  Memory Clock rate:                             %.0f Mhz\n",
            memoryClock * 1e-3f);
        int memBusWidth;
        getCudaAttribute<int>(&memBusWidth,
            CU_DEVICE_ATTRIBUTE_GLOBAL_MEMORY_BUS_WIDTH, dev);
        appGeneral("  Memory Bus Width:                              %d-bit\n",
            memBusWidth);
        int L2CacheSize;
        getCudaAttribute<int>(&L2CacheSize, CU_DEVICE_ATTRIBUTE_L2_CACHE_SIZE, dev);

        if (L2CacheSize) 
        {
            appGeneral("  L2 Cache Size:                                 %d bytes\n",
                L2CacheSize);
        }

#endif

        appGeneral(
            "  Maximum Texture Dimension Size (x,y,z)         1D=(%d), 2D=(%d, "
            "%d), 3D=(%d, %d, %d)\n",
            deviceProp.maxTexture1D, deviceProp.maxTexture2D[0],
            deviceProp.maxTexture2D[1], deviceProp.maxTexture3D[0],
            deviceProp.maxTexture3D[1], deviceProp.maxTexture3D[2]);
        appGeneral(
            "  Maximum Layered 1D Texture Size, (num) layers  1D=(%d), %d layers\n",
            deviceProp.maxTexture1DLayered[0], deviceProp.maxTexture1DLayered[1]);
        appGeneral(
            "  Maximum Layered 2D Texture Size, (num) layers  2D=(%d, %d), %d "
            "layers\n",
            deviceProp.maxTexture2DLayered[0], deviceProp.maxTexture2DLayered[1],
            deviceProp.maxTexture2DLayered[2]);

        appGeneral("  Total amount of constant memory:               %lu bytes\n",
            deviceProp.totalConstMem);
        appGeneral("  Total amount of shared memory per block:       %lu bytes\n",
            deviceProp.sharedMemPerBlock);
        appGeneral("  Total number of registers available per block: %d\n",
            deviceProp.regsPerBlock);
        appGeneral("  Warp size:                                     %d\n",
            deviceProp.warpSize);
        appGeneral("  Maximum number of threads per multiprocessor:  %d\n",
            deviceProp.maxThreadsPerMultiProcessor);
        appGeneral("  Maximum number of threads per block:           %d\n",
            deviceProp.maxThreadsPerBlock);
        appGeneral("  Max dimension size of a thread block (x,y,z): (%d, %d, %d)\n",
            deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1],
            deviceProp.maxThreadsDim[2]);
        appGeneral("  Max dimension size of a grid size    (x,y,z): (%d, %d, %d)\n",
            deviceProp.maxGridSize[0], deviceProp.maxGridSize[1],
            deviceProp.maxGridSize[2]);
        appGeneral("  Maximum memory pitch:                          %lu bytes\n",
            deviceProp.memPitch);
        appGeneral("  Texture alignment:                             %lu bytes\n",
            deviceProp.textureAlignment);
        appGeneral(
            "  Concurrent copy and kernel execution:          %s with %d copy "
            "engine(s)\n",
            (deviceProp.deviceOverlap ? "Yes" : "No"), deviceProp.asyncEngineCount);
        appGeneral("  Run time limit on kernels:                     %s\n",
            deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No");
        appGeneral("  Integrated GPU sharing Host Memory:            %s\n",
            deviceProp.integrated ? "Yes" : "No");
        appGeneral("  Support host page-locked memory mapping:       %s\n",
            deviceProp.canMapHostMemory ? "Yes" : "No");
        appGeneral("  Alignment requirement for Surfaces:            %s\n",
            deviceProp.surfaceAlignment ? "Yes" : "No");
        appGeneral("  Device has ECC support:                        %s\n",
            deviceProp.ECCEnabled ? "Enabled" : "Disabled");
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        appGeneral("  CUDA Device Driver Mode (TCC or WDDM):         %s\n",
            deviceProp.tccDriver ? "TCC (Tesla Compute Cluster Driver)"
            : "WDDM (Windows Display Driver Model)");
#endif
        appGeneral("  Device supports Unified Addressing (UVA):      %s\n",
            deviceProp.unifiedAddressing ? "Yes" : "No");
        appGeneral("  Device supports Compute Preemption:            %s\n",
            deviceProp.computePreemptionSupported ? "Yes" : "No");
        appGeneral("  Supports Cooperative Kernel Launch:            %s\n",
            deviceProp.cooperativeLaunch ? "Yes" : "No");
        appGeneral("  Supports MultiDevice Co-op Kernel Launch:      %s\n",
            deviceProp.cooperativeMultiDeviceLaunch ? "Yes" : "No");
        appGeneral("  Device PCI Domain ID / Bus ID / location ID:   %d / %d / %d\n",
            deviceProp.pciDomainID, deviceProp.pciBusID, deviceProp.pciDeviceID);

        const char *sComputeMode[] = {
            "Default (multiple host threads can use ::cudaSetDevice() with device "
            "simultaneously)",
            "Exclusive (only one host thread in one process is able to use "
            "::cudaSetDevice() with this device)",
            "Prohibited (no host thread can use ::cudaSetDevice() with this "
            "device)",
            "Exclusive Process (many threads in one process is able to use "
            "::cudaSetDevice() with this device)",
            "Unknown",
            NULL };
        appGeneral("  Compute Mode:\n");
        appGeneral("     < %s >\n", sComputeMode[deviceProp.computeMode]);
    }

    // If there are 2 or more GPUs, query to determine whether RDMA is supported
    if (deviceCount >= 2) {
        cudaDeviceProp prop[64];
        int gpuid[64];  // we want to find the first two GPUs that can support P2P
        int gpu_p2p_count = 0;

        for (int i = 0; i < deviceCount; i++) 
        {
            checkCudaErrors(cudaGetDeviceProperties(&prop[i], i));

            // Only boards based on Fermi or later can support P2P
            if ((prop[i].major >= 2)
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
                // on Windows (64-bit), the Tesla Compute Cluster driver for windows
                // must be enabled to support this
                && prop[i].tccDriver
#endif
                ) {
                // This is an array of P2P capable GPUs
                gpuid[gpu_p2p_count++] = i;
            }
        }

        // Show all the combinations of support P2P GPUs
        int can_access_peer;

        if (gpu_p2p_count >= 2) {
            for (int i = 0; i < gpu_p2p_count; i++) 
            {
                for (int j = 0; j < gpu_p2p_count; j++) 
                {
                    if (gpuid[i] == gpuid[j]) 
                    {
                        continue;
                    }
                    checkCudaErrors(
                        cudaDeviceCanAccessPeer(&can_access_peer, gpuid[i], gpuid[j]));
                    appGeneral("> Peer access from %s (GPU%d) -> %s (GPU%d) : %s\n",
                        prop[gpuid[i]].name, gpuid[i], prop[gpuid[j]].name, gpuid[j],
                        can_access_peer ? "Yes" : "No");
                }
            }
        }
    }

    // csv masterlog info
    // *****************************
    // exe and CUDA driver name
    appGeneral("\n");
    std::string sProfileString = "deviceQuery, CUDA Driver = CUDART";
    char cTemp[16];

    // driver version
    sProfileString += ", CUDA Driver Version = ";
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    appSprintf(cTemp, 10, "%d.%d", driverVersion / 1000, (driverVersion % 100) / 10);
#else
    appSprintf(cTemp, sizeof(cTemp), "%d.%d", driverVersion / 1000,
        (driverVersion % 100) / 10);
#endif
    sProfileString += cTemp;

    // Runtime version
    sProfileString += ", CUDA Runtime Version = ";
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    appSprintf(cTemp, 10, "%d.%d", runtimeVersion / 1000, (runtimeVersion % 100) / 10);
#else
    appSprintf(cTemp, sizeof(cTemp), "%d.%d", runtimeVersion / 1000,
        (runtimeVersion % 100) / 10);
#endif
    sProfileString += cTemp;

    // Device count
    sProfileString += ", NumDevs = ";
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    appSprintf(cTemp, 10, "%d", deviceCount);
#else
    appSprintf(cTemp, sizeof(cTemp), "%d", deviceCount);
#endif
    sProfileString += cTemp;
    sProfileString += "\n";
    appGeneral("%s", sProfileString.c_str());

    appGeneral("Result = PASS\n");

}

void CCudaHelper::MemoryQuery()
{
    size_t availableMemory, totalMemory, usedMemory;
    cudaMemGetInfo(&availableMemory, &totalMemory);
    usedMemory = totalMemory - availableMemory;
    appGeneral(_T("Device Memory: used %llu, available %llu, total %llu\n"), usedMemory, availableMemory, totalMemory);
}

void CCudaHelper::DebugFunction()
{
    appGeneral(_T("mult: %d, %d, %d\n"), _HC_MultX, _HC_MultY, _HC_MultZ);
    appGeneral(_T("l: %d, %d, %d, %d\n"), _HC_Lx, _HC_Ly, _HC_Lz, _HC_Lt);

    _kernelDebugFunction << <1,1 >> > ();

    preparethread;
    _kernelDebugFunctionForSites << <block, threads >> > ();
}

void CCudaHelper::CopyConstants() const
{
    checkCudaErrors(cudaMemcpyToSymbol(_constIntegers, m_ConstIntegers, sizeof(UINT) * kContentLength));
    checkCudaErrors(cudaMemcpyToSymbol(_constFloats, m_ConstFloats, sizeof(Real) * kContentLength));
}

void CCudaHelper::CopyRandomPointer(const CRandom* r) const
{
    checkCudaErrors(cudaMemcpyToSymbol(__r, &r, sizeof(CRandom*)));
}

void CCudaHelper::CreateGammaMatrix() const
{
    gammaMatrix* pChiralGamma;
    deviceSU3* pSU3;

    //create pointer
    checkCudaErrors(cudaMalloc((void**)&pChiralGamma, sizeof(gammaMatrix) * EGM_MAX));
    checkCudaErrors(cudaMalloc((void**)&pSU3, sizeof(deviceSU3) * 9));

    //craete content
    _kernelCreateMatrix << <1, 1 >> > (pChiralGamma, pSU3);

    //copy to constant
    checkCudaErrors(cudaMemcpyToSymbol(__chiralGamma, pChiralGamma, sizeof(gammaMatrix) * EGM_MAX));
    checkCudaErrors(cudaMemcpyToSymbol(__SU3Generators, pSU3, sizeof(deviceSU3) * 9));

    //free pointers (already copy to constant, no need)
    checkCudaErrors(cudaFree(pChiralGamma));
    checkCudaErrors(cudaFree(pSU3));
}

void CCudaHelper::SetDeviceIndex(class CIndexData* pIdx) const
{
    checkCudaErrors(cudaMalloc((void**)&m_pDevicePtrIndexData, sizeof(CIndexData)));
    checkCudaErrors(cudaMemcpy(m_pDevicePtrIndexData, pIdx, sizeof(CIndexData), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpyToSymbol(__idx, &m_pDevicePtrIndexData, sizeof(CIndexData*)));
}

void CCudaHelper::SetFieldPointers()
{
    for (BYTE i = 0; i < _kMaxFieldCount; ++i)
    {
        CField * pField = appGetLattice()->GetFieldById(i);
        if (NULL == pField)
        {
            m_deviceFieldPointers[i] = NULL;
        }
        else
        {
            const UINT uiSize = pField->GetClass()->GetSize();
            CField* pDeviceField = NULL;
            checkCudaErrors(cudaMalloc((void**)&pDeviceField, static_cast<size_t>(uiSize)));
            checkCudaErrors(cudaMemcpy(pDeviceField, pField, static_cast<size_t>(uiSize), cudaMemcpyHostToDevice));
            m_deviceFieldPointers[i] = pDeviceField;
        }

        CFieldBoundary * pBoundaryField = appGetLattice()->GetBoundaryFieldById(i);
        if (NULL == pBoundaryField)
        {
            m_deviceBoundaryFieldPointers[i] = NULL;
        }
        else
        {
            const UINT uiSize = pBoundaryField->GetClass()->GetSize();
            CFieldBoundary* pDeviceBoundaryField = NULL;
            checkCudaErrors(cudaMalloc((void**)&pDeviceBoundaryField, static_cast<size_t>(uiSize)));
            checkCudaErrors(cudaMemcpy(pDeviceBoundaryField, pBoundaryField, static_cast<size_t>(uiSize), cudaMemcpyHostToDevice));
            m_deviceBoundaryFieldPointers[i] = pDeviceBoundaryField;
        }
    }

    checkCudaErrors(cudaMemcpyToSymbol(__fieldPointers, m_deviceFieldPointers, sizeof(CField*) * _kMaxFieldCount));
    checkCudaErrors(cudaMemcpyToSymbol(__boundaryFieldPointers, m_deviceBoundaryFieldPointers, sizeof(CFieldBoundary*) * _kMaxFieldCount));
}

TArray<UINT> CCudaHelper::GetMaxThreadCountAndThreadPerblock()
{
    TArray<UINT> ret;

    INT deviceCount = 0;
    const cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

    if (error_id != cudaSuccess)
    {
        appCrucial("cudaGetDeviceCount returned %d\n-> %s\n",
            static_cast<INT>(error_id), cudaGetErrorString(error_id));
        appCrucial("Result = FAIL\n");
        _FAIL_EXIT;
    }

    if (0 == deviceCount)
    {
        appCrucial(_T("This program need GPU but you do NOT have a GPU.\n"));
        _FAIL_EXIT;
    }

    cudaSetDevice(0);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);

    //We need to constrain it further for shared memeory per block
    ret.AddItem(deviceProp.maxThreadsPerBlock);
    if (CCommonData::m_uiMaxThreadPerBlock > 0 && ret[0] > CCommonData::m_uiMaxThreadPerBlock)
    {
        ret[0] = CCommonData::m_uiMaxThreadPerBlock;
    }
#if _CLG_USE_LAUNCH_BOUND
    ret[0] = ret[0] > _CLG_LAUNCH_MAX_THREAD ? _CLG_LAUNCH_MAX_THREAD : ret[0];
#endif
    ret.AddItem(deviceProp.maxThreadsDim[0]);
    ret.AddItem(deviceProp.maxThreadsDim[1]);
    ret.AddItem(deviceProp.maxThreadsDim[2]);

    return ret;
}

/**
* The buffer size is NOT thread count of a block, but thread count of a grid
*/
void CCudaHelper::AllocateTemeraryBuffers(UINT uiThreadCount)
{
    m_uiThreadCount = uiThreadCount;
    m_uiReducePower = GetReduceDim((uiThreadCount + 1) >> 1);
    checkCudaErrors(cudaMalloc((void**)&m_pRealBufferThreadCount, sizeof(Real)* uiThreadCount));
    checkCudaErrors(cudaMalloc((void**)&m_pComplexBufferThreadCount, sizeof(CLGComplex)* uiThreadCount));
}

/**
* When block is not (1,1,1), some times, this will be excuted when only few block is finished.
* Must wait until all blocks finished.
*/
CLGComplex CCudaHelper::ThreadBufferSum(CLGComplex * pDeviceBuffer)
{
    return ReduceComplexWithThreadCount(pDeviceBuffer);
}

Real CCudaHelper::ThreadBufferSum(Real * pDeviceBuffer)
{
    return ReduceRealWithThreadCount(pDeviceBuffer);
}

void CCudaHelper::ThreadBufferZero(CLGComplex * pDeviceBuffer, CLGComplex cInitial) const
{
    preparethread;
    _kernelThreadBufferZeroComplex<<<block, threads>>>(pDeviceBuffer, cInitial);
}

void CCudaHelper::ThreadBufferZero(Real * pDeviceBuffer, Real fInitial) const
{
    preparethread;
    _kernelThreadBufferZeroReal << <block, threads >> >(pDeviceBuffer, fInitial);
}

Real CCudaHelper::ReduceReal(Real* deviceBuffer, UINT uiLength)
{
    const UINT iRequiredDim = (uiLength + 1) >> 1;
    const UINT iPower = GetReduceDim(iRequiredDim);
    for (UINT i = 0; i <= iPower; ++i)
    {
        UINT iJump = 1 << i;
        UINT iThreadNeeded = 1 << (iPower - i);
        UINT iBlock = iThreadNeeded > _HC_ThreadConstraint ? iThreadNeeded / _HC_ThreadConstraint : 1;
        UINT iThread = iThreadNeeded > _HC_ThreadConstraint ? _HC_ThreadConstraint : iThreadNeeded;
        _kernelReduceReal << <iBlock, iThread >> > (deviceBuffer, iJump, uiLength);
    }
    Real result[1];
    cudaMemcpy(result, deviceBuffer, sizeof(Real), cudaMemcpyDeviceToHost);
    return result[0];
}

Real CCudaHelper::ReduceRealWithThreadCount(Real* deviceBuffer)
{
    for (UINT i = 0; i <= m_uiReducePower; ++i)
    {
        UINT iJump = 1 << i;
        UINT iThreadNeeded = 1 << (m_uiReducePower - i);
        UINT iBlock = iThreadNeeded > _HC_ThreadConstraint ? iThreadNeeded / _HC_ThreadConstraint : 1;
        UINT iThread = iThreadNeeded > _HC_ThreadConstraint ? _HC_ThreadConstraint : iThreadNeeded;
        _kernelReduceReal << <iBlock, iThread >> > (deviceBuffer, iJump, m_uiThreadCount);
    }
    Real result[1];
    cudaMemcpy(result, deviceBuffer, sizeof(Real), cudaMemcpyDeviceToHost);
    return result[0];
}

CLGComplex CCudaHelper::ReduceComplex(CLGComplex* deviceBuffer, UINT uiLength)
{
    const UINT iRequiredDim = (uiLength + 1) >> 1;
    const UINT iPower = GetReduceDim(iRequiredDim);
    for (UINT i = 0; i <= iPower; ++i)
    {
        UINT iJump = 1 << i;
        UINT iThreadNeeded = 1 << (iPower - i);
        UINT iBlock = iThreadNeeded > _HC_ThreadConstraint ? iThreadNeeded / _HC_ThreadConstraint : 1;
        UINT iThread = iThreadNeeded > _HC_ThreadConstraint ? _HC_ThreadConstraint : iThreadNeeded;
        _kernelReduceComp << <iBlock, iThread >> > (deviceBuffer, iJump, uiLength);
    }
    CLGComplex result[1];
    cudaMemcpy(result, deviceBuffer, sizeof(CLGComplex), cudaMemcpyDeviceToHost);
    return result[0];
}

CLGComplex CCudaHelper::ReduceComplexWithThreadCount(CLGComplex* deviceBuffer)
{
    for (UINT i = 0; i <= m_uiReducePower; ++i)
    {
        UINT iJump = 1 << i;
        UINT iThreadNeeded = 1 << (m_uiReducePower - i);
        UINT iBlock = iThreadNeeded > _HC_ThreadConstraint ? iThreadNeeded / _HC_ThreadConstraint : 1;
        UINT iThread = iThreadNeeded > _HC_ThreadConstraint ? _HC_ThreadConstraint : iThreadNeeded;
        _kernelReduceComp << <iBlock, iThread >> > (deviceBuffer, iJump, m_uiThreadCount);
    }
    CLGComplex result[1];
    cudaMemcpy(result, deviceBuffer, sizeof(CLGComplex), cudaMemcpyDeviceToHost);
    return result[0];
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================