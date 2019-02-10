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
__constant__ CIndex* __idx;
__constant__ gammaMatrixSet* __diracGamma;
__constant__ gammaMatrixSet* __chiralGamma;
__constant__ deviceSU3 __SU3Generators[9];

/**
* The construction is on device
*/
__global__ void _kernelCreateMatrix(gammaMatrixSet** ppPtrDirac, gammaMatrixSet** ppPtrChiral, deviceSU3* pGenerator)
{
    (*ppPtrDirac) = new gammaMatrixSet(EGMS_Dirac);
    (*ppPtrChiral) = new gammaMatrixSet(EGMS_Chiral);

    for (int i = 0; i < 9; ++i)
    {
        pGenerator[i] = deviceSU3::makeSU3Generator(i);
    }
}

__global__ void _kernelDebugFunction()
{
    deviceSU3 a1 = deviceSU3::makeSU3Random(0);
    deviceSU3Vector v1 = deviceSU3Vector::makeZeroSU3Vector();

    deviceWilsonVectorSU3 v2;
    //v2.DebugPrint();
    //a1.DebugPrint();
    printf("deviceSu3: %d, deviceSU3Vector: %d, deviceWilsonVectorSU3: %d\n", 
        sizeof(deviceSU3), 
        sizeof(deviceSU3Vector), 
        sizeof(deviceWilsonVectorSU3));

    printf("deviceSU3Vector: %d, deviceWilsonVectorSU3: %d\n",
        sizeof(deviceSU3Vector),
        sizeof(deviceWilsonVectorSU3));

    printf("deviceWilsonVectorSU3: %d\n",
        sizeof(deviceWilsonVectorSU3));

    printf("SIndex: %d\n",
        sizeof(SIndex));

    for (UINT i = 0; i < 9; ++i)
    {
        __SU3Generators[i].DebugPrint();
    }
}

struct complex_plus_for_thrust
{
    __host__ __device__ _Complex operator()(const _Complex &lhs, const _Complex &rhs) const { return _cuCaddf(lhs, rhs); }
};

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
        exit(EXIT_FAILURE);
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
    _kernelDebugFunction << <1,1 >> > ();
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
    gammaMatrixSet** ppDiracGamma;
    gammaMatrixSet** ppChiralGamma;
    deviceSU3* pSU3;

    //create pointer
    checkCudaErrors(cudaMalloc((void**)&ppDiracGamma, sizeof(gammaMatrixSet*)));
    checkCudaErrors(cudaMalloc((void**)&ppChiralGamma, sizeof(gammaMatrixSet*)));
    checkCudaErrors(cudaMalloc((void**)&pSU3, sizeof(deviceSU3) * 9));

    //craete content
    _kernelCreateMatrix << <1, 1 >> > (ppDiracGamma, ppChiralGamma, pSU3);

    //copy to constant
    checkCudaErrors(cudaMemcpyToSymbol(__diracGamma, ppDiracGamma, sizeof(gammaMatrixSet*)));
    checkCudaErrors(cudaMemcpyToSymbol(__chiralGamma, ppChiralGamma, sizeof(gammaMatrixSet*)));
    checkCudaErrors(cudaMemcpyToSymbol(__SU3Generators, pSU3, sizeof(deviceSU3) * 9));

    //free pointers (already copy to constant, no need)
    checkCudaErrors(cudaFree(ppDiracGamma));
    checkCudaErrors(cudaFree(ppChiralGamma));
    checkCudaErrors(cudaFree(pSU3));
}

void CCudaHelper::SetDeviceIndex(class CIndex** ppIdx) const
{
    checkCudaErrors(cudaMemcpyToSymbol(__idx, ppIdx, sizeof(CIndex*)));
}

TArray<UINT> CCudaHelper::GetMaxThreadCountAndThreadPerblock()
{
    TArray<UINT> ret;

    INT deviceCount = 0;
    cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

    if (error_id != cudaSuccess)
    {
        appCrucial("cudaGetDeviceCount returned %d\n-> %s\n",
            static_cast<INT>(error_id), cudaGetErrorString(error_id));
        appCrucial("Result = FAIL\n");
        exit(EXIT_FAILURE);
    }

    if (0 == deviceCount)
    {
        appCrucial(_T("This program need GPU but you do NOT have a GPU.\n"));
        exit(EXIT_FAILURE);
    }

    cudaSetDevice(0);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);

    //We need to constrain it further for shared memeory per block
    ret.AddItem(deviceProp.maxThreadsPerBlock > kSharedLength ? kSharedLength : deviceProp.maxThreadsPerBlock);
    ret.AddItem(deviceProp.maxThreadsDim[0]);
    ret.AddItem(deviceProp.maxThreadsDim[1]);
    ret.AddItem(deviceProp.maxThreadsDim[2]);

    return ret;
}

void CCudaHelper::AllocateTemeraryBuffers(UINT uiThreadCount)
{
    m_uiThreadCount = uiThreadCount;
    checkCudaErrors(cudaMalloc((void**)&m_pRealBufferThreadCount, sizeof(Real)* uiThreadCount));
    checkCudaErrors(cudaMalloc((void**)&m_pComplexBufferThreadCount, sizeof(_Complex)* uiThreadCount));
    //checkCudaErrors(cudaMalloc((void**)&m_pIndexBuffer, sizeof(SIndex)* kMaxPlaqutteCache));
}

/**
* When block is not (1,1,1), some times, this will be excuted when only few block is finished.
* Must wait until all blocks finished.
*/
_Complex CCudaHelper::ThreadBufferSum(_Complex * pDeviceBuffer)
{
    //checkCudaErrors(cudaDeviceSynchronize());
    thrust::device_ptr<_Complex> dp(pDeviceBuffer);
    thrust::device_vector<_Complex> d_x(dp, dp + m_uiThreadCount);
    return thrust::reduce(d_x.begin(), d_x.end(), _make_cuComplex(0, 0), complex_plus_for_thrust());
}

Real CCudaHelper::ThreadBufferSum(Real * pDeviceBuffer)
{
    //checkCudaErrors(cudaDeviceSynchronize());
    thrust::device_ptr<Real> dp(pDeviceBuffer);
    thrust::device_vector<Real> d_x(dp, dp + m_uiThreadCount);
    return thrust::reduce(d_x.begin(), d_x.end(), (Real)0, thrust::plus<Real>());
}

////////////////////////////////////////////////////////////////////////////////
// export C interface
extern "C"
void computeGold(char *reference, char *idata, const unsigned int len);
extern "C"
void computeGold2(int2 *reference, int2 *idata, const unsigned int len);

////////////////////////////////////////////////////////////////////////////////
//! Compute reference data set
//! Each element is multiplied with the number of threads / array length
//! @param reference  reference data, computed but preallocated
//! @param idata      input data as provided to device
//! @param len        number of elements in reference / idata
////////////////////////////////////////////////////////////////////////////////
void
computeGold(char *reference, char *idata, const unsigned int len)
{
    for (unsigned int i = 0; i < len; ++i)
        reference[i] = idata[i] - 10;
}

////////////////////////////////////////////////////////////////////////////////
//! Compute reference data set for int2 version
//! Each element is multiplied with the number of threads / array length
//! @param reference  reference data, computed but preallocated
//! @param idata      input data as provided to device
//! @param len        number of elements in reference / idata
////////////////////////////////////////////////////////////////////////////////
void
computeGold2(int2 *reference, int2 *idata, const unsigned int len)
{
    for (unsigned int i = 0; i < len; ++i)
    {
        reference[i].x = idata[i].x - idata[i].y;
        reference[i].y = idata[i].y;
    }
}



__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================