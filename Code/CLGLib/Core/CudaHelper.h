//=============================================================================
// FILENAME : CudaHelper.h
// 
// DESCRIPTION:
// This is the file for some common CUDA usage
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CUDAHELPER_H_
#define _CUDAHELPER_H_

__BEGIN_NAMESPACE

extern "C" bool runCudaTest(const int argc, const char **argv,
    char *data, int2 *data_int2, unsigned int len);

#define __thread_id ((threadIdx.x + blockIdx.x * blockDim.x) * blockDim.y * gridDim.y * blockDim.z * gridDim.z + (threadIdx.y + blockIdx.y * blockDim.y) * blockDim.z * gridDim.z + (threadIdx.z + blockIdx.z * blockDim.z))

enum { kSharedLength = 1024, };

enum { kContentLength = 1024,};

extern __constant__ UINT _constIntegers[kContentLength];
extern __constant__ Real _constFloats[kContentLength];

extern __constant__ class CRandom* __r;
extern __constant__ class CIndex* __idx;

extern __constant__ class gammaMatrixSet* __diracGamma;
extern __constant__ class gammaMatrixSet* __chiralGamma;

extern __constant__ struct deviceSU3* __SU3Generators[8];

enum EConstIntId
{
    ECI_Dim,
    ECI_Dir,
    ECI_Lx,
    ECI_Ly,
    ECI_Lz,
    ECI_Lt,
    ECI_Volumn,
    ECI_PlaqutteCount,
    ECI_LinkCount,
    ECI_MultX,
    ECI_MultY,
    ECI_MultZ,
    ECI_DecompX, //number of blocks
    ECI_DecompY,
    ECI_DecompZ,
    ECI_DecompLx, //threads per block
    ECI_DecompLy,
    ECI_DecompLz,
    ECI_ThreadCountPerBlock, //thread per block
    ECI_ThreadCount, //thread per grid (just lx*ly*lz)
    ECI_RandomSeed,
    ECI_ExponentPrecision,
    ECI_ActionListLength,
    ECI_MeasureListLength,
    ECI_SUN,

    ECI_ForceDWORD = 0x7fffffff,
};

enum EConstFloatId
{
    ECF_InverseSqrtLink16,
};

class CLGAPI CCudaHelper
{
public:
    CCudaHelper()
    {
        memset(m_ConstIntegers, 0, sizeof(UINT) * kContentLength);
        memset(m_ConstFloats, 0, sizeof(Real) * kContentLength);
    }
    ~CCudaHelper();

    static void DeviceQuery();
    static void MemoryQuery();
    void CopyConstants() const;
    void CopyRandomPointer(const class CRandom* r) const;
    void SetDeviceIndex(class CIndex** ppIdx) const;
    //we never need gamma matrix on host, so this is purely hiden in device
    void CreateGammaMatrix() const;

    /**ret[0] = max thread count, ret[1,2,3] = max thread for x,y,z per block*/
    static TArray<UINT> GetMaxThreadCountAndThreadPerblock();


    UINT m_ConstIntegers[kContentLength];
    Real m_ConstFloats[kContentLength];

    #pragma region global temperary buffers

    //make sure this is called after thread is partioned
    void AllocateTemeraryBuffers(UINT uiThreadCount)
    {
        m_uiThreadCount = uiThreadCount;
        checkCudaErrors(cudaMalloc((void**)&m_pRealBufferThreadCount, sizeof(Real)* uiThreadCount));
        checkCudaErrors(cudaMalloc((void**)&m_pComplexBufferThreadCount, sizeof(_Complex)* uiThreadCount));
    }

    void ReleaseTemeraryBuffers()
    {
        checkCudaErrors(cudaFree(m_pRealBufferThreadCount));
        checkCudaErrors(cudaFree(m_pComplexBufferThreadCount));
    }

    _Complex ThreadBufferSum(_Complex * pDeviceBuffer);
    Real ThreadBufferSum(Real * pDeviceBuffer);

    _Complex * m_pComplexBufferThreadCount;
    Real * m_pRealBufferThreadCount;
    UINT m_uiThreadCount;

    #pragma endregion
};

__END_NAMESPACE


#endif //#ifndef _CUDAHELPER_H_

//=============================================================================
// END OF FILE
//=============================================================================