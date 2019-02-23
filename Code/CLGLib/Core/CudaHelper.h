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

//((threadIdx.x + blockIdx.x * blockDim.x) * blockDim.y * gridDim.y * blockDim.z * gridDim.z + (threadIdx.y + blockIdx.y * blockDim.y) * blockDim.z * gridDim.z + (threadIdx.z + blockIdx.z * blockDim.z))
//#define __thread_id ((threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z))

enum { kContentLength = 1024, };

extern __constant__ UINT _constIntegers[kContentLength];
extern __constant__ Real _constFloats[kContentLength];

extern __constant__ class CRandom* __r;
extern __constant__ class CIndex* __idx;

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
    CHARGECONJG,
    EGM_MAX,
    )

//extern __constant__ struct gammaMatrix __diracGamma[EGM_MAX];
extern __constant__ struct gammaMatrix __chiralGamma[EGM_MAX];

extern __constant__ struct deviceSU3 __SU3Generators[9];

enum EConstIntId
{
    ECI_Dim,
    ECI_Dir,
    ECI_Lx,
    ECI_Ly,
    ECI_Lz,
    ECI_Lt,
    ECI_Volumn,
    ECI_Volumn_xyz,
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
    ECI_ThreadCountPerBlock, //thread per block
    ECI_RandomSeed,
    ECI_ExponentPrecision,
    ECI_ActionListLength,
    ECI_FermionFieldLength,
    ECI_MeasureListLength,
    ECI_SUN,
    ECI_ThreadConstaint,
    ECI_SummationDecompose,

    ECI_ForceDWORD = 0x7fffffff,
};

enum EConstFloatId
{
    ECF_InverseSqrtLink16, //This is not using... remember to remove it
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

    static Real ReduceReal(Real* deviceBuffer, UINT uiLength);
    Real ReduceRealWithThreadCount(Real* deviceBuffer);
    static _Complex ReduceComplex(_Complex* deviceBuffer, UINT uiLength);
    _Complex ReduceComplexWithThreadCount(_Complex* deviceBuffer);

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
    void AllocateTemeraryBuffers(UINT uiThreadCount);

    void ReleaseTemeraryBuffers()
    {
        checkCudaErrors(cudaFree(m_pRealBufferThreadCount));
        checkCudaErrors(cudaFree(m_pComplexBufferThreadCount));
        //checkCudaErrors(cudaFree(m_pIndexBuffer));
    }

    _Complex ThreadBufferSum(_Complex * pDeviceBuffer);
    Real ThreadBufferSum(Real * pDeviceBuffer);

    //struct SIndex* m_pIndexBuffer;
    _Complex * m_pComplexBufferThreadCount;
    Real * m_pRealBufferThreadCount;

    //thread per grid ( = volumn)
    UINT m_uiThreadCount;
    UINT m_uiReducePower;

    #pragma endregion
};

__END_NAMESPACE


#endif //#ifndef _CUDAHELPER_H_

//=============================================================================
// END OF FILE
//=============================================================================