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

enum 
{
    kContentLength = 256, _kMaxFieldCount = 8,
};

extern __constant__ UINT _constIntegers[kContentLength];
extern __constant__ Real _constFloats[kContentLength];

/**
* Note that, the pointers are copyied here. So the virtual functions should not be used!
*/
extern __constant__ class CField* __fieldPointers[_kMaxFieldCount];
extern __constant__ class CFieldBoundary* __boundaryFieldPointers[_kMaxFieldCount];

extern __constant__ class CRandom* __r;
extern __constant__ class CIndexData* __idx;

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
    ECI_Volume,
    ECI_Volume_xyz,
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
    ECI_ThreadConstaintX,
    ECI_ThreadConstaintY,
    ECI_ThreadConstaintZ,
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
        : m_pDevicePtrIndexData(NULL)
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
    static CLGComplex ReduceComplex(CLGComplex* deviceBuffer, UINT uiLength);
    CLGComplex ReduceComplexWithThreadCount(CLGComplex* deviceBuffer);

    void CopyConstants() const;
    void CopyRandomPointer(const class CRandom* r) const;
    void SetDeviceIndex(class CIndexData* ppIdx) const;

    class CIndexData* m_pDevicePtrIndexData;

    //we never need gamma matrix on host, so this is purely hiden in device
    void CreateGammaMatrix() const;

    void SetFieldPointers();

    /**ret[0] = max thread count, ret[1,2,3] = max thread for x,y,z per block*/
    static TArray<UINT> GetMaxThreadCountAndThreadPerblock();


    UINT m_ConstIntegers[kContentLength];
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
            checkCudaErrors(cudaFree(m_pDevicePtrIndexData));
        }

        checkCudaErrors(cudaFree(m_pRealBufferThreadCount));
        checkCudaErrors(cudaFree(m_pComplexBufferThreadCount));
        //checkCudaErrors(cudaFree(m_pIndexBuffer));

        for (UINT i = 0; i < _kMaxFieldCount; ++i)
        {
            if (NULL != m_deviceFieldPointers[i])
            {
                checkCudaErrors(cudaFree(m_deviceFieldPointers[i]));
                m_deviceFieldPointers[i] = NULL;
            }
            if (NULL != m_deviceBoundaryFieldPointers[i])
            {
                checkCudaErrors(cudaFree(m_deviceBoundaryFieldPointers[i]));
                m_deviceBoundaryFieldPointers[i] = NULL;
            }
        }
    }

    void ThreadBufferZero(CLGComplex * pDeviceBuffer, CLGComplex cInitial = _make_cuComplex(F(0.0),F(0.0))) const;
    void ThreadBufferZero(Real * pDeviceBuffer, Real fInitial = F(0.0)) const;

    CLGComplex ThreadBufferSum(CLGComplex * pDeviceBuffer);
    Real ThreadBufferSum(Real * pDeviceBuffer);

    //struct SIndex* m_pIndexBuffer;
    CLGComplex * m_pComplexBufferThreadCount;
    Real * m_pRealBufferThreadCount;
    class CField * m_deviceFieldPointers[_kMaxFieldCount];
    class CFieldBoundary * m_deviceBoundaryFieldPointers[_kMaxFieldCount];

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