//=============================================================================
// FILENAME : CudaBuffer.h
// 
// DESCRIPTION:
// Use this to avoid memory fragment
//
// REVISION:
//  [02/20/2019 nbale]
//=============================================================================

#ifndef _CUDABUFFER_H_
#define _CUDABUFFER_H_

#define __cudaMalloc(a, b) GetBuffer()->CudaMalloc(a, b)
#define __cudaFree(a) GetBuffer()->CudaFree(a)

__BEGIN_NAMESPACE

extern CLGAPI void appGeneral(const TCHAR *format, ...);

class CLGAPI CCudaBuffer
{
public:
    CCudaBuffer() 
        : m_pDevicePtr(NULL)
        , m_bUseBuffer(FALSE)
        , m_ulTotal(0)
        , m_ulFree(0)
    {

    }

    ~CCudaBuffer()
    {
        if (NULL != m_pDevicePtr)
        {
            checkCudaErrors(cudaFree(m_pDevicePtr));
        }
    }

    void Initial(FLOAT fSize)
    {
        m_ulTotal = static_cast<ULONGLONG>(fSize * 1024.f) * (1ull << 20);
        checkCudaErrors(cudaMalloc((void**)&m_pDevicePtr, m_ulTotal));
        m_ulFree = m_ulTotal;
        m_bUseBuffer = TRUE;
        appGeneral(_T("Total %d MB Allocated.\n"), static_cast<UINT>(m_ulTotal / (1ull << 20)));
    }

    inline cudaError_t CudaMalloc(void** pPtr, size_t size)
    {
        if (!m_bUseBuffer)
        {
            return cudaMalloc(pPtr, size);
        }

        if (m_ulFree < size)
        {
            checkCudaErrors(cudaFree(m_pDevicePtr));
            _FAIL_EXIT;
            //return cudaErrorMemoryAllocation;
        }

        (*pPtr) = m_pDevicePtr + (m_ulTotal - m_ulFree);
        m_ulFree -= size;
#if _CLG_DEBUG
        appGeneral(_T("Allocated %d MB, Free memory %d MB Left.\n"), 
            static_cast<UINT>(size / (1ull << 20)), 
            static_cast<UINT>(m_ulFree / (1ull << 20)));
#endif
        return cudaSuccess;
    }

    inline cudaError_t CudaFree(void* pPtr)
    {
        if (!m_bUseBuffer)
        {
            return cudaFree(pPtr);
        }
        return cudaSuccess;
    }

    BYTE* m_pDevicePtr;
    ULONGLONG m_ulTotal;
    ULONGLONG m_ulFree;
    UBOOL m_bUseBuffer;
};

inline CCudaBuffer* GetBuffer();

__END_NAMESPACE


#endif //#ifndef _CUDABUFFER_H_

//=============================================================================
// END OF FILE
//=============================================================================