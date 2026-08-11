//=============================================================================
// FILENAME : CudaBuffer.h
// 
// DESCRIPTION:
// Use this to avoid memory fragment
//
// REVISION:
//  [mm/dd/yy]
//  [02/20/2019 nbale]
//=============================================================================
#pragma once

#include "Tools/Data/TArray.h"
#include "Tools/Data/THashMap.h"

#ifndef _CUDABUFFER_H_
#define _CUDABUFFER_H_

#define __cudaMalloc(pptr, size) GetBuffer()->CudaMalloc(pptr, size, _T(__FILE__), __LINE__)
#define __cudaFree(ptr) GetBuffer()->CudaFree(ptr)

#define appSafeCudaFree(a) if (NULL != a) { checkCudaErrors(cudaFree(a)); a = NULL; }

__BEGIN_NAMESPACE

extern CLGAPI void appGeneral(const TCHAR *format, ...);
inline void appFlushLog();



class CLGAPI CCudaBuffer
{
public:

    struct CLGAPI SCudaMemRecord
    {
        char m_sFilename[1024];
        size_t m_iLinenumber;
        size_t m_uiSize;
        void* m_ptr;
    };

    CCudaBuffer() 
        : m_pDevicePtr(NULL)
        , m_ulAllocated(0)
        , m_ulTotal(0)
        , m_ulFree(0)
        , m_bUseBuffer(FALSE)
    {
        
    }

    ~CCudaBuffer()
    {
        if (NULL != m_pDevicePtr)
        {
            checkCudaErrors(cudaFree(m_pDevicePtr));
        }
        appDetailed(_T("======== Unfreed Memory (start) ==========.\n"));
        TArray<ULONGLONG> allKeys = m_mapRecords.GetAllKeys();
        for (INT i = 0; i < allKeys.Num(); ++i)
        {
            const SCudaMemRecord& record = m_mapRecords[allKeys[i]];
            appDetailed(_T("Unfreed %d MB, Allocated at %s(%d).\n"),
                static_cast<UINT>(record.m_uiSize / (1ull << 20)),
                record.m_sFilename,
                static_cast<INT>(record.m_iLinenumber)
            );
        }
        appDetailed(_T("======== Unfreed Memory (end) ==========.\n"));
    }

    void Initial(FLOAT fSize)
    {
        //(1ull << 20) = 1M
        //fSize is in GB
        m_ulTotal = static_cast<ULONGLONG>(1024.f * fSize) * (1ull << 20);
        checkCudaErrors(cudaMalloc((void**)&m_pDevicePtr, m_ulTotal));
        m_ulFree = m_ulTotal;
        m_bUseBuffer = TRUE;

        appGeneral(_T("Total %d MB Allocated.\n"), static_cast<UINT>(m_ulTotal / (1ull << 20)));
    }

    inline cudaError_t CudaMalloc(void** pPtr, size_t size, const char* filename, size_t linenumber)
    {
        if (!m_bUseBuffer || size <= (1U << 12))
        {
            cudaError_t err = cudaMalloc(pPtr, size);
            if (cudaSuccess == err)
            {
                SCudaMemRecord record;
                appStrcpy(record.m_sFilename, 1023, filename);
                record.m_sFilename[1023] = 0;
                record.m_iLinenumber = linenumber;
                record.m_uiSize = size;
                record.m_ptr = *pPtr;
                m_mapRecords[reinterpret_cast<ULONGLONG>(*pPtr)] = record;

                m_ulAllocated += size;
                appDetailed(_T("Allocated %d MB, Allocated at %s(%d), Total Allocated %d MB.\n"),
                    static_cast<UINT>(size / (1ull << 20)),
                    record.m_sFilename,
                    static_cast<INT>(record.m_iLinenumber),
                    static_cast<UINT>(m_ulAllocated / (1ull << 20)));
            }

            return err;
        }

        if (m_ulFree < size)
        {
            checkCudaErrors(cudaFree(m_pDevicePtr));
            _FAIL_EXIT;
            //return cudaErrorMemoryAllocation;
        }

        (*pPtr) = m_pDevicePtr + (m_ulTotal - m_ulFree);
        m_ulFree -= size;

        SCudaMemRecord record;
        appStrcpy(record.m_sFilename, 1023, filename);
        record.m_sFilename[1023] = 0;
        record.m_iLinenumber = linenumber;
        record.m_uiSize = size;
        record.m_ptr = *pPtr;
        m_mapRecords[reinterpret_cast<ULONGLONG>(*pPtr)] = record;

        appDetailed(_T("Allocated %d MB, Allocated at %s(%d), Free memory %d MB Left.\n"),
            static_cast<UINT>(size / (1ull << 20)),
            record.m_sFilename,
            static_cast<INT>(record.m_iLinenumber),
            static_cast<UINT>(m_ulFree / (1ull << 20)));

        return cudaSuccess;
    }

    inline cudaError_t CudaFree(void* pPtr)
    {
        if (!m_bUseBuffer)
        {
            if (m_mapRecords.Exist(reinterpret_cast<ULONGLONG>(pPtr)))
            {
                const SCudaMemRecord& record = m_mapRecords[reinterpret_cast<ULONGLONG>(pPtr)];
                appDetailed(_T("Free %d MB, Allocated at %s(%d), Total Allocated %d MB.\n"),
                    static_cast<UINT>(record.m_uiSize / (1ull << 20)),
                    record.m_sFilename,
                    static_cast<INT>(record.m_iLinenumber),
                    static_cast<UINT>((m_ulAllocated - record.m_uiSize) / (1ull << 20)));
                m_ulAllocated -= record.m_uiSize;
                m_mapRecords.RemoveKey(reinterpret_cast<ULONGLONG>(pPtr));
            }
            else
            {
                appGeneral(_T("Warning: Free a pointer not allocated by CudaBuffer: %p\n"), pPtr);
            }

            return cudaFree(pPtr);
        }
        else
        {
            if (m_mapRecords.Exist(reinterpret_cast<ULONGLONG>(pPtr)))
            {
                const SCudaMemRecord& record = m_mapRecords[reinterpret_cast<ULONGLONG>(pPtr)];
                appDetailed(_T("Free %d MB, Allocated at %s(%d), Total Allocated %d MB.\n"),
                    static_cast<UINT>(record.m_uiSize / (1ull << 20)),
                    record.m_sFilename,
                    static_cast<INT>(record.m_iLinenumber),
                    static_cast<UINT>((m_ulAllocated - record.m_uiSize) / (1ull << 20)));
                m_ulAllocated -= record.m_uiSize;
                m_mapRecords.RemoveKey(reinterpret_cast<ULONGLONG>(pPtr));

                if (record.m_uiSize <= (1UL << 12))
                {
                    return cudaFree(pPtr);
                }
            }
        }
        return cudaSuccess;
    }

    BYTE* m_pDevicePtr;
    ULONGLONG m_ulAllocated;
    ULONGLONG m_ulTotal;
    ULONGLONG m_ulFree;
    UBOOL m_bUseBuffer;
    THashMap<ULONGLONG, SCudaMemRecord> m_mapRecords;
};

inline CCudaBuffer* GetBuffer();

//an interface
class CLGAPI CRegisteredBufferCache
{
public:

    CRegisteredBufferCache()
    {

    }

    virtual ~CRegisteredBufferCache()
    {

    }
};

__END_NAMESPACE


#endif //#ifndef _CUDABUFFER_H_

//=============================================================================
// END OF FILE
//=============================================================================