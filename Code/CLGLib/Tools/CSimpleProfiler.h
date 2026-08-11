//=============================================================================
// FILENAME : CSimpleProfiler.h
// 
// DESCRIPTION:
// record number of calls of a function and the excution time
//
// REVISION:
//  [mm/dd/yy]
//  [01/07/2025 nbale]
//=============================================================================

#ifndef _SIMPLEPROFILER_H_
#define _SIMPLEPROFILER_H_

#if _CLG_PROFILER
#define _RECORD(name) CSimpleProfileRecord re(#name)
#define _RECORD2(name, var) CSimpleProfileRecord var(#name)
#else
#define _RECORD(name)
#define _RECORD2(name, var)
#endif

__BEGIN_NAMESPACE

class CLGAPI CSimpleProfiler
{
public:

    struct CLGAPI SSimpleProfilerRecord
    {
        ULONGLONG m_uiCalled;
        LONGLONG m_uiTimeNanoseconds;

        //to put it int tarray, i need to implement ==, but this is useles
        inline UBOOL operator==(const SSimpleProfilerRecord& Other) const
        {
            return FALSE;
        }
    };

    CSimpleProfiler()
    {

    }

    ~CSimpleProfiler()
    {
        
    }

    void AddRecord(const CCString& sName, LONGLONG uiTime)
    {
        if (m_NameHash.Exist(sName))
        {
            INT idx = m_NameHash[sName];
            m_lstRecoreds[idx].m_uiCalled = m_lstRecoreds[idx].m_uiCalled + 1;
            m_lstRecoreds[idx].m_uiTimeNanoseconds = m_lstRecoreds[idx].m_uiTimeNanoseconds + uiTime;
        }
        else
        {
            m_NameHash[sName] = m_lstRecoreds.Num();
            m_lstRecoreds.AddItem(SSimpleProfilerRecord{
                1,
                uiTime
                });
        }
    }

    void Dump() const;

    void Clear()
    {
        m_NameHash.RemoveAll();
        m_lstRecoreds.RemoveAll();
    }

    THashMap<CCString, INT> m_NameHash;
    TArray<SSimpleProfilerRecord> m_lstRecoreds;
};

extern CLGAPI CSimpleProfiler GSimpleProfiler;

class CLGAPI CSimpleProfileRecord
{
public:

    CSimpleProfileRecord(const TCHAR* name);
    ~CSimpleProfileRecord();

private:

    CCString m_sName;
    LONGLONG m_uiStart;
};

void inline appDumpProfiler()
{
    GSimpleProfiler.Dump();
}

void inline appClearProfiler()
{
    GSimpleProfiler.Clear();
}

__END_NAMESPACE

#endif //#ifndef _SIMPLEPROFILER_H_


//=============================================================================
// END OF FILE
//=============================================================================
