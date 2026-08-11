//=============================================================================
// FILENAME : Tracer.h
// 
// DESCRIPTION:
// This is class for messages
//
// REVISION:
//  [mm/dd/yy]
//  [12/2/2018 nbale]
//=============================================================================

#ifndef _TRACER_H_
#define _TRACER_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM( EVerboseLevel,

    CRUCIAL,
    WARNING,
    GENERAL,
    DETAILED,
    PARANOIAC,

    ForceDWORD = 0x7fffffff,

    )

enum 
{
    _kTraceBuffSize = 4096,
};

__DEFINE_ENUM(EVerboseColor,
    EVC_RESET,
    EVC_RED,
    EVC_GREEN,
    EVC_YELLOW,
    EVC_BLUE,
    EVC_MAGENTA,
    EVC_CYAN,
    EVC_WHITE,
    EVC_Gray,
    EVC_Max
    )

class CLGAPI CTracer
{
public:
    CTracer(void)
        : m_eLevel(CRUCIAL)
        , m_pStream(NULL)
        , m_pStdStream(NULL)
        , m_uiFatalCount(0)
    {
        Initial(CRUCIAL);
    }

    ~CTracer(void)
    {
        if (NULL != m_pStream)
        {
            m_pStream->flush();
        }
        appSafeDelete(m_pStream);
    }

    inline void SetVerboseLevel(EVerboseLevel eLevel) { m_eLevel = eLevel; }

    //I10 (multi-GPU-improve1.md 3.10): --mg-config must keep stdout limited to
    //the single "rankCount gx gy gz gt" line, so all tracer output is routed to
    //stderr in that mode.
    inline void SetLogToStdErr()
    {
        appSafeDelete(m_pStdStream);
        if (NULL != m_pStream)
        {
            m_pStream->flush();
            appSafeDelete(m_pStream);
        }
        m_pStdStream = new OSTREAM(CERR.rdbuf());
        if (NULL == m_pStdStream)
        {
            printf(_T("ERROR: CTracer: no output stream."));
            exit(EXIT_FAILURE);
        }
    }

    /** Number of CRUCIAL messages printed so far (worker fatal-status summary). */
    inline UINT GetFatalCount() const { return m_uiFatalCount; }

    inline void SetOutStream(const CCString& filename = _T("stdout"))
    {
        appSafeDelete(m_pStdStream);
        if (NULL != m_pStream)
        {
            m_pStream->flush();
            appSafeDelete(m_pStream);
        }

        m_pStdStream = new OSTREAM(COUT.rdbuf());
        UBOOL bShowHasFile = FALSE;
        if (filename == _T("stdout"))
        {
            m_pStream = NULL;
        }
        else if (filename == _T("timestamp"))
        {
            CCString sRealFile;
            sRealFile.Format(_T("%d.log"), appGetTimeStamp());
            m_pStream = new OFSTREAM(sRealFile);
            bShowHasFile = TRUE;
        }
        else if (filename == _T("datetime"))
        {
            CCString sRealFile;
            static TCHAR datetime[256];
            appGetTimeNow(datetime, 256);
            sRealFile.Format(_T("%s.log"), datetime);
            m_pStream = new OFSTREAM(sRealFile);
            bShowHasFile = TRUE;
        }
        else
        {
            m_pStream = new OFSTREAM(filename);
            bShowHasFile = TRUE;
        }

        if (NULL == m_pStdStream || (bShowHasFile && NULL == m_pStream))
        {
            printf(_T("ERROR: CTracer: no output stream."));
            if (NULL != m_pStream)
            {
                m_pStream->flush();
            }
            exit(EXIT_FAILURE);
        }
    }

    inline void Initial(EVerboseLevel eLevel, const CCString& filename = _T("stdout"))
    {
        m_eLevel = eLevel;
        m_pStdStream = new OSTREAM(COUT.rdbuf());
        UBOOL bShowHasFile = FALSE;
        if (filename == _T("stdout"))
        {
            m_pStream = NULL;
        }
        else if (filename == _T("timestamp"))
        {
            CCString sRealFile;
            sRealFile.Format(_T("%d.log"), appGetTimeStamp());
            m_pStream = new OFSTREAM(sRealFile);
            bShowHasFile = TRUE;
        }
        else if (filename == _T("datetime"))
        {
            CCString sRealFile;
            static TCHAR datetime[256];
            appGetTimeNow(datetime, 256);
            sRealFile.Format(_T("%s.log"), datetime);
            m_pStream = new OFSTREAM(sRealFile);
            bShowHasFile = TRUE;
        }
        else
        {
            m_pStream = new OFSTREAM(filename);
            bShowHasFile = TRUE;
        }

        if (NULL == m_pStdStream || (bShowHasFile && NULL == m_pStream))
        {
            printf(_T("ERROR: CTracer: no output stream."));
            if (NULL != m_pStream)
            {
                m_pStream->flush();
            }
            exit(EXIT_FAILURE);
        }

        m_sErrors = _T("");
    }

    inline void Print(EVerboseLevel level, const TCHAR *format, va_list& arg)
    {
        if ((level <= m_eLevel))
        {
            //appAssert(NULL != m_pStdStream);
            if (NULL == m_pStdStream)
            {
                //Maybe the first initial is not entered?
            }
            if (CRUCIAL == level)
            {
                *m_pStdStream << _T("\033[31;1m");
                ++m_uiFatalCount;
            }
            else if (WARNING == level)
            {
                *m_pStdStream << _T("\033[33m");
            }
            else if (PARANOIAC == level)
            {
                *m_pStdStream << _T("\033[90m");
            }
            UBOOL bLogData = (m_lstLogDate.Num() > 0) ? m_lstLogDate[m_lstLogDate.Num() - 1] : TRUE;
            if (bLogData)
            {
                static TCHAR timeBuffer[256];
                if (level <= GENERAL)
                {
                    appGetTimeNow(timeBuffer, 256);
                    *m_pStdStream << _T("[") << timeBuffer << "|" << m_sTraceHeader.c_str() << _T("]");
                    if (NULL != m_pStream)
                    {
                        *m_pStream << _T("[") << timeBuffer << "|" << m_sTraceHeader.c_str() << _T("]");
                    }
                }
            }
            appVsnprintf(m_cBuff, _kTraceBuffSize - 1, format, arg);
            appAssert(NULL != m_pStdStream);
            *m_pStdStream << m_cBuff;
            if (CRUCIAL == level)
            {
                m_sErrors = m_sErrors + _T("\n") + m_cBuff;
            }
            if (CRUCIAL == level || WARNING == level || PARANOIAC == level)
            {
                *m_pStdStream << _T("\033[0m");
            }
            if (NULL != m_pStream)
            {
                *m_pStream << m_cBuff;
            }
            if (level <= DETAILED)
            {
                Flush();
            }
            else
            {
#if _CLG_DEBUG
                Flush();
#endif
            }
        }
    }

    inline void Flush() const
    {
        if (NULL != m_pStream)
        {
            m_pStream->flush();
        }
        m_pStdStream->flush();
    }

    inline void PrintAllErrors() const
    {
        *m_pStdStream << m_sErrors;
        if (NULL != m_pStream)
        {
            *m_pStream << m_sErrors;
#ifdef _CLG_DEBUG
            * m_pStream << std::flush;
#endif
        }
    }

    inline void PushLogDate(UBOOL bLog) { m_lstLogDate.PushBack(bLog); }
    inline void PopLogDate() 
    {
        if (m_lstLogDate.Num() > 0)
        {
            m_lstLogDate.Pop();
        }
    }
    inline void SetLogHeader(const CCString& sHeader) { m_sTraceHeader = sHeader; }

    //static inline CCString CapsuleURL(const CCString& url, const CCString& title)
    //{
    //    return _T("\\e]8;;") + url + _T("\\e\\\\") + title + _T("\\e]8;;\\e\\\\");
    //}

    //static CCString CapsuleTextFile(const CCString& name, const CCString& filename, INT line);

private:

    EVerboseLevel m_eLevel;
    OSTREAM * m_pStream;
    OSTREAM * m_pStdStream;
    UINT m_uiFatalCount;
    TCHAR m_cBuff[_kTraceBuffSize];
    TArray<UBOOL> m_lstLogDate;
    CCString m_sTraceHeader;
    CCString m_sErrors;
};

extern CLGAPI void appInitialTracer(EVerboseLevel eLevel, const CCString& filename = _T("stdout"));
extern CLGAPI void appVOut(EVerboseLevel eLevel, const TCHAR *format, ...);
extern CLGAPI void _appCrucial(const TCHAR *format, ...);
extern CLGAPI void _appWarning(const TCHAR* format, ...);
extern CLGAPI void appGeneral(const TCHAR *format, ...);
extern CLGAPI void appDetailed(const TCHAR *format, ...);
extern CLGAPI void appParanoiac(const TCHAR *format, ...);

extern CLGAPI CCString appDressColor(EVerboseColor eColor, const TCHAR* content);

#ifdef _CLG_DEBUG
#   define appCrucial(...) {char ___msg[1024];appSprintf(___msg, 1024, __VA_ARGS__);_appCrucial(_T("%s(%d): Error: %s\n"), _T(__FILE__), __LINE__, ___msg);}
#   define appWarning(...) {char ___msg[1024];appSprintf(___msg, 1024, __VA_ARGS__);_appWarning(_T("%s(%d): Warning: %s\n"), _T(__FILE__), __LINE__, ___msg);}
#else
#   define appCrucial(...) {_appCrucial(__VA_ARGS__);}
#   define appWarning(...) {_appWarning(__VA_ARGS__);}
#endif

extern CLGAPI CTracer GTracer;

inline void appSetTracer(EVerboseLevel eLevel, const CCString& filename)
{
    GTracer.SetVerboseLevel(eLevel);
    GTracer.SetOutStream(filename);
}

inline void appFlushLog()
{
    GTracer.Flush();
}

inline void appPopLogDate()
{
    GTracer.PopLogDate();
}

inline void appPushLogDate(UBOOL bLog)
{
    GTracer.PushLogDate(bLog);
}

inline void appSetLogHeader(const CCString& sHeader)
{
    GTracer.SetLogHeader(sHeader);
}

inline void appPrintAllErrors()
{
    GTracer.PrintAllErrors();
}

__END_NAMESPACE

#endif //_TRACER_H_

//=============================================================================
// END OF FILE
//=============================================================================
