//=============================================================================
// FILENAME : Tracer.h
// 
// DESCRIPTION:
// This is class for messages
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#ifndef _TRACER_H_
#define _TRACER_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM( EVerboseLevel,

    CRUCIAL,
    GENERAL,
    DETAILED,
    PARANOIAC,

    ForceDWORD = 0x7fffffff,

    )

enum 
{
    _kTraceBuffSize = 1024,
};

class CLGAPI CTracer
{
public:
    CTracer(void)
        : m_eLevel(CRUCIAL)
        , m_pStream(NULL)
        , m_pStdStream(NULL)
        , m_bLogDate(TRUE)
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
    }

    inline void Print(EVerboseLevel level, const TCHAR *format, va_list& arg)
    {
        if ((level <= m_eLevel))
        {
            //assert(NULL != m_pStdStream);
            if (NULL == m_pStdStream)
            {
                //Maybe the first initial is not entered?
            }
            if (m_bLogDate)
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
            *m_pStdStream << m_cBuff;
            if (NULL != m_pStream)
            {
                *m_pStream << m_cBuff;
#ifdef _CLG_DEBUG
                *m_pStream << std::flush;
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
    }

    inline void SetLogDate(UBOOL bLog) { m_bLogDate = bLog; }
    inline void SetLogHeader(const CCString& sHeader) { m_sTraceHeader = sHeader; }

private:

    EVerboseLevel m_eLevel;
    OSTREAM * m_pStream;
    OSTREAM * m_pStdStream;
    TCHAR m_cBuff[_kTraceBuffSize];
    UBOOL m_bLogDate;
    CCString m_sTraceHeader;
};

extern CLGAPI void appInitialTracer(EVerboseLevel eLevel, const CCString& filename = _T("stdout"));
extern CLGAPI void appVOut(EVerboseLevel eLevel, const TCHAR *format, ...);
extern CLGAPI void _appCrucial(const TCHAR *format, ...);
extern CLGAPI void appGeneral(const TCHAR *format, ...);
extern CLGAPI void appDetailed(const TCHAR *format, ...);
extern CLGAPI void appParanoiac(const TCHAR *format, ...);

#ifdef _CLG_DEBUG
#   define appCrucial(...) {char ___msg[1024];appSprintf(___msg, 1024, __VA_ARGS__);_appCrucial(_T("%s(%d): Error: %s\n"), _T(__FILE__), __LINE__, ___msg);}
#else
#   define appCrucial(...) {_appCrucial(__VA_ARGS__);}
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

inline void appSetLogDate(UBOOL bLog)
{
    GTracer.SetLogDate(bLog);
}

inline void appSetLogHeader(const CCString& sHeader)
{
    GTracer.SetLogHeader(sHeader);
}

__END_NAMESPACE

#endif //_TRACER_H_

//=============================================================================
// END OF FILE
//=============================================================================
