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

enum EVerboseLevel
{
    CRUCIAL,
    GENERAL,
    DETAILED,
    PARANOIAC,

    ForceDWORD = 0x7fffffff,
};

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
    {

    }

    ~CTracer(void)
    {
        appSafeDelete(m_pStream);
    }

    inline void Initial(EVerboseLevel eLevel, const STRING& filename = _T("stdout"))
    {
        m_eLevel = eLevel;
        if (filename == _T("stdout"))
        {
            m_pStream = new OSTREAM(COUT.rdbuf());
        }
        else
        {
            m_pStream = new OFSTREAM(filename.c_str());
        }

        if (!m_pStream)
        {
            appFailMessage(_T("ERROR: CTracer: no output stream."));
            exit(EXIT_FAILURE);
        }
    }

    inline void Print(EVerboseLevel level, const TCHAR *format, va_list& arg)
    {
        if ((level <= m_eLevel))
        {
            if (!m_pStream) 
            {
                appFailMessage(_T("ERROR: CTracer: no output stream."));
                exit(EXIT_FAILURE);
            }

            appEnterCriticalSection();

            appVsnprintf(m_cBuff, _kTraceBuffSize - 1, format, arg);

            *m_pStream << m_cBuff;
#ifdef _CLG_DEBUG
            *m_pStream << std::flush;
#endif

            appLeaveCriticalSection();
            if (!m_pStream->good()) 
            {
                appFailMessage(_T("CTracer: output failed."));
                exit(EXIT_FAILURE);
            }
        }
    }

private:

    EVerboseLevel m_eLevel;
    OSTREAM * m_pStream;
    TCHAR m_cBuff[_kTraceBuffSize];
};

extern CLGAPI void appInitialTracer(EVerboseLevel eLevel, const STRING& filename = _T("stdout"));
extern CLGAPI void appVOut(EVerboseLevel eLevel, const TCHAR *format, ...);
extern CLGAPI void _appCrucial(const TCHAR *format, ...);
extern CLGAPI void appGeneral(const TCHAR *format, ...);
extern CLGAPI void appDetailed(const TCHAR *format, ...);
extern CLGAPI void appParanoiac(const TCHAR *format, ...);

#ifdef _CLG_DEBUG
#   define appCrucial(...) {char ___msg[1024];appSprintf(___msg, sizeof(___msg), __VA_ARGS__);appTrace(_T("%s(%d): Error: %s\n"), _T(__FILE__), __LINE__, ___msg);_appCrucial(_T("%s(%d): Error: %s\n"), _T(__FILE__), __LINE__, ___msg);}
#else
#   define appCrucial(...) {_appCrucial(__VA_ARGS__);}
#endif

extern CLGAPI CTracer GTracer;

__END_NAMESPACE

#endif //_TRACER_H_

//=============================================================================
// END OF FILE
//=============================================================================
