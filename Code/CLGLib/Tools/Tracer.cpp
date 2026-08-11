//=============================================================================
// FILENAME : Tracer.cpp
// 
// DESCRIPTION:
// This is class for messages
//
// REVISION:
//  [mm/dd/yy]
//  [12/2/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

CLGAPI CTracer GTracer;

/**
*
*
*/
CLGAPI void appInitialTracer(EVerboseLevel eLevel, const CCString& filename)
{
    GTracer.Initial(eLevel, filename);
}

/**
*
*
*/
namespace
{
    //P5-2.2: only the log rank prints (root on single-GPU builds), so an -nN
    //run produces one log instead of N interleaved copies. CRUCIAL errors are
    //also suppressed on non-root ranks (the root prints the same fatal
    //message; _FAIL_EXIT aborts every rank anyway).
    inline UBOOL _clgShouldLog()
    {
#if _CLG_MULTI_GPU
        return _CLG_IS_LOG_RANK;
#else
        return TRUE;
#endif
    }
}

CLGAPI void appVOut(EVerboseLevel level, const TCHAR *format, ...)
{
    if (!_clgShouldLog())
    {
        return;
    }
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(level, format, arg);
        va_end(arg);
    }
}

/**
*
*
*/
CLGAPI void _appCrucial(const TCHAR *format, ...)
{
    if (!_clgShouldLog())
    {
        return;
    }
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(CRUCIAL, format, arg);
        GTracer.Flush();
        va_end(arg);
    }
}

CLGAPI void _appWarning(const TCHAR* format, ...)
{
    if (!_clgShouldLog())
    {
        return;
    }
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(WARNING, format, arg);
        GTracer.Flush();
        va_end(arg);
    }
}

/**
*
*
*/
CLGAPI void appGeneral(const TCHAR *format, ...)
{
    if (!_clgShouldLog())
    {
        return;
    }
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(GENERAL, format, arg);
        va_end(arg);
    }
}

CLGAPI CCString appDressColor(EVerboseColor eColor, const TCHAR* content)
{
    static const TCHAR* colors[static_cast<INT>(EVC_Max)] = {
        _T("\033[0m"),
        _T("\033[31m"),
        _T("\033[32m"),
        _T("\033[33m"),
        _T("\033[34m"),
        _T("\033[35m"),
        _T("\033[36m"),
        _T("\033[37m"),
        _T("\033[90m")
    };

    static TCHAR tmp[static_cast<INT>(_kTraceBuffSize)];
    appSprintf(tmp, static_cast<INT>(_kTraceBuffSize), _T("%s%s\033[0m"), colors[static_cast<INT>(eColor)], content);
    return CCString(tmp);
}

/**
*
*
*/
CLGAPI void appDetailed(const TCHAR *format, ...)
{
    if (!_clgShouldLog())
    {
        return;
    }
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(DETAILED, format, arg);
        va_end(arg);
    }
}

/**
*
*
*/
CLGAPI void appParanoiac(const TCHAR *format, ...)
{
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(PARANOIAC, format, arg);
        va_end(arg);
    }
}

//CCString CTracer::CapsuleTextFile(const CCString& name, const CCString& filename, INT line)
//{
//    static TCHAR path[CCString::_CLG_MAX_PATH];
//    appGetPath(path, CCString::_CLG_MAX_PATH);
//    const CCString strpath(path);
//#if _CLG_WIN
//    const std::string hyperlinkstart = u8"\xE0\xA4\xA7";
//    const std::string hyperlinkend = u8"\xE0\xA4\x88";
//    return CCString(hyperlinkstart.c_str()) + _T("file:///") + strpath + filename + _T("#n") + appToString(line) + CCString(hyperlinkend.c_str());
//#else
//    return _T("\\e]8;;file:///") + strpath + filename + _T("#n") + appToString(line) + _T("\\e\\\\") + name + _T("\\e]8;;\\e\\\\");
//#endif
//}

__END_NAMESPACE

//====================================================================
//====================================================================
