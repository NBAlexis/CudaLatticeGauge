//=============================================================================
// FILENAME : Tracer.cpp
// 
// DESCRIPTION:
// This is class for messages
//
// REVISION:
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
CLGAPI void appVOut(EVerboseLevel level, const TCHAR *format, ...)
{
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
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(CRUCIAL, format, arg);
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
    va_list arg;
    {
        va_start(arg, format);
        GTracer.Print(GENERAL, format, arg);
        va_end(arg);
    }
}

/**
*
*
*/
CLGAPI void appDetailed(const TCHAR *format, ...)
{
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
