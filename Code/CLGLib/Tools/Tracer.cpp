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

__END_NAMESPACE

//====================================================================
//====================================================================
