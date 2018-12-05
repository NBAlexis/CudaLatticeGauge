//=============================================================================
// FILENAME : WinFunction.cpp
// 
// DESCRIPTION:
// This is the system functions for MS-VC platform
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region Trace functions


CRITICAL_SECTION s_cs;

/**
*
*
*/
CLGAPI void appInitCriticalSection()
{
    ::InitializeCriticalSection(&s_cs);
}

/**
*
*
*/
CLGAPI void appUninitCriticalSection()
{
    ::DeleteCriticalSection(&s_cs);
}

/**
*
*
*/
CLGAPI void appEnterCriticalSection()
{
    ::EnterCriticalSection(&s_cs);
}

/**
*
*
*/
CLGAPI void appLeaveCriticalSection()
{
    ::LeaveCriticalSection(&s_cs);
}


/**
*
*
*/
CLGAPI void appTraceA(const ANSICHAR* Fmt, ...)
{
    static char buf[1024] = "";
    appEnterCriticalSection();
    va_list ArgPtr;
    va_start(ArgPtr, Fmt);
    _vsnprintf_s(buf, 1023, Fmt, ArgPtr);
    va_end(ArgPtr);
    OutputDebugStringA(buf);
    appLeaveCriticalSection();
}

/**
*
*
*/
CLGAPI void appTraceW(const UNICHAR* Fmt, ...)
{
    static wchar_t buf[1024] = L"";
    appEnterCriticalSection();
    va_list ArgPtr;
    va_start(ArgPtr, Fmt);
    _vsnwprintf_s(buf, 1023, Fmt, ArgPtr);
    va_end(ArgPtr);
    OutputDebugStringW(buf);
    appLeaveCriticalSection();
}

#pragma endregion Trace functions

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
