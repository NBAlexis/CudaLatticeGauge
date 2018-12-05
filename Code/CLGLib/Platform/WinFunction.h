//=============================================================================
// FILENAME : WinFunction.h
// 
// DESCRIPTION:
// This is the system functions for MS-VC platform
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#ifndef _WINFUNCTION_H_
#define _WINFUNCTION_H_

#ifdef _MSC_VER
#pragma region String functions
#endif

#   define appStrstr	_tcsstr
#   define appStrcpy	_tcscpy_s
#   define appStrcat	_tcscat_s
#   define appStrlen	_tcslen
#   define appStrcmp	_tcscmp
#   define appStricmp	_tcsicmp
#   define appStrncmp   _tcsncmp
#   define appStrupr	_tcsupr
#   define appStrlwr	_tcslwr
#   define appStrchr	_tcschr
#   define appStrrchr	_tcsrchr
#   define appSprintf	_stprintf_s
#   define appVsprintf	_vstprintf_s
#   define appVsnprintf	_vsntprintf_s

#   define appToUpper	_totupper
#   define appIsSpace	_istspace
#   define appIsDigit	_istdigit
#   define appStoi		_tstoi
#   define appStof		_tstof

#ifdef _MSC_VER
#pragma endregion String functions
#endif

#ifdef _MSC_VER
#pragma region Math functions
#endif

#ifdef _MSC_VER
#pragma endregion Math functions
#endif

#ifdef _MSC_VER
#pragma region Trace and Debug
#endif

__BEGIN_NAMESPACE

CLGAPI void appInitCriticalSection();
CLGAPI void appUninitCriticalSection();
CLGAPI void appEnterCriticalSection();
CLGAPI void appLeaveCriticalSection();

//This trace is for simple debug usage, for release, use tracer instead.
CLGAPI void CDECL appTraceA(const ANSICHAR* fmt, ...);
CLGAPI void CDECL appTraceW(const UNICHAR* fmt, ...);

__END_NAMESPACE

#if _CLG_DEBUG
#   ifdef _CLG_UNICODE
#       define appTrace	appTraceW
#   else
#       define appTrace	appTraceA
#   endif
#else
#   define appTrace {}
#endif

#define appBreak()				DebugBreak()
#if defined _CLG_DEBUG
#   define appAssert(exp) {if(!(exp)){appTrace(_T("%s(%d): Assert failed: %s\n"), _T(__FILE__), __LINE__, _T(#exp)); appBreak();}}
#   define appVerify(exp) appAssert(exp)
#else
#   define appAssert(exp) { (void)(exp); }
#   define appVerify(exp) { (void)(exp); }
#endif

#ifdef _CLG_DEBUG
#	define appFailMessage(msg) {appTrace(_T("%s(%d): Error: %s\n"), _T(__FILE__), __LINE__, (msg)); appBreak();}
#else
#	define appFailMessage(msg) {}
#endif


#ifdef _MSC_VER
#pragma endregion Trace and Debug
#endif

#ifdef _MSC_VER
#pragma region Time functions
#endif

__BEGIN_NAMESPACE

FORCEINLINE DWORD appGetCycles() 
{
    LARGE_INTEGER ret;
    QueryPerformanceFrequency(&ret); 
    return ret.LowPart;
}

FORCEINLINE void appStartTimer(DWORD& timer) { timer -= appGetCycles(); }
FORCEINLINE void appEndTimer(DWORD& timer) { timer += appGetCycles(); }
FORCEINLINE FLOAT appGetTime() { return static_cast<FLOAT>(appGetCycles()) * 0.001f; }

__END_NAMESPACE

#define appGetSystemTime()		::GetTickCount()
#define appInterlockedIncrement(n)	InterlockedIncrement(CAST(LONG*,(n)))
#define appInterlockedDecrement(n)	InterlockedDecrement(CAST(LONG*,(n)))

#ifdef _MSC_VER
#pragma endregion Time functions
#endif

#endif //#ifndef _WINFUNCTION_H_

//=============================================================================
// END OF FILE
//=============================================================================