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

//aligned alloca
extern "C" void* __cdecl _alloca(SIZE_T);
#define appAlloca(size) ((0 == size) ? 0 : _alloca((size+7)&~7))

#define appInterlockedIncrement(n)	InterlockedIncrement((LONG*)(n))
#define appInterlockedDecrement(n)	InterlockedDecrement((LONG*)(n))

#pragma region String functions

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


#ifdef UNICODE
#define TCHAR_TO_ANSI(str) winToANSI((ANSICHAR*)appAlloca(winGetSizeANSI(str)),str,winGetSizeANSI(str))
#define ANSI_TO_TCHAR(str) winToUNICODE((TCHAR*)appAlloca(sizeof(UNICHAR)*winGetSizeUNICODE(str)),str,winGetSizeUNICODE(str))
#else
#define TCHAR_TO_ANSI(str) str
#define ANSI_TO_TCHAR(str) str
#endif

#pragma endregion String functions

#pragma region Math functions



#pragma endregion Math functions



#pragma region Trace and Debug


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


#pragma endregion Trace and Debug


#pragma region Time functions

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

#pragma endregion Time functions

#endif //#ifndef _WINFUNCTION_H_

//=============================================================================
// END OF FILE
//=============================================================================