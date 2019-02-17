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
#   define appStrupr	_tcsupr_s
#   define appStrlwr	_tcslwr_s
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



#pragma region Time functions

__BEGIN_NAMESPACE

FORCEINLINE UINT appGetTimeStamp(void)
{
    return static_cast<UINT>(time(0));
}

FORCEINLINE void appGetTimeNow(TCHAR* outchar, UINT buffSize)
{
    time_t now = time(0);
    ctime_s(outchar, buffSize, &now);
}

FORCEINLINE void appGetTimeUtc(TCHAR* outchar, UINT buffSize)
{
    time_t now = time(0);
    tm gmtm;
    gmtime_s(&gmtm, &now);
    asctime_s(outchar, static_cast<size_t>(buffSize), &gmtm);
}

CLGAPI INT appUnicodeToAnsi(ANSICHAR* mbstr, const UNICHAR* wcstr, INT bufsize);
CLGAPI INT appAnsiToUnicode(UNICHAR* wcstr, const ANSICHAR* mbstr, INT bufsize);

__END_NAMESPACE

#pragma endregion Time functions

#endif //#ifndef _WINFUNCTION_H_

//=============================================================================
// END OF FILE
//=============================================================================