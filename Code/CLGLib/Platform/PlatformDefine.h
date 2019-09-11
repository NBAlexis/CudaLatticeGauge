//=============================================================================
// FILENAME : PlatformDefine.h
// 
// DESCRIPTION:
// This is the data-type, defination file for MS-VC platform
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#ifndef _PLATFORMDEFINE_H_
#define _PLATFORMDEFINE_H_

#pragma region Warnings

//No using when compile with nvcc...
#pragma warning(disable : 4819) /* The file contains a character that cannot be represented in the current code page (936) (In Cuda XML file)*/
#pragma warning(disable : 4324) /* 'struct_name' : structure was padded due to __declspec(align())*/
#pragma warning(disable : 4505) /* local unreferenced function removed */
#pragma warning(disable : 4100) /* unreferenced formal parameter */
#pragma warning(disable : 4201) /* nonstandard extension used: nameless struct/union */

//It dosen't work for intellisense warnings.
#pragma warning(disable : 26451) /* arithmetic overflow */
#pragma warning(disable : 26495) /* always initialize a member variable */

#pragma endregion Warnings

#pragma region Type definations

#undef FALSE
#undef TRUE
#undef NULL

#define FALSE   0
#define TRUE    1
#define NULL    0

// Undo any Windows defines.
#undef BYTE
#undef WORD
#undef DWORD
#undef INT
#undef Real
#undef VOID

#define VOID    void

//Type Define
//NOTE: We assume int is int32, and long long is int64, which is true both for MSVC and GCC
//We do NOT use long, long is int32 in MSVC and int64 in GCC.

typedef double                    DOUBLE;
typedef unsigned long long        QWORD;
typedef long long                SQWORD;
typedef long long               LONGLONG;
typedef unsigned long long      ULONGLONG;

// Unsigned base types.
typedef unsigned char        BYTE;        // 8-bit  unsigned.
typedef unsigned short        WORD;        // 16-bit unsigned.
typedef unsigned int        UINT;        // 32-bit unsigned.
typedef unsigned long        DWORD;        // 32-bit unsigned.
typedef unsigned long long    QWORD;        // 64-bit unsigned.

                                        // Signed base types.
typedef    signed char            SBYTE;        // 8-bit  signed.
typedef signed short        SWORD;        // 16-bit signed.
typedef signed int          INT;        // 32-bit signed.
typedef signed long long    SQWORD;        // 64-bit signed.
typedef size_t              SIZE_T;

                                        // Character types.
typedef char                ANSICHAR;    // An ANSI character.
                                        //typedef unsigned short      UNICHAR;    // A UNICODE character.
typedef short               UNICHAR;
typedef unsigned char        ANSICHARU;    // An ANSI character.
typedef unsigned short      UNICHARU;    // A UNICODE character.

                                        // Other base types.
typedef signed int           UBOOL;        // Boolean 0 (FALSE) or 1 (TRUE).
typedef float                FLOAT;        // 32-bit IEEE floating point.
typedef double                DOUBLE;        // 64-bit IEEE double.
                                        //TODO even undef SIZE_T not work..
                                        //typedef unsigned long       SIZE_T;     // Corresponds to C SIZE_T.

#ifdef _CLG_X64
typedef unsigned long long    PTRINT;        // Integer large enough to hold a pointer.
#else
typedef unsigned int        PTRINT;        // Integer large enough to hold a pointer.
#endif

#undef MAXBYTE
#undef MAXWORD
#undef MAXDWORD
#undef MAXINT

enum { MAXBYTE = 0xff };
enum { MAXWORD = 0xffffU };
enum { MAXDWORD = 0xffffffffU };
enum { MAXSBYTE = 0x7f };
enum { MAXSWORD = 0x7fff };
enum { MAXINT = 0x7fffffff };
enum { INDEX_NONE = -1 };
enum { UNICODE_BOM = 0xfeff };

#if _CLG_UNICODE

#ifndef _TCHAR_DEFINED
typedef UNICHAR  TCHAR;
typedef UNICHARU TCHARU;
#else
#undef TCHAR
#undef TCHARU
typedef UNICHAR  TCHAR;
typedef UNICHARU TCHARU;
#endif

#ifndef _TEXT_DEFINED
#define _TEXT_DEFINED
#undef TEXT
#undef _T
#define TEXT(s) L ## s
#define _T(s) TEXT(s)
#endif

#else

#ifndef _TCHAR_DEFINED
typedef ANSICHAR  TCHAR;
typedef ANSICHARU TCHARU;
#endif
#undef TEXT
#undef _T
#define TEXT(s) s
#define _T(s) TEXT(s)

#endif

#pragma endregion Type definations

#define __REMINDER_STR(s)        #s
#define __REMINDER_STRINT(s)    __REMINDER_STR(s)
#define REMINDER(prefix, msg)    message( __FILE__ "(" __REMINDER_STRINT(__LINE__) ") : " prefix msg )
#define TODO(msg)                REMINDER("TODO: ", #msg)

__BEGIN_NAMESPACE

FORCEINLINE UINT appGetTimeStamp(void)
{
    return static_cast<UINT>(time(0));
}

FORCEINLINE void appGetTimeNow(TCHAR* outchar, UINT buffSize)
{
    time_t now = time(0);
#if _CLG_WIN
    //ctime_s(outchar, buffSize, &now);
    tm now_tm;
    gmtime_s(&now_tm, &now);
    strftime(outchar, buffSize, _T("%d-%m-%Y %H-%M-%S"), &now_tm);
#else
    tm now_tm = *localtime(&now);
    strftime(outchar, buffSize, "%d-%m-%Y %H-%M-%S", &now_tm);
#endif
}

FORCEINLINE void appGetTimeUtc(TCHAR* outchar, UINT buffSize)
{
    time_t now = time(0);
    tm gmtm;
#if _CLG_WIN
    gmtime_s(&gmtm, &now);
    asctime_s(outchar, static_cast<size_t>(buffSize), &gmtm);
#else
    gmtm = *gmtime(&now);
    strftime(outchar, buffSize, "%A %c", &gmtm);
#endif
}

__END_NAMESPACE

#endif //#ifndef _PLATFORMDEFINE_H_



//=============================================================================
// END OF FILE
//=============================================================================