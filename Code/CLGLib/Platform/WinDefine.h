//=============================================================================
// FILENAME : WinDefine.h
// 
// DESCRIPTION:
// This is the data-type, defination file for MS-VC platform
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#ifndef _WINDEFINE_H_
#define _WINDEFINE_H_

#pragma region Warnings

//No using when compile with nvcc...
#pragma warning(disable : 4996) /* Something wrong with the code from bridge++ */
#pragma warning(disable : 4251) /* needs to have dll-interface to be used by clients of class (When using STL) */
#pragma warning(disable : 4819) /* The file contains a character that cannot be represented in the current code page (936) (In Cuda XML file)*/

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

typedef double					DOUBLE;
typedef unsigned __int64		QWORD;
typedef __int64					SQWORD;

// Unsigned base types.
typedef unsigned char		BYTE;		// 8-bit  unsigned.
typedef unsigned short		WORD;		// 16-bit unsigned.
typedef unsigned __int32	UINT;		// 32-bit unsigned.
typedef unsigned long		DWORD;		// 32-bit unsigned.
typedef unsigned __int64	QWORD;		// 64-bit unsigned.

                                        // Signed base types.
typedef	signed char			SBYTE;		// 8-bit  signed.
typedef signed short		SWORD;		// 16-bit signed.
typedef signed __int32  	INT;		// 32-bit signed.
typedef signed __int64		SQWORD;		// 64-bit signed.

                                        // Character types.
typedef char				ANSICHAR;	// An ANSI character.
                                        //typedef unsigned short      UNICHAR;	// A UNICODE character.
typedef wchar_t             UNICHAR;
typedef unsigned char		ANSICHARU;	// An ANSI character.
typedef unsigned short      UNICHARU;	// A UNICODE character.

                                        // Other base types.
typedef signed __int32		UBOOL;		// Boolean 0 (FALSE) or 1 (TRUE).
typedef float				FLOAT;		// 32-bit IEEE floating point.
typedef double				DOUBLE;		// 64-bit IEEE double.
                                        //TODO even undef SIZE_T not work..
                                        //typedef unsigned long       SIZE_T;     // Corresponds to C SIZE_T.

#ifdef _WIN64
typedef unsigned __int64	PTRINT;		// Integer large enough to hold a pointer.
#else
typedef unsigned long		PTRINT;		// Integer large enough to hold a pointer.
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

#ifdef _CLG_UNICODE

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

#define __REMINDER_STR(s)		#s
#define __REMINDER_STRINT(s)    __REMINDER_STR(s)
#define REMINDER(prefix, msg)	message( __FILE__ "(" __REMINDER_STRINT(__LINE__) ") : " prefix msg )
#define TODO(msg)				REMINDER("TODO: ", #msg)


#endif //#ifndef _WINDEFINE_H_

//=============================================================================
// END OF FILE
//=============================================================================