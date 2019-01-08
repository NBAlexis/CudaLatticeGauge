//=============================================================================
// FILENAME : CLGDefine.h
// 
// DESCRIPTION:
// This is the file for some common definations
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#ifndef _CLGDEFINE_H_
#define _CLGDEFINE_H_

#pragma region Namespace

#ifdef  __NAMESPACE
#undef  __NAMESPACE
#endif
#define __NAMESPACE				CLGLib

#ifdef  __GVERSION
#undef  __GVERSION
#endif
#define __GVERSION              (0.001)

#ifdef  __BEGIN_NAMESPACE
#undef  __BEGIN_NAMESPACE
#endif
#define __BEGIN_NAMESPACE		namespace __NAMESPACE{

#ifdef  __END_NAMESPACE
#undef  __END_NAMESPACE
#endif
#define __END_NAMESPACE			}

#ifdef  __USE_NAMESPACE
#undef  __USE_NAMESPACE
#endif
#define __USE_NAMESPACE			using namespace __NAMESPACE;

#pragma endregion Namespace


#pragma region Function call

#if _CLG_WIN
# define __DLL_IMPORT			__declspec(dllimport)
# define CLGAPIPRIVATE
# define __DLL_EXPORT			__declspec(dllexport)
# define __IMPORT_LIB(libname)	comment(lib, libname)
# undef FORCEINLINE
# undef CDECL
# define FORCEINLINE             __forceinline
# define CDECL                   __cdecl

# define SUPPORTS_PRAGMA_PACK 1
# define __PACK_PUSH				pack(push, 8)
# define __PACK_POP				pack(pop)
#else
# define __DLL_IMPORT			
# define CCPRIVATE
# define __DLL_EXPORT			
# define __IMPORT_LIB(libname)	
# undef FORCEINLINE
# undef CDECL
# define FORCEINLINE             __forceinline
# define CDECL                   __cdecl

# define SUPPORTS_PRAGMA_PACK 0
# define __PACK_PUSH			
# define __PACK_POP				
#endif

#pragma endregion

#pragma region Helpers

#define ARRAY_COUNT( aarray ) \
    ( sizeof(aarray) / sizeof((aarray)[0]) )

#define appSafeDelete(p)		{if(p){delete p; p=NULL;}}
#define appSafeDeleteArray(p)	{if(p){delete[] p; p=NULL;}}

#define UN_USE(a) (void)a

#pragma endregion

#endif //#ifndef _CLGDEFINE_H_

//=============================================================================
// END OF FILE
//=============================================================================