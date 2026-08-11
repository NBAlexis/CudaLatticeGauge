//=============================================================================
// FILENAME : CLGDefine.h
// 
// DESCRIPTION:
// This is the file for some common definations
//
// REVISION:
//  [mm/dd/yy]
//  [12/2/2018 nbale]
//=============================================================================
#pragma once

#ifndef _CLGDEFINE_H_
#define _CLGDEFINE_H_

#pragma region Namespace

#ifdef  __NAMESPACE
#undef  __NAMESPACE
#endif
#define __NAMESPACE                CLGLib

#ifdef  __GVERSION
#undef  __GVERSION
#endif
#define __GVERSION      (1)

#ifdef  __GVERSION_S
#undef  __GVERSION_S
#endif
#define __GVERSION_S    (3)

#ifdef  __BEGIN_NAMESPACE
#undef  __BEGIN_NAMESPACE
#endif
#define __BEGIN_NAMESPACE        namespace __NAMESPACE{

#ifdef  __END_NAMESPACE
#undef  __END_NAMESPACE
#endif
#define __END_NAMESPACE            }

#ifdef  __USE_NAMESPACE
#undef  __USE_NAMESPACE
#endif
#define __USE_NAMESPACE            using namespace __NAMESPACE;

#pragma endregion Namespace


#pragma region Function call

#if _CLG_WIN
# define __DLL_IMPORT            __declspec(dllimport)
# define CLGAPIPRIVATE
# define __DLL_EXPORT            __declspec(dllexport)
# define __IMPORT_LIB(libname)    comment(lib, libname)
# undef FORCEINLINE
# undef CDECL
# define FORCEINLINE             __forceinline
# define FORCE_NOT_INLINE        __declspec(noinline)
# define CDECL                   __cdecl

# define SUPPORTS_PRAGMA_PACK 1
# define __PACK_PUSH                pack(push, 8)
# define __PACK_POP                pack(pop)
#else
# define __DLL_IMPORT            
# define CLGAPIPRIVATE
# define __DLL_EXPORT            __attribute__((visibility("default")))    
# define __IMPORT_LIB(libname)    
# undef FORCEINLINE
# undef CDECL
# define FORCEINLINE             inline
# define FORCE_NOT_INLINE        __attribute__((noinline))
# define CDECL                   __cdecl

# define SUPPORTS_PRAGMA_PACK 0
# define __PACK_PUSH            
# define __PACK_POP                
#endif

#pragma endregion

#pragma region Helpers

#define ARRAY_COUNT( aarray ) \
    ( sizeof(aarray) / sizeof((aarray)[0]) )

#define appSafeFree(p)        {if(p){free(p); p=NULL;}}
#define appSafeDelete(p)        {if(p){delete p; p=NULL;}}
#define appSafeDeleteArray(p)    {if(p){delete[] p; p=NULL;}}

#define UN_USE(a) (void)a

//aligned alloca
//extern "C" void* __cdecl _alloca(size_t);
#define appAlloca(size) ((0 == size) ? 0 : alloca((size+7)&~7))

#if _CLG_ENABLE_RELEASE_ASSERT
    #define appAssert(cond) assert(cond)
    #define appDebugBreak assert(0)
#else
    #if _CLG_DEBUG
        #define appAssert(cond) assert(cond)
        #define appDebugBreak assert(0)
    #else
        #define appAssert(cond)
        #define appDebugBreak
    #endif
#endif

#pragma endregion

//from https://github.com/netcan/recipes/blob/master/cpp/metaproggramming/MetaMacro.hpp
/*************************************************************************
    > File Name: MetaMacro.hpp
    > Author: Netcan
    > Descripton: MetaMacro
    > Blog: https://netcan.github.io/
    > Mail: 1469709759@qq.com
    > Created Time: 2020-08-30 08:56
************************************************************************/
#define PP_THIRD_ARG(a, b, c, ...) c
#define VA_OPT_SUPPORTED_I(...) PP_THIRD_ARG(__VA_OPT__(, ), 1, 0, )
#define VA_OPT_SUPPORTED VA_OPT_SUPPORTED_I(?)

// Traditional MSVC requires a special EXPAND phase
#if (defined(_MSC_VER) && !defined(_MSVC_TRADITIONAL)) ||                  \
    (defined(_MSVC_TRADITIONAL) && _MSVC_TRADITIONAL)

#define GET_ARG_COUNT(...)                                                 \
    INTERNAL_EXPAND_ARGS_PRIVATE(INTERNAL_ARGS_AUGMENTER(__VA_ARGS__))

#define INTERNAL_ARGS_AUGMENTER(...) unused, __VA_ARGS__
#define INTERNAL_EXPAND(x) x
#define INTERNAL_EXPAND_ARGS_PRIVATE(...)                                  \
    INTERNAL_EXPAND(INTERNAL_GET_ARG_COUNT_PRIVATE(                        \
        __VA_ARGS__, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90, 89, 88,  \
        87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74, 73, 72,    \
        71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56,    \
        55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40,    \
        39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24,    \
        23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7,   \
        6, 5, 4, 3, 2, 1, 0))

#else  // Other compilers

#if VA_OPT_SUPPORTED  // Standardized in C++20
#define GET_ARG_COUNT(...)                                                 \
    INTERNAL_GET_ARG_COUNT_PRIVATE(                                        \
        unused __VA_OPT__(, ) __VA_ARGS__, 100, 99, 98, 97, 96, 95, 94,    \
        93, 92, 91, 90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78,    \
        77, 76, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62,    \
        61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46,    \
        45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30,    \
        29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14,    \
        13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)
#elif defined(__GNUC__)  // Extension in GCC/Clang
#define GET_ARG_COUNT(...)                                                 \
    INTERNAL_GET_ARG_COUNT_PRIVATE(                                        \
        unused, ##__VA_ARGS__, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91,    \
        90, 89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75,    \
        74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59,    \
        58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43,    \
        42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27,    \
        26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11,    \
        10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)
#else  // GET_ARG_COUNT() may return 1 here
#define GET_ARG_COUNT(...)                                                 \
    INTERNAL_GET_ARG_COUNT_PRIVATE(                                        \
        unused, __VA_ARGS__, 100, 99, 98, 97, 96, 95, 94, 93, 92, 91, 90,  \
        89, 88, 87, 86, 85, 84, 83, 82, 81, 80, 79, 78, 77, 76, 75, 74,    \
        73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58,    \
        57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42,    \
        41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26,    \
        25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, \
        8, 7, 6, 5, 4, 3, 2, 1, 0)
#endif

#endif

#define INTERNAL_GET_ARG_COUNT_PRIVATE(                                    \
    e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15,  \
    e16, e17, e18, e19, e20, e21, e22, e23, e24, e25, e26, e27, e28, e29,  \
    e30, e31, e32, e33, e34, e35, e36, e37, e38, e39, e40, e41, e42, e43,  \
    e44, e45, e46, e47, e48, e49, e50, e51, e52, e53, e54, e55, e56, e57,  \
    e58, e59, e60, e61, e62, e63, e64, e65, e66, e67, e68, e69, e70, e71,  \
    e72, e73, e74, e75, e76, e77, e78, e79, e80, e81, e82, e83, e84, e85,  \
    e86, e87, e88, e89, e90, e91, e92, e93, e94, e95, e96, e97, e98, e99,  \
    e100, count, ...)                                                      \
    count


#define EXPAND(...) __VA_ARGS__
#define TMPARG(...) <__VA_ARGS__>
#define ARG_SPLIT(args) EXPAND args


//===================== The following is for SUN gauge only =====================
#if _CLG_SU4_GAUGE
#define _DEF_F2_4(imp, b, ...) imp(4, b);
#else
#define _DEF_F2_4(imp, b, ...)
#endif

#if _CLG_SU5_GAUGE
#define _DEF_F2_5(imp, b, ...) imp(5, b); EXPAND(_DEF_F2_4(imp, __VA_ARGS__))
#else
#define _DEF_F2_5(imp, b, ...) EXPAND(_DEF_F2_4(imp, __VA_ARGS__))
#endif

#if _CLG_SU6_GAUGE
#define _DEF_F2_6(imp, b, ...) imp(6, b); EXPAND(_DEF_F2_5(imp, __VA_ARGS__))
#else
#define _DEF_F2_6(imp, b, ...) EXPAND(_DEF_F2_5(imp, __VA_ARGS__))
#endif

#if _CLG_SU7_GAUGE
#define _DEF_F2_7(imp, b, ...) imp(7, b); EXPAND(_DEF_F2_6(imp, __VA_ARGS__))
#else
#define _DEF_F2_7(imp, b, ...) EXPAND(_DEF_F2_6(imp, __VA_ARGS__))
#endif

#if _CLG_SU8_GAUGE
#define _DEF_F2_8(imp, b, ...) imp(8, b); EXPAND(_DEF_F2_7(imp, __VA_ARGS__))
#else
#define _DEF_F2_8(imp, b, ...) EXPAND(_DEF_F2_7(imp, __VA_ARGS__))
#endif

#if _CLG_SU9_GAUGE
#define _DEF_F2_9(imp, b, ...) imp(9, b); EXPAND(_DEF_F2_8(imp, __VA_ARGS__))
#else
#define _DEF_F2_9(imp, b, ...) EXPAND(_DEF_F2_8(imp, __VA_ARGS__))
#endif

#if _CLG_SU10_GAUGE
#define _DEF_F2_10(imp, b, ...) imp(10, b); EXPAND(_DEF_F2_9(imp, __VA_ARGS__))
#else
#define _DEF_F2_10(imp, b, ...) EXPAND(_DEF_F2_9(imp, __VA_ARGS__))
#endif

#if _CLG_SU11_GAUGE
#define _DEF_F2_11(imp, b, ...) imp(11, b); EXPAND(_DEF_F2_10(imp, __VA_ARGS__))
#else
#define _DEF_F2_11(imp, b, ...) EXPAND(_DEF_F2_10(imp, __VA_ARGS__))
#endif

#if _CLG_SU12_GAUGE
#define _DEF_F2_12(imp, b, ...) imp(12, b); EXPAND(_DEF_F2_11(imp, __VA_ARGS__))
#else
#define _DEF_F2_12(imp, b, ...) EXPAND(_DEF_F2_11(imp, __VA_ARGS__))
#endif

#if _CLG_SU13_GAUGE
#define _DEF_F2_13(imp, b, ...) imp(13, b); EXPAND(_DEF_F2_12(imp, __VA_ARGS__))
#else
#define _DEF_F2_13(imp, b, ...) EXPAND(_DEF_F2_12(imp, __VA_ARGS__))
#endif

#if _CLG_SU14_GAUGE
#define _DEF_F2_14(imp, b, ...) imp(14, b); EXPAND(_DEF_F2_13(imp, __VA_ARGS__))
#else
#define _DEF_F2_14(imp, b, ...) EXPAND(_DEF_F2_13(imp, __VA_ARGS__))
#endif

#if _CLG_SU15_GAUGE
#define _DEF_F2_15(imp, b, ...) imp(15, b); EXPAND(_DEF_F2_14(imp, __VA_ARGS__))
#else
#define _DEF_F2_15(imp, b, ...) EXPAND(_DEF_F2_14(imp, __VA_ARGS__))
#endif

#if _CLG_SU16_GAUGE
#define _DEF_F2_16(imp, b, ...) imp(16, b); EXPAND(_DEF_F2_15(imp, __VA_ARGS__))
#else
#define _DEF_F2_16(imp, b, ...) EXPAND(_DEF_F2_15(imp, __VA_ARGS__))
#endif

#define _DEF_F2_N_IMPL(n, args) _DEF_F2_##n args
#define _DEF_F2_N(n, ...)       _DEF_F2_N_IMPL(n, (__VA_ARGS__))

//use it like #define _DEF_F2_TO(n, imp) _DEF_F2_N(n, imp, 4, 3, 2)


#if _CLG_SU4_GAUGE
#define _DEF_F_4(imp) imp(4);
#else
#define _DEF_F_4(imp)
#endif

#if _CLG_SU5_GAUGE
#define _DEF_F_5(imp) imp(5); _DEF_F_4(imp)
#else
#define _DEF_F_5(imp) _DEF_F_4(imp)
#endif

#if _CLG_SU6_GAUGE
#define _DEF_F_6(imp) imp(6); _DEF_F_5(imp)
#else
#define _DEF_F_6(imp) _DEF_F_5(imp)
#endif

#if _CLG_SU7_GAUGE
#define _DEF_F_7(imp) imp(7); _DEF_F_6(imp)
#else
#define _DEF_F_7(imp) _DEF_F_6(imp)
#endif

#if _CLG_SU8_GAUGE
#define _DEF_F_8(imp) imp(8); _DEF_F_7(imp)
#else
#define _DEF_F_8(imp) _DEF_F_7(imp)
#endif

#if _CLG_SU9_GAUGE
#define _DEF_F_9(imp) imp(9); _DEF_F_8(imp)
#else
#define _DEF_F_9(imp) _DEF_F_8(imp)
#endif

#if _CLG_SU10_GAUGE
#define _DEF_F_10(imp) imp(10); _DEF_F_9(imp)
#else
#define _DEF_F_10(imp) _DEF_F_9(imp)
#endif

#if _CLG_SU11_GAUGE
#define _DEF_F_11(imp) imp(11); _DEF_F_10(imp)
#else
#define _DEF_F_11(imp) _DEF_F_10(imp)
#endif

#if _CLG_SU12_GAUGE
#define _DEF_F_12(imp) imp(12); _DEF_F_11(imp)
#else
#define _DEF_F_12(imp) _DEF_F_11(imp)
#endif

#if _CLG_SU13_GAUGE
#define _DEF_F_13(imp) imp(13); _DEF_F_12(imp)
#else
#define _DEF_F_13(imp) _DEF_F_12(imp)
#endif

#if _CLG_SU14_GAUGE
#define _DEF_F_14(imp) imp(14); _DEF_F_13(imp)
#else
#define _DEF_F_14(imp) _DEF_F_13(imp)
#endif

#if _CLG_SU15_GAUGE
#define _DEF_F_15(imp) imp(15); _DEF_F_14(imp)
#else
#define _DEF_F_15(imp) _DEF_F_14(imp)
#endif

#if _CLG_SU16_GAUGE
#define _DEF_F_16(imp) imp(16); _DEF_F_15(imp)
#else
#define _DEF_F_16(imp) _DEF_F_15(imp)
#endif

#define _DEF_F_N(n, imp) _DEF_F_##n(imp)

// ==================== Z_N dispatch chain (parallel to _DEF_F_N) ====================
#if _CLG_Z2_GAUGE
#define _DEF_F_Z2(imp) imp(2, Z2);
#else
#define _DEF_F_Z2(imp)
#endif

#if _CLG_Z3_GAUGE
#define _DEF_F_Z3(imp) imp(3, Z3); _DEF_F_Z2(imp)
#else
#define _DEF_F_Z3(imp) _DEF_F_Z2(imp)
#endif

#if _CLG_Z4_GAUGE
#define _DEF_F_Z4(imp) imp(4, Z4); _DEF_F_Z3(imp)
#else
#define _DEF_F_Z4(imp) _DEF_F_Z3(imp)
#endif

#if _CLG_Z5_GAUGE
#define _DEF_F_Z5(imp) imp(5, Z5); _DEF_F_Z4(imp)
#else
#define _DEF_F_Z5(imp) _DEF_F_Z4(imp)
#endif

#if _CLG_Z6_GAUGE
#define _DEF_F_Z6(imp) imp(6, Z6); _DEF_F_Z5(imp)
#else
#define _DEF_F_Z6(imp) _DEF_F_Z5(imp)
#endif

#define _DEF_F_ZN(n, imp) _DEF_F_Z##n(imp)

// ==================== D_N dispatch chain (parallel to Z_N) ====================
#if _CLG_D3_GAUGE
#define _DEF_F_D3(imp) imp(3, D3);
#else
#define _DEF_F_D3(imp)
#endif

#if _CLG_D4_GAUGE
#define _DEF_F_D4(imp) imp(4, D4); _DEF_F_D3(imp)
#else
#define _DEF_F_D4(imp) _DEF_F_D3(imp)
#endif

#if _CLG_D5_GAUGE
#define _DEF_F_D5(imp) imp(5, D5); _DEF_F_D4(imp)
#else
#define _DEF_F_D5(imp) _DEF_F_D4(imp)
#endif

#if _CLG_D6_GAUGE
#define _DEF_F_D6(imp) imp(6, D6); _DEF_F_D5(imp)
#else
#define _DEF_F_D6(imp) _DEF_F_D5(imp)
#endif

#if _CLG_D7_GAUGE
#define _DEF_F_D7(imp) imp(7, D7); _DEF_F_D6(imp)
#else
#define _DEF_F_D7(imp) _DEF_F_D6(imp)
#endif

#if _CLG_D8_GAUGE
#define _DEF_F_D8(imp) imp(8, D8); _DEF_F_D7(imp)
#else
#define _DEF_F_D8(imp) _DEF_F_D7(imp)
#endif

#define _DEF_F_DN(n, imp) _DEF_F_D##n(imp)


#if _CLG_SL2C_GAUGE
#define _DEF_F2SLNC_2(imp, b, ...) imp(2, b);
#define _DEF_FSLNC_2(imp) imp(2)
#else
#define _DEF_F2SLNC_2(imp, b, ...)
#define _DEF_FSLNC_2(imp)
#endif

#if _CLG_SL3C_GAUGE
#define _DEF_F2SLNC_3(imp, b, ...) imp(3, b); EXPAND(_DEF_F2SLNC_2(imp, __VA_ARGS__))
#define _DEF_FSLNC_3(imp) imp(3); _DEF_FSLNC_2(imp)
#else
#define _DEF_F2SLNC_3(imp, b, ...) EXPAND(_DEF_F2SLNC_2(imp, __VA_ARGS__))
#define _DEF_FSLNC_3(imp) _DEF_FSLNC_2(imp)
#endif

#if _CLG_SL4C_GAUGE
#define _DEF_F2SLNC_4(imp, b, ...) imp(4, b); EXPAND(_DEF_F2SLNC_3(imp, __VA_ARGS__))
#define _DEF_FSLNC_4(imp) imp(4); _DEF_FSLNC_3(imp)
#else
#define _DEF_F2SLNC_4(imp, b, ...) EXPAND(_DEF_F2SLNC_3(imp, __VA_ARGS__))
#define _DEF_FSLNC_4(imp) _DEF_FSLNC_3(imp)
#endif

#if _CLG_SL5C_GAUGE
#define _DEF_F2SLNC_5(imp, b, ...) imp(5, b); EXPAND(_DEF_F2SLNC_4(imp, __VA_ARGS__))
#define _DEF_FSLNC_5(imp) imp(5); _DEF_FSLNC_4(imp)
#else
#define _DEF_F2SLNC_5(imp, b, ...) EXPAND(_DEF_F2SLNC_4(imp, __VA_ARGS__))
#define _DEF_FSLNC_5(imp) _DEF_FSLNC_4(imp)
#endif

#if _CLG_SL6C_GAUGE
#define _DEF_F2SLNC_6(imp, b, ...) imp(6, b); EXPAND(_DEF_F2SLNC_5(imp, __VA_ARGS__))
#define _DEF_FSLNC_6(imp) imp(6); _DEF_FSLNC_5(imp)
#else
#define _DEF_F2SLNC_6(imp, b, ...) EXPAND(_DEF_F2SLNC_5(imp, __VA_ARGS__))
#define _DEF_FSLNC_6(imp) _DEF_FSLNC_5(imp)
#endif

#if _CLG_SL7C_GAUGE
#define _DEF_F2SLNC_7(imp, b, ...) imp(7, b); EXPAND(_DEF_F2SLNC_6(imp, __VA_ARGS__))
#define _DEF_FSLNC_7(imp) imp(7); _DEF_FSLNC_6(imp)
#else
#define _DEF_F2SLNC_7(imp, b, ...) EXPAND(_DEF_F2SLNC_6(imp, __VA_ARGS__))
#define _DEF_FSLNC_7(imp) _DEF_FSLNC_6(imp)
#endif

#if _CLG_SL8C_GAUGE
#define _DEF_F2SLNC_8(imp, b, ...) imp(8, b); EXPAND(_DEF_F2SLNC_7(imp, __VA_ARGS__))
#define _DEF_FSLNC_8(imp) imp(8); _DEF_FSLNC_7(imp)
#else
#define _DEF_F2SLNC_8(imp, b, ...) EXPAND(_DEF_F2SLNC_7(imp, __VA_ARGS__))
#define _DEF_FSLNC_8(imp) _DEF_FSLNC_7(imp)
#endif

#if _CLG_SL9C_GAUGE
#define _DEF_F2SLNC_9(imp, b, ...) imp(9, b); EXPAND(_DEF_F2SLNC_8(imp, __VA_ARGS__))
#define _DEF_FSLNC_9(imp) imp(9); _DEF_FSLNC_8(imp)
#else
#define _DEF_F2SLNC_9(imp, b, ...) EXPAND(_DEF_F2SLNC_8(imp, __VA_ARGS__))
#define _DEF_FSLNC_9(imp) _DEF_FSLNC_8(imp)
#endif

#if _CLG_SL10C_GAUGE
#define _DEF_F2SLNC_10(imp, b, ...) imp(10, b); EXPAND(_DEF_F2SLNC_9(imp, __VA_ARGS__))
#define _DEF_FSLNC_10(imp) imp(10); _DEF_FSLNC_9(imp)
#else
#define _DEF_F2SLNC_10(imp, b, ...) EXPAND(_DEF_F2SLNC_9(imp, __VA_ARGS__))
#define _DEF_FSLNC_10(imp) _DEF_FSLNC_9(imp)
#endif

#if _CLG_SL11C_GAUGE
#define _DEF_F2SLNC_11(imp, b, ...) imp(11, b); EXPAND(_DEF_F2SLNC_10(imp, __VA_ARGS__))
#define _DEF_FSLNC_11(imp) imp(11); _DEF_FSLNC_10(imp)
#else
#define _DEF_F2SLNC_11(imp, b, ...) EXPAND(_DEF_F2SLNC_10(imp, __VA_ARGS__))
#define _DEF_FSLNC_11(imp) _DEF_FSLNC_10(imp)
#endif

#if _CLG_SL12C_GAUGE
#define _DEF_F2SLNC_12(imp, b, ...) imp(12, b); EXPAND(_DEF_F2SLNC_11(imp, __VA_ARGS__))
#define _DEF_FSLNC_12(imp) imp(12); _DEF_FSLNC_11(imp)
#else
#define _DEF_F2SLNC_12(imp, b, ...) EXPAND(_DEF_F2SLNC_11(imp, __VA_ARGS__))
#define _DEF_FSLNC_12(imp) _DEF_FSLNC_11(imp)
#endif

#if _CLG_SL13C_GAUGE
#define _DEF_F2SLNC_13(imp, b, ...) imp(13, b); EXPAND(_DEF_F2SLNC_12(imp, __VA_ARGS__))
#define _DEF_FSLNC_13(imp) imp(13); _DEF_FSLNC_12(imp)
#else
#define _DEF_F2SLNC_13(imp, b, ...) EXPAND(_DEF_F2SLNC_12(imp, __VA_ARGS__))
#define _DEF_FSLNC_13(imp) _DEF_FSLNC_12(imp)
#endif

#if _CLG_SL14C_GAUGE
#define _DEF_F2SLNC_14(imp, b, ...) imp(14, b); EXPAND(_DEF_F2SLNC_13(imp, __VA_ARGS__))
#define _DEF_FSLNC_14(imp) imp(14); _DEF_FSLNC_13(imp)
#else
#define _DEF_F2SLNC_14(imp, b, ...) EXPAND(_DEF_F2SLNC_13(imp, __VA_ARGS__))
#define _DEF_FSLNC_14(imp) _DEF_FSLNC_13(imp)
#endif

#if _CLG_SL15C_GAUGE
#define _DEF_F2SLNC_15(imp, b, ...) imp(15, b); EXPAND(_DEF_F2SLNC_14(imp, __VA_ARGS__))
#define _DEF_FSLNC_15(imp) imp(15); _DEF_FSLNC_14(imp)
#else
#define _DEF_F2SLNC_15(imp, b, ...) EXPAND(_DEF_F2SLNC_14(imp, __VA_ARGS__))
#define _DEF_FSLNC_15(imp) _DEF_FSLNC_14(imp)
#endif

#if _CLG_SL16C_GAUGE
#define _DEF_F2SLNC_16(imp, b, ...) imp(16, b); EXPAND(_DEF_F2SLNC_15(imp, __VA_ARGS__))
#define _DEF_FSLNC_16(imp) imp(16); _DEF_FSLNC_15(imp)
#else
#define _DEF_F2SLNC_16(imp, b, ...) EXPAND(_DEF_F2SLNC_15(imp, __VA_ARGS__))
#define _DEF_FSLNC_16(imp) _DEF_FSLNC_15(imp)
#endif

#define _DEF_F2SLNC_N_IMPL(n, args) _DEF_F2SLNC_##n args
#define _DEF_F2SLNC_N(n, ...)       _DEF_F2SLNC_N_IMPL(n, (__VA_ARGS__))
#define _DEF_FSLNC_N(n, imp) _DEF_FSLNC_##n(imp)

#if _CLG_U2_GAUGE
#define _DEF_F2UN_2(imp, b, ...) imp(2, b);
#define _DEF_FUN_2(imp) imp(2)
#else
#define _DEF_F2UN_2(imp, b, ...)
#define _DEF_FUN_2(imp)
#endif

#if _CLG_U3_GAUGE
#define _DEF_F2UN_3(imp, b, ...) imp(3, b); EXPAND(_DEF_F2UN_2(imp, __VA_ARGS__))
#define _DEF_FUN_3(imp) imp(3); _DEF_FUN_2(imp)
#else
#define _DEF_F2UN_3(imp, b, ...) EXPAND(_DEF_F2UN_2(imp, __VA_ARGS__))
#define _DEF_FUN_3(imp) _DEF_FUN_2(imp)
#endif

#if _CLG_U4_GAUGE
#define _DEF_F2UN_4(imp, b, ...) imp(4, b); EXPAND(_DEF_F2UN_3(imp, __VA_ARGS__))
#define _DEF_FUN_4(imp) imp(4); _DEF_FUN_3(imp)
#else
#define _DEF_F2UN_4(imp, b, ...) EXPAND(_DEF_F2UN_3(imp, __VA_ARGS__))
#define _DEF_FUN_4(imp) _DEF_FUN_3(imp)
#endif

#if _CLG_U5_GAUGE
#define _DEF_F2UN_5(imp, b, ...) imp(5, b); EXPAND(_DEF_F2UN_4(imp, __VA_ARGS__))
#define _DEF_FUN_5(imp) imp(5); _DEF_FUN_4(imp)
#else
#define _DEF_F2UN_5(imp, b, ...) EXPAND(_DEF_F2UN_4(imp, __VA_ARGS__))
#define _DEF_FUN_5(imp) _DEF_FUN_4(imp)
#endif

#if _CLG_U6_GAUGE
#define _DEF_F2UN_6(imp, b, ...) imp(6, b); EXPAND(_DEF_F2UN_5(imp, __VA_ARGS__))
#define _DEF_FUN_6(imp) imp(6); _DEF_FUN_5(imp)
#else
#define _DEF_F2UN_6(imp, b, ...) EXPAND(_DEF_F2UN_5(imp, __VA_ARGS__))
#define _DEF_FUN_6(imp) _DEF_FUN_5(imp)
#endif

#if _CLG_U7_GAUGE
#define _DEF_F2UN_7(imp, b, ...) imp(7, b); EXPAND(_DEF_F2UN_6(imp, __VA_ARGS__))
#define _DEF_FUN_7(imp) imp(7); _DEF_FUN_6(imp)
#else
#define _DEF_F2UN_7(imp, b, ...) EXPAND(_DEF_F2UN_6(imp, __VA_ARGS__))
#define _DEF_FUN_7(imp) _DEF_FUN_6(imp)
#endif

#if _CLG_U8_GAUGE
#define _DEF_F2UN_8(imp, b, ...) imp(8, b); EXPAND(_DEF_F2UN_7(imp, __VA_ARGS__))
#define _DEF_FUN_8(imp) imp(8); _DEF_FUN_7(imp)
#else
#define _DEF_F2UN_8(imp, b, ...) EXPAND(_DEF_F2UN_7(imp, __VA_ARGS__))
#define _DEF_FUN_8(imp) _DEF_FUN_7(imp)
#endif

#if _CLG_U9_GAUGE
#define _DEF_F2UN_9(imp, b, ...) imp(9, b); EXPAND(_DEF_F2UN_8(imp, __VA_ARGS__))
#define _DEF_FUN_9(imp) imp(9); _DEF_FUN_8(imp)
#else
#define _DEF_F2UN_9(imp, b, ...) EXPAND(_DEF_F2UN_8(imp, __VA_ARGS__))
#define _DEF_FUN_9(imp) _DEF_FUN_8(imp)
#endif

#if _CLG_U10_GAUGE
#define _DEF_F2UN_10(imp, b, ...) imp(10, b); EXPAND(_DEF_F2UN_9(imp, __VA_ARGS__))
#define _DEF_FUN_10(imp) imp(10); _DEF_FUN_9(imp)
#else
#define _DEF_F2UN_10(imp, b, ...) EXPAND(_DEF_F2UN_9(imp, __VA_ARGS__))
#define _DEF_FUN_10(imp) _DEF_FUN_9(imp)
#endif

#if _CLG_U11_GAUGE
#define _DEF_F2UN_11(imp, b, ...) imp(11, b); EXPAND(_DEF_F2UN_10(imp, __VA_ARGS__))
#define _DEF_FUN_11(imp) imp(11); _DEF_FUN_10(imp)
#else
#define _DEF_F2UN_11(imp, b, ...) EXPAND(_DEF_F2UN_10(imp, __VA_ARGS__))
#define _DEF_FUN_11(imp) _DEF_FUN_10(imp)
#endif

#if _CLG_U12_GAUGE
#define _DEF_F2UN_12(imp, b, ...) imp(12, b); EXPAND(_DEF_F2UN_11(imp, __VA_ARGS__))
#define _DEF_FUN_12(imp) imp(12); _DEF_FUN_11(imp)
#else
#define _DEF_F2UN_12(imp, b, ...) EXPAND(_DEF_F2UN_11(imp, __VA_ARGS__))
#define _DEF_FUN_12(imp) _DEF_FUN_11(imp)
#endif

#if _CLG_U13_GAUGE
#define _DEF_F2UN_13(imp, b, ...) imp(13, b); EXPAND(_DEF_F2UN_12(imp, __VA_ARGS__))
#define _DEF_FUN_13(imp) imp(13); _DEF_FUN_12(imp)
#else
#define _DEF_F2UN_13(imp, b, ...) EXPAND(_DEF_F2UN_12(imp, __VA_ARGS__))
#define _DEF_FUN_13(imp) _DEF_FUN_12(imp)
#endif

#if _CLG_U14_GAUGE
#define _DEF_F2UN_14(imp, b, ...) imp(14, b); EXPAND(_DEF_F2UN_13(imp, __VA_ARGS__))
#define _DEF_FUN_14(imp) imp(14); _DEF_FUN_13(imp)
#else
#define _DEF_F2UN_14(imp, b, ...) EXPAND(_DEF_F2UN_13(imp, __VA_ARGS__))
#define _DEF_FUN_14(imp) _DEF_FUN_13(imp)
#endif

#if _CLG_U15_GAUGE
#define _DEF_F2UN_15(imp, b, ...) imp(15, b); EXPAND(_DEF_F2UN_14(imp, __VA_ARGS__))
#define _DEF_FUN_15(imp) imp(15); _DEF_FUN_14(imp)
#else
#define _DEF_F2UN_15(imp, b, ...) EXPAND(_DEF_F2UN_14(imp, __VA_ARGS__))
#define _DEF_FUN_15(imp) _DEF_FUN_14(imp)
#endif

#if _CLG_U16_GAUGE
#define _DEF_F2UN_16(imp, b, ...) imp(16, b); EXPAND(_DEF_F2UN_15(imp, __VA_ARGS__))
#define _DEF_FUN_16(imp) imp(16); _DEF_FUN_15(imp)
#else
#define _DEF_F2UN_16(imp, b, ...) EXPAND(_DEF_F2UN_15(imp, __VA_ARGS__))
#define _DEF_FUN_16(imp) _DEF_FUN_15(imp)
#endif

#define _DEF_F2UN_N_IMPL(n, args) _DEF_F2UN_##n args
#define _DEF_F2UN_N(n, ...)       _DEF_F2UN_N_IMPL(n, (__VA_ARGS__))

#define _DEF_FUN_N(n, imp) _DEF_FUN_##n(imp)


#if _CLG_O2_GAUGE
#define _DEF_F2ON_2(imp, b, ...) imp(2, b);
#define _DEF_FON_2(imp) imp(2)
#else
#define _DEF_F2ON_2(imp, b, ...)
#define _DEF_FON_2(imp)
#endif

#if _CLG_O3_GAUGE
#define _DEF_F2ON_3(imp, b, ...) imp(3, b); EXPAND(_DEF_F2ON_2(imp, __VA_ARGS__))
#define _DEF_FON_3(imp) imp(3); _DEF_FON_2(imp)
#else
#define _DEF_F2ON_3(imp, b, ...) EXPAND(_DEF_F2ON_2(imp, __VA_ARGS__))
#define _DEF_FON_3(imp) _DEF_FON_2(imp)
#endif

#if _CLG_O4_GAUGE
#define _DEF_F2ON_4(imp, b, ...) imp(4, b); EXPAND(_DEF_F2ON_3(imp, __VA_ARGS__))
#define _DEF_FON_4(imp) imp(4); _DEF_FON_3(imp)
#else
#define _DEF_F2ON_4(imp, b, ...) EXPAND(_DEF_F2ON_3(imp, __VA_ARGS__))
#define _DEF_FON_4(imp) _DEF_FON_3(imp)
#endif

#if _CLG_O5_GAUGE
#define _DEF_F2ON_5(imp, b, ...) imp(5, b); EXPAND(_DEF_F2ON_4(imp, __VA_ARGS__))
#define _DEF_FON_5(imp) imp(5); _DEF_FON_4(imp)
#else
#define _DEF_F2ON_5(imp, b, ...) EXPAND(_DEF_F2ON_4(imp, __VA_ARGS__))
#define _DEF_FON_5(imp) _DEF_FON_4(imp)
#endif

#if _CLG_O6_GAUGE
#define _DEF_F2ON_6(imp, b, ...) imp(6, b); EXPAND(_DEF_F2ON_5(imp, __VA_ARGS__))
#define _DEF_FON_6(imp) imp(6); _DEF_FON_5(imp)
#else
#define _DEF_F2ON_6(imp, b, ...) EXPAND(_DEF_F2ON_5(imp, __VA_ARGS__))
#define _DEF_FON_6(imp) _DEF_FON_5(imp)
#endif

#if _CLG_O7_GAUGE
#define _DEF_F2ON_7(imp, b, ...) imp(7, b); EXPAND(_DEF_F2ON_6(imp, __VA_ARGS__))
#define _DEF_FON_7(imp) imp(7); _DEF_FON_6(imp)
#else
#define _DEF_F2ON_7(imp, b, ...) EXPAND(_DEF_F2ON_6(imp, __VA_ARGS__))
#define _DEF_FON_7(imp) _DEF_FON_6(imp)
#endif

#if _CLG_O8_GAUGE
#define _DEF_F2ON_8(imp, b, ...) imp(8, b); EXPAND(_DEF_F2ON_7(imp, __VA_ARGS__))
#define _DEF_FON_8(imp) imp(8); _DEF_FON_7(imp)
#else
#define _DEF_F2ON_8(imp, b, ...) EXPAND(_DEF_F2ON_7(imp, __VA_ARGS__))
#define _DEF_FON_8(imp) _DEF_FON_7(imp)
#endif

#if _CLG_O9_GAUGE
#define _DEF_F2ON_9(imp, b, ...) imp(9, b); EXPAND(_DEF_F2ON_8(imp, __VA_ARGS__))
#define _DEF_FON_9(imp) imp(9); _DEF_FON_8(imp)
#else
#define _DEF_F2ON_9(imp, b, ...) EXPAND(_DEF_F2ON_8(imp, __VA_ARGS__))
#define _DEF_FON_9(imp) _DEF_FON_8(imp)
#endif

#if _CLG_O10_GAUGE
#define _DEF_F2ON_10(imp, b, ...) imp(10, b); EXPAND(_DEF_F2ON_9(imp, __VA_ARGS__))
#define _DEF_FON_10(imp) imp(10); _DEF_FON_9(imp)
#else
#define _DEF_F2ON_10(imp, b, ...) EXPAND(_DEF_F2ON_9(imp, __VA_ARGS__))
#define _DEF_FON_10(imp) _DEF_FON_9(imp)
#endif

#if _CLG_O11_GAUGE
#define _DEF_F2ON_11(imp, b, ...) imp(11, b); EXPAND(_DEF_F2ON_10(imp, __VA_ARGS__))
#define _DEF_FON_11(imp) imp(11); _DEF_FON_10(imp)
#else
#define _DEF_F2ON_11(imp, b, ...) EXPAND(_DEF_F2ON_10(imp, __VA_ARGS__))
#define _DEF_FON_11(imp) _DEF_FON_10(imp)
#endif

#if _CLG_O12_GAUGE
#define _DEF_F2ON_12(imp, b, ...) imp(12, b); EXPAND(_DEF_F2ON_11(imp, __VA_ARGS__))
#define _DEF_FON_12(imp) imp(12); _DEF_FON_11(imp)
#else
#define _DEF_F2ON_12(imp, b, ...) EXPAND(_DEF_F2ON_11(imp, __VA_ARGS__))
#define _DEF_FON_12(imp) _DEF_FON_11(imp)
#endif

#if _CLG_O13_GAUGE
#define _DEF_F2ON_13(imp, b, ...) imp(13, b); EXPAND(_DEF_F2ON_12(imp, __VA_ARGS__))
#define _DEF_FON_13(imp) imp(13); _DEF_FON_12(imp)
#else
#define _DEF_F2ON_13(imp, b, ...) EXPAND(_DEF_F2ON_12(imp, __VA_ARGS__))
#define _DEF_FON_13(imp) _DEF_FON_12(imp)
#endif

#if _CLG_O14_GAUGE
#define _DEF_F2ON_14(imp, b, ...) imp(14, b); EXPAND(_DEF_F2ON_13(imp, __VA_ARGS__))
#define _DEF_FON_14(imp) imp(14); _DEF_FON_13(imp)
#else
#define _DEF_F2ON_14(imp, b, ...) EXPAND(_DEF_F2ON_13(imp, __VA_ARGS__))
#define _DEF_FON_14(imp) _DEF_FON_13(imp)
#endif

#if _CLG_O15_GAUGE
#define _DEF_F2ON_15(imp, b, ...) imp(15, b); EXPAND(_DEF_F2ON_14(imp, __VA_ARGS__))
#define _DEF_FON_15(imp) imp(15); _DEF_FON_14(imp)
#else
#define _DEF_F2ON_15(imp, b, ...) EXPAND(_DEF_F2ON_14(imp, __VA_ARGS__))
#define _DEF_FON_15(imp) _DEF_FON_14(imp)
#endif

#if _CLG_O16_GAUGE
#define _DEF_F2ON_16(imp, b, ...) imp(16, b); EXPAND(_DEF_F2ON_15(imp, __VA_ARGS__))
#define _DEF_FON_16(imp) imp(16); _DEF_FON_15(imp)
#else
#define _DEF_F2ON_16(imp, b, ...) EXPAND(_DEF_F2ON_15(imp, __VA_ARGS__))
#define _DEF_FON_16(imp) _DEF_FON_15(imp)
#endif

#define _DEF_F2ON_N_IMPL(n, args) _DEF_F2ON_##n args
#define _DEF_F2ON_N(n, ...)       _DEF_F2ON_N_IMPL(n, (__VA_ARGS__))

#define _DEF_FON_N(n, imp) _DEF_FON_##n(imp)


#if _CLG_SO2_GAUGE
#define _DEF_F2SON_2(imp, b, ...) imp(2, b);
#define _DEF_FSON_2(imp) imp(2)
#else
#define _DEF_F2SON_2(imp, b, ...)
#define _DEF_FSON_2(imp)
#endif

#if _CLG_SO3_GAUGE
#define _DEF_F2SON_3(imp, b, ...) imp(3, b); EXPAND(_DEF_F2SON_2(imp, __VA_ARGS__))
#define _DEF_FSON_3(imp) imp(3); _DEF_FSON_2(imp)
#else
#define _DEF_F2SON_3(imp, b, ...) EXPAND(_DEF_F2SON_2(imp, __VA_ARGS__))
#define _DEF_FSON_3(imp) _DEF_FSON_2(imp)
#endif

#if _CLG_SO4_GAUGE
#define _DEF_F2SON_4(imp, b, ...) imp(4, b); EXPAND(_DEF_F2SON_3(imp, __VA_ARGS__))
#define _DEF_FSON_4(imp) imp(4); _DEF_FSON_3(imp)
#else
#define _DEF_F2SON_4(imp, b, ...) EXPAND(_DEF_F2SON_3(imp, __VA_ARGS__))
#define _DEF_FSON_4(imp) _DEF_FSON_3(imp)
#endif

#if _CLG_SO5_GAUGE
#define _DEF_F2SON_5(imp, b, ...) imp(5, b); EXPAND(_DEF_F2SON_4(imp, __VA_ARGS__))
#define _DEF_FSON_5(imp) imp(5); _DEF_FSON_4(imp)
#else
#define _DEF_F2SON_5(imp, b, ...) EXPAND(_DEF_F2SON_4(imp, __VA_ARGS__))
#define _DEF_FSON_5(imp) _DEF_FSON_4(imp)
#endif

#if _CLG_SO6_GAUGE
#define _DEF_F2SON_6(imp, b, ...) imp(6, b); EXPAND(_DEF_F2SON_5(imp, __VA_ARGS__))
#define _DEF_FSON_6(imp) imp(6); _DEF_FSON_5(imp)
#else
#define _DEF_F2SON_6(imp, b, ...) EXPAND(_DEF_F2SON_5(imp, __VA_ARGS__))
#define _DEF_FSON_6(imp) _DEF_FSON_5(imp)
#endif

#if _CLG_SO7_GAUGE
#define _DEF_F2SON_7(imp, b, ...) imp(7, b); EXPAND(_DEF_F2SON_6(imp, __VA_ARGS__))
#define _DEF_FSON_7(imp) imp(7); _DEF_FSON_6(imp)
#else
#define _DEF_F2SON_7(imp, b, ...) EXPAND(_DEF_F2SON_6(imp, __VA_ARGS__))
#define _DEF_FSON_7(imp) _DEF_FSON_6(imp)
#endif

#if _CLG_SO8_GAUGE
#define _DEF_F2SON_8(imp, b, ...) imp(8, b); EXPAND(_DEF_F2SON_7(imp, __VA_ARGS__))
#define _DEF_FSON_8(imp) imp(8); _DEF_FSON_7(imp)
#else
#define _DEF_F2SON_8(imp, b, ...) EXPAND(_DEF_F2SON_7(imp, __VA_ARGS__))
#define _DEF_FSON_8(imp) _DEF_FSON_7(imp)
#endif

#if _CLG_SO9_GAUGE
#define _DEF_F2SON_9(imp, b, ...) imp(9, b); EXPAND(_DEF_F2SON_8(imp, __VA_ARGS__))
#define _DEF_FSON_9(imp) imp(9); _DEF_FSON_8(imp)
#else
#define _DEF_F2SON_9(imp, b, ...) EXPAND(_DEF_F2SON_8(imp, __VA_ARGS__))
#define _DEF_FSON_9(imp) _DEF_FSON_8(imp)
#endif

#if _CLG_SO10_GAUGE
#define _DEF_F2SON_10(imp, b, ...) imp(10, b); EXPAND(_DEF_F2SON_9(imp, __VA_ARGS__))
#define _DEF_FSON_10(imp) imp(10); _DEF_FSON_9(imp)
#else
#define _DEF_F2SON_10(imp, b, ...) EXPAND(_DEF_F2SON_9(imp, __VA_ARGS__))
#define _DEF_FSON_10(imp) _DEF_FSON_9(imp)
#endif

#if _CLG_SO11_GAUGE
#define _DEF_F2SON_11(imp, b, ...) imp(11, b); EXPAND(_DEF_F2SON_10(imp, __VA_ARGS__))
#define _DEF_FSON_11(imp) imp(11); _DEF_FSON_10(imp)
#else
#define _DEF_F2SON_11(imp, b, ...) EXPAND(_DEF_F2SON_10(imp, __VA_ARGS__))
#define _DEF_FSON_11(imp) _DEF_FSON_10(imp)
#endif

#if _CLG_SO12_GAUGE
#define _DEF_F2SON_12(imp, b, ...) imp(12, b); EXPAND(_DEF_F2SON_11(imp, __VA_ARGS__))
#define _DEF_FSON_12(imp) imp(12); _DEF_FSON_11(imp)
#else
#define _DEF_F2SON_12(imp, b, ...) EXPAND(_DEF_F2SON_11(imp, __VA_ARGS__))
#define _DEF_FSON_12(imp) _DEF_FSON_11(imp)
#endif

#if _CLG_SO13_GAUGE
#define _DEF_F2SON_13(imp, b, ...) imp(13, b); EXPAND(_DEF_F2SON_12(imp, __VA_ARGS__))
#define _DEF_FSON_13(imp) imp(13); _DEF_FSON_12(imp)
#else
#define _DEF_F2SON_13(imp, b, ...) EXPAND(_DEF_F2SON_12(imp, __VA_ARGS__))
#define _DEF_FSON_13(imp) _DEF_FSON_12(imp)
#endif

#if _CLG_SO14_GAUGE
#define _DEF_F2SON_14(imp, b, ...) imp(14, b); EXPAND(_DEF_F2SON_13(imp, __VA_ARGS__))
#define _DEF_FSON_14(imp) imp(14); _DEF_FSON_13(imp)
#else
#define _DEF_F2SON_14(imp, b, ...) EXPAND(_DEF_F2SON_13(imp, __VA_ARGS__))
#define _DEF_FSON_14(imp) _DEF_FSON_13(imp)
#endif

#if _CLG_SO15_GAUGE
#define _DEF_F2SON_15(imp, b, ...) imp(15, b); EXPAND(_DEF_F2SON_14(imp, __VA_ARGS__))
#define _DEF_FSON_15(imp) imp(15); _DEF_FSON_14(imp)
#else
#define _DEF_F2SON_15(imp, b, ...) EXPAND(_DEF_F2SON_14(imp, __VA_ARGS__))
#define _DEF_FSON_15(imp) _DEF_FSON_14(imp)
#endif

#if _CLG_SO16_GAUGE
#define _DEF_F2SON_16(imp, b, ...) imp(16, b); EXPAND(_DEF_F2SON_15(imp, __VA_ARGS__))
#define _DEF_FSON_16(imp) imp(16); _DEF_FSON_15(imp)
#else
#define _DEF_F2SON_16(imp, b, ...) EXPAND(_DEF_F2SON_15(imp, __VA_ARGS__))
#define _DEF_FSON_16(imp) _DEF_FSON_15(imp)
#endif

#define _DEF_F2SON_N_IMPL(n, args) _DEF_F2SON_##n args
#define _DEF_F2SON_N(n, ...)       _DEF_F2SON_N_IMPL(n, (__VA_ARGS__))

#define _DEF_FSON_N(n, imp) _DEF_FSON_##n(imp)


#endif //#ifndef _CLGDEFINE_H_

//=============================================================================
// END OF FILE
//=============================================================================