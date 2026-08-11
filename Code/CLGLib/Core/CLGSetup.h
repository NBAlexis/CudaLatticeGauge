//=============================================================================
// FILENAME : CLGSetup.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION:
//  [mm/dd/yy]
//  [12/2/2018 nbale]
//=============================================================================
#pragma once

#ifndef _CLGSETUP_H_
#define _CLGSETUP_H_

//No support for unicode, the unicode is not well supported in CUDA
//#define _CLG_UNICODE 1

#ifdef DEBUG
#define _CLG_DEBUG 1
#endif

//Note: Important!
//This is the tag for windows, msvc specific
//Ignore _MSC_VER, which is just for Visual Studio IDE specific, and should be harmless
#ifdef WIN64
#define _CLG_WIN 1
#elif defined(_WIN64)
#define _CLG_WIN 1
#endif

//Since cuda only support x64, this is reductent. 
//This is just for seperate those code depending on whether 32 or 64
#define _CLG_X64 1

//_CLG_DOUBLEFLOAT = 0 or 1. 
//Note that single float is rarely the problem for accuracy, but much much faster
#ifndef _CLG_DOUBLEFLOAT
#if _CLG_DEBUG
#define _CLG_DOUBLEFLOAT 0
#else
#define _CLG_DOUBLEFLOAT 0
#endif
#endif

//_CLG_USE_LAUNCH_BOUND = 0 or 1.
//NOTE: If the regcount required is out-numbered, sometimes, there is NO error message!
//So, either be sure to build with _CLG_USE_LAUNCH_BOUND = 1, or reduce the thread count
//reduce the thread count is expansive, so _CLG_USE_LAUNCH_BOUND = 1 is recommanded
//It's better to complie using the maximum thread per block of the device of the computer.
#if _CLG_DEBUG
#define _CLG_USE_LAUNCH_BOUND 0
#else
#define _CLG_USE_LAUNCH_BOUND 0
#endif

#ifndef _CLG_LAUNCH_MAX_THREAD
#if _CLG_USE_LAUNCH_BOUND
#if _CLG_DEBUG
    #define _CLG_LAUNCH_MAX_THREAD 512U
    #define _CLG_LAUNCH_MAX_THREADHALF 256U
#else
    #define _CLG_LAUNCH_MAX_THREAD 1024U
    #define _CLG_LAUNCH_MAX_THREADHALF 512U
#endif
#else
#if _CLG_DEBUG
    #define _CLG_LAUNCH_MAX_THREAD 256U
    #define _CLG_LAUNCH_MAX_THREADHALF 128U
#else
    #define _CLG_LAUNCH_MAX_THREAD 1024U
    #define _CLG_LAUNCH_MAX_THREADHALF 512U
#endif
#endif
#else
    #define _CLG_LAUNCH_MAX_THREADHALF (_CLG_LAUNCH_MAX_THREAD >> 1U)
#endif

#ifndef _CLG_PROFILER
#define _CLG_PROFILER 0
#endif

//Synchronize the device after each kernel launch (in launchKernel and _CHECKCUDA) to catch errors early
//NOTE: this serializes all kernel launches, use it for debug only
#ifndef _CLG_CHECKSYNCHRONIZE
#define _CLG_CHECKSYNCHRONIZE 0
#endif

#ifndef _CLG_ENABLE_RELEASE_ASSERT
#define _CLG_ENABLE_RELEASE_ASSERT 0
#endif

#ifndef _CLG_CHECKCUDAERRORS
#define _CLG_CHECKCUDAERRORS 1
#endif

#ifndef _CLG_ASSUME_SQUARE_LATTICE
#define _CLG_ASSUME_SQUARE_LATTICE 1
#endif

//_CLG_MULTI_GPU = 0 or 1.
//Multi-GPU support: one MPI process per GPU, with the lattice split into sub-lattices.
//When 0, everything degenerates to the original single-GPU code path with no behaviour change.
//Defined by the build system: VS uses the Debug_MG/Release_MG configurations,
//CMake uses -DCLG_MULTI_GPU=1. Do NOT include <mpi.h> outside of _CLG_MULTI_GPU guards,
//so that a machine without MPI installed can still build the single-GPU configurations.
#ifndef _CLG_MULTI_GPU
#define _CLG_MULTI_GPU 0
#endif

#define _CLG_LAUNCH_KERNEL 1
#define _CLG_PADDING 0

//Compile-time array dimensions shared by CudaHelper.h and CHaloManager.h.
//Defined here (before both are included) to avoid a circular include dependency:
//  CudaHelper.h needs appGetHaloManager() -> must include CHaloManager.h first,
//  CHaloManager.h needs kMaxFieldCount -> must see CudaHelper.h first.
//Moving the constants here breaks the cycle.
enum
{
    kContentLength  = 128,
    kMaxFieldCount  = 32,
};

//do not change this, unless the SUN, SUNVectors are changed acoordingly
#define _MAX_SUN 16

//not yet implemented
//less than _MAX_SUN
#define _CLG_SU2_GAUGE 0

#define _CLG_SU4_GAUGE 1
#define _CLG_SU4_GAUGER 1
#define _CLG_SU4_KS 0
#define _CLG_SU4_BOSON 1

#define _CLG_SU5_GAUGE 0
#define _CLG_SU5_GAUGER 0
#define _CLG_SU5_KS 0
#define _CLG_SU5_BOSON 0

#define _CLG_SU6_GAUGE 0
#define _CLG_SU6_GAUGER 0
#define _CLG_SU6_KS 0
#define _CLG_SU6_BOSON 0

#define _CLG_SU7_GAUGE 0
#define _CLG_SU7_GAUGER 0
#define _CLG_SU7_KS 0
#define _CLG_SU7_BOSON 0

#define _CLG_SU8_GAUGE 1
#define _CLG_SU8_GAUGER 0
#define _CLG_SU8_KS 0
#define _CLG_SU8_BOSON 1

#define _CLG_SU9_GAUGE 0
#define _CLG_SU10_GAUGE 0
#define _CLG_SU11_GAUGE 0
#define _CLG_SU12_GAUGE 0
#define _CLG_SU13_GAUGE 0
#define _CLG_SU14_GAUGE 0
#define _CLG_SU15_GAUGE 0
#define _CLG_SU16_GAUGE 0

// Z_N gauge group support
#define _MAX_ZN 6

#define _CLG_Z2_GAUGE 1
#define _CLG_Z3_GAUGE 1
#define _CLG_Z4_GAUGE 0
#define _CLG_Z5_GAUGE 0
#define _CLG_Z6_GAUGE 0

// D_N dihedral gauge group support (D3, D4, D8)
#define _MAX_DN 8

#define _CLG_D3_GAUGE 1
#define _CLG_D4_GAUGE 1
#define _CLG_D8_GAUGE 0


#define _MAX_SLNC 16
#define _CLG_SL2C_GAUGE 0
#define _CLG_SL3C_GAUGE 1
#define _CLG_SL4C_GAUGE 0
#define _CLG_SL5C_GAUGE 0
#define _CLG_SL6C_GAUGE 0
#define _CLG_SL7C_GAUGE 0
#define _CLG_SL8C_GAUGE 0
#define _CLG_SL9C_GAUGE 0
#define _CLG_SL10C_GAUGE 0
#define _CLG_SL11C_GAUGE 0
#define _CLG_SL12C_GAUGE 0
#define _CLG_SL13C_GAUGE 0
#define _CLG_SL14C_GAUGE 0
#define _CLG_SL15C_GAUGE 0
#define _CLG_SL16C_GAUGE 0

#define _MAX_UN 16
#define _CLG_U2_GAUGE 0
#define _CLG_U3_GAUGE 1
#define _CLG_U4_GAUGE 0
#define _CLG_U5_GAUGE 0
#define _CLG_U6_GAUGE 0
#define _CLG_U7_GAUGE 0
#define _CLG_U8_GAUGE 0
#define _CLG_U9_GAUGE 0
#define _CLG_U10_GAUGE 0
#define _CLG_U11_GAUGE 0
#define _CLG_U12_GAUGE 0
#define _CLG_U13_GAUGE 0
#define _CLG_U14_GAUGE 0
#define _CLG_U15_GAUGE 0
#define _CLG_U16_GAUGE 0

#define _MAX_ON 16
#define _CLG_O2_GAUGE 0
#define _CLG_O3_GAUGE 1
#define _CLG_O4_GAUGE 0
#define _CLG_O5_GAUGE 0
#define _CLG_O6_GAUGE 0
#define _CLG_O7_GAUGE 0
#define _CLG_O8_GAUGE 0
#define _CLG_O9_GAUGE 0
#define _CLG_O10_GAUGE 0
#define _CLG_O11_GAUGE 0
#define _CLG_O12_GAUGE 0
#define _CLG_O13_GAUGE 0
#define _CLG_O14_GAUGE 0
#define _CLG_O15_GAUGE 0
#define _CLG_O16_GAUGE 0

#define _MAX_SON 16
#define _CLG_SO2_GAUGE 0
#define _CLG_SO3_GAUGE 1
#define _CLG_SO4_GAUGE 0
#define _CLG_SO5_GAUGE 0
#define _CLG_SO6_GAUGE 0
#define _CLG_SO7_GAUGE 0
#define _CLG_SO8_GAUGE 0
#define _CLG_SO9_GAUGE 0
#define _CLG_SO10_GAUGE 0
#define _CLG_SO11_GAUGE 0
#define _CLG_SO12_GAUGE 0
#define _CLG_SO13_GAUGE 0
#define _CLG_SO14_GAUGE 0
#define _CLG_SO15_GAUGE 0
#define _CLG_SO16_GAUGE 0

#endif //#ifndef _CLGSETUP_H_

//=============================================================================
// END OF FILE
//=============================================================================