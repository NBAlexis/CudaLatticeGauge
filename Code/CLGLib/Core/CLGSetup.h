//=============================================================================
// FILENAME : CLGSetup.h
// 
// DESCRIPTION:
// This is the file for building options
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

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
#define _CLG_WIN 1

//_CLG_DOUBLEFLOAT = 0 or 1. 
//Note that single float is rarely the problem for accuracy, but much much faster
#define _CLG_DOUBLEFLOAT 0

//_CLG_USE_LAUNCH_BOUND = 0 or 1.
//NOTE: For release version, the if the regcount required is out-numbered, there is NO error message!
//So, either be sure to build with _CLG_USE_LAUNCH_BOUND = 1, or reduce the thread count
//reduce the thread count is expansive, so _CLG_USE_LAUNCH_BOUND = 1 is recommanded
//It's better to complie using the maximum thread per block of the device of the computer.
#if _CLG_DEBUG
#define _CLG_USE_LAUNCH_BOUND 1
#else
#define _CLG_USE_LAUNCH_BOUND 1
#endif
#define _CLG_LAUNCH_MAX_THREAD 1024


#endif //#ifndef _CLGSETUP_H_

//=============================================================================
// END OF FILE
//=============================================================================