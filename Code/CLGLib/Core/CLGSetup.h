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
//It's better to complie using the maximum thread per block of the device of the computer.
#define _CLG_USE_LAUNCH_BOUND 1
#define _CLG_LAUNCH_MAX_THREAD 1024

//It is best to compile both 0 and 1 to see the difference.
//Sometimes, the compiler will try to optimize the instractions but just make things slower.
//Then, _CLG_USE_INTRINSC_FLOAT = 1 should be used to prevent compiler from 'optimize' it.
#define _CLG_USE_INTRINSC_FLOAT 0

#endif //#ifndef _CLGSETUP_H_

//=============================================================================
// END OF FILE
//=============================================================================