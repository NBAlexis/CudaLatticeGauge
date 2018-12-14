//=============================================================================
// FILENAME : PlatfomrIncs.h
// 
// DESCRIPTION:
// This is the file for all system include files
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#ifndef _PLATFORMINCS_H_
#define _PLATFORMINCS_H_

#ifdef _CLG_WIN

#include <iostream>
#include <fstream>
#include <string>
#include <stack>
#include <cstdarg>

#include <limits.h>
#include <windows.h>
#include <tchar.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <sstream>
#include <assert.h>
#include <time.h>

#ifdef _CLG_UNICODE

using ISTREAM = std::wistream;
using OSTREAM = std::wostream;
using OFSTREAM = std::wofstream;
using IFSTREAM = std::wifstream;
using ISTRINGSTREAM = std::wistringstream;
#define COUT std::wcout
#define CIN std::wcin

#else

using ISTREAM = std::istream;
using OSTREAM = std::ostream;
using OFSTREAM = std::ofstream;
using IFSTREAM = std::ifstream;
using ISTRINGSTREAM = std::istringstream;
#define COUT std::cout
#define CIN std::cin

#endif

#endif

#endif //#ifndef _PLATFORMINCS_H_

//=============================================================================
// END OF FILE
//=============================================================================