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

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stack>
#include <cstdarg>

#include <limits.h>
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

#include <algorithm> //c++14
#include <atomic> //replace interlock
#include <chrono> //for timer

#if _CLG_UNICODE

using ISTREAM = std::wistream;
using OSTREAM = std::wostream;
using OFSTREAM = std::wofstream;
using IFSTREAM = std::wifstream;
using ISTRINGSTREAM = std::wistringstream;
using STDSTRING = std::wstring;
#define COUT std::wcout
#define CIN std::wcin

#else

using ISTREAM = std::istream;
using OSTREAM = std::ostream;
using OFSTREAM = std::ofstream;
using IFSTREAM = std::ifstream;
using ISTRINGSTREAM = std::istringstream;
using STDSTRING = std::string;
#define COUT std::cout
#define CIN std::cin

#endif



#endif //#ifndef _PLATFORMINCS_H_

//=============================================================================
// END OF FILE
//=============================================================================