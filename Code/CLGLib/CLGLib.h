//=============================================================================
// FILENAME : CLGLib.h
// 
// DESCRIPTION:
// This is the one header file for all
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#ifndef _CLGLIB_H_
#define _CLGLIB_H_

#define _CLG_PUBLIC 1

#include "Core/CLGSetup.h"
#include "Core/CLGDefine.h"

#if defined(_CLG_WIN)
#   if !defined(CLGAPI)
#	    define __LIB_TITLE__	"CLGLib"
#       ifdef _CLG_PRIVATE
#           define CLGAPI __DLL_EXPORT
#       else
#   	    define CLGAPI __DLL_IMPORT
#       endif
#       ifndef _CLG_PRIVATE
#   	    ifdef _CLG_DEBUG
#	    	    define __LIB_FILE__	__LIB_TITLE__ "_d.lib"
#	        else
#		        define __LIB_FILE__ __LIB_TITLE__ ".lib"
#	        endif
#	        pragma __IMPORT_LIB(__LIB_FILE__)
#	        pragma message("linking with " __LIB_FILE__ "...")
#	        undef __LIB_FILE__
#	        undef __LIB_TITLE__
#       endif
#   endif
#else
#	define CLGAPI  
#endif

#if _CLG_WIN
#include "Platform/PlatformIncs.h"
#include "Platform/WinDefine.h"
#include "Platform/WinFunction.h"
#endif

#include "Tools/Tracer.h"
#include "Tools/Timer.h"
#include "Tools/CYAMLParser.h"

#include "Core/CudaHelperFunctions.h"
#include "Core/CudaHelper.h"

#include "Tools/Math/SU3.h"

#include "Data/CCommonData.h"
#include "Data/Boundary/CBoundaryCondition.h"
#include "Data/Boundary/CBoundaryConditionTorusSquare.h"
#include "Data/Lattice/CIndex.h"
#include "Data/Lattice/CIndexSquare.h"
#include "Data/Lattice/CLatticeData.h"
#include "Data/Field/CField.h"

#include "Data/Field/CFieldGauge.h"
#include "Data/Field/CFieldGaugeSU3.h"
#include "Data/Field/CFieldGaugeZ2.h"

#include "Data/Field/CFieldBoson.h"

#include "Data/Field/CFieldSpin.h"

#include "Data/Field/CFieldFermion.h"

#include "Core/CLGLibManager.h"

#ifndef _CLG_PRIVATE
__USE_NAMESPACE
#endif

#endif //#ifndef _CLGLIB_H_

//=============================================================================
// END OF FILE
//=============================================================================
