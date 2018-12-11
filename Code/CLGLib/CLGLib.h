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

#include "Tools/Data/CLinkedList.h"
#include "Tools/Data/TemplateFunctions.h"
#include "Tools/Data/MemStack.h"
#include "Tools/Data/CBitFlag.h"
//using these class to avoid warnings by STL...
#include "Tools/Data/TArray.h"
#include "Tools/Data/CCString.h"
#include "Tools/Data/THashMap.h"

#include "Tools/Tracer.h"
#include "Tools/Timer.h"
#include "Tools/CYAMLParser.h"

#include "Core/CudaHelperFunctions.h"
#include "Core/CudaHelper.h"

//====================================
//define some common function before decompose threads
#define preparethread \
CLatticeData* pLattice = CLatticeData::GetInstance(); \
dim3 block(pLattice->m_uiLatticeDecompose[0], pLattice->m_uiLatticeDecompose[1], pLattice->m_uiLatticeDecompose[2]); \
dim3 threads(pLattice->m_uiLatticeDecompose[3], pLattice->m_uiLatticeDecompose[4], pLattice->m_uiLatticeDecompose[5]);


#define preparethreadsimple \
dim3 block(pLattice->m_uiLatticeDecompose[0], pLattice->m_uiLatticeDecompose[1], pLattice->m_uiLatticeDecompose[2]); \
dim3 threads(pLattice->m_uiLatticeDecompose[3], pLattice->m_uiLatticeDecompose[4], pLattice->m_uiLatticeDecompose[5]);


#define intokernal \
UINT coord[4]; \
coord[0] = threadIdx.x + blockIdx.x * blockDim.x; \
coord[1] = threadIdx.y + blockIdx.y * blockDim.y; \
coord[2] = threadIdx.z + blockIdx.z * blockDim.z; \
UINT uiTLength = pLattice->m_uiTLength; \
UINT uiDir = pLattice->m_uiDir;


#include "Tools/Math/CudaComplexFunction.h"
#include "Tools/Math/Random.h"
#include "Tools/Math/SU3.h"
#include "Tools/Math/cuComplexI.h"
#include "Tools/Math/GammaMatrix.h"

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

#include "Data/Action/CAction.h"
#include "Data/Action/CActionGaugePlaquette.h"

#include "Update/CUpdator.h"
#include "Update/Continous/CIntegrator.h"
#include "Update/Continous/CIntegratorLeapFrog.h"
#include "Update/Continous/CIntegratorOmelyan.h"
#include "Update/Continous/CHMC.h"

#include "Core/CLGLibManager.h"

#ifndef _CLG_PRIVATE
__USE_NAMESPACE
#endif

#endif //#ifndef _CLGLIB_H_

//=============================================================================
// END OF FILE
//=============================================================================
