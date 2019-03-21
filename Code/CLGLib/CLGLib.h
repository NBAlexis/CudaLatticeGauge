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


#include "Platform/PlatformIncs.h"
#include "Platform/PlatformDefine.h"
#include "Tools/Data/STDStringFunctions.h"

#include "Core/CudaHelperFunctions.h"
#include "Core/CCudaBuffer.h"
#include "Core/CLGFloat.h"

#include "Tools/Data/CLinkedList.h"
#include "Tools/Data/TemplateFunctions.h"
#include "Tools/Data/MemStack.h"
#include "Tools/Data/CBitFlag.h"
//using these class to avoid warnings by STL...
#include "Tools/Data/TArray.h"
#include "Tools/Data/CCString.h"
#include "Tools/Data/THashMap.h"
#include "Tools/EnumGather.h"

#include "Platform/CFile.h"

#include "Tools/Tracer.h"
#include "Tools/Timer.h"
#include "Tools/CYAMLParser.h"

#include "Core/CBase.h"
#include "Core/CudaHelper.h"


//====================================
//define some common function before decompose threads
#define preparethread \
dim3 block(_HC_DecompX, _HC_DecompY, _HC_DecompZ); \
dim3 threads(_HC_DecompLx, _HC_DecompLy, _HC_DecompLz);

#define preparethreadE(element_count) \
dim3 block(_HC_DecompX * element_count, _HC_DecompY, _HC_DecompZ); \
dim3 threads(_HC_DecompLx, _HC_DecompLy, _HC_DecompLz);


#define intokernal \
UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z)); 


#define intokernalE(element_count)\
UINT blockIdxX = blockIdx.x / element_count;\
UINT elementIdx = blockIdx.x % element_count; \
UINT uiSiteIndex = ((threadIdx.x + blockIdxX * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z));

#define intokernaldir \
UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z)); \
BYTE uiDir = static_cast<BYTE>(_DC_Dir);


#include "Tools/Math/CudaComplexFunction.h"
#include "Tools/Math/Random.h"
#include "Tools/Math/Vectors.h" //vectors.h must ealier than gamma matrix
#include "Tools/Math/SU3.h"
#include "Tools/Math/GammaMatrix.h" //gamma matrix must later than cuComplexI

#include "Tools/Math/CLinearAlgebraHelper.h"

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

//#include "Data/Field/CFieldBoson.h"

//#include "Data/Field/CFieldSpin.h"

#include "Data/Field/CFieldFermion.h"
#include "Data/Field/CFieldFermionWilsonSquareSU3.h"

#include "Data/Action/CAction.h"
#include "Data/Action/CActionGaugePlaquette.h"
#include "Data/Action/CActionFermionWilsonNf2.h"

#include "SparseLinearAlgebra/CSLASolver.h"
#include "SparseLinearAlgebra/CSolverBiCGstab.h"
#include "SparseLinearAlgebra/CSolverGMRES.h"
#include "SparseLinearAlgebra/CSolverGCR.h"
#include "SparseLinearAlgebra/CSolverGCRODR.h"

#include "Measurement/CMeasure.h"
#include "Measurement/CMeasurePlaqutteEnergy.h"
#include "Measurement/CMeasureMesonCorrelator.h"
#include "Measurement/CMeasurementManager.h"
#include "Measurement/GaugeSmearing/CGaugeSmearing.h"
#include "Measurement/GaugeSmearing/CGaugeSmearingAPEProj.h"
#include "Measurement/GaugeSmearing/CGaugeSmearingAPEStout.h"

#include "Update/CUpdator.h"
#include "Update/Continous/CIntegrator.h"
#include "Update/Continous/CIntegratorLeapFrog.h"
#include "Update/Continous/CIntegratorOmelyan.h"
#include "Update/Continous/CIntegratorForceGradient.h"
#include "Update/Continous/CIntegratorNestedLeapFrog.h"
#include "Update/Continous/CIntegratorNestedOmelyan.h"
#include "Update/Continous/CIntegratorNestedForceGradient.h"
#include "Update/Continous/CHMC.h"

#include "Core/CLGLibManager.h"

#ifndef _CLG_PRIVATE
__USE_NAMESPACE
#endif

#endif //#ifndef _CLGLIB_H_

//=============================================================================
// END OF FILE
//=============================================================================
