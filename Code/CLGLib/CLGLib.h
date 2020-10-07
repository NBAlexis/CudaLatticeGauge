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
#        define __LIB_TITLE__    "CLGLib"
#       ifdef _CLG_PRIVATE
#           define CLGAPI __DLL_EXPORT
#       else
#           define CLGAPI __DLL_IMPORT
#       endif
#       ifndef _CLG_PRIVATE
#           ifdef _CLG_DEBUG
#                define __LIB_FILE__    __LIB_TITLE__ "_d.lib"
#            else
#                define __LIB_FILE__ __LIB_TITLE__ ".lib"
#            endif
#            pragma __IMPORT_LIB(__LIB_FILE__)
#            pragma message("linking with " __LIB_FILE__ "...")
#            undef __LIB_FILE__
#            undef __LIB_TITLE__
#       endif
#   endif
#else
#    define CLGAPI  
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
const dim3 block(_HC_DecompX, _HC_DecompY, _HC_DecompZ); \
const dim3 threads(_HC_DecompLx, _HC_DecompLy, _HC_DecompLz);

#define preparethread_even \
UINT uiDecompX = ((_HC_DecompX & 1) == 0) ? (_HC_DecompX >> 1) : _HC_DecompX; \
UINT uiDecompLX = ((_HC_DecompX & 1) == 0) ? _HC_DecompLx : (_HC_DecompLx >> 1); \
const dim3 block(uiDecompX, _HC_DecompY, _HC_DecompZ); \
const dim3 threads(uiDecompLX, _HC_DecompLy, _HC_DecompLz);

#define preparethreadE(element_count) \
const dim3 block(_HC_DecompX * element_count, _HC_DecompY, _HC_DecompZ); \
const dim3 threads(_HC_DecompLx, _HC_DecompLy, _HC_DecompLz);


#define intokernal \
const UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z)); 


/**
 * Even to Site Index transform
 * The uiSiteIndex of even sites maybe not even, so need this transform
 * When run with evenSites, for convinient to decompose threads, we assume Nx * Ny is even
 * so evenIndex is 0 -> (#sites / 2) - 1
 * If sSite4 of (2 * evenIndex) is even, the uiSiteIndex is 2 * evenIndex
 * If sSite4 of (2 * evenIndex) is odd, the uiSiteIndex is 2 * evenIndex - 1
 *
 * For preparethread_even, threadIdx.x + blockIdx.x * blockDim.x = 0 -> Nx * Ny / 2,
 * so the total number is correctly 0 -> (#sites / 2) - 1
 */
#define intokernal_even \
UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z)) << 1; \
SSmallInt4 _sSite4 = __deviceSiteIndexToInt4(uiSiteIndex); \
if (_sSite4.IsOdd()) \
{ \
    uiSiteIndex = uiSiteIndex + 1; \
} 

#define intokernalInt4_even \
UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z)) << 1; \
SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex); \
if (sSite4.IsOdd()) \
{ \
    uiSiteIndex = uiSiteIndex + 1; \
    sSite4 = __deviceSiteIndexToInt4(uiSiteIndex); \
} 


#define intokernal_odd \
UINT uiSiteIndex = (((threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z)) << 1) + 1; \
SSmallInt4 _sSite4 = __deviceSiteIndexToInt4(uiSiteIndex); \
if (!_sSite4.IsOdd()) \
{ \
    uiSiteIndex = uiSiteIndex - 1; \
} 

#define intokernalInt4_odd \
UINT uiSiteIndex = (((threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z)) << 1) + 1; \
SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex); \
if (!sSite4.IsOdd()) \
{ \
    uiSiteIndex = uiSiteIndex - 1; \
    sSite4 = __deviceSiteIndexToInt4(uiSiteIndex); \
}

#define intokernalInt4 \
SSmallInt4 sSite4; \
const UINT _ixy = (threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.x = static_cast<SBYTE> (_ixy / _DC_Ly); \
sSite4.y = static_cast<SBYTE> (_ixy % _DC_Ly); \
sSite4.z = static_cast<SBYTE>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.w = static_cast<SBYTE>(threadIdx.z + blockIdx.z * blockDim.z); \
const UINT uiSiteIndex = _ixy * _DC_GridDimZT + sSite4.z * _DC_Lt + sSite4.w; 


#define intokernalOnlyInt4 \
SSmallInt4 sSite4; \
UINT _ixy = (threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.x = static_cast<SBYTE> (_ixy / _DC_Ly); \
sSite4.y = static_cast<SBYTE> (_ixy % _DC_Ly); \
sSite4.z = static_cast<SBYTE>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.w = static_cast<SBYTE>(threadIdx.z + blockIdx.z * blockDim.z); 

#define _QUICK_AXPY_BLOCK 2

#define intokernalE(element_count)\
UINT blockIdxX = blockIdx.x / element_count;\
UINT elementIdx = blockIdx.x % element_count; \
const UINT uiSiteIndex = ((threadIdx.x + blockIdxX * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z));

#define intokernaldir \
const UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z)); \
const BYTE uiDir = static_cast<BYTE>(_DC_Dir);


#define preparethread_S \
const dim3 block(m_pHDecomp[0], m_pHDecomp[1], m_pHDecomp[2]); \
const dim3 threads(m_pHDecomp[3], m_pHDecomp[4], m_pHDecomp[5]);

#define intokernalInt4_S \
SSmallInt4 sSite4; \
sSite4.x = static_cast<SBYTE>(threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.y = static_cast<SBYTE>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.z = static_cast<SBYTE>(threadIdx.z + blockIdx.z * blockDim.z); \
sSite4.w = uiT; \
const UINT uiSiteIndex = sSite4.x * _DC_MultX + sSite4.y * _DC_MultY + sSite4.z * _DC_Lt + sSite4.w; \
const UINT uiSiteIndex3D = (sSite4.x * _DC_Ly + sSite4.y) * _DC_Lz + sSite4.z;

#define intokernalInt4_S_Only3D \
SSmallInt4 sSite4; \
sSite4.x = static_cast<SBYTE>(threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.y = static_cast<SBYTE>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.z = static_cast<SBYTE>(threadIdx.z + blockIdx.z * blockDim.z); \
sSite4.w = uiT; \
const UINT uiSiteIndex3D = (sSite4.x * _DC_Ly + sSite4.y) * _DC_Lz + sSite4.z;

#define cudaSafeFree(ptr) if (NULL != ptr) { checkCudaErrors(cudaFree(ptr)); ptr = NULL; }

#include "Tools/Math/CudaComplexFunction.h"
#include "Tools/Math/Random.h"
#include "Tools/Math/Vectors.h" //vectors.h must ealier than gamma matrix
#include "Tools/Math/SU3.h"
#include "Tools/Math/GammaMatrix.h" //gamma matrix must later than cuComplexI

#include "Tools/Math/CLinearAlgebraHelper.h"
#include "Tools/Math/CLGFFT.h"
#include "Tools/Math/CRationalApproximation.h"

#include "Data/CCommonData.h"
#include "Data/Boundary/CBoundaryCondition.h"
#include "Data/Boundary/CBoundaryConditionTorusSquare.h"
#include "Data/Boundary/CBoundaryConditionPeriodicAndDirichletSquare.h"
#include "Data/Boundary/CBoundaryConditionProjectivePlaneSquare.h"
#include "Data/Lattice/CIndex.h"
#include "Data/Lattice/CIndexSquare.h"
#include "Data/Lattice/CIndexData.h"
#include "Data/Lattice/CLatticeData.h"

//=======================================================
//Field
#include "Data/Field/CField.h"
#include "Data/Field/BoundaryField/CFieldBoundary.h"

#include "Data/Field/CFieldGauge.h"
#include "Data/Field/BoundaryField/CFieldBoundaryGaugeSU3.h"
#include "Data/Field/CFieldGaugeSU3.h"
#include "Data/Field/CFieldGaugeSU3D.h"
#include "Data/Field/CFieldGaugeZ2.h"
#include "Data/Field/CFieldGaugeSU3TreeImproved.h"

//#include "Data/Field/CFieldBoson.h"

//#include "Data/Field/CFieldSpin.h"

#include "Data/Field/CFieldFermion.h"
#include "Data/Field/BoundaryField/CFieldBoundaryWilsonSquareSU3.h"
#include "Data/Field/CFieldFermionWilsonSquareSU3.h"
#include "Data/Field/CFieldFermionWilsonSquareSU3D.h"
#include "Data/Field/CFieldFermionWilsonSquareSU3DR.h"
#include "Data/Field/CFieldFermionWilsonSquareSU3DEven.h"
#include "Data/Field/CFieldFermionWilsonSquareSU3DREven.h"
#include "Data/Field/CFieldFermionWilsonSquareSU3Acc.h"
#include "Data/Field/CFieldFermionWilsonSquareSU3DRigidAcc.h"
#include "Data/Field/CFieldFermionWilsonSquareSU3Boost.h"
#include "Data/Field/CFieldFermionKSSU3.h"
#include "Data/Field/CFieldFermionKSSU3R.h"
#include "Data/Field/CFieldFermionKSSU3D.h"
#include "Data/Field/CFieldFermionKSSU3DR.h"
#include "Data/Field/CFieldFermionKSSU3EM.h"

//=====================================================

#include "Data/Action/CAction.h"
#include "Data/Action/CActionGaugePlaquette.h"
#include "Data/Action/CActionFermionWilsonNf2.h"
#include "Data/Action/CActionGaugePlaquetteRotating.h"
#include "Data/Action/CActionFermionKS.h"
#include "Data/Action/CActionGaugePlaquetteAcceleration.h"
#include "Data/Action/CActionGaugePlaquetteBoost.h"
#include "Data/Action/CActionGaugePlaquetteRigidAcc.h"

#include "SparseLinearAlgebra/CSLASolver.h"
#include "SparseLinearAlgebra/CSolverBiCGstab.h"
#include "SparseLinearAlgebra/CSolverGMRES.h"
#include "SparseLinearAlgebra/CSolverGCR.h"
#include "SparseLinearAlgebra/CSolverGCRODR.h"
#include "SparseLinearAlgebra/CSolverGMRESMDR.h"
#include "SparseLinearAlgebra/CSolverTFQMR.h"
#include "SparseLinearAlgebra/CMultiShiftSolver.h"
#include "SparseLinearAlgebra/CMultiShiftGMRES.h"
#include "SparseLinearAlgebra/CMultiShiftFOM.h"
#include "SparseLinearAlgebra/CMultiShiftBiCGStab.h"
#include "SparseLinearAlgebra/CMultiShiftNested.h"

#include "Measurement/CMeasure.h"
#include "Measurement/CMeasurePlaqutteEnergy.h"
#include "Measurement/CMeasureAction.h"
#include "Measurement/CMeasureMesonCorrelator.h"
#include "Measurement/CMeasureAMomentumJG.h"
#include "Measurement/CMeasureAMomentumJF.h"
#include "Measurement/CMeasureChargeAndCurrents.h"
#include "Measurement/CMeasureTopologicChargeXY.h"
#include "Measurement/CMeasurePolyakovXY.h"
#include "Measurement/CMeasurePolyakov.h"
#include "Measurement/CMeasureWilsonLoop.h"
#include "Measurement/CMeasureChiralCondensate.h"
#include "Measurement/CMeasureAMomentumStochastic.h"
#include "Measurement/CMeasureMesonCorrelatorStaggered.h"
#include "Measurement/CMeasureMesonCorrelatorStaggeredSimple.h"
#include "Measurement/CMeasureChiralCondensateKS.h"
#include "Measurement/CMeasureConnectedChiralSusceptibilityKS.h"

#include "Measurement/CMeasurementManager.h"
#include "Measurement/GaugeSmearing/CGaugeSmearing.h"
#include "Measurement/GaugeSmearing/CGaugeSmearingAPEProj.h"
#include "Measurement/GaugeSmearing/CGaugeSmearingAPEStout.h"

#include "GaugeFixing/CGaugeFixing.h"
#include "GaugeFixing/CGaugeFixingLandauCornell.h"
#include "GaugeFixing/CGaugeFixingCoulombCornell.h"
#include "GaugeFixing/CGaugeFixingLandauLosAlamos.h"
#include "GaugeFixing/CGaugeFixingCoulombLosAlamos.h"
#include "GaugeFixing/CGaugeFixingRandom.h"

#include "Update/CUpdator.h"
#include "Update/Continous/CIntegrator.h"
#include "Update/Continous/CIntegratorLeapFrog.h"
#include "Update/Continous/CIntegratorOmelyan.h"
#include "Update/Continous/CIntegratorForceGradient.h"
#include "Update/Continous/CIntegratorNestedLeapFrog.h"
#include "Update/Continous/CIntegratorNestedOmelyan.h"
#include "Update/Continous/CIntegratorNested11Stage.h"
#include "Update/Continous/CIntegratorNestedForceGradient.h"
#include "Update/Continous/CIntegratorMultiLevelNestedOmelyan.h"
#include "Update/Continous/CIntegratorMultiLevelNestedForceGradient.h"
#include "Update/Continous/CHMC.h"

#include "Core/CLGLibManager.h"

#ifndef _CLG_PRIVATE
__USE_NAMESPACE
#endif

#endif //#ifndef _CLGLIB_H_

//=============================================================================
// END OF FILE
//=============================================================================
