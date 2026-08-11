//=============================================================================
// FILENAME : CLGLib_Private.h
// 
// DESCRIPTION:
// This is the header file for pre-compile header
//
// REVISION:
//  [mm/dd/yy]
//  [12/2/2018 nbale]
//=============================================================================

#ifndef _CLGLIB_PRIVATE_H_
#define _CLGLIB_PRIVATE_H_

#include "Core/CLGSetup.h"
#include "Core/CLGDefine.h"

#if !defined(CLGAPI)
#    define CLGAPI __DLL_EXPORT
#endif


#include "Platform/PlatformIncs.h"
#include "Platform/PlatformDefine.h"
#include "Tools/Data/STDStringFunctions.h"

#include "Core/CudaHelperFunctions.h"
#include "Core/CLGFloat.h"

#include "Tools/Data/CLinkedList.h"
#include "Tools/Data/TemplateFunctions.h"
#include "Tools/Data/MemStack.h"
#include "Tools/Data/CBitFlag.h"
//using these class to avoid warnings by STL...
#include "Tools/Data/TArray.h"
#include "Tools/Data/CCString.h"
#include "Tools/Data/CUUID.h"
#include "Tools/Data/CLGMD5.h"
#include "Tools/Data/THashMap.h"
#include "Tools/EnumGather.h"

#include "Platform/CFile.h"

#include "Tools/Tracer.h"
#include "Tools/Timer.h"
#include "Tools/CYAMLParser.h"
#include "Tools/CSimpleProfiler.h"
#include "Core/CCudaBuffer.h"

#include "Core/CBase.h"
//Multi-GPU: must precede CudaHelper.h, whose launch guard helpers use appGetHaloManager().
#include "Core/Distributed/CLGComm.h"
#include "Core/Distributed/CHaloManager.h"
#include "Core/Distributed/CLGHaloLayout.h"
#include "Core/CudaHelper.h"

#include "Tools/Math/CudaComplexFunction.h"
#include "Tools/Math/Random.h"
#include "Tools/Math/Vectors.h" //vectors.h must ealier than gamma matrix
#include "Tools/Math/VectorsN.h"
#include "Tools/Math/VectorsN2.h"
#include "Tools/Math/SU3.h"
#include "Tools/Math/SU2.h"
#include "Tools/Math/SUN.h"
#include "Tools/Math/ZN.h"
#include "Tools/Math/DN.h"
#include "Tools/Math/SLNC.h"
#include "Tools/Math/UN.h"
#include "Tools/Math/ON.h"
#include "Tools/Math/SON.h"
#include "Tools/Math/GammaMatrix.h" //gamma matrix must later than cuComplexI

#include "Tools/Math/CLinearAlgebraHelper.h"
#include "Tools/Math/CLGFFT.h"
#include "Tools/Math/CRationalApproximation.h"

#include "Data/CCommonData.h"
//Multi-GPU halo runtime wrappers: need the _HC_* const-integer macros above.
#include "Core/Distributed/CLGHaloLayoutRuntime.h"
#include "Data/Boundary/CBoundaryCondition.h"
#include "Data/Boundary/CBoundaryConditionTorusSquare.h"
#include "Data/Boundary/CBoundaryConditionPeriodicAndDirichletSquare.h"
#include "Data/Boundary/CBoundaryConditionProjectivePlaneSquare.h"
#include "Data/Lattice/CIndex.h"
#include "Data/Lattice/CIndexSquare.h"
#include "Data/Lattice/CIndexData.h"
#include "Data/Lattice/CLatticeData.h"

//====================================
//define some common function before decompose threads

#define preparethread \
UINT block = _HC_DecompBlock; \
UINT threads = _HC_DecompThread; 

#define preparethreadHalf \
UINT block = _HC_DecompBlockHalf; \
UINT threads = _HC_DecompThreadHalf; 

#define preparethreadDir \
UINT block = _HC_DecompBlockDir; \
UINT threads = _HC_DecompThreadDir; 

#define intokernal \
const UINT uiSiteIndex = threadIdx.x + blockIdx.x * blockDim.x; \
if (uiSiteIndex >= _DC_Volume) \
{ \
    return; \
}


//To be discarded
#define intokernaldir \
intokernal \
const UINT uiDir = _DC_Dir;

//To be discarded
#define intokernalInt4dirC \
intokernal \
const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex); \
const UINT uiDir = _DC_Dir;

#define intokernalInt4 \
intokernal \
const SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

#define intokernalInt4EO \
intokernalInt4 \
if (_UBOOLXOR((pEtaTable[uiSiteIndex] >> 4U) & 1U, bEven)) \
{ \
    return; \
}

// Works for even lattice only
// NOTE: This is to SKIP even if bEven =1, or skip odd if bEven = 0
// 
// loop over volume/2
// at first, siteindex = 1, 3, 5, 7, 9, ... when bEven = 1, 
//       and siteindex = 0, 2, 4, 6, 8, ... when bEven = 0
// then, check whether this site is odd, if mask = isodd ^ needodd 
// (why not iseven ^ needeven? because if bEven=1, we are on odd sites, and D act on even (neighbour) sites, and give results on odd sites)
// so, bEven = need odd
// so:
// mask is:
//           need odd   0  1
// is odd       0       0  1
//              1       1  0
// if mask = 1, it means we are wrong, need to move to neighbor site
// If bEven,  we need uiSiteIndex - 1
// If !bEven, we need uiSiteIndex + 1
// why? see below, denote e=even, o=odd,
// 
// e o e o   0  1  2  3
// o e o e   4  5  6  7
// e o e o   8  9 10 11
// o e o e  12 13 14 15
// 
// when  bEven siteindex = 1,3,5,7,9,11,13,15; we need, 1,3,4,6,9,11,12,14
// when !bEven siteindex = 0,2,4,6,8,10,12,14; we need, 0,2,5,7,8,10,13,15   
// 
// So update uiSiteIndex = uiSiteIndex + (1 - 2 * bEven)
// 
// Do we need to check out of bound again? Note that, only when !bEven, which means we need even site, and we add one 
// However, if !bEven, before shift, the largest one is volume - 2, which gives volume -1 if there was a shift, still in bound
// update eta = pEtaTable[uiSiteIndex]
#define intokernalEOHalf \
UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) << 1U) | bEven; \
if (uiSiteIndex >= _DC_Volume) \
{ \
    return; \
} \
BYTE eta = pEtaTable[uiSiteIndex]; \
const BYTE mask = ((eta >> 4U) & 1U) ^ bEven; \
uiSiteIndex += mask * (1U - (bEven << 1U)); \
eta = pEtaTable[uiSiteIndex]; 


#define intokernalInt4NoConstant \
intokernal \
SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

//skip even
#define intokernalEO \
intokernal \
if (_UBOOLXOR((pEtaTable[uiSiteIndex] >> 4U) & 1U, bEven)) \
{ \
    return; \
}

//To be removed
#define intokernaldirEO \
intokernaldir \
if (_UBOOLXOR((pEtaTable[uiSiteIndex] >> 4U) & 1U, bEven)) \
{ \
    return; \
}

#if _CLG_ASSUME_SQUARE_LATTICE
#define intokernalDir_NoDir \
const UINT uiLinkIndex = threadIdx.x + blockIdx.x * blockDim.x; \
const UINT uiSiteIndex = (uiLinkIndex >> 2U); \
if (uiSiteIndex >= _DC_Volume) \
{ \
    return; \
}

#define intokernalDir \
intokernalDir_NoDir \
const UINT dir = (uiLinkIndex & 3U); 

#else
#define intokernalDir_NoDir \
const UINT uiLinkIndex = threadIdx.x + blockIdx.x * blockDim.x; \
const UINT uiSiteIndex = uiLinkIndex / _DC_Dir; \
if (uiSiteIndex >= _DC_Volume) \
{ \
    return; \
}

#define intokernalDir \
intokernalDir_NoDir \
const UINT dir = uiLinkIndex % _DC_Dir; 
#endif

#define intokernalDirInt4 \
intokernalDir \
SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

//skip even
#define intokernalDirEO \
intokernalDir \
if (_UBOOLXOR((pEtaTable[uiSiteIndex] >> 4U) & 1U, bEven)) \
{ \
    return; \
}

#define intokernalDirEO_NoDir \
intokernalDir_NoDir \
if (_UBOOLXOR((pEtaTable[uiSiteIndex] >> 4U) & 1U, bEven)) \
{ \
    return; \
}

#define intokernalDirInt4EO \
intokernalDirEO \
SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);

//make sure every site in one block
#define preparethreadE(element_count) \
UINT block, threads; \
appBlockThreadsE(_HC_Volume, element_count, block, threads); 

//make sure every link in one block
#define preparethreadEDir(element_count) \
UINT block, threads; \
appBlockThreadsE(_HC_Volume * _HC_Dir, element_count, block, threads); 

#define preparethreadEVar(element_count, blockvar, threadvar) \
UINT blockvar, threadvar; \
appBlockThreadsE(_HC_Volume, element_count, blockvar, threadvar); 

#define intokernalE(element_count)\
const UINT uiTotalIndex = threadIdx.x + blockIdx.x * blockDim.x; \
const UINT uiSiteIndex = uiTotalIndex / element_count; \
if (uiSiteIndex >= _DC_Volume) \
{ \
    return; \
} \
const BYTE elementIdx = uiTotalIndex % element_count; 

#if _CLG_ASSUME_SQUARE_LATTICE

#define intokernalEDir_NoDir(element_count)\
const UINT uiTotalIndex = threadIdx.x + blockIdx.x * blockDim.x; \
const UINT uiLinkIndex = uiTotalIndex / element_count; \
const UINT uiSiteIndex = uiLinkIndex >> 2U; \
if (uiSiteIndex >= _DC_Volume) \
{ \
    return; \
} \
const BYTE elementIdx = uiTotalIndex % element_count; 

#define intokernalEDir(element_count)\
intokernalEDir_NoDir(element_count) \
const BYTE dir = uiLinkIndex & 3U;

#else

#define intokernalEDir(element_count)\
const UINT uiTotalIndex = threadIdx.x + blockIdx.x * blockDim.x; \
const UINT uiLinkIndex = uiTotalIndex / element_count; \
const UINT uiSiteIndex = uiLinkIndex / _DC_Dir; \
if (uiSiteIndex >= _DC_Volume) \
{ \
    return; \
} \
const UINT elementIdx = uiTotalIndex % element_count; \
const UINT dir = uiLinkIndex % _DC_Dir;

#define intokernalEDir_NoDir(element_count)\
const UINT uiTotalIndex = threadIdx.x + blockIdx.x * blockDim.x; \
const UINT uiLinkIndex = uiTotalIndex / element_count; \
const UINT uiSiteIndex = uiLinkIndex / _DC_Dir; \
if (uiSiteIndex >= _DC_Volume) \
{ \
    return; \
} \
const UINT elementIdx = uiTotalIndex % element_count; 

#endif

//================= 3D threads, the x, y, z are x, y, and z =====================

#define preparethread_S \
const dim3 block3d(_HC_DecompX3D, _HC_DecompY3D, _HC_DecompZ3D); \
const dim3 threads3d(_HC_DecompLx3D, _HC_DecompLy3D, _HC_DecompLz3D);

#define intokernalInt4_S(uiT) \
SSmallInt4 sSite4; \
sSite4.x = static_cast<SCHAR>(threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.y = static_cast<SCHAR>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.z = static_cast<SCHAR>(threadIdx.z + blockIdx.z * blockDim.z); \
sSite4.w = uiT; \
const UINT uiSiteIndex = sSite4.x * _DC_MultX + sSite4.y * _DC_MultY + sSite4.z * _DC_MultZ + sSite4.w; \
const UINT uiSiteIndex3D = (sSite4.x * _DC_Ly + sSite4.y) * _DC_Lz + sSite4.z;

#define intokernalInt4_S_Only3D(uiT) \
SSmallInt4 sSite4; \
sSite4.x = static_cast<SCHAR>(threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.y = static_cast<SCHAR>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.z = static_cast<SCHAR>(threadIdx.z + blockIdx.z * blockDim.z); \
sSite4.w = uiT; \
const UINT uiSiteIndex3D = (sSite4.x * _DC_Ly + sSite4.y) * _DC_Lz + sSite4.z;


#define preparethread_Sxyt \
const dim3 block3dxyt(_HC_DecompX3DXYT, _HC_DecompY3DXYT, _HC_DecompZ3DXYT); \
const dim3 threads3dxyt(_HC_DecompLx3DXYT, _HC_DecompLy3DXYT, _HC_DecompLz3DXYT);

#define intokernalInt4_Sxyt(uiZ) \
SSmallInt4 sSite4; \
sSite4.x = static_cast<SCHAR>(threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.y = static_cast<SCHAR>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.w = static_cast<SCHAR>(threadIdx.z + blockIdx.z * blockDim.z); \
sSite4.z = uiZ; \
const UINT uiSiteIndex = sSite4.x * _DC_MultX + sSite4.y * _DC_MultY + sSite4.z * _DC_MultZ + sSite4.w; \
const UINT uiSiteIndex3DXYT = (sSite4.x * _DC_Ly + sSite4.y) * _DC_Lt + sSite4.w;

#define intokernalInt4_Sxyt_Only3D(uiZ) \
SSmallInt4 sSite4; \
sSite4.x = static_cast<SCHAR>(threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.y = static_cast<SCHAR>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.w = static_cast<SCHAR>(threadIdx.z + blockIdx.z * blockDim.z); \
sSite4.z = uiZ; \
const UINT uiSiteIndex3DXYT = (sSite4.x * _DC_Ly + sSite4.y) * _DC_Lt + sSite4.w;

#define preparethread_Sxzt \
const dim3 block3dxzt(_HC_DecompX3DXZT, _HC_DecompY3DXZT, _HC_DecompZ3DXZT); \
const dim3 threads3dxzt(_HC_DecompLx3DXZT, _HC_DecompLy3DXZT, _HC_DecompLz3DXZT);

#define intokernalInt4_Sxzt(uiY) \
SSmallInt4 sSite4; \
sSite4.x = static_cast<SCHAR>(threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.z = static_cast<SCHAR>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.w = static_cast<SCHAR>(threadIdx.z + blockIdx.z * blockDim.z); \
sSite4.y = uiY; \
const UINT uiSiteIndex = sSite4.x * _DC_MultX + sSite4.y * _DC_MultY + sSite4.z * _DC_MultZ + sSite4.w; \
const UINT uiSiteIndex3DXZT = (sSite4.x * _DC_Lz + sSite4.z) * _DC_Lt + sSite4.w;

#define intokernalInt4_Sxzt_Only3D(uiY) \
SSmallInt4 sSite4; \
sSite4.x = static_cast<SCHAR>(threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.z = static_cast<SCHAR>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.w = static_cast<SCHAR>(threadIdx.z + blockIdx.z * blockDim.z); \
sSite4.y = uiY; \
const UINT uiSiteIndex3DXZT = (sSite4.x * _DC_Lz + sSite4.z) * _DC_Lt + sSite4.w;

#define preparethread_Syzt \
const dim3 block3dyzt(_HC_DecompX3DYZT, _HC_DecompY3DYZT, _HC_DecompZ3DYZT); \
const dim3 threads3dyzt(_HC_DecompLx3DYZT, _HC_DecompLy3DYZT, _HC_DecompLz3DYZT);

#define intokernalInt4_Syzt(uiX) \
SSmallInt4 sSite4; \
sSite4.x = uiX; \
sSite4.y = static_cast<SCHAR>(threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.z = static_cast<SCHAR>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.w = static_cast<SCHAR>(threadIdx.z + blockIdx.z * blockDim.z); \
const UINT uiSiteIndex = sSite4.x * _DC_MultX + sSite4.y * _DC_MultY + sSite4.z * _DC_MultZ + sSite4.w; \
const UINT uiSiteIndex3DYZT = (sSite4.y * _DC_Lz + sSite4.z) * _DC_Lt + sSite4.w;

#define intokernalInt4_Syzt_Only3D(uiX) \
SSmallInt4 sSite4; \
sSite4.x = uiX; \
sSite4.y = static_cast<SCHAR>(threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.z = static_cast<SCHAR>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.w = static_cast<SCHAR>(threadIdx.z + blockIdx.z * blockDim.z); \
const UINT uiSiteIndex3DYZT = (sSite4.y * _DC_Lz + sSite4.z) * _DC_Lt + sSite4.w;

 //=======================================================
 //Field
#include "Data/Field/CField.h"
#include "Data/Field/CFieldCommonKernel.h"
#include "Data/Field/BoundaryField/CFieldBoundary.h"

#include "Data/Field/CFieldGauge.h"
#include "Data/Field/CFieldFermion.h"
#include "Data/Field/CFieldBoson.h"
#include "Data/Field/CFieldTensor2.h"

#include "SparseLinearAlgebra/CSLASolver.h"
#include "SparseLinearAlgebra/CMultiShiftSolver.h"

#include "Data/Field/Gauge/CFieldGaugeLink.h"
#include "Data/Field/Gauge/CFieldGaugeSU3_12.h"
#include "Data/Field/Staggered/CFieldFermionKS.h"
#include "Data/Field/Boson/CFieldBosonVNKernel.h"
#include "Data/Field/Boson/CFieldBosonT.h"

#include "Data/Field/Tensor2/CFieldTensor2Kernel.h"
#include "Data/Field/Tensor2/CFieldTensor2T.h"


//=====================================================

#include "Data/Action/CAction.h"

#include "Measurement/CMeasureData.h"
#include "Measurement/CMeasure.h"
#include "Measurement/CMeasurementManager.h"

#include "GaugeSmearing/CGaugeSmearing.h"

#include "GaugeFixing/CGaugeFixing.h"

#include "Tools/Math/DeviceInlineSU3.h"

#include "Update/CUpdator.h"
#include "Update/Continous/CIntegrator.h"
#include "Update/Continous/CHMC.h"
#include "Update/Discrete/CHeatbath.h"

#include "Core/CLGLibManager.h"

#endif //#ifndef _CLGLIB_PRIVATE_H_


//=============================================================================
// END OF FILE
//=============================================================================
