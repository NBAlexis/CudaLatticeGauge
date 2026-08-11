//=============================================================================
// FILENAME : CLGLib.h
// 
// DESCRIPTION:
// This is the one header file for all
//
// REVISION:
//  [mm/dd/yy]
//  [12/2/2018 nbale]
//=============================================================================

#ifndef _CLGLIB_H_
#define _CLGLIB_H_

#define _CLG_PUBLIC 1

#include "Core/CLGSetup.h"
#include "Core/CLGDefine.h"

#if defined(_CLG_WIN)
#   if !defined(CLGAPI)
#       define __LIB_TITLE__    "CLGLib"
#       define CLGAPI __DLL_IMPORT
#       ifdef _CLG_DEBUG
#           define __LIB_FILE__    __LIB_TITLE__ "_d.lib"
#       else
#           define __LIB_FILE__ __LIB_TITLE__ ".lib"
#       endif
#       pragma __IMPORT_LIB(__LIB_FILE__)
#       pragma message("linking with " __LIB_FILE__ "...")
#       undef __LIB_FILE__
#       undef __LIB_TITLE__
#   endif
#else
#    define CLGAPI  
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

#include "Data/Field/BoundaryField/CFieldBoundaryZero.h"
#include "Data/Field/BoundaryField/CFieldBoundaryOne.h"

#include "Data/Field/Gauge/CFieldGaugeLink.h"
#include "Data/Field/Gauge/CFieldGaugeSU3_12.h"
#include "Data/Field/Gauge/CFieldGaugeLinkDirichlet.h"
#include "Data/Field/Gauge/CFieldGaugeU1Real.h"
#include "Data/Field/Gauge/CFieldGaugeZ2.h"
#include "Data/Field/Gauge/CFieldGaugeSU3TreeImproved.h"
#include "Data/Field/Gauge/CFieldGaugeSU3OneLoopImproved.h"

//#include "Data/Field/CFieldSpin.h"

#include "Data/Field/Staggered/CFieldFermionKS.h"
#include "Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3.h"
#include "Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3D.h"
#include "Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3DR.h"
#include "Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3Acc.h"
#include "Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3DRigidAcc.h"
#include "Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3Boost.h"
#include "Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3Gamma.h"
#include "Data/Field/WilsonDirac/CFieldFermionWilsonSquareCloverSU3.h"

#include "Data/Field/Staggered/CFieldFermionKST.h"
#include "Data/Field/Staggered/CFieldFermionKSTD.h"
#include "Data/Field/Staggered/CFieldFermionKSSU3Gamma.h"
#include "Data/Field/Staggered/CFieldFermionKSSU3GammaEM.h"
#include "Data/Field/Staggered/CFieldFermionKSSU3Acc.h"
#include "Data/Field/Staggered/CFieldFermionKSSU3RigidAcc.h"
#include "Data/Field/Staggered/CFieldFermionKSTR.h"
#include "Data/Field/Staggered/CFieldFermionKSSU3REM.h"
#include "Data/Field/Staggered/CFieldFermionKSHISQ.h"
#include "Data/Field/Staggered/CFieldFermionKSHISQAnisotropic.h"
#include "Data/Field/Staggered/CFieldFermionKSHISQAnisotropicKernel.h"
#include "Data/Field/Staggered/CFieldFermionKSHISQWithPhase.h"

#include "Data/Field/Boson/CFieldBosonVN.h"
#include "Data/Field/Boson/CFieldBosonVNRotation.h"
#include "Data/Field/Boson/CFieldBosonReal.h"

#include "Data/Field/Tensor2/CFieldTensor2Kernel.h"
#include "Data/Field/Tensor2/CFieldTensor2Fields.h"
#include "Data/Field/Tensor2/CFieldTensor2Z3.h"

//=====================================================

#include "Data/Action/CAction.h"
#include "Data/Action/CActionGaugePlaquette.h"
#include "Data/Action/CActionGaugePlaquetteAnisotropic.h"
#include "Data/Action/CActionGaugePlaquettePSU3.h"
#include "Data/Action/CActionGaugePlaquettePSU3WithBoundary.h"
#include "Data/Action/CActionGaugePlaquetteRotatingT.h"
#include "Data/Action/CActionFermionKS.h"
#include "Data/Action/CActionFermionKSImprove.h"
#include "Data/Action/CActionGaugePlaquetteAcceleration.h"
#include "Data/Action/CActionGaugePlaquetteBoost.h"
#include "Data/Action/CActionGaugePlaquetteRigidAcc.h"
#include "Data/Action/CActionGaugePlaquetteBetaGradient.h"
#include "Data/Action/CActionGaugePlaquetteAtGradient.h"
#include "Data/Action/CActionGaugePlaquetteCylinder.h"
#include "Data/Action/CActionGaugePlaquetteRotatingT3D.h"
#include "Data/Action/CActionPhi4.h"
#include "Data/Action/CActionTemperatureDistribution.h"
#include "Data/Action/CActionFermionKSCombined.h"
#include "Data/Action/CActionFermionKSImproveCombined.h"
#include "Data/Action/CActionFermionHISQCombined.h"

#include "SparseLinearAlgebra/CSolverBiCGstab.h"
#include "SparseLinearAlgebra/CSolverGMRES.h"
#include "SparseLinearAlgebra/CSolverGCR.h"
#include "SparseLinearAlgebra/CSolverGCRODR.h"
#include "SparseLinearAlgebra/CSolverGMRESMDR.h"
#include "SparseLinearAlgebra/CSolverTFQMR.h"
#include "SparseLinearAlgebra/CSolverCG.h"
#include "SparseLinearAlgebra/CSolverDeflatedCG.h"
#include "SparseLinearAlgebra/CMultiShiftGMRES.h"
#include "SparseLinearAlgebra/CMultiShiftFOM.h"
#include "SparseLinearAlgebra/CMultiShiftCG.h"
#include "SparseLinearAlgebra/CMultiShiftBiCGStab.h"
#include "SparseLinearAlgebra/CMultiShiftNested.h"

#include "Measurement/CMeasureData.h"
#include "Measurement/CMeasure.h"
#include "Measurement/CMeasurePlaqutteEnergy.h"
#include "Measurement/CMeasureAction.h"
#include "Measurement/CMeasureRotatingAction.h"
#include "Measurement/CMeasureMesonCorrelator.h"
#include "Measurement/CMeasureAMomentumJG.h"
#include "Measurement/CMeasureChargeAndCurrents.h"
#include "Measurement/CMeasureTopologicChargeXY.h"
#include "Measurement/CMeasurePolyakovXY.h"
#include "Measurement/CMeasureWilsonLoop.h"
#include "Measurement/CMeasureWilsonLoopXY.h"
#include "Measurement/CMeasureChiralCondensate.h"
#include "Measurement/CMeasureAMomentumStochastic.h"
#include "Measurement/CMeasureMesonCorrelatorStaggered.h"
#include "Measurement/CMeasureMesonCorrelatorStaggeredSimple2.h"
#include "Measurement/CMeasureChiralCondensateKS.h"
#include "Measurement/CMeasureAngularMomentumKS.h"
#include "Measurement/CMeasureConnectedChiralSusceptibilityKS.h"
#include "Measurement/CMeasureBerryPhase.h"
#include "Measurement/CMeasurePandChiralTalor.h"
#include "Measurement/CMeasureWilsonLoopWithPath.h"
#include "Measurement/CMeasureAngularMomentumKSREM.h"
#include "Measurement/CMeasurePolyakovXY3D.h"
#include "Measurement/CMeasurePandChiralTalorKS.h"
#include "Measurement/CMeasureBosonCond.h"
#include "Measurement/CMeasureBosonValue.h"

#include "Measurement/CMeasurementManager.h"

#include "GaugeSmearing/CGaugeSmearing.h"
#include "GaugeSmearing/CGaugeSmearingAPEProj.h"
#include "GaugeSmearing/CGaugeSmearingAPEStout.h"
#include "GaugeSmearing/CGaugeSmearingASQTAD.h"
#include "GaugeSmearing/CGaugeSmearingHISQ.h"
#include "GaugeSmearing/CGaugeSmearingHISQWithPhase.h"
#include "GaugeSmearing/CGaugeSmearingStout.h"

#include "GaugeFixing/CGaugeFixing.h"
#include "GaugeFixing/CGaugeFixingLandauCornell.h"
#include "GaugeFixing/CGaugeFixingCoulombCornell.h"
#include "GaugeFixing/CGaugeFixingLandauLosAlamos.h"
#include "GaugeFixing/CGaugeFixingCoulombLosAlamos.h"
#include "GaugeFixing/CGaugeFixingRandom.h"
#include "GaugeFixing/CGaugeFixingMAG.h"
#include "GaugeFixing/CGaugeFixingMCGDirect.h"
#include "GaugeFixing/CGaugeFixingMCGIndirect.h"

#include "Update/CUpdator.h"
#include "Update/CStapleCache.h"
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
#include "Update/Discrete/CHeatbath.h"

#include "Core/CLGLibManager.h"

#ifndef _CLG_PRIVATE
__USE_NAMESPACE
#endif

#endif //#ifndef _CLGLIB_H_

//=============================================================================
// END OF FILE
//=============================================================================
