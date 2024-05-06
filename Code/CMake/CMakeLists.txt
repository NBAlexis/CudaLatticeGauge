cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

if (DEFINED NVCCROOT)
    set(CMAKE_CUDA_COMPILER ${NVCCROOT})
    MESSAGE("CMAKE_CUDA_COMPILER = ${CMAKE_CUDA_COMPILER}")
endif()

set(CUDA_CMP "compute_61")
set(CUDA_SM "sm_61")
if (DEFINED CUDASM)
    set(CUDA_CMP "compute_${CUDASM}")
    set(CUDA_SM "sm_${CUDASM}")
    if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set (CMAKE_CUDA_ARCHITECTURES ${CUDASM})
    endif(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
endif()

project(CLG LANGUAGES CXX CUDA)

set(CMAKE_GENERATOR_PLATFORM x64)

# We start from CMAKE_SOURCE_DIR which should be /Code/CMake
set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_SOURCE_DIR}/../../Bin/UbuntuDebug)
set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_CURRENT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH  ${CMAKE_CURRENT_BINARY_DIR})
# This is our code file dir
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/..)
# Flags
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -Wall -Wno-unknown-pragmas -Wno-strict-overflow -Wno-class-memaccess")
add_definitions(-D_UBUNTU)
# to enable double float, add the following line:
if (DEFINED CLGDOUBLE)
  add_definitions(-D_CLG_DOUBLEFLOAT=1)
  MESSAGE("Note: double float is enabled, arch is ${CUDA_CMP} and ${CUDA_SM}.")
else()
  add_definitions(-D_CLG_DOUBLEFLOAT=0)
  MESSAGE("Note: double float is disabled, arch is ${CUDA_CMP} and ${CUDA_SM}.")
endif()
MESSAGE("CMAKE_CUDA_FLAGS flag = ${CMAKE_CUDA_FLAGS}")
MESSAGE("CMAKE_CXX_FLAGS flag = ${CMAKE_CXX_FLAGS}")

include_directories(${PROJECT_SOURCE_DIR}/CLGLib)
add_library(CLGLib STATIC
    ${PROJECT_SOURCE_DIR}/CLGLib/CLGLib.h
    ${PROJECT_SOURCE_DIR}/CLGLib/CLGLib_Private.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/CBase.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/CCudaBuffer.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/CLGDefine.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/CLGFloat.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/CLGLibManager.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/CLGSetup.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/CudaHelper.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/CudaHelperFunctions.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CAction.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionFermionKS.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionFermionWilsonNf2.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquette.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteAcceleration.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteBetaGradient.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteBoost.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquettePolar.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteRigidAcc.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteRotating.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteRotating3D.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteRotatingU1.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteRotatingU1_3D.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryCondition.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionPeriodicAndDirichletSquare.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionProjectivePlaneSquare.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionTorusSquare.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/CCommonData.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Boson/CFieldBosonU1.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundary.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundaryGaugeSU3.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundaryGaugeU1.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundaryWilsonSquareSU3.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CField.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldBoson.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldFermion.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldSpin.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGauge.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeSU3.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeSU3D.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeSU3TreeImproved.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeU1.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeU1Real.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeZ2.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKS.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3Acc.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3D.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3DR.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3EM.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3Gamma.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3GammaEM.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3P4.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3R.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3REM.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3RigidAcc.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSU1.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSU1R.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3Acc.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3Boost.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3D.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3DR.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3DRigidAcc.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3Gamma.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CIndex.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CIndexData.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CIndexSquare.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CLatticeData.h
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixing.h
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingCoulombCornell.h
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingCoulombLosAlamos.h
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingLandauCornell.h
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingLandauLosAlamos.h
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingRandom.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasure.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAction.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAMomentumStochastic.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAngularMomentumKS.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAngularMomentumKSREM.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureBerryPhase.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureChargeAndCurrents.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureChiralCondensate.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureChiralCondensateKS.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureConnectedChiralSusceptibilityKS.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurementManager.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelator.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAMomentumJG.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelatorStaggered.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelatorStaggeredSimple.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelatorStaggeredSimple2.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePandChiralTalor.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePandChiralTalorKS.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePlaqutteEnergy.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePolyakov.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePolyakovU1XY.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePolyakovXY.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePolyakovXY3D.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureTopologicChargeXY.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureWilsonLoop.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureWilsonLoopWithPath.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureWilsonLoopXY.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearing.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearingAPEProj.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearingAPEStout.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Platform/CFile.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Platform/PlatformIncs.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Platform/PlatformDefine.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftFOM.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftGMRES.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftNested.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftSolver.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftBiCGStab.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSLASolver.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverBiCGstab.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGCR.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGCRODR.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGMRES.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGMRESMDR.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverTFQMR.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/CYAMLParser.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/CBitFlag.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/CCString.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/CLGMD5.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/CLinkedList.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/STDStringFunctions.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/TArray.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/TemplateFunctions.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/THashMap.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/MemStack.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/EnumGather.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/CLGFFT.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/CLinearAlgebraHelper.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/CRationalApproximation.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/GammaMatrix.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/CudaComplexFunction.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/Random.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/SU3.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/Vectors.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Tracer.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Timer.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CHMC.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegrator.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorForceGradient.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorLeapFrog.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorMultiLevelNestedForceGradient.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorMultiLevelNestedOmelyan.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNested11Stage.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNestedLeapFrog.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNestedOmelyan.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNestedForceGradient.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorOmelyan.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/CUpdator.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteAcceleration.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteBetaGradient.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteBoost.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquettePolar.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteRigidAcc.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteRotating3D.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteRotatingU1.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteRotatingU1_3D.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionProjectivePlaneSquare.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Boson/CFieldBosonU1.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundaryGaugeU1.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGauge.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeSU3.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeSU3D.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeSU3TreeImproved.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeU1.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeU1Real.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3Acc.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3D.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3DR.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3EM.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3Gamma.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3GammaEM.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3P4.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3R.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3REM.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3RigidAcc.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSU1.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSU1R.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3Acc.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3Boost.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3D.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3DR.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3DRigidAcc.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3Gamma.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingCoulombCornell.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingCoulombLosAlamos.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingLandauLosAlamos.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingRandom.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAMomentumStochastic.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAngularMomentumKS.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAngularMomentumKSREM.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureChiralCondensate.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingLandauCornell.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasure.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureBerryPhase.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/CudaHelper.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteRotating.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionPeriodicAndDirichletSquare.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundaryGaugeSU3.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundaryWilsonSquareSU3.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CIndexData.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CLatticeData.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureChargeAndCurrents.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureChiralCondensateKS.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureConnectedChiralSusceptibilityKS.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelator.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAMomentumJG.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelatorStaggered.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelatorStaggeredSimple.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelatorStaggeredSimple2.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePandChiralTalor.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePandChiralTalorKS.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePolyakov.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePolyakovU1XY.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePolyakovXY.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePolyakovXY3D.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureTopologicChargeXY.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureWilsonLoop.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureWilsonLoopWithPath.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureWilsonLoopXY.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearing.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearingAPEProj.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearingAPEStout.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/CLGFFT.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/CLinearAlgebraHelper.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/Random.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionTorusSquare.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CIndexSquare.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CAction.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionFermionKS.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldBoson.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldFermion.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAction.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftFOM.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftGMRES.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftBiCGStab.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftNested.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGMRESMDR.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverTFQMR.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/HostRandom.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorForceGradient.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorMultiLevelNestedForceGradient.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorMultiLevelNestedOmelyan.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNested11Stage.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNestedForceGradient.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNestedLeapFrog.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNestedOmelyan.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/CUpdator.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/CLGLib_Private.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/CLGLibManager.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGCR.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGMRES.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorOmelyan.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionFermionWilsonNf2.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CField.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurementManager.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePlaqutteEnergy.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Platform/CFile.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/CYAMLParser.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverBiCGstab.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGCRODR.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/CCString.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Tracer.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CHMC.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegrator.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquette.cpp
    )

# Request that CLGLib be built with -std=c++14
# As this is a public compile feature anything that links to 
# CLGLib will also build with -std=c++14
target_compile_features(CLGLib PUBLIC cxx_std_14)
set_property(TARGET CLGLib PROPERTY CXX_STANDARD 14)
 
# We need to explicitly state that we need all CUDA files in the 
# CLGLib library to be built with -dc as the member functions 
# could be called by other libraries and executables
set_target_properties( CLGLib
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(CLGLib -lcurand)
target_link_libraries(CLGLib -lcufft)

# To enable the double, the minimum arch is 6.0
target_compile_options(CLGLib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=${CUDA_CMP},code=${CUDA_SM}>)




# ==================== 
# CLGTest 
# =================

include_directories(${PROJECT_SOURCE_DIR}/CLGTest)
add_executable(CLGTest 
    ${PROJECT_SOURCE_DIR}/CLGTest/CLGTest.h
    ${PROJECT_SOURCE_DIR}/CLGTest/CLGTest.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestAngularMomentum.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestBoundary.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestConfigurationFileIO.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestFermionMatrix.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestFermionUpdator.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestCommon.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestKSUpdator.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestMultiShiftSolver.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestRandom.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestRotationAcc.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestSolver.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestUpdator.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestGaugeFixing.cpp
    )

target_compile_features(CLGTest PUBLIC cxx_std_14)
target_link_libraries(CLGTest CLGLib)



# ==================== 
# ConstAcc 
# =================

include_directories(${PROJECT_SOURCE_DIR}/Applications/ConstAcc)
add_executable(ConstAcc 
    ${PROJECT_SOURCE_DIR}/Applications/ConstAcc/ConstAcc.h
    ${PROJECT_SOURCE_DIR}/Applications/ConstAcc/ConstAcc.cpp
    ${PROJECT_SOURCE_DIR}/Applications/ConstAcc/Measure.cpp
    ${PROJECT_SOURCE_DIR}/Applications/ConstAcc/Simulate.cpp
    )

target_compile_features(ConstAcc PUBLIC cxx_std_14)
target_link_libraries(ConstAcc CLGLib)




# ==================== 
# StaggeredRotation 
# =================

include_directories(${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation)
add_executable(StaggeredRotation 
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/StaggeredRotation.h
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/GaugeFixingJob.cpp
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/Measure.cpp
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/Measure3D.cpp
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/MeasureEM.cpp
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/MeasureREM.cpp
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/MeasureU1.cpp
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/Simulate3D.cpp
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/SimulateBeta.cpp
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/SimulateEM.cpp
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/SimulateREM.cpp
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/SimulateU1.cpp
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/SimulateU13D.cpp
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/StaggeredRotation.cpp
    ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/Simulate.cpp
    )

target_compile_features(StaggeredRotation PUBLIC cxx_std_14)
target_link_libraries(StaggeredRotation CLGLib)




# ==================== 
# BetaGradient 
# =================

include_directories(${PROJECT_SOURCE_DIR}/Applications/BetaGradient)
add_executable(BetaGradient 
    ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/BetaGradient.h
    ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/GaugeFixingJob.cpp
    ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/Measure.cpp
    ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/BetaGradient.cpp
    ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/MeasureBetaScan.cpp
    ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/Simulate.cpp
    ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/SimulateBetaScan.cpp
    )

target_compile_features(BetaGradient PUBLIC cxx_std_14)
target_link_libraries(BetaGradient CLGLib)




# ==================== 
# ElectricChemical 
# =================

include_directories(${PROJECT_SOURCE_DIR}/Applications/ElectricChemical)
add_executable(ElectricChemical 
    ${PROJECT_SOURCE_DIR}/Applications/ElectricChemical/ElectricChemical.h
    ${PROJECT_SOURCE_DIR}/Applications/ElectricChemical/GaugeFixingJob.cpp
    ${PROJECT_SOURCE_DIR}/Applications/ElectricChemical/Measure.cpp
    ${PROJECT_SOURCE_DIR}/Applications/ElectricChemical/ElectricChemical.cpp
    ${PROJECT_SOURCE_DIR}/Applications/ElectricChemical/MeasureImageChemical.cpp
    ${PROJECT_SOURCE_DIR}/Applications/ElectricChemical/Simulate.cpp
    ${PROJECT_SOURCE_DIR}/Applications/ElectricChemical/SimulateRW.cpp
    )

target_compile_features(ElectricChemical PUBLIC cxx_std_14)
target_link_libraries(ElectricChemical CLGLib)

