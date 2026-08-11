cmake_minimum_required(VERSION 3.17 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(CLG LANGUAGES CXX)

set(CLG_BACKEND "CUDA" CACHE STRING "CLG backend (CUDA, HIP, DTK)")
set_property(CACHE CLG_BACKEND PROPERTY STRINGS "CUDA" "HIP" "DTK")
set(BACKEND "${CLG_BACKEND}")

set(CLG_GPU_ARCH "60" CACHE STRING "CLG GPU architecture")

if(BACKEND STREQUAL "CUDA")
    set(CMAKE_CUDA_ARCHITECTURES ${CLG_GPU_ARCH})
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3")
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
elseif(BACKEND STREQUAL "HIP")
    set(CMAKE_HIP_ARCHITECTURES ${CLG_GPU_ARCH})
    set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -O3 -Wno-return-type")
    enable_language(HIP)
    find_package(HIP REQUIRED)
    find_package(hipfft REQUIRED)
    find_package(hiprand REQUIRED)
elseif(BACKEND STREQUAL "DTK")
    set(AMDGPU_TARGETS ${CLG_GPU_ARCH})
    set(CMAKE_HIP_ARCHITECTURES ${CLG_GPU_ARCH})
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 -Wno-return-type --offload-arch=${CLG_GPU_ARCH}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fgpu-rdc --hip-link -fuse-ld=lld -flto=thin -Xoffload-linker --whole-archive")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fgpu-rdc --hip-link -fuse-ld=lld -flto=thin -Xoffload-linker --whole-archive")
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
endif()

# We start from CMAKE_SOURCE_DIR which should be /Code/CMake
option(CLG_DEBUG "CLG debug mode" OFF)

if(CLG_DEBUG)
    set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_SOURCE_DIR}/../../Bin/UbuntuDebug)
    message("Note: This is debug mode.")
else()
    set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_SOURCE_DIR}/../../Bin/Ubuntu)
endif()
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

# This is our code file dir
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/..)
if(CLG_DEBUG)
    add_definitions(-DDEBUG=1)
    add_definitions(-D_DEBUG=1)
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -Wno-unknown-pragmas -Wno-strict-overflow -Wno-maybe-uninitialized -Wno-return-type -flto=thin -fvisibility=hidden -fvisibility-inlines-hidden")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -Wno-unknown-pragmas -Wno-strict-overflow -Wno-maybe-uninitialized -Wno-class-memaccess -fvisibility=hidden -fvisibility-inlines-hidden")
    endif()
else()
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -Wall -Wno-unknown-pragmas -Wno-strict-overflow -Wno-return-type -flto=thin -fvisibility=hidden -fvisibility-inlines-hidden")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -Wall -Wno-unknown-pragmas -Wno-strict-overflow -Wno-class-memaccess -fvisibility=hidden -fvisibility-inlines-hidden")
    endif()
endif()
if(BACKEND STREQUAL "DTK")
    add_definitions(-D_CLG_DTK=1)
endif()
# to enable double float, add the following line:
if(DEFINED CLG_DOUBLE)
    add_definitions(-D_CLG_DOUBLEFLOAT=1)
    message("Note: double float is enabled.")
else()
    add_definitions(-D_CLG_DOUBLEFLOAT=0)
    message("Note: double float is disabled.")
endif()

# to enable profiler (synchronize in _RECORD to measure kernel time), add the following line:
if(DEFINED CLG_PROFILER)
    add_definitions(-D_CLG_PROFILER=1)
    message("Note: profiler is enabled.")
else()
    add_definitions(-D_CLG_PROFILER=0)
    message("Note: profiler is disabled.")
endif()

# to synchronize after each kernel launch (for debug), add the following line:
if(DEFINED CLG_CHECKSYNCHRONIZE)
    add_definitions(-D_CLG_CHECKSYNCHRONIZE=1)
    message("Note: synchronize after each kernel launch is enabled.")
else()
    add_definitions(-D_CLG_CHECKSYNCHRONIZE=0)
    message("Note: synchronize after each kernel launch is disabled.")
endif()

# to enable multi-GPU (one MPI process per GPU, lattice split into sub-lattices),
# use -DCLG_MULTI_GPU=1. MPI is only required in this case, so a machine without
# MPI can still configure and build the single-GPU library.
if(DEFINED CLG_MULTI_GPU)
    find_package(MPI REQUIRED)
    add_definitions(-D_CLG_MULTI_GPU=1)
    include_directories(${MPI_CXX_INCLUDE_DIRS})
    message("Note: multi-GPU is enabled. MPI: ${MPI_CXX_COMPILER}")
else()
    add_definitions(-D_CLG_MULTI_GPU=0)
    message("Note: multi-GPU is disabled.")
endif()

message("CMAKE_CUDA_FLAGS = ${CMAKE_CUDA_FLAGS}")
message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message("CMAKE_CUDA_ARCHITECTURES = ${CMAKE_CUDA_ARCHITECTURES}")
message("CMAKE_HIP_ARCHITECTURES = ${CMAKE_HIP_ARCHITECTURES}")

message("C++ Compiler Path: ${CMAKE_CXX_COMPILER}")
message("C++ Compiler ID: ${CMAKE_CXX_COMPILER_ID}")
message("C++ Compiler Version: ${CMAKE_CXX_COMPILER_VERSION}")

message("CUDA Compiler: ${CMAKE_CUDA_COMPILER}")
message("HIP Compiler: ${CMAKE_HIP_COMPILER}")

include_directories(${PROJECT_SOURCE_DIR}/CLGLib)
add_library(CLGCPULib STATIC
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/Version.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CAction.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionDiscreteGauge.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionFermionHISQCombined.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionFermionKS.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionFermionKSCombined.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionFermionKSImprove.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionFermionKSImproveCombined.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionPhi4.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryCondition.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Boson/CFieldBosonVN.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Boson/CFieldBosonVNRotation.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundaryOne.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundaryZero.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldBoson.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldFermion.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldGauge.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeLink.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeLinkDirichlet.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeSU3OneLoopImproved.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeSU3TreeImproved.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKS.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSHISQ.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSHISQAnisotropic.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3Gamma.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3GammaEM.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3REM.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKST.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSTD.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSTR.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareCloverSU3.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CLatticeData.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeSmearing/CGaugeSmearingHISQWithPhase.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAction.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureBosonCond.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureRotatingAction.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftFOM.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftCG.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftGMRES.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftBiCGStab.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CMultiShiftNested.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGMRESMDR.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverTFQMR.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverCG.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverDeflatedCG.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/CSimpleProfiler.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/CRationalApproximation.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/HostRandom.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorForceGradient.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorMultiLevelNestedForceGradient.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorMultiLevelNestedOmelyan.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNested11Stage.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNestedForceGradient.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNestedLeapFrog.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNestedOmelyan.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/CStapleCache.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/CUpdator.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/CLGLib_Private.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/CLGLibManager.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/Distributed/CLGComm.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/Distributed/CHaloBufferHandle.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/Distributed/CHaloManager.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGCR.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGMRES.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorOmelyan.cpp
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
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteAnisotropic.cpp
    )
target_precompile_headers(CLGCPULib PRIVATE ${PROJECT_SOURCE_DIR}/CLGLib/CLGLib_Private.h)
set(CLGLibSource
    ${PROJECT_SOURCE_DIR}/CLGLib/CLGLib_Private.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionDiscreteDNPlaquette.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionDiscreteZNPlaquette.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteAcceleration.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteAtGradient.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteBetaGradient.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteCylinder.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteBoost.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquettePSU3.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquettePSU3WithBoundary.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteRigidAcc.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteRotatingT.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquetteRotatingT3D.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionTemperatureDistribution.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionProjectivePlaneSquare.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Boson/CFieldBosonReal.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Boson/CFieldBosonVNKernel.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldCommonKernel.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeKernel.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeSU3_12.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Gauge/CFieldGaugeU1Real.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSHISQAnisotropicKernel.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSHISQWithPhase.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3Acc.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3RigidAcc.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSTKernel.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSTKernelGamma.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Staggered/CFieldFermionKSTKernelR.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Tensor2/CFieldTensor2.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Tensor2/CFieldTensor2Kernel.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/Tensor2/CFieldTensor2Z3.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonKernel.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3Acc.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3Boost.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3D.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3DR.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3DRigidAcc.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3EM.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3Gamma.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixing.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingCoulombCornell.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingCoulombLosAlamos.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingLandauLosAlamos.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingMAG.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingMCGDirect.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingMCGIndirect.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingRandom.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeSmearing/CGaugeSmearing.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeSmearing/CGaugeSmearingAPEProj.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeSmearing/CGaugeSmearingAPEStout.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeSmearing/CGaugeSmearingASQTAD.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeSmearing/CGaugeSmearingHISQ.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeSmearing/CGaugeSmearingStout.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAMomentumStochastic.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAngularMomentumKS.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAngularMomentumKSREM.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureBosonValue.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureChiralCondensate.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/GaugeFixing/CGaugeFixingLandauCornell.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasure.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureBerryPhase.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/CudaHelper.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionPeriodicAndDirichletSquare.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CIndexData.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureChargeAndCurrents.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureChiralCondensateKS.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureConnectedChiralSusceptibilityKS.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelator.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureAMomentumJG.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelatorStaggered.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelatorStaggeredSimple2.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePandChiralTalor.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePandChiralTalorKS.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePolyakovXY.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePolyakovXY3D.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureTopologicChargeXY.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureWilsonLoop.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureWilsonLoopWithPath.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureWilsonLoopXY.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/CLGFFT.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/CLinearAlgebraHelper.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/Random.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Discrete/CHeatbath.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionTorusSquare.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CIndexSquare.cu
    )
if(BACKEND STREQUAL "CUDA")
    if((DEFINED CLG_DOUBLE) AND (CLG_GPU_ARCH VERSION_LESS 60))
        message(FATAL_ERROR "CLG_GPU_ARCH must be >= 60 for double float atomicAdd")
    endif()
    add_library(CLGLib STATIC ${CLGLibSource})
    target_precompile_headers(CLGLib PRIVATE ${PROJECT_SOURCE_DIR}/CLGLib/CLGLib_Private.h)
    target_link_libraries(CLGCPULib CUDA::cudart CUDA::cufft CUDA::curand)
    target_link_libraries(CLGLib CUDA::cudart CUDA::cufft CUDA::curand)
elseif(BACKEND STREQUAL "HIP")
    set_source_files_properties(${CLGLibSource} PROPERTIES LANGUAGE HIP)
    add_library(CLGLib STATIC ${CLGLibSource})
    target_precompile_headers(CLGLib PRIVATE ${PROJECT_SOURCE_DIR}/CLGLib/CLGLib_Private.h)
    target_link_libraries(CLGLib hip::hipfft hip::hiprand)
    target_link_libraries(CLGCPULib hip::hipfft hip::hiprand)
elseif(BACKEND STREQUAL "DTK")
    add_library(CLGLib STATIC ${CLGLibSource})
    target_precompile_headers(CLGLib PRIVATE ${PROJECT_SOURCE_DIR}/CLGLib/CLGLib_Private.h)
    target_link_libraries(CLGCPULib CUDA::cudart CUDA::cufft CUDA::curand)
    target_link_libraries(CLGLib CUDA::cudart CUDA::cufft CUDA::curand)
endif()
set_target_properties(CLGLib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

# Multi-GPU: link MPI into both libraries. Only reached when -DCLG_MULTI_GPU=1,
# so single-GPU builds never require MPI to be present.
if(DEFINED CLG_MULTI_GPU)
    target_link_libraries(CLGCPULib MPI::MPI_CXX)
    target_link_libraries(CLGLib MPI::MPI_CXX)
endif()


# ==================== 
# CLGTest 
# =================

include_directories(${PROJECT_SOURCE_DIR}/CLGTest)
add_executable(CLGTest 
    ${PROJECT_SOURCE_DIR}/CLGTest/CLGTest.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestAngularMomentum.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestBoson.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestBoundary.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestConfigurationFileIO.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestMG.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestImprove1.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestImprove2.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestImprove5.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestImprove6.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestImprove8.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestFermionMatrix.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestFermionUpdator.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestCommon.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestHmcDiagnostics.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestKSUpdator.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestMultiShiftSolver.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestNumpyExport.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestRandom.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestRotationAcc.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestSolver.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestSUN.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestUpdator.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestGaugeFixing.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestPSU3WithBoundary.cpp
    )

target_link_libraries(CLGTest CLGLib CLGCPULib)

# Improve-1 I5 negative compile check (multi-GPU-improve1.md 3.3 gate).
# Opt-in, NOT part of ALL: compiles an intentionally type-mismatched
# _LAUNCH_KERNEL call and EXPECTS failure. Run manually on a
# _CLG_LAUNCH_KERNEL=0 (raw) build: `make CheckRawLaunchTypeFail`.
add_custom_target(CheckRawLaunchTypeFail
    COMMAND bash ${PROJECT_SOURCE_DIR}/Tools/check_raw_launch_typefail.sh
            ${CMAKE_BINARY_DIR} ${PROJECT_SOURCE_DIR}
    VERBATIM)




# ==================== 
# RotatingReproduce 
# =================

if (CLG_RotatingReproduce)
    include_directories(${PROJECT_SOURCE_DIR}/Applications/RotatingReproduce)
    add_executable(RotatingReproduce 
        ${PROJECT_SOURCE_DIR}/Applications/RotatingReproduce/AngularMomentumTest.cpp
        ${PROJECT_SOURCE_DIR}/Applications/RotatingReproduce/GaugeFixingJob.cpp
        ${PROJECT_SOURCE_DIR}/Applications/RotatingReproduce/PolykovDistrib.cpp
        ${PROJECT_SOURCE_DIR}/Applications/RotatingReproduce/RotatingReproduce.cpp
        ${PROJECT_SOURCE_DIR}/Applications/RotatingReproduce/RotationTest.cpp
        ${PROJECT_SOURCE_DIR}/Applications/RotatingReproduce/ThermalTest.cpp
            )

    target_link_libraries(RotatingReproduce CLGLib CLGCPULib)
else()
    message("project RotatingReproduce not built, to enable it using -DCLG_RotatingReproduce=1")
endif()




# ==================== 
# MatchingRho 
# =================

if (CLG_MatchingRho)
    include_directories(${PROJECT_SOURCE_DIR}/Applications/MatchingRho)
    add_executable(MatchingRho 
        ${PROJECT_SOURCE_DIR}/Applications/MatchingRho/MatchingRho.cpp
            )

    target_link_libraries(MatchingRho CLGLib CLGCPULib)
else()
    message("project MatchingRho not built, to enable it using -DCLG_MatchingRho=1")
endif()




# ==================== 
# ConstAcc 
# =================

if (CLG_ConstAcc)
    include_directories(${PROJECT_SOURCE_DIR}/Applications/ConstAcc)
    add_executable(ConstAcc 
        ${PROJECT_SOURCE_DIR}/Applications/ConstAcc/ConstAcc.cpp
        ${PROJECT_SOURCE_DIR}/Applications/ConstAcc/Measure.cpp
        ${PROJECT_SOURCE_DIR}/Applications/ConstAcc/Simulate.cpp
            )

    target_link_libraries(ConstAcc CLGLib CLGCPULib)
else()
    message("project ConstAcc not built, to enable it using -DCLG_ConstAcc=1")
endif()




# ==================== 
# StaggeredSpectrum 
# =================

if (CLG_StaggeredSpectrum)
    include_directories(${PROJECT_SOURCE_DIR}/Applications/StaggeredSpectrum)
    add_executable(StaggeredSpectrum 
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredSpectrum/GaugeFixing.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredSpectrum/Measurement.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredSpectrum/StaggeredSpectrum.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredSpectrum/Simulate.cpp
            )

    target_link_libraries(StaggeredSpectrum CLGLib CLGCPULib)
else()
    message("project StaggeredSpectrum not built, to enable it using -DCLG_StaggeredSpectrum=1")
endif()




# ==================== 
# StaggeredRotation 
# =================

if (CLG_StaggeredRotation)
    include_directories(${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation)
    add_executable(StaggeredRotation 
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/GaugeFixingJob.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/Measure.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/Measure3D.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/MeasureREM.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/MeasureU1.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/Simulate3D.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/SimulateBeta.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/SimulateREM.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/SimulateU1.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/SimulateU13D.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/StaggeredRotation.cpp
        ${PROJECT_SOURCE_DIR}/Applications/StaggeredRotation/Simulate.cpp
            )

    target_link_libraries(StaggeredRotation CLGLib CLGCPULib)
else()
    message("project StaggeredRotation not built, to enable it using -DCLG_StaggeredRotation=1")
endif()




# ==================== 
# BetaGradient 
# =================

if (CLG_BetaGradient)
    include_directories(${PROJECT_SOURCE_DIR}/Applications/BetaGradient)
    add_executable(BetaGradient 
        ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/GaugeFixingJob.cpp
        ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/MCGGaugeFixingJob.cpp
        ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/Measure.cpp
        ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/BetaGradient.cpp
        ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/MeasureBetaScan.cpp
        ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/MeasureTemperatureDistri.cpp
        ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/Simulate.cpp
        ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/SimulateAt.cpp
        ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/SimulateBetaScan.cpp
        ${PROJECT_SOURCE_DIR}/Applications/BetaGradient/SimulateTemperatureDistri.cpp
            )

    target_link_libraries(BetaGradient CLGLib CLGCPULib)
else()
    message("project BetaGradient not built, to enable it using -DCLG_BetaGradient=1")
endif()




# ==================== 
# ElectricChemical 
# =================

if (CLG_ElectricChemical)
    include_directories(${PROJECT_SOURCE_DIR}/Applications/ElectricChemical)
    add_executable(ElectricChemical 
        ${PROJECT_SOURCE_DIR}/Applications/ElectricChemical/GaugeFixingJob.cpp
        ${PROJECT_SOURCE_DIR}/Applications/ElectricChemical/Measure.cpp
        ${PROJECT_SOURCE_DIR}/Applications/ElectricChemical/ElectricChemical.cpp
        ${PROJECT_SOURCE_DIR}/Applications/ElectricChemical/MeasureImageChemical.cpp
        ${PROJECT_SOURCE_DIR}/Applications/ElectricChemical/Simulate.cpp
        ${PROJECT_SOURCE_DIR}/Applications/ElectricChemical/SimulateRW.cpp
            )

    target_link_libraries(ElectricChemical CLGLib CLGCPULib)
else()
    message("project ElectricChemical not built, to enable it using -DCLG_ElectricChemical=1")
endif()




# ==================== 
# HISQ 
# =================

if (CLG_HISQ)
    include_directories(${PROJECT_SOURCE_DIR}/Applications/HISQ)
    add_executable(HISQ 
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/CompareHISQ.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/DebugSU3_12.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/GaugeFixing.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/Measure.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/HISQ.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/Simulate.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/SimulateTest.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/Task1EO.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/Task1Full.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/Task2EO.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/Task2Full.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/Task3EO.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/Task4Action.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/Task5Plaq.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/Task6FermionAction.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/Task7Mom.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/TaskRotD.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/TaskRotFinalForce.cpp
        ${PROJECT_SOURCE_DIR}/Applications/HISQ/TaskRotGaugeForce.cpp
            )

    target_link_libraries(HISQ CLGLib CLGCPULib)
else()
    message("project HISQ not built, to enable it using -DCLG_HISQ=1")
endif()




# ==================== 
# WilsonDiracElectricMagnetic 
# =================

if (CLG_WilsonDiracElectricMagnetic)
    include_directories(${PROJECT_SOURCE_DIR}/Applications/WilsonDiracElectricMagnetic)
    add_executable(WilsonDiracElectricMagnetic 
        ${PROJECT_SOURCE_DIR}/Applications/WilsonDiracElectricMagnetic/WilsonDiracElectricMagnetic.cpp
            )

    target_link_libraries(WilsonDiracElectricMagnetic CLGLib CLGCPULib)
else()
    message("project WilsonDiracElectricMagnetic not built, to enable it using -DCLG_WilsonDiracElectricMagnetic=1")
endif()




# ==================== 
# RotationImprovedFermion 
# =================

if (CLG_RotationImprovedFermion)
    include_directories(${PROJECT_SOURCE_DIR}/Applications/RotationImprovedFermion)
    add_executable(RotationImprovedFermion 
        ${PROJECT_SOURCE_DIR}/Applications/RotationImprovedFermion/Measure.cpp
        ${PROJECT_SOURCE_DIR}/Applications/RotationImprovedFermion/RotationImprovedFermion.cpp
        ${PROJECT_SOURCE_DIR}/Applications/RotationImprovedFermion/Simulate.cpp
        ${PROJECT_SOURCE_DIR}/Applications/RotationImprovedFermion/SimulateWithRotation.cpp
            )

    target_link_libraries(RotationImprovedFermion CLGLib CLGCPULib)
else()
    message("project RotationImprovedFermion not built, to enable it using -DCLG_RotationImprovedFermion=1")
endif()

if (CLG_PSU3WithBoundary)
    include_directories(${PROJECT_SOURCE_DIR}/Applications/PSU3WithBoundary)
    add_executable(PSU3WithBoundary 
        ${PROJECT_SOURCE_DIR}/Applications/PSU3WithBoundary/PSU3WithBoundary.cpp
        ${PROJECT_SOURCE_DIR}/Applications/PSU3WithBoundary/Simulate.cpp
        ${PROJECT_SOURCE_DIR}/Applications/PSU3WithBoundary/Measure.cpp
            )

    target_link_libraries(PSU3WithBoundary CLGLib CLGCPULib)
else()
    message("project PSU3WithBoundary not built, to enable it using -DCLG_PSU3WithBoundary=1")
endif()

if (CLG_DEBUG)
  message("Will output to Bin/UbuntuDebug")
else()
  message("Will output to Bin/Ubuntu")
endif()
