cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(CLG LANGUAGES CXX CUDA)

set(CMAKE_GENERATOR_PLATFORM x64)

# We start from CMAKE_SOURCE_DIR which should be /Code/CMake
set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_SOURCE_DIR}/../../Bin/UbuntuDebug)
set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_CURRENT_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH  ${CMAKE_CURRENT_BINARY_DIR})
# This is our code file dir
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/..)
add_definitions(-D_UBUNTU)

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
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionFermion.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionFermionWilsonNf2.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGauge.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquette.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryCondition.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionPeriodicAndDirichletSquare.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionTorusSquare.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/CCommonData.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundary.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundaryGaugeSU3.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundaryWilsonSquareSU3.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CField.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldScalar.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldFermion.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldFermionWilsonSquareSU3.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldGauge.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldGaugeSU3.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldGaugeZ2.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldSpin.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CIndex.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CIndexSquare.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CLatticeData.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasure.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurementManager.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelator.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasurePlaqutteEnergy.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearing.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearingAPEProj.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearingAPEStout.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Platform/CFile.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Platform/PlatformIncs.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Platform/PlatformDefine.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSLASolver.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverBiCGstab.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGCR.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGCRODR.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGMRES.h
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGMRESMDR.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/CYAMLParser.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/CBitFlag.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/CCString.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/CLinkedList.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/STDStringFunctions.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/TArray.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/TemplateFunctions.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/THashMap.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Data/MemStack.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/EnumGather.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/CLinearAlgebraHelper.h
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
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNestedLeapFrog.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNestedOmelyan.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorNestedForceGradient.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorOmelyan.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/CUpdator.h
    ${PROJECT_SOURCE_DIR}/CLGLib/Core/CudaHelper.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryCondition.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionPeriodicAndDirichletSquare.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundaryGaugeSU3.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/BoundaryField/CFieldBoundaryWilsonSquareSU3.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldFermionWilsonSquareSU3.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldGauge.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CLatticeData.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CIndex.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelator.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearingAPEProj.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearingAPEStout.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/CLinearAlgebraHelper.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Tools/Math/Random.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionTorusSquare.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Field/CFieldGaugeSU3.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/Data/Lattice/CIndexSquare.cu
    ${PROJECT_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGMRESMDR.cpp
    ${PROJECT_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorForceGradient.cpp
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
 
# We need to explicitly state that we need all CUDA files in the 
# CLGLib library to be built with -dc as the member functions 
# could be called by other libraries and executables
set_target_properties( CLGLib
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_link_libraries(CLGLib -lcurand)




# ==================== 
# CLGTest 
# =================

include_directories(${PROJECT_SOURCE_DIR}/CLGTest)
add_executable(CLGTest 
    ${PROJECT_SOURCE_DIR}/CLGTest/CLGTest.h
    ${PROJECT_SOURCE_DIR}/CLGTest/CLGTest.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestBoundary.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestConfigurationFileIO.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestFermionUpdator.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestGammaMatrixSpeed.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestOperators.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestRandom.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestSolver.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestSU3Generator.cpp
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestUpdator.cpp
    )

target_compile_features(CLGTest PUBLIC cxx_std_14)
target_link_libraries(CLGTest CLGLib)
