cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(CLG LANGUAGES CXX CUDA)

set(CMAKE_GENERATOR_PLATFORM x64)

# We start from CMAKE_SOURCE_DIR which should be /Code/CMake# First, we change it to /Code
set(CMAKE_SOURCE_DIR  ${CMAKE_SOURCE_DIR}/..)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/../Bin/UbuntuDebug)
set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR})
set(CMAKE_CURRENT_BINARY_DIR  ${CMAKE_BINARY_DIR})

# This is our code file dir
set(${PROJECT_SOURCE_DIR} ${CMAKE_SOURCE_DIR})
add_definitions(-D_UBUNTU)

include_directories(${PROJECT_SOURCE_DIR}/CLGLib)
add_subdirectory(${PROJECT_SOURCE_DIR}/CLGLib)

add_library(cudadevrt STATIC IMPORTED)
add_library(cudart_static STATIC IMPORTED)
add_library(cudrand_static STATIC IMPORTED)

add_library(CLGLib STATIC
    ${CMAKE_SOURCE_DIR}/CLGLib/CLGLib.h
    ${CMAKE_SOURCE_DIR}/CLGLib/CLGLib_Private.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Core/CBase.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Core/CCudaBuffer.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Core/CLGDefine.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Core/CLGFloat.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Core/CLGLibManager.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Core/CLGSetup.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Core/CudaHelper.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Core/CudaHelperFunctions.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Action/CAction.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Action/CActionFermion.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Action/CActionFermionWilsonNf2.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Action/CActionGauge.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquette.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryCondition.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionDirichlet.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionTorusSquare.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/CCommonData.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Field/CField.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Field/CFieldScalar.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Field/CFieldFermion.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Field/CFieldFermionWilsonSquareSU3.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Field/CFieldGauge.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Field/CFieldGaugeSU3.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Field/CFieldGaugeZ2.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Field/CFieldSpin.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Lattice/CIndex.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Lattice/CIndexSquare.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Lattice/CLatticeData.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Measurement/CMeasure.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Measurement/CMeasurementManager.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelator.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Measurement/CMeasurePlaqutteEnergy.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearing.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearingAPEProj.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearingAPEStout.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Platform/CFile.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Platform/PlatformIncs.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Platform/PlatformDefine.h
    ${CMAKE_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSLASolver.h
    ${CMAKE_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverBiCGstab.h
    ${CMAKE_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGCR.h
    ${CMAKE_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGMRES.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/CYAMLParser.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Data/CBitFlag.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Data/CCString.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Data/CLinkedList.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Data/STDStringFunctions.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Data/TArray.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Data/TemplateFunctions.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Data/THashMap.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Data/MemStack.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/EnumGather.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Math/GammaMatrix.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Math/CudaComplexFunction.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Math/Random.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Math/SU3.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Math/Vectors.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Tracer.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Timer.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Update/Continous/CHMC.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Update/Continous/CIntegrator.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorLeapFrog.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorOmelyan.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Update/CUpdator.h
    ${CMAKE_SOURCE_DIR}/CLGLib/Core/CudaHelper.cu
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryCondition.cu
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Field/CFieldFermionWilsonSquareSU3.cu
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Field/CFieldGauge.cu
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Lattice/CLatticeData.cu
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Lattice/CIndex.cu
    ${CMAKE_SOURCE_DIR}/CLGLib/Measurement/CMeasureMesonCorrelator.cu
    ${CMAKE_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearingAPEProj.cu
    ${CMAKE_SOURCE_DIR}/CLGLib/Measurement/GaugeSmearing/CGaugeSmearingAPEStout.cu
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Math/Random.cu
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Boundary/CBoundaryConditionTorusSquare.cu
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Field/CFieldGaugeSU3.cu
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Lattice/CIndexSquare.cu
    ${CMAKE_SOURCE_DIR}/CLGLib/Update/CUpdator.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/CLGLib_Private.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/Core/CLGLibManager.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGCR.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverGMRES.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/Update/Continous/CIntegratorOmelyan.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Action/CActionFermionWilsonNf2.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Field/CField.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/Measurement/CMeasurementManager.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/Measurement/CMeasurePlaqutteEnergy.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/Platform/CFile.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/CYAMLParser.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/SparseLinearAlgebra/CSolverBiCGstab.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Data/CCString.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/Tools/Tracer.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/Update/Continous/CHMC.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/Update/Continous/CIntegrator.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib/Data/Action/CActionGaugePlaquette.cpp
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

target_link_libraries(CLGLib libcurand_static.a libcudadevrt.a libcudart_static.a)




# ==================== 
# CLGTest 
# =================

include_directories(${PROJECT_SOURCE_DIR}/CLGTest)
add_subdirectory(${PROJECT_SOURCE_DIR}/CLGTest)

add_executable(CLGTest 
    ${CMAKE_SOURCE_DIR}/CLGTest/CLGTest.h
    ${CMAKE_SOURCE_DIR}/CLGTest/CLGTest.cpp
    ${CMAKE_SOURCE_DIR}/CLGTest/Tests/TestConfigurationFileIO.cpp
    ${CMAKE_SOURCE_DIR}/CLGTest/Tests/TestFermionUpdator.cpp
    ${CMAKE_SOURCE_DIR}/CLGTest/Tests/TestGammaMatrixSpeed.cpp
    ${CMAKE_SOURCE_DIR}/CLGTest/Tests/TestOperators.cpp
    ${CMAKE_SOURCE_DIR}/CLGTest/Tests/TestRandom.cpp
    ${CMAKE_SOURCE_DIR}/CLGTest/Tests/TestSolver.cpp
    ${CMAKE_SOURCE_DIR}/CLGTest/Tests/TestSU3Generator.cpp
    ${CMAKE_SOURCE_DIR}/CLGTest/Tests/TestUpdator.cpp
    )

target_compile_features(CLGTest PUBLIC cxx_std_14)
target_link_libraries(CLGTest CLGLib)
