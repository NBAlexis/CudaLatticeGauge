cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(CLGLib LANGUAGES CXX CUDA)

set(CMAKE_GENERATOR_PLATFORM x64)

# We start from CMAKE_SOURCE_DIR which should be /Code/CMake# First, we change it to /Code/CLGLib
set(CMAKE_SOURCE_DIR  ${CMAKE_SOURCE_DIR}/../CLGLib)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/../../Bin/UbuntuDebug)
set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR})
set(CMAKE_CURRENT_BINARY_DIR  ${CMAKE_BINARY_DIR})

# This is our code file dir
set(${PROJECT_SOURCE_DIR} ${CMAKE_SOURCE_DIR})
include_directories("${PROJECT_SOURCE_DIR}")

message("CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message("CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message("CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")

include_directories(${CMAKE_SOURCE_DIR})

add_definitions(-D_UBUNTU)

add_library(CLGLib STATIC
    ${CMAKE_SOURCE_DIR}/CLGLib.h
    ${CMAKE_SOURCE_DIR}/CLGLib_Private.h
    ${CMAKE_SOURCE_DIR}/Core/CBase.h
    ${CMAKE_SOURCE_DIR}/Core/CCudaBuffer.h
    ${CMAKE_SOURCE_DIR}/Core/CLGDefine.h
    ${CMAKE_SOURCE_DIR}/Core/CLGFloat.h
    ${CMAKE_SOURCE_DIR}/Core/CLGLibManager.h
    ${CMAKE_SOURCE_DIR}/Core/CLGSetup.h
    ${CMAKE_SOURCE_DIR}/Core/CudaHelper.h
    ${CMAKE_SOURCE_DIR}/Core/CudaHelperFunctions.h
    ${CMAKE_SOURCE_DIR}/Data/Action/CAction.h
    ${CMAKE_SOURCE_DIR}/Data/Action/CActionFermion.h
    ${CMAKE_SOURCE_DIR}/Data/Action/CActionFermionWilsonNf2.h
    ${CMAKE_SOURCE_DIR}/Data/Action/CActionGauge.h
    ${CMAKE_SOURCE_DIR}/Data/Action/CActionGaugePlaquette.h
    ${CMAKE_SOURCE_DIR}/Data/Boundary/CBoundaryCondition.h
    ${CMAKE_SOURCE_DIR}/Data/Boundary/CBoundaryConditionDirichlet.h
    ${CMAKE_SOURCE_DIR}/Data/Boundary/CBoundaryConditionTorusSquare.h
    ${CMAKE_SOURCE_DIR}/Data/CCommonData.h
    ${CMAKE_SOURCE_DIR}/Data/Field/CField.h
    ${CMAKE_SOURCE_DIR}/Data/Field/CFieldScalar.h
    ${CMAKE_SOURCE_DIR}/Data/Field/CFieldFermion.h
    ${CMAKE_SOURCE_DIR}/Data/Field/CFieldFermionWilsonSquareSU3.h
    ${CMAKE_SOURCE_DIR}/Data/Field/CFieldGauge.h
    ${CMAKE_SOURCE_DIR}/Data/Field/CFieldGaugeSU3.h
    ${CMAKE_SOURCE_DIR}/Data/Field/CFieldGaugeZ2.h
    ${CMAKE_SOURCE_DIR}/Data/Field/CFieldSpin.h
    ${CMAKE_SOURCE_DIR}/Data/Lattice/CIndex.h
    ${CMAKE_SOURCE_DIR}/Data/Lattice/CIndexSquare.h
    ${CMAKE_SOURCE_DIR}/Data/Lattice/CLatticeData.h
    ${CMAKE_SOURCE_DIR}/Measurement/CMeasure.h
    ${CMAKE_SOURCE_DIR}/Measurement/CMeasurementManager.h
    ${CMAKE_SOURCE_DIR}/Measurement/CMeasureMesonCorrelator.h
    ${CMAKE_SOURCE_DIR}/Measurement/CMeasurePlaqutteEnergy.h
    ${CMAKE_SOURCE_DIR}/Measurement/GaugeSmearing/CGaugeSmearing.h
    ${CMAKE_SOURCE_DIR}/Measurement/GaugeSmearing/CGaugeSmearingAPEProj.h
    ${CMAKE_SOURCE_DIR}/Measurement/GaugeSmearing/CGaugeSmearingAPEStout.h
    ${CMAKE_SOURCE_DIR}/Platform/CFile.h
    ${CMAKE_SOURCE_DIR}/Platform/PlatformIncs.h
    ${CMAKE_SOURCE_DIR}/Platform/PlatformDefine.h
    ${CMAKE_SOURCE_DIR}/SparseLinearAlgebra/CSLASolver.h
    ${CMAKE_SOURCE_DIR}/SparseLinearAlgebra/CSolverBiCGstab.h
    ${CMAKE_SOURCE_DIR}/SparseLinearAlgebra/CSolverGCR.h
    ${CMAKE_SOURCE_DIR}/SparseLinearAlgebra/CSolverGMRES.h
    ${CMAKE_SOURCE_DIR}/Tools/CYAMLParser.h
    ${CMAKE_SOURCE_DIR}/Tools/Data/CBitFlag.h
    ${CMAKE_SOURCE_DIR}/Tools/Data/CCString.h
    ${CMAKE_SOURCE_DIR}/Tools/Data/CLinkedList.h
    ${CMAKE_SOURCE_DIR}/Tools/Data/STDStringFunctions.h
    ${CMAKE_SOURCE_DIR}/Tools/Data/TArray.h
    ${CMAKE_SOURCE_DIR}/Tools/Data/TemplateFunctions.h
    ${CMAKE_SOURCE_DIR}/Tools/Data/THashMap.h
    ${CMAKE_SOURCE_DIR}/Tools/Data/MemStack.h
    ${CMAKE_SOURCE_DIR}/Tools/EnumGather.h
    ${CMAKE_SOURCE_DIR}/Tools/Math/GammaMatrix.h
    ${CMAKE_SOURCE_DIR}/Tools/Math/CudaComplexFunction.h
    ${CMAKE_SOURCE_DIR}/Tools/Math/Random.h
    ${CMAKE_SOURCE_DIR}/Tools/Math/SU3.h
    ${CMAKE_SOURCE_DIR}/Tools/Math/Vectors.h
    ${CMAKE_SOURCE_DIR}/Tools/Tracer.h
    ${CMAKE_SOURCE_DIR}/Tools/Timer.h
    ${CMAKE_SOURCE_DIR}/Update/Continous/CHMC.h
    ${CMAKE_SOURCE_DIR}/Update/Continous/CIntegrator.h
    ${CMAKE_SOURCE_DIR}/Update/Continous/CIntegratorLeapFrog.h
    ${CMAKE_SOURCE_DIR}/Update/Continous/CIntegratorOmelyan.h
    ${CMAKE_SOURCE_DIR}/Update/CUpdator.h
    ${CMAKE_SOURCE_DIR}/Core/CudaHelper.cu
    ${CMAKE_SOURCE_DIR}/Data/Boundary/CBoundaryCondition.cu
    ${CMAKE_SOURCE_DIR}/Data/Field/CFieldFermionWilsonSquareSU3.cu
    ${CMAKE_SOURCE_DIR}/Data/Field/CFieldGauge.cu
    ${CMAKE_SOURCE_DIR}/Data/Lattice/CLatticeData.cu
    ${CMAKE_SOURCE_DIR}/Data/Lattice/CIndex.cu
    ${CMAKE_SOURCE_DIR}/Measurement/CMeasureMesonCorrelator.cu
    ${CMAKE_SOURCE_DIR}/Measurement/GaugeSmearing/CGaugeSmearingAPEProj.cu
    ${CMAKE_SOURCE_DIR}/Measurement/GaugeSmearing/CGaugeSmearingAPEStout.cu
    ${CMAKE_SOURCE_DIR}/Tools/Math/Random.cu
    ${CMAKE_SOURCE_DIR}/Data/Boundary/CBoundaryConditionTorusSquare.cu
    ${CMAKE_SOURCE_DIR}/Data/Field/CFieldGaugeSU3.cu
    ${CMAKE_SOURCE_DIR}/Data/Lattice/CIndexSquare.cu
    ${CMAKE_SOURCE_DIR}/Update/CUpdator.cpp
    ${CMAKE_SOURCE_DIR}/CLGLib_Private.cpp
    ${CMAKE_SOURCE_DIR}/Core/CLGLibManager.cpp
    ${CMAKE_SOURCE_DIR}/SparseLinearAlgebra/CSolverGCR.cpp
    ${CMAKE_SOURCE_DIR}/SparseLinearAlgebra/CSolverGMRES.cpp
    ${CMAKE_SOURCE_DIR}/Update/Continous/CIntegratorOmelyan.cpp
    ${CMAKE_SOURCE_DIR}/Data/Action/CActionFermionWilsonNf2.cpp
    ${CMAKE_SOURCE_DIR}/Data/Field/CField.cpp
    ${CMAKE_SOURCE_DIR}/Measurement/CMeasurementManager.cpp
    ${CMAKE_SOURCE_DIR}/Measurement/CMeasurePlaqutteEnergy.cpp
    ${CMAKE_SOURCE_DIR}/Platform/CFile.cpp
    ${CMAKE_SOURCE_DIR}/Tools/CYAMLParser.cpp
    ${CMAKE_SOURCE_DIR}/SparseLinearAlgebra/CSolverBiCGstab.cpp
    ${CMAKE_SOURCE_DIR}/Tools/Data/CCString.cpp
    ${CMAKE_SOURCE_DIR}/Tools/Tracer.cpp
    ${CMAKE_SOURCE_DIR}/Update/Continous/CHMC.cpp
    ${CMAKE_SOURCE_DIR}/Update/Continous/CIntegrator.cpp
    ${CMAKE_SOURCE_DIR}/Data/Action/CActionGaugePlaquette.cpp
    )

# Request that CLGLib be built with -std=c++14
# As this is a public compile feature anything that links to 
# CLGLib will also build with -std=c++14
target_compile_features(CLGLib PUBLIC cxx_std_14)
 
# We need to explicitly state that we need all CUDA files in the 
# CLGLib library to be built with -dc as the member functions 
# could be called by other libraries and executables
set_target_properties( CLGLib
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)/n/ntarget_link_libraries(CLGLib libcurand.so libcudadevrt.a libcudart.so)




# ==================== 
# CLGTest 
# =================

project(CLGTest LANGUAGES CXX)

set(CMAKE_GENERATOR_PLATFORM x64)

# We start from CMAKE_SOURCE_DIR which should be /Code/CMake# First, we change it to /Code/CLGTest
set(CMAKE_SOURCE_DIR  ${CMAKE_SOURCE_DIR}/../CLGTest)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/../../Bin/UbuntuDebug)
set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR})
set(CMAKE_CURRENT_BINARY_DIR  ${CMAKE_BINARY_DIR})

# This is our code file dir
set(${PROJECT_SOURCE_DIR} ${CMAKE_SOURCE_DIR})
include_directories("${PROJECT_SOURCE_DIR}")

message("CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message("CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message("CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")

include_directories(${CMAKE_SOURCE_DIR}/../CLGLib)
include_directories(${CMAKE_SOURCE_DIR})

add_definitions(-D_UBUNTU)

add_executable(CLGTest 
    ${CMAKE_SOURCE_DIR}/CLGTest.h
    ${CMAKE_SOURCE_DIR}/CLGTest.cpp
    ${CMAKE_SOURCE_DIR}/Tests/TestConfigurationFileIO.cpp
    ${CMAKE_SOURCE_DIR}/Tests/TestFermionUpdator.cpp
    ${CMAKE_SOURCE_DIR}/Tests/TestGammaMatrixSpeed.cpp
    ${CMAKE_SOURCE_DIR}/Tests/TestOperators.cpp
    ${CMAKE_SOURCE_DIR}/Tests/TestRandom.cpp
    ${CMAKE_SOURCE_DIR}/Tests/TestSolver.cpp
    ${CMAKE_SOURCE_DIR}/Tests/TestSU3Generator.cpp
    ${CMAKE_SOURCE_DIR}/Tests/TestUpdator.cpp
    )

target_compile_features(CLGTest PUBLIC cxx_std_14)
target_link_libraries(CLGTest CLGLib)
