cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

message("CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message("CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message("CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")

project(CLGLib LANGUAGES CXX CUDA)

include(CLGLib)

add_library(CLGLib STATIC
    CLGLib.h
    CLGLib_Private.h
    Core/CBase.h
    Core/CCudaBuffer.h
    Core/CLGDefine.h
    Core/CLGFloat.h
    Core/CLGLibManager.h
    Core/CLGSetup.h
    Core/CudaHelper.h
    Core/CudaHelperFunctions.h
    Data/Action/CAction.h
    Data/Action/CActionFermion.h
    Data/Action/CActionFermionWilsonNf2.h
    Data/Action/CActionGauge.h
    Data/Action/CActionGaugePlaquette.h
    Data/Boundary/CBoundaryCondition.h
    Data/Boundary/CBoundaryConditionDirichlet.h
    Data/Boundary/CBoundaryConditionTorusSquare.h
    Data/CCommonData.h
    Data/Field/CField.h
    Data/Field/CFieldScalar.h
    Data/Field/CFieldFermion.h
    Data/Field/CFieldFermionWilsonSquareSU3.h
    Data/Field/CFieldGauge.h
    Data/Field/CFieldGaugeSU3.h
    Data/Field/CFieldGaugeZ2.h
    Data/Field/CFieldSpin.h
    Data/Lattice/CIndex.h
    Data/Lattice/CIndexSquare.h
    Data/Lattice/CLatticeData.h
    Measurement/CMeasure.h
    Measurement/CMeasurementManager.h
    Measurement/CMeasureMesonCorrelator.h
    Measurement/CMeasurePlaqutteEnergy.h
    Measurement/GaugeSmearing/CGaugeSmearing.h
    Measurement/GaugeSmearing/CGaugeSmearingAPEProj.h
    Measurement/GaugeSmearing/CGaugeSmearingAPEStout.h
    Platform/CFile.h
    Platform/PlatformIncs.h
    Platform/WinDefine.h
    Platform/WinFunction.h
    SparseLinearAlgebra/CSLASolver.h
    SparseLinearAlgebra/CSolverBiCGstab.h
    SparseLinearAlgebra/CSolverGCR.h
    SparseLinearAlgebra/CSolverGMRES.h
    Tools/CYAMLParser.h
    Tools/Data/CBitFlag.h
    Tools/Data/CCString.h
    Tools/Data/CLinkedList.h
    Tools/Data/TArray.h
    Tools/Data/TemplateFunctions.h
    Tools/Data/THashMap.h
    Tools/Data/MemStack.h
    Tools/EnumGather.h
    Tools/Math/GammaMatrix.h
    Tools/Math/CudaComplexFunction.h
    Tools/Math/Random.h
    Tools/Math/SU3.h
    Tools/Math/Vectors.h
    Tools/Tracer.h
    Tools/Timer.h
    Update/Continous/CHMC.h
    Update/Continous/CIntegrator.h
    Update/Continous/CIntegratorLeapFrog.h
    Update/Continous/CIntegratorOmelyan.h
    Update/CUpdator.h
    Core/CudaHelper.cu
    Data/Boundary/CBoundaryCondition.cu
    Data/Field/CFieldFermionWilsonSquareSU3.cu
    Data/Field/CFieldGauge.cu
    Data/Lattice/CLatticeData.cu
    Data/Lattice/CIndex.cu
    Measurement/CMeasureMesonCorrelator.cu
    Measurement/GaugeSmearing/CGaugeSmearingAPEProj.cu
    Measurement/GaugeSmearing/CGaugeSmearingAPEStout.cu
    Tools/Math/Random.cu
    Data/Boundary/CBoundaryConditionTorusSquare.cu
    Data/Field/CFieldGaugeSU3.cu
    Data/Lattice/CIndexSquare.cu
    Update/CUpdator.cpp
    CLGLib_Private.cpp
    Core/CLGLibManager.cpp
    SparseLinearAlgebra/CSolverGCR.cpp
    SparseLinearAlgebra/CSolverGMRES.cpp
    Update/Continous/CIntegratorOmelyan.cpp
    Data/Action/CActionFermionWilsonNf2.cpp
    Data/Field/CField.cpp
    Measurement/CMeasurementManager.cpp
    Measurement/CMeasurePlaqutteEnergy.cpp
    Platform/CFileWin.cpp
    Platform/WinFunction.cpp
    Tools/CYAMLParser.cpp
    SparseLinearAlgebra/CSolverBiCGstab.cpp
    Tools/Data/CCString.cpp
    Tools/Tracer.cpp
    Update/Continous/CHMC.cpp
    Update/Continous/CIntegrator.cpp
    Data/Action/CActionGaugePlaquette.cpp
    )

# Request that CLGLib be built with -std=c++11
# As this is a public compile feature anything that links to 
# CLGLib will also build with -std=c++11
target_compile_features(CLGLib PUBLIC cxx_std_11)
 
# We need to explicitly state that we need all CUDA files in the 
# CLGLib library to be built with -dc as the member functions 
# could be called by other libraries and executables
set_target_properties( CLGLib
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)