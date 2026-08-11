cmake_minimum_required(VERSION 3.17 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(AgentDriver LANGUAGES CXX CUDA)

set(CMAKE_CUDA_STANDARD 17)
set(CLG_GPU_ARCH "86" CACHE STRING "CLG GPU architecture")
set(CMAKE_CUDA_ARCHITECTURES ${CLG_GPU_ARCH})

set(CLG_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
set(CLG_BIN ${CLG_ROOT}/Bin/Ubuntu)

if(CLG_DEBUG)
    set(CLG_BIN ${CLG_ROOT}/Bin/UbuntuDebug)
endif()

set(CMAKE_CURRENT_BINARY_DIR ${CLG_BIN})
set(EXECUTABLE_OUTPUT_PATH ${CLG_BIN})

if(CLG_DOUBLE)
    add_definitions(-D_CLG_DOUBLEFLOAT=1)
    message("Note: double float is enabled.")
else()
    add_definitions(-D_CLG_DOUBLEFLOAT=0)
    message("Note: double float is disabled.")
endif()

include_directories(
    ${CLG_ROOT}/Code/CLGLib
    ${CMAKE_CURRENT_SOURCE_DIR}
)

find_package(CUDAToolkit REQUIRED)

add_executable(AgentDriver
    AgentDriver.cpp
    Simulate.cpp
    GaugeFixing.cpp
    Measure.cpp
)

set_target_properties(AgentDriver PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

target_link_directories(AgentDriver PRIVATE ${CLG_BIN})

if(EXISTS ${CLG_BIN}/libCLGLib.a AND EXISTS ${CLG_BIN}/libCLGCPULib.a)
    target_link_libraries(AgentDriver
        -Wl,--whole-archive CLGCPULib CLGLib -Wl,--no-whole-archive
        CUDA::cudart
        CUDA::cufft
        CUDA::curand)
else()
    message(FATAL_ERROR "CLGLib not found in ${CLG_BIN}. Please build CLGLib first using Code/CMake.")
endif()
