using System;
using System.Collections.Generic;
using System.IO;


namespace CLGMakeWriter
{

    enum EArch
    {
        EArcSM52,
        EArcSM61,
        EArcSM70,
        EArcSM75,
    }

    class CMakeWritter
    {
        public bool m_bDouble = false;
        public EArch m_eArch = EArch.EArcSM61;
        public bool m_bDebug = true;
        public bool m_bWinOrUbuntu = true;

        public bool m_bHasWilsonDiracRotation = true;
        public bool m_bHasWilsonDiracMatching = false;

        //turn it off, recently we will not use this
        public bool m_bHasCompresser = false;

        //turn it off, we have sign problem here
        public bool m_bHasConstAcc = false;

        public bool m_bHasStaggeredSpectrum = true;
        
        readonly static string[] FileSurfix = { "_DebugMSVC.txt", "_ReleaseMSVC.txt", "_DebugGCC.txt", "_ReleaseGCC.txt" };
        readonly static string[] ArchNames = 
        {
            //GTX 970M
            "compute_52",

            //GTX 1060, 1070
            "compute_61",

            //V100
            "compute_70",

            //RTX2060-2080
            "compute_75",
        };

        readonly static string[] CodeNames =
        {
            //GTX 970M
            "sm_52",

            //GTX 1060, 1070
            "sm_61",

            //V100
            "sm_70",

            //RTX2060-2080
            "sm_75",
        };

        public void WritteTheFile(string sSolDir, CProjFile projFile, Dictionary<string, CProjFile> excutables)
        {
            string sContent = sCmakeContentOrignal1;

            #region Add CLGLib

            sContent += "add_library(CLGCPULib STATIC\n    ";
            foreach (string sFileName in projFile.m_lstAllCppFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/CLGLib/" + sFileName + "\n    ";
            }
            sContent += ")\n";

            sContent += "target_precompile_headers(CLGCPULib PRIVATE ${PROJECT_SOURCE_DIR}/CLGLib/CLGLib_Private.h)\nset(CLGLibSource\n    ";
            foreach (string sFileName in projFile.m_lstAllCuFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/CLGLib/" + sFileName + "\n    ";
            }

            sContent += sCmakeContentOrignal2;

            #endregion

            #region Add CLGTest

            CProjFile slgTest = excutables["CLGTest"];

            sContent += "\n\n\n# ==================== \n# CLGTest \n# =================\n\n";

            sContent += "include_directories(${PROJECT_SOURCE_DIR}/CLGTest)\n";

            sContent += "add_executable(CLGTest \n    ";
            foreach (string sFileName in slgTest.m_lstAllCppFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/CLGTest/" + sFileName + "\n    ";
            }
            sContent += ")\n\n";
            sContent += "target_link_libraries(CLGTest CLGLib CLGCPULib)\n\n";

            #endregion

            #region Add Applications



            //if (m_bHasCompresser)
            //{
            //    sContent += AddApplication(excutables["ConfigurationCompresser"]);
            //}

            sContent += AddApplication(excutables["RotatingReproduce"]);
            sContent += AddApplication(excutables["MatchingRho"]);
            sContent += AddApplication(excutables["ConstAcc"]);
            sContent += AddApplication(excutables["StaggeredSpectrum"]);
            sContent += AddApplication(excutables["StaggeredRotation"]);
            sContent += AddApplication(excutables["BetaGradient"]);
            sContent += AddApplication(excutables["ElectricChemical"]);
            sContent += AddApplication(excutables["HISQ"]);
            sContent += AddApplication(excutables["WilsonDiracElectricMagnetic"]);
            sContent += AddApplication(excutables["RotationImprovedFermion"]);

            //sContent += AddApplication(excutables["CLGExample"]);

            #endregion

            sContent += "if (CLG_DEBUG)\n";
            sContent += "  message(\"Will output to Bin/UbuntuDebug\")\n";
            sContent += "else()\n";
            sContent += "  message(\"Will output to Bin/Ubuntu\")\n";
            sContent += "endif()\n";

            sContent = sContent.Replace("\r\n", "\n");
            sContent = sContent.Replace("\n\r", "\n");
            sContent = sContent.Replace("\r", "\n");
            sContent = sContent.Replace("\\", "/");

            File.WriteAllText(sSolDir + "CMake/CMakeLists.txt", sContent);

            //File.WriteAllText(sSolDir +"CMake/CMakeLists" 
            //    + (m_bWinOrUbuntu ? (m_bDebug ? FileSurfix[0] : FileSurfix[1]) : (m_bDebug ? FileSurfix[2] : FileSurfix[3])), sContent);1
        }

        protected string AddApplication(CProjFile addProj)
        {
            string sAppName = addProj.m_sName;

            var sRet = string.Format("\n\n\n# ==================== \n# {0} \n# =================\n\n", sAppName);
            sRet += string.Format("if (CLG_{0})\n", sAppName);
            sRet += string.Format("    include_directories({1}/Applications/{0})\n", sAppName, "${PROJECT_SOURCE_DIR}");
            sRet += string.Format("    add_executable({0} \n        ", sAppName);

            foreach (string sFileName in addProj.m_lstAllCppFiles)
            {
                sRet += string.Format("{1}/Applications/{0}/" + sFileName + "\n        ", sAppName, "${PROJECT_SOURCE_DIR}");
            }

            sRet += "    )\n\n";
            sRet += string.Format("    target_link_libraries({0} CLGLib CLGCPULib)\n", sAppName);
            sRet += "else()\n";
            sRet += string.Format("    message(\"project {0} not built, to enable it using -DCLG_{0}=1\")\n", sAppName);
            sRet += "endif()\n\n";
            return sRet;
        }

        public string sCmakeContentOrignal1 = @"cmake_minimum_required(VERSION 3.17 FATAL_ERROR)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(CLG LANGUAGES CXX)

set(CLG_BACKEND ""CUDA"" CACHE STRING ""CLG backend (CUDA, HIP, DTK)"")
set_property(CACHE CLG_BACKEND PROPERTY STRINGS ""CUDA"" ""HIP"" ""DTK"")
set(BACKEND ""${CLG_BACKEND}"")

set(CLG_GPU_ARCH ""60"" CACHE STRING ""CLG GPU architecture"")

if(BACKEND STREQUAL ""CUDA"")
    set(CMAKE_CUDA_ARCHITECTURES ${CLG_GPU_ARCH})
    set(CMAKE_CUDA_FLAGS ""${CMAKE_CUDA_FLAGS} -O3"")
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
elseif(BACKEND STREQUAL ""HIP"")
    set(CMAKE_HIP_ARCHITECTURES ${CLG_GPU_ARCH})
    set(CMAKE_HIP_FLAGS ""${CMAKE_HIP_FLAGS} -O3 -Wno-return-type"")
    enable_language(HIP)
    find_package(HIP REQUIRED)
    find_package(hipfft REQUIRED)
    find_package(hiprand REQUIRED)
elseif(BACKEND STREQUAL ""DTK"")
    set(AMDGPU_TARGETS ${CLG_GPU_ARCH})
    set(CMAKE_HIP_ARCHITECTURES ${CLG_GPU_ARCH})
    set(CMAKE_CUDA_FLAGS ""${CMAKE_CUDA_FLAGS} -O3 -Wno-return-type --offload-arch=${CLG_GPU_ARCH}"")
    set(CMAKE_EXE_LINKER_FLAGS ""${CMAKE_EXE_LINKER_FLAGS} -fgpu-rdc --hip-link -fuse-ld=lld -flto=thin -Xoffload-linker --whole-archive"")
    set(CMAKE_SHARED_LINKER_FLAGS ""${CMAKE_EXE_LINKER_FLAGS} -fgpu-rdc --hip-link -fuse-ld=lld -flto=thin -Xoffload-linker --whole-archive"")
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
endif()

# We start from CMAKE_SOURCE_DIR which should be /Code/CMake
option(CLG_DEBUG ""CLG debug mode"" OFF)

if(CLG_DEBUG)
    set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_SOURCE_DIR}/../../Bin/UbuntuDebug)
    message(""Note: This is debug mode."")
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
    if(CMAKE_CXX_COMPILER_ID STREQUAL ""Clang"")
        set(CMAKE_CXX_FLAGS ""${CMAKE_CXX_FLAGS} -g -O0 -Wall -Wno-unknown-pragmas -Wno-strict-overflow -Wno-maybe-uninitialized -Wno-return-type -flto=thin -fvisibility=hidden -fvisibility-inlines-hidden"")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL ""GNU"")
        set(CMAKE_CXX_FLAGS ""${CMAKE_CXX_FLAGS} -g -O0 -Wall -Wno-unknown-pragmas -Wno-strict-overflow -Wno-maybe-uninitialized -Wno-class-memaccess -fvisibility=hidden -fvisibility-inlines-hidden"")
    endif()
else()
    if(CMAKE_CXX_COMPILER_ID STREQUAL ""Clang"")
        set(CMAKE_CXX_FLAGS ""${CMAKE_CXX_FLAGS} -Ofast -Wall -Wno-unknown-pragmas -Wno-strict-overflow -Wno-return-type -flto=thin -fvisibility=hidden -fvisibility-inlines-hidden"")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL ""GNU"")
        set(CMAKE_CXX_FLAGS ""${CMAKE_CXX_FLAGS} -Ofast -Wall -Wno-unknown-pragmas -Wno-strict-overflow -Wno-class-memaccess -fvisibility=hidden -fvisibility-inlines-hidden"")
    endif()
endif()
if(BACKEND STREQUAL ""DTK"")
    add_definitions(-D_CLG_DTK=1)
endif()
# to enable double float, add the following line:
if(DEFINED CLG_DOUBLE)
    add_definitions(-D_CLG_DOUBLEFLOAT=1)
    message(""Note: double float is enabled."")
else()
    add_definitions(-D_CLG_DOUBLEFLOAT=0)
    message(""Note: double float is disabled."")
endif()

# to enable profiler (synchronize in _RECORD to measure kernel time), add the following line:
if(DEFINED CLG_PROFILER)
    add_definitions(-D_CLG_PROFILER=1)
    message(""Note: profiler is enabled."")
else()
    add_definitions(-D_CLG_PROFILER=0)
    message(""Note: profiler is disabled."")
endif()

# to synchronize after each kernel launch (for debug), add the following line:
if(DEFINED CLG_CHECKSYNCHRONIZE)
    add_definitions(-D_CLG_CHECKSYNCHRONIZE=1)
    message(""Note: synchronize after each kernel launch is enabled."")
else()
    add_definitions(-D_CLG_CHECKSYNCHRONIZE=0)
    message(""Note: synchronize after each kernel launch is disabled."")
endif()

message(""CMAKE_CUDA_FLAGS = ${CMAKE_CUDA_FLAGS}"")
message(""CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}"")
message(""CMAKE_CUDA_ARCHITECTURES = ${CMAKE_CUDA_ARCHITECTURES}"")
message(""CMAKE_HIP_ARCHITECTURES = ${CMAKE_HIP_ARCHITECTURES}"")

message(""C++ Compiler Path: ${CMAKE_CXX_COMPILER}"")
message(""C++ Compiler ID: ${CMAKE_CXX_COMPILER_ID}"")
message(""C++ Compiler Version: ${CMAKE_CXX_COMPILER_VERSION}"")

message(""CUDA Compiler: ${CMAKE_CUDA_COMPILER}"")
message(""HIP Compiler: ${CMAKE_HIP_COMPILER}"")

include_directories(${PROJECT_SOURCE_DIR}/CLGLib)
";
        public string sCmakeContentOrignal2 = @")
if(BACKEND STREQUAL ""CUDA"")
    if((DEFINED CLG_DOUBLE) AND (CLG_GPU_ARCH VERSION_LESS 60))
        message(FATAL_ERROR ""CLG_GPU_ARCH must be >= 60 for double float atomicAdd"")
    endif()
    add_library(CLGLib STATIC ${CLGLibSource})
    target_precompile_headers(CLGLib PRIVATE ${PROJECT_SOURCE_DIR}/CLGLib/CLGLib_Private.h)
    target_link_libraries(CLGCPULib CUDA::cudart CUDA::cufft CUDA::curand)
    target_link_libraries(CLGLib CUDA::cudart CUDA::cufft CUDA::curand)
elseif(BACKEND STREQUAL ""HIP"")
    set_source_files_properties(${CLGLibSource} PROPERTIES LANGUAGE HIP)
    add_library(CLGLib STATIC ${CLGLibSource})
    target_precompile_headers(CLGLib PRIVATE ${PROJECT_SOURCE_DIR}/CLGLib/CLGLib_Private.h)
    target_link_libraries(CLGLib hip::hipfft hip::hiprand)
    target_link_libraries(CLGCPULib hip::hipfft hip::hiprand)
elseif(BACKEND STREQUAL ""DTK"")
    add_library(CLGLib STATIC ${CLGLibSource})
    target_precompile_headers(CLGLib PRIVATE ${PROJECT_SOURCE_DIR}/CLGLib/CLGLib_Private.h)
    target_link_libraries(CLGCPULib CUDA::cudart CUDA::cufft CUDA::curand)
    target_link_libraries(CLGLib CUDA::cudart CUDA::cufft CUDA::curand)
endif()
set_target_properties(CLGLib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)";
    }
}