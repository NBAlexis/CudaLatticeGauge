
using System.Collections.Generic;
using System.IO;

namespace CLGMakeWriter
{

    enum EArch
    {
        EArcSM52,
        EArcSM61,
        EArcSM70,
    }

    class CMakeWritter
    {
        public bool m_bDouble = false;
        public EArch m_eArch = EArch.EArcSM61;
        public bool m_bDebug = true;
        public bool m_bWinOrUbuntu = true;
        
        readonly static string[] FileSurfix = { "_DebugMSVC.txt", "_ReleaseMSVC.txt", "_DebugGCC.txt", "_ReleaseGCC.txt" };
        readonly static string[] ArchNames = 
        {
            //GTX 970M
            "arch=compute_52,code=sm_52",

            //GTX 1060, 1070
            "arch=compute_61,code=sm_61",

            //V100
            "arch=compute_70,code=sm_70",
        };

        public void WritteTheFile(string sSolDir, CProjFile projFile, Dictionary<string, CProjFile> excutables)
        {
            string sContent = "cmake_minimum_required(VERSION 3.8 FATAL_ERROR)\n\n";

            sContent += "project(CLG LANGUAGES CXX CUDA)\n\n";
            sContent += "set(CMAKE_GENERATOR_PLATFORM x64)\n\n";

            sContent += "# We start from CMAKE_SOURCE_DIR which should be /Code/CMake\n";
            sContent += "set(CMAKE_CURRENT_BINARY_DIR ${CMAKE_SOURCE_DIR}/../../Bin/UbuntuDebug)\n";
            sContent += "set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_CURRENT_BINARY_DIR})\n";
            sContent += "set(LIBRARY_OUTPUT_PATH  ${CMAKE_CURRENT_BINARY_DIR})\n";

            sContent += "# This is our code file dir\n";
            sContent += "set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/..)\n";

            sContent += "# Flags\n";
            sContent += "set(CMAKE_CUDA_FLAGS \"${CMAKE_CUDA_FLAGS} -O3\")\n";
            sContent += "set(CMAKE_CXX_FLAGS \"${CMAKE_CXX_FLAGS} -Ofast -Wall -Wno-unknown-pragmas -Wno-strict-overflow\")\n";
            sContent += "add_definitions(-D_UBUNTU)\n";
            sContent += "# to enable double float, add the following line:\n";
            if (m_bDouble)
            {
                sContent += "add_definitions(-D_CLG_DOUBLEFLOAT=1)\n";
                sContent += string.Format("MESSAGE(\"Note: double float is enabled, arch is {0}.\")\n", ArchNames[(int)m_eArch]);
            }
            else
            {
                sContent += "add_definitions(-D_CLG_DOUBLEFLOAT=0)\n";
                sContent += string.Format("MESSAGE(\"Note: double float NOT is enabled, arch is {0}.\")\n", ArchNames[(int)m_eArch]);
            }
            
            sContent += "MESSAGE(\"CMAKE_CUDA_FLAGS flag = ${CMAKE_CUDA_FLAGS}\")\n";
            sContent += "MESSAGE(\"CMAKE_CXX_FLAGS flag = ${CMAKE_CXX_FLAGS}\")\n\n";

            #region Add CLGLib

            sContent += "include_directories(${PROJECT_SOURCE_DIR}/CLGLib)\n";

            sContent += "add_library(CLGLib STATIC\n    ";
            foreach (string sFileName in projFile.m_lstAllHeaderFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/CLGLib/" + sFileName + "\n    ";
            }
            foreach (string sFileName in projFile.m_lstAllCuFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/CLGLib/" + sFileName + "\n    ";
            }
            foreach (string sFileName in projFile.m_lstAllCppFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/CLGLib/" + sFileName + "\n    ";
            }
            sContent += ")\n\n";
 
            sContent += @"# Request that CLGLib be built with -std=c++14
# As this is a public compile feature anything that links to 
# CLGLib will also build with -std=c++14
target_compile_features(CLGLib PUBLIC cxx_std_14)
 
# We need to explicitly state that we need all CUDA files in the 
# CLGLib library to be built with -dc as the member functions 
# could be called by other libraries and executables
set_target_properties( CLGLib
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)";

            sContent += "\n\ntarget_link_libraries(CLGLib -lcurand)\n\n";


            sContent += "# To enable the double, the minimum arch is 6.0\n";
            sContent += string.Format("target_compile_options(CLGLib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode {0}>)\n\n", ArchNames[(int)m_eArch]);

            #endregion

            #region Add CLGTest

            CProjFile slgTest = excutables["CLGTest"];

            sContent += "\n\n\n# ==================== \n# CLGTest \n# =================\n\n";

            sContent += "include_directories(${PROJECT_SOURCE_DIR}/CLGTest)\n";

            sContent += "add_executable(CLGTest \n    ";
            foreach (string sFileName in slgTest.m_lstAllHeaderFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/CLGTest/" + sFileName + "\n    ";
            }
            foreach (string sFileName in slgTest.m_lstAllCppFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/CLGTest/" + sFileName + "\n    ";
            }
            sContent += ")\n\n";

            sContent += "target_compile_features(CLGTest PUBLIC cxx_std_14)\n";
            sContent += "target_link_libraries(CLGTest CLGLib)\n";

            #endregion

            #region Add Application/RotatingReproduce

            CProjFile rotatingProj = excutables["RotatingReproduce"];

            sContent += "\n\n\n# ==================== \n# RotatingReproduce \n# =================\n\n";

            sContent += "include_directories(${PROJECT_SOURCE_DIR}/Applications/RotatingReproduce)\n";
            //sContent += "add_subdirectory(${PROJECT_SOURCE_DIR}/CLGTest)\n\n";

            sContent += "add_executable(RotatingReproduce \n    ";
            foreach (string sFileName in rotatingProj.m_lstAllHeaderFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/Applications/RotatingReproduce/" + sFileName + "\n    ";
            }
            foreach (string sFileName in rotatingProj.m_lstAllCppFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/Applications/RotatingReproduce/" + sFileName + "\n    ";
            }
            sContent += ")\n\n";

            sContent += "target_compile_features(RotatingReproduce PUBLIC cxx_std_14)\n";
            sContent += "target_link_libraries(RotatingReproduce CLGLib)\n";

            #endregion

            sContent = sContent.Replace("\r\n", "\n");
            sContent = sContent.Replace("\n\r", "\n");
            sContent = sContent.Replace("\r", "\n");
            sContent = sContent.Replace("\\", "/");

            File.WriteAllText(sSolDir + "CMake/CMakeLists.txt", sContent);

            //File.WriteAllText(sSolDir +"CMake/CMakeLists" 
            //    + (m_bWinOrUbuntu ? (m_bDebug ? FileSurfix[0] : FileSurfix[1]) : (m_bDebug ? FileSurfix[2] : FileSurfix[3])), sContent);1
        }
    }
}