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
        public bool m_bDouble = true;
        public EArch m_eArch = EArch.EArcSM61;
        public bool m_bDebug = true;
        public bool m_bWinOrUbuntu = true;

        //turn it off, recently we will not use this
        public bool m_bHasCompresser = false;

        //turn it off, we have sign problem here
        public bool m_bHasConstAcc = false;
        
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
            string sContent = "cmake_minimum_required(VERSION 3.8 FATAL_ERROR)\n\n";

            sContent += "if (DEFINED NVCCROOT)\n";
            sContent += "    set(CMAKE_CUDA_COMPILER ${NVCCROOT})\n";
            sContent += "    MESSAGE(\"CMAKE_CUDA_COMPILER = ${CMAKE_CUDA_COMPILER}\")\n";
            sContent += "endif()\n\n";

            sContent += string.Format("set(CUDA_CMP \"{0}\")\n", ArchNames[(int)m_eArch]);
            sContent += string.Format("set(CUDA_SM \"{0}\")\n", CodeNames[(int)m_eArch]);

            sContent += "if (DEFINED CUDASM)\n";
            sContent += "    set(CUDA_CMP \"compute_${CUDASM}\")\n";
            sContent += "    set(CUDA_SM \"sm_${CUDASM}\")\n";
            sContent += "endif()\n\n";

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
            }
            else
            {
                sContent += "add_definitions(-D_CLG_DOUBLEFLOAT=0)\n";
            }
            sContent += "MESSAGE(\"Note: double float is enabled, arch is ${CUDA_CMP} and ${CUDA_SM}.\")\n";

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

            sContent += "\n\ntarget_link_libraries(CLGLib -lcurand)\n";
            sContent += "target_link_libraries(CLGLib -lcufft)\n";

            sContent += "\n# To enable the double, the minimum arch is 6.0\n";
            sContent += "target_compile_options(CLGLib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=${CUDA_CMP},code=${CUDA_SM}>)\n\n";

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

            #region Add Applications

            sContent += AddApplication(excutables["RotatingReproduce"]);
            sContent += AddApplication(excutables["MatchingRho"]);

            if (m_bHasCompresser)
            {
                sContent += AddApplication(excutables["ConfigurationCompresser"]);
            }

            if (m_bHasConstAcc)
            {
                sContent += AddApplication(excutables["ConstAcc"]);
            }

            sContent += AddApplication(excutables["StaggeredSpectrum"]);
            sContent += AddApplication(excutables["StaggeredRotation"]);

            #endregion

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
            sRet += string.Format("include_directories({1}/Applications/{0})\n", sAppName, "${PROJECT_SOURCE_DIR}");
            sRet += string.Format("add_executable({0} \n    ", sAppName);
            foreach (string sFileName in addProj.m_lstAllHeaderFiles)
            {
                sRet += string.Format("{1}/Applications/{0}/" + sFileName + "\n    ", sAppName, "${PROJECT_SOURCE_DIR}");
            }

            foreach (string sFileName in addProj.m_lstAllCppFiles)
            {
                sRet += string.Format("{1}/Applications/{0}/" + sFileName + "\n    ", sAppName, "${PROJECT_SOURCE_DIR}");
            }

            sRet += ")\n\n";

            sRet += string.Format("target_compile_features({0} PUBLIC cxx_std_14)\n", sAppName);
            sRet += string.Format("target_link_libraries({0} CLGLib)\n\n", sAppName);
            return sRet;
        }
    }
}