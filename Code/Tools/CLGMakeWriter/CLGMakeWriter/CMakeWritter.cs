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
                sContent += string.Format("MESSAGE(\"Note: double float is enabled, arch is {0} and {1}.\")\n", ArchNames[(int)m_eArch], CodeNames[(int)m_eArch]);
            }
            else
            {
                sContent += "add_definitions(-D_CLG_DOUBLEFLOAT=0)\n";
                sContent += string.Format("MESSAGE(\"Note: double float NOT is enabled, arch is {0} and {1}.\")\n", ArchNames[(int)m_eArch], CodeNames[(int)m_eArch]);
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

            sContent += "\n\ntarget_link_libraries(CLGLib -lcurand)\n";
            sContent += "target_link_libraries(CLGLib -lcufft)\n";

            sContent += "\n# To enable the double, the minimum arch is 6.0\n";
            sContent += string.Format("target_compile_options(CLGLib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch={0},code={1}>)\n\n", ArchNames[(int)m_eArch], CodeNames[(int)m_eArch]);

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

            #region Add Matching

            CProjFile matchingProj = excutables["MatchingRho"];

            sContent += "\n\n\n# ==================== \n# MatchingRho \n# =================\n\n";

            sContent += "include_directories(${PROJECT_SOURCE_DIR}/Applications/MatchingRho)\n";
            //sContent += "add_subdirectory(${PROJECT_SOURCE_DIR}/CLGTest)\n\n";

            sContent += "add_executable(MatchingRho \n    ";
            foreach (string sFileName in matchingProj.m_lstAllHeaderFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/Applications/MatchingRho/" + sFileName + "\n    ";
            }
            foreach (string sFileName in matchingProj.m_lstAllCppFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/Applications/MatchingRho/" + sFileName + "\n    ";
            }
            sContent += ")\n\n";

            sContent += "target_compile_features(MatchingRho PUBLIC cxx_std_14)\n";
            sContent += "target_link_libraries(MatchingRho CLGLib)\n";

            #endregion

            #region Add Compresser

            if (m_bHasCompresser)
            {
                CProjFile compresserProj = excutables["ConfigurationCompresser"];

                sContent += "\n\n\n# ==================== \n# ConfigurationCompresser \n# =================\n\n";
                sContent += "include_directories(${PROJECT_SOURCE_DIR}/Applications/ConfigurationCompresser)\n";
                sContent += "add_executable(ConfigurationCompresser \n    ";
                foreach (string sFileName in compresserProj.m_lstAllHeaderFiles)
                {
                    sContent += "${PROJECT_SOURCE_DIR}/Applications/ConfigurationCompresser/" + sFileName + "\n    ";
                }

                foreach (string sFileName in compresserProj.m_lstAllCppFiles)
                {
                    sContent += "${PROJECT_SOURCE_DIR}/Applications/ConfigurationCompresser/" + sFileName + "\n    ";
                }

                sContent += ")\n\n";

                sContent += "target_compile_features(ConfigurationCompresser PUBLIC cxx_std_14)\n";
                sContent += "target_link_libraries(ConfigurationCompresser CLGLib)\n";
            }

            #endregion

            #region Add ConstAcc

            if (m_bHasConstAcc)
            {
                CProjFile constaccProj = excutables["ConstAcc"];

                sContent += "\n\n\n# ==================== \n# ConstAcc \n# =================\n\n";
                sContent += "include_directories(${PROJECT_SOURCE_DIR}/Applications/ConstAcc)\n";
                sContent += "add_executable(ConstAcc \n    ";
                foreach (string sFileName in constaccProj.m_lstAllHeaderFiles)
                {
                    sContent += "${PROJECT_SOURCE_DIR}/Applications/ConstAcc/" + sFileName + "\n    ";
                }
                foreach (string sFileName in constaccProj.m_lstAllCppFiles)
                {
                    sContent += "${PROJECT_SOURCE_DIR}/Applications/ConstAcc/" + sFileName + "\n    ";
                }
                sContent += ")\n\n";

                sContent += "target_compile_features(ConstAcc PUBLIC cxx_std_14)\n";
                sContent += "target_link_libraries(ConstAcc CLGLib)\n";
            }

            #endregion

            #region Staggered Spectrum

            CProjFile staggeredSpectrumProj = excutables["StaggeredSpectrum"];

            sContent += "\n\n\n# ==================== \n# Staggered Spectrum \n# =================\n\n";
            sContent += "include_directories(${PROJECT_SOURCE_DIR}/Applications/StaggeredSpectrum)\n";
            sContent += "add_executable(StaggeredSpectrum \n    ";
            foreach (string sFileName in staggeredSpectrumProj.m_lstAllHeaderFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/Applications/StaggeredSpectrum/" + sFileName + "\n    ";
            }
            foreach (string sFileName in staggeredSpectrumProj.m_lstAllCppFiles)
            {
                sContent += "${PROJECT_SOURCE_DIR}/Applications/StaggeredSpectrum/" + sFileName + "\n    ";
            }
            sContent += ")\n\n";

            sContent += "target_compile_features(StaggeredSpectrum PUBLIC cxx_std_14)\n";
            sContent += "target_link_libraries(StaggeredSpectrum CLGLib)\n";

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