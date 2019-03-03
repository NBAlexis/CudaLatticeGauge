
using System.Collections.Generic;
using System.IO;

namespace CLGMakeWriter
{
    class CMakeWritter
    {
        public bool m_bDebug = true;
        public bool m_bWinOrUbuntu = true;
        readonly static string[] FileSurfix = { "_DebugMSVC.txt", "_ReleaseMSVC.txt", "_DebugGCC.txt", "_ReleaseGCC.txt" };

        public void WritteTheFile(string sSolDir, CProjFile projFile, Dictionary<string, CProjFile> excutables)
        {
            string sContent = "cmake_minimum_required(VERSION 3.8 FATAL_ERROR)\n\n";

            sContent += "project(CLG LANGUAGES CXX CUDA)\n\n";
            sContent += "set(CMAKE_GENERATOR_PLATFORM x64)\n\n";

            sContent += "# We start from CMAKE_SOURCE_DIR which should be /Code/";
            sContent += "set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/../Bin/UbuntuDebug)\n";
            sContent += "set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_BINARY_DIR})\n";
            sContent += "set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR})\n";
            sContent += "set(CMAKE_CURRENT_BINARY_DIR  ${CMAKE_BINARY_DIR})\n\n";

            sContent += "# This is our code file dir\n";
            sContent += "set(${PROJECT_SOURCE_DIR} ${CMAKE_SOURCE_DIR})\n";

            sContent += "add_definitions(-D_UBUNTU)\n\n";

            #region Add CLGLib

            sContent += "include_directories(${PROJECT_SOURCE_DIR}/CLGLib)\n";
            sContent += "add_subdirectory(${PROJECT_SOURCE_DIR}/CLGLib)\n\n";

            sContent += "add_library(cudadevrt STATIC IMPORTED)\n";
            sContent += "add_library(cudart_static STATIC IMPORTED)\n";
            sContent += "add_library(cudrand_static STATIC IMPORTED)\n\n";

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

            sContent += "\n\ntarget_link_libraries(CLGLib cudadevrt cudart_static curand_static)\n\n";

            #endregion

            #region Add CLGTest

            CProjFile slgTest = excutables["CLGTest"];

            sContent += "\n\n\n# ==================== \n# CLGTest \n# =================\n\n";

            sContent += "include_directories(${PROJECT_SOURCE_DIR}/CLGTest)\n";
            sContent += "add_subdirectory(${PROJECT_SOURCE_DIR}/CLGTest)\n\n";

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

            sContent = sContent.Replace("\r\n", "\n");
            sContent = sContent.Replace("\n\r", "\n");
            sContent = sContent.Replace("\r", "\n");
            sContent = sContent.Replace("\\", "/");

            File.WriteAllText(sSolDir + "CMakeLists.txt", sContent);

            //File.WriteAllText(sSolDir +"CMake/CMakeLists" 
            //    + (m_bWinOrUbuntu ? (m_bDebug ? FileSurfix[0] : FileSurfix[1]) : (m_bDebug ? FileSurfix[2] : FileSurfix[3])), sContent);1
        }
    }
}