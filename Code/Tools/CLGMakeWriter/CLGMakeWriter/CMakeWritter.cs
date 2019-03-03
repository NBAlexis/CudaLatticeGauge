
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

            #region Add CLGLib

            sContent += "project(CLGLib LANGUAGES CXX CUDA)\n\n";

            sContent += "set(CMAKE_GENERATOR_PLATFORM x64)\n\n";

            sContent += "# We start from CMAKE_SOURCE_DIR which should be /Code/CMake";
            sContent += "# First, we change it to /Code/CLGLib\n";
            sContent += "set(CMAKE_SOURCE_DIR  ${CMAKE_SOURCE_DIR}/../CLGLib)\n";

            sContent += "set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/../../Bin/UbuntuDebug)\n";
            sContent += "set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_BINARY_DIR})\n";
            sContent += "set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR})\n";
            sContent += "set(CMAKE_CURRENT_BINARY_DIR  ${CMAKE_BINARY_DIR})\n\n";
            
            sContent += "# This is our code file dir\n";
            sContent += "set(${PROJECT_SOURCE_DIR} ${CMAKE_SOURCE_DIR})\n";
            sContent += "include_directories(\"${PROJECT_SOURCE_DIR}\")\n\n";

            sContent += "message(\"CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}\")\n";
            sContent += "message(\"CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}\")\n";
            sContent += "message(\"CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}\")\n\n";


            sContent += "include_directories(${CMAKE_SOURCE_DIR})\n\n";

            sContent += "add_definitions(-D_UBUNTU)\n\n";

            sContent += "add_library(CLGLib STATIC\n    ";
            foreach (string sFileName in projFile.m_lstAllHeaderFiles)
            {
                sContent += "${CMAKE_SOURCE_DIR}/" + sFileName + "\n    ";
            }
            foreach (string sFileName in projFile.m_lstAllCuFiles)
            {
                sContent += "${CMAKE_SOURCE_DIR}/" + sFileName + "\n    ";
            }
            foreach (string sFileName in projFile.m_lstAllCppFiles)
            {
                sContent += "${CMAKE_SOURCE_DIR}/" + sFileName + "\n    ";
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

            sContent += "\n\ntarget_link_libraries(CLGLib libcurand.so libcudadevrt.a libcudart.so)\n\n";

            #endregion

            #region Add CLGTest

            CProjFile slgTest = excutables["CLGTest"];

            sContent += "\n\n\n# ==================== \n# CLGTest \n# =================\n\nproject(CLGTest LANGUAGES CXX)\n\n";

            sContent += "set(CMAKE_GENERATOR_PLATFORM x64)\n\n";

            sContent += "# We start from CMAKE_SOURCE_DIR which should be /Code/CMake";
            sContent += "# First, we change it to /Code/CLGTest\n";
            sContent += "set(CMAKE_SOURCE_DIR  ${CMAKE_SOURCE_DIR}/../CLGTest)\n";

            sContent += "set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/../../Bin/UbuntuDebug)\n";
            sContent += "set(EXECUTABLE_OUTPUT_PATH  ${CMAKE_BINARY_DIR})\n";
            sContent += "set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR})\n";
            sContent += "set(CMAKE_CURRENT_BINARY_DIR  ${CMAKE_BINARY_DIR})\n\n";

            sContent += "# This is our code file dir\n";
            sContent += "set(${PROJECT_SOURCE_DIR} ${CMAKE_SOURCE_DIR})\n";
            sContent += "include_directories(\"${PROJECT_SOURCE_DIR}\")\n\n";

            sContent += "message(\"CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}\")\n";
            sContent += "message(\"CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}\")\n";
            sContent += "message(\"CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}\")\n\n";


            sContent += "include_directories(${CMAKE_SOURCE_DIR}/../CLGLib)\n";
            sContent += "include_directories(${CMAKE_SOURCE_DIR})\n\n";

            sContent += "add_definitions(-D_UBUNTU)\n\n";

            sContent += "add_executable(CLGTest \n    ";
            foreach (string sFileName in slgTest.m_lstAllHeaderFiles)
            {
                sContent += "${CMAKE_SOURCE_DIR}/" + sFileName + "\n    ";
            }
            foreach (string sFileName in slgTest.m_lstAllCppFiles)
            {
                sContent += "${CMAKE_SOURCE_DIR}/" + sFileName + "\n    ";
            }
            sContent += ")\n\n";

            sContent += "target_compile_features(CLGTest PUBLIC cxx_std_14)\n";
            sContent += "target_link_libraries(CLGTest CLGLib)\n";

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