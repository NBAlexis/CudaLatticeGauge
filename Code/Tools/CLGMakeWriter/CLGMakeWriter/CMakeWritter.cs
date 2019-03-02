
using System.IO;

namespace CLGMakeWriter
{
    class CMakeWritter
    {
        public bool m_bDebug = true;
        public bool m_bWinOrUbuntu = true;
        readonly static string[] FileSurfix = { "_DebugMSVC.txt", "_ReleaseMSVC.txt", "_DebugGCC.txt", "_ReleaseGCC.txt" };

        public void WritteTheFile(string sSolDir, CProjFile projFile)
        {
            string sContent = "cmake_minimum_required(VERSION 3.8 FATAL_ERROR)\n\n";

            //CMAKE_BINARY_DIR, CMAKE_SOURCE_DIR, CMAKE_CURRENT_BINARY_DIR all = /Code/CMake
            sContent += "message(\"CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}\")\n";
            sContent += "message(\"CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}\")\n";
            sContent += "message(\"CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}\")\n\n";

            sContent += "project(CLGLib LANGUAGES CXX CUDA)\n\n";
            sContent += "include(CLGLib)\n\n";
            sContent += "add_library(CLGLib STATIC\n    ";
            foreach (string sFileName in projFile.m_lstAllHeaderFiles)
            {
                sContent += sFileName + "\n    ";
            }
            foreach (string sFileName in projFile.m_lstAllCuFiles)
            {
                sContent += sFileName + "\n    ";
            }
            foreach (string sFileName in projFile.m_lstAllCppFiles)
            {
                sContent += sFileName + "\n    ";
            }
            sContent += ")\n\n";
 
            sContent += @"# Request that CLGLib be built with -std=c++11
# As this is a public compile feature anything that links to 
# CLGLib will also build with -std=c++11
target_compile_features(CLGLib PUBLIC cxx_std_11)
 
# We need to explicitly state that we need all CUDA files in the 
# CLGLib library to be built with -dc as the member functions 
# could be called by other libraries and executables
set_target_properties( CLGLib
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)";

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