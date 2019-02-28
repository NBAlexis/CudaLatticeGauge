
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
            string sContent = "cmake_minimum_required (VERSION 2.6)\n";
            
            sContent += "message(\"CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}\")\n";
            sContent += "message(\"CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}\")\n";
            sContent += "message(\"CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}\")\n";

            sContent += "project(CLGLib)\n";

            File.WriteAllText(sSolDir +"CMake/CMakeLists" 
                + (m_bWinOrUbuntu ? (m_bDebug ? FileSurfix[0] : FileSurfix[1]) : (m_bDebug ? FileSurfix[2] : FileSurfix[3])), sContent);
        }
    }
}