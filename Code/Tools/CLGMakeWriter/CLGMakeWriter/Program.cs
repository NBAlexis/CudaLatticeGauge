using System;
using System.Collections.Generic;
using System.IO;


namespace CLGMakeWriter
{
    /// <summary>
    /// This is for prepare to make on Linux
    /// This program should write CMakelist.txt, or even make with Bazel
    /// </summary>
    class Program
    {
        static void Main(string[] args)
        {
            string sSolFileName = System.AppDomain.CurrentDomain.BaseDirectory + "../../../../../../CudaLatticeGauge.sln";
            //string textSol;


            string projSolPath = Path.Combine(new[]{ System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../"});
            string projSolFilePath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../CudaLatticeGauge.sln" });

            string projCLGLibPath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../CLGLib" });
            string projCLGLibFilePath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../CLGLib/CLGLib.vcxproj" });
            string projCLGLibFilterPath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../CLGLib/CLGLib.vcxproj.filters" });

            string projCLGTestPath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../CLGTest" });
            string projCLGTestFilePath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../CLGTest/CLGTest.vcxproj" });
            string projCLGTestFilterPath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../CLGTest/CLGTest.vcxproj.filters" });

            string projRotatingPath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../Applications/RotatingReproduce" });
            string projRotatingFilePath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../Applications/RotatingReproduce/RotatingReproduce.vcxproj" });
            //string projRotatingFilterPath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../Applications/RotatingReproduce/CLGTest.vcxproj.filters" });

            string projMatchingRhoPath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../Applications/MatchingRho" });
            string projMatchingRhoFilePath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../Applications/MatchingRho/MatchingRho.vcxproj" });

            string projConfigurationCompresserPath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../Applications/ConfigurationCompresser" });
            string projConfigurationCompresserFilePath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../Applications/ConfigurationCompresser/ConfigurationCompresser.vcxproj" });

            if (!File.Exists(projSolFilePath)
             || !File.Exists(projCLGLibFilePath)
             || !File.Exists(projCLGLibFilterPath)
             || !File.Exists(projCLGTestFilePath)
             || !File.Exists(projCLGTestFilterPath))
            {
                Console.WriteLine("Failed...path incorrect...");
                return;
            }

            if (!File.Exists(projRotatingFilePath))
            {
                Console.WriteLine("Failed...cannot find rotating project...");
                return;
            }


            List<string> applicationProjFiles = new List<string>();
            List<string> applicationProjDirs = new List<string>();

            //gather applications

            CProjFile clgLibProj = new CProjFile(projCLGLibFilePath, projCLGLibPath);
            int[] versionNumber = clgLibProj.FindVersionNumber();
            Console.WriteLine(versionNumber[0].ToString() + "." + versionNumber[1].ToString());

            CMakeWritter writer = new CMakeWritter();
            Dictionary<string, CProjFile> apps = new Dictionary<string, CProjFile>();
            apps.Add("CLGTest", new CProjFile(projCLGTestFilePath, projCLGTestPath));
            apps.Add("RotatingReproduce", new CProjFile(projRotatingFilePath, projRotatingPath));
            apps.Add("MatchingRho", new CProjFile(projMatchingRhoFilePath, projMatchingRhoPath));
            apps.Add("ConfigurationCompresser", new CProjFile(projConfigurationCompresserFilePath, projConfigurationCompresserPath));
            
            writer.WritteTheFile(projSolPath, clgLibProj, apps);

            Console.WriteLine("work done, press enter to exit...");
            string byebye = Console.ReadLine();
        }
    }
}
