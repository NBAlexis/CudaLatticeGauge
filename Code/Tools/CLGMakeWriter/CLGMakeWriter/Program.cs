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

            if (!File.Exists(projSolFilePath)
             || !File.Exists(projCLGLibFilePath)
             || !File.Exists(projCLGLibFilterPath)
             || !File.Exists(projCLGTestFilePath)
             || !File.Exists(projCLGTestFilterPath))
            {
                Console.WriteLine("Failed...path incorrect...");
                return;
            }

            //gather applications

            CProjFile clgLibProj = new CProjFile("CLGLib", projCLGLibFilePath, projCLGLibPath);
            int[] versionNumber = clgLibProj.FindVersionNumber();
            Console.WriteLine(versionNumber[0].ToString() + "." + versionNumber[1].ToString());

            CMakeWritter writer = new CMakeWritter();
            Dictionary<string, CProjFile> apps = new Dictionary<string, CProjFile>();
            apps.Add("CLGTest", new CProjFile("CLGTest", projCLGTestFilePath, projCLGTestPath));
            apps.Add("RotatingReproduce", new CProjFile("RotatingReproduce"));
            apps.Add("MatchingRho", new CProjFile("MatchingRho"));
            apps.Add("ConfigurationCompresser", new CProjFile("ConfigurationCompresser"));
            apps.Add("ConstAcc", new CProjFile("ConstAcc"));
            apps.Add("StaggeredSpectrum", new CProjFile("StaggeredSpectrum"));
            apps.Add("StaggeredRotation", new CProjFile("StaggeredRotation"));

            writer.WritteTheFile(projSolPath, clgLibProj, apps);

            Console.WriteLine("work done, press enter to exit...");
            string byebye = Console.ReadLine();
        }
    }
}
