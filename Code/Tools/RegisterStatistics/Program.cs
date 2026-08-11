using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Xml.Linq;

namespace RegisterStatistics
{
    public class CuobjdumpParser
    {

        public static string FindCuobjdumpPath()
        {
            var cudaPath = Environment.GetEnvironmentVariable("CUDA_PATH");
            if (string.IsNullOrEmpty(cudaPath))
            {
                throw new FileNotFoundException("CUDA_PATH not found");
            }
            var cuobjdumpPath = Path.Combine(cudaPath, "bin", "cuobjdump.exe");
            if (!File.Exists(cuobjdumpPath))
            {
                throw new FileNotFoundException($"cuobjdump.exe not found: {cudaPath}");
            }
            return cuobjdumpPath;
        }

        public static string RunCuobjdump(string cuobjdumpPath, string targetExePath)
        {
            if (!File.Exists(targetExePath))
            {
                throw new FileNotFoundException($"Target file not exist {targetExePath}");
            }

            var processInfo = new ProcessStartInfo
            {
                FileName = cuobjdumpPath,
                Arguments = $"--dump-resource-usage --gpu-architecture sm_86 \"{targetExePath}\"", 
                RedirectStandardOutput = true,
                RedirectStandardError = true,
                UseShellExecute = false,
                CreateNoWindow = true,
                StandardOutputEncoding = Encoding.UTF8,
                StandardErrorEncoding = Encoding.UTF8
            };

            var output = new StringBuilder();
            using (var process = new Process())
            {
                process.StartInfo = processInfo;
                process.OutputDataReceived += (sender, e) =>
                {
                    if (!string.IsNullOrEmpty(e.Data))
                    {
                        output.AppendLine(e.Data);
                    }
                };
                process.ErrorDataReceived += (sender, e) =>
                {
                    if (!string.IsNullOrEmpty(e.Data))
                    {
                        output.AppendLine($"[ERROR] {e.Data}");
                    }
                };

                process.Start();
                process.BeginOutputReadLine();
                process.BeginErrorReadLine();

                if (!process.WaitForExit(5000))
                {
                    throw new TimeoutException("cuobjdump timeout");
                }

                if (process.ExitCode != 0)
                {
                    throw new Exception($"cuobjdump return error code: {process.ExitCode}\n{output}");
                }

                return output.ToString();
            }
        }


        public static void ParseRegisterUsage(string cuobjdumpOutput, string sOutput)
        {
            const int MaxNameDisplayLength = 80;
            const int RegistersColumnWidth = 10;

            int totalWidth = 4 + MaxNameDisplayLength + RegistersColumnWidth + 2;
            string separatorLine = new string('=', totalWidth);
            string middleLine = new string('-', totalWidth);

            var kernelList = new List<(string Name, int Registers)>();

            var regex = new Regex(
                @"Function\s+([^:\r\n]+):[\s\S]*?REG:(\d+)",
                RegexOptions.Multiline | RegexOptions.Compiled
            );

            var matches = regex.Matches(cuobjdumpOutput);

            foreach (Match match in matches)
            {
                if (match.Success &&
                    match.Groups.Count > 2 &&
                    int.TryParse(match.Groups[2].Value, out int registers))
                {
                    string originalName = match.Groups[1].Value.Trim();
                    //Console.WriteLine(originalName);
                    string cleanedName = LightweightDemangler.Demangle(originalName); 
                    //Console.WriteLine(cleanedName);
                    kernelList.Add((cleanedName, registers));
                }
            }

            // desend
            var sortedKernels = kernelList
                .OrderByDescending(k => k.Registers)
                .ThenBy(k => k.Name)
                .ToList();

            string headerFormat = $"{{0,-4}} {{1,-{MaxNameDisplayLength}}} {{2}}";
            string rowFormat = $"{{0,-4}} {{1,-{MaxNameDisplayLength}}} {{2}}";

            var sb = new StringBuilder();
            sb.AppendLine("Kernel Register usages:");
            sb.AppendLine(separatorLine);
            sb.AppendLine(string.Format(headerFormat, "Rank", "Kernel Name", "Registers"));
            sb.AppendLine(middleLine);

            for (int i = 0; i < sortedKernels.Count; i++)
            {
                sb.AppendLine(string.Format(rowFormat,
                    $"[{i + 1}]",
                    TruncateString(sortedKernels[i].Name, MaxNameDisplayLength),
                    sortedKernels[i].Registers));
            }
            sb.AppendLine(separatorLine + "\n");

            Console.WriteLine(sb.ToString());

            File.WriteAllText(sOutput, sb.ToString() + "\n" + cuobjdumpOutput);
            OpenFileWithDefaultProgram(sOutput);
        }

        private static string TruncateString(string value, int maxLength)
        {
            if (string.IsNullOrEmpty(value)) return value;
            return value.Length <= maxLength ?
                   value :
                   value.Substring(0, maxLength - 3) + "...";
        }

        public static void Main(string[] args)
        {
            try
            {
                string sOutputFile = "../../../../../../Bin/Debug/CLGLib.dll.txt";
                string targetExe = Path.Combine(Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory,
                    "../../../../../../Bin/Debug/CLGLib_d.dll" }));
                if (Release)
                {
                    sOutputFile = "../../../../../../Bin/Release/CLGLib.dll.txt";
                    targetExe = Path.Combine(Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory,
                    "../../../../../../Bin/Release/CLGLib.dll" }));
                }

                //Console.WriteLine(Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../Bin/" }));

                string cuobjdumpPath = FindCuobjdumpPath();
                Console.WriteLine($"Found cuobjdump: {cuobjdumpPath}");

                string output = RunCuobjdump(cuobjdumpPath, targetExe);

                //Console.WriteLine(output);
                ParseRegisterUsage(output, sOutputFile);
            }
            catch (Exception ex)
            {
                Console.WriteLine($"Failed: {ex.Message}");
            }
        }


        private static void OpenFileWithDefaultProgram(string filePath)
        {
            try
            {
                Process.Start(new ProcessStartInfo
                {
                    FileName = filePath,
                    UseShellExecute = true 
                });
            }
            catch (Win32Exception ex) 
            {
                Console.WriteLine($"Default program failed {ex.NativeErrorCode}: {ex.Message}");

                try
                {
                    Process.Start(new ProcessStartInfo
                    {
                        FileName = "cmd.exe", 
                        Arguments = $"/c code \"{filePath}\"", 
                        WindowStyle = ProcessWindowStyle.Hidden, 
                        UseShellExecute = false
                    });
                }
                catch (Exception codeEx)
                {
                    Console.WriteLine($"VSCode faield: {codeEx.Message}");
                    if (codeEx is Win32Exception win32Ex)
                    {
                        Console.WriteLine($"error code: {win32Ex.NativeErrorCode}");
                        Console.WriteLine($"detailed reason: {win32Ex.Message}");
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine($"unknown error: {ex.Message}");
            }
        }

        static public bool Release = false;
    }
}