using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection.Metadata;
using System.Runtime.InteropServices;
using System.Text.RegularExpressions;

namespace Precommit
{
    class Program
    {
        static string logfile = "pre-commit.log";

        static int Main(string[] args)
        {
            //Console.WriteLine("Current path dir: " + Directory.GetCurrentDirectory());
            //Console.WriteLine("Current path env: " + Environment.CurrentDirectory);
            //Console.WriteLine("Current path dom: " + AppDomain.CurrentDomain.BaseDirectory);
            string textContent = File.ReadAllText(AppDomain.CurrentDomain.BaseDirectory + "../../Code/CLGLib/Core/Version.cpp");
            int version = int.Parse(Regex.Match(textContent, @"(?<a>appVersion\(\)[\s]*\{[\s]*return[\s]+)([\d]+)(?<b>;[\s]*\})").Groups[1].Value) + 1;
            textContent = Regex.Replace(textContent, @"(?<a>appVersion\(\)[\s]*\{[\s]*return[\s]+)([\d]+)(?<b>;[\s]*\})", "${a}" + version.ToString() + "${b}");
            File.WriteAllText(AppDomain.CurrentDomain.BaseDirectory + "../../Code/CLGLib/Core/Version.cpp", textContent);

            Console.WriteLine("Time: " + DateTime.Now.ToString() + " Version: " + version.ToString());

            return 0;
        }

    }
}
