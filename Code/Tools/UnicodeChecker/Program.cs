using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;

namespace UnicodeChecker
{
    class Program
    {
        static void Main(string[] args)
        {
            string sCodeFolder = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../Code/" });
            string sBinFolder = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../Bin/" });

            List<string> sProblemFiles = new List<string>();

            //All code file
            DirectoryInfo codeFolder = new DirectoryInfo(sCodeFolder);
            FileInfo[] allheader = codeFolder.GetFiles("*.h", SearchOption.AllDirectories);
            for (int i = 0; i < allheader.Length; ++i)
            {
                string sFileName = allheader[i].FullName;
                string sText = File.ReadAllText(sFileName);
                if (CheckExistUnicode(sText))
                {
                    Console.WriteLine("==== {0} is NOT OK", allheader[i].Name);
                    sProblemFiles.Add(sFileName);
                }
                else
                {
                    Console.WriteLine("==== {0} is OK", allheader[i].Name);
                }
            }

            allheader = codeFolder.GetFiles("*.cpp", SearchOption.AllDirectories);
            for (int i = 0; i < allheader.Length; ++i)
            {
                string sFileName = allheader[i].FullName;
                string sText = File.ReadAllText(sFileName);
                if (CheckExistUnicode(sText))
                {
                    Console.WriteLine("==== {0} is NOT OK", allheader[i].Name);
                    sProblemFiles.Add(sFileName);
                }
                else
                {
                    Console.WriteLine("==== {0} is OK", allheader[i].Name);
                }
            }

            allheader = codeFolder.GetFiles("*.cuh", SearchOption.AllDirectories);
            for (int i = 0; i < allheader.Length; ++i)
            {
                string sFileName = allheader[i].FullName;
                string sText = File.ReadAllText(sFileName);
                if (CheckExistUnicode(sText))
                {
                    Console.WriteLine("==== {0} is NOT OK", allheader[i].Name);
                    sProblemFiles.Add(sFileName);
                }
                else
                {
                    Console.WriteLine("==== {0} is OK", allheader[i].Name);
                }
            }

            allheader = codeFolder.GetFiles("*.cu", SearchOption.AllDirectories);
            for (int i = 0; i < allheader.Length; ++i)
            {
                string sFileName = allheader[i].FullName;
                string sText = File.ReadAllText(sFileName);
                if (CheckExistUnicode(sText))
                {
                    Console.WriteLine("==== {0} is NOT OK", allheader[i].Name);
                    sProblemFiles.Add(sFileName);
                }
                else
                {
                    Console.WriteLine("==== {0} is OK", allheader[i].Name);
                }
            }

            allheader = codeFolder.GetFiles("*.c", SearchOption.AllDirectories);
            for (int i = 0; i < allheader.Length; ++i)
            {
                string sFileName = allheader[i].FullName;
                string sText = File.ReadAllText(sFileName);
                if (CheckExistUnicode(sText))
                {
                    Console.WriteLine("==== {0} is NOT OK", allheader[i].Name);
                    sProblemFiles.Add(sFileName);
                }
                else
                {
                    Console.WriteLine("==== {0} is OK", allheader[i].Name);
                }
            }

            DirectoryInfo binFolder = new DirectoryInfo(sBinFolder);
            allheader = binFolder.GetFiles("*.yaml", SearchOption.AllDirectories);
            for (int i = 0; i < allheader.Length; ++i)
            {
                string sFileName = allheader[i].FullName;
                string sText = File.ReadAllText(sFileName);
                if (CheckExistUnicode(sText))
                {
                    Console.WriteLine("==== {0} is NOT OK", allheader[i].Name);
                    sProblemFiles.Add(sFileName);
                }
                else
                {
                    Console.WriteLine("==== {0} is OK", allheader[i].Name);
                }
            }

            if (sProblemFiles.Count > 0)
            {
                for (int i = 0; i < sProblemFiles.Count; ++i)
                {
                    Console.WriteLine("{0} is NOT OK", sProblemFiles[i]);
                    ShowLineNumberOfUnicode(sProblemFiles[i]);
                }
            }
            else
            {
                Console.WriteLine("==== All files are OK");
            }

            //Console.WriteLine(sCodeFolder);
            //Console.WriteLine(sBinFolder);


            Console.WriteLine("work done, press enter to exit...");
            string byebye = Console.ReadLine();
        }

        public static bool CheckExistUnicode(string strInput)
        {
            int i = strInput.Length;
            if (i == 0)
                return false;

            int j = System.Text.Encoding.Default.GetBytes(strInput).Length;

            if (i != j)
                return true;
            else
                return false;
        }

        public static void ShowLineNumberOfUnicode(string sFileName)
        {
            string fileContent = File.ReadAllText(sFileName);
            string[] lines = fileContent.Split('\n');
            for (int i = 0; i < lines.Length; ++i)
            {
                if (CheckExistUnicode(lines[i]))
                {
                    Console.WriteLine("Line{0}: {1}", i + 1, lines[i]);
                }
            }
        }
    }
}
