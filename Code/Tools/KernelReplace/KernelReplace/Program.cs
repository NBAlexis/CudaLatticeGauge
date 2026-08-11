using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using static System.Net.Mime.MediaTypeNames;

namespace KernelReplace
{
    internal class Program
    {
        static void Main(string[] args)
        {
            string sCodeFolder = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../Code/" });
            string sBinFolder = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, "../../../../../../Bin/" });

            List<string> sProblemFiles = new List<string>();

            //All code file
            DirectoryInfo codeFolder = new DirectoryInfo(sCodeFolder);
            FileInfo[] allheader = codeFolder.GetFiles("*.cu", SearchOption.AllDirectories);
            for (int i = 0; i < allheader.Length; ++i)
            {
                //string sFileName = allheader[i].FullName;
                //string sText = File.ReadAllText(sFileName);
                Console.WriteLine(allheader[i]);

                string sFileName = allheader[i].FullName;
                string sText = File.ReadAllText(sFileName);

                string sResult = Regex.Replace(sText, 
                    @"([a-zA-Z_][\w]*)[\s]*<[\s]*<[\s]*<[\s]*([a-zA-Z\d_]+)[\s]*,[\s]*([a-zA-Z\d_]+)[\s]*>[\s]*>[\s]*>[\s]*\(",
                    "_LAUNCH_KERNEL($1, $2, $3, ");

                Console.WriteLine(sResult);

                //File.WriteAllText(sFileName, sResult);
            }
        }
    }
}
