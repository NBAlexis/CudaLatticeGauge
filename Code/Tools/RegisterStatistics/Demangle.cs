using System;
using System.Collections.Generic;
using System.Text.RegularExpressions;

namespace RegisterStatistics
{

    public static class LightweightDemangler
    {
        static List<string> m_lstNameSpaces = new List<string>();

        public static string FindNameSpace(string s)
        {
            m_lstNameSpaces.Clear();

            Match match = Regex.Match(s, @"^_ZN(\d+)");

            if (match.Success)
            {
                string prefixNumber = match.Groups[1].Value;
                int numeral = int.Parse(prefixNumber);
                string sNSName = s.Substring(3 + prefixNumber.Length, numeral);
                if (!m_lstNameSpaces.Contains(sNSName))
                {
                    m_lstNameSpaces.Add(sNSName);
                }
                s = s.Substring(3 + prefixNumber.Length + numeral);
                return s;
            }

            return s;
        }

        public static string FindAllFunctionOrStructureNames(string s)
        {
            Regex regex = new Regex(@"^(\d+)(.+?)[IE]");
            Match match = regex.Match(s);  // 关键修改点

            try
            {
                if (match.Success)
                {
                    string startNumber = match.Groups[1].Value;
                    int numeral = int.Parse(startNumber);
                    if ("I" == s.Substring(startNumber.Length + numeral, 1))
                    {
                        string sFunctionName = s.Substring(startNumber.Length, numeral);
                        string sTemplate;
                        string sArgument;
                        if (ParseTemplate(s.Substring(startNumber.Length + numeral + 1), out sTemplate, out sArgument))
                        {
                            return sFunctionName + ReplaceTemplateArgs(sTemplate) + "(" + ParseArgument(sArgument) + ")";
                        }
                        return sFunctionName + s.Substring(startNumber.Length + numeral + 1);
                    }
                    else if ("E" == s.Substring(startNumber.Length + numeral, 1))
                    {
                        string sFunctionName = s.Substring(startNumber.Length, numeral);
                        return sFunctionName + "(" + ParseArgument(s.Substring(startNumber.Length + numeral + 1)) + ")";
                    }
                    return s;
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.StackTrace);
            }

            return s;
        }

        public static bool ParseTemplate(string s, out string result, out string arg)
        {
            result = "<";
            int depth = 1;
            while (depth > 0 && s.Length > 0)
            {
                if (s.StartsWith("NS_"))
                {
                    s = s.Substring(3);
                    continue;
                }

                if (s.StartsWith("I"))
                {
                    s = s.Substring(1);
                    result = result + "<";
                    depth = depth + 1;
                    continue;
                }

                if (s.StartsWith("EE"))
                {
                    s = s.Substring(2);
                    result = result + ">";
                    depth = depth - 1;
                    continue;
                }

                if (s.StartsWith("E"))
                {
                    s = s.Substring(1);
                    result = result + ", ";
                    continue;
                }

                result = result + s.Substring(0, 1);
                s = s.Substring(1);
            }

            if (0 == depth && s.StartsWith("Ev"))
            {
                arg = s.Substring(2);
                return true;
            }
            if (0 == depth && s.StartsWith("v"))
            {
                arg = s.Substring(1);
                return true;
            }
            arg = "";
            return false;
        }
        public static string ParseArgument(string s)
        {
            string result = "";
            bool inPointer = false;
            bool inConstantPointer = false;
            bool inConstant = false;
            bool inConstReference = false;
            List<string> parsed = new List<string>();
            while (s.Length > 0)
            {
                if (s.StartsWith("NS_"))
                {
                    s = s.Substring(3);
                    continue;
                }

                if (s[0] == 'i')
                {
                    result = result + "INT, ";
                    s = s.Substring(1);
                    parsed.Add("INT");
                    continue;
                }
                if (s[0] == 'h')
                {
                    result = result + "BYTE, ";
                    s = s.Substring(1);
                    parsed.Add("BYTE");
                    continue;
                }
                if (s[0] == 'd')
                {
                    result = result + "DOUBLE, ";
                    s = s.Substring(1);
                    parsed.Add("DOUBLE");
                    continue;
                }
                if (s[0] == 'j')
                {
                    result = result + "UINT, ";
                    s = s.Substring(1);
                    parsed.Add("UINT");
                    continue;
                }

                if (s.StartsWith("S0_")
                 || s.StartsWith("S1_")
                 || s.StartsWith("S2_")
                 || s.StartsWith("S3_")
                 || s.StartsWith("S4_")
                 || s.StartsWith("S5_"))
                {
                    s = s.Substring(3);
                    result = result + parsed[parsed.Count - 1] + ", ";
                    parsed.Add(parsed[parsed.Count - 1]);
                    continue;
                }

                if (s.StartsWith("T0_"))
                {
                    s = s.Substring(3);
                    result = result + "T0, ";
                    parsed.Add("T0");
                    continue;
                }

                if (s.StartsWith("T_"))
                {
                    s = s.Substring(2);
                    result = result + "T, ";
                    parsed.Add("T");
                    continue;
                }

                if (s.StartsWith("KT_"))
                {
                    s = s.Substring(3);
                    result = result + "const T, ";
                    parsed.Add("T");
                    continue;
                }
                if (s.StartsWith("PKT_"))
                {
                    s = s.Substring(4);
                    result = result + "const T*, ";
                    parsed.Add("T");
                    continue;
                }

                if (s.StartsWith("PKi"))
                {
                    s = s.Substring(3);
                    result = result + "const INT*, ";
                    parsed.Add("INT");
                    continue;
                }
                if (s.StartsWith("PKh"))
                {
                    s = s.Substring(3);
                    result = result + "const BYTE*, ";
                    parsed.Add("BYTE");
                    continue;
                }
                if (s.StartsWith("PKd"))
                {
                    s = s.Substring(3);
                    result = result + "const DOUBLE*, ";
                    parsed.Add("DOUBLE");
                    continue;
                }
                if (s.StartsWith("PKj"))
                {
                    s = s.Substring(3);
                    result = result + "const UINT*, ";
                    parsed.Add("UINT");
                    continue;
                }

                if (s.StartsWith("PKS0_")
                 || s.StartsWith("PKS1_")
                 || s.StartsWith("PKS2_")
                 || s.StartsWith("PKS3_")
                 || s.StartsWith("PKS4_"))
                {
                    s = s.Substring(5);
                    result = result + "const " + parsed[parsed.Count - 1] + "*, ";
                    parsed.Add(parsed[parsed.Count - 1]);
                    continue;
                }

                if (s.StartsWith("PS0_")
                 || s.StartsWith("PS1_")
                 || s.StartsWith("PS2_")
                 || s.StartsWith("PS3_")
                 || s.StartsWith("PS4_"))
                {
                    s = s.Substring(4);
                    result = result + parsed[parsed.Count - 1] + "*, ";
                    parsed.Add(parsed[parsed.Count - 1]);
                    continue;
                }

                if (s.StartsWith("RK"))
                {
                    if (!inPointer && !inConstantPointer && !inConstant && !inConstReference)
                    {
                        s = s.Substring(2);
                        inConstReference = true;
                        continue;
                    }
                }
                if (s.StartsWith("PK"))
                {
                    if (!inPointer && !inConstantPointer && !inConstant && !inConstReference)
                    {
                        s = s.Substring(2);
                        inConstantPointer = true;
                        continue;
                    }
                }
                if (s.StartsWith("P"))
                {
                    if (!inPointer && !inConstantPointer && !inConstant && !inConstReference)
                    {
                        s = s.Substring(1);
                        inPointer = true;
                        continue;
                    }
                }
                if (s.StartsWith("K"))
                {
                    if (!inPointer && !inConstantPointer && !inConstant && !inConstReference)
                    {
                        s = s.Substring(1);
                        inConstant = true;
                        continue;
                    }
                }

                var match = Regex.Match(s, @"^(\d+)");
                if (match.Success)
                {
                    string sL = match.Groups[1].Value;
                    int iL = int.Parse(sL);
                    string sVariable = s.Substring(sL.Length, iL);
                    
                    if (s.Length == iL + sL.Length)
                    {
                        s = s.Substring(sL.Length + iL);
                        if (inPointer)
                        {
                            result = result + sVariable + "*, ";
                            inPointer = false;
                        }
                        else if (inConstantPointer)
                        {
                            result = result + "const " + sVariable + "*, ";
                            inConstantPointer = false;
                        }
                        else if (inConstant)
                        {
                            result = result + "const " + sVariable + ", ";
                            inConstant = false;
                        }
                        else if (inConstReference)
                        {
                            result = result + "const " + sVariable + "&, ";
                            inConstReference = false;
                        }
                        else
                        {
                            result = result + sVariable + ", ";
                        }
                        
                        parsed.Add(sVariable);
                        continue;
                    }
                    else if ("E" == s.Substring(iL + sL.Length, 1))
                    {
                        s = s.Substring(sL.Length + iL + 1);
                        if (inPointer)
                        {
                            result = result + sVariable + "*, ";
                            inPointer = false;
                        }
                        else if (inConstantPointer)
                        {
                            result = result + "const " + sVariable + "*, ";
                            inConstantPointer = false;
                        }
                        else if (inConstant)
                        {
                            result = result + "const " + sVariable + ", ";
                            inConstant = false;
                        }
                        else if (inConstReference)
                        {
                            result = result + "const " + sVariable + "&, ";
                            inConstReference = false;
                        }
                        else
                        {
                            result = result + sVariable + ", ";
                        }
                        parsed.Add(sVariable);
                        continue;
                    }
                    else
                    {
                        s = s.Substring(sL.Length + iL);
                        if (inPointer)
                        {
                            result = result + sVariable + "*, ";
                            inPointer = false;
                        }
                        else if (inConstantPointer)
                        {
                            result = result + "const " + sVariable + "*, ";
                            inConstantPointer = false;
                        }
                        else if (inConstant)
                        {
                            result = result + "const " + sVariable + ", ";
                            inConstant = false;
                        }
                        else if (inConstReference)
                        {
                            result = result + "const " + sVariable + "&, ";
                            inConstReference = false;
                        }
                        else
                        {
                            result = result + sVariable + ", ";
                        }
                        parsed.Add(sVariable);
                        continue;
                    }
                }

                //not recongnized
                break;
            }

            result += s;

            if (result.EndsWith(", "))
            {
                result = result.Substring(0, result.Length - 2);
            }
            return result;
        }

        public static string ReplaceTemplateArgs(string s)
        {
            for (int i = 1; i < 100; ++i)
            {
                s = s.Replace("Li" + i.ToString(), i.ToString());
            }

            var match = Regex.Match(s, @"[<, ](\d+)([a-zA-Z_][\w]*)[,<>]");
            while (match.Success)
            {
                string sL = match.Groups[1].Value;
                int iL = int.Parse(sL);
                s = s.Replace(match.Groups[1].Value + match.Groups[2].Value, match.Groups[2].Value);
                match = Regex.Match(s, @"[<, ](\d+)([a-zA-Z_][\w]*)[,<>]");
            }

            return s;
        }

        public static string RemoveHeadNS(string s)
        {
            for (int i = 0; i < m_lstNameSpaces.Count; ++i)
            {
                string shead = m_lstNameSpaces[i].Length.ToString() + m_lstNameSpaces[i];
                if (s.StartsWith(shead))
                {
                    return s.Substring(shead.Length);
                }
            }
            return s;
        }

        #region Basic Structure

        /// <summary>
        /// N70_INTERNAL_cd27788a_33_CActionGaugePlaquetteRotatingT_cu_203f0e8a_331606CLGLib17_deviceHiShifted2EhNS0_10SSmallInt4ERKNS0_6SIndexE
        /// </summary>
        /// <param name="sin"></param>
        /// <returns></returns>
        private static string RemoveInternal(string sin)
        {
            return Regex.Replace(sin, @"_INTERNAL_[a-f0-9]+_\d+_", "");
        }

        /// <summary>
        /// N6CLGLib3SubE -> CLGLib::Sub
        /// namespace::functionname
        /// </summary>
        /// <param name="encoded"></param>
        /// <returns></returns>
        private static string ParseNameSpace(string sin)
        {
            return Regex.Replace(sin, @"N\d+([A-Za-z]+)[\d]+([A-Za-z_][\w]+)E", "$1::$2(");
        }

        public static string DemangleBasicStructure(string mangled)
        {
            mangled = FindNameSpace(mangled);
            return mangled;
        }

        #endregion

        public static string Demangle(string mangled)
        {
            if (!mangled.StartsWith("_Z"))
            {
                //not Itanium
                return mangled;
            }

            string demangled = FindNameSpace(mangled);
            demangled = RemoveHeadNS(demangled);
            demangled = FindAllFunctionOrStructureNames(demangled);

            return demangled;
        }


        // 测试代码
        //public static void Main()
        //{
        //    string[] tests = {
        //    "_ZN6CLGLib10_kernelAddE",
        //    "_ZN6CLGLib10_kernelAddIdEE",
        //    "_ZN6CLGLib10_kernelAddIdEEvPT_PKS1_j",
        //    "_ZN6CLGLib21_kernelBosonValueRealIdEEvPKT_Pd"
        //};

        //    foreach (var test in tests)
        //    {
        //        Console.WriteLine($"{test} → {Demangle(test)}");
        //    }
        //}
    }
}
   