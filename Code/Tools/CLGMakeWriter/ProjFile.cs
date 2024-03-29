using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;
using System.Xml;


namespace CLGMakeWriter
{
    class CProjFile
    {
        public CProjFile(string sName, string sFileName, string sProjPath)
        {
            m_sName = sName;
            m_sProjectDir = sProjPath + "/";
            m_sContent = File.ReadAllText(sFileName);
            XmlDocument doc = new XmlDocument();
            doc.LoadXml(m_sContent);

            XmlNodeList allClInclude = doc.GetElementsByTagName("ClInclude");
            for (int i = 0; i < allClInclude.Count; ++i)
            {
                m_lstAllHeaderFiles.Add(allClInclude[i].Attributes.GetNamedItem("Include").InnerText);
            }

            XmlNodeList allClCompile = doc.GetElementsByTagName("ClCompile");
            for (int i = 0; i < allClCompile.Count; ++i)
            {
                if (null == allClCompile[i].Attributes.GetNamedItem("Include"))
                {
                    //those are compile settings
                }
                else
                {
                    //m_lstAllCppFiles.Add(Path.GetFullPath(Path.Combine(sProjPath, allClCompile[i].Attributes.GetNamedItem("Include").InnerText)));
                    m_lstAllCppFiles.Add(allClCompile[i].Attributes.GetNamedItem("Include").InnerText);
                }
            }

            XmlNodeList allCuCompile = doc.GetElementsByTagName("CudaCompile");
            for (int i = 0; i < allCuCompile.Count; ++i)
            {
                if (null == allCuCompile[i].Attributes.GetNamedItem("Include"))
                {
                    //those are nvcc settings
                    ///<seealso cref="m_sContent">
                }
                else
                {
                    m_lstAllCuFiles.Add(allCuCompile[i].Attributes.GetNamedItem("Include").InnerText);
                }
            }
        }


        public CProjFile(string sAppName)
        {
            m_sName = sAppName;

            string sProjPath = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, 
                string.Format("../../Code/Applications/{0}", sAppName)
            });
            string sFileName = Path.Combine(new[] { System.AppDomain.CurrentDomain.BaseDirectory, 
                string.Format("../../Code/Applications/{0}/{0}.vcxproj", sAppName)
            });

            m_sProjectDir = sProjPath + "/";
            m_sContent = File.ReadAllText(sFileName);
            XmlDocument doc = new XmlDocument();
            doc.LoadXml(m_sContent);

            XmlNodeList allClInclude = doc.GetElementsByTagName("ClInclude");
            for (int i = 0; i < allClInclude.Count; ++i)
            {
                m_lstAllHeaderFiles.Add(allClInclude[i].Attributes.GetNamedItem("Include").InnerText);
            }

            XmlNodeList allClCompile = doc.GetElementsByTagName("ClCompile");
            for (int i = 0; i < allClCompile.Count; ++i)
            {
                if (null == allClCompile[i].Attributes.GetNamedItem("Include"))
                {
                    //those are compile settings
                }
                else
                {
                    //m_lstAllCppFiles.Add(Path.GetFullPath(Path.Combine(sProjPath, allClCompile[i].Attributes.GetNamedItem("Include").InnerText)));
                    m_lstAllCppFiles.Add(allClCompile[i].Attributes.GetNamedItem("Include").InnerText);
                }
            }

            XmlNodeList allCuCompile = doc.GetElementsByTagName("CudaCompile");
            for (int i = 0; i < allCuCompile.Count; ++i)
            {
                if (null == allCuCompile[i].Attributes.GetNamedItem("Include"))
                {
                    //those are nvcc settings
                    ///<seealso cref="m_sContent">
                }
                else
                {
                    m_lstAllCuFiles.Add(allCuCompile[i].Attributes.GetNamedItem("Include").InnerText);
                }
            }
        }

        public int[] FindVersionNumber()
        {
            int[] ret = {0,0};

            foreach (string sHeader in m_lstAllHeaderFiles)
            {
                if (sHeader.Contains("CLGDefine.h"))
                {
                    string sFullPath = Path.GetFullPath(Path.Combine(m_sProjectDir, sHeader));
                    string sFileContent = File.ReadAllText(sFullPath);
                    Match allMatch = Regex.Match(sFileContent, @"__GVERSION[\s]+\(([\d])+\)");
                    if (allMatch.Success && allMatch.Groups.Count > 1)
                    {
                        int iMajor = 1;
                        int.TryParse(allMatch.Groups[1].ToString(), out iMajor);
                        ret[0] = iMajor;
                    }

                    allMatch = Regex.Match(sFileContent, @"__GVERSION_S[\s]+\(([\d])+\)");
                    if (allMatch.Success && allMatch.Groups.Count > 1)
                    {
                        int iMinor = 0;
                        int.TryParse(allMatch.Groups[1].ToString(), out iMinor);
                        ret[1] = iMinor;
                    }
                }
            }

            return ret;
        }

        public readonly List<string> m_lstAllHeaderFiles = new List<string>();
        public readonly List<string> m_lstAllCuFiles = new List<string>();
        public readonly List<string> m_lstAllCppFiles = new List<string>();

        public string m_sName;
        public string m_sContent;
        public string m_sProjectDir;
    }
}

