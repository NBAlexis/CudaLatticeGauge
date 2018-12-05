//=============================================================================
// FILENAME : WinFunction.h
// 
// DESCRIPTION:
// This is the class read YAML paramters
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#ifndef _YAMLPARSER_H_
#define _YAMLPARSER_H_

__BEGIN_NAMESPACE

//! Parameter manager with YAML parser.

/*!
This is a simple parser to read parameters from a file
prepared with YAML format.
Only simple cases were checked.
[17 Jul 2012 H.Matsufuru]
*/

//TODO(YAML only support strings now...)

class CLGAPI CParameters
{
public:

    CParameters() {;}
    ~CParameters() {;}

    void SetStringVaule(const STRING& key, const STRING& value)
    {
        m_pStrings.insert_or_assign(key, value);
    }

    void SetStringVectorVaule(const STRING& key, std::vector<STRING>& value)
    {
        m_pStringVector.insert_or_assign(key, value);
    }

    void SetParameterVaule(const STRING& key, const CParameters& value)
    {
        m_pParameters.insert_or_assign(key, value);
    }

    CParameters& GetParameter(const STRING& key)
    {
        std::map<STRING, CParameters>::iterator p = m_pParameters.find(key);
        if (p != m_pParameters.end())
        {
            return p->second;
        }

        appCrucial("key '%s' not found.\n", key.c_str());
        return *this;
    }

    BOOL FetchStringValue(const STRING& key, STRING& value) const
    {
        std::map<STRING, STRING>::const_iterator p = m_pStrings.find(key);
        if (p != m_pStrings.end())
        {
            value = p->second;
            return TRUE;
        }

        return FALSE;
    }

    BOOL FetchStringVectorValue(const STRING& key, std::vector<STRING>& value) const
    {
        std::map<STRING, std::vector<STRING>>::const_iterator p = m_pStringVector.find(key);
        if (p != m_pStringVector.end())
        {
            value = p->second;
            return TRUE;
        }

        return FALSE;
    }

    BOOL FetchParameterValue(const STRING& key, CParameters& value) const
    {
        std::map<STRING, CParameters>::const_iterator p = m_pParameters.find(key);
        if (p != m_pParameters.end())
        {
            value = p->second;
            return TRUE;
        }

        return FALSE;
    }

    BOOL Exist(const STRING& key) const
    {
        if (m_pStrings.find(key) != m_pStrings.end())
        {
            return TRUE;
        }
        if (m_pStringVector.find(key) != m_pStringVector.end())
        {
            return TRUE;
        }
        if (m_pParameters.find(key) != m_pParameters.end())
        {
            return TRUE;
        }
        return FALSE;
    }

    void Dump(const STRING& indent = "") const;

private:

    // scalar
    std::map<STRING, STRING>                m_pStrings;
    // array
    std::map<STRING, std::vector<STRING>>   m_pStringVector;
    // map
    std::map<STRING, CParameters>           m_pParameters;
};

class CLGAPI CYAMLParser
{
public:
    CYAMLParser() {}

    //! read parameters from file.
    static void ParseFile(const STRING& params_file, CParameters& params);

    //These are functions for parse
    static void Parse(ISTREAM& iss, CParameters& params)
    {
        INT result = ParseStream(iss, params);

        if (result != EXIT_SUCCESS) 
        {
            appCrucial("YAML: parse failed.\n");
            exit(EXIT_FAILURE);
        }
    }

    static INT ParseStream(ISTREAM& iss, CParameters& params);
    static INT ParseLine(TCHAR *buf, STRING& key, STRING& value);
    static INT ParseVector(TCHAR *buf, std::vector<STRING>& vec);
};

__END_NAMESPACE

#endif //#ifndef _YAMLPARSER_H_

//=============================================================================
// END OF FILE
//=============================================================================

