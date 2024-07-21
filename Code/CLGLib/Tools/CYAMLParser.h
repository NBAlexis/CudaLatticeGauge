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

#define _FetchFunction(typen) UBOOL FetchValue##typen(const CCString& key, typen& value) const \
{ \
    CCString v; \
    if (!FetchStringValue(key, v)) \
    { \
        return FALSE; \
    } \
    value = appStrTo##typen(v); \
    return TRUE; \
}


#define _FetchFunctionArray(typen) UBOOL FetchValueArray##typen(const CCString& key, TArray<typen>& value) const \
{ \
    value.RemoveAll();   \
    TArray<CCString> vs; \
    if (!FetchStringVectorValue(key, vs)) \
    { \
        return FALSE; \
    } \
    for (INT i = 0; i < vs.Num(); ++i) \
    { \
        value.AddItem(appStrTo##typen(vs[i])); \
    } \
    return TRUE; \
}



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

    CParameters() : m_iLine(-1) { }
    CParameters(const CCString& sName, const CCString& sFileName, INT iLine)
        : m_sName(sName) 
        , m_sFileName(sFileName)
        , m_iLine(iLine)
    {
    }

    CParameters(const CParameters& other)
        : m_sName(other.m_sName)
        , m_sFileName(other.m_sFileName)
        , m_iLine(other.m_iLine)
    { 
        m_pStrings = other.m_pStrings;
        m_pStringVector = other.m_pStringVector;
        m_pParameters = other.m_pParameters;
    }

    ~CParameters() {}

    void SetStringVaule(const CCString& key, const CCString& value)
    {
        m_pStrings[key] = value;
    }

    void SetStringVectorVaule(const CCString& key, const TArray<CCString>& value)
    {
        m_pStringVector[key] = value;
    }

    void SetParameterVaule(const CCString& key, const CParameters& value)
    {
        m_pParameters[key] = value;
    }

    CParameters& GetParameter(const CCString& key)
    {
        if (m_pParameters.Exist(key))
        {
            return m_pParameters[key];
        }

        appCrucial(_T("key '%s' not found.\n"), key.c_str());
        return *this;
    }

    _FetchFunction(INT)

    _FetchFunction(Real)

    _FetchFunction(FLOAT)

    _FetchFunction(DOUBLE)

    UBOOL FetchStringValue(const CCString& key, CCString& value) const
    {
        if (m_pStrings.Lookup(key, value))
        {
            return TRUE;
        }

        return FALSE;
    }

    _FetchFunctionArray(INT)

    _FetchFunctionArray(BYTE)

    _FetchFunctionArray(UINT)

    _FetchFunctionArray(Real)

    _FetchFunctionArray(FLOAT)

    _FetchFunctionArray(DOUBLE)

    UBOOL FetchStringVectorValue(const CCString& key, TArray<CCString>& value) const
    {
        if (m_pStringVector.Lookup(key, value))
        {
            return TRUE;
        }
        return FALSE;
    }

    UBOOL FetchParameterValue(const CCString& key, CParameters& value) const
    {
        if (m_pParameters.Lookup(key, value))
        {
            return TRUE;
        }

        return FALSE;
    }

    UBOOL Exist(const CCString& key) const
    {
        if (m_pStrings.Exist(key))
        {
            return TRUE;
        }
        if (m_pStringVector.Exist(key))
        {
            return TRUE;
        }
        if (m_pParameters.Exist(key))
        {
            return TRUE;
        }
        return FALSE;
    }

    void Dump(const CCString& indent = "") const;

    inline void Copy(const CParameters& other)
    {
        m_sName = other.m_sName;
        m_sFileName = other.m_sFileName;
        m_iLine = other.m_iLine;
        m_pStrings = other.m_pStrings;
        m_pStringVector = other.m_pStringVector;
        m_pParameters = other.m_pParameters;
    }

    inline CParameters& operator=(const CParameters& other)
    {
        Copy(other);
        return *this;
    }

    CCString GetName() const { return m_sName; }
    CCString GetLocation() const 
    { 
        static TCHAR path[CCString::_CLG_MAX_PATH];
        appGetPath(path, CCString::_CLG_MAX_PATH);
        const CCString strpath(path);
#if _CLG_WIN
        CCString ret = strpath + _T("/") + m_sFileName + _T("(") + appToString(m_iLine) + _T(")");
#else
        CCString ret = _T(" file:///") + strpath + _T("/") + m_sFileName + _T("  line:") + appToString(m_iLine);
#endif
        ret = ret.Replace(_T("\\"), _T("/"));
        return ret;
    }

    void SetFile(const CCString& path, INT iLine) { m_sFileName = path; m_iLine = iLine; }
    UBOOL ValidFile() const { return m_iLine >= 0; }

    void RemoveAll()
    {
        m_sName = _T("");
        m_sFileName = _T("");
        m_iLine = -1;

        m_pStrings.RemoveAll();
        m_pStringVector.RemoveAll();
        m_pParameters.RemoveAll();
    }

private:

    CCString m_sName;
    CCString m_sFileName;
    INT m_iLine;

    // scalar
    THashMap<CCString, CCString>               m_pStrings;
    // array
    THashMap<CCString, TArray<CCString>>       m_pStringVector;
    // map
    THashMap<CCString, CParameters>            m_pParameters;
};

class CLGAPI CYAMLParser
{
public:
    CYAMLParser() {}

    //! read parameters from file.
    static void ParseFile(const CCString& params_file, CParameters& params);

    //These are functions for parse
    static void Parse(const CCString& sName, ISTREAM& iss, CParameters& params)
    {
        const INT result = ParseStream(sName, sName, 0, iss, params);

        if (result != EXIT_SUCCESS) 
        {
            appCrucial(_T("YAML: parse failed.\n"));
            exit(EXIT_FAILURE);
        }
    }

    static INT ParseStream(const CCString &sName, const CCString& sFileName, INT iStartLine, ISTREAM& iss, CParameters& params);
    static INT ParseLine(TCHAR *buf, CCString& key, CCString& value);
    static INT ParseVector(TCHAR *buf, TArray<CCString>& vec);
};

__END_NAMESPACE

#endif //#ifndef _YAMLPARSER_H_

//=============================================================================
// END OF FILE
//=============================================================================

