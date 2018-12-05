//=============================================================================
// FILENAME : CYaml.cpp
// 
// DESCRIPTION:
// This is the class read YAML paramters
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region CParameters

void CParameters::Dump(const STRING& indent) const
{
    appGeneral("%sScalar<string>:\n", indent.c_str());
    if (m_pStrings.size() == 0)
    {
        appGeneral("%s  (none)\n", indent.c_str());
    }
    else 
    {
        for (std::map<STRING, STRING>::const_iterator p = m_pStrings.begin(); p != m_pStrings.end(); ++p)
        {
            appGeneral("%s  key = %s\tvalue = %s\n", indent.c_str(), p->first.c_str(), p->second.c_str());
        }
    }

    appGeneral("%sVector<string>:\n", indent.c_str());
    if (m_pStringVector.size() == 0)
    {
        appGeneral("%s  (none)\n", indent.c_str());
    }
    else 
    {
        for (std::map<STRING, std::vector<STRING> >::const_iterator p = m_pStringVector.begin(); p != m_pStringVector.end(); ++p)
        {
            appGeneral("%s  key = %s\tvalue = [ ", indent.c_str(), p->first.c_str());
            for (size_t i = 0; i < p->second.size(); ++i) 
            {
                appGeneral("%s, ", p->second[i].c_str());
            }
            appGeneral("]\n");
        }
    }

    appGeneral("%sParameters:\n", indent.c_str());
    if (m_pParameters.size() == 0)
    {
        appGeneral("%s  (none)\n", indent.c_str());
    }
    else {
        for (std::map<STRING, CParameters>::const_iterator p = m_pParameters.begin(); p != m_pParameters.end(); ++p)
        {
            appGeneral("%s  key = %s, value:\n", indent.c_str(), p->first.c_str());
            p->second.Dump(indent + "    ");
        }
    }
}

#pragma endregion CParameters

#pragma region CYAMLParser

INT CYAMLParser::ParseStream(ISTREAM& iss, CParameters& params)
{
    INT retv = EXIT_SUCCESS;

    const size_t buf_size = 1024;
    TCHAR buf[buf_size];

    typedef std::pair<STRING, CParameters*> env_t;
    typedef std::pair<INT, env_t> level_t;
    typedef std::stack<level_t> stack_t;

    stack_t levels;

    CParameters *current_params = &params;
    INT current_indent = 0;

    BOOL expect_map = FALSE;

    while (iss.getline(buf, buf_size))
    {
        STRING key, value;

        INT indent = ParseLine(buf, key, value);

        if (indent < 0) 
        {
            appParanoiac(_T("CYAMLParser: empty line. skip.\n"));
            continue;
        }

        // level up/down
        if (indent > current_indent) 
        {
            if (!expect_map) 
            {
                appGeneral(_T("CYAMLParser: Error: unexpected nest level.\n"));
                retv = EXIT_FAILURE;
                continue;
            }

            // start new level
            current_params = new CParameters;
            current_indent = indent;

            expect_map = FALSE;
        }
        else 
        {
            if (expect_map) 
            {
                // open key in the previous line actually correspond to empty value
                level_t lv = levels.top();
                levels.pop();

                STRING     key_s = lv.second.first;
                CParameters *stored_params = lv.second.second;

                stored_params->SetStringVaule(key_s, STRING()); //NULL STRING
                current_params = stored_params;
                current_indent = lv.first;

                expect_map = FALSE;
            }

            if (indent < current_indent) 
            {
                while (indent < current_indent)
                {
                    level_t lv = levels.top();
                    levels.pop();

                    STRING     key_s = lv.second.first;
                    CParameters *stored_params = lv.second.second;
                    
                    stored_params->SetParameterVaule(key_s, *current_params);
                    appSafeDelete(current_params);

                    // restore upper level
                    current_params = stored_params;
                    current_indent = lv.first;
                }
            }
            else 
            {  // indent == current_indent
                    //
            }
        }

        // store key-value
        if (value.length() > 0) 
        {
            if (value[0] == '[') 
            {
                memset(buf, '\0', buf_size);
                appStrcpy(buf, buf_size, value.c_str());
                //value.copy(buf, buf_size);  // reuse buffer

                std::vector<STRING> v;

                INT nvalues = ParseVector(buf, v);

                if (nvalues < 0) 
                {
                    appGeneral("YAMLParser: ERROR: parse_vector failed.\n");
                    continue;
                }

                // store key - vector value.
                current_params->SetStringVectorVaule(key, v);
            }
            else 
            {
                // store key - scalar value.
                current_params->SetStringVaule(key, value);
            }
        }
        else 
        {
            // key for a map in subsequent lines
            expect_map = TRUE;

            // push current environment to stack
            levels.push(level_t(indent, env_t(key, current_params)));
        }
    }

    while (current_indent > 0)
    {
        level_t lv = levels.top();
        levels.pop();

        STRING     key = lv.second.first;
        CParameters *stored_params = lv.second.second;

        stored_params->SetParameterVaule(key, *current_params);
        appSafeDelete(current_params);

        // restore upper level
        current_params = stored_params;
        current_indent = lv.first;
    }

    return retv;
}

INT CYAMLParser::ParseLine(TCHAR *buf, STRING& key, STRING& value)
{
    // N.B. buf modified on exit.

    const TCHAR delim = _T(':');

    // remove comments
    if (TCHAR *q = appStrchr(buf, _T('#'))) { *q = _T('\0'); }

    // remove trailing spaces
    TCHAR *s = buf + appStrlen(buf) - 1;
    while (*s == _T(' '))
    {
        *s-- = _T('\0');
    }

    // find indent
    INT indent = 0;

    TCHAR *p = buf;
    while (*p == _T(' '))
    {
        ++p;
        ++indent;
    }

    // find key-value separator
    TCHAR *q = appStrchr(buf, delim);

    if (!q) 
    {
        key = STRING();
        value = STRING();

        return -1;
    }

    // find key
    TCHAR *r = q;

    *r = _T('\0');
    --r;
    while (r >= p && *r == _T(' '))
    {
        *r-- = _T('\0');
    }

    key = STRING(p);

    // find value
    ++q;
    while (*q == _T(' '))
    {
        ++q;
    }

    value = STRING(q);

    // return indent
    return indent;
}

INT CYAMLParser::ParseVector(TCHAR *buf, std::vector<STRING>& vec)
{
    // N.B. buf modified on exit.
    const TCHAR sep = _T(',');

    INT count = 0;

    if ((buf[0] != _T('[')) || (buf[appStrlen(buf) - 1] != _T(']'))) {
        return -1;
    }

    buf[appStrlen(buf) - 1] = _T('\0');

    TCHAR *p = buf + 1;

    while (p && *p)
    {
        // skip heading spaces
        while (*p == _T(' '))
        {
            ++p;
        }

        // find separator
        TCHAR *q = appStrchr(p, sep);

        if (q) 
        {
            // eliminate separator
            *q = _T('\0');

            // eliminate trailing spaces of the item
            TCHAR *r = q - 1;
            while (*r == _T(' '))
            {
                *r-- = _T('\0');
            }

            vec.push_back(STRING(p));
            ++count;

            // go to next item
            p = q + 1;
        }
        else 
        {
            // separator is not found; maybe last item in the sequence.

            // eliminate trailing spaces;
            TCHAR *r = p + appStrlen(p) - 1;
            while (r >= p && *r == _T(' '))
            {
                *r-- = _T('\0');
            }

            if (appStrlen(p) > 0) {
                vec.push_back(STRING(p));
                ++count;
            }
            else 
            {
                // discard
            }

            p = q;
        }
    }

    return count;
}

void CYAMLParser::ParseFile(const STRING& params_file, CParameters& params)
{
    INT  filesize = 0;
    TCHAR *buf = 0;

    IFSTREAM fin(params_file.c_str());
    if (!fin) 
    {
        appCrucial("Error at YAML: unable to read parameter file: %s.\n", params_file.c_str());
    }

    fin.seekg(0, std::ios::end);
    filesize = (INT)fin.tellg();
    fin.seekg(0, std::ios::beg);

    INT padding = 8 - (filesize % 8);

    appParanoiac("YAML::read_params: filesize = %d, padding = %d\n", filesize, padding);

    filesize += padding;


    buf = new TCHAR[filesize];
    memset(buf, 0, filesize);

    fin.read(buf, filesize - padding);


    ISTRINGSTREAM iss(buf);
    Parse(iss, params);

    appSafeDeleteArray(buf);
}

#pragma endregion CYAMLParser

__END_NAMESPACE


//====================================================================
//============================================================END=====
