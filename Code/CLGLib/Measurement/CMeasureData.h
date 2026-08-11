//=============================================================================
// FILENAME : CMeasureData.h
// 
// DESCRIPTION:
// Type-erased container for measurement result arrays, plus CSV / NumPy export
// helpers. All data arrays are owned by CMeasurementManager and stored in a
// CMemStack.
//
// REVISION:
//  [06/14/2026]
//=============================================================================
#pragma once

#ifndef _CMEASUREDATA_H_
#define _CMEASUREDATA_H_

#include "Core/CLGDefine.h"
#include "Core/CLGFloat.h"
#include "Tools/Data/CCString.h"
#include "Tools/Data/TArray.h"
#include "Tools/Data/MemStack.h"
#include "Tools/Data/THashMap.h"
#include "Tools/CYAMLParser.h"
#include "Tools/Tracer.h"

#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <stdint.h>

__BEGIN_NAMESPACE

#if !_CLG_WIN

inline void strerror_s(TCHAR* buffer, size_t bufferSize, INT error)
{
    strcpy(buffer, strerror(error));
}

inline void _gcvt_s(TCHAR* buff, UINT uiBuffLength, Real fVaule, UINT uiDigit)
{
    static TCHAR tmpBuff[10];
    appSprintf(tmpBuff, 10, _T("%s.%df"), _T("%"), uiDigit);
    appSprintf(buff, uiBuffLength, tmpBuff, fVaule);
}

#endif

enum { kFileDigital = 20, };

template <class T>
void WriteRealArray(const CCString& sFileName, const TArray<T>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kFileDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.Num(); ++i)
    {
        _gcvt_s(str, 50, static_cast<DOUBLE>(lst[i]), iDigital);
        CCString sReal = CCString(str);
        sReal = sReal.Replace(_T("e"), _T("*^"));
        file << _T(" ");
        file << sReal;
        if (i != lst.GetCount() - 1)
        {
            file << _T(",");
        }
    }
    file.flush();
    file.close();
}

template <class T>
void WriteRealArray2(const CCString& sFileName, const TArray<TArray<T>>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kFileDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        for (INT j = 0; j < lst[i].GetCount(); ++j)
        {
            _gcvt_s(str, 50, lst[i][j], iDigital);
            CCString sReal = CCString(str);
            sReal = sReal.Replace(_T("e"), _T("*^"));
            file << _T(" ");
            file << sReal;
            if (j != lst[i].GetCount() - 1)
            {
                file << _T(",");
            }
        }
        file << _T("\n");
    }
    file.flush();
    file.close();
}

inline void WriteComplexArray(const CCString& sFileName, const TArray<cuComplex>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kFileDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.Num(); ++i)
    {
        _gcvt_s(str, 50, lst[i].x, iDigital);
        CCString sReal = CCString(str);
        sReal = sReal.Replace(_T("e"), _T("*^"));
        file << _T(" ");
        file << sReal;

        _gcvt_s(str, 50, lst[i].y, iDigital);
        CCString sImg = CCString(str);
        sImg = sImg.Replace(_T("e"), _T("*^"));
        CCString sMid = _T(" + ");
        if (sImg.Left(1) == _T("-"))
        {
            sImg = sImg.Right(sImg.GetLength() - 1);
            sMid = _T(" - ");
        }
        file << sMid;
        file << sImg;
        if (i == lst.GetCount() - 1)
        {
            file << _T(" I");
        }
        else
        {
            file << _T(" I,");
        }
    }
    file.flush();
    file.close();
}

inline void WriteComplexArray(const CCString& sFileName, const TArray<cuDoubleComplex>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kFileDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.Num(); ++i)
    {
        _gcvt_s(str, 50, lst[i].x, iDigital);
        CCString sReal = CCString(str);
        sReal = sReal.Replace(_T("e"), _T("*^"));
        file << _T(" ");
        file << sReal;

        _gcvt_s(str, 50, lst[i].y, iDigital);
        CCString sImg = CCString(str);
        sImg = sImg.Replace(_T("e"), _T("*^"));
        CCString sMid = _T(" + ");
        if (sImg.Left(1) == _T("-"))
        {
            sImg = sImg.Right(sImg.GetLength() - 1);
            sMid = _T(" - ");
        }
        file << sMid;
        file << sImg;
        if (i == lst.GetCount() - 1)
        {
            file << _T(" I");
        }
        else
        {
            file << _T(" I,");
        }
    }
    file.flush();
    file.close();
}

inline void WriteComplexArray2(const CCString& sFileName, const TArray<TArray<cuComplex>>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kFileDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        for (INT j = 0; j < lst[i].GetCount(); ++j)
        {
            _gcvt_s(str, 50, lst[i][j].x, iDigital);
            CCString sReal = CCString(str);
            sReal = sReal.Replace(_T("e"), _T("*^"));
            _gcvt_s(str, 50, lst[i][j].y, iDigital);
            CCString sImg = CCString(str);
            sImg = sImg.Replace(_T("e"), _T("*^"));
            CCString sMid = _T(" + ");
            if (sImg.Left(1) == _T("-"))
            {
                sImg = sImg.Right(sImg.GetLength() - 1);
                sMid = _T(" - ");
            }
            file << _T(" ");
            file << sReal;
            file << sMid;
            file << sImg;
            if (j == lst[i].GetCount() - 1)
            {
                file << _T(" I");
            }
            else
            {
                file << _T(" I,");
            }
        }
        file << _T("\n");
    }
    file.flush();
    file.close();
}

inline void WriteComplexArray2(const CCString& sFileName, const TArray<TArray<cuDoubleComplex>>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kFileDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        for (INT j = 0; j < lst[i].GetCount(); ++j)
        {
            _gcvt_s(str, 50, lst[i][j].x, iDigital);
            CCString sReal = CCString(str);
            sReal = sReal.Replace(_T("e"), _T("*^"));
            _gcvt_s(str, 50, lst[i][j].y, iDigital);
            CCString sImg = CCString(str);
            sImg = sImg.Replace(_T("e"), _T("*^"));
            CCString sMid = _T(" + ");
            if (sImg.Left(1) == _T("-"))
            {
                sImg = sImg.Right(sImg.GetLength() - 1);
                sMid = _T(" - ");
            }
            file << _T(" ");
            file << sReal;
            file << sMid;
            file << sImg;
            if (j == lst[i].GetCount() - 1)
            {
                file << _T(" I");
            }
            else
            {
                file << _T(" I,");
            }
        }
        file << _T("\n");
    }
    file.flush();
    file.close();
}

inline void WriteComplexArray2Simple(const CCString& sFileName, const TArray<TArray<CLGComplex>>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kFileDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        for (INT j = 0; j < lst[i].GetCount(); ++j)
        {
            _gcvt_s(str, 50, lst[i][j].x, iDigital);
            CCString sReal = CCString(str);
            sReal = sReal.Replace(_T("e"), _T("*^"));
            file << _T(" ");
            file << sReal;

            if (abs(lst[i][j].y) > F(0.0001) * lst[i][j].x)
            {
                _gcvt_s(str, 50, lst[i][j].y, iDigital);
                CCString sImg = CCString(str);
                sImg = sImg.Replace(_T("e"), _T("*^"));
                CCString sMid = _T(" + ");
                if (sImg.Left(1) == _T("-"))
                {
                    sImg = sImg.Right(sImg.GetLength() - 1);
                    sMid = _T(" - ");
                }

                file << sMid;
                file << sImg;
                if (j == lst[i].GetCount() - 1)
                {
                    file << _T(" I");
                }
                else
                {
                    file << _T(" I,");
                }
            }
            else
            {
                if (j != lst[i].GetCount() - 1)
                {
                    file << _T(",");
                }
            }
        }
        file << _T("\n");
    }
    file.flush();
    file.close();
}

template<typename T>
struct TNpyElementType { using Type = T; };

template<typename T, typename A>
struct TNpyElementType<TArray<T, A>> { using Type = typename TNpyElementType<T>::Type; };

template<typename T>
inline CCString _NpyDtype()
{
    appCrucial(_T("SaveAsNumpyFile: unsupported data type for NumPy output\n"));
    return _T("<f8");
}

template<> inline CCString _NpyDtype<FLOAT>() { return _T("<f4"); }
template<> inline CCString _NpyDtype<DOUBLE>() { return _T("<f8"); }
template<> inline CCString _NpyDtype<INT>() { return sizeof(INT) == 4 ? _T("<i4") : _T("<i8"); }
template<> inline CCString _NpyDtype<UINT>() { return sizeof(UINT) == 4 ? _T("<u4") : _T("<u8"); }
template<> inline CCString _NpyDtype<cuComplex>() { return _T("<c8"); }
template<> inline CCString _NpyDtype<cuDoubleComplex>() { return _T("<c16"); }

template<typename T>
inline void _SaveAsNumpyFileImpl(
    const CCString& sFileName,
    const T* pData,
    INT nTotalElements,
    const TArray<INT>& shape)
{
    std::string sHeader = "{'descr': '";
    sHeader += _NpyDtype<T>().c_str();
    sHeader += "', 'fortran_order': False, 'shape': (";

    for (INT i = 0; i < shape.Num(); ++i)
    {
        sHeader += appToString(shape[i]).c_str();
        if (i < shape.Num() - 1)
        {
            sHeader += ", ";
        }
    }
    if (shape.Num() == 1)
    {
        sHeader += ",";
    }
    sHeader += "), ";
    sHeader += "}";

    // NPY format: the header is padded with spaces so that the total of
    // magic+version+len+header is a multiple of 64, and the *last* byte must
    // be '\n'. Pad first, then terminate -- a trailing newline followed by
    // spaces makes numpy fail with "Cannot parse header".
    size_t nPrefix = 6 + 2 + 2;
    size_t nTotal = nPrefix + sHeader.size() + 1;
    size_t nPadding = (64 - (nTotal % 64)) % 64;
    sHeader.append(nPadding, ' ');
    sHeader += "\n";

    std::ofstream file(sFileName.c_str(), std::ios::binary);

    if (file.fail())
    {
        appCrucial(_T("SaveAsNumpyFile: failed to open %s\n"), sFileName.c_str());
        return;
    }

    file.write("\x93NUMPY", 6);
    unsigned char version[2] = {0x01, 0x00};
    file.write(reinterpret_cast<const char*>(version), 2);
    unsigned short usHeaderLen = static_cast<unsigned short>(sHeader.size());
    file.write(reinterpret_cast<const char*>(&usHeaderLen), 2);
    file.write(sHeader.c_str(), static_cast<std::streamsize>(sHeader.size()));
    file.write(reinterpret_cast<const char*>(pData),
        static_cast<std::streamsize>(nTotalElements * sizeof(T)));
    file.close();
}

inline void _NpyAppend(std::vector<BYTE>& v, const void* pData, size_t uiSize)
{
    if (0 == uiSize)
    {
        return;
    }
    const BYTE* pBytes = reinterpret_cast<const BYTE*>(pData);
    v.insert(v.end(), pBytes, pBytes + uiSize);
}

inline void _NpyAppendU16(std::vector<BYTE>& v, UINT uiValue)
{
    v.push_back(static_cast<BYTE>(uiValue & 0xff));
    v.push_back(static_cast<BYTE>((uiValue >> 8) & 0xff));
}

inline void _NpyAppendU32(std::vector<BYTE>& v, UINT uiValue)
{
    v.push_back(static_cast<BYTE>(uiValue & 0xff));
    v.push_back(static_cast<BYTE>((uiValue >> 8) & 0xff));
    v.push_back(static_cast<BYTE>((uiValue >> 16) & 0xff));
    v.push_back(static_cast<BYTE>((uiValue >> 24) & 0xff));
}

inline std::string _NpyBuildHeader(const CCString& sDtype, const TArray<INT>& shape)
{
    std::string sHeader = "{'descr': '";
    sHeader += sDtype.c_str();
    sHeader += "', 'fortran_order': False, 'shape': (";

    for (INT i = 0; i < shape.Num(); ++i)
    {
        sHeader += appToString(shape[i]).c_str();
        if (i < shape.Num() - 1)
        {
            sHeader += ", ";
        }
    }
    if (shape.Num() == 1)
    {
        sHeader += ",";
    }
    sHeader += "), }";
    sHeader += "\n";

    const size_t nPrefix = 6 + 2 + 2;
    const size_t nTotal = nPrefix + sHeader.size();
    const size_t nPadding = (64 - (nTotal % 64)) % 64;
    sHeader.insert(sHeader.size() - 1, nPadding, ' ');
    return sHeader;
}

inline std::vector<BYTE> _NpyBuildBytes(
    const void* pData,
    size_t uiDataBytes,
    const TArray<INT>& shape,
    const CCString& sDtype)
{
    std::vector<BYTE> ret;
    const std::string sHeader = _NpyBuildHeader(sDtype, shape);
    ret.reserve(6 + 2 + 2 + sHeader.size() + uiDataBytes);
    _NpyAppend(ret, "\x93NUMPY", 6);
    ret.push_back(0x01);
    ret.push_back(0x00);
    _NpyAppendU16(ret, static_cast<UINT>(sHeader.size()));
    _NpyAppend(ret, sHeader.c_str(), sHeader.size());
    _NpyAppend(ret, pData, uiDataBytes);
    return ret;
}

template<typename T>
inline std::vector<BYTE> _NpyBuildBytes(
    const T* pData,
    INT nTotalElements,
    const TArray<INT>& shape)
{
    return _NpyBuildBytes(pData, static_cast<size_t>(nTotalElements) * sizeof(T), shape, _NpyDtype<T>());
}

inline std::string _NpyPickleQuote(const CCString& sValue)
{
    std::string ret = "'";
    const std::string value = sValue.c_str();
    for (size_t i = 0; i < value.size(); ++i)
    {
        const unsigned char c = static_cast<unsigned char>(value[i]);
        switch (c)
        {
        case '\\':
            ret += "\\\\";
            break;
        case '\'':
            ret += "\\'";
            break;
        case '\n':
            ret += "\\n";
            break;
        case '\r':
            ret += "\\r";
            break;
        case '\t':
            ret += "\\t";
            break;
        default:
            if (c < 32 || c >= 127)
            {
                ret.push_back('\\');
                ret.push_back(static_cast<char>('0' + ((c >> 6) & 7)));
                ret.push_back(static_cast<char>('0' + ((c >> 3) & 7)));
                ret.push_back(static_cast<char>('0' + (c & 7)));
            }
            else
            {
                ret.push_back(static_cast<char>(c));
            }
            break;
        }
    }
    ret += "'";
    return ret;
}

inline std::vector<BYTE> _NpyBuildPickleDictPayload(const THashMap<CCString, CCString>& metadata)
{
    std::string payload = "(d";
    TArray<CCString> keys = metadata.GetAllKeys();
    for (INT i = 0; i < keys.Num(); ++i)
    {
        payload += "S";
        payload += _NpyPickleQuote(keys[i]);
        payload += "\nS";
        payload += _NpyPickleQuote(metadata.GetAt(keys[i]));
        payload += "\ns";
    }
    payload += ".";

    std::vector<BYTE> ret;
    _NpyAppend(ret, payload.c_str(), payload.size());
    return ret;
}

inline THashMap<CCString, CCString> _NpyMetadataFromParameters(const CParameters& metadata)
{
    THashMap<CCString, CCString> ret;
    TArray<CCString> keys = metadata.GetAllStringKeys();
    for (INT i = 0; i < keys.Num(); ++i)
    {
        CCString value;
        if (metadata.FetchStringValue(keys[i], value))
        {
            ret.SetAt(keys[i], value);
        }
    }
    return ret;
}

inline std::vector<BYTE> _NpyBuildMetadataBytes(const THashMap<CCString, CCString>& metadata)
{
    TArray<INT> shape;
    const std::vector<BYTE> payload = _NpyBuildPickleDictPayload(metadata);
    return _NpyBuildBytes(payload.data(), payload.size(), shape, _T("|O"));
}

inline UINT _NpyCrc32(const BYTE* pData, size_t uiSize)
{
    UINT crc = 0xffffffffU;
    for (size_t i = 0; i < uiSize; ++i)
    {
        crc ^= pData[i];
        for (UINT j = 0; j < 8; ++j)
        {
            crc = (crc >> 1) ^ (0xedb88320U & (0U - (crc & 1U)));
        }
    }
    return crc ^ 0xffffffffU;
}

struct SNumpyZipArray
{
    CCString m_sName;
    std::vector<BYTE> m_data;

    UBOOL operator == (const SNumpyZipArray& other) const
    {
        return m_sName == other.m_sName;
    }
};

inline void _NpyWriteZipEntryLocal(std::ofstream& file, const SNumpyZipArray& entry, UINT uiCrc)
{
    const std::string sName = entry.m_sName.c_str();
    std::vector<BYTE> vHeader;
    vHeader.reserve(30 + sName.size());
    _NpyAppendU32(vHeader, 0x04034b50U);
    _NpyAppendU16(vHeader, 20);
    _NpyAppendU16(vHeader, 0);
    _NpyAppendU16(vHeader, 0);
    _NpyAppendU16(vHeader, 0);
    _NpyAppendU16(vHeader, 0);
    _NpyAppendU32(vHeader, uiCrc);
    _NpyAppendU32(vHeader, static_cast<UINT>(entry.m_data.size()));
    _NpyAppendU32(vHeader, static_cast<UINT>(entry.m_data.size()));
    _NpyAppendU16(vHeader, static_cast<UINT>(sName.size()));
    _NpyAppendU16(vHeader, 0);
    file.write(reinterpret_cast<const char*>(vHeader.data()), static_cast<std::streamsize>(vHeader.size()));
    file.write(sName.c_str(), static_cast<std::streamsize>(sName.size()));
    file.write(reinterpret_cast<const char*>(entry.m_data.data()), static_cast<std::streamsize>(entry.m_data.size()));
}

inline void _NpyWriteZipCentral(
    std::vector<BYTE>& central,
    const SNumpyZipArray& entry,
    UINT uiCrc,
    UINT uiLocalOffset)
{
    const std::string sName = entry.m_sName.c_str();
    _NpyAppendU32(central, 0x02014b50U);
    _NpyAppendU16(central, 20);
    _NpyAppendU16(central, 20);
    _NpyAppendU16(central, 0);
    _NpyAppendU16(central, 0);
    _NpyAppendU16(central, 0);
    _NpyAppendU16(central, 0);
    _NpyAppendU32(central, uiCrc);
    _NpyAppendU32(central, static_cast<UINT>(entry.m_data.size()));
    _NpyAppendU32(central, static_cast<UINT>(entry.m_data.size()));
    _NpyAppendU16(central, static_cast<UINT>(sName.size()));
    _NpyAppendU16(central, 0);
    _NpyAppendU16(central, 0);
    _NpyAppendU16(central, 0);
    _NpyAppendU16(central, 0);
    _NpyAppendU32(central, 0);
    _NpyAppendU32(central, uiLocalOffset);
    _NpyAppend(central, sName.c_str(), sName.size());
}

inline void _SaveAsNumpyZipFileImpl(const CCString& sFileName, const TArray<SNumpyZipArray>& entries)
{
    std::ofstream file(sFileName.c_str(), std::ios::binary);
    if (file.fail())
    {
        appCrucial(_T("SaveAsNumpyZipFile: failed to open %s\n"), sFileName.c_str());
        return;
    }

    std::vector<BYTE> central;
    TArray<UINT> crcs;
    TArray<UINT> offsets;

    for (INT i = 0; i < entries.Num(); ++i)
    {
        const UINT uiOffset = static_cast<UINT>(file.tellp());
        const UINT uiCrc = _NpyCrc32(entries[i].m_data.data(), entries[i].m_data.size());
        offsets.AddItem(uiOffset);
        crcs.AddItem(uiCrc);
        _NpyWriteZipEntryLocal(file, entries[i], uiCrc);
    }

    const UINT uiCentralOffset = static_cast<UINT>(file.tellp());
    for (INT i = 0; i < entries.Num(); ++i)
    {
        _NpyWriteZipCentral(central, entries[i], crcs[i], offsets[i]);
    }
    file.write(reinterpret_cast<const char*>(central.data()), static_cast<std::streamsize>(central.size()));

    std::vector<BYTE> end;
    _NpyAppendU32(end, 0x06054b50U);
    _NpyAppendU16(end, 0);
    _NpyAppendU16(end, 0);
    _NpyAppendU16(end, static_cast<UINT>(entries.Num()));
    _NpyAppendU16(end, static_cast<UINT>(entries.Num()));
    _NpyAppendU32(end, static_cast<UINT>(central.size()));
    _NpyAppendU32(end, uiCentralOffset);
    _NpyAppendU16(end, 0);
    file.write(reinterpret_cast<const char*>(end.data()), static_cast<std::streamsize>(end.size()));
    file.close();
}

inline CCString _NpyZipMemberName(const CCString& sName)
{
    std::string name = sName.c_str();
    if (name.size() < 4 || name.substr(name.size() - 4) != ".npy")
    {
        name += ".npy";
    }
    return CCString(name.c_str());
}

inline INT _NpyTotalSize(const TArray<INT>& shape)
{
    INT nTotal = 1;
    for (INT i = 0; i < shape.Num(); ++i)
    {
        nTotal *= shape[i];
    }
    return nTotal;
}

template<typename T>
void _NpyGetShape(const TArray<T>& arr, TArray<INT>& shape)
{
    shape.AddItem(arr.Num());
}

template<typename T, typename A>
void _NpyGetShape(const TArray<TArray<T, A>>& arr, TArray<INT>& shape)
{
    shape.AddItem(arr.Num());
    if (arr.Num() > 0)
    {
        _NpyGetShape(arr[0], shape);
    }
}

// Flatten nested TArray to a linear buffer of leaf elements.
// E = leaf element type (e.g. CLGComplex for TArray<TArray<CLGComplex>>)
// T = current level type
template<typename E, typename T>
void _NpyFlattenTo(const TArray<T>& arr, E* pOut, INT& idx)
{
    for (INT i = 0; i < arr.Num(); ++i)
    {
        pOut[idx++] = arr[i];  // T == E at the leaf level
    }
}

template<typename E, typename T, typename A>
void _NpyFlattenTo(const TArray<TArray<T, A>>& arr, E* pOut, INT& idx)
{
    for (INT i = 0; i < arr.Num(); ++i)
    {
        _NpyFlattenTo<E>(arr[i], pOut, idx);
    }
}

// Deprecated: old _NpyFlatten kept for backwards compat with rank-0/1 only.
// DO NOT USE for nested arrays (rank >= 2) — the pOut type never reaches the leaf.
template<typename T>
void _NpyFlatten(const TArray<T>& arr, T* pOut, INT& idx)
{
    for (INT i = 0; i < arr.Num(); ++i)
    {
        pOut[idx++] = arr[i];
    }
}

template<typename T, typename A>
void _NpyFlatten(const TArray<TArray<T, A>>& arr, T* pOut, INT& idx)
{
    for (INT i = 0; i < arr.Num(); ++i)
    {
        _NpyFlatten(arr[i], pOut, idx);
    }
}


template<typename T>
inline void SaveAsNumpyFile(const CCString& sFileName, const T* pData, const TArray<INT>& shape)
{
    _SaveAsNumpyFileImpl(sFileName, pData, _NpyTotalSize(shape), shape);
}

template<typename T>
inline SNumpyZipArray MakeNumpyZipArray(const CCString& sName, const T* pData, const TArray<INT>& shape)
{
    SNumpyZipArray ret;
    ret.m_sName = _NpyZipMemberName(sName);
    ret.m_data = _NpyBuildBytes(pData, _NpyTotalSize(shape), shape);
    return ret;
}

inline void SaveAsNumpyZipFile(
    const CCString& sFileName,
    const THashMap<CCString, CCString>& metadata,
    const TArray<SNumpyZipArray>& arrays)
{
    TArray<SNumpyZipArray> entries;
    SNumpyZipArray metaEntry;
    metaEntry.m_sName = _T("metadata.npy");
    metaEntry.m_data = _NpyBuildMetadataBytes(metadata);
    entries.AddItem(metaEntry);
    for (INT i = 0; i < arrays.Num(); ++i)
    {
        SNumpyZipArray entry = arrays[i];
        entry.m_sName = _NpyZipMemberName(entry.m_sName);
        entries.AddItem(entry);
    }
    _SaveAsNumpyZipFileImpl(sFileName, entries);
}

inline void SaveAsNumpyZipFile(
    const CCString& sFileName,
    const CParameters& metadata,
    const TArray<SNumpyZipArray>& arrays)
{
    SaveAsNumpyZipFile(sFileName, _NpyMetadataFromParameters(metadata), arrays);
}

template<typename T>
inline void SaveAsNumpyZipFile(
    const CCString& sFileName,
    const THashMap<CCString, CCString>& metadata,
    const T* pData,
    const TArray<INT>& shape,
    const CCString& sArrayName = _T("array"))
{
    TArray<SNumpyZipArray> arrays;
    arrays.AddItem(MakeNumpyZipArray(sArrayName, pData, shape));
    SaveAsNumpyZipFile(sFileName, metadata, arrays);
}

template<typename T>
inline void SaveAsNumpyZipFile(
    const CCString& sFileName,
    const CParameters& metadata,
    const T* pData,
    const TArray<INT>& shape,
    const CCString& sArrayName = _T("array"))
{
    SaveAsNumpyZipFile(sFileName, _NpyMetadataFromParameters(metadata), pData, shape, sArrayName);
}

template<typename T>
inline void SaveAsNumpyFile(
    const CCString& sFileName,
    const T* pData,
    const TArray<INT>& shape,
    const THashMap<CCString, CCString>& extra)
{
    SaveAsNumpyZipFile(sFileName, extra, pData, shape);
}

template<typename T>
inline void SaveAsNumpyFile(
    const CCString& sFileName,
    const T* pData,
    const TArray<INT>& shape,
    const CParameters& extra)
{
    SaveAsNumpyZipFile(sFileName, extra, pData, shape);
}

template<typename T>
void SaveAsNumpyFile(const CCString& sFileName, const TArray<T>& arr)
{
    TArray<INT> shape;
    _NpyGetShape(arr, shape);
    INT nTotal = _NpyTotalSize(shape);

    CCString sShape;
    for (INT i = 0; i < shape.Num(); ++i)
    {
        sShape += appToString(shape[i]);
        if (i < shape.Num() - 1)
        {
            sShape += _T("x");
        }
    }
    appGeneral(_T("[CMeasureData] Writing NPY %s shape=(%s) total=%d\n"),
        sFileName.c_str(), sShape.c_str(), nTotal);

    using ElementType = typename TNpyElementType<T>::Type;
    TArray<ElementType> flat;
    flat.SetSize(nTotal);
    INT idx = 0;
    _NpyFlattenTo<ElementType>(arr, flat.GetData(), idx);
    _SaveAsNumpyFileImpl<ElementType>(sFileName, flat.GetData(), nTotal, shape);
}

template<typename T>
void SaveAsNumpyFile(
    const CCString& sFileName,
    const TArray<T>& arr,
    const THashMap<CCString, CCString>& extra)
{
    TArray<INT> shape;
    _NpyGetShape(arr, shape);
    INT nTotal = _NpyTotalSize(shape);
    using ElementType = typename TNpyElementType<T>::Type;
    TArray<ElementType> flat;
    flat.SetSize(nTotal);
    INT idx = 0;
    _NpyFlattenTo<ElementType>(arr, flat.GetData(), idx);
    SaveAsNumpyZipFile<ElementType>(sFileName, extra, flat.GetData(), shape);
}

template<typename T>
void SaveAsNumpyFile(
    const CCString& sFileName,
    const TArray<T>& arr,
    const CParameters& extra)
{
    TArray<INT> shape;
    _NpyGetShape(arr, shape);
    INT nTotal = _NpyTotalSize(shape);
    using ElementType = typename TNpyElementType<T>::Type;
    TArray<ElementType> flat;
    flat.SetSize(nTotal);
    INT idx = 0;
    _NpyFlattenTo<ElementType>(arr, flat.GetData(), idx);
    SaveAsNumpyZipFile<ElementType>(sFileName, extra, flat.GetData(), shape);
}

//=============================================================================
// Type-erased measurement data container.
//=============================================================================
class CLGAPI CMeasureData
{
public:

    CMeasureData() {}

    virtual ~CMeasureData() {}

    virtual void Reset() = 0;

    virtual INT Num() const = 0;

    virtual void Export(const CCString& sPath) const = 0;
};

//=============================================================================
// Compile-time rank of nested TArray.
//=============================================================================
template<typename T>
struct TMeasureRank { enum { Value = 0 }; };

template<typename T, typename A>
struct TMeasureRank<TArray<T, A>> { enum { Value = 1 + TMeasureRank<T>::Value }; };

//=============================================================================
// Helper dispatch for CSV export according to scalar type.
//=============================================================================
inline void WriteScalarCsv(const CCString& sFileName, const TArray<FLOAT>& arr)
{
    appGeneral(_T("[CMeasureData] Writing CSV %s (%d real values)\n"), sFileName.c_str(), arr.Num());
    WriteRealArray(sFileName, arr);
}
inline void WriteScalarCsv(const CCString& sFileName, const TArray<DOUBLE>& arr)
{
    appGeneral(_T("[CMeasureData] Writing CSV %s (%d double values)\n"), sFileName.c_str(), arr.Num());
    WriteRealArray(sFileName, arr);
}
inline void WriteScalarCsv(const CCString& sFileName, const TArray<cuComplex>& arr)
{
    appGeneral(_T("[CMeasureData] Writing CSV %s (%d complex values)\n"), sFileName.c_str(), arr.Num());
    WriteComplexArray(sFileName, arr);
}
inline void WriteScalarCsv(const CCString& sFileName, const TArray<cuDoubleComplex>& arr)
{
    appGeneral(_T("[CMeasureData] Writing CSV %s (%d complex values)\n"), sFileName.c_str(), arr.Num());
    WriteComplexArray(sFileName, arr);
}
inline void WriteScalarCsv(const CCString& sFileName, const TArray<INT>& arr)
{
    appGeneral(_T("[CMeasureData] Writing CSV %s (%d int values)\n"), sFileName.c_str(), arr.Num());
    WriteRealArray(sFileName, arr);
}
inline void WriteScalarCsv(const CCString& sFileName, const TArray<UINT>& arr)
{
    appGeneral(_T("[CMeasureData] Writing CSV %s (%d uint values)\n"), sFileName.c_str(), arr.Num());
    WriteRealArray(sFileName, arr);
}

inline void WriteScalarCsv2(const CCString& sFileName, const TArray<TArray<FLOAT>>& arr)
{
    appGeneral(_T("[CMeasureData] Writing 2D CSV %s (%d rows)\n"), sFileName.c_str(), arr.Num());
    WriteRealArray2(sFileName, arr);
}
inline void WriteScalarCsv2(const CCString& sFileName, const TArray<TArray<DOUBLE>>& arr)
{
    appGeneral(_T("[CMeasureData] Writing 2D CSV %s (%d rows)\n"), sFileName.c_str(), arr.Num());
    WriteRealArray2(sFileName, arr);
}
inline void WriteScalarCsv2(const CCString& sFileName, const TArray<TArray<cuComplex>>& arr)
{
    appGeneral(_T("[CMeasureData] Writing 2D CSV %s (%d rows)\n"), sFileName.c_str(), arr.Num());
    WriteComplexArray2(sFileName, arr);
}
inline void WriteScalarCsv2(const CCString& sFileName, const TArray<TArray<cuDoubleComplex>>& arr)
{
    appGeneral(_T("[CMeasureData] Writing 2D CSV %s (%d rows)\n"), sFileName.c_str(), arr.Num());
    WriteComplexArray2(sFileName, arr);
}
inline void WriteScalarCsv2(const CCString& sFileName, const TArray<TArray<INT>>& arr)
{
    appGeneral(_T("[CMeasureData] Writing 2D CSV %s (%d rows)\n"), sFileName.c_str(), arr.Num());
    WriteRealArray2(sFileName, arr);
}
inline void WriteScalarCsv2(const CCString& sFileName, const TArray<TArray<UINT>>& arr)
{
    appGeneral(_T("[CMeasureData] Writing 2D CSV %s (%d rows)\n"), sFileName.c_str(), arr.Num());
    WriteRealArray2(sFileName, arr);
}

//=============================================================================
// Exporter by rank.
// Rank is the nested depth of T inside CMeasureDataT<T>::m_lstData.
//   Rank 0 -> m_lstData is TArray<T>              (1D array)
//   Rank 1 -> m_lstData is TArray<TArray<T>>      (2D array)
//   Rank >= 2 -> 3D or higher, only NPY is exported.
//=============================================================================
template<typename T, INT Rank>
struct TMeasureExporterByRank
{
    static void Export(const TArray<T>& arr, const CCString& sPath)
    {
        appGeneral(_T("[CMeasureData] Export rank-%d NPY to %s.npy (%d top-level elements)\n"),
            Rank, sPath.c_str(), arr.Num());
        SaveAsNumpyFile(sPath + _T(".npy"), arr);
    }
};

template<typename T>
struct TMeasureExporterByRank<T, 0>
{
    static void Export(const TArray<T>& arr, const CCString& sPath)
    {
        appGeneral(_T("[CMeasureData] Export rank-0 (1D scalar) to %s (%d values)\n"),
            sPath.c_str(), arr.Num());
        WriteScalarCsv(sPath + _T(".csv"), arr);
        SaveAsNumpyFile(sPath + _T(".npy"), arr);
    }
};

template<typename T>
struct TMeasureExporterByRank<T, 1>
{
    static void Export(const TArray<T>& arr, const CCString& sPath)
    {
        appGeneral(_T("[CMeasureData] Export rank-1 (2D scalar) to %s (%d rows)\n"),
            sPath.c_str(), arr.Num());
        WriteScalarCsv2(sPath + _T(".csv"), arr);
        SaveAsNumpyFile(sPath + _T(".npy"), arr);
    }
};

template<typename T>
class CMeasureDataT : public CMeasureData
{
public:

    CMeasureDataT() {}

    ~CMeasureDataT() {}

    void Reset() override
    {
        m_lstData.RemoveAll();
    }

    INT Num() const override
    {
        return m_lstData.Num();
    }

    void Export(const CCString& sPath) const override
    {
        appGeneral(_T("[CMeasureDataT] Exporting data (count=%d) to prefix %s\n"),
            m_lstData.Num(), sPath.c_str());
        TMeasureExporterByRank<T, TMeasureRank<T>::Value>::Export(m_lstData, sPath);
    }

    TArray<T>& Get()
    {
        return m_lstData;
    }

    const TArray<T>& Get() const
    {
        return m_lstData;
    }

private:

    TArray<T> m_lstData;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================
