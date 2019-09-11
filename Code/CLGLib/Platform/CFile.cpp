//=============================================================================
// FILENAME : CFileWin.cpp
// 
// DESCRIPTION:
//  Implement of CFileSytem in Windows
//
// REVISION:
//  [01/31/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#if !_CLG_WIN
void fopen_s(FILE** fp, const TCHAR* fileName, const TCHAR* fileMode)
{
    (*fp) = fopen(fileName, fileMode);
}
#endif

BYTE* CFileSystem::ReadAllBytes(const TCHAR* sFilename, UINT& size)
{
    FILE* fp = NULL;
    fopen_s(&fp, sFilename, _T("rb"));

    if (NULL == fp)
    {
        return NULL;
    }
    fseek(fp, 0, SEEK_END);
    const UINT fileSize = static_cast<UINT>(ftell(fp));
    size = fileSize;
    fseek(fp, 0, SEEK_SET);
    BYTE* ret = (BYTE*)malloc(fileSize);
    if (NULL == ret)
    {
        return NULL;
    }
#if defined(__GNUC__)
    size_t notusing __attribute__((unused)) = fread(ret, 1, fileSize, fp);
#else
    fread(ret, 1, fileSize, fp);
#endif
    fclose(fp);
    return ret;
}

UBOOL CFileSystem::WriteAllBytes(const TCHAR* sFilename, BYTE* data, UINT uiSize)
{
    FILE* fp = NULL;
    fopen_s(&fp, sFilename, _T("wb"));
    if (NULL == fp)
    {
        return FALSE;
    }

    fwrite(data, 1, uiSize, fp);
    fflush(fp);
    fclose(fp);
    return TRUE;
}

CCString CFileSystem::ReadAllText(const TCHAR* sFilename)
{
    FILE* fp = NULL;
    fopen_s(&fp, sFilename, _T("rb"));
    if (NULL == fp)
    {
        return CCString("");
    }
    fseek(fp, 0, SEEK_END);
    const UINT fileSize = static_cast<UINT>(ftell(fp));
    fseek(fp, 0, SEEK_SET);
    BYTE* file_data = (BYTE*)malloc(fileSize + sizeof(TCHAR));

    //cast to void dose not help
#if defined(__GNUC__)
    size_t notusing __attribute__((unused)) = fread(file_data, 1, fileSize, fp);
#else
    fread(file_data, 1, fileSize, fp);
#endif
    fclose(fp);

    file_data[fileSize] = 0;
#if _CLG_UNICODE
    file_data[fileSize + 1] = 0;
#endif

    BYTE* string_start = NULL;
    INT string_len = 0;
    const UBOOL ansi_to_unicode = FALSE;
    UBOOL unicode_to_ansi = FALSE;

#if __BIG_ENDIAN
    UBOOL change_byte_order = TRUE;
#else
    UBOOL change_byte_order = FALSE;
#endif

    if (0xff == file_data[0] && 0xfe == file_data[1])
    {
#if !_CLG_UNICODE
        unicode_to_ansi = TRUE;
#endif
            
        string_start = file_data + 2;
        string_len = (fileSize - 2) / 2;
    }
    else if (0xfe == file_data[0] && 0xff == file_data[1])
    {
#if !_CLG_UNICODE
        unicode_to_ansi = TRUE;
#endif
        string_start = file_data + 2;
        string_len = (fileSize - 2) / 2;
        change_byte_order = TRUE;
    }
    else
    {
#if _CLG_UNICODE
        ansi_to_unicode = TRUE;
        string_start = file_data;
        string_len = fileSize;
#else
        string_start = file_data;
        string_len = fileSize;
#endif
    }

    if (change_byte_order)
    {
        appGeneral(_T("We are changing byte order for config file: %s"), sFilename);
        UNICHAR* p = (UNICHAR*)string_start;
        for (INT i = 0; i < string_len; i++, p++)
        {
            BYTE* a = (BYTE*)p;
            BYTE* b = a + 1;
            const BYTE c = *a;
            *a = *b;
            *b = c;
        }
    }

    CCString s;
    if (ansi_to_unicode)
    {
        //Not supported yet
#if 0
        INT buf_len = (string_len + 1) * 2;
        UNICHAR* new_data = (UNICHAR*)(malloc(buf_len));
        INT ret_len = appAnsiToUnicode(new_data, (ANSICHAR*)string_start/*with terminator*/, buf_len/*out buffer bytes*/);
        memset(new_data, 0, buf_len);
        assert(ret_len <= buf_len);
        free(file_data);
        file_data = (BYTE*)new_data;
        string_start = (BYTE*)new_data;
        ret_len /= 2; //bytes => count
        string_len = ret_len - 1; //do not include terminator

        s = (TCHAR*)string_start;
        free(new_data);
#endif
    }
    else if (unicode_to_ansi)
    {
        //Not supported yet
#if 0
        INT buf_len = string_len * 2 + 1;
        ANSICHAR* new_data = (ANSICHAR*)(malloc(buf_len)); //double the output buffer for pure multi-byte system
        memset(new_data, 0, buf_len);
        INT ret_len = appUnicodeToAnsi(new_data, (UNICHAR*)string_start/*with terminator*/, buf_len/*out buffer bytes*/);
        assert(ret_len <= buf_len);
        free(file_data);
        file_data = (BYTE*)new_data;
        string_start = (BYTE*)new_data;
        string_len = ret_len - 1; //do not include terminatror

        s = (TCHAR*)string_start;
        free(new_data);
#endif
    }
    else
    {
        s = (TCHAR*)file_data;
        free(file_data);
    }

    return s;
}

UBOOL CFileSystem::WriteAllText(const TCHAR* sFilename, const CCString& data)
{
    return WriteAllBytes(sFilename, (BYTE*)data.c_str(), data.GetLength() * sizeof(TCHAR));
}

//UBOOL CFileSystem::SetDefaultDirectory(const TCHAR* Filename)
//{
//    return SetCurrentDirectory(Filename) != 0;
//}

//CCString CFileSystem::GetDefaultDirectory()
//{
//#if _CLG_UNICODE
//    TCHAR Buffer[kMaxPathNameLength] = TEXT("");
//    GetCurrentDirectoryW(ARRAY_COUNT(Buffer), Buffer);
//    return CCString(Buffer);
//#endif
//    UNICHAR Buffer[kMaxPathNameLength] = L"";
//    GetCurrentDirectoryW(ARRAY_COUNT(Buffer), Buffer);
//    ANSICHAR retBuff[kMaxPathNameLength];
//    appUnicodeToAnsi(retBuff, Buffer, kMaxPathNameLength);
//    return CCString(retBuff);
//}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================

