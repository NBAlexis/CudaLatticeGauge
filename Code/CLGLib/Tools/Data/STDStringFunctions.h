//=============================================================================
// FILENAME : STDStringFunctions.h
// 
// DESCRIPTION:
//
// Use C++11 std:: string functions instead of <windows.h>, prepare for the Ubuntu build
//
// REVISION:
//  [03/02/2019 nbale]
//=============================================================================

#ifndef _STLSTRINGFUNCTIONS_H_
#define _STLSTRINGFUNCTIONS_H_

#   define appStrpbrk   std::strpbrk
#   define appStrstr    std::strstr
#   define appStrcpy    strcpy_s
#   define appStrcat    std::strcat
#   define appStrlen    std::strlen
#   define appStrcmp    std::strcmp
#   define appStrncmp   std::strncmp
#   define appStrchr    std::strchr
#   define appStrrchr    std::strrchr
#   define appSprintf    sprintf_s
#   define appVsprintf    vsprintf_s
#   define appVsnprintf    vsnprintf_s

#   define appIsSpace    ::isspace
#   define appIsDigit    ::isdigit

__BEGIN_NAMESPACE

#if !_CLG_WIN

inline INT sprintf_s(TCHAR* buffer, size_t sizeOfBuffer, const TCHAR* format, ...)
{
    va_list ap;
    va_start(ap, format);
    INT result = std::vsnprintf(buffer, sizeOfBuffer, format, ap);
    va_end(ap);
    return result;
}

inline INT vsprintf_s(TCHAR* buffer, size_t sizeOfBuffer, const TCHAR* format, va_list ap)
{
    return std::vsnprintf(buffer, sizeOfBuffer, format, ap);
}

inline INT vsnprintf_s(TCHAR* buffer, size_t sizeOfBuffer, const TCHAR* format, va_list ap)
{
    return std::vsnprintf(buffer, sizeOfBuffer, format, ap);
}

inline TCHAR* strcpy_s(TCHAR* dest, size_t length, const TCHAR* source)
{
    size_t lengthOfdest = static_cast<size_t>(sizeof(dest) / sizeof(TCHAR));
    if (lengthOfdest < length)
        length = lengthOfdest;
    return std::strncpy(dest, source, length);
}

#endif

inline TCHAR appToUpper(TCHAR in)
{
    return static_cast<TCHAR>(::toupper(in));
}

inline TCHAR appToLower(TCHAR in)
{
    return static_cast<TCHAR>(::tolower(in));
}

inline INT appStoI(const TCHAR* str, INT iBase = 10)
{
    try { return std::stoi(STDSTRING(str), 0, iBase); }
    catch (...)
    {
        return 0;
    }
}

inline UINT appStoUI(const TCHAR* str, INT iBase = 10)
{
    //ulong is uint in MSVC, but not in GCC
    try { return static_cast<UINT>(std::stoul(STDSTRING(str), 0, iBase)); }
    catch (...)
    {
        return 0;
    }
}

inline LONGLONG appStoLL(const TCHAR* str, INT iBase = 10)
{
    try { return std::stoll(STDSTRING(str), 0, iBase); }
    catch (...)
    {
        return 0;
    }
}

inline ULONGLONG appStoULL(const TCHAR* str, INT iBase = 10)
{
    try { return std::stoull(STDSTRING(str), 0, iBase); }
    catch (...)
    {
        return 0;
    }
}

inline FLOAT appStoF(const TCHAR* str)
{
    
    try { return std::stof(STDSTRING(str)); }
    catch (...)
    {
        return 0.0f;
    }
}

inline DOUBLE appStoD(const TCHAR* str)
{
    
    try { return std::stod(STDSTRING(str)); }
    catch (...)
    {
        return 0.0;
    }
}

inline void appStrupr(TCHAR* buff, size_t bufferL, const TCHAR* s)
{
    STDSTRING str(s);
    std::transform(str.begin(), str.end(), str.begin(), appToUpper);
    appStrcpy(buff, bufferL, str.c_str());
}

inline void appStrlwr(TCHAR* buff, size_t bufferL, const TCHAR* s)
{
    STDSTRING str(s);
    std::transform(str.begin(), str.end(), str.begin(), appToLower);
    appStrcpy(buff, bufferL, str.c_str());
}

inline INT appStricmp(const TCHAR* a, const TCHAR* b)
{
    STDSTRING stra(a);
    std::transform(stra.begin(), stra.end(), stra.begin(), appToUpper);
    STDSTRING strb(b);
    std::transform(strb.begin(), strb.end(), strb.begin(), appToUpper);

    return appStrcmp(stra.c_str(), strb.c_str());
}

inline TCHAR* appStrInc(const TCHAR* a)
{
    return const_cast<TCHAR*>(a + 1);
}

inline void appStrRev(TCHAR* s)
{
    std::string str(s);
    std::reverse(str.begin(), str.end());
    appStrcpy(s, appStrlen(s), str.c_str());
}

__END_NAMESPACE

#endif //#ifndef _STLSTRINGFUNCTIONS_H_

//=============================================================================
// END OF FILE
//=============================================================================