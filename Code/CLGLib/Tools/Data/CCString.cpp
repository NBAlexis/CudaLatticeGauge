//=============================================================================
// FILENAME : CCString.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [3/13/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

//============================================================
//	String
//============================================================
static TCHAR _NullChar = _T('\0');
static INT _NullString[] = {-1,0,0,0};
static CCStringData* _EmptyStringData = (CCStringData*)&_NullString;
CLGAPI const TCHAR* __GEmptyString = (const TCHAR*)((BYTE*)&_NullString+sizeof(CCStringData));

/**
*
*
*/
CCString::CCString(const CCString& stringSrc)
{
    assert(stringSrc.GetData()->m_nRefs != 0);
    if (stringSrc.GetData()->m_nRefs >= 0)
    {
        assert(stringSrc.GetData() != _EmptyStringData);
        m_pchData = stringSrc.m_pchData;
        //appInterlockedIncrement(&GetData()->m_nRefs);
        GetData()->m_nRefs++;
    }
    else
    {
        Init();
        *this = stringSrc.m_pchData;
    }
}

/**
*
*
*/
CCString::CCString(const TCHAR* lpsz)
{
    Init();
    SIZE_T nLen = (SIZE_T)(__SafeStrlen(lpsz));
    if (nLen != 0)
    {
        AllocBuffer((INT)(nLen));
        memcpy(m_pchData, lpsz, nLen*sizeof(TCHAR));
    }
}

/**
*
*
*/
CCString::CCString(TCHAR ch, INT nLength)
{
    Init();
    if (nLength >= 1)
    {
        AllocBuffer(nLength);
        for (INT i = 0; i < nLength; ++i) //use appMemset instead?
            m_pchData[i] = ch;
    }
}

/**
*
*
*/
CCString::CCString(const TCHAR* lpch, INT nLength)
{
    Init();
    if (nLength != 0)
    {
        AllocBuffer(nLength);
        memcpy(m_pchData, lpch, nLength*sizeof(TCHAR));
    }
}

/**
* free any attached data
*
*/
CCString::~CCString()
{
    if (GetData() != _EmptyStringData)
    {
        GetData()->m_nRefs--;
        //if (appInterlockedDecrement(&GetData()->m_nRefs) <= 0)
        if (GetData()->m_nRefs <= 0)
            FreeData(GetData());
    }
}

/**
*
*
*/
void CCString::Release()
{
    if (GetData() != _EmptyStringData)
    {
        assert(GetData()->m_nRefs != 0);
        GetData()->m_nRefs--;
        //if (appInterlockedDecrement(&GetData()->m_nRefs) <= 0)
        if (GetData()->m_nRefs <= 0)
            FreeData(GetData());
        Init();
    }
}

/**
*
*
*/
void CCString::Release(CCStringData* pData)
{
    if (pData != _EmptyStringData)
    {
        assert(pData->m_nRefs != 0);
        pData->m_nRefs--;
        //if (appInterlockedDecrement(&pData->m_nRefs) <= 0)
        if (pData->m_nRefs <= 0)
            FreeData(pData);
    }
}

/**
*
*
*/
void CCString::Empty()
{
    if (GetData()->m_nDataLength == 0)
        return;
    if (GetData()->m_nRefs >= 0)
        Release();
    else
        *this = &_NullChar;
    assert(GetData()->m_nDataLength == 0);
    assert(GetData()->m_nRefs < 0 || GetData()->m_nAllocLength == 0);
}

/**
*
*
*/
const CCString& CCString::operator=(const CCString& stringSrc)
{
    if (m_pchData != stringSrc.m_pchData)
    {
        if ((GetData()->m_nRefs < 0 && GetData() != _EmptyStringData) ||
            stringSrc.GetData()->m_nRefs < 0)
        {
            // actual copy necessary since one of the strings is locked
            AssignCopy(stringSrc.GetData()->m_nDataLength, stringSrc.m_pchData);
        }
        else
        {
            // can just copy references around
            Release();
            assert(stringSrc.GetData() != _EmptyStringData);
            m_pchData = stringSrc.m_pchData;
            //appInterlockedIncrement(&GetData()->m_nRefs);
            GetData()->m_nRefs++;
        }
    }
    return *this;
}

//////////////////////////////////////////////////////////////////////////////
// concatenation
// NOTE: "operator+" is done as friend functions for simplicity
//      There are three variants:
//          String + String
// and for ? = TCHAR, const TCHAR*
//          String + ?
//          ? + String

/**
*
*
*/
CCString CLGAPI operator+(const CCString& string1, const CCString& string2)
{
    CCString s;
    s.ConcatCopy(string1.GetData()->m_nDataLength, string1.m_pchData,
        string2.GetData()->m_nDataLength, string2.m_pchData);
    return s;
}

/**
*
*
*/
CCString CLGAPI operator+(const CCString& string, const TCHAR* lpsz)
{
    CCString s;
    s.ConcatCopy(string.GetData()->m_nDataLength, string.m_pchData,
        __SafeStrlen(lpsz), lpsz);
    return s;
}

/**
*
*
*/
CCString CLGAPI operator+(const TCHAR* lpsz, const CCString& string)
{
    CCString s;
    s.ConcatCopy(__SafeStrlen(lpsz), lpsz, string.GetData()->m_nDataLength,
        string.m_pchData);
    return s;
}

/**
*
*
*/
CCString CLGAPI operator+(const CCString& string1, TCHAR ch)
{
    CCString s;
    s.ConcatCopy(string1.GetData()->m_nDataLength, string1.m_pchData, 1, &ch);
    return s;
}

/**
*
*
*/
CCString CLGAPI operator+(TCHAR ch, const CCString& string)
{
    CCString s;
    s.ConcatCopy(1, &ch, string.GetData()->m_nDataLength, string.m_pchData);
    return s;
}

///////////////////////////////////////////////////////////////////////////////
// Advanced direct buffer access

/**
* Do not use this please
*
*/
TCHAR* CCString::GetBuffer(INT nMinBufLength)
{
    assert(nMinBufLength >= 0);

    if (GetData()->m_nRefs > 1 || nMinBufLength > GetData()->m_nAllocLength)
    {
#ifdef _DEBUG
        // give a warning in case locked string becomes unlocked
        if (GetData() != _EmptyStringData && GetData()->m_nRefs < 0)
            appGeneral(_T("Warning: GetBuffer on locked FString creates unlocked FString!\n"));
#endif
        // we have to grow the buffer
        CCStringData* pOldData = GetData();
        INT nOldLen = GetData()->m_nDataLength;   // AllocBuffer will tromp it
        if (nMinBufLength < nOldLen)
            nMinBufLength = nOldLen;
        AllocBuffer(nMinBufLength);
        memcpy(m_pchData, pOldData->Data(), (nOldLen + 1) * sizeof(TCHAR));
        GetData()->m_nDataLength = nOldLen;
        CCString::Release(pOldData);
    }
    assert(GetData()->m_nRefs <= 1);

    // return a pointer to the character storage for this string
    assert(m_pchData != NULL);
    return m_pchData;
}

/**
*
*
*/
void CCString::ReleaseBuffer(INT nNewLength)
{
    CopyBeforeWrite();  // just in case GetBuffer was not called

    if (nNewLength == -1)
        nNewLength = (INT)appStrlen(m_pchData); // zero terminated

    assert(nNewLength <= GetData()->m_nAllocLength);
    GetData()->m_nDataLength = nNewLength;
    m_pchData[nNewLength] = _T('\0');
}

/**
*
*
*/
TCHAR* CCString::GetBufferSetLength(INT nNewLength)
{
    assert(nNewLength >= 0);

    GetBuffer(nNewLength);
    GetData()->m_nDataLength = nNewLength;
    m_pchData[nNewLength] = _T('\0');
    return m_pchData;
}

/**
*
*
*/
void CCString::FreeExtra()
{
    assert(GetData()->m_nDataLength <= GetData()->m_nAllocLength);
    if (GetData()->m_nDataLength != GetData()->m_nAllocLength)
    {
        CCStringData* pOldData = GetData();
        AllocBuffer(GetData()->m_nDataLength);
        memcpy(m_pchData, pOldData->Data(), pOldData->m_nDataLength * sizeof(TCHAR));
        assert(_T('\0') == m_pchData[GetData()->m_nDataLength]);
        CCString::Release(pOldData);
    }
    assert(GetData() != NULL);
}

/**
*
*
*/
void CCString::UnlockBuffer()
{
    assert(GetData()->m_nRefs == -1);
    if (GetData() != _EmptyStringData)
        GetData()->m_nRefs = 1;
}

/**
* find first non-space character
*
*/
void CCString::TrimLeft()
{
    CopyBeforeWrite();
    const TCHAR* lpsz = m_pchData;

    while (appIsSpace(*lpsz))
        lpsz = appStrInc(lpsz);

    if (lpsz != m_pchData)
    {
        // fix up data and length
        INT nDataLength = GetData()->m_nDataLength - (INT)(lpsz - m_pchData);
        memmove(m_pchData, lpsz, (nDataLength + 1) * sizeof(TCHAR));
        GetData()->m_nDataLength = nDataLength;
    }
}

/**
*
*
*/
void __cdecl CCString::Format(const TCHAR* lpszFormat, ...)
{
    va_list argList;
    va_start(argList, lpszFormat);
    FormatV(lpszFormat, argList);
    va_end(argList);
}

CCString CCString::FormatS(const TCHAR* lpszFormat, ...)
{
    CCString ret;
    va_list argList;
    va_start(argList, lpszFormat);
    ret = FormatVS(lpszFormat, argList);
    va_end(argList);
    return ret;
}

CCString CCString::FormatVS(const TCHAR* lpszFormat, va_list argList)
{
    CCString ret;
    ret.FormatV(lpszFormat, argList);
    return ret;
}

#define TCHAR_ARG   TCHAR
#define DOUBLE_ARG  DOUBLE

/**
*
*
*/
void CCString::FormatV(const TCHAR* lpszFormat, va_list argList)
{
    va_list argListSave = argList;

    // make a guess at the maximum length of the resulting string
    INT nMaxLen = 0;
    for (const TCHAR* lpsz = lpszFormat; *lpsz != _T('\0'); lpsz = appStrInc(lpsz))
    {
        // handle '%' character, but watch out for '%%'
        if (*lpsz != _T('%') || *(lpsz = appStrInc(lpsz)) == _T('%'))
        {
            nMaxLen += (INT)appStrlen(lpsz);
            continue;
        }

        INT nItemLen = 0;

        // handle '%' character with format
        INT nWidth = 0;
        for (; *lpsz != _T('\0'); lpsz = appStrInc(lpsz))
        {
            // check for valid flags
            if (*lpsz == _T('#'))
                nMaxLen += 2;   // for '0x'
            else if (*lpsz == _T('*'))
                nWidth = va_arg(argList, INT);
            else if (*lpsz == _T('-') || *lpsz == _T('+') || *lpsz == _T('0') ||
                *lpsz == _T(' '))
                ;
            else // hit non-flag character
                break;
        }
        // get width and skip it
        if (nWidth == 0)
        {
            // width indicated by
            if (appIsDigit(*lpsz))
            {
                nWidth = appStoI(lpsz);
                for (; *lpsz != _T('\0') && appIsDigit(*lpsz); lpsz = appStrInc(lpsz))
                    ;
            }
        }
        assert(nWidth >= 0);

        INT nPrecision = 0;
        if (*lpsz == _T('.'))
        {
            // skip past '.' separator (width.precision)
            lpsz = appStrInc(lpsz);

            // get precision and skip it
            if (*lpsz == _T('*'))
            {
                nPrecision = va_arg(argList, INT);
                lpsz = appStrInc(lpsz);
            }
            else
            {
                nPrecision = appStoI(lpsz);
                for (; *lpsz != _T('\0') && appIsDigit(*lpsz); lpsz = appStrInc(lpsz))
                    ;
            }
            assert(nPrecision >= 0);
        }

        // now should be on specifier
        switch (*lpsz)
        {
            // single characters
        case _T('c'):
        case _T('C'):
            nItemLen = 2;
            va_arg(argList, TCHAR_ARG);
            break;

            // strings
        case _T('s'):
        case _T('S'):
            {
                const TCHAR* pstrNextArg = va_arg(argList, const TCHAR*);
                if (pstrNextArg == NULL)
                    nItemLen = 6;  // "(null)"
                else
                {
                    nItemLen = (INT)appStrlen(pstrNextArg);
                    nItemLen = appMax(1, nItemLen);
                }
            }
            break;
        }

        // adjust nItemLen for strings
        if (nItemLen != 0)
        {
            if (nPrecision != 0)
                nItemLen = appMin(nItemLen, nPrecision);
            nItemLen = appMax(nItemLen, nWidth);
        }
        else
        {
            switch (*lpsz)
            {
                // integers
            case _T('d'):
            case _T('i'):
            case _T('u'):
            case _T('x'):
            case _T('X'):
            case _T('o'):
                va_arg(argList, INT);
                nItemLen = 32;
                nItemLen = appMax(nItemLen, nWidth+nPrecision);
                break;

            case _T('e'):
            case _T('g'):
            case _T('G'):

                va_arg(argList, DOUBLE_ARG);
                nItemLen = 128;
                nItemLen = appMax(nItemLen, nWidth+nPrecision);
                break;

            case _T('f'):
                {
                    DOUBLE f;
                    TCHAR* pszTemp;

                    // 312 == strlen("-1+(309 zeroes).")
                    // 309 zeroes == max precision of a double
                    // 6 == adjustment in case precision is not specified,
                    //   which means that the precision defaults to 6
                    DWORD nLength = appMax(nWidth, 312 + nPrecision + 6);
                    pszTemp = (TCHAR*)appAlloca(nLength);

                    f = va_arg(argList, DOUBLE);
                    appSprintf( pszTemp, nLength, _T( "%*.*f" ), nWidth, nPrecision+6, f );
                    nItemLen = (INT)appStrlen(pszTemp);
                }
                break;

            case _T('p'):
                va_arg(argList, void*);
                nItemLen = 32;
                nItemLen = appMax(nItemLen, nWidth+nPrecision);
                break;

                // no output
            case _T('n'):
                va_arg(argList, INT*);
                break;

            default:
                assert(FALSE);  // unknown formatting option
            }
        }

        // adjust nMaxLen for output nItemLen
        nMaxLen += nItemLen;
    }

    GetBuffer(nMaxLen);
    appVsprintf(m_pchData, GetAllocLength(), lpszFormat, argListSave);
    ReleaseBuffer();

    va_end(argListSave);
}

//////////////////////////////////////////////////////////////////////////////
// Advanced manipulation

/**
*
*
*/
INT CCString::Delete(INT nIndex, INT nCount /* = 1 */)
{
    if (nIndex < 0)
        nIndex = 0;
    INT nNewLength = GetData()->m_nDataLength;
    if (nCount > 0 && nIndex < nNewLength)
    {
        CopyBeforeWrite();
        INT nBytesToCopy = nNewLength - (nIndex + nCount) + 1;

        memcpy(m_pchData + nIndex,
            m_pchData + nIndex + nCount, nBytesToCopy * sizeof(TCHAR));
        GetData()->m_nDataLength = nNewLength - nCount;
    }

    return nNewLength;
}

/**
*
*
*/
INT CCString::Insert(INT nIndex, TCHAR ch)
{
    CopyBeforeWrite();

    if (nIndex < 0)
        nIndex = 0;

    INT nNewLength = GetData()->m_nDataLength;
    if (nIndex > nNewLength)
        nIndex = nNewLength;
    ++nNewLength;

    if (GetData()->m_nAllocLength < nNewLength)
    {
        CCStringData* pOldData = GetData();
        TCHAR* pstr = m_pchData;
        AllocBuffer(nNewLength);
        memcpy(m_pchData, pstr, (pOldData->m_nDataLength + 1) * sizeof(TCHAR));
        CCString::Release(pOldData);
    }

    // move existing bytes down
    //memcpy(m_pchData + nIndex + 1,
    //	m_pchData + nIndex, (nNewLength-nIndex)*sizeof(TCHAR));
    memmove(m_pchData + nIndex + 1, m_pchData + nIndex, (nNewLength - nIndex) * sizeof(TCHAR));
    m_pchData[nIndex] = ch;
    GetData()->m_nDataLength = nNewLength;

    return nNewLength;
}

/**
*
*
*/
INT CCString::Insert(INT nIndex, const TCHAR* pstr)
{
    if (nIndex < 0)
        nIndex = 0;

    INT nInsertLength = __SafeStrlen(pstr);
    INT nNewLength = GetData()->m_nDataLength;
    if (nInsertLength > 0)
    {
        CopyBeforeWrite();
        if (nIndex > nNewLength)
            nIndex = nNewLength;
        nNewLength += nInsertLength;

        if (GetData()->m_nAllocLength < nNewLength)
        {
            CCStringData* pOldData = GetData();
            TCHAR* pstrd = m_pchData;
            AllocBuffer(nNewLength);
            memcpy(m_pchData, pstrd, (pOldData->m_nDataLength + 1) * sizeof(TCHAR));
            CCString::Release(pOldData);
        }

        // move existing bytes down
        //memcpy(m_pchData + nIndex + nInsertLength,
        //	m_pchData + nIndex,
        //	(nNewLength-nIndex-nInsertLength+1)*sizeof(TCHAR));
        memmove(m_pchData + nIndex + nInsertLength, m_pchData + nIndex, (nNewLength - nIndex - nInsertLength + 1) * sizeof(TCHAR));
        memcpy(m_pchData + nIndex, pstr, nInsertLength * sizeof(TCHAR));
        GetData()->m_nDataLength = nNewLength;
    }

    return nNewLength;
}

/**
*
*
*/
INT CCString::Replace(TCHAR chOld, TCHAR chNew)
{
    INT nCount = 0;

    // short-circuit the nop case
    if (chOld != chNew)
    {
        // otherwise modify each character that matches in the string
        CopyBeforeWrite();
        TCHAR* psz = m_pchData;
        TCHAR* pszEnd = psz + GetData()->m_nDataLength;
        while (psz < pszEnd)
        {
            // replace instances of the specified character only
            if (*psz == chOld)
            {
                *psz = chNew;
                ++nCount;
            }
            psz = appStrInc(psz);
        }
    }
    return nCount;
}

/**
*
*
*/
INT CCString::Replace(const TCHAR* lpszOld, const TCHAR* lpszNew)
{
    // can't have empty or NULL lpszOld
    INT nSourceLen = __SafeStrlen(lpszOld);
    if (nSourceLen == 0)
        return 0;
    INT nReplacementLen = __SafeStrlen(lpszNew);

    // loop once to figure out the size of the result string
    INT nCount = 0;
    TCHAR* lpszStart = m_pchData;
    TCHAR* lpszEnd = m_pchData + GetData()->m_nDataLength;
    TCHAR* lpszTarget;
    while (lpszStart < lpszEnd)
    {
        while ((lpszTarget = appStrstr(lpszStart, lpszOld)) != NULL)
        {
            ++nCount;
            lpszStart = lpszTarget + nSourceLen;
        }
        lpszStart += appStrlen(lpszStart) + 1; //appStrlen different from lstrlen?
    }

    // if any changes were made, make them
    if (nCount > 0)
    {
        CopyBeforeWrite();

        // if the buffer is too small, just
        //   allocate a new buffer (slow but sure)
        INT nOldLength = GetData()->m_nDataLength;
        INT nNewLength =  nOldLength + (nReplacementLen - nSourceLen) * nCount;
        if (GetData()->m_nAllocLength < nNewLength || GetData()->m_nRefs > 1)
        {
            CCStringData* pOldData = GetData();
            TCHAR* pstr = m_pchData;
            AllocBuffer(nNewLength);
            memcpy(m_pchData, pstr, pOldData->m_nDataLength * sizeof(TCHAR));
            CCString::Release(pOldData);
        }
        // else, we just do it in-place
        lpszStart = m_pchData;
        lpszEnd = m_pchData + GetData()->m_nDataLength;

        // loop again to actually do the work
        while (lpszStart < lpszEnd)
        {
            while ( (lpszTarget = appStrstr(lpszStart, lpszOld)) != NULL)
            {
                INT nBalance = nOldLength - (INT)(lpszTarget - m_pchData) + nSourceLen;
                memmove(lpszTarget + nReplacementLen, lpszTarget + nSourceLen,
                    nBalance * sizeof(TCHAR));
                memcpy(lpszTarget, lpszNew, nReplacementLen * sizeof(TCHAR));
                lpszStart = lpszTarget + nReplacementLen;
                lpszStart[nBalance] = _T('\0');
                nOldLength += (nReplacementLen - nSourceLen);
            }
            lpszStart += appStrlen(lpszStart) + 1; //appStrlen different from lstrlen?
        }
        assert(m_pchData[nNewLength] == _T('\0'));
        GetData()->m_nDataLength = nNewLength;
    }

    return nCount;
}

/**
*
*
*/
INT CCString::Remove(TCHAR chRemove)
{
    CopyBeforeWrite();

    TCHAR* pstrSource = m_pchData;
    TCHAR* pstrDest = m_pchData;
    TCHAR* pstrEnd = m_pchData + GetData()->m_nDataLength;

    while (pstrSource < pstrEnd)
    {
        if (*pstrSource != chRemove)
        {
            *pstrDest = *pstrSource;
            pstrDest = appStrInc(pstrDest);
        }
        pstrSource = appStrInc(pstrSource);
    }
    *pstrDest = _T('\0');
    INT nCount = (INT)(pstrSource - pstrDest);
    GetData()->m_nDataLength -= nCount;

    return nCount;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
