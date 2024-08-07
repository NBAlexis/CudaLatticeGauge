//=============================================================================
// FILENAME : CLGTest.h
// 
// DESCRIPTION:
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#include "CLGLib.h"

#define ___REGIST_TEST(functionname, catogary, paramName, showName, tag) \
struct STestSuits##functionname##paramName : public TestList \
{ \
    STestSuits##functionname##paramName(testfunction pf) \
    { \
        m_uiIndex = 0; \
        m_uiTag = tag; \
        m_pfTest = pf; \
        m_sCatogary = _T(#catogary); \
        m_sParamName = _T(#paramName); \
        m_sShowName = _T(#showName); \
        Link(_testSuits); \
    } \
}; \
static STestSuits##functionname##paramName registTest##functionname##paramName(functionname);

#define __REGIST_TEST(functionname, catogary, paramName, showName) ___REGIST_TEST(functionname, catogary, paramName, showName, 0)

typedef UINT (*testfunction)(CParameters& sParamName);

enum
{
    _TEST_BOUND = 0x01,
    _TEST_DOUBLE = 0x02,
    _TEST_RELEASE = 0x04,
    _TEST_SINGLE = 0x08,
    _TEST_NOCHECK = 0x10,
};

struct STestSuits
{
    UINT m_uiIndex;
    UINT m_uiTag;
    testfunction m_pfTest;
    const TCHAR* m_sCatogary;
    const TCHAR* m_sParamName;
    const TCHAR* m_sShowName;

    UBOOL OnlyBound() const
    {
        return 0 != (m_uiTag & static_cast<UINT>(_TEST_BOUND));
    }

    UBOOL OnlyDouble() const
    {
        return 0 != (m_uiTag & static_cast<UINT>(_TEST_DOUBLE));
    }

    UBOOL OnlyRelease() const
    {
        return 0 != (m_uiTag & static_cast<UINT>(_TEST_RELEASE));
    }

    UBOOL OnlySingle() const
    {
        return 0 != (m_uiTag & static_cast<UINT>(_TEST_SINGLE));
    }

    UBOOL NoCheck() const
    {
        return 0 != (m_uiTag & static_cast<UINT>(_TEST_NOCHECK));
    }

    CCString GetName() const
    {
        CCString sTag;
        if (OnlyBound())
        {
            sTag = sTag + _T("B");
        }
        if (OnlyDouble())
        {
            sTag = sTag + _T("D");
        }
        if (OnlyRelease())
        {
            sTag = sTag + _T("R");
        }
        if (OnlySingle())
        {
            sTag = sTag + _T("S");
        }
        if (NoCheck())
        {
            sTag = sTag + _T("V");
        }

        if (sTag.GetLength() > 0)
        {
            return sTag + _T("_") + m_sShowName;
        }
        return m_sShowName;
    }
};

typedef TSimpleDoubleLinkedList<STestSuits> TestList;

extern TestList* _testSuits;
extern CCString _bug;

inline void LastProbem(const CCString& problem)
{
    _bug = _bug + problem;
}

extern CCString _msg;
inline void AddMsg(const CCString& msg)
{
    _msg = _msg + msg;
}

//=============================================================================
// Switches which turn on or turn off tests
//=============================================================================

//=============================================================================
// Common used functions
//=============================================================================
extern UINT TestUpdateCommon(CParameters& sParam);

//=============================================================================
// END OF FILE
//=============================================================================
