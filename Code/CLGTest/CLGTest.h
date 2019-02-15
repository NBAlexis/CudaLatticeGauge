//=============================================================================
// FILENAME : CLGTest.h
// 
// DESCRIPTION:
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#include "CLGLib.h"

#define __REGIST_TEST(functionname, catogary, paramName) \
struct STestSuits##paramName : public TestList \
{ \
    STestSuits##paramName(testfunction pf) \
    { \
        m_uiIndex = 0; \
        m_pfTest = pf; \
        m_sCatogary = _T(#catogary); \
        m_sParamName = _T(#paramName); \
        Link(_testSuits); \
    } \
}; \
static STestSuits##paramName registTest##paramName(functionname); 


typedef UINT (*testfunction)(CParameters& sParamName);

struct STestSuits
{
    UINT m_uiIndex;
    testfunction m_pfTest;
    const TCHAR* m_sCatogary;
    const TCHAR* m_sParamName;
};

typedef TSimpleDoubleLinkedList<STestSuits> TestList;

extern TestList* _testSuits;

//=============================================================================
// END OF FILE
//=============================================================================
