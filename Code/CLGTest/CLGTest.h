//=============================================================================
// FILENAME : CLGTest.h
// 
// DESCRIPTION:
//
// REVISION:
//  [mm/dd/yy]
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
        m_bSupportMultiGPU = \
            (0 != ((tag) & static_cast<UINT>(_TEST_MULTIGPU))) ? TRUE : FALSE; \
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
    _TEST_SINGLE = 0x08,        // single-PRECISION only (float build)
    _TEST_NOCHECK = 0x10,
    _TEST_MULTIGPU = 0x20,      // explicit opt-in: this test may run under a multi-rank grid
};

struct STestSuits
{
    UINT m_uiIndex;
    UINT m_uiTag;
    testfunction m_pfTest;
    const TCHAR* m_sCatogary;
    const TCHAR* m_sParamName;
    const TCHAR* m_sShowName;
    UBOOL m_bSupportMultiGPU;   // default FALSE (single-GPU only); TRUE only via the explicit _TEST_MULTIGPU tag

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

    /** Explicitly registered as multi-rank capable (--mg-config / --mg-worker accept only these). */
    UBOOL IsMultiGPUTagged() const
    {
        return 0 != (m_uiTag & static_cast<UINT>(_TEST_MULTIGPU));
    }

    /** Test may run under a multi-rank grid (only _TEST_MULTIGPU opt-in tests). */
    UBOOL SupportMultiGPU() const
    {
        return m_bSupportMultiGPU;
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
        if (IsMultiGPUTagged())
        {
            sTag = sTag + _T("M");
        }

        if (sTag.GetLength() > 0)
        {
            return sTag + _T("_") + m_sShowName;
        }
        return m_sShowName;
    }
};

// Global multi-GPU mode switch. Always FALSE except inside the --mg-worker run
// mode (multi-GPU-improve1.md 3.10): the interactive menu is strictly
// single-process/single-GPU (the 'g' toggle is removed) and batch modes keep
// the identity GpuGrid=[1,1,1,1]. When TRUE, RunTest injects the worker
// --gpu-grid / --device-per-node values into the test's parameter block before
// appInitialCLG for every _TEST_MULTIGPU-tagged test.
extern UBOOL g_bCLGMultiGPU;

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
extern UINT TestHeatbath(CParameters& sParam);

//=============================================================================
// END OF FILE
//=============================================================================
