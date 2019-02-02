//=============================================================================
// FILENAME : CLGTest.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#include "CLGTest.h"

TestList* _testSuits;

UINT RunTest(CParameters&params, TestList* pTest)
{
    appGeneral("=========== Testing:%s \n", pTest->m_sParamName);

    CParameters paramForTheTest = params.GetParameter(pTest->m_sParamName);

    //Initial
    if (!appInitialCLG(paramForTheTest))
    {
        return 1;
    }

    //Do the work
    CTimer timer;
    timer.Start();
    UINT uiErrors = (*pTest->m_pfTest)(paramForTheTest);
    timer.Stop();
    appGeneral(_T("=========== Finished, errors: %d, cost: %f(ms)\n"), uiErrors, timer.Elapsed());

    //Final
    appQuitCLG();

    return uiErrors;
}

void ListAllTests(const TArray<TestList*>& allTests)
{
    for (INT i = 0; i <= (allTests.Num() + 1) / 4; ++i)
    {
        for (INT j = 0; j < 4; ++j)
        {
            INT indexOfTest = i * 4 + j;
            if (indexOfTest < allTests.Num())
            {
                TCHAR names[256];
                sprintf_s(names, _T("%d - %s    "), indexOfTest + 1, allTests[indexOfTest]->m_sParamName);
                COUT << names;
            }
            else if (indexOfTest == allTests.Num())
            {
                COUT << "0 - All";
            }
        }
        COUT << std::endl;
    }
}

int main(int argc, char * argv[])
{
    //Load settings
    CParameters params;
    CYAMLParser::ParseFile(_T("TestSuit.yaml"), params);
    appSetupLog(params);
    UINT uiTotalError = 0;

    TArray<TestList*> allTests;
    for (TestList* pTest = _testSuits; NULL != pTest; pTest = pTest->m_pNext)
    {
        if (params.Exist(pTest->m_sParamName))
        {
            allTests.AddItem(pTest);
        }
    }

    UBOOL bCorrectInput = FALSE;
    INT inputNumber = -1;
    while (!bCorrectInput)
    {
        ListAllTests(allTests);
        inputNumber = -1;
        CIN >> inputNumber;
        if (inputNumber >= 0 && inputNumber < allTests.Num() + 1)
        {
            bCorrectInput = TRUE;
        }
    }

    UINT uiError = 0;
    if (0 == inputNumber)
    {
        for (INT i = 0; i < allTests.Num(); ++i)
        {
            uiError += RunTest(params, allTests[i]);
        }
    }
    else
    {
        uiError += RunTest(params, allTests[inputNumber - 1]);
    }

    return uiError;
}

//=============================================================================
// END OF FILE
//=============================================================================
