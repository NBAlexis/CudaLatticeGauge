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

UINT RunTest(CParameters&params, const TestList* pTest)
{
    appQuitCLG();

    appGeneral("\n=========== Testing:%s \n", pTest->m_sParamName);
    CParameters paramForTheTest = params.GetParameter(pTest->m_sParamName);
    appGeneral(_T("============= Parameters =============\n"));
    paramForTheTest.Dump(_T(""));
    //Initial
    if (!appInitialCLG(paramForTheTest))
    {
        return 1;
    }

    //Do the work
    CTimer timer;
    timer.Start();
    const UINT uiErrors = (*pTest->m_pfTest)(paramForTheTest);
    timer.Stop();
    appGeneral(_T("=========== Finished %s, errors: %d, cost: %f(ms)\n ------------- End --------------\n\n"), pTest->m_sParamName, uiErrors, timer.Elapsed());

    //Final
    //appQuitCLG();

    return uiErrors;
}

void ListAllTests(const THashMap<CCString, TArray<TestList*>*>& category)
{
    TArray<CCString> sKeys = category.GetAllKeys();
    UINT uiIdx = 0;
    for (INT k = 0; k < sKeys.Num(); ++k)
    {
        COUT << _T("============== ") << sKeys[k] << _T(" ==============\n");
        TArray<TestList*>* lst = category.GetAt(sKeys[k]); //category[] only work with non-const THashMap
        for (INT i = 0; i <= lst->Num() / 3; ++i)
        {
            for (INT j = 0; j < 4; ++j)
            {
                const INT indexOfTest = i * 4 + j;
                if (indexOfTest < lst->Num())
                {
                    TCHAR names[256];
                    ++uiIdx;
                    lst->GetAt(indexOfTest)->m_uiIndex = uiIdx;
                    appSprintf(names, 256, _T("%d - %s"), uiIdx, appStrWithLen(lst->GetAt(indexOfTest)->GetName(), 25).c_str());
                    COUT << names;
                }
            }
            COUT << std::endl;
        }
    }
}

TestList* GetTest(const TArray<TestList*>& alltest, UINT idx)
{
    for (INT i = 0; i < alltest.Num(); ++i)
    {
        if (idx == alltest[i]->m_uiIndex)
        {
            return alltest[i];
        }
    }
    return NULL;
}

void DeleteAllLists(THashMap<CCString, TArray<TestList*>*>& category)
{
    //delete the lists
    TArray<CCString> sKeys = category.GetAllKeys();
    for (INT i = 0; i < sKeys.Num(); ++i)
    {
        appSafeDelete(category[sKeys[i]]);
    }
}

int main(int argc, char * argv[])
{
    //Load settings
    CParameters params;
#if _CLG_DEBUG
    CYAMLParser::ParseFile(_T("TestSuit.yaml"), params);
    CYAMLParser::ParseFile(_T("TestSuit_Common.yaml"), params);
    CYAMLParser::ParseFile(_T("TestSuit_Random.yaml"), params);
    CYAMLParser::ParseFile(_T("TestSuit_Boundary.yaml"), params);
    CYAMLParser::ParseFile(_T("TestSuit_FileIO.yaml"), params);
    CYAMLParser::ParseFile(_T("TestSuit_FermionMatrix.yaml"), params);
    CYAMLParser::ParseFile(_T("TestSuit_FermionUpdator.yaml"), params);
    CYAMLParser::ParseFile(_T("TestSuit_FermionUpdatorKS.yaml"), params);
    CYAMLParser::ParseFile(_T("TestSuit_GaugeFixing.yaml"), params);
    CYAMLParser::ParseFile(_T("TestSuit_Rotation.yaml"), params);
    CYAMLParser::ParseFile(_T("TestSuit_MiscUpdate.yaml"), params);
    CYAMLParser::ParseFile(_T("TestSuit_Solver.yaml"), params);
    CYAMLParser::ParseFile(_T("TestSuit_Updator.yaml"), params);
    CYAMLParser::ParseFile(_T("TestSuit_Boson.yaml"), params);
    //CYAMLParser::ParseFile(_T("TestSuit_EvenOdd.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/TestSuit.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Common.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Random.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Boundary.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_FileIO.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_FermionMatrix.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_FermionUpdator.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_FermionUpdatorKS.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_GaugeFixing.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Rotation.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_MiscUpdate.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Solver.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Updator.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Boson.yaml"), params);
    //CYAMLParser::ParseFile(_T("../Debug/TestSuit_EvenOdd.yaml"), params);
#endif
    appSetupLog(params);

    BYTE realByte[8];
    Real testReal = F(-1.2345);
    memset(realByte, 0, 8);
    memcpy(realByte, &testReal, sizeof(Real));
    CCString sRealByte;
    for (UINT i = 0; i < 8; ++i)
    {
        sRealByte += appToString(realByte[i]) + _T(", ");
    }
    appGeneral(_T(" \n ================= sizeof(Real) : %d and -1.2345 is %s =============\n"), sizeof(Real), sRealByte.c_str());
    
    TArray<TestList*> allTests;
    THashMap<CCString, TArray<TestList*>*> category;
    //UINT uiIndex = 0;
    for (TestList* pTest = _testSuits; NULL != pTest; pTest = pTest->m_pNext)
    {
        if (params.Exist(pTest->m_sParamName))
        {
            //pTest->m_uiIndex = uiIndex;
            //++uiIndex;
            allTests.AddItem(pTest);
            CCString sCategory = pTest->m_sCatogary;
            if (category.Exist(sCategory))
            {
                category[sCategory]->AddItem(pTest);
            }
            else
            {
                TArray<TestList*>* newList = new TArray<TestList*>();
                newList->AddItem(pTest);
                category.SetAt(sCategory, newList);
            }
        }
    }

    //INT inputNumber = -1;
    ListAllTests(category);
    while (TRUE)
    {
        COUT << _T("============== CLG v") << GetCLGVersion().c_str() << _T("==============\nq - Quit,  l - List all,  r - Run all,  c - Check all,  d - Device info,  i - lattice info\n");
        //ListAllTests(category);
        //inputNumber = -1;
        std::string name;
        std::getline(std::cin, name);
        CCString sRes(name.c_str());
        INT number = appStrToINT(sRes);
        UBOOL bExcuted = FALSE;
        sRes.MakeLower();
        if (sRes == _T("q"))
        {
            break;
        }

        if (sRes == _T("l"))
        {
            ListAllTests(category);
            bExcuted = TRUE;
        }
        else if (sRes == _T("d"))
        {
            CCudaHelper::DeviceQuery();
            bExcuted = TRUE;
        }
        else if (sRes == _T("i"))
        {
            appGeneral(_T("\n=======================================\n") + appGetLattice()->GetInfos(_T("")) + _T("\n=======================================\n"));
            bExcuted = TRUE;
        }
        else if (appToString(number) == sRes)
        {
            if (number > 0 && number <= allTests.Num())
            {
                TestList* pTest = GetTest(allTests, static_cast<UINT>(number));
                UBOOL bPass = FALSE;
#if _CLG_DEBUG
                if (pTest->OnlyRelease())
                {
                    bPass = TRUE;
                    appGeneral(_T("Cannot run this test in DEBUG mode\n"));
                }
#endif

#if !_CLG_USE_LAUNCH_BOUND
                if (pTest->OnlyBound())
                {
                    bPass = TRUE;
                    appGeneral(_T("Can ONLY run this test in BOUND mode\n"));
                }
#endif

#if !_CLG_DOUBLEFLOAT
                if (pTest->OnlyDouble())
                {
                    bPass = TRUE;
                    appGeneral(_T("Can ONLY this test in DOUBLE mode\n"));
                }
#else
                if (pTest->OnlySingle())
                {
                    bPass = TRUE;
                    appGeneral(_T("Can ONLY this test in FLOAT mode\n"));
                }
#endif
                if (!bPass)
                {
                    RunTest(params, pTest);
                    bExcuted = TRUE;
                }
            }
        }
        else if (sRes == _T("r"))
        {
            CTimer timer;
            timer.Start();
            UINT uiError = 0;
            UINT uiPassed = 0;
            TArray<CCString> problemTest;
            for (INT i = 0; i < allTests.Num(); ++i)
            {
                UINT uiThisError = RunTest(params, allTests[i]);
                if (0 == uiThisError)
                {
                    ++uiPassed;
                }
                else
                {
                    uiError += uiThisError;
                    problemTest.AddItem(allTests[i]->m_sParamName);
                }
            }
            timer.Stop();
            appGeneral(_T("\nRun all test with %d(success) / %d(total) (with %d errors) and %f secs\n\n\n================\n"), uiPassed, allTests.Num(), uiError, timer.Elapsed() * 0.001f);
            for (INT i = 0; i < problemTest.Num(); ++i)
            {
                appGeneral(_T("problem test suits: %s\n"), problemTest[i].c_str());
            }
            bExcuted = TRUE;
        }
        else if (sRes == _T("c"))
        {
            CTimer timer;
            timer.Start();
            UINT uiError = 0;
            UINT uiPassed = 0;
            UINT uiTotal = 0;
            TArray<CCString> problemTest;
            for (INT i = 0; i < allTests.Num(); ++i)
            {
                if (allTests[i]->IsCheck())
                {
                    ++uiTotal;
                    UINT uiThisError = RunTest(params, allTests[i]);
                    if (0 == uiThisError)
                    {
                        ++uiPassed;
                    }
                    else
                    {
                        uiError += uiThisError;
                        problemTest.AddItem(allTests[i]->m_sParamName);
                    }
                }
            }
            timer.Stop();
            appGeneral(_T("\nRun all test with %d(success) / %d(total) (with %d errors) and %f secs\n\n\n================\n"), uiPassed, uiTotal, uiError, timer.Elapsed() * 0.001f);
            for (INT i = 0; i < problemTest.Num(); ++i)
            {
                appGeneral(_T("problem test suits: %s\n"), problemTest[i].c_str());
            }
            bExcuted = TRUE;
        }
        else if (sRes == _T("class"))
        {
            GClassGather.TraceAllClass();
            appGeneral(_T("\n================================\n"));
        }
        else if (0 == sRes.Find(_T("print ")) && sRes.GetLength() > 6)
        {
            CCString sfieldid = sRes.Right(sRes.GetLength() - 6);
            BYTE byNum = appStrToBYTE(sfieldid);
            CField* field = appGetLattice()->GetFieldById(byNum);
            if (NULL != field)
            {
                field->DebugPrintMe();
            }
            else
            {
                appGeneral(_T("Print field, but field not found, the command was: %s\n"), sRes.c_str());
            }
        }
        else
        {
            TArray<CCString> keys = category.GetAllKeys();
            for (INT i = 0; i < keys.Num(); ++i)
            {
                CCString sInput = sRes;
                CCString sKey = keys[i];
                sKey.MakeLower();
                if (sInput == sKey)
                {
                    CTimer timer;
                    timer.Start();
                    UINT uiError = 0;
                    UINT uiPassed = 0;
                    TArray<CCString> unpassed;
                    for (INT j = 0; j < category[keys[i]]->Num(); ++j)
                    {
                        const TestList* pTest = category[keys[i]]->GetAt(j);
#if _CLG_DEBUG
                        if (pTest->OnlyRelease())
                            continue;
#endif

#if !_CLG_USE_LAUNCH_BOUND
                        if (pTest->OnlyBound())
                            continue;
#endif

#if !_CLG_DOUBLEFLOAT
                        if (pTest->OnlyDouble())
                            continue;
#else
                        if (pTest->OnlySingle())
                            continue;
#endif

                        UINT uiThisError = RunTest(params, pTest);
                        if (0 == uiThisError)
                        {
                            ++uiPassed;
                        }
                        else
                        {
                            uiError += uiThisError;
                            unpassed.AddItem(category[keys[i]]->GetAt(j)->m_sParamName);
                        }
                    }
                    timer.Stop();
                    appGeneral(_T("Run all %s test with %d(success) / %d(total) (with %d errors) and %f secs\n\n\n================\n"), 
                        keys[i].c_str(), uiPassed, category[keys[i]]->Num(),
                        uiError, timer.Elapsed() * 0.001f);
                    for (INT unpassidx = 0; unpassidx < unpassed.Num(); ++unpassidx)
                    {
                        appGeneral(_T("Faield:%s\n"), unpassed[unpassidx].c_str());
                    }
                    break;
                    //bExcuted = TRUE;
                }
            }
        }

        if (!bExcuted)
        {
            COUT << _T("Input commond:") << name << _T(" not kown") << std::endl;
        }
    }
    DeleteAllLists(category);
    appQuitCLG();

    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
