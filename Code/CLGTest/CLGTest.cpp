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
CCString _bug;
CCString _msg;

UINT RunTest(CParameters&params, const TestList* pTest)
{
    _bug = _T("");
    _msg = _T("");

    appQuitCLG();

    appGeneral("\n=========== Testing:%s \n", pTest->m_sParamName);
    CParameters paramForTheTest = params.GetParameter(pTest->m_sParamName);
    appGeneral(_T("============= Parameters %s =============\n"), paramForTheTest.GetLocation().c_str());
    paramForTheTest.Dump(_T(""));
    //Initial
    if (!appInitialCLG(paramForTheTest))
    {
        _bug = _T("Initial CLG Failed");
        return 1;
    }

    //Do the work
    CTimer timer;
    timer.Start();
    const UINT uiErrors = (*pTest->m_pfTest)(paramForTheTest);
    timer.Stop();
    Real fCost = timer.Elapsed();
    appGeneral(_T("=========== Finished %s, errors: %d, cost: %f(ms)\n ======== Param: %s \n ------------- End --------------\n\n"), 
        pTest->m_sParamName, 
        uiErrors, 
        fCost,
        paramForTheTest.GetLocation().c_str());

    CCString sFinishMsg;
    sFinishMsg.Format(_T(" cost:%f ms"), fCost);
    AddMsg(sFinishMsg);

#if _CLG_WIN
    OutputDebugString(_T("Param Name: "));
    OutputDebugString(paramForTheTest.GetName().c_str());
    OutputDebugString(_T(", Error: "));
    OutputDebugString(appToString(uiErrors).c_str());
    OutputDebugString(_T(", Cost: "));
    OutputDebugString(appToString(timer.Elapsed()).c_str());
    OutputDebugString(_T(", Last probem: "));
    OutputDebugString(_bug.c_str());
    OutputDebugString(_T(", double click following:\n"));
    OutputDebugString(paramForTheTest.GetLocation().c_str());
    OutputDebugString(_T("\n"));
#endif
    //std::cerr << paramForTheTest.GetLocation().c_str() << std::endl;
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

void LoadParams(CParameters& params)
{
    params.RemoveAll();
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
    CYAMLParser::ParseFile(_T("TestSuit_SUN.yaml"), params);
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
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_SUN.yaml"), params);
    //CYAMLParser::ParseFile(_T("../Debug/TestSuit_EvenOdd.yaml"), params);
#endif
}

void Print(const CCString& command)
{
    TArray <CCString> sArgs = appGetStringList(command, _T(' '), EGSLF_IgnorTabSpaceInSide | EGSLF_IgnorEmety);
    if (sArgs.Num() > 1)
    {
        if (sArgs[1] == _T("class"))
        {
            GClassGather.TraceAllClass();
            appGeneral(_T("\n================================\n"));
            return;
        }
        else if (sArgs[1] == _T("device"))
        {
            CCudaHelper::DeviceQuery();
        }
        else if (sArgs[1] == _T("lattice"))
        {
            appGeneral(_T("\n=======================================\n") + appGetLattice()->GetInfos(_T("")) + _T("\n=======================================\n"));
        }
        else if (sArgs[1] == _T("field"))
        {
            if (sArgs.Num() > 2)
            {
                BYTE byNum = appStrToBYTE(sArgs[2]);
                CField* field = appGetLattice()->GetFieldById(byNum);
                if (NULL != field)
                {
                    field->DebugPrintMe();
                }
                else
                {
                    appGeneral(_T("Print field, but field not found, the command was: %s\n"), command.c_str());
                }
                return;
            }
        }
        else if (sArgs[1] == _T("bond"))
        {
            if (sArgs.Num() > 2)
            {
                BYTE byNum = appStrToBYTE(sArgs[2]);
                CIndexData::DebugLinkDirichletOrDagger(byNum);
                return;
            }
        }
    }

    appGeneral(_T("not recognizd cmd: %s\n supported: class, device, lattice, field [d], position [d], bond [d]\n"), command.c_str());
}

int main(int argc, char * argv[])
{
    //Load settings
    CParameters params;
    LoadParams(params);
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
        COUT << _T("============== CLG v") << GetCLGVersion().c_str() << _T(": (") << appVersion() << _T(") ==============\nq - Quit,  l - List all,  r - Run all,  p - reload params, print - print info\n");
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

        if (sRes == _T("p"))
        {
            LoadParams(params);
            ListAllTests(category);
            bExcuted = TRUE;
        }

        if (sRes == _T("l"))
        {
            ListAllTests(category);
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
            TArray<CCString> skippedname;
            TArray<CCString> unpassedname;
            TArray<CCString> unpassedbug;
            TArray<CCString> unpassedmsg;
            TArray<CCString> passedname;
            TArray<CCString> passedmsg;
            for (INT i = 0; i < allTests.Num(); ++i)
            {
                const TestList* pTest = allTests[i];
#if _CLG_DEBUG
                if (pTest->OnlyRelease())
                {
                    skippedname.AddItem(pTest->m_sParamName);
                    continue;
                }
#endif

#if !_CLG_USE_LAUNCH_BOUND
                if (pTest->OnlyBound())
                {
                    skippedname.AddItem(pTest->m_sParamName);
                    continue;
                }
#endif

#if !_CLG_DOUBLEFLOAT
                if (pTest->OnlyDouble())
                {
                    skippedname.AddItem(pTest->m_sParamName);
                    continue;
                }
#else
                if (pTest->OnlySingle())
                {
                    skippedname.AddItem(pTest->m_sParamName);
                    continue;
                }
#endif
                if (pTest->NoCheck())
                {
                    skippedname.AddItem(pTest->m_sParamName);
                    continue;
                }

                UINT uiThisError = RunTest(params, pTest);
                if (0 == uiThisError)
                {
                    ++uiPassed;
                    passedname.AddItem(pTest->m_sParamName);
                    unpassedmsg.AddItem(_msg);
                }
                else
                {
                    uiError += uiThisError;
                    unpassedname.AddItem(pTest->m_sParamName);
                    unpassedbug.AddItem(_bug);
                    passedmsg.AddItem(_msg);
                }
            }
            timer.Stop();
            appGeneral(_T("Run all test with %d(success) / %d(total, %d skipped) (with %d errors) and %f secs\n\n\n================\n"),
                uiPassed, allTests.Num() - skippedname.Num(), skippedname.Num(),
                uiError, timer.Elapsed() * 0.001f);

            for (INT unpassidx = 0; unpassidx < skippedname.Num(); ++unpassidx)
            {
                appGeneral(_T("Skipped:%s\n"), skippedname[unpassidx].c_str());
            }

            for (INT unpassidx = 0; unpassidx < passedname.Num(); ++unpassidx)
            {
                appGeneral(_T("Success:%s, %s\n"), passedname[unpassidx].c_str(), unpassedmsg[unpassidx].c_str());
            }

            for (INT unpassidx = 0; unpassidx < unpassedname.Num(); ++unpassidx)
            {
                appGeneral(_T("Faield:%s, Bug:%s, Msg:%s\n"), unpassedname[unpassidx].c_str(), unpassedbug[unpassidx].c_str(), passedmsg[unpassidx].c_str());
            }
            bExcuted = TRUE;
        }
        else if (0 == sRes.Find(_T("print")))
        {
            Print(sRes);
            bExcuted = TRUE;
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
                    TArray<CCString> skippedname;
                    TArray<CCString> unpassedname;
                    TArray<CCString> unpassedbug;
                    TArray<CCString> unpassedmsg;
                    TArray<CCString> passedname;
                    TArray<CCString> passedmsg;
                    for (INT j = 0; j < category[keys[i]]->Num(); ++j)
                    {
                        const TestList* pTest = category[keys[i]]->GetAt(j);
#if _CLG_DEBUG
                        if (pTest->OnlyRelease())
                        {
                            skippedname.AddItem(pTest->m_sParamName);
                            continue;
                        }
#endif

#if !_CLG_USE_LAUNCH_BOUND
                        if (pTest->OnlyBound())
                        {
                            skippedname.AddItem(pTest->m_sParamName);
                            continue;
                        }
#endif

#if !_CLG_DOUBLEFLOAT
                        if (pTest->OnlyDouble())
                        {
                            skippedname.AddItem(pTest->m_sParamName);
                            continue;
                        }
#else
                        if (pTest->OnlySingle())
                        {
                            skippedname.AddItem(pTest->m_sParamName);
                            continue;
                        }
#endif
                        if (pTest->NoCheck())
                        {
                            skippedname.AddItem(pTest->m_sParamName);
                            continue;
                        }

                        UINT uiThisError = RunTest(params, pTest);
                        if (0 == uiThisError)
                        {
                            ++uiPassed;
                            passedname.AddItem(pTest->m_sParamName);
                            unpassedmsg.AddItem(_msg);
                        }
                        else
                        {
                            uiError += uiThisError;
                            unpassedname.AddItem(pTest->m_sParamName);
                            unpassedbug.AddItem(_bug);
                            passedmsg.AddItem(_msg);
                        }
                    }
                    timer.Stop();
                    appGeneral(_T("Run all %s test with %d(success) / %d(total %d skipped) (with %d errors) and %f secs\n\n\n================\n"), 
                        keys[i].c_str(), uiPassed, category[keys[i]]->Num() - skippedname.Num(), skippedname.Num(),
                        uiError, timer.Elapsed() * 0.001f);

                    for (INT unpassidx = 0; unpassidx < skippedname.Num(); ++unpassidx)
                    {
                        appGeneral(_T("Skipped:%s\n"), skippedname[unpassidx].c_str());
                    }

                    for (INT unpassidx = 0; unpassidx < passedname.Num(); ++unpassidx)
                    {
                        appGeneral(_T("Success:%s, %s\n"), passedname[unpassidx].c_str(), unpassedmsg[unpassidx].c_str());
                    }

                    for (INT unpassidx = 0; unpassidx < unpassedname.Num(); ++unpassidx)
                    {
                        appGeneral(_T("Faield:%s, Bug:%s, Msg:%s\n"), unpassedname[unpassidx].c_str(), unpassedbug[unpassidx].c_str(), passedmsg[unpassidx].c_str());
                    }
                    break;
                }
            }
            bExcuted = TRUE;
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
