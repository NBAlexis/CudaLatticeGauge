//=============================================================================
// FILENAME : CLGTest.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [mm/dd/yy]
//  [12/2/2018 nbale]
//=============================================================================

#include "CLGTest.h"

#if _CLG_MULTI_GPU
#include <mpi.h>
#endif

TestList* _testSuits;
CCString _bug;
CCString _msg;

// Global multi-GPU mode switch. Always FALSE except inside the --mg-worker run
// mode (multi-GPU-improve1.md 3.10): the interactive menu is strictly
// single-process/single-GPU (the 'g' toggle is removed) and the plain batch
// modes keep the identity GpuGrid=[1,1,1,1]. When TRUE, RunTest injects the
// worker --gpu-grid / --device-per-node values into the test's parameter block
// before appInitialCLG (only _TEST_MULTIGPU-tagged tests reach that path).
UBOOL g_bCLGMultiGPU = FALSE;

// Explicit grid from the worker --gpu-grid argument; g_uiDevicePerNode is the
// worker --device-per-node argument (the runtest GpuCount).
static UINT g_MultiGPUGrid[4] = { 0, 0, 0, 0 };
static UBOOL g_bGridSpecified = FALSE;
static UINT g_uiDevicePerNode = 1;

// Set by RunTest when the test was skipped (single-only under a grid, or the
// RequireSplit guard); read by the --mg-worker collective summary.
static UBOOL g_bLastTestSkipped = FALSE;

#if _CLG_MULTI_GPU
static void _clgTestFinalizeMPI()
{
    INT bFinalized = 0;
    MPI_Finalized(&bFinalized);
    if (!bFinalized)
    {
        MPI_Finalize();
    }
}
#endif

// MPI world size, initialising MPI on first use (idempotent with CLGComm::Initial).
static UINT GetMPIWorldSize()
{
#if _CLG_MULTI_GPU
    INT bInit = 0;
    MPI_Initialized(&bInit);
    if (!bInit)
    {
        MPI_Init(NULL, NULL);
        atexit(_clgTestFinalizeMPI);
    }
    INT iSize = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &iSize);
    return static_cast<UINT>(iSize);
#else
    return 1;
#endif
}

// MPI world rank, initialising MPI on first use (via GetMPIWorldSize).
static UINT GetMPIWorldRank()
{
#if _CLG_MULTI_GPU
    GetMPIWorldSize();
    INT iRank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &iRank);
    return static_cast<UINT>(iRank);
#else
    return 0;
#endif
}

// Inject the worker --gpu-grid / --device-per-node values into the test's
// parameter block before appInitialCLG (multi-GPU-improve1.md 3.10 mode 4).
// DevicePerNode is no longer forced to 1: it is exactly the runtest GpuCount
// argument, so CCLGLibManager maps rank -> device as rank % GpuCount.
static void ConfigureGridForTest(CParameters& paramForTheTest)
{
    if (!g_bGridSpecified)
    {
        //Worker argument validation guarantees a grid; without it there is
        //nothing to inject (identity [1,1,1,1] is the library default).
        return;
    }
    TArray<CCString> gridStr;
    for (UINT i = 0; i < 4; ++i)
    {
        gridStr.AddItem(appToString(g_MultiGPUGrid[i]));
    }
    paramForTheTest.SetStringVectorVaule(_T("GpuGrid"), gridStr);
    paramForTheTest.SetStringVaule(_T("DevicePerNode"), appToString(g_uiDevicePerNode));
    appGeneral(_T("[CLGTest] mg-worker: grid [%d,%d,%d,%d], device-per-node %d\n"),
        g_MultiGPUGrid[0], g_MultiGPUGrid[1], g_MultiGPUGrid[2], g_MultiGPUGrid[3], g_uiDevicePerNode);
}

// Parse a "[x, y, z, t]" or "x,y,z,t" worker --gpu-grid value into a 4-factor
// grid. Returns FALSE on a malformed argument.
static UBOOL ParseGridArgument(const CCString& sArg, UINT* pOutGrid)
{
    CCString s = sArg;
    s.TrimLeft(_T("[ "));
    s.TrimRight(_T("] "));
    TArray<CCString> parts = appGetStringList(s, _T(','), EGSLF_IgnorEmety);
    if (parts.Num() < 4)
    {
        return FALSE;
    }
    for (INT i = 0; i < 4; ++i)
    {
        CCString sPart = parts[i];
        sPart.TrimLeft(_T(" "));
        sPart.TrimRight(_T(" "));
        pOutGrid[i] = static_cast<UINT>(appStrToINT(sPart));
        if (pOutGrid[i] < 1)
        {
            return FALSE;
        }
    }
    return TRUE;
}

UINT RunTest(CParameters&params, const TestList* pTest)
{
    _bug = _T("");
    _msg = _T("");
    g_bLastTestSkipped = FALSE;

    appQuitCLG();

    //Retire any GPU work left by the previous test before the next one starts.
    //Single-card oversubscribe (multiple ranks on one device) intermittently
    //stalls a rank's kernel launch right after a test teardown; synchronising
    //here makes the boundary deterministic.
#if _CLG_MULTI_GPU
    if (g_bCLGMultiGPU)
    {
        cudaDeviceSynchronize();
        cudaGetLastError();
    }
#endif

    //Single-GPU-only test under multi-GPU mode: skip instead of running with an
    //invalid (unconfigured) grid.
    if (g_bCLGMultiGPU && GetMPIWorldSize() > 1 && !pTest->SupportMultiGPU())
    {
        appGeneral(_T("Skipped %s: single-GPU only test under multi-GPU mode.\n"), pTest->m_sParamName);
        g_bLastTestSkipped = TRUE;
        return 0;
    }

    //Test-boundary synchronisation: a previous test's collective MPI calls
    //(e.g. gauge fixing gather/fix/scatter) must be fully retired on every rank
    //before this test issues its own collectives, otherwise a rank that is late
    //leaving the previous scatter can match this test's gather and hang (seen
    //intermittently as TestMGCompressedSaveLoad hanging after
    //TestMGGaugeFixingConsistency on -n 2).
#if _CLG_MULTI_GPU
    if (g_bCLGMultiGPU && GetMPIWorldSize() > 1)
    {
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif

    appGeneral("\n=========== Testing:%s \n", pTest->m_sParamName);
    CParameters paramForTheTest = params.GetParameter(pTest->m_sParamName);
    if (!paramForTheTest.Exist(_T("MaxThreadPerBlock")))
    {
#if _CLG_DEBUG
        paramForTheTest.SetStringVaule(_T("MaxThreadPerBlock"), _T("128"));
#else
        paramForTheTest.SetStringVaule(_T("MaxThreadPerBlock"), _T("256"));
#endif
    }

    //Multi-GPU mode: configure the process grid before Initial so the test runs
    //under the requested decomposition without touching the yaml by hand.
    if (g_bCLGMultiGPU && pTest->SupportMultiGPU())
    {
        ConfigureGridForTest(paramForTheTest);
    }

    //Opt-in guard (yaml "RequireSplit : 1"): the test's lattice is only safe to
    //build when the process grid splits its oversized dimensions. The baked
    //mapping table stores coordinates in 8-bit SCHARs, so a local dimension
    //above 119 (local + 8 must fit a SCHAR) wraps and dead-loops the bake
    //kernel. Skip instead of initialising a lattice the index build cannot
    //represent (e.g. the global-T > 127 coordinate test under -n 1 or in the
    //single-GPU suite).
    INT iRequireSplit = 0;
    paramForTheTest.FetchValueINT(_T("RequireSplit"), iRequireSplit);
    if (iRequireSplit > 0)
    {
        UINT uiLatticeChk[4] = { 8, 8, 8, 8 };
        TArray<INT> latticeChk;
        if (paramForTheTest.FetchValueArrayINT(_T("LatticeLength"), latticeChk) && latticeChk.Num() >= 4)
        {
            for (UINT i = 0; i < 4; ++i)
            {
                uiLatticeChk[i] = static_cast<UINT>(latticeChk[i]);
            }
        }
        UINT uiGridChk[4] = { 1, 1, 1, 1 };
        TArray<INT> gridChk;
        if (paramForTheTest.FetchValueArrayINT(_T("GpuGrid"), gridChk) && gridChk.Num() >= 4)
        {
            for (UINT i = 0; i < 4; ++i)
            {
                uiGridChk[i] = static_cast<UINT>(gridChk[i] > 1 ? gridChk[i] : 1);
            }
        }
        UBOOL bLocalSafe = TRUE;
        for (UINT i = 0; i < 4; ++i)
        {
            const UINT uiLocal = (uiLatticeChk[i] + uiGridChk[i] - 1) / uiGridChk[i];
            if (uiLocal + 8 > 127)
            {
                bLocalSafe = FALSE;
            }
        }
        if (!bLocalSafe)
        {
            appGeneral(_T("Skipped %s: lattice needs a split grid to keep every local dimension within the 8-bit baked-coordinate bound (run with mpiexec -n >= 2 so the oversized direction decomposes).\n"), pTest->m_sParamName);
            g_bLastTestSkipped = TRUE;
            return 0;
        }
    }

    appGeneral(_T("============= Parameters %s =============\n"), paramForTheTest.GetLocation().c_str());
    paramForTheTest.Dump(_T(""));
    //Initial
    if (!appInitialCLG(paramForTheTest))
    {
        _bug = _T("Initial CLG Failed");
        return 1;
    }

    //Do the work
    appClearProfiler();
    CTimer timer;
    timer.Start();
    const UINT uiErrors = (*pTest->m_pfTest)(paramForTheTest);
    timer.Stop();
    Real fCost = timer.Elapsed();
    appGeneral(_T("=========== Finished %s, errors: %s, cost: %f(ms)\n ======== Param: %s \n ------------- End --------------\n\n"), 
        pTest->m_sParamName, 
        (0 == uiErrors) ? appDressColor(EVC_GREEN, appToString(uiErrors).c_str()).c_str() : appDressColor(EVC_RED, appToString(uiErrors).c_str()).c_str(),
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

    appDumpProfiler();

    return uiErrors;
}

void ListAllTests(const THashMap<CCString, TArray<TestList*>*>& category)
{
    COUT << _T("============== CLG v") << GetCLGVersion().c_str() << _T(": (") << appVersion() << _T(") test mode: ")
        << (g_bCLGMultiGPU ? _T("Multi-GPU") : _T("Single-GPU")) << _T(" ==============\n");
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

TestList* GetTestByName(const TArray<TestList*>& alltest, const CCString& name)
{
    for (INT i = 0; i < alltest.Num(); ++i)
    {
        CCString sShow = alltest[i]->m_sShowName;
        sShow.MakeLower();
        if (name == sShow)
        {
            return alltest[i];
        }
    }
    for (INT i = 0; i < alltest.Num(); ++i)
    {
        CCString sParam = alltest[i]->m_sParamName;
        sParam.MakeLower();
        if (name == sParam)
        {
            return alltest[i];
        }
    }
    return NULL;
}

UINT RunTestBatch(CParameters& params, const TArray<TestList*>& tests, const CCString& label)
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
    for (INT i = 0; i < tests.Num(); ++i)
    {
        const TestList* pTest = tests[i];
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
        if (g_bCLGMultiGPU && GetMPIWorldSize() > 1 && !pTest->SupportMultiGPU())
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
    appGeneral(_T("Run %s test with %d(success) / %d(total, %d skipped) (with %d errors) and %f secs\n\n\n================\n"),
        label.c_str(), uiPassed, tests.Num() - skippedname.Num(), skippedname.Num(),
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
        appGeneral(_T("Failed:%s, Bug:%s, Msg:%s\n"), unpassedname[unpassidx].c_str(), unpassedbug[unpassidx].c_str(), passedmsg[unpassidx].c_str());
    }
    return uiError;
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
    CYAMLParser::ParseFile(_T("../Debug/TestSuit.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Common.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Random.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Boundary.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_FileIO.yaml"), params);
    //P5-3.1: the multi-GPU tests live in their own group file. Parsed after
    //TestSuit_FileIO.yaml so params do not collide (segment names are unique).
    //Improve-1: the I1 known-failure group file was dissolved at I4 -- every
    //block moved here as its fix landed (multi-GPU-improve1.md I4).
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_MG.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_FermionMatrix.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_FermionUpdator.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_FermionUpdatorKS.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_GaugeFixing.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Rotation.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_MiscUpdate.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Solver.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Updator.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_HmcDiagnostics.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_Boson.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_SUN.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_NumpyExport.yaml"), params);
    CYAMLParser::ParseFile(_T("../Debug/TestSuit_PSU3WithBoundary.yaml"), params);
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

//3.10 mode 3 (--mg-config TestName): print the test's multi-rank driver
//metadata as exactly one ASCII line "rankCount gx gy gz gt" on stdout. No CUDA
//initialisation, no MPI collective; diagnostics go to stderr (the caller has
//already routed the tracer there). Unknown tests, tests not registered with
//_TEST_MULTIGPU, and missing/invalid metadata all return non-zero and never
//print a fake config.
static INT RunMGConfigQuery(const TArray<TestList*>& allTests, CParameters& params, const CCString& sTestName)
{
    CCString sName(sTestName);
    sName.MakeLower();
    const TestList* pTest = GetTestByName(allTests, sName);
    if (NULL == pTest)
    {
        fprintf(stderr, "[CLGTest] --mg-config: unknown test: %s\n", sTestName.c_str());
        return 1;
    }
    if (!pTest->IsMultiGPUTagged())
    {
        fprintf(stderr, "[CLGTest] --mg-config: %s is not registered with _TEST_MULTIGPU\n", pTest->m_sParamName);
        return 1;
    }

    CParameters paramForTheTest = params.GetParameter(pTest->m_sParamName);
    INT iRankCount = 0;
    if (!paramForTheTest.FetchValueINT(_T("MultiGPUTestRankCount"), iRankCount) || iRankCount < 1)
    {
        fprintf(stderr, "[CLGTest] --mg-config: %s has no valid MultiGPUTestRankCount\n", pTest->m_sParamName);
        return 1;
    }
    TArray<INT> grid;
    if (!paramForTheTest.FetchValueArrayINT(_T("MultiGPUTestGrid"), grid) || grid.Num() < 4)
    {
        fprintf(stderr, "[CLGTest] --mg-config: %s has no valid MultiGPUTestGrid\n", pTest->m_sParamName);
        return 1;
    }
    INT iProduct = 1;
    for (INT i = 0; i < 4; ++i)
    {
        if (grid[i] < 1)
        {
            fprintf(stderr, "[CLGTest] --mg-config: %s MultiGPUTestGrid factors must be >= 1\n", pTest->m_sParamName);
            return 1;
        }
        iProduct *= grid[i];
    }
    if (iProduct != iRankCount)
    {
        fprintf(stderr, "[CLGTest] --mg-config: %s MultiGPUTestGrid product %d != MultiGPUTestRankCount %d\n",
            pTest->m_sParamName, iProduct, iRankCount);
        return 1;
    }

    printf("%d %d %d %d %d\n", iRankCount, grid[0], grid[1], grid[2], grid[3]);
    return 0;
}

//3.10 mode 4 (TestName --mg-worker --gpu-grid x,y,z,t --device-per-node G):
//validate the launch, then run exactly one _TEST_MULTIGPU-tagged test under
//the given grid. Every validation failure fails collectively -- all ranks
//return the same non-zero code before any test kernel runs.
static INT RunMGWorker(CParameters& params, const TestList* pTest, const UINT* pGrid, UINT uiDevicePerNode)
{
#if _CLG_MULTI_GPU
    INT iFail = 0;
    const UINT uiWorld = GetMPIWorldSize();
    const UINT uiRank = GetMPIWorldRank();

    const UINT uiProduct = pGrid[0] * pGrid[1] * pGrid[2] * pGrid[3];
    if (uiProduct != uiWorld)
    {
        fprintf(stderr, "[CLGTest] rank %d: --mg-worker: MPI world size %d != grid [%d,%d,%d,%d] product %d\n",
            uiRank, uiWorld, pGrid[0], pGrid[1], pGrid[2], pGrid[3], uiProduct);
        iFail = 2;
    }

    //Single-node only this round: refuse a multi-node world outright instead
    //of silently applying a global "rank % G" mapping across nodes.
    MPI_Comm nodeComm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &nodeComm);
    INT iNodeSize = 1;
    MPI_Comm_size(nodeComm, &iNodeSize);
    MPI_Comm_free(&nodeComm);
    if (static_cast<UINT>(iNodeSize) != uiWorld)
    {
        fprintf(stderr, "[CLGTest] rank %d: --mg-worker: multi-node worlds are not supported by this driver (node ranks %d of world %d)\n",
            uiRank, iNodeSize, uiWorld);
        iFail = 2;
    }

    //The visible device count is verified here with the CUDA API (the script
    //side must not parse nvidia-smi); a shortage fails before any test kernel.
    INT iDeviceCount = 0;
    if (cudaSuccess != cudaGetDeviceCount(&iDeviceCount) || iDeviceCount < 0)
    {
        iDeviceCount = 0;
    }
    if (static_cast<UINT>(iDeviceCount) < uiDevicePerNode)
    {
        fprintf(stderr, "[CLGTest] rank %d: --mg-worker: --device-per-node %d but only %d visible device(s) on this node\n",
            uiRank, uiDevicePerNode, iDeviceCount);
        iFail = 3;
    }

    INT iFailAll = 0;
    MPI_Allreduce(&iFail, &iFailAll, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (0 != iFailAll)
    {
        return iFailAll;
    }

    for (UINT i = 0; i < 4; ++i)
    {
        g_MultiGPUGrid[i] = pGrid[i];
    }
    g_bGridSpecified = TRUE;
    g_uiDevicePerNode = uiDevicePerNode;
    g_bCLGMultiGPU = TRUE;

    //rank -> device mapping is rank % GpuCount (CCLGLibManager applies it as
    //rank % DevicePerNode); rank 0 reports the full mapping once so the
    //diagnostic log is verifiable without N duplicated summaries.
    if (0 == uiRank)
    {
        for (UINT r = 0; r < uiWorld; ++r)
        {
            appGeneral(_T("[CLGTest] mg-worker mapping: rank %d -> device %d\n"), r, r % uiDevicePerNode);
        }
        if (uiWorld > uiDevicePerNode)
        {
            appGeneral(_T("[CLGTest] mg-worker: single-GPU oversubscription: %d ranks on %d device(s); MPI decomposition, halo, collectives and indexing are validated, real cross-device behaviour is NOT.\n"),
                uiWorld, uiDevicePerNode);
        }
    }

    const UINT uiErrors = RunTest(params, pTest);

    //Collective end-of-test summary: error count, fatal (printed CRUCIAL)
    //count and skip status are reduced over all ranks; only rank 0 prints the
    //summary, and every rank returns the same exit code. The exit code follows
    //the batch-mode rule: errors only. The fatal count is diagnostic —
    //expected-fatal tests (launch guards, handle-registry misuse) deliberately
    //trigger CRUCIALs and still finish with errors=0.
    INT iLocal[3];
    iLocal[0] = static_cast<INT>(uiErrors);
    iLocal[1] = static_cast<INT>(GTracer.GetFatalCount());
    iLocal[2] = g_bLastTestSkipped ? 1 : 0;
    INT iSum[3] = { 0, 0, 0 };
    MPI_Allreduce(iLocal, iSum, 3, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (0 == uiRank)
    {
        appGeneral(_T("[CLGTest] mg-worker summary: test=%s ranks=%d grid=[%d,%d,%d,%d] devicePerNode=%d errors=%d fatals=%d skippedRanks=%d oversubscription=%s\n"),
            pTest->m_sParamName, uiWorld, pGrid[0], pGrid[1], pGrid[2], pGrid[3], uiDevicePerNode,
            iSum[0], iSum[1], iSum[2], (uiWorld > uiDevicePerNode) ? _T("yes") : _T("no"));
    }
    return (iSum[0] > 0) ? 1 : 0;
#else
    fprintf(stderr, "[CLGTest] --mg-worker requires a multi-GPU (-DCLG_MULTI_GPU=1) build.\n");
    return 2;
#endif
}

int main(int argc, char * argv[])
{
    //multi-GPU-improve1.md 3.10 mode 3: "--mg-config TestName" prints exactly
    //one ASCII line ("rankCount gx gy gz gt") on stdout, so every tracer/log
    //line is routed to stderr for the whole process in this mode.
    const UBOOL bMGConfigMode = (argc > 1 && 0 == strcmp(argv[1], "--mg-config"));
    if (bMGConfigMode)
    {
        GTracer.SetLogToStdErr();
    }

    //Load settings
    CParameters params;
    LoadParams(params);
    if (!bMGConfigMode)
    {
        appSetupLog(params);
    }

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

    BYTE* testCudaBuffer = NULL;
    CCudaBuffer testBuffer;
    
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

    //3.10 mode 3: metadata query. The tracer was routed to stderr at startup,
    //so stdout carries only the single config line this prints.
    if (bMGConfigMode)
    {
        INT iRet = 1;
        if (argc > 2)
        {
            iRet = RunMGConfigQuery(allTests, params, CCString(argv[2]));
        }
        else
        {
            fprintf(stderr, "[CLGTest] usage: CLGTest --mg-config TestName\n");
        }
        DeleteAllLists(category);
        return iRet;
    }

    //3.10: the interactive menu and the plain batch modes are strictly
    //single-process. Under an MPI world with more than one rank only the
    //--mg-worker form is legal; refuse before any rank reaches a stdin read
    //or runs a test with an unconfigured (identity) grid.
    UBOOL bWorkerRequested = FALSE;
    for (INT i = 1; i < argc; ++i)
    {
        CCString sToken(argv[i]);
        if (0 == sToken.CompareNoCase(_T("--mg-worker")))
        {
            bWorkerRequested = TRUE;
        }
    }
    if (!bWorkerRequested && GetMPIWorldSize() > 1)
    {
        if (0 == GetMPIWorldRank())
        {
            fprintf(stderr, "[CLGTest] interactive and plain batch modes are single-process only; run multi-rank tests via runtest.sh TestName GpuCount (mpiexec ... --mg-worker).\n");
        }
        DeleteAllLists(category);
        return 2;
    }

    // CLI modes (multi-GPU-improve1.md 3.10):
    //   CLGTest <name|category|all>  -- single-process batch, identity grid
    //   CLGTest TestName --mg-worker --gpu-grid gx,gy,gz,gt --device-per-node G
    if (argc > 1)
    {
        CCString sArg(argv[1]);
        UBOOL bHasGridArg = FALSE;
        CCString sGridValue;
        UINT uiWorkerDevicePerNode = 0;
        for (INT i = 2; i < argc; ++i)
        {
            CCString sToken(argv[i]);
            if (0 == sToken.CompareNoCase(_T("--mg-worker")))
            {
                continue;
            }
            if (0 == sToken.CompareNoCase(_T("--gpu-grid")) && i + 1 < argc)
            {
                //Collect "[2, 2, 1, 1]" even when the shell split it into
                //separate argv tokens (no quotes around the brackets).
                sGridValue = CCString(argv[++i]);
                while ((sGridValue.Find(_T("[")) >= 0) && (sGridValue.Find(_T("]")) < 0) && (i + 1 < argc))
                {
                    sGridValue += _T(" ") + CCString(argv[++i]);
                }
                bHasGridArg = TRUE;
                continue;
            }
            if (0 == sToken.CompareNoCase(_T("--device-per-node")) && i + 1 < argc)
            {
                uiWorkerDevicePerNode = static_cast<UINT>(appStrToINT(CCString(argv[++i])));
                continue;
            }
            appCrucial(_T("Unknown argument: %s\n"), argv[i]);
            DeleteAllLists(category);
            appQuitCLG();
            return 2;
        }

        if (bWorkerRequested)
        {
            CCString sTestName(sArg);
            sTestName.MakeLower();
            TestList* pTest = GetTestByName(allTests, sTestName);
            if (NULL == pTest || !pTest->IsMultiGPUTagged())
            {
                if (0 == GetMPIWorldRank())
                {
                    fprintf(stderr, "[CLGTest] --mg-worker: %s must name a test registered with _TEST_MULTIGPU (not a category, not all)\n", sArg.c_str());
                }
                DeleteAllLists(category);
                return 2;
            }
            UINT uiWorkerGrid[4] = { 1, 1, 1, 1 };
            if (!bHasGridArg || !ParseGridArgument(sGridValue, uiWorkerGrid) || uiWorkerDevicePerNode < 1)
            {
                if (0 == GetMPIWorldRank())
                {
                    fprintf(stderr, "[CLGTest] usage: CLGTest TestName --mg-worker --gpu-grid gx,gy,gz,gt --device-per-node G (G >= 1)\n");
                }
                DeleteAllLists(category);
                return 2;
            }
            INT iRet = RunMGWorker(params, pTest, uiWorkerGrid, uiWorkerDevicePerNode);
            DeleteAllLists(category);
            appQuitCLG();
            return iRet;
        }

        sArg.MakeLower();

        if (sArg == _T("all"))
        {
            UINT errors = RunTestBatch(params, allTests, _T("all"));
            DeleteAllLists(category);
            appQuitCLG();
            return (errors > 0) ? 1 : 0;
        }

        // Check category
        TArray<CCString> keys = category.GetAllKeys();
        for (INT i = 0; i < keys.Num(); ++i)
        {
            CCString sKey = keys[i];
            sKey.MakeLower();
            if (sArg == sKey)
            {
                UINT errors = RunTestBatch(params, *category[keys[i]], keys[i]);
                DeleteAllLists(category);
                appQuitCLG();
                return (errors > 0) ? 1 : 0;
            }
        }

        // Check single test by name
        TestList* pTest = GetTestByName(allTests, sArg);
        if (NULL != pTest)
        {
            UINT errors = RunTest(params, pTest);
            DeleteAllLists(category);
            appQuitCLG();
            return (errors > 0) ? 1 : 0;
        }

        appCrucial(_T("Unknown test or category: %s\n"), argv[1]);
        DeleteAllLists(category);
        appQuitCLG();
        return 1;
    }

    ListAllTests(category);
    while (TRUE)
    {
        COUT << _T("============== CLG v") << GetCLGVersion().c_str() << _T(": (") << appVersion() << _T(") test mode: ")
            << _T("Single-GPU") << _T(" ==============\nq - Quit,  l - List all,  r - Run all,  p - reload params, print - print info\ntest memory: allocate n (in MB)/free\n");
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
                RunTest(params, pTest);
                bExcuted = TRUE;
            }
        }
        else if (sRes == _T("r"))
        {
            RunTestBatch(params, allTests, _T("all"));
            bExcuted = TRUE;
        }
        else if (0 == sRes.Find(_T("print")))
        {
            Print(sRes);
            bExcuted = TRUE;
        }
        else if (0 == sRes.Find(_T("allocate")))
        {
            TArray <CCString> sArgs = appGetStringList(sRes, _T(' '), EGSLF_IgnorTabSpaceInSide | EGSLF_IgnorEmety);
            if (2 == sArgs.Num())
            {
                ULONGLONG numberInMB = static_cast<ULONGLONG>(appStrToINT(sArgs[1]));
                if (numberInMB > 0 && numberInMB < 100000)
                {
                    if (NULL != testCudaBuffer)
                    {
                        checkCudaErrors(testBuffer.CudaFree(testCudaBuffer));
                        testCudaBuffer = NULL;
                    }
                    checkCudaErrors(testBuffer.CudaMalloc((void**)&testCudaBuffer, numberInMB * (1 << 20), _T(__FILE__), __LINE__));
                }
            }
            bExcuted = TRUE;
        }
        else if (0 == sRes.Find(_T("free")))
        {
            if (NULL != testCudaBuffer)
            {
                checkCudaErrors(testBuffer.CudaFree(testCudaBuffer));
                testCudaBuffer = NULL;
            }
            bExcuted = TRUE;
        }
        else
        {
            TArray<CCString> keys = category.GetAllKeys();
            for (INT i = 0; i < keys.Num(); ++i)
            {
                CCString sKey = keys[i];
                sKey.MakeLower();
                if (sRes == sKey)
                {
                    RunTestBatch(params, *category[keys[i]], keys[i]);
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
