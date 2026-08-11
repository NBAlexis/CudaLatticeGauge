//=============================================================================
// FILENAME : DebugSU3_12.cpp
//
// DESCRIPTION:
//   Test job for deviceSU3_12 compact representation
//
// REVISION:
//  [05/10/2026 nbale]
//=============================================================================

#include "HISQ.h"

INT DebugSU3_12Job(CParameters& params)
{
    appSetupLog(params);

    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    appGeneral(_T("\n===== Running deviceSU3_12 Debug Tests =====\n\n"));

    // This launches _kernelDebugFunction which contains all SU3_12 tests
    CCudaHelper::DebugFunction();

    appGeneral(_T("\n===== deviceSU3_12 Debug Tests Done =====\n\n"));

    appQuitCLG();

    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
