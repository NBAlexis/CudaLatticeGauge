//=============================================================================
// FILENAME : Simulate.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [12/08/2024 nbale]
//=============================================================================

#include "HISQ.h"

INT SimulateTest(CParameters& params)
{
#pragma region Parameters

    appSetupLog(params);

    //Number of trajectories, calculate h_diff, but accept all
    INT iValue = 0;
    params.FetchValueINT(_T("Warmup"), iValue);
    UINT uiBeforeMetropolis = static_cast<UINT>(iValue);

    //Number of trajectories, calculate h_diff, but accept all
    iValue = 1;
    params.FetchValueINT(_T("Metropolis"), iValue);
    UINT uiMetropolis = static_cast<UINT>(iValue);

    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    //Measure
    appGetLattice()->m_pMeasurements->Reset();

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
    appGetLattice()->m_pUpdator->Update(uiBeforeMetropolis, FALSE);

    appClearProfiler();
    appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, _T("HISQTEST"));
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
    appGetLattice()->m_pUpdator->Update(uiMetropolis, FALSE);

    appGeneral(_T("\n========= Finished! ==========\n\n"));

    appDumpProfiler();

    appQuitCLG();

    return 0;
}
