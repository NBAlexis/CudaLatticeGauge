//=============================================================================
// FILENAME : CLGExample.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [04/18/2024 nbale]
//=============================================================================

#include "CLGLib.h"


int main(int argc, char * argv[])
{
    CParameters params;
#if _CLG_DEBUG && _CLG_WIN
    CYAMLParser::ParseFile(_T("CLGExample.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/CLGExample.yaml"), params);
#endif
    appSetupLog(params);
    appInitialCLG(params);

    //warm up for 5 trajectories
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
    appGetLattice()->m_pUpdator->Update(5, FALSE);

    //update for 10 trajectories and save
    appClearProfiler();
    appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, _T("test"));
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
    appGetLattice()->m_pUpdator->UpdateUntileAccept(20, FALSE);

    appDumpProfiler();

    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
