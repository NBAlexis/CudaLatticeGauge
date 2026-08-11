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
    CYAMLParser::ParseFile(_T("WilsonDiracElectricMagnetic.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/WilsonDiracElectricMagnetic.yaml"), params);
#endif
    CCString sJob;
    if (params.FetchStringValue(_T("Job"), sJob))
    {
        CParameters jobparam;
        if (params.FetchParameterValue(sJob, jobparam))
        {
            appSetupLog(jobparam);
            appInitialCLG(jobparam);

            INT iWarmUp = 5;
            jobparam.FetchValueINT(_T("Warmup"), iWarmUp);
            INT iConfigurationNumber = 20;
            jobparam.FetchValueINT(_T("ConfigurationNumber"), iConfigurationNumber);
            CCString sFileName = _T("Configuration");
            jobparam.FetchStringValue(_T("FileName"), sFileName);

            //warm up for 5 trajectories
            appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
            appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
            appGetLattice()->m_pUpdator->Update(iWarmUp, FALSE);

            //update for 10 trajectories and save
            appClearProfiler();
            appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, sFileName);
            appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
            appGetLattice()->m_pUpdator->UpdateUntileAccept(iConfigurationNumber, FALSE);

            appDumpProfiler();
        }
    }

    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
