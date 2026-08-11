//=============================================================================
// FILENAME : PSU3WithBoundary.cpp
//
// DESCRIPTION:
// Main entry of the PSU3WithBoundary application. Dispatches on WorkJob:
//   EPSU3J_Simulate -> Simulate (HMC production with tensor2 companion saves)
//   EPSU3J_Measure  -> Measure  (load configurations, measure observables)
//
// Usage:
//   PSU3WithBoundary [alternative.yaml]
// Without an argument the YAML is read from ../Debug/PSU3WithBoundary.yaml.
//
// REVISION:
//  [08/11/26]
//=============================================================================
#include "PSU3WithBoundary.h"

int main(int argc, char* argv[])
{
    CParameters params;
    CCString sJob = _T("EPSU3J_Simulate");

    if (2 == argc)
    {
        CCString sAlternativeYAML = argv[1];
        CYAMLParser::ParseFile(sAlternativeYAML, params);
        params.FetchStringValue(_T("WorkJob"), sJob);
    }
    else
    {
#if _CLG_DEBUG && _CLG_WIN
        CYAMLParser::ParseFile(_T("PSU3WithBoundary.yaml"), params);
#else
        CYAMLParser::ParseFile(_T("../Debug/PSU3WithBoundary.yaml"), params);
#endif
        params.FetchStringValue(_T("WorkJob"), sJob);
    }

    EPSU3Job eJob = __STRING_TO_ENUM(EPSU3Job, sJob);

    INT res = 0;
    switch (eJob)
    {
    case EPSU3J_Simulate:
        {
            CParameters workingParam = params.GetParameter(_T("JobSimulate"));
            res = Simulate(workingParam);
        }
        break;
    case EPSU3J_Measure:
        {
            CParameters workingParam = params.GetParameter(_T("JobMeasure"));
            res = Measure(workingParam);
        }
        break;
    default:
        appCrucial(_T("PSU3WithBoundary: unknown WorkJob %s\n"), sJob.c_str());
        res = 1;
        break;
    }

    return res;
}

//=============================================================================
// END OF FILE
//=============================================================================
