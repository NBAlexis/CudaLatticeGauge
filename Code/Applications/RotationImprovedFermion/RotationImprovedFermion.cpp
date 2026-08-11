//=============================================================================
// FILENAME : RotationImprovedFermion.cpp
// 
// DESCRIPTION:
//
// REVISION:
// [mm/dd/yy]
// [08/22/2025 nbale]
//=============================================================================

#include "RotationImprovedFermion.h"


int main(int argc, char * argv[])
{
    CParameters params;
    CCString sJob = _T("ERIF_SimulateHISQ");

    if (2 == argc)
    {
        CCString sAlternativeYAML = argv[1];
        CYAMLParser::ParseFile(sAlternativeYAML, params);
        params.FetchStringValue(_T("WorkJob"), sJob);
    }
    else
    {
#if _CLG_DEBUG && _CLG_WIN
        CYAMLParser::ParseFile(_T("RotationImprovedFermion.yaml"), params);
#else
        CYAMLParser::ParseFile(_T("../Debug/RotationImprovedFermion.yaml"), params);
#endif

        params.FetchStringValue(_T("WorkJob"), sJob);
    }
    ERotationImproveFermionJob eJob = __STRING_TO_ENUM(ERotationImproveFermionJob, sJob);

    INT res = 0;
    switch (eJob)
    {
    case ERIF_SimulateHISQ:
        {
            CParameters workingParam1 = params.GetParameter(_T("HISQ"));
            res = Simulate(workingParam1);
        }
        break;
    case ERIF_SimulateHISQRotation:
        {
            CParameters workingParam2 = params.GetParameter(_T("HISQRotation"));
            res = Simulate(workingParam2);
        }
        break;
    case ERIF_SimulateStoutlink:
        {
            CParameters workingParam3 = params.GetParameter(_T("Stoutlink"));
            res = Simulate(workingParam3);
        }
        break;
    case ERIF_SimulateStoutlinkRotation:
        {
            CParameters workingParam4 = params.GetParameter(_T("StoulinkRotation"));
            res = SimulateRotation(workingParam4);
        }
        break;
    case ERIF_MeasureHISQ:
        {
            CParameters workingParam5 = params.GetParameter(_T("MeasureHISQ"));
            res = MeasureHISQ(workingParam5);
        }
        break;
    case ERIF_MeasureHISQRotation:
        {
            CParameters workingParam6 = params.GetParameter(_T("MeasureHISQRotation"));
            res = MeasureHISQRotation(workingParam6);
        }
        break;
    }
    return res;
}

//=============================================================================
// END OF FILE
//=============================================================================
