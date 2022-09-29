//=============================================================================
// FILENAME : ElectricChemical.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [09/29/2022 nbale]
//=============================================================================

#include "ElectricChemical.h"


int main(int argc, char * argv[])
{
    CParameters params;
#if _CLG_DEBUG
    CYAMLParser::ParseFile(_T("ElectricChemical.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/ElectricChemical.yaml"), params);
#endif

    CCString sJob = _T("EBGJ_Simulate");

    params.FetchStringValue(_T("WorkJob"), sJob);
    EElectricChemical eJob = __STRING_TO_ENUM(EElectricChemical, sJob);

    INT res = 0;
    switch (eJob)
    {
    case EEC_Simulate:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobSimulate"));
            res = Simulate(workingParam1);
        }
        break;
    case EEC_Measure:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobSimulateBetaScan"));
            res = Measurement(workingParam1);
        }
        break;
    default:
        break;
    }

    return res;
}

//=============================================================================
// END OF FILE
//=============================================================================
