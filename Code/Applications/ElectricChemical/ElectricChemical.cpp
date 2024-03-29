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

    CCString sJob = _T("EEC_Simulate");

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
            CParameters workingParam1 = params.GetParameter(_T("JobMeasure"));
            res = Measurement(workingParam1);
        }
        break;
    case EEC_GaugeFixing:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobGaugeFixing"));
            res = GaugeFixing(workingParam1);
        }
        break;
    case EEC_SimulateRW:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobSimulateRW"));
            res = SimulateRW(workingParam1);
        }
        break;
    case EEC_MeasureRW:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobMeasureImageChemical"));
            res = MeasureRW(workingParam1);
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
