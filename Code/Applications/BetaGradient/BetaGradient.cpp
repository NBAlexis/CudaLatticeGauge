//=============================================================================
// FILENAME : BetaGradient.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [08/18/2022 nbale]
//=============================================================================

#include "BetaGradient.h"


int main(int argc, char * argv[])
{
    CParameters params;
#if _CLG_DEBUG
    CYAMLParser::ParseFile(_T("StaggeredRotation.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/BetaGradient.yaml"), params);
#endif

    CCString sJob = _T("EBGJ_Simulate");

    params.FetchStringValue(_T("WorkJob"), sJob);
    EBetaGradientJob eJob = __STRING_TO_ENUM(EBetaGradientJob, sJob);

    INT res = 0;
    switch (eJob)
    {
    case EBGJ_Simulate:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobSimulate"));
            res = Simulate(workingParam1);
        }
        break;
    case EBGJ_SimulateQ:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobSimulateQ"));
            res = Simulate(workingParam1);
        }
        break;
    case EBGJ_SimulateScanBeta:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobSimulateBetaScan"));
            res = SimulateBetaScan(workingParam1);
        }
        break;
    case EBGJ_SimulateScanQ:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobSimulateBetaScanQ"));
            res = SimulateBetaScan(workingParam1);
        }
        break;
    case EBGJ_Measure:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobMeasure"));
            res = Measurement(workingParam1);
        }
        break;
    case EBGJ_MeasureScanBeta:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobMeasureBetaScan"));
            res = MeasurementBetaScan(workingParam1);
        }
        break;
    case EBGJ_GaugeFixing:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobGaugeFixing"));
            appCrucial(_T("Gauge fixing not implemented yet!\n"));
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
