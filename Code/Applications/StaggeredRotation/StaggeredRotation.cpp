//=============================================================================
// FILENAME : StaggeredRotation.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [09/24/2020 nbale]
//=============================================================================

#include "StaggeredRotation.h"


int main(int argc, char * argv[])
{
    CParameters params;
#if _CLG_DEBUG
    CYAMLParser::ParseFile(_T("StaggeredRotation.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/StaggeredRotation.yaml"), params);
#endif

    CCString sJob = _T("ESR_Simulate");

    params.FetchStringValue(_T("RotationJob"), sJob);
    EStaggeredRotationJob eJob = __STRING_TO_ENUM(EStaggeredRotationJob, sJob);

    INT res = 0;
    switch (eJob)
    {
    case ESR_Simulate:
    {
        CParameters workingParam1 = params.GetParameter(_T("JobSimulate"));
        res = SimulateStaggeredRotation(workingParam1);
    }
    break;
    case ESR_SimulateEM:
    {
        CParameters workingParam1 = params.GetParameter(_T("JobSimulateEM"));
        res = SimulateStaggeredEM(workingParam1);
    }
    break;
    case ESR_Measure:
    {
        CParameters workingParam1 = params.GetParameter(_T("JobMeasure"));
        res = Measurement(workingParam1);
    }
    break;
    case ESR_MeasureEM:
    {
        CParameters workingParam1 = params.GetParameter(_T("JobMeasureEM"));
        res = MeasurementEM(workingParam1);
    }
    break;
    case ESR_GaugeFixing:
    {
        CParameters workingParam1 = params.GetParameter(_T("JobGaugeFixing"));
        res = GaugeFixing(workingParam1);
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
