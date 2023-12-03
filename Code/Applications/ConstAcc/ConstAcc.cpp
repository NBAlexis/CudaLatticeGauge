//=============================================================================
// FILENAME : ConfigurationCompresser.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [04/10/2020 nbale]
//=============================================================================

#include "ConstAcc.h"


int main(int argc, char * argv[])
{
    CParameters params;
#if _CLG_DEBUG
    CYAMLParser::ParseFile(_T("ConstAcc.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/ConstAcc.yaml"), params);
#endif

    CCString sJob = _T("EAJ_Simulate");

    params.FetchStringValue(_T("AccJob"), sJob);
    EAccJob eJob = __STRING_TO_ENUM(EAccJob, sJob);

    INT res = 0;
    switch (eJob)
    {
    case EAJ_Simulate:
    {
        CParameters workingParam = params.GetParameter(_T("JobSimulate"));
        res = SimulateAcc(workingParam);
    }
    break;
    case EAJ_SimulateQ:
    {
        CParameters workingParam = params.GetParameter(_T("JobSimulateQ"));
        res = SimulateAcc(workingParam);
    }
    break;
    case EAJ_Measure:
    {
        CParameters workingParam = params.GetParameter(_T("JobMeasure"));
        res = Measurement(workingParam);
    }
    break;
    }

    return res;
}

//=============================================================================
// END OF FILE
//=============================================================================
