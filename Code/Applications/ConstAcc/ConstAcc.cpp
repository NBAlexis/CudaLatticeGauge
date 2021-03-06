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
        CParameters workingParam1 = params.GetParameter(_T("JobSimulate"));
        res = SimulateAcc(workingParam1);
    }
    break;
    case EAJ_SimulateBoost:
    {
        CParameters workingParam1 = params.GetParameter(_T("JobSimulateBoost"));
        res = SimulateAcc(workingParam1);
    }
    break;
    case EAJ_MeasureGauge:
    {
        //CParameters workingParam2 = params.GetParameter(_T("JobThermal"));

    }
    break;
    }

    return res;
}

//=============================================================================
// END OF FILE
//=============================================================================
