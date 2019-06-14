//=============================================================================
// FILENAME : RotatingReproduce.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [02/26/2019 nbale]
//=============================================================================

#include "RotatingReproduce.h"

int main(int argc, char * argv[])
{
    CParameters params;
#if _CLG_DEBUG
    CYAMLParser::ParseFile(_T("RotatingReproduce.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/RotatingReproduce.yaml"), params);
#endif

    CCString sJob = _T("ERJ_AngularMomentum");

    params.FetchStringValue(_T("RotatingJob"), sJob);
    ERotatingJob eJob = __STRING_TO_ENUM(ERotatingJob, sJob);

    INT res = 0;
    switch (eJob)
    {
    case ERJ_AngularMomentum:
    {
        CParameters workingParam1 = params.GetParameter(_T("JobAugular"));
        res = TestAngularMomentum(workingParam1);
    }
        break;
    case ERJ_Thermal:
    {
        CParameters workingParam2 = params.GetParameter(_T("JobThermal"));
        res = TestThermal(workingParam2);
    }
        break;
    }

    return res;
}

//=============================================================================
// END OF FILE
//=============================================================================
