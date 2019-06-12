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
        res = TestAngularMomentum(params);
        break;
    case ERJ_Thermal:
        res = TestThermal(params);
        break;
    }

    return res;
}

//=============================================================================
// END OF FILE
//=============================================================================
