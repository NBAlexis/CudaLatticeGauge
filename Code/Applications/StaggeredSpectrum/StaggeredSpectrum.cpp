//=============================================================================
// FILENAME : StaggeredSpectrum.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [08/21/2020 nbale]
//=============================================================================

#include "StaggeredSpectrum.h"


int main(int argc, char * argv[])
{
    CParameters params;
#if _CLG_DEBUG
    CYAMLParser::ParseFile(_T("StaggeredSpectrum.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/StaggeredSpectrum.yaml"), params);
#endif

    CCString sJob = _T("ESSJ_Simulate");

    params.FetchStringValue(_T("StaggeredSpectrumJob"), sJob);
    const EStaggeredSpectrumJob eJob = __STRING_TO_ENUM(EStaggeredSpectrumJob, sJob);

    INT res = 0;
    switch (eJob)
    {
    case ESSJ_Simulate:
    {
        CParameters workingParam1 = params.GetParameter(_T("JobSimulate"));
        res = SimulateStaggered(workingParam1);
    }
    break;
    case ESSJ_Measure:
    {
        CParameters workingParam2 = params.GetParameter(_T("JobMeasure"));
        res = Measurement(workingParam2);
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
