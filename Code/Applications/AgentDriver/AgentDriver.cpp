//=============================================================================
// FILENAME : AgentDriver.cpp
// 
// DESCRIPTION:
// Unified YAML-driven entry point for agent use.
// Reads a YAML file and dispatches to Simulate / GaugeFixing / Measure.
//
// REVISION:
//  [06/14/2026]
//=============================================================================

#include "AgentDriver.h"
#include "AgentDriver.h"

__BEGIN_NAMESPACE

enum EAgentDriverTask
{
    EADT_Simulate = 0,
    EADT_GaugeFixing,
    EADT_Measure,

    EADT_Count,
};

static EAgentDriverTask __TaskNameToEnum(const CCString& sName)
{
    if (sName == _T("Simulate")) return EADT_Simulate;
    if (sName == _T("GaugeFixing")) return EADT_GaugeFixing;
    if (sName == _T("Measure")) return EADT_Measure;
    return EADT_Count;
}

__END_NAMESPACE

__USE_NAMESPACE

INT main(INT argc, TCHAR* argv[])
{
    CParameters params;
    if (argc < 2)
    {
        appGeneral(_T("Usage: AgentDriver <yaml-file>\n"));
        return 1;
    }

    CCString sYamlFile = argv[1];

    // Default to GENERAL verbosity so the driver is always observable.
    // appSetupLog() may override this from the YAML file.
    appSetTracer(GENERAL, _T("stdout"));

    appGeneral(_T("\n==================================================\n"));
    appGeneral(_T("AgentDriver started\n"));
    appGeneral(_T("YAML file: %s\n"), sYamlFile.c_str());
    appGeneral(_T("==================================================\n"));

    CYAMLParser::ParseFile(sYamlFile, params);
    appGeneral(_T("YAML parsed successfully.\n"));

    CCString sTask = _T("Simulate");
    params.FetchStringValue(_T("Task"), sTask);
    const EAgentDriverTask eTask = __TaskNameToEnum(sTask);
    appGeneral(_T("Task: %s\n"), sTask.c_str());

    INT iRet = 0;
    switch (eTask)
    {
    case EADT_Simulate:
        iRet = RunAgentSimulate(params);
        break;
    case EADT_GaugeFixing:
        iRet = RunAgentGaugeFixing(params);
        break;
    case EADT_Measure:
        iRet = RunAgentMeasure(params);
        break;
    default:
        appCrucial(_T("Unknown Task: %s, expected Simulate/GaugeFixing/Measure\n"), sTask.c_str());
        iRet = 1;
        break;
    }

    appGeneral(_T("\n==================================================\n"));
    appGeneral(_T("AgentDriver finished with exit code %d\n"), iRet);
    appGeneral(_T("==================================================\n\n"));

    return iRet;
}

//=============================================================================
// END OF FILE
//=============================================================================
