//=============================================================================
// FILENAME : HISQ.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [12/08/2024 nbale]
//=============================================================================

#include "HISQ.h"


int main(int argc, char * argv[])
{
    CParameters params;
    if (argc > 1)
    {
        CCString sYaml(argv[1]);
        CYAMLParser::ParseFile(sYaml, params);
    }
    else
    {
#if _CLG_DEBUG && _CLG_WIN
        CYAMLParser::ParseFile(_T("HISQ.yaml"), params);
#else
        CYAMLParser::ParseFile(_T("../Debug/HISQ.yaml"), params);
#endif
    }

    CCString sJob = _T("EHJ_SimulateNf2p1");

    if (argc > 2)
    {
        sJob = argv[2];
    }
    else
    {
        params.FetchStringValue(_T("WorkJob"), sJob);
    }
    EHISQJob eJob = __STRING_TO_ENUM(EHISQJob, sJob);

    INT res = 0;
    switch (eJob)
    {
    case EHJ_SimulateTest:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobSimulateTest"));
            res = SimulateTest(workingParam1);
        }
        break;
    case EHJ_SimulateNf2p1:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobSimulateNf2p1"));
            res = SimulateNf2p1(workingParam1);
        }
        break;
    case EHJ_Measure:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobMeasure"));
            res = Measurement(workingParam1);
        }
        break;
    case EHJ_GaugeFixing:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobGaugeFixing"));
            res = GaugeFixing(workingParam1);
        }
        break;
    case EHJ_Task1Full:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobTask1Full"));
            if (argc > 3) workingParam1.SetStringVaule(_T("CfgFile"), argv[3]);
            if (argc > 4) workingParam1.SetStringVaule(_T("VFile"), argv[4]);
            if (argc > 5) workingParam1.SetStringVaule(_T("OutputFile"), argv[5]);
            res = Task1FullJob(workingParam1);
        }
        break;
    case EHJ_Task1EO:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobTask1EO"));
            if (argc > 3) workingParam1.SetStringVaule(_T("CfgFile"), argv[3]);
            if (argc > 4) workingParam1.SetStringVaule(_T("VFile"), argv[4]);
            if (argc > 5) workingParam1.SetStringVaule(_T("OutputFile"), argv[5]);
            res = Task1EOJob(workingParam1);
        }
        break;
    case EHJ_Task2Full:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobTask2Full"));
            if (argc > 3) workingParam1.SetStringVaule(_T("CfgFile"), argv[3]);
            if (argc > 4) workingParam1.SetStringVaule(_T("VFile"), argv[4]);
            if (argc > 5) workingParam1.SetStringVaule(_T("OutputFile"), argv[5]);
            res = Task2FullJob(workingParam1);
        }
        break;
    case EHJ_Task2EO:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobTask2EO"));
            if (argc > 3) workingParam1.SetStringVaule(_T("CfgFile"), argv[3]);
            if (argc > 4) workingParam1.SetStringVaule(_T("VFile"), argv[4]);
            if (argc > 5) workingParam1.SetStringVaule(_T("OutputFile"), argv[5]);
            res = Task2EOJob(workingParam1);
        }
        break;
    case EHJ_Task3EO:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobTask3EO"));
            if (argc > 3) workingParam1.SetStringVaule(_T("CfgFile"), argv[3]);
            if (argc > 4) workingParam1.SetStringVaule(_T("VFile"), argv[4]);
            if (argc > 5) workingParam1.SetStringVaule(_T("OutputFile"), argv[5]);
            res = Task3EOJob(workingParam1);
        }
        break;
    case EHJ_Task4Action:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobTask4Action"));
            if (argc > 3) workingParam1.SetStringVaule(_T("CfgFile"), argv[3]);
            if (argc > 4) workingParam1.SetStringVaule(_T("OutputFile"), argv[4]);
            res = Task4ActionJob(workingParam1);
        }
        break;
    case EHJ_Task5Plaq:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobTask5Plaq"));
            if (argc > 3) workingParam1.SetStringVaule(_T("CfgFile"), argv[3]);
            if (argc > 4) workingParam1.SetStringVaule(_T("OutputFile"), argv[4]);
            res = Task5PlaqJob(workingParam1);
        }
        break;
    case EHJ_Task6FermionAction:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobTask6FermionAction"));
            if (argc > 3) workingParam1.SetStringVaule(_T("CfgFile"), argv[3]);
            if (argc > 4) workingParam1.SetStringVaule(_T("V1File"), argv[4]);
            if (argc > 5) workingParam1.SetStringVaule(_T("V2File"), argv[5]);
            if (argc > 6) workingParam1.SetStringVaule(_T("V3File"), argv[6]);
            if (argc > 7) workingParam1.SetStringVaule(_T("OutputFile"), argv[7]);
            res = Task6FermionActionJob(workingParam1);
        }
        break;
    case EHJ_Task7Mom:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobTask7Mom"));
            if (argc > 3) workingParam1.SetStringVaule(_T("CfgFile"), argv[3]);
            if (argc > 4) workingParam1.SetStringVaule(_T("MomFile"), argv[4]);
            if (argc > 5) workingParam1.SetStringVaule(_T("OutputFile"), argv[5]);
            res = Task7MomJob(workingParam1);
        }
        break;
    case EHJ_TaskRotD:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobTaskRotD"));
            if (argc > 3) workingParam1.SetStringVaule(_T("CfgFile"), argv[3]);
            if (argc > 4) workingParam1.SetStringVaule(_T("VFile"), argv[4]);
            if (argc > 5) workingParam1.SetStringVaule(_T("OutputFile"), argv[5]);
            res = TaskRotDJob(workingParam1);
        }
        break;
    case EHJ_TaskRotFinalForce:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobTaskRotFinalForce"));
            if (argc > 3) workingParam1.SetStringVaule(_T("CfgFile"), argv[3]);
            if (argc > 4) workingParam1.SetStringVaule(_T("VFile"), argv[4]);
            if (argc > 5) workingParam1.SetStringVaule(_T("OutputFile"), argv[5]);
            res = TaskRotFinalForceJob(workingParam1);
        }
        break;
    case EHJ_TaskRotGaugeForce:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobTaskRotGaugeForce"));
            if (argc > 3) workingParam1.SetStringVaule(_T("CfgFile"), argv[3]);
            if (argc > 4) workingParam1.SetStringVaule(_T("OutputFile"), argv[4]);
            res = TaskRotGaugeForceJob(workingParam1);
        }
        break;
    case EHJ_SimulateTestCacheFull:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobSimulateTestCacheFull"));
            res = SimulateTest(workingParam1);
        }
        break;
    case EHJ_SimulateTestCacheMedian:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobSimulateTestCacheMedian"));
            res = SimulateTest(workingParam1);
        }
        break;
    case EHJ_SimulateTestCacheNone:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobSimulateTestCacheNone"));
            res = SimulateTest(workingParam1);
        }
        break;
    case EHJ_DebugSU3_12:
        {
            CParameters workingParam1 = params.GetParameter(_T("JobDebugSU3_12"));
            res = DebugSU3_12Job(workingParam1);
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
