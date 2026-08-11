//=============================================================================
// FILENAME : HISQ.h
// 
// DESCRIPTION:
//   
//
// REVISION:
//  [mm/dd/yy]
//  [12/08/2024 nbale]
//=============================================================================
#pragma once
#define _CLG_PROFILER 1
#include "CLGLib.h"

__DEFINE_ENUM(EHISQJob,
    EHJ_SimulateTest,
    EHJ_SimulateNf2,
    EHJ_SimulateNf2MassPrecondition,
    EHJ_SimulateNf2p1,
    EHJ_Measure,
    EHJ_MeasureWilson, //To put the definition here because "EHMJ_Wilson" needs another YAML parameter for gauge smearing
    EHJ_GaugeFixing,
    EHJ_Task1Full,
    EHJ_Task1EO,
    EHJ_Task2Full,
    EHJ_Task2EO,
    EHJ_Task3EO,
    EHJ_Task4Action,
    EHJ_Task5Plaq,
    EHJ_Task6FermionAction,
    EHJ_Task7Mom,
    EHJ_TaskRotD,
    EHJ_TaskRotFinalForce,
    EHJ_TaskRotGaugeForce,
    EHJ_SimulateTestCacheFull,
    EHJ_SimulateTestCacheMedian,
    EHJ_SimulateTestCacheNone,
    EHJ_DebugSU3_12,
    )


extern INT SimulateTest(CParameters& params);
extern INT SimulateNf2p1(CParameters& params);
extern INT Measurement(CParameters& params);
extern INT GaugeFixing(CParameters& params);
extern INT Task1FullJob(CParameters& params);
extern INT Task1EOJob(CParameters& params);
extern INT Task2FullJob(CParameters& params);
extern INT Task2EOJob(CParameters& params);
extern INT Task3EOJob(CParameters& params);
extern INT Task4ActionJob(CParameters& params);
extern INT Task5PlaqJob(CParameters& params);
extern INT Task6FermionActionJob(CParameters& params);
extern INT Task7MomJob(CParameters& params);
extern INT TaskRotDJob(CParameters& params);
extern INT TaskRotFinalForceJob(CParameters& params);
extern INT TaskRotGaugeForceJob(CParameters& params);
extern INT DebugSU3_12Job(CParameters& params);

// Shared helpers for comparison tasks (defined in CompareHISQ.cpp)
struct SCompareCtx
{
    CCString sCfgFile;
    CCString sVFile;
    CCString sOutFile;
    CFieldGaugeSU3* pGauge;
    CFieldFermionHISQSU3* pFerm;
    CGaugeSmearingHISQSU3* pSmear;
    TArray<const CFieldGauge*> effGauges;
    TArray<const CFieldGauge*> rawGauges;
};

extern INT SetupCompare(CParameters& params, SCompareCtx& ctx, const CCString& sTag);
extern INT SetupGaugeOnly(CParameters& params, SCompareCtx& ctx, const CCString& sTag);


//=============================================================================
// END OF FILE
//=============================================================================
