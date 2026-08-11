//=============================================================================
// FILENAME : BetaGradientJob.h
// 
// DESCRIPTION:
//   
//
// REVISION:
//  [mm/dd/yy]
//  [08/17/2022 nbale]
//=============================================================================
#pragma once

#include "CLGLib.h"

__DEFINE_ENUM(EBetaGradientJob,
    EBGJ_Simulate,
    EBGJ_SimulateScanBeta,
    EBGJ_Measure,
    EBGJ_MeasureScanBeta,
    EBGJ_SimulateQ,
    EBGJ_SimulateScanQ,
    EBGJ_GaugeFixing,
    EBGJ_MCGGaugeFixing,
    EBGJ_Simulate2,
    EBGJ_Measure2,
    EBGJ_SimulateAtGradient,
    EBGJ_SimulateTemperatureDistri,
    EBGJ_MeasureTemperatureDistri,
    )



#define _CLG_EXPORT_ANGULAR(measureName, lstName, variableName, fileIdxHead) \
CCString sFileNameWrite##lstName = _T("%s_angular"); \
CCString sFileNameWrite##lstName##In = _T("%s_angular"); \
CCString sFileNameWrite##lstName##All = _T("%s_angular"); \
sFileNameWrite##lstName = sFileNameWrite##lstName + _T(#lstName) + _T("_Nt%d_") + _T(#fileIdxHead) + _T("%d.csv"); \
sFileNameWrite##lstName##In = sFileNameWrite##lstName##In + _T(#lstName) + _T("_Nt%d_In_") + _T(#fileIdxHead) + _T("%d.csv"); \
sFileNameWrite##lstName##All = sFileNameWrite##lstName##All + _T(#lstName) + _T("_Nt%d_All_") + _T(#fileIdxHead) + _T("%d.csv"); \
sFileNameWrite##lstName.Format(sFileNameWrite##lstName, sCSVSavePrefix.c_str(), _HC_Lt, variableName); \
sFileNameWrite##lstName##In.Format(sFileNameWrite##lstName##In, sCSVSavePrefix.c_str(), _HC_Lt, variableName); \
sFileNameWrite##lstName##All.Format(sFileNameWrite##lstName##All, sCSVSavePrefix.c_str(), _HC_Lt, variableName); \
TArray<TArray<Real>> lstName##OverR; \
TArray<Real> lstName##In; \
TArray<Real> lstName##All; \
for (UINT j = 0; j < (iEndN - iStartN + 1); ++j) \
{ \
    TArray<Real> thisConfiguration; \
    for (INT i = 0; i < measureName->m_lstR.Num(); ++i) \
    { \
        thisConfiguration.AddItem(measureName->m_lst##lstName[j * measureName->m_lstR.Num() + i]); \
    } \
    lstName##OverR.AddItem(thisConfiguration); \
    lstName##In.AddItem(measureName->m_lst##lstName##Inner[j]); \
    lstName##All.AddItem(measureName->m_lst##lstName##All[j]); \
} \
WriteRealArray2(sFileNameWrite##lstName, lstName##OverR); \
WriteRealArray(sFileNameWrite##lstName##In, lstName##In); \
WriteRealArray(sFileNameWrite##lstName##All, lstName##All);

extern INT Simulate(CParameters& params);
extern INT Measurement(CParameters& params);
extern INT SimulateBetaScan(CParameters& params);
extern INT MeasurementBetaScan(CParameters& params);
extern INT GaugeFixing(CParameters& params);
extern INT MCGGaugeFixing(CParameters& params);
extern INT SimulateAtGradient(CParameters& params);
extern INT SimulateTempDist(CParameters& params);
extern INT MeasurementTemperatureDistri(CParameters & params);

extern TArray<TArray<INT>> ParseWilsonPath(const CCString& sFileName);
extern TArray<SCHAR> GetOnePath(const TArray<INT>&dirs, SCHAR mu, SCHAR nu);

//=============================================================================
// END OF FILE
//=============================================================================
