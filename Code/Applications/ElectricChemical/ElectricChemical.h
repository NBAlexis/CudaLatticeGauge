//=============================================================================
// FILENAME : ElectricChemical.h
// 
// DESCRIPTION:
//   
//
// REVISION:
//  [mm/dd/yy]
//  [09/29/2022 nbale]
//=============================================================================
#pragma once

#include "CLGLib.h"

__DEFINE_ENUM(EElectricChemical,
    EEC_Simulate,
    EEC_Measure,
    EEC_GaugeFixing,
    EEC_SimulateRW,
    EEC_MeasureRW,
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
WriteComplexArray2(sFileNameWrite##lstName, lstName##OverR); \
WriteComplexArray(sFileNameWrite##lstName##In, lstName##In); \
WriteComplexArray(sFileNameWrite##lstName##All, lstName##All);

inline void AppendStringFile(const CCString& sFileName, const CCString& sContent)
{
    appGetFileSystem()->AppendAllText(sFileName, sContent);
}


extern INT Simulate(CParameters& params);
extern INT Measurement(CParameters& params);
extern INT GaugeFixing(CParameters& params);
extern INT SimulateRW(CParameters& params);
extern INT MeasureRW(CParameters& params);

//=============================================================================
// END OF FILE
//=============================================================================
