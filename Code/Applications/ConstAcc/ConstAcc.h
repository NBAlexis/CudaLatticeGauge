//=============================================================================
// FILENAME : ConstAcc.h
// 
// DESCRIPTION:
//   This is to compress the configuration files to half the size
//
// REVISION:
//  [04/10/2020 nbale]
//=============================================================================

#include "CLGLib.h"

__DEFINE_ENUM(EAccJob,
    EAJ_Simulate,
    EAJ_SimulateQ,
    EAJ_SimulateQMidCenter,
    EAJ_Measure,
    EAJ_MeasureMidCenter,
    )

#define _CLG_EXPORT_CHIRAL(measureName, lstName) \
CCString sFileNameWrite##measureName##lstName = _T("%s_%s_condensate"); \
CCString sFileNameWrite##measureName##lstName##ZSlice = _T("%s_%s_condensateZSlice"); \
sFileNameWrite##measureName##lstName = sFileNameWrite##measureName##lstName + _T(#measureName) + _T(#lstName) + _T(".csv"); \
sFileNameWrite##measureName##lstName##ZSlice = sFileNameWrite##measureName##lstName##ZSlice + _T(#measureName) + _T(#lstName) + _T(".csv"); \
sFileNameWrite##measureName##lstName.Format(sFileNameWrite##measureName##lstName, sCSVSavePrefix.c_str(), PrefixList[uiOmega].c_str()); \
sFileNameWrite##measureName##lstName##ZSlice.Format(sFileNameWrite##measureName##lstName##ZSlice, sCSVSavePrefix.c_str(), PrefixList[uiOmega].c_str()); \
TArray<CLGComplex> lstName##measureName; \
TArray<TArray<CLGComplex>> lstName##measureName##ZSlice; \
for (UINT j = 0; j < (iEndN - iStartN + 1); ++j) \
{ \
    lstName##measureName.AddItem(measureName->m_lstCondAll[lstName][j]); \
    if (measureName->m_bMeasureZSlice) \
    { \
        TArray<CLGComplex> thisConfiguration##measureName##lstName##ZSlice; \
        for (UINT i = 0; i < _HC_Lz; ++i) \
        { \
            thisConfiguration##measureName##lstName##ZSlice.AddItem(measureName->m_lstCondZSlice[lstName][j * _HC_Lz + i]); \
        } \
        lstName##measureName##ZSlice.AddItem(thisConfiguration##measureName##lstName##ZSlice); \
    } \
} \
WriteComplexArray(sFileNameWrite##measureName##lstName, lstName##measureName); \
WriteComplexArray2(sFileNameWrite##measureName##lstName##ZSlice, lstName##measureName##ZSlice); 

extern INT SimulateAcc(CParameters& params);
extern INT Measurement(CParameters& params);

//=============================================================================
// END OF FILE
//=============================================================================
