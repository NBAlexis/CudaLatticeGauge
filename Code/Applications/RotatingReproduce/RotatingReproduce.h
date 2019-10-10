//=============================================================================
// FILENAME : RotatingReproduce.h
// 
// DESCRIPTION:
// Reproduce the results of PRL
//
// REVISION:
//  [02/26/2019 nbale]
//=============================================================================

#include "CLGLib.h"

__DEFINE_ENUM(ERotatingJob,
    ERJ_AngularMomentum,
    ERJ_Thermal,
    ERJ_PolyakovDist,
    ERJ_GaugeFixing,
    )


extern INT TestAngularMomentum(CParameters& params);
extern INT TestThermal(CParameters& params);
extern INT MeasurePolyakovDist(CParameters& params);
extern INT GaugeFixing(CParameters& params);


#define _CLG_EXPORT_ANGULAR(measureName, lstName) \
CCString sFileNameWrite##lstName = _T("%s_angular"); \
CCString sFileNameWrite##lstName##In = _T("%s_angular"); \
sFileNameWrite##lstName = sFileNameWrite##lstName + _T(#lstName) + _T("_Nt%d_O%d.csv"); \
sFileNameWrite##lstName##In = sFileNameWrite##lstName##In + _T(#lstName) + _T("_Nt%d_In.csv"); \
sFileNameWrite##lstName.Format(sFileNameWrite##lstName, sCSVSavePrefix.c_str(), _HC_Lt, uiOmega); \
sFileNameWrite##lstName##In.Format(sFileNameWrite##lstName##In, sCSVSavePrefix.c_str(), _HC_Lt); \
TArray<TArray<Real>> lstName##OverR; \
TArray<Real> lstName##In; \
for (UINT j = 0; j < (iEndN - iStartN + 1); ++j) \
{ \
    TArray<Real> thisConfiguration; \
    for (INT i = 0; i < measureName->m_lstR.Num(); ++i) \
    { \
        thisConfiguration.AddItem(measureName->m_lst##lstName[j * measureName->m_lstR.Num() + i]); \
    } \
    lstName##OverR.AddItem(thisConfiguration); \
    lstName##In.AddItem(measureName->m_lst##lstName##Inner[j]); \
} \
WriteStringFile(sFileNameWrite##lstName, ExportRealArray2(lstName##OverR)); \
CCString sHead##lstName; \
sHead##lstName.Format("\n\n\n(* Omega = %d *)\n\n\n", uiOmega); \
AppendStringFile(sFileNameWrite##lstName##In, sHead##lstName); \
AppendStringFile(sFileNameWrite##lstName##In, ExportRealArray(lstName##In));


//=============================================================================
// END OF FILE
//=============================================================================
