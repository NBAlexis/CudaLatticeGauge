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
CCString sFileNameWrite##lstName##All = _T("%s_angular"); \
sFileNameWrite##lstName = sFileNameWrite##lstName + _T(#lstName) + _T("_Nt%d_O%d.csv"); \
sFileNameWrite##lstName##In = sFileNameWrite##lstName##In + _T(#lstName) + _T("_Nt%d_In_O%d.csv"); \
sFileNameWrite##lstName##All = sFileNameWrite##lstName##All + _T(#lstName) + _T("_Nt%d_All_O%d.csv"); \
sFileNameWrite##lstName.Format(sFileNameWrite##lstName, sCSVSavePrefix.c_str(), _HC_Lt, uiOmega); \
sFileNameWrite##lstName##In.Format(sFileNameWrite##lstName##In, sCSVSavePrefix.c_str(), _HC_Lt, uiOmega); \
sFileNameWrite##lstName##All.Format(sFileNameWrite##lstName##All, sCSVSavePrefix.c_str(), _HC_Lt, uiOmega); \
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
WriteStringFileRealArray2(sFileNameWrite##lstName, lstName##OverR); \
WriteStringFileRealArray(sFileNameWrite##lstName##In, lstName##In); \
WriteStringFileRealArray(sFileNameWrite##lstName##All, lstName##All);


enum ECondList
{
    Chiral,
    Gamma1,
    Gamma2,
    Gamma3,
    Gamma4,
    Gamma5,
    Gamma45,
    GammaX,
    GammaY,
};


#define _CLG_EXPORT_CHIRAL(measureName, lstName) \
CCString sFileNameWrite##lstName = _T("%s_condensate"); \
CCString sFileNameWrite##lstName##All = _T("%s_condensate"); \
sFileNameWrite##lstName = sFileNameWrite##lstName + _T(#lstName) + _T("_Nt%d_O%d.csv"); \
sFileNameWrite##lstName##All = sFileNameWrite##lstName##All + _T(#lstName) + _T("_Nt%d_All_O%d.csv"); \
sFileNameWrite##lstName.Format(sFileNameWrite##lstName, sCSVSavePrefix.c_str(), _HC_Lt, uiOmega); \
sFileNameWrite##lstName##All.Format(sFileNameWrite##lstName##All, sCSVSavePrefix.c_str(), _HC_Lt, uiOmega); \
TArray<TArray<CLGComplex>> lstName##OverR; \
TArray<CLGComplex> lstName##All; \
for (UINT j = 0; j < (iEndN - iStartN + 1); ++j) \
{ \
    TArray<CLGComplex> thisConfiguration; \
    for (INT i = 0; i < measureName->m_lstR.Num(); ++i) \
    { \
        thisConfiguration.AddItem(measureName->m_lstCond[lstName][j * measureName->m_lstR.Num() + i]); \
    } \
    lstName##OverR.AddItem(thisConfiguration); \
    lstName##All.AddItem(measureName->m_lstCondAll[lstName][j]); \
} \
appGeneral(_T("Saving to %s and %s"), sFileNameWrite##lstName, sFileNameWrite##lstName##All); \
WriteStringFileComplexArray2(sFileNameWrite##lstName, lstName##OverR); \
WriteStringFileComplexArray(sFileNameWrite##lstName##All, lstName##All); 


//=============================================================================
// END OF FILE
//=============================================================================
