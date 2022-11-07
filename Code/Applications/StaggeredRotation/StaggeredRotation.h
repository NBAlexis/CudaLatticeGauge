//=============================================================================
// FILENAME : StaggeredRotation.h
// 
// DESCRIPTION:
//   This is to compress the configuration files to half the size
//
// REVISION:
//  [09/24/2020 nbale]
//=============================================================================
#pragma once
#include "CLGLib.h"

__DEFINE_ENUM(EStaggeredRotationJob,
    ESR_Simulate,
    ESR_SimulateEM,
    ESR_SimulateU1,
    ESR_SimulateREM,
    ESR_Measure,
    ESR_MeasureEM,
    ESR_MeasureU1,
    ESR_MeasureREM,
    ESR_SimulateQuenched,
    ESR_SimulateQuenchedU1,
    ESR_SimulateQuenched3D,
    ESR_SimulateQuenchedU13D,
    ESR_Measure3D,
    ESR_MeasureU13D,
    ESR_GaugeFixing,
    )


#define _CLG_EXPORT_CHIRAL(measureName, lstName, variableName, fileIdxHead) \
CCString sFileNameWrite##measureName##lstName = _T("%s_condensate"); \
CCString sFileNameWrite##measureName##lstName##All = _T("%s_condensate"); \
CCString sFileNameWrite##measureName##lstName##In = _T("%s_condensate"); \
sFileNameWrite##measureName##lstName = sFileNameWrite##measureName##lstName + _T(#measureName) + _T(#lstName) + _T("_Nt%d_") + _T(#fileIdxHead) + _T("%d.csv"); \
sFileNameWrite##measureName##lstName##All = sFileNameWrite##measureName##lstName##All + _T(#measureName) + _T(#lstName) + _T("_Nt%d_All_") + _T(#fileIdxHead) + _T("%d.csv"); \
sFileNameWrite##measureName##lstName##In = sFileNameWrite##measureName##lstName##In + _T(#measureName) + _T(#lstName) + _T("_Nt%d_In_") + _T(#fileIdxHead) + _T("%d.csv"); \
sFileNameWrite##measureName##lstName.Format(sFileNameWrite##measureName##lstName, sCSVSavePrefix.c_str(), _HC_Lt, variableName); \
sFileNameWrite##measureName##lstName##All.Format(sFileNameWrite##measureName##lstName##All, sCSVSavePrefix.c_str(), _HC_Lt, variableName); \
sFileNameWrite##measureName##lstName##In.Format(sFileNameWrite##measureName##lstName##In, sCSVSavePrefix.c_str(), _HC_Lt, variableName); \
TArray<TArray<CLGComplex>> lstName##measureName##OverR; \
TArray<CLGComplex> lstName##measureName##All; \
TArray<CLGComplex> lstName##measureName##In; \
for (UINT j = 0; j < (iEndN - iStartN + 1); ++j) \
{ \
    TArray<CLGComplex> thisConfiguration; \
    for (INT i = 0; i < measureName->m_lstR.Num(); ++i) \
    { \
        thisConfiguration.AddItem(measureName->m_lstCond[lstName][j * measureName->m_lstR.Num() + i]); \
    } \
    lstName##measureName##OverR.AddItem(thisConfiguration); \
    lstName##measureName##All.AddItem(measureName->m_lstCondAll[lstName][j]); \
    lstName##measureName##In.AddItem(measureName->m_lstCondIn[lstName][j]); \
} \
WriteStringFileComplexArray2(sFileNameWrite##measureName##lstName, lstName##measureName##OverR); \
WriteStringFileComplexArray(sFileNameWrite##measureName##lstName##All, lstName##measureName##All); \
WriteStringFileComplexArray(sFileNameWrite##measureName##lstName##In, lstName##measureName##In); 



#define _CLG_EXPORT_CHIRALZSLICE(measureName, lstName, variableName, fileIdxHead) \
CCString sFileNameWrite##measureName##lstName##ZSlice = _T("%s_condensate"); \
sFileNameWrite##measureName##lstName##ZSlice = sFileNameWrite##measureName##lstName##ZSlice + _T(#measureName) + _T(#lstName) + _T("ZSlice_Nt%d_") + _T(#fileIdxHead) + _T("%d.csv"); \
sFileNameWrite##measureName##lstName##ZSlice.Format(sFileNameWrite##measureName##lstName##ZSlice, sCSVSavePrefix.c_str(), _HC_Lt, variableName); \
TArray<TArray<CLGComplex>> lstName##measureName##ZSlice; \
for (UINT j = 0; j < (iEndN - iStartN + 1); ++j) \
{ \
    TArray<CLGComplex> thisConfiguration; \
    for (UINT i = 0; i < _HC_Lz; ++i) \
    { \
        thisConfiguration.AddItem(measureName->m_lstCondZSlice[lstName][j * _HC_Lz + i]); \
    } \
    lstName##measureName##ZSlice.AddItem(thisConfiguration); \
} \
WriteStringFileComplexArray2(sFileNameWrite##measureName##lstName##ZSlice, lstName##measureName##ZSlice); 


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
WriteStringFileRealArray2(sFileNameWrite##lstName, lstName##OverR); \
WriteStringFileRealArray(sFileNameWrite##lstName##In, lstName##In); \
WriteStringFileRealArray(sFileNameWrite##lstName##All, lstName##All);


#if !_CLG_WIN
inline void strerror_s(TCHAR* buffer, size_t bufferSize, INT error)
{
    strcpy(buffer, strerror(error));
}
#endif

enum { kExportDigital = 20, };

#if !_CLG_WIN

inline void _gcvt_s(TCHAR* buff, UINT uiBuffLength, Real fVaule, UINT uiDigit)
{
    static TCHAR tmpBuff[10];
    appSprintf(tmpBuff, 10, _T("%s.%df"), _T("%"), uiDigit);
    appSprintf(buff, uiBuffLength, tmpBuff, fVaule);
}

#endif

inline void WriteStringFile(const CCString& sFileName, const CCString& sContent)
{
    appGetFileSystem()->WriteAllText(sFileName, sContent);
}

template <class T>
void WriteStringFileRealArray(const CCString& sFileName, const TArray<T>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kExportDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.Num(); ++i)
    {
        _gcvt_s(str, 50, static_cast<DOUBLE>(lst[i]), iDigital);
        CCString sReal = CCString(str);
        sReal = sReal.Replace(_T("e"), _T("*^"));
        file << _T(" ");
        file << sReal;
        if (i != lst.GetCount() - 1)
        {
            file << _T(",");
        }
    }
    file.flush();
    file.close();
}

inline void WriteStringFileRealArray2(const CCString& sFileName, const TArray<TArray<Real>>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kExportDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        for (INT j = 0; j < lst[i].GetCount(); ++j)
        {
            _gcvt_s(str, 50, lst[i][j], iDigital);
            CCString sReal = CCString(str);
            sReal = sReal.Replace(_T("e"), _T("*^"));
            file << _T(" ");
            file << sReal;
            if (j != lst[i].GetCount() - 1)
            {
                file << _T(",");
            }
        }
        file << _T("\n");
    }
    file.flush();
    file.close();
}

inline void WriteStringFileComplexArray(const CCString& sFileName, const TArray<CLGComplex>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kExportDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.Num(); ++i)
    {
        _gcvt_s(str, 50, lst[i].x, iDigital);
        CCString sReal = CCString(str);
        sReal = sReal.Replace(_T("e"), _T("*^"));
        _gcvt_s(str, 50, lst[i].y, iDigital);
        CCString sImg = CCString(str);
        sImg = sImg.Replace(_T("e"), _T("*^"));
        CCString sMid = _T(" + ");
        if (sImg.Left(1) == _T("-"))
        {
            sImg = sImg.Right(sImg.GetLength() - 1);
            sMid = _T(" - ");
        }

        file << _T(" ");
        file << sReal;
        file << sMid;
        file << sImg;
        if (i == lst.GetCount() - 1)
        {
            file << _T(" I");
        }
        else
        {
            file << _T(" I,");
        }
    }
    file.flush();
    file.close();
}

inline void WriteStringFileComplexArray2(const CCString& sFileName, const TArray<TArray<CLGComplex>>& lst, UBOOL bAppend = FALSE)
{
    const INT iDigital = static_cast<INT>(kExportDigital);
    std::ofstream file;
    if (!bAppend)
    {
        file.open(sFileName.c_str(), std::ios::out);
    }
    else
    {
        file.open(sFileName.c_str(), std::ios::app | std::ios::out);
    }

    if (file.fail())
    {
        static TCHAR errorMsg[256];
        strerror_s(errorMsg, 256, errno);
        appCrucial(_T("Saving %s failed! Because %s\n"), sFileName.c_str(), errorMsg);
    }

    TCHAR str[50];
    for (INT i = 0; i < lst.GetCount(); ++i)
    {
        for (INT j = 0; j < lst[i].GetCount(); ++j)
        {
            _gcvt_s(str, 50, lst[i][j].x, iDigital);
            CCString sReal = CCString(str);
            sReal = sReal.Replace(_T("e"), _T("*^"));
            _gcvt_s(str, 50, lst[i][j].y, iDigital);
            CCString sImg = CCString(str);
            sImg = sImg.Replace(_T("e"), _T("*^"));
            CCString sMid = _T(" + ");
            if (sImg.Left(1) == _T("-"))
            {
                sImg = sImg.Right(sImg.GetLength() - 1);
                sMid = _T(" - ");
            }
            file << _T(" ");
            file << sReal;
            file << sMid;
            file << sImg;
            if (j == lst[i].GetCount() - 1)
            {
                file << _T(" I");
            }
            else
            {
                file << _T(" I,");
            }
        }
        file << _T("\n");
    }
    file.flush();
    file.close();
}

inline void AppendStringFile(const CCString& sFileName, const CCString& sContent)
{
    appGetFileSystem()->AppendAllText(sFileName, sContent);
}


extern INT SimulateStaggeredBeta(CParameters& params);
extern INT SimulateStaggeredRotation(CParameters& params);
extern INT SimulateStaggeredRotation3D(CParameters& params);
extern INT SimulateStaggeredEM(CParameters& params);
extern INT SimulateStaggeredRotationU1(CParameters& params);
extern INT SimulateStaggeredRotationU13D(CParameters& params);
extern INT SimulateStaggeredRotationEM(CParameters& params);
extern INT Measurement(CParameters& params);
extern INT MeasurementEM(CParameters& params);
extern INT MeasurementU1(CParameters& params);
extern INT MeasurementREM(CParameters& params);
extern INT Measurement3D(CParameters& params);
extern INT GaugeFixing(CParameters& params);

//=============================================================================
// END OF FILE
//=============================================================================
