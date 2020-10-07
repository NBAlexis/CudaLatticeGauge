//=============================================================================
// FILENAME : Measurement.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [09/25/2020 nbale]
//=============================================================================

#include "StaggeredSpectrum.h"

enum { kExportDigital = 20, };

#if !_CLG_WIN

void _gcvt_s(TCHAR* buff, UINT uiBuffLength, Real fVaule, UINT uiDigit)
{
    static TCHAR tmpBuff[10];
    appSprintf(tmpBuff, 10, _T("%s.%df"), _T("%"), uiDigit);
    appSprintf(buff, uiBuffLength, tmpBuff, fVaule);
}

#endif

void WriteStringFile(const CCString& sFileName, const CCString& sContent)
{
    appGetFileSystem()->WriteAllText(sFileName, sContent);
}

void WriteStringFileRealArray(const CCString& sFileName, const TArray<Real>& lst, UBOOL bAppend = FALSE)
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
    TCHAR str[50];
    for (INT i = 0; i < lst.Num(); ++i)
    {
        _gcvt_s(str, 50, lst[i], iDigital);
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

void WriteStringFileDoubleColumn(const CCString& sFileName, const TArray<Real>& lst1, const TArray<CLGComplex>& lst2, UBOOL bAppend = FALSE)
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
    TCHAR strleft[50];
    TCHAR strreal[50];
    TCHAR strimg[50];
    for (INT i = 0; i < lst1.Num(); ++i)
    {
        _gcvt_s(strleft, 50, lst1[i], iDigital);
        _gcvt_s(strreal, 50, lst2[i].x, iDigital);
        _gcvt_s(strimg, 50, lst2[i].y, iDigital);
        CCString s1 = CCString(strleft);
        CCString s2 = CCString(strreal);
        CCString s3 = CCString(strimg);
        s1 = s1.Replace(_T("e"), _T("*^"));
        s2 = s2.Replace(_T("e"), _T("*^"));
        s3 = s3.Replace(_T("e"), _T("*^"));
        
        file << s1;
        file << _T(", ");
        file << s2;

        if (s3.Left(1) == _T("-"))
        {
            s3 = s3.Right(s3.GetLength() - 1);
            file << _T(" - ");
        }
        else
        {
            file << _T(" + ");
        }
        file << s3;
        file << _T(" I");
    }
    file.flush();
    file.close();
}

void WriteStringFileRealArray2(const CCString& sFileName, const TArray<TArray<Real>>& lst, UBOOL bAppend = FALSE)
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

void WriteStringFileComplexArray(const CCString& sFileName, const TArray<CLGComplex>& lst, UBOOL bAppend = FALSE)
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

void WriteStringFileComplexArray2(const CCString& sFileName, const TArray<TArray<CLGComplex>>& lst, UBOOL bAppend = FALSE)
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

void AppendStringFile(const CCString& sFileName, const CCString& sContent)
{
    appGetFileSystem()->AppendAllText(sFileName, sContent);
}

INT Measurement(CParameters& params)
{

#pragma region read parameters

    appSetupLog(params);

    INT iVaule = 1;
    params.FetchValueINT(_T("StartN"), iVaule);
    UINT iStartN = static_cast<UINT>(iVaule);

    iVaule = 200;
    params.FetchValueINT(_T("EndN"), iVaule);
    UINT iEndN = static_cast<UINT>(iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("DoSmearing"), iVaule);
    UBOOL bDoSmearing = (0 != iVaule);

    //iVaule = 1;
    //params.FetchValueINT(_T("CheckGaugeFixing"), iVaule);
    //UBOOL bCheckGaugeFixing = 0 != iVaule;

    //iVaule = 0;
    //params.FetchValueINT(_T("UseZ4"), iVaule);
    //UBOOL bZ4 = 0 != iVaule;

    CCString sValue = _T("ESSM_Polyakov");
    params.FetchStringValue(_T("MeasureType"), sValue);
    EStaggeredSpectrumMeasure eJob = __STRING_TO_ENUM(EStaggeredSpectrumMeasure, sValue);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    CCString sCSVSavePrefix;
    params.FetchStringValue(_T("CSVSavePrefix"), sCSVSavePrefix);
    appGeneral(_T("csv save prefix: %s\n"), sCSVSavePrefix.c_str());

    CCString sSubFolderPrefix;
    params.FetchStringValue(_T("SubFolderPrefix"), sSubFolderPrefix);
    appGeneral(_T("sub folder prefix: %s\n"), sSubFolderPrefix.c_str());

    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

#pragma endregion

    const UINT uiNewLine = (iEndN - iStartN + 1) / 5;
    CFieldGaugeSU3* pStaple = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField->GetCopy());
    CMeasureWilsonLoop* pPL = dynamic_cast<CMeasureWilsonLoop*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureMesonCorrelatorStaggered* pMC = dynamic_cast<CMeasureMesonCorrelatorStaggered*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CMeasureMesonCorrelatorStaggeredSimple* pMCSimple = dynamic_cast<CMeasureMesonCorrelatorStaggeredSimple*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
    pPL->Reset();
    pMC->Reset();

#pragma region Measure

    appGeneral(_T("(*\n"));
    appSetLogDate(FALSE);
    for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
    {
        CCString sFileName;
        sFileName.Format(_T("%sMatching_%d.con"), sSavePrefix.c_str(), uiN);
        appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, EFFT_CLGBin);

        switch (eJob)
        {
        case ESSM_Polyakov:
        {
            if (bDoSmearing)
            {
                appGetLattice()->m_pGaugeField->CalculateOnlyStaple(pStaple);
                appGetLattice()->m_pGaugeSmearing->GaugeSmearing(appGetLattice()->m_pGaugeField, pStaple);
            }

            pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
            if (uiN == iStartN)
            {
                TArray<Real> lstRadius;
                for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                {
                    lstRadius.AddItem(_hostsqrt(static_cast<Real>(pPL->m_lstR[i])));
                }
                CCString sRadiousFile;
                sRadiousFile.Format(_T("%s_VR_R.csv"), sCSVSavePrefix.c_str());
                WriteStringFileRealArray(sRadiousFile, lstRadius);
            }
        }
        break;
        case ESSM_Correlator:
        {
            pMC->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
        }
        break;
        case ESSM_CorrelatorSimple:
        {
            pMCSimple->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
        }
        break;
        default:
            break;
        }

        if ((iEndN - uiN + 1) % uiNewLine == 0)
        {
            appSetLogDate(TRUE);
            appGeneral(_T("\n="));
            appSetLogDate(FALSE);
        }
        else
        {
            appSetLogDate(FALSE);
            appGeneral(_T("="));
        }

    }
    appGeneral(_T("\n*)\n"));

#pragma endregion

    switch (eJob)
    {
    case ESSM_Polyakov:
    {
        //Write result to file
        CCString sCSVFile;
        sCSVFile.Format(_T("%s_VR.csv"), sCSVSavePrefix.c_str());
        TArray<TArray<CLGComplex>> vrs;
        for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
        {
            TArray<CLGComplex> thisConfiguration;
            for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
            {
                for (UINT t = 0; t < _HC_Lt / 2; ++t)
                {
                    thisConfiguration.AddItem(pPL->m_lstC[j][i][t]);
                }
            }
            vrs.AddItem(thisConfiguration);
        }
        WriteStringFileComplexArray2(sCSVFile, vrs);
    }
    break;
    case ESSM_Correlator:
    {
        for (INT ty = 0; ty < CMeasureMesonCorrelatorStaggered::_kMesonCorrelatorType; ++ty)
        {
            CCString sCSVFile;
            sCSVFile.Format(_T("%s_meson%d.csv"), sCSVSavePrefix.c_str(), ty);
            TArray<TArray<Real>> res;
            for (INT conf = 0; conf < pMC->m_lstResults.Num(); ++conf)
            {
                TArray<Real> oneConf;
                for (INT t = 0; t < _HC_Lti - 1; ++t)
                {
                    oneConf.AddItem(pMC->m_lstResults[conf][ty][t].x);
                }
                res.AddItem(oneConf);
            }
            WriteStringFileRealArray2(sCSVFile, res);
        }
    }
    break;
    case ESSM_CorrelatorSimple:
    {
        for (INT ty = 0; ty < CMeasureMesonCorrelatorStaggeredSimple::_kMesonCorrelatorTypeSimple; ++ty)
        {
            CCString sCSVFile;
            sCSVFile.Format(_T("%s_mesonsimple%d.csv"), sCSVSavePrefix.c_str(), ty);
            TArray<TArray<Real>> res;
            for (INT conf = 0; conf < pMCSimple->m_lstResults.Num(); ++conf)
            {
                TArray<Real> oneConf;
                for (INT t = 0; t < _HC_Lti - 1; ++t)
                {
                    oneConf.AddItem(pMCSimple->m_lstResults[conf][ty][t]);
                }
                res.AddItem(oneConf);
            }
            WriteStringFileRealArray2(sCSVFile, res);
        }
    }
    break;
    default:
        break;
    }

    appGeneral(_T("\n"));

    appGeneral(_T("\n(*"));
    appSetLogDate(TRUE);

    appGeneral(_T("\n=====================================\n========= finished! ==========\n*)"));

    appSafeDelete(pStaple);

    appQuitCLG();

    return 0;
}


