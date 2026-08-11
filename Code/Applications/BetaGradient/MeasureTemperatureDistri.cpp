//=============================================================================
// FILENAME : MeasurementTemperatureDistri.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [mm/dd/yy]
//  [11/07/2024 nbale]
//=============================================================================

#include "BetaGradient.h"

__DEFINE_ENUM(ETemperatureMeasureJob,
    ETMJ_Polyakov,
    ETMJ_BosonValue,
    ETMJ_Plaquette,
    )

INT MeasurementTemperatureDistri(CParameters& params)
{

#pragma region read parameters

    appSetupLog(params);

    INT iVaule = 0;
    iVaule = 1;
    params.FetchValueINT(_T("StartN"), iVaule);
    UINT iStartN = static_cast<UINT>(iVaule);

    iVaule = 200;
    params.FetchValueINT(_T("EndN"), iVaule);
    UINT iEndN = static_cast<UINT>(iVaule);

    //iVaule = 10;
    //params.FetchValueINT(_T("StochasticFieldCount"), iVaule);
    //UINT iFieldCount = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("AlsoCheckMD5"), iVaule);
    UBOOL bCheckMd5 = (0 != iVaule);

    //iVaule = 0;
    //params.FetchValueINT(_T("UseZ4"), iVaule);
    //UBOOL bZ4 = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("SubFolder"), iVaule);
    UBOOL bSubFolder = 0 != iVaule;

    CCString sSubFolderName;
    params.FetchStringValue(_T("SubFolderPrefix"), sSubFolderName);

    CCString sValue = _T("ETMJ_Polyakov");
    params.FetchStringValue(_T("DistributionJob"), sValue);
    ETemperatureMeasureJob eJob = __STRING_TO_ENUM(ETemperatureMeasureJob, sValue);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    CCString sCSVSavePrefix;
    params.FetchStringValue(_T("CSVSavePrefix"), sCSVSavePrefix);
    appGeneral(_T("csv save prefix: %s\n"), sCSVSavePrefix.c_str());

    iVaule = 0;
    params.FetchValueINT(_T("ListStart"), iVaule);
    UINT iListStart = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("ListEnd"), iVaule);
    const INT iListEnd = iVaule;

    TArray<DOUBLE> BetaList;
    params.FetchValueArrayDOUBLE(_T("BetaList"), BetaList);

    TArray<CCString> PrefixList;
    params.FetchStringVectorValue(_T("PrefixList"), PrefixList);

    if (BetaList.Num() != PrefixList.Num() || BetaList.Num() < 1)
    {
        appCrucial(_T("BetaList and PrefixList not correctly set!\n"));
        return 0;
    }

    /*
    iVaule = 0;
    params.FetchValueINT(_T("SaveFermionFile"), iVaule);
    UBOOL bSaveFermion = 0 != iVaule;

    iVaule = 1;
    params.FetchValueINT(_T("FermionFileIndexStart"), iVaule);
    UINT uiSaveFermionStart = static_cast<UINT>(iVaule);

    CCString sFermionHead;
    params.FetchStringValue(_T("FermionFileHead"), sFermionHead);
    appGeneral(_T("FermionFileHead: %s\n"), sFermionHead.c_str());

    iVaule = 0;
    params.FetchValueINT(_T("LoadFermion"), iVaule);
    //const UINT uiLoadFermion = iVaule;

    CCString sLoadFermionFile;
    params.FetchStringValue(_T("LoadFermionFile"), sLoadFermionFile);
    appGeneral(_T("Load Fermion File Name: %s\n"), sLoadFermionFile.c_str());

    CCString sLoadFermionHead;
    params.FetchStringValue(_T("LoadFermionHead"), sLoadFermionHead);
    appGeneral(_T("Load Fermion File Head: %s\n"), sLoadFermionHead.c_str());
    */

    CCString sLoadType = _T("EFFT_CLGBin");
    EFieldFileType eLoadType = EFFT_CLGBin;
    if (params.FetchStringValue(_T("LoadType"), sLoadType))
    {
        eLoadType = __STRING_TO_ENUM(EFieldFileType, sLoadType);
    }
    appGeneral(_T("load type: %s\n"), __ENUM_TO_STRING(EFieldFileType, eLoadType).c_str());

    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

#pragma endregion


    UINT uiNewLine = (iEndN - iStartN + 1) / 5;
    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureBosonValueReal* pBV = dynamic_cast<CMeasureBosonValueReal*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CMeasureWilsonLoopWithPath* pWilsonPath = dynamic_cast<CMeasureWilsonLoopWithPath*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
    TArray<TArray<SCHAR>> wilsonPaths;
    for (INT i = 0; i < 4; ++i)
    {
        for (INT j = i + 1; j < 4; ++j)
        {
            TArray<SCHAR> onepath;
            onepath.AddItem(static_cast<SCHAR>(i + 1));
            onepath.AddItem(static_cast<SCHAR>(j + 1));
            onepath.AddItem(static_cast<SCHAR>(-i - 1));
            onepath.AddItem(static_cast<SCHAR>(-j - 1));
            wilsonPaths.AddItem(onepath);
        }
    }
    pWilsonPath->SetPath(wilsonPaths);

    appPushLogDate(FALSE);

    for (INT uiOmega = iListStart; uiOmega < BetaList.Num() && uiOmega < iListEnd; ++uiOmega)
    {
        appGeneral(_T("(* ==== Beta(%f) ========= *)\n"), BetaList[uiOmega]);
        pPL->Reset();
        pBV->Reset();
        pWilsonPath->Reset();

#pragma region Measure

        appGeneral(_T("(*"));
        for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
        {
            CCString sFileName;
            CCString sTxtFileName;
            if (bSubFolder)
            {
                sFileName.Format(_T("%s/%s/%sGradient_%s_%dG.con"), sSubFolderName.c_str(), PrefixList[uiOmega].c_str(), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
                sTxtFileName.Format(_T("%s/%s/%sGradient_%s_%dG.txt"), sSubFolderName.c_str(), PrefixList[uiOmega].c_str(), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
            }
            else
            {
                sFileName.Format(_T("%sGradient_%s_%dG.con"), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
                sTxtFileName.Format(_T("%sGradient_%s_%d.txt"), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
            }

            //appGeneral(_T("checking %s ..."), sFileName);
            if (bCheckMd5)
            {
                UINT uiSize = 0;
                BYTE* fileContent = appGetFileSystem()->ReadAllBytes(sFileName, uiSize);
                CCString sMD5 = "MD5 : " + CLGMD5Hash(fileContent, uiSize);
                CCString sFileContent = appGetFileSystem()->ReadAllText(sTxtFileName);
                if (sFileContent.Find(sMD5) >= 0)
                {
                    appGeneral(_T("-"));
                }
                else if (sFileContent.Find("MD5 : ") >= 0)
                {
                    appCrucial(_T("MD5 Found and NOT good %s \n"), sFileName.c_str());
                }
            }
            
            appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName, eLoadType);

            if (bSubFolder)
            {
                sFileName.Format(_T("%s/%s/%sGradient_%s_%dB.con"), sSubFolderName.c_str(), PrefixList[uiOmega].c_str(), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
            }
            else
            {
                sFileName.Format(_T("%sGradient_%s_%dB.con"), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
            }

            //appGeneral(_T("checking %s ..."), sFileName);
            if (bCheckMd5)
            {
                UINT uiSize = 0;
                BYTE* fileContent = appGetFileSystem()->ReadAllBytes(sFileName, uiSize);
                CCString sMD5 = "MD5 : " + CLGMD5Hash(fileContent, uiSize);
                CCString sFileContent = appGetFileSystem()->ReadAllText(sTxtFileName);
                if (sFileContent.Find(sMD5) >= 0)
                {
                    appGeneral(_T("-"));
                }
                else if (sFileContent.Find("MD5 : ") >= 0)
                {
                    appCrucial(_T("MD5 Found and NOT good %s \n"), sFileName.c_str());
                }
            }

            appGetLattice()->m_pBosonField[0]->InitialFieldWithFile(sFileName, eLoadType);

            switch (eJob)
            {
                case ETMJ_Polyakov:
                    {
                        pPL->OnConfigurationAccepted(_FIELDS, NULL);
                    }
                    break;
                case ETMJ_BosonValue:
                    {
                        pBV->OnConfigurationAccepted(_FIELDS, NULL);
                    }
                    break;
                case ETMJ_Plaquette:
                    {
                        pWilsonPath->OnConfigurationAccepted(_FIELDS, NULL);
                    }
                    break;
                default:
                    break;
            }

            if (uiNewLine > 0 && ((iEndN - uiN + 1) % uiNewLine == 0))
            {
                appPushLogDate(TRUE);
                appGeneral(_T("\n="));
                appPopLogDate();
            }
            else
            {
                appPushLogDate(FALSE);
                appGeneral(_T("="));
                appPopLogDate();
            }
            
        }
        appGeneral(_T("\n*)\n"));

#pragma endregion

        switch (eJob)
        {
            case ETMJ_Polyakov:
                {
                    pPL->Export(sCSVSavePrefix, iStartN, iEndN, PrefixList[uiOmega], uiOmega, iListStart);
                }
                break;
            case ETMJ_BosonValue:
                {
                    CCString sFileNameWrite1;
                    CCString sFileNameWrite2;
                    sFileNameWrite1.Format(_T("%s_%s_bosonvalueC.csv"), sCSVSavePrefix.c_str(), PrefixList[uiOmega].c_str());
                    sFileNameWrite2.Format(_T("%s_%s_bosonvalueC_ZSlice.csv"), sCSVSavePrefix.c_str(), PrefixList[uiOmega].c_str());

                    WriteComplexArray(sFileNameWrite1, pBV->m_lstEveryConfigurationC);
                    WriteComplexArray2(sFileNameWrite2, pBV->m_lstEveryConfigurationZsliceC);

                    sFileNameWrite1.Format(_T("%s_%s_bosonvalueR.csv"), sCSVSavePrefix.c_str(), PrefixList[uiOmega].c_str());
                    sFileNameWrite2.Format(_T("%s_%s_bosonvalueR_ZSlice.csv"), sCSVSavePrefix.c_str(), PrefixList[uiOmega].c_str());

                    WriteRealArray(sFileNameWrite1, pBV->m_lstEveryConfigurationR);
                    WriteRealArray2(sFileNameWrite2, pBV->m_lstEveryConfigurationZsliceR);
                }
                break;
            case ETMJ_Plaquette:
                {
                    CCString sCSVFile;
                    sCSVFile.Format(_T("%s_%s_wilsonloops.csv"), sCSVSavePrefix.c_str(), PrefixList[uiOmega].c_str());
                    WriteComplexArray2(sCSVFile, pWilsonPath->m_lstV);
                }
                break;
            default:
                break;
        }

        appGeneral(_T("\n"));
    }

    appGeneral(_T("\n(*"));
    appPopLogDate();

    appGeneral(_T("\n=====================================\n========= finished! ==========\n*)"));

    appQuitCLG();

    return 0;
}


