//=============================================================================
// FILENAME : MeasureBetaScan.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [08/17/2022 nbale]
//=============================================================================

#include "BetaGradient.h"

__DEFINE_ENUM(EBetaScanMeasureJob,
    EBSMJ_Polyakov,
    EBSMJ_Chiral,
    EBSMJ_Wilson,
    EBSMJ_DoubleToFloat,
    )

INT MeasurementBetaScan(CParameters& params)
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

    iVaule = 10;
    params.FetchValueINT(_T("StochasticFieldCount"), iVaule);
    UINT iFieldCount = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("AlsoCheckMD5"), iVaule);
    UBOOL bCheckMd5 = (0 != iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("UseZ4"), iVaule);
    UBOOL bZ4 = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("SubFolder"), iVaule);
    UBOOL bSubFolder = 0 != iVaule;

    CCString sSubFolderName;
    params.FetchStringValue(_T("SubFolderPrefix"), sSubFolderName);

    CCString sValue = _T("EBSMJ_Polyakov");
    params.FetchStringValue(_T("DistributionJob"), sValue);
    EBetaScanMeasureJob eJob = __STRING_TO_ENUM(EBetaScanMeasureJob, sValue);

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
    CMeasureWilsonLoop* pWL = dynamic_cast<CMeasureWilsonLoop*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
    CFieldGaugeSU3* pStaple = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField->GetCopy());

    CMeasureChiralCondensateKS* pCCLight = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    //CMeasureChiralCondensateKS* pCCHeavy = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));

    CActionGaugePlaquette* pAG = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->m_pActionList.Num() > 0 ? appGetLattice()->m_pActionList[0] : NULL);

    CFieldFermionKSSU3* pF1Light = NULL;
    CFieldFermionKSSU3* pF2Light = NULL;
    //CFieldFermionKSSU3* pF1Heavy = NULL;
    //CFieldFermionKSSU3* pF2Heavy = NULL;

    if (EBSMJ_Chiral == eJob)
    {
        pF1Light = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2));
        pF2Light = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2));
        //pF1Heavy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(3));
        //pF2Heavy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(3));
    }

    appSetLogDate(FALSE);

    for (INT uiOmega = iListStart; uiOmega < BetaList.Num() && uiOmega < iListEnd; ++uiOmega)
    {
        if (NULL != pAG)
        {
            pAG->SetBeta(BetaList[uiOmega]);
        }
        appGeneral(_T("(* ==== Beta(%f) ========= *)\n"), BetaList[uiOmega]);
        pPL->Reset();
        pWL->Reset();
        pCCLight->Reset();
        //pCCHeavy->Reset();

        pCCLight->SetFieldCount(iFieldCount);
        //pCCHeavy->SetFieldCount(iFieldCount);

#pragma region Measure

        appGeneral(_T("(*"));
        for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
        {
            CCString sFileName;
            CCString sTxtFileName;
            if (bSubFolder)
            {
                sFileName.Format(_T("%s/%s/%sBetaScan_%s_%d.con"), sSubFolderName.c_str(), PrefixList[uiOmega].c_str(), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
                sTxtFileName.Format(_T("%s/%s/%sBetaScan_%s_%d.txt"), sSubFolderName.c_str(), PrefixList[uiOmega].c_str(), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
            }
            else
            {
                sFileName.Format(_T("%sBetaScan_%s_%d.con"), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
                sTxtFileName.Format(_T("%sBetaScan_%s_%d.txt"), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
            }


            //appGeneral(_T("checking %s ..."), sFileName);
            if (bCheckMd5)
            {
                UINT uiSize = 0;
                BYTE* fileContent = appGetFileSystem()->ReadAllBytes(sFileName, uiSize);
                CCString sMD5 = "MD5 : " + CLGMD5Hash(fileContent, uiSize);
                CCString sMD5old = "MD5 : " + CLGMD5Hash_OLD(fileContent, uiSize);
                CCString sFileContent = appGetFileSystem()->ReadAllText(sTxtFileName);
                if (sFileContent.Find(sMD5) >= 0)
                {
                    appGeneral(_T("-"));
                }
                else if (sFileContent.Find(sMD5old) >= 0)
                {
                    appGeneral(_T("-"));
                    sFileContent = sFileContent.Replace(sMD5old, sMD5);
                    appGetFileSystem()->WriteAllText(sTxtFileName, sFileContent);
                }
                else if (sFileContent.Find("MD5 : ") >= 0)
                {
                    appCrucial(_T("MD5 Found and NOT good %s \n"), sFileName.c_str());
                }
                else
                {
                    appGeneral(_T("+"));
                    sFileContent = sFileContent + "\n" + sMD5 + "\n";
                    appGetFileSystem()->WriteAllText(sTxtFileName, sFileContent);
                }
            }
            
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, eLoadType);

            switch (eJob)
            {
                case EBSMJ_Polyakov:
                {
                    pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                }
                break;
                case EBSMJ_Wilson:
                {
                    appGetLattice()->m_pGaugeField->CalculateOnlyStaple(pStaple);
                    appGetLattice()->m_pGaugeSmearing->GaugeSmearing(appGetLattice()->m_pGaugeField, pStaple);
                    pWL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                    if (uiN == iStartN)
                    {
                        TArray<Real> lstRadius;
                        for (INT i = 0; i < pWL->m_lstR.Num(); ++i)
                        {
                            lstRadius.AddItem(_hostsqrt(static_cast<Real>(pWL->m_lstR[i])));
                        }
                        CCString sRadiousFile;
                        sRadiousFile.Format(_T("%s_VR_R.csv"), sCSVSavePrefix.c_str());
                        WriteStringFileRealArray(sRadiousFile, lstRadius);
                    }
                }
                break;
                case EBSMJ_Chiral:
                {
                    for (UINT i = 0; i < iFieldCount; ++i)
                    {
                        if (bZ4)
                        {
                            pF1Light->InitialField(EFIT_RandomZ4);
                        }
                        else
                        {
                            pF1Light->InitialField(EFIT_RandomGaussian);
                        }
                        pF1Light->FixBoundary();
                        pF1Light->CopyTo(pF2Light);
                        pF1Light->InverseD(appGetLattice()->m_pGaugeField);
                        pF1Light->FixBoundary();
                        if (bSaveFermion)
                        {
                            CCString sFermionFile = "";
                            sFermionFile.Format(_T("%s_Light_%s_%d_F%d"), sFermionHead.c_str(), PrefixList[uiOmega].c_str(), uiN, uiSaveFermionStart + i);
                            CCString sMD51 = pF1Light->SaveToFile(sFermionFile + _T("_F1.con"));
                            CCString sMD52 = pF2Light->SaveToFile(sFermionFile + _T("_F2.con"));
                            CCString sFileContent = "";
                            sFileContent = _T("Stochastic Fermion File for ") + sFileName;
                            if (bZ4)
                            {
                                sFileContent = sFileContent + _T("\nZ4\n");
                            }
                            else
                            {
                                sFileContent = sFileContent + _T("\nGaussian\n");
                            }
                            sFileContent = sFileContent + _T("MD51: ") + sMD51 + _T("\n");
                            sFileContent = sFileContent + _T("MD52: ") + sMD52 + _T("\n");
                            appGetFileSystem()->WriteAllText(sFermionFile + _T(".txt"), sFileContent);
                        }

                        pCCLight->OnConfigurationAcceptedZ4(
                            appGetLattice()->m_pGaugeField,
                            NULL,
                            pF2Light,
                            pF1Light,
                            0 == i,
                            iFieldCount == i + 1);

                        /*
                        if (bZ4)
                        {
                            pF1Heavy->InitialField(EFIT_RandomZ4);
                        }
                        else
                        {
                            pF1Heavy->InitialField(EFIT_RandomGaussian);
                        }
                        pF1Heavy->FixBoundary();
                        pF1Heavy->CopyTo(pF2Heavy);
                        pF1Heavy->InverseD(appGetLattice()->m_pGaugeField);
                        pF1Heavy->FixBoundary();
                        if (bSaveFermion)
                        {
                            CCString sFermionFile = "";
                            sFermionFile.Format(_T("%s_Heavy_%s_%d_F%d"), sFermionHead.c_str(), PrefixList[uiOmega].c_str(), uiN, uiSaveFermionStart + i);
                            CCString sMD51 = pF1Heavy->SaveToFile(sFermionFile + _T("_F1.con"));
                            CCString sMD52 = pF2Heavy->SaveToFile(sFermionFile + _T("_F2.con"));
                            CCString sFileContent = "";
                            sFileContent = _T("Stochastic Fermion File for ") + sFileName;
                            if (bZ4)
                            {
                                sFileContent = sFileContent + _T("\nZ4\n");
                            }
                            else
                            {
                                sFileContent = sFileContent + _T("\nGaussian\n");
                            }
                            sFileContent = sFileContent + _T("MD51: ") + sMD51 + _T("\n");
                            sFileContent = sFileContent + _T("MD52: ") + sMD52 + _T("\n");
                            appGetFileSystem()->WriteAllText(sFermionFile + _T(".txt"), sFileContent);
                        }

                        pCCHeavy->OnConfigurationAcceptedZ4(
                            appGetLattice()->m_pGaugeField,
                            NULL,
                            pF2Heavy,
                            pF1Heavy,
                            0 == i,
                            iFieldCount == i + 1);
                        */
                    }
                }
                break;
                default:
                    break;
            }

            if (uiNewLine > 0 && ((iEndN - uiN + 1) % uiNewLine == 0))
            {
                appSetLogDate(TRUE);
                appGeneral(_T("\n="));
                appSetLogDate(FALSE);
            }
            else
            {
                appSetLogDate(FALSE);
                appGeneral(_T("="));
                appSetLogDate(TRUE);
            }
            
        }
        appGeneral(_T("\n*)\n"));

#pragma endregion

        switch (eJob)
        {
            case EBSMJ_Polyakov:
            {
                CCString sFileNameWrite1;
                CCString sFileNameWrite2;
                sFileNameWrite1.Format(_T("%s_%s_polyakov.csv"), sCSVSavePrefix.c_str(), PrefixList[uiOmega].c_str());
                sFileNameWrite2.Format(_T("%s_%s_polyakov_ZSlice.csv"), sCSVSavePrefix.c_str(), PrefixList[uiOmega].c_str());

                //extract result
                TArray<CLGComplex> polyOut;
                TArray<TArray<CLGComplex>> polyakovOmgZSlice;
                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    polyOut.AddItem(pPL->m_lstLoop[j]);

                    if (pPL->m_bMeasureZSlice)
                    {
                        TArray<CLGComplex> thisConfigurationZSlice;
                        for (UINT i = 0; i < _HC_Lz; ++i)
                        {
                            thisConfigurationZSlice.AddItem(pPL->m_lstPZSlice[j * _HC_Lz + i]);
                        }
                        polyakovOmgZSlice.AddItem(thisConfigurationZSlice);
                    }
                }
                WriteStringFileComplexArray(sFileNameWrite1, polyOut);
                if (pPL->m_bMeasureZSlice)
                {
                    WriteStringFileComplexArray2(sFileNameWrite2, polyakovOmgZSlice);
                }
            }
            break;
            case EBSMJ_Chiral:
            {
                _CLG_EXPORT_CHIRAL(pCCLight, ChiralKS);
                if (pCCLight->m_bMeasureConnect)
                {
                    _CLG_EXPORT_CHIRAL(pCCLight, ConnectSusp);
                }
                
                _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma3);
                _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma4);

                /*
                _CLG_EXPORT_CHIRAL_SCAN(pCCHeavy, ChiralKS);
                if (pCCHeavy->m_bMeasureConnect)
                {
                    _CLG_EXPORT_CHIRAL_SCAN(pCCHeavy, ConnectSusp);
                }
                
                _CLG_EXPORT_CHIRAL_SCAN(pCCHeavy, CMTKSGamma3);
                _CLG_EXPORT_CHIRAL_SCAN(pCCHeavy, CMTKSGamma4);
                */
            }
            break;
            case EBSMJ_Wilson:
            {
                CCString sCSVFile;
                sCSVFile.Format(_T("%s_VR_Nt%d_O%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiOmega);
                TArray<TArray<CLGComplex>> vrs;
                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    TArray<CLGComplex> thisConfiguration;
                    for (INT i = 0; i < pWL->m_lstR.Num(); ++i)
                    {
                        for (UINT t = 0; t < _HC_Lt / 2; ++t)
                        {
                            thisConfiguration.AddItem(pWL->m_lstC[j][i][t]);
                        }
                    }
                    vrs.AddItem(thisConfiguration);
                }
                WriteStringFileComplexArray2(sCSVFile, vrs);
            }
            break;
            default:
                break;
        }

        appGeneral(_T("\n"));
    }

    appGeneral(_T("\n(*"));
    appSetLogDate(TRUE);

    appGeneral(_T("\n=====================================\n========= finished! ==========\n*)"));
    if (NULL != pF1Light)
    {
        pF1Light->Return();
        pF2Light->Return();
        //pF1Heavy->Return();
        //pF2Heavy->Return();
    }

    appSafeDelete(pStaple);

    appQuitCLG();

    return 0;
}


