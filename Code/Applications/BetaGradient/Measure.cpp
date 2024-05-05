//=============================================================================
// FILENAME : Measure.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [08/17/2022 nbale]
//=============================================================================

#include "BetaGradient.h"

__DEFINE_ENUM(EGradientMeasureJob,
    EGMJ_Polyakov,
    EGMJ_Chiral,
    EGMJ_DoubleToFloat,
    )

INT Measurement(CParameters& params)
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

    CCString sValue = _T("EGMJ_Polyakov");
    params.FetchStringValue(_T("DistributionJob"), sValue);
    EGradientMeasureJob eJob = __STRING_TO_ENUM(EGradientMeasureJob, sValue);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    CCString sCSVSavePrefix;
    params.FetchStringValue(_T("CSVSavePrefix"), sCSVSavePrefix);
    appGeneral(_T("csv save prefix: %s\n"), sCSVSavePrefix.c_str());

    TArray<CCString> PrefixList;
    params.FetchStringVectorValue(_T("PrefixList"), PrefixList);

    iVaule = 0;
    params.FetchValueINT(_T("BetaStride"), iVaule);
    const INT iBetaStride = iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("ListStart"), iVaule);
    const INT iListStart = iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("ListEnd"), iVaule);
    const INT iListEnd = iVaule;

    TArray<DOUBLE> fMiddleBeta;
    params.FetchValueArrayDOUBLE(_T("MiddleBetaList"), fMiddleBeta);
    if (fMiddleBeta.Num() != PrefixList.Num() || PrefixList.Num() < 1)
    {
        appCrucial(_T("sSavePrefix and fMiddleBeta not corrected!\n"));
        return 0;
    }

    DOUBLE fDeltaBeta = 3.0;
    params.FetchValueDOUBLE(_T("DeltaBeta"), fDeltaBeta);

    iVaule = 0;
    params.FetchValueINT(_T("SubFolder"), iVaule);
    UBOOL bSubFolder = 0 != iVaule;

    CCString sSubFolderName;
    params.FetchStringValue(_T("SubFolderPrefix"), sSubFolderName);

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
    CMeasureChiralCondensateKS* pCCLight = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    //CMeasureChiralCondensateKS* pCCHeavy = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));

    CActionGaugePlaquetteGradient* pGaugeGradient = NULL;
    if (appGetLattice()->m_pActionList.Num() > 0)
    {
        pGaugeGradient = dynamic_cast<CActionGaugePlaquetteGradient*>(appGetLattice()->GetActionById(1));
        TArray<DOUBLE> betaArray;
        if (0 == iBetaStride)
        {
            for (INT i = 0; i < _HC_Lzi; ++i)
            {
                betaArray.AddItem(sin((-1.0 + 2.0 * i / static_cast<DOUBLE>(_HC_Lzi)) * PI) * fDeltaBeta + fMiddleBeta[0]);
            }
        }
        else
        {
            if (15 != _HC_Lzi)
            {
                appCrucial(_T("Lz must be 11"));
            }
            DOUBLE upper = fMiddleBeta[0] + fDeltaBeta;
            DOUBLE lower = fMiddleBeta[0] - fDeltaBeta;

            if (1 == iBetaStride)
            {
                // ----- . ++++++
                for (INT i = 0; i < 7; ++i)
                {
                    betaArray.AddItem(lower);
                }

                betaArray.AddItem(fMiddleBeta[0]);

                for (INT i = 0; i < 7; ++i)
                {
                    betaArray.AddItem(upper);
                }
            }
            else if (2 == iBetaStride)
            {
                // ---- ... +++++
                for (INT i = 0; i < 6; ++i)
                {
                    betaArray.AddItem(lower);
                }

                for (INT i = -1; i <= 1; ++i)
                {
                    betaArray.AddItem(fMiddleBeta[0] + i * 0.5 * fDeltaBeta);
                }

                for (INT i = 0; i < 6; ++i)
                {
                    betaArray.AddItem(upper);
                }
            }
            else if (3 == iBetaStride)
            {
                // --- ..... +++
                for (INT i = 0; i < 5; ++i)
                {
                    betaArray.AddItem(lower);
                }

                for (INT i = -2; i <= 2; ++i)
                {
                    betaArray.AddItem(fMiddleBeta[0] + i * 0.333333333333333 * fDeltaBeta);
                }

                for (INT i = 0; i < 5; ++i)
                {
                    betaArray.AddItem(upper);
                }
            }
            else if (4 == iBetaStride)
            {
                // --- ..... +++
                for (INT i = 0; i < 4; ++i)
                {
                    betaArray.AddItem(lower);
                }

                for (INT i = -3; i <= 3; ++i)
                {
                    betaArray.AddItem(fMiddleBeta[0] + i * 0.25 * fDeltaBeta);
                }

                for (INT i = 0; i < 4; ++i)
                {
                    betaArray.AddItem(upper);
                }
            }
        }
        pGaugeGradient->SetBeta(betaArray);
    }

    CFieldFermionKSSU3* pF1Light = NULL;
    CFieldFermionKSSU3* pF2Light = NULL;
    //CFieldFermionKSSU3* pF1Heavy = NULL;
    //CFieldFermionKSSU3* pF2Heavy = NULL;

    if (EGMJ_Chiral == eJob)
    {
        pF1Light = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2));
        pF2Light = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2));
        //pF1Heavy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(3));
        //pF2Heavy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(3));
    }

    appPushLogDate(FALSE);
    appGeneral(_T("(* ==== Start Measure ========= *)\n"));
    pPL->Reset();
    pCCLight->Reset();
    //pCCHeavy->Reset();

    pCCLight->SetFieldCount(iFieldCount);
    //pCCHeavy->SetFieldCount(iFieldCount);

#pragma region Measure

    for (INT uiOmega = iListStart; uiOmega < PrefixList.Num() && uiOmega < iListEnd; ++uiOmega)
    {
        if (NULL != pGaugeGradient)
        {
            TArray<DOUBLE> betaArray;
            if (0 == iBetaStride)
            {
                for (INT i = 0; i < _HC_Lzi; ++i)
                {
                    betaArray.AddItem(sin((-1.0 + 2.0 * i / static_cast<DOUBLE>(_HC_Lzi)) * PI) * fDeltaBeta + fMiddleBeta[uiOmega]);
                }
            }
            else
            {
                if (15 != _HC_Lzi)
                {
                    appCrucial(_T("Lz must be 11"));
                }
                DOUBLE upper = fMiddleBeta[uiOmega] + fDeltaBeta;
                DOUBLE lower = fMiddleBeta[uiOmega] - fDeltaBeta;

                if (1 == iBetaStride)
                {
                    // ----- . ++++++
                    for (INT i = 0; i < 7; ++i)
                    {
                        betaArray.AddItem(lower);
                    }

                    betaArray.AddItem(fMiddleBeta[uiOmega]);

                    for (INT i = 0; i < 7; ++i)
                    {
                        betaArray.AddItem(upper);
                    }
                }
                else if (2 == iBetaStride)
                {
                    // ---- ... +++++
                    for (INT i = 0; i < 6; ++i)
                    {
                        betaArray.AddItem(lower);
                    }

                    for (INT i = -1; i <= 1; ++i)
                    {
                        betaArray.AddItem(fMiddleBeta[uiOmega] + i * 0.5 * fDeltaBeta);
                    }

                    for (INT i = 0; i < 6; ++i)
                    {
                        betaArray.AddItem(upper);
                    }
                }
                else if (3 == iBetaStride)
                {
                    // --- ..... +++
                    for (INT i = 0; i < 5; ++i)
                    {
                        betaArray.AddItem(lower);
                    }

                    for (INT i = -2; i <= 2; ++i)
                    {
                        betaArray.AddItem(fMiddleBeta[uiOmega] + i * 0.333333333333333 * fDeltaBeta);
                    }

                    for (INT i = 0; i < 5; ++i)
                    {
                        betaArray.AddItem(upper);
                    }
                }
                else if (4 == iBetaStride)
                {
                    // --- ..... +++
                    for (INT i = 0; i < 4; ++i)
                    {
                        betaArray.AddItem(lower);
                    }

                    for (INT i = -3; i <= 3; ++i)
                    {
                        betaArray.AddItem(fMiddleBeta[uiOmega] + i * 0.25 * fDeltaBeta);
                    }

                    for (INT i = 0; i < 4; ++i)
                    {
                        betaArray.AddItem(upper);
                    }
                }
            }
            pGaugeGradient->SetBeta(betaArray);
        }
        appGeneral(_T("(* ==== Beta(%s) ========= *)\n"), PrefixList[uiOmega].c_str());
        pPL->Reset();
        pCCLight->Reset();
        pCCLight->SetFieldCount(iFieldCount);

        appGeneral(_T("(*"));
        for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
        {
            CCString sFileName;
            CCString sTxtFileName;
            if (bSubFolder)
            {
                sFileName.Format(_T("%s/%s/%sGradient_%s_%d.con"), sSubFolderName.c_str(), PrefixList[uiOmega].c_str(), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
                sTxtFileName.Format(_T("%s/%s/%sGradient_%s_%d.txt"), sSubFolderName.c_str(), PrefixList[uiOmega].c_str(), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
            }
            else
            {
                sFileName.Format(_T("%sGradient_%s_%d.con"), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
                sTxtFileName.Format(_T("%sGradient_%s_%d.txt"), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
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

            appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName, eLoadType);

            switch (eJob)
            {
            case EGMJ_Polyakov:
            {
                pPL->OnConfigurationAccepted(_FIELDS, NULL);
            }
            break;
            case EGMJ_Chiral:
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
                    pF1Light->InverseD(_FIELDS);
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
                        _FIELDS,
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
                        sFermionFile.Format(_T("%s_Heavy_%d_F%d"), sFermionHead.c_str(), uiN, uiSaveFermionStart + i);
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
                appPushLogDate(TRUE);
                appGeneral(_T("\n="));
                appPopLogDate();
            }
            else
            {
                appGeneral(_T("="));
            }

        }
        appGeneral(_T("\n*)\n"));

#pragma endregion

        switch (eJob)
        {
        case EGMJ_Polyakov:
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
            WriteStringFileComplexArray2(sFileNameWrite2, polyakovOmgZSlice);
        }
        break;
        case EGMJ_Chiral:
        {
            _CLG_EXPORT_CHIRAL(pCCLight, ChiralKS);
            if (pCCLight->m_bMeasureConnect)
            {
                _CLG_EXPORT_CHIRAL(pCCLight, ConnectSusp);
            }

            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma3);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma4);

            /*
            _CLG_EXPORT_CHIRAL(pCCHeavy, ChiralKS);
            if (pCCHeavy->m_bMeasureConnect)
            {
                _CLG_EXPORT_CHIRAL(pCCHeavy, ConnectSusp);
            }

            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma3);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma4);
            */
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
    if (NULL != pF1Light)
    {
        pF1Light->Return();
        pF2Light->Return();
        //pF1Heavy->Return();
        //pF2Heavy->Return();
    }

    appQuitCLG();

    return 0;
}


