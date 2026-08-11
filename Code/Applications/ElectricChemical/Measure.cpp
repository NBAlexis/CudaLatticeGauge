//=============================================================================
// FILENAME : Measure.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [09/29/2022 nbale]
//=============================================================================

#include "ElectricChemical.h"

__DEFINE_ENUM(EGradientMeasureJob,
    EGMJ_Polyakov,
    EGMJ_Chiral,
    EGMJ_Meson,
    EGMJ_BerryPhase,
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

    iVaule = 0;
    params.FetchValueINT(_T("ListStart"), iVaule);
    const INT iListStart = iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("ListEnd"), iVaule);
    const INT iListEnd = iVaule;

    TArray<Real> lstElectric;
    params.FetchValueArrayReal(_T("Electric"), lstElectric);

    TArray<Real> lstChemical;
    params.FetchValueArrayReal(_T("Chemical"), lstChemical);

    TArray<Real> lstMagnetic;
    params.FetchValueArrayReal(_T("Magnetic"), lstMagnetic);

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
    CMeasureChiralCondensateKS* pCCHeavy = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
    CMeasureMesonCorrelatorStaggeredSimple2* pMeson = dynamic_cast<CMeasureMesonCorrelatorStaggeredSimple2*>(appGetLattice()->m_pMeasurements->GetMeasureById(4));
    CMeasureBerryPhase* pBerryPhaseU = dynamic_cast<CMeasureBerryPhase*>(appGetLattice()->m_pMeasurements->GetMeasureById(5));
    CMeasureBerryPhase* pBerryPhaseD = dynamic_cast<CMeasureBerryPhase*>(appGetLattice()->m_pMeasurements->GetMeasureById(6));

    CFieldFermionKSSU3GammaEM* pU = NULL;
    CFieldFermionKSSU3GammaEM* pD = NULL;

    CFieldFermionKSSU3GammaEM* pF1Light = NULL;
    CFieldFermionKSSU3GammaEM* pF2Light = NULL;
    CFieldFermionKSSU3GammaEM* pF1Heavy = NULL;
    CFieldFermionKSSU3GammaEM* pF2Heavy = NULL;

    pU = dynamic_cast<CFieldFermionKSSU3GammaEM*>(appGetLattice()->GetFieldById(2));
    pD = dynamic_cast<CFieldFermionKSSU3GammaEM*>(appGetLattice()->GetFieldById(3));
    CFieldGaugeU1Real* pU1 = dynamic_cast<CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(4));

    if (EGMJ_Chiral == eJob)
    {
        pF1Light = dynamic_cast<CFieldFermionKSSU3GammaEM*>(appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));
        pF2Light = dynamic_cast<CFieldFermionKSSU3GammaEM*>(appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));
        pF1Heavy = dynamic_cast<CFieldFermionKSSU3GammaEM*>(appGetLattice()->GetPooledFieldById(3, _T(__FILE__), __LINE__));
        pF2Heavy = dynamic_cast<CFieldFermionKSSU3GammaEM*>(appGetLattice()->GetPooledFieldById(3, _T(__FILE__), __LINE__));
    }

    appPushLogDate(FALSE);
    appGeneral(_T("(* ==== Start Measure ========= *)\n"));
    pPL->Reset();
    pCCLight->Reset();
    pCCHeavy->Reset();

    pCCLight->SetFieldCount(iFieldCount);
    pCCHeavy->SetFieldCount(iFieldCount);

#pragma region Measure

    for (INT uiOmega = iListStart; uiOmega < lstElectric.Num() && uiOmega < iListEnd; ++uiOmega)
    {
        appGeneral(_T("\n========= Electric =%f Chemical = %f Magnetic = %f ==========\n"), lstElectric[uiOmega], lstChemical[uiOmega], lstMagnetic[uiOmega]);

        pU->m_fCoeffGamma54 = lstChemical[uiOmega];
        //pU->UpdatePooledParamters();
        pD->m_fCoeffGamma54 = lstChemical[uiOmega];
        //pD->UpdatePooledParamters();
        EU1RealType eEzType = pU1->m_eE;
        EU1RealType eBzType = pU1->m_eB;
        pU1->InitialU1Real(EURT_None, eEzType, eBzType, F(0.0), lstElectric[uiOmega], lstMagnetic[uiOmega], FALSE);

        pPL->Reset();
        pCCLight->Reset();
        pCCLight->SetFieldCount(iFieldCount);
        pCCHeavy->Reset();
        pCCHeavy->SetFieldCount(iFieldCount);

        pMeson->Reset();

        pBerryPhaseU->Reset();
        pBerryPhaseD->Reset();

        appGeneral(_T("(*"));
        for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
        {
            CCString sFileName;
            CCString sTxtFileName;
            if (bSubFolder)
            {
                sFileName.Format(_T("%s/EC%d/%sEC_%d_%d.con"), sSubFolderName.c_str(), uiOmega, sSavePrefix.c_str(), uiOmega, uiN);
                sTxtFileName.Format(_T("%s/EC%d/%sEC_%d_%d.txt"), sSubFolderName.c_str(), uiOmega, sSavePrefix.c_str(), uiOmega, uiN);
            }
            else
            {
                sFileName.Format(_T("%sEC_%d_%d.con"), sSavePrefix.c_str(), uiOmega, uiN);
                sTxtFileName.Format(_T("%sEC_%d_%d.txt"), sSavePrefix.c_str(), uiOmega, uiN);
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
                    pF1Light->FixBoundary(EFB_Field);
                    pF1Light->CopyTo(pF2Light);
                    pF1Light->InverseD(_FIELDS);
                    pF1Light->FixBoundary(EFB_Field);
                    if (bSaveFermion)
                    {
                        CCString sFermionFile = "";
                        sFermionFile.Format(_T("%s_Light_%d_%d_F%d"), sFermionHead.c_str(), uiOmega, uiN, uiSaveFermionStart + i);
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


                    if (bZ4)
                    {
                        pF1Heavy->InitialField(EFIT_RandomZ4);
                    }
                    else
                    {
                        pF1Heavy->InitialField(EFIT_RandomGaussian);
                    }
                    pF1Heavy->FixBoundary(EFB_Field);
                    pF1Heavy->CopyTo(pF2Heavy);
                    pF1Heavy->InverseD(_FIELDS);
                    pF1Heavy->FixBoundary(EFB_Field);
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
                        _FIELDS,
                        NULL,
                        pF2Heavy,
                        pF1Heavy,
                        0 == i,
                        iFieldCount == i + 1);

                }
            }
            break;
            case EGMJ_Meson:
                {
                    pMeson->OnConfigurationAccepted(_FIELDS, NULL);
                }
                break;
            case EGMJ_BerryPhase:
                {
                    pBerryPhaseU->OnConfigurationAccepted(_FIELDS, NULL);
                    pBerryPhaseD->OnConfigurationAccepted(_FIELDS, NULL);
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
            pPL->Export(sCSVSavePrefix, iStartN, iEndN, uiOmega, iListStart);
        }
        break;
        case EGMJ_Chiral:
        {
            _CLG_EXPORT_CHIRAL(pCCLight, ChiralKS, uiOmega);
            if (pCCLight->m_bMeasureConnect)
            {
                _CLG_EXPORT_CHIRAL(pCCLight, ConnectSusp, uiOmega);
            }

            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma1, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma2, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma3, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma4, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma5, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma51, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma52, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma53, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma54, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSSigma12, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSSigma13, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSSigma14, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSSigma23, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSSigma24, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSSigma34, uiOmega);

            _CLG_EXPORT_CHIRAL(pCCHeavy, ChiralKS, uiOmega);
            if (pCCHeavy->m_bMeasureConnect)
            {
                _CLG_EXPORT_CHIRAL(pCCHeavy, ConnectSusp, uiOmega);
            }

            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma1, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma2, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma3, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma4, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma5, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma51, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma52, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma53, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma54, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSSigma12, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSSigma13, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSSigma14, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSSigma23, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSSigma24, uiOmega);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSSigma34, uiOmega);

        }
        break;
        case EGMJ_Meson:
        {
            static const TCHAR* heads1[4] =
            {
                "uu",
                "ud",
                "du",
                "dd"
            };

            for (INT i = 0; i < 16; ++i)
            {
                for (INT j = 0; j < 4; ++j)
                {
                    TArray<TArray<DOUBLE>> onemesonconfig;
                    TArray<TArray<DOUBLE>> onemesonconfigX;
                    TArray<TArray<DOUBLE>> onemesonconfigY;
                    TArray<TArray<DOUBLE>> onemesonconfigZ;

                    for (UINT k = 0; k < (iEndN - iStartN + 1); ++k)
                    {
                        TArray<DOUBLE> onemeson_oneconfig;
                        TArray<DOUBLE> onemeson_oneconfigX;
                        TArray<DOUBLE> onemeson_oneconfigY;
                        TArray<DOUBLE> onemeson_oneconfigZ;

                        for (INT uiT = 0; uiT < _HC_Lti; ++uiT)
                        {
                            onemeson_oneconfig.AddItem(pMeson->m_lstResults[k][i * 4 + j][uiT]);
                        }
                        for (INT uiX = 0; uiX < _HC_Lxi; ++uiX)
                        {
                            onemeson_oneconfigX.AddItem(pMeson->m_lstResultsX[k][i * 4 + j][uiX]);
                        }
                        for (INT uiY = 0; uiY < _HC_Lyi; ++uiY)
                        {
                            onemeson_oneconfigY.AddItem(pMeson->m_lstResultsY[k][i * 4 + j][uiY]);
                        }
                        for (INT uiZ = 0; uiZ < _HC_Lzi; ++uiZ)
                        {
                            onemeson_oneconfigZ.AddItem(pMeson->m_lstResultsZ[k][i * 4 + j][uiZ]);
                        }
                        onemesonconfig.AddItem(onemeson_oneconfig);
                        onemesonconfigX.AddItem(onemeson_oneconfigX);
                        onemesonconfigY.AddItem(onemeson_oneconfigY);
                        onemesonconfigZ.AddItem(onemeson_oneconfigZ);
                    }

                    CCString sFileNameMeson;
                    CCString sFileNameMesonX;
                    CCString sFileNameMesonY;
                    CCString sFileNameMesonZ;
                    sFileNameMeson.Format(_T("%sT_%s%d_%d.csv"), sCSVSavePrefix.c_str(), heads1[j], i, uiOmega);
                    sFileNameMesonX.Format(_T("%sX_%s%d_%d.csv"), sCSVSavePrefix.c_str(), heads1[j], i, uiOmega);
                    sFileNameMesonY.Format(_T("%sY_%s%d_%d.csv"), sCSVSavePrefix.c_str(), heads1[j], i, uiOmega);
                    sFileNameMesonZ.Format(_T("%sZ_%s%d_%d.csv"), sCSVSavePrefix.c_str(), heads1[j], i, uiOmega);
                    WriteRealArray2(sFileNameMeson, onemesonconfig);
                    WriteRealArray2(sFileNameMesonX, onemesonconfigX);
                    WriteRealArray2(sFileNameMesonY, onemesonconfigY);
                    WriteRealArray2(sFileNameMesonZ, onemesonconfigZ);
                }
            }
        }
        break;
        case EGMJ_BerryPhase:
        {
            CCString sFileNameMeson;
            sFileNameMeson.Format(_T("%s_berryphaseU_%d.csv"), sCSVSavePrefix.c_str(), uiOmega);
            WriteRealArray2(sFileNameMeson, pBerryPhaseU->m_lstData);
            sFileNameMeson.Format(_T("%s_berryphaseUXY_%d.csv"), sCSVSavePrefix.c_str(), uiOmega);
            WriteRealArray2(sFileNameMeson, pBerryPhaseU->m_lstDataXY);
            sFileNameMeson.Format(_T("%s_berryphaseUXZ_%d.csv"), sCSVSavePrefix.c_str(), uiOmega);
            WriteRealArray2(sFileNameMeson, pBerryPhaseU->m_lstDataXZ);
            sFileNameMeson.Format(_T("%s_berryphaseUYZ_%d.csv"), sCSVSavePrefix.c_str(), uiOmega);
            WriteRealArray2(sFileNameMeson, pBerryPhaseU->m_lstDataYZ);

            sFileNameMeson.Format(_T("%s_berryphaseD_%d.csv"), sCSVSavePrefix.c_str(), uiOmega);
            WriteRealArray2(sFileNameMeson, pBerryPhaseD->m_lstData);
            sFileNameMeson.Format(_T("%s_berryphaseDXY_%d.csv"), sCSVSavePrefix.c_str(), uiOmega);
            WriteRealArray2(sFileNameMeson, pBerryPhaseD->m_lstDataXY);
            sFileNameMeson.Format(_T("%s_berryphaseDXZ_%d.csv"), sCSVSavePrefix.c_str(), uiOmega);
            WriteRealArray2(sFileNameMeson, pBerryPhaseD->m_lstDataXZ);
            sFileNameMeson.Format(_T("%s_berryphaseDYZ_%d.csv"), sCSVSavePrefix.c_str(), uiOmega);
            WriteRealArray2(sFileNameMeson, pBerryPhaseD->m_lstDataYZ);

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
        pF1Heavy->Return();
        pF2Heavy->Return();
    }

    appQuitCLG();

    return 0;
}


