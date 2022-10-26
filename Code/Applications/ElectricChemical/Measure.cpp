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
        pF1Light = dynamic_cast<CFieldFermionKSSU3GammaEM*>(appGetLattice()->GetPooledFieldById(2));
        pF2Light = dynamic_cast<CFieldFermionKSSU3GammaEM*>(appGetLattice()->GetPooledFieldById(2));
        pF1Heavy = dynamic_cast<CFieldFermionKSSU3GammaEM*>(appGetLattice()->GetPooledFieldById(3));
        pF2Heavy = dynamic_cast<CFieldFermionKSSU3GammaEM*>(appGetLattice()->GetPooledFieldById(3));
    }

    appSetLogDate(FALSE);
    appGeneral(_T("(* ==== Start Measure ========= *)\n"));
    pPL->Reset();
    pCCLight->Reset();
    pCCHeavy->Reset();

    pCCLight->SetFieldCount(iFieldCount);
    pCCHeavy->SetFieldCount(iFieldCount);

#pragma region Measure

    for (INT uiOmega = iListStart; uiOmega < lstElectric.Num() && uiOmega < iListEnd; ++uiOmega)
    {
        appGeneral(_T("\n========= Electric =%f Chemical = %f  ==========\n"), lstElectric[uiOmega], lstChemical[uiOmega]);

        pU->m_fCoeffGamma54 = lstChemical[uiOmega];
        pD->UpdatePooledParamters();
        pD->m_fCoeffGamma54 = lstChemical[uiOmega];
        pD->UpdatePooledParamters();
        pU1->InitialU1Real(EURT_None, EURT_E_t, EURT_None, F(0.0), lstElectric[uiOmega], F(0.0));

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
                sFileName.Format(_T("%s\\EC%d\\%sEC_%d_%d.con"), sSubFolderName.c_str(), uiOmega, sSavePrefix.c_str(), uiOmega, uiN);
                sTxtFileName.Format(_T("%s\\EC%d\\%sEC_%d_%d.txt"), sSubFolderName.c_str(), uiOmega, sSavePrefix.c_str(), uiOmega, uiN);
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

            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, eLoadType);

            switch (eJob)
            {
            case EGMJ_Polyakov:
            {
                pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
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
                    pF1Light->InverseD(appGetLattice()->m_pGaugeField);
                    pF1Light->FixBoundary();
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
                        appGetLattice()->m_pGaugeField,
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
            sFileNameWrite1.Format(_T("%s_%d_polyakov.csv"), sCSVSavePrefix.c_str(), uiOmega);
            sFileNameWrite2.Format(_T("%s_%d_polyakov_ZSlice.csv"), sCSVSavePrefix.c_str(), uiOmega);

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

            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma1);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma2);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma3);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma4);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma5);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma51);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma52);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma53);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma54);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSSigma12);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSSigma13);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSSigma14);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSSigma23);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSSigma24);
            _CLG_EXPORT_CHIRAL(pCCLight, CMTKSSigma34);

            _CLG_EXPORT_CHIRAL(pCCHeavy, ChiralKS);
            if (pCCHeavy->m_bMeasureConnect)
            {
                _CLG_EXPORT_CHIRAL(pCCHeavy, ConnectSusp);
            }

            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma1);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma2);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma3);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma4);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma5);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma51);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma52);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma53);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma54);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSSigma12);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSSigma13);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSSigma14);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSSigma23);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSSigma24);
            _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSSigma34);

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
        pF1Heavy->Return();
        pF2Heavy->Return();
    }

    appQuitCLG();

    return 0;
}


