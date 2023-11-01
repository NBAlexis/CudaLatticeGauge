//=============================================================================
// FILENAME : Measure.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [09/29/2022 nbale]
//=============================================================================

#include "ElectricChemical.h"

__DEFINE_ENUM(EGradientMeasureRWJob,
    EGMJRW_Polyakov,
    EGMJRW_Chiral,
    EGMJRW_ChiralDiagnal,
    )

INT MeasureRW(CParameters& params)
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

    CCString sValue = _T("EGMJRW_Chiral");
    params.FetchStringValue(_T("DistributionJob"), sValue);
    EGradientMeasureRWJob eJob = __STRING_TO_ENUM(EGradientMeasureRWJob, sValue);

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

    iVaule = 0;
    params.FetchValueINT(_T("FreeFermion"), iVaule);
    UBOOL bFreeFermion = (0 != iVaule);

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

    CFieldFermionKSSU3Gamma* pU = NULL;

    CFieldFermionKSSU3Gamma* pF1Light = NULL;
    CFieldFermionKSSU3Gamma* pF2Light = NULL;

    pU = dynamic_cast<CFieldFermionKSSU3Gamma*>(appGetLattice()->GetFieldById(2));

    if (EGMJRW_Chiral == eJob || EGMJRW_ChiralDiagnal == eJob)
    {
        pF1Light = dynamic_cast<CFieldFermionKSSU3Gamma*>(appGetLattice()->GetPooledFieldById(2));
        pF2Light = dynamic_cast<CFieldFermionKSSU3Gamma*>(appGetLattice()->GetPooledFieldById(2));
    }

    appPushLogDate(FALSE);
    appGeneral(_T("(* ==== Start Measure ========= *)\n"));
    pPL->Reset();
    pCCLight->Reset();

    pCCLight->SetFieldCount(iFieldCount);

#pragma region Measure

    for (INT uiOmega = iListStart; uiOmega < lstChemical.Num() && uiOmega < iListEnd; ++uiOmega)
    {
        appGeneral(_T("\n========= Chemical = %f  ==========\n"), lstChemical[uiOmega]);

        pU->m_fCoeffGamma4 = lstChemical[uiOmega];
        pU->UpdatePooledParamters();

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

            if (bFreeFermion)
            {
                appGetLattice()->m_pGaugeField->InitialField(EFIT_Identity);
            }
            else
            {
                appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, eLoadType);
            }
            

            switch (eJob)
            {
            case EGMJRW_Polyakov:
            {
                pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
            }
            break;
            case EGMJRW_Chiral:
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

                }
            }
            break;
            case EGMJRW_ChiralDiagnal:
                {
                    //appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField);
                    //appGeneral(_T("debug gamma4 : %f, %f\n"), pF1Light->m_fCoeffGamma4, pF2Light->m_fCoeffGamma4);
                    CCString sFileDiagnal;
                    sFileDiagnal.Format(_T("%s_diagnal_Nt%d_IC%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiOmega);
                    TArray<TArray<CLGComplex>> lightdiagnal = pCCLight->ExportDiagnal(appGetLattice()->m_pGaugeField, pF1Light, pF2Light);
                    WriteStringFileComplexArray2(sFileDiagnal, lightdiagnal);
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
        case EGMJRW_Polyakov:
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
        case EGMJRW_Chiral:
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
    }

    appQuitCLG();

    return 0;
}


