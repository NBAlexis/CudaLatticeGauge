//=============================================================================
// FILENAME : MeasureU1.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [10/05/2021 nbale]
//=============================================================================

#include "StaggeredRotation.h"

__DEFINE_ENUM(EDistributionJobKSU1,
    EDJKSU1_Polyakov,
    EDJKSU1_Chiral,
    )


INT MeasurementU1(CParameters& params)
{

#pragma region read parameters

    appSetupLog(params);

    INT iVaule = 0;
    params.FetchValueINT(_T("StartOmega"), iVaule);
    UINT iStartOmega = static_cast<UINT>(iVaule);

    iVaule = 10;
    params.FetchValueINT(_T("EndOmega"), iVaule);
    UINT iEndOmega = static_cast<UINT>(iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("StartN"), iVaule);
    UINT iStartN = static_cast<UINT>(iVaule);

    //iVaule = 1;
    //params.FetchValueINT(_T("FermionMomentum"), iVaule);
    //UBOOL bJF = 0 != iVaule;

    iVaule = 200;
    params.FetchValueINT(_T("EndN"), iVaule);
    UINT iEndN = static_cast<UINT>(iVaule);

    iVaule = 10;
    params.FetchValueINT(_T("StochasticFieldCount"), iVaule);
    UINT iFieldCount = static_cast<UINT>(iVaule);

    //iVaule = 0;
    //params.FetchValueINT(_T("AlsoCheckMD5"), iVaule);
    //UBOOL bCheckMd5 = (0 != iVaule);

    //iVaule = 0;
    //params.FetchValueINT(_T("MeasureCCS"), iVaule);
    //UBOOL bMeasureCCS = (0 != iVaule);

    //iVaule = 1;
    //params.FetchValueINT(_T("CheckGaugeFixing"), iVaule);
    //UBOOL bCheckGaugeFixing = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("UseZ4"), iVaule);
    UBOOL bZ4 = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("SubFolder"), iVaule);
    UBOOL bSubFolder = 0 != iVaule;

    CCString sValue = _T("EDJKSU1_Chiral");
    params.FetchStringValue(_T("EDistributionJobKSU1"), sValue);
    EDistributionJobKSU1 eJob = __STRING_TO_ENUM(EDistributionJobKSU1, sValue);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    CCString sCSVSavePrefix;
    params.FetchStringValue(_T("CSVSavePrefix"), sCSVSavePrefix);
    appGeneral(_T("csv save prefix: %s\n"), sCSVSavePrefix.c_str());

    CCString sSubFolderPrefix;
    params.FetchStringValue(_T("SubFolderPrefix"), sSubFolderPrefix);
    appGeneral(_T("sub folder prefix: %s\n"), sSubFolderPrefix.c_str());

    Real fBeta = F(0.0);
    params.FetchValueReal(_T("GaugeBate"), fBeta);

    Real fOmega = F(0.1);
    params.FetchValueReal(_T("OmegaRange"), fOmega);
    fOmega = fOmega / iEndOmega;

    iVaule = 0;
    params.FetchValueINT(_T("SaveFermionFile"), iVaule);
    UBOOL bSaveFermion = 0 != iVaule;

    iVaule = 1;
    params.FetchValueINT(_T("FermionFileIndexStart"), iVaule);
    UINT uiSaveFermionStart = static_cast<UINT>(iVaule);

    CCString sFermionHead;
    params.FetchStringValue(_T("FermionFileHead"), sFermionHead);
    appGeneral(_T("FermionFileHead: %s\n"), sFermionHead.c_str());

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

    UINT uiMaxL = (_HC_Lx + 1) / 2 - 1;
    uiMaxL = uiMaxL * uiMaxL;
    TArray<TArray<CCString>> r_omega_idx;
    for (UINT i = 0; i < iEndOmega * iEndOmega * uiMaxL; ++i)
    {
        TArray<CCString> newlst;
        r_omega_idx.AddItem(newlst);
    }
    TArray<Real> lstR;
    TArray<TArray<CLGComplex>> lstPolyIn;
    TArray<TArray<CLGComplex>> lstPolyOut;
    TArray<TArray<Real>> lstPolyInAbs;
    TArray<TArray<Real>> lstPolyOutAbs;
    //TArray<TArray<CLGComplex>> lstPolyInZ;
    //TArray<TArray<CLGComplex>> lstPolyOutZ;

    //============================================
    // Attention here heavy light is wrong
    // "heavy" is light, "light" is heavy
    //============================================
    
    //CCommonData::m_fBeta = fBeta;
    UINT uiNewLine = (iEndN - iStartN + 1) / 5;
    CMeasurePolyakovU1XY* pPL = dynamic_cast<CMeasurePolyakovU1XY*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
    CMeasureChiralCondensateKS* pCCLight = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureChiralCondensateKS* pCCHeavy = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    //CMeasureAngularMomentumKS* pFALight = dynamic_cast<CMeasureAngularMomentumKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(4));
    //CMeasureAngularMomentumKS* pFAHeavy = dynamic_cast<CMeasureAngularMomentumKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(5));
    //CMeasureAMomentumJG* pJG = dynamic_cast<CMeasureAMomentumJG*>(appGetLattice()->m_pMeasurements->GetMeasureById(6));
    //CMeasureConnectedSusceptibilityKS* pCCSLight = dynamic_cast<CMeasureConnectedSusceptibilityKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(7));
    //CMeasureConnectedSusceptibilityKS* pCCSHeavy = dynamic_cast<CMeasureConnectedSusceptibilityKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(8));
    //CMeasureWilsonLoopXY* pWilson = dynamic_cast<CMeasureWilsonLoopXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(9));

    //CMeasureAction* pPE = dynamic_cast<CMeasureAction*>(appGetLattice()->m_pMeasurements->GetMeasureById(6));
    //CActionFermionWilsonNf2* pAF = dynamic_cast<CActionFermionWilsonNf2*>(appGetLattice()->m_pActionList[1]);

    CActionGaugePlaquetteRotatingU1* pAG = dynamic_cast<CActionGaugePlaquetteRotatingU1*>(appGetLattice()->m_pActionList.Num() > 0 ? appGetLattice()->m_pActionList[0] : NULL);
    CFieldFermionKSU1R* pLight = dynamic_cast<CFieldFermionKSU1R*>(appGetLattice()->GetFieldById(2));
    CFieldFermionKSU1R* pHeavy = dynamic_cast<CFieldFermionKSU1R*>(appGetLattice()->GetFieldById(3));

    CFieldFermionKSU1* pF1Light = NULL;
    CFieldFermionKSU1* pF2Light = NULL;
    CFieldFermionKSU1* pF1Heavy = NULL;
    CFieldFermionKSU1* pF2Heavy = NULL;

    if (EDJKSU1_Chiral == eJob)
    {
        pF1Light = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(2));
        pF2Light = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(2));
        pF1Heavy = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(3));
        pF2Heavy = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(3));
    }

    appPushLogDate(FALSE);
    //CFieldGaugeSU3* pStaple = NULL;
    //if (EDJKS_VR == eJob)
    //{
    //    pStaple = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField->GetCopy());
    //}

    for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
    {
        //CCommonData::m_fOmega = fOmega * uiOmega;
        if (NULL != pAG)
        {
            pAG->SetGaugeOmega(fOmega * uiOmega);
        }
        if (NULL != pLight)
        {
            pLight->SetFermionOmega(fOmega * uiOmega);
        }
        if (NULL != pHeavy)
        {
            pHeavy->SetFermionOmega(fOmega * uiOmega);
        }
        appGeneral(_T("(* ==== Omega(%f) ========= *)\n"), fOmega * uiOmega);
        pPL->Reset();
        //pJG->Reset();
        pCCLight->Reset();
        pCCHeavy->Reset();
        //pFALight->Reset();
        //pFAHeavy->Reset();
        //pCCSLight->Reset();
        //pCCSHeavy->Reset();
        //pWilson->Reset();

        pCCLight->SetFieldCount(iFieldCount);
        pCCHeavy->SetFieldCount(iFieldCount);
        //pFALight->SetFieldCount(iFieldCount);
        //pFAHeavy->SetFieldCount(iFieldCount);

#pragma region Measure

        appGeneral(_T("(*"));
        for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
        {
            CCString sFileName;
            CCString sTxtFileName;
            if (bSubFolder)
            {
                sFileName.Format(_T("%s/O%d/%sR_Nt%d_O%d_%d.con"), sSubFolderPrefix.c_str(), uiOmega, sSavePrefix.c_str(), _HC_Lt, uiOmega, uiN);
                sTxtFileName.Format(_T("%s/O%d/%sR_Nt%d_O%d_%d.txt"), sSubFolderPrefix.c_str(), uiOmega, sSavePrefix.c_str(), _HC_Lt, uiOmega, uiN);
            }
            else
            {
                sFileName.Format(_T("%sR_Nt%d_O%d_%d.con"), sSavePrefix.c_str(), _HC_Lt, uiOmega, uiN);
                sTxtFileName.Format(_T("%sR_Nt%d_O%d_%d.txt"), sSavePrefix.c_str(), _HC_Lt, uiOmega, uiN);
            }
            //appGeneral(_T("checking %s ..."), sFileName);
            //if (EDJKS_CheckMD5 == eJob || bCheckMd5)
            //{
            //    UINT uiSize = 0;
            //    BYTE* fileContent = appGetFileSystem()->ReadAllBytes(sFileName, uiSize);
            //    CCString sMD5 = "MD5 : " + CLGMD5Hash(fileContent, uiSize);
            //    CCString sMD5old = "MD5 : " + CLGMD5Hash_OLD(fileContent, uiSize);
            //    CCString sFileContent = appGetFileSystem()->ReadAllText(sTxtFileName);
            //    if (sFileContent.Find(sMD5) >= 0)
            //    {
            //        if (EDJKS_CheckMD5 == eJob)
            //        {
            //            appGeneral(_T("MD5 Found and good %s \n"), sFileName.c_str());
            //        }
            //        else
            //        {
            //            appGeneral(_T("-"));
            //        }
            //    }
            //    else if (sFileContent.Find(sMD5old) >= 0)
            //    {
            //        if (EDJKS_CheckMD5 == eJob)
            //        {
            //            appGeneral(_T("MD5 Found and good but use old bad MD5 %s \n"), sFileName.c_str());
            //        }
            //        else
            //        {
            //            appGeneral(_T("-"));
            //        }
            //        sFileContent = sFileContent.Replace(sMD5old, sMD5);
            //        appGetFileSystem()->WriteAllText(sTxtFileName, sFileContent);
            //    }
            //    else if (sFileContent.Find("MD5 : ") >= 0)
            //    {
            //        appCrucial(_T("MD5 Found and NOT good %s \n"), sFileName.c_str());
            //    }
            //    else
            //    {
            //        if (EDJKS_CheckMD5 == eJob)
            //        {
            //            appGeneral(_T("MD5 Not Found so add to it %s \n"), sFileName.c_str());
            //        }
            //        else
            //        {
            //            appGeneral(_T("+"));
            //        }
            //        sFileContent = sFileContent + "\n" + sMD5 + "\n";
            //        appGetFileSystem()->WriteAllText(sTxtFileName, sFileContent);
            //    }

            //    if (EDJKS_CheckMD5 == eJob)
            //    {
            //        continue;
            //    }
            //}
            
            appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName, eLoadType);

            switch (eJob)
            {
                case EDJKSU1_Polyakov:
                {
                    pPL->OnConfigurationAccepted(_FIELDS, NULL);
                }
                break;
                case EDJKSU1_Chiral:
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
                            sFermionFile.Format(_T("%s_Light_Nt%d_O%d_%d_F%d"), sFermionHead.c_str(), _HC_Lt, uiOmega, uiN, uiSaveFermionStart + i);
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

                        //pFALight->OnConfigurationAcceptedZ4(
                        //    appGetLattice()->m_pGaugeField,
                        //    NULL,
                        //    pF2Light,
                        //    pF1Light,
                        //    0 == i,
                        //    iFieldCount == i + 1);

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
                        pF1Heavy->InverseD(_FIELDS);
                        pF1Heavy->FixBoundary();
                        if (bSaveFermion)
                        {
                            CCString sFermionFile = "";
                            sFermionFile.Format(_T("%s_Heavy_Nt%d_O%d_%d_F%d"), sFermionHead.c_str(), _HC_Lt, uiOmega, uiN, uiSaveFermionStart + i);
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

                        //pFAHeavy->OnConfigurationAcceptedZ4(
                        //    appGetLattice()->m_pGaugeField,
                        //    NULL,
                        //    pF2Heavy,
                        //    pF1Heavy,
                        //    0 == i,
                        //    iFieldCount == i + 1);
                    }

                    //if (bMeasureCCS)
                    //{
                    //    pCCSLight->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                    //    pCCSHeavy->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                    //}
                }
                break;
                //case EDJKS_AngularMomentum:
                //{
                //    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField);
                //    pJG->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                //}
                //break;
                //case EDJKS_ChiralAndFermionMomentum:
                //{
                //    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField);
                //    pJG->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                //    for (UINT i = 0; i < iFieldCount; ++i)
                //    {
                //        if (bZ4)
                //        {
                //            pF1Light->InitialField(EFIT_RandomZ4);
                //        }
                //        else
                //        {
                //            pF1Light->InitialField(EFIT_RandomGaussian);
                //        }
                //        pF1Light->FixBoundary();
                //        pF1Light->CopyTo(pF2Light);
                //        pF1Light->InverseD(appGetLattice()->m_pGaugeField);
                //        pF1Light->FixBoundary();
                //        if (bSaveFermion)
                //        {
                //            CCString sFermionFile = "";
                //            sFermionFile.Format(_T("%s_Light_Nt%d_O%d_%d_F%d"), sFermionHead.c_str(), _HC_Lt, uiOmega, uiN, uiSaveFermionStart + i);
                //            CCString sMD51 = pF1Light->SaveToFile(sFermionFile + _T("_F1.con"));
                //            CCString sMD52 = pF2Light->SaveToFile(sFermionFile + _T("_F2.con"));
                //            CCString sFileContent = "";
                //            sFileContent = _T("Stochastic Fermion File for ") + sFileName;
                //            if (bZ4)
                //            {
                //                sFileContent = sFileContent + _T("\nZ4\n");
                //            }
                //            else
                //            {
                //                sFileContent = sFileContent + _T("\nGaussian\n");
                //            }
                //            sFileContent = sFileContent + _T("MD51: ") + sMD51 + _T("\n");
                //            sFileContent = sFileContent + _T("MD52: ") + sMD52 + _T("\n");
                //            appGetFileSystem()->WriteAllText(sFermionFile + _T(".txt"), sFileContent);
                //        }

                //        pCCLight->OnConfigurationAcceptedZ4(
                //            appGetLattice()->m_pGaugeField,
                //            NULL,
                //            pF2Light,
                //            pF1Light,
                //            0 == i,
                //            iFieldCount == i + 1);

                //        pFALight->OnConfigurationAcceptedZ4(
                //            appGetLattice()->m_pGaugeField,
                //            NULL,
                //            pF2Light,
                //            pF1Light,
                //            0 == i,
                //            iFieldCount == i + 1);

                //        if (bZ4)
                //        {
                //            pF1Heavy->InitialField(EFIT_RandomZ4);
                //        }
                //        else
                //        {
                //            pF1Heavy->InitialField(EFIT_RandomGaussian);
                //        }
                //        pF1Heavy->FixBoundary();
                //        pF1Heavy->CopyTo(pF2Heavy);
                //        pF1Heavy->InverseD(appGetLattice()->m_pGaugeField);
                //        pF1Heavy->FixBoundary();
                //        if (bSaveFermion)
                //        {
                //            CCString sFermionFile = "";
                //            sFermionFile.Format(_T("%s_Heavy_Nt%d_O%d_%d_F%d"), sFermionHead.c_str(), _HC_Lt, uiOmega, uiN, uiSaveFermionStart + i);
                //            CCString sMD51 = pF1Heavy->SaveToFile(sFermionFile + _T("_F1.con"));
                //            CCString sMD52 = pF1Heavy->SaveToFile(sFermionFile + _T("_F2.con"));
                //            CCString sFileContent = "";
                //            sFileContent = _T("Stochastic Fermion File for ") + sFileName;
                //            if (bZ4)
                //            {
                //                sFileContent = sFileContent + _T("\nZ4\n");
                //            }
                //            else
                //            {
                //                sFileContent = sFileContent + _T("\nGaussian\n");
                //            }
                //            sFileContent = sFileContent + _T("MD51: ") + sMD51 + _T("\n");
                //            sFileContent = sFileContent + _T("MD52: ") + sMD52 + _T("\n");
                //            appGetFileSystem()->WriteAllText(sFermionFile + _T(".txt"), sFileContent);
                //        }

                //        pCCHeavy->OnConfigurationAcceptedZ4(
                //            appGetLattice()->m_pGaugeField,
                //            NULL,
                //            pF2Heavy,
                //            pF1Heavy,
                //            0 == i,
                //            iFieldCount == i + 1);

                //        pFAHeavy->OnConfigurationAcceptedZ4(
                //            appGetLattice()->m_pGaugeField,
                //            NULL,
                //            pF2Heavy,
                //            pF1Heavy,
                //            0 == i,
                //            iFieldCount == i + 1);
                //    }

                //    if (bMeasureCCS)
                //    {
                //        pCCSLight->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                //        pCCSHeavy->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                //    }
                //}
                //break;
                //case EDJKS_PlaqutteEnergy:
                //{
                //    //pPE->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                //}
                //break;
                //case EDJKS_VR:
                //    {
                //        appGetLattice()->m_pGaugeField->CalculateOnlyStaple(pStaple);
                //        appGetLattice()->m_pGaugeSmearing->GaugeSmearing(appGetLattice()->m_pGaugeField, pStaple);
                //        pWilson->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                //        if (uiN == iStartN)
                //        {
                //            TArray<Real> lstRadius;
                //            for (INT i = 0; i < pWilson->m_lstR.Num(); ++i)
                //            {
                //                lstRadius.AddItem(_hostsqrt(static_cast<Real>(pWilson->m_lstR[i])));
                //            }
                //            CCString sRadiousFile;
                //            sRadiousFile.Format(_T("%s_VR_R.csv"), sCSVSavePrefix.c_str());
                //            WriteStringFileRealArray(sRadiousFile, lstRadius);
                //        }
                //    }
                //break;
                default:
                    break;
            }

            if ((iEndN - uiN + 1) % uiNewLine == 0)
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

        //if (EDJKS_CheckMD5 == eJob)
        //{
        //    continue;
        //}

#pragma endregion

        switch (eJob)
        {
            case EDJKSU1_Polyakov:
            {
                CCString sFileNameWrite1;
                CCString sFileNameWrite2;
                CCString sFileNameWrite3;
                sFileNameWrite1.Format(_T("%s_polyakov_Nt%d_R.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
                sFileNameWrite2.Format(_T("%s_polyakov_Nt%d_O%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiOmega);
                sFileNameWrite3.Format(_T("%s_polyakovabs_Nt%d_O%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiOmega);

                //extract result
                assert(static_cast<INT>(iEndN - iStartN + 1)* pPL->m_lstR.Num() == pPL->m_lstP.Num());

                if (uiOmega == iStartOmega)
                {
                    for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                    {
                        lstR.AddItem(F(0.5) * _hostsqrt(static_cast<Real>(pPL->m_lstR[i])));
                    }
                    WriteStringFileRealArray(sFileNameWrite1, lstR);
                }

                TArray<TArray<CLGComplex>> polyakovOmgR;
                TArray<CLGComplex> polyIn;
                TArray<CLGComplex> polyOut;

                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    TArray<CLGComplex> thisConfiguration;
                    for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                    {
                        thisConfiguration.AddItem(pPL->m_lstP[j * pPL->m_lstR.Num() + i]);
                    }
                    polyakovOmgR.AddItem(thisConfiguration);
                    polyIn.AddItem(pPL->m_lstLoopInner[j]);
                    polyOut.AddItem(pPL->m_lstLoop[j]);
                }
                lstPolyIn.AddItem(polyIn);
                lstPolyOut.AddItem(polyOut);
                WriteStringFileComplexArray2(sFileNameWrite2, polyakovOmgR);

                TArray<TArray<Real>> polyakovOmgRAbs;
                TArray<Real> polyInAbs;
                TArray<Real> polyOutAbs;

                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    TArray<Real> thisConfiguration;
                    for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                    {
                        thisConfiguration.AddItem(pPL->m_lstPAbs[j * pPL->m_lstR.Num() + i]);
                    }
                    polyakovOmgRAbs.AddItem(thisConfiguration);
                    polyInAbs.AddItem(pPL->m_lstLoopAbsInner[j]);
                    polyOutAbs.AddItem(pPL->m_lstLoopAbs[j]);
                }
                lstPolyInAbs.AddItem(polyInAbs);
                lstPolyOutAbs.AddItem(polyOutAbs);
                WriteStringFileRealArray2(sFileNameWrite3, polyakovOmgRAbs);
            }
            break;
            case EDJKSU1_Chiral:
            {
                _CLG_EXPORT_CHIRAL(pCCLight, ChiralKS, uiOmega, O);
                _CLG_EXPORT_CHIRAL(pCCLight, ConnectSusp, uiOmega, O);
                _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma3, uiOmega, O);
                _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma4, uiOmega, O);
                _CLG_EXPORT_CHIRAL(pCCHeavy, ChiralKS, uiOmega, O);
                _CLG_EXPORT_CHIRAL(pCCHeavy, ConnectSusp, uiOmega, O);
                _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma3, uiOmega, O);
                _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma4, uiOmega, O);

                if (uiOmega == iStartOmega)
                {
                    TArray<Real> lstRadius;
                    for (INT i = 0; i < pCCLight->m_lstR.Num(); ++i)
                    {
                        lstRadius.AddItem(F(0.5)* _hostsqrt(static_cast<Real>(pCCLight->m_lstR[i])));
                    }
                    CCString sRadiousFile;
                    sRadiousFile.Format(_T("%s_condensateR.csv"), sCSVSavePrefix.c_str());
                    WriteStringFileRealArray(sRadiousFile, lstRadius);
                }

                //if (bMeasureCCS)
                //{
                //    CCString sFileNameWriteCCS;
                //    sFileNameWriteCCS.Format(_T("%s_lightCCS_Nt%d_O%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiOmega);
                //    WriteStringFileComplexArray(sFileNameWriteCCS, pCCSLight->m_lstResults);
                //    sFileNameWriteCCS.Format(_T("%s_heavyCCS_Nt%d_O%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiOmega);
                //    WriteStringFileComplexArray(sFileNameWriteCCS, pCCSHeavy->m_lstResults);
                //}
            }
            break;
            //case EDJKS_AngularMomentum:
            //{
            //    _CLG_EXPORT_ANGULAR(pJG, JG, uiOmega, O);
            //    _CLG_EXPORT_ANGULAR(pJG, JGS, uiOmega, O);
            //    _CLG_EXPORT_ANGULAR(pJG, JGChen, uiOmega, O);
            //    _CLG_EXPORT_ANGULAR(pJG, JGSurf, uiOmega, O);
            //    _CLG_EXPORT_ANGULAR(pJG, JGPot, uiOmega, O);
            //}
            //break;
            //case EDJKS_ChiralAndFermionMomentum:
            //{
            //    _CLG_EXPORT_ANGULAR(pJG, JG, uiOmega, O);
            //    _CLG_EXPORT_ANGULAR(pJG, JGS, uiOmega, O);
            //    _CLG_EXPORT_ANGULAR(pJG, JGChen, uiOmega, O);
            //    _CLG_EXPORT_ANGULAR(pJG, JGSurf, uiOmega, O);
            //    _CLG_EXPORT_ANGULAR(pJG, JGPot, uiOmega, O);

            //    _CLG_EXPORT_CHIRAL(pCCLight, ChiralKS, uiOmega, O);
            //    _CLG_EXPORT_CHIRAL(pCCLight, ConnectSusp, uiOmega, O);
            //    _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma3, uiOmega, O);
            //    _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma4, uiOmega, O);
            //    _CLG_EXPORT_CHIRAL(pCCHeavy, ChiralKS, uiOmega, O);
            //    _CLG_EXPORT_CHIRAL(pCCHeavy, ConnectSusp, uiOmega, O);
            //    _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma3, uiOmega, O);
            //    _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma4, uiOmega, O);

            //    _CLG_EXPORT_CHIRAL(pFALight, OrbitalKS, uiOmega, O);
            //    _CLG_EXPORT_CHIRAL(pFALight, SpinKS, uiOmega, O);
            //    _CLG_EXPORT_CHIRAL(pFALight, PotentialKS, uiOmega, O);
            //    _CLG_EXPORT_CHIRAL(pFAHeavy, OrbitalKS, uiOmega, O);
            //    _CLG_EXPORT_CHIRAL(pFAHeavy, SpinKS, uiOmega, O);
            //    _CLG_EXPORT_CHIRAL(pFAHeavy, PotentialKS, uiOmega, O);

            //    if (uiOmega == iStartOmega)
            //    {
            //        TArray<Real> lstRadius;
            //        for (INT i = 0; i < pCCLight->m_lstR.Num(); ++i)
            //        {
            //            lstRadius.AddItem(F(0.5) * _hostsqrt(static_cast<Real>(pCCLight->m_lstR[i])));
            //        }
            //        CCString sRadiousFile;
            //        sRadiousFile.Format(_T("%s_condensateR.csv"), sCSVSavePrefix.c_str());
            //        WriteStringFileRealArray(sRadiousFile, lstRadius);
            //    }

            //    if (bMeasureCCS)
            //    {
            //        CCString sFileNameWriteCCS;
            //        sFileNameWriteCCS.Format(_T("%s_lightCCS_Nt%d_O%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiOmega);
            //        WriteStringFileComplexArray(sFileNameWriteCCS, pCCSLight->m_lstResults);
            //        sFileNameWriteCCS.Format(_T("%s_heavyCCS_Nt%d_O%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiOmega);
            //        WriteStringFileComplexArray(sFileNameWriteCCS, pCCSHeavy->m_lstResults);
            //    }
            //}
            //break;
            //case EDJKS_PlaqutteEnergy:
            //{
            //    //CCString sFileName;
            //    //sFileName.Format(_T("%s_plaqutte.csv"), sCSVSavePrefix.c_str());
            //    //WriteStringFileRealArray(sFileName, pPE->m_lstData);
            //}
            //break;
            //case EDJKS_VR:
            //    {
            //        CCString sCSVFile;
            //        sCSVFile.Format(_T("%s_VR_Nt%d_O%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiOmega);
            //        TArray<TArray<CLGComplex>> vrs;
            //        for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
            //        {
            //            TArray<CLGComplex> thisConfiguration;
            //            for (INT i = 0; i < pWilson->m_lstR.Num(); ++i)
            //            {
            //                for (UINT t = 0; t < (_HC_Lt - 1); ++t)
            //                {
            //                    thisConfiguration.AddItem(pWilson->m_lstC[j][i][t]);
            //                }
            //            }
            //            vrs.AddItem(thisConfiguration);
            //        }
            //        WriteStringFileComplexArray2(sCSVFile, vrs);
            //    }
            //break;
            default:
                break;
        }

        appGeneral(_T("\n"));
    }

    //if (NULL != pStaple)
    //{
    //    appSafeDelete(pStaple);
    //}

    //if (EDJKS_CheckMD5 == eJob)
    //{
    //    if (NULL != pF1Light)
    //    {
    //        pF1Light->Return();
    //        pF2Light->Return();
    //        pF1Heavy->Return();
    //        pF2Heavy->Return();
    //    }

    //    appQuitCLG();
    //    return 0;
    //}

    switch (eJob)
    {
        case EDJKSU1_Polyakov:
        {
            CCString sFileNameWrite1;
            CCString sFileNameWrite2;
            sFileNameWrite1.Format(_T("%s_polyakov_Nt%d_In.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            sFileNameWrite2.Format(_T("%s_polyakov_Nt%d_Out.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            WriteStringFileComplexArray2(sFileNameWrite1, lstPolyIn);
            WriteStringFileComplexArray2(sFileNameWrite2, lstPolyOut);

            sFileNameWrite1.Format(_T("%s_polyakovabs_Nt%d_In.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            sFileNameWrite2.Format(_T("%s_polyakovabs_Nt%d_Out.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            WriteStringFileRealArray2(sFileNameWrite1, lstPolyInAbs);
            WriteStringFileRealArray2(sFileNameWrite2, lstPolyOutAbs);

        }
        break;
        default:
            break;
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


