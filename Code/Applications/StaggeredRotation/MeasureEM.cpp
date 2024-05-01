//=============================================================================
// FILENAME : MeasureEM.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [10/06/2020 nbale]
//=============================================================================

#include "StaggeredRotation.h"



__DEFINE_ENUM(EDistributionJobKSEM,
    EDJKSEM_Polyakov,
    EDJKSEM_Chiral,
    EDJKSEM_BerryPhase,
    EDJKSEM_Meson,
    )


INT MeasurementEM(CParameters& params)
{

#pragma region read parameters

    appSetupLog(params);

    INT iVaule = 0;
    params.FetchValueINT(_T("EMStart"), iVaule);
    UINT iStartEM = static_cast<UINT>(iVaule);

    iVaule = 10;
    params.FetchValueINT(_T("EMEnd"), iVaule);
    UINT iEMEnd = static_cast<UINT>(iVaule);

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
    params.FetchValueINT(_T("UseZ4"), iVaule);
    UBOOL bZ4 = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("SubFolder"), iVaule);
    UBOOL bSubFolder = 0 != iVaule;

    CCString sValue = _T("EDJKSEM_Polyakov");
    params.FetchStringValue(_T("DistributionJob"), sValue);
    EDistributionJobKSEM eJob = __STRING_TO_ENUM(EDistributionJobKSEM, sValue);

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

    TArray<Real> lstE;
    params.FetchValueArrayReal(_T("ElectricList"), lstE);

    TArray<Real> lstM;
    params.FetchValueArrayReal(_T("MagneticList"), lstM);

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
    const UINT uiLoadFermion = iVaule;

    CCString sLoadFermionFile;
    params.FetchStringValue(_T("LoadFermionFile"), sLoadFermionFile);
    appGeneral(_T("Load Fermion File Name: %s\n"), sLoadFermionFile.c_str());

    CCString sLoadFermionHead;
    params.FetchStringValue(_T("LoadFermionHead"), sLoadFermionHead);
    appGeneral(_T("Load Fermion File Head: %s\n"), sLoadFermionHead.c_str());

    if (0 == lstE.Num() && 0 != lstM.Num())
    {
        for (INT i = 0; i < lstM.Num(); ++i)
        {
            lstE.AddItem(F(0.0));
        }
    }

    if (0 == lstM.Num() && 0 != lstE.Num())
    {
        for (INT i = 0; i < lstE.Num(); ++i)
        {
            lstM.AddItem(F(0.0));
        }
    }

    if ((0 == lstE.Num() && 0 == lstM.Num()) || lstE.Num() != lstM.Num())
    {
        appCrucial(_T("The length of the list ElectricList or MagneticList is wrong!\n"));
        return 1;
    }

    if (iEMEnd < static_cast<UINT>(lstE.Num()))
    {
        iEMEnd = static_cast<UINT>(lstE.Num());
    }

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
    for (UINT i = 0; i < iEMEnd * iEMEnd * uiMaxL; ++i)
    {
        TArray<CCString> newlst;
        r_omega_idx.AddItem(newlst);
    }
    TArray<Real> lstR;
    TArray<TArray<CLGComplex>> lstPolyIn;
    TArray<TArray<CLGComplex>> lstPolyOut;
    TArray<TArray<CLGComplex>> lstPolyInZ;
    TArray<TArray<CLGComplex>> lstPolyOutZ;

    CCommonData::m_fBeta = fBeta;
    UINT uiNewLine = (iEndN - iStartN + 1) / 5;
    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureChiralCondensateKS* pCCLight = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CMeasureChiralCondensateKS* pCCHeavy = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
    CMeasureAngularMomentumKS* pFALight = dynamic_cast<CMeasureAngularMomentumKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(4));
    CMeasureAngularMomentumKS* pFAHeavy = dynamic_cast<CMeasureAngularMomentumKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(5));
    CMeasureAMomentumJG* pJG = dynamic_cast<CMeasureAMomentumJG*>(appGetLattice()->m_pMeasurements->GetMeasureById(6));
    CMeasureBerryPhase* pBPu = dynamic_cast<CMeasureBerryPhase*>(appGetLattice()->m_pMeasurements->GetMeasureById(7));
    CMeasureBerryPhase* pBPd = dynamic_cast<CMeasureBerryPhase*>(appGetLattice()->m_pMeasurements->GetMeasureById(8));
    CMeasureMesonCorrelatorStaggeredSimple2* pMeson = dynamic_cast<CMeasureMesonCorrelatorStaggeredSimple2*>(appGetLattice()->m_pMeasurements->GetMeasureById(9));
    //CMeasureAction* pPE = dynamic_cast<CMeasureAction*>(appGetLattice()->m_pMeasurements->GetMeasureById(6));
    //CActionFermionWilsonNf2* pAF = dynamic_cast<CActionFermionWilsonNf2*>(appGetLattice()->m_pActionList[1]);

    //CActionGaugePlaquetteRotating* pAG = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->m_pActionList.Num() > 0 ? appGetLattice()->m_pActionList[0] : NULL);

    CFieldFermionKSSU3* pF1Light = NULL;
    CFieldFermionKSSU3* pF2Light = NULL;
    CFieldFermionKSSU3* pF1Heavy = NULL;
    CFieldFermionKSSU3* pF2Heavy = NULL;

    if (EDJKSEM_Chiral == eJob)
    {
        pF1Light = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2));
        pF2Light = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2));
        pF1Heavy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(3));
        pF2Heavy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(3));
    }

    appPushLogDate(FALSE);


    for (UINT uiEM = iStartEM; uiEM < iEMEnd; ++uiEM)
    {
        appGeneral(_T("(* ==== Electric(%f) Magnetic(%f) ========= *)\n"), lstE[uiEM], lstM[uiEM]);
        pPL->Reset();
        pJG->Reset();
        pCCLight->Reset();
        pCCHeavy->Reset();
        pFALight->Reset();
        pFAHeavy->Reset();
        pBPu->Reset();
        pBPd->Reset();
        pCCLight->SetFieldCount(iFieldCount);
        pCCHeavy->SetFieldCount(iFieldCount);
        pFALight->SetFieldCount(iFieldCount);
        pFAHeavy->SetFieldCount(iFieldCount);

        pMeson->Reset();

        CCommonData::m_fBz = lstM[uiEM];
        CCommonData::m_fEz = lstE[uiEM];

#pragma region Measure

        appGeneral(_T("(*"));
        for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
        {
            CCString sFileName;
            CCString sTxtFileName;
            if (bSubFolder)
            {
                sFileName.Format(_T("%s/EM%d/%sR_Nt%d_EM%d_%d.con"), sSubFolderPrefix.c_str(), uiEM, sSavePrefix.c_str(), _HC_Lt, uiEM, uiN);
                sTxtFileName.Format(_T("%s/EM%d/%sR_Nt%d_EM%d_%d.txt"), sSubFolderPrefix.c_str(), uiEM, sSavePrefix.c_str(), _HC_Lt, uiEM, uiN);
            }
            else
            {
                sFileName.Format(_T("%sR_Nt%d_EM%d_%d.con"), sSavePrefix.c_str(), _HC_Lt, uiEM, uiN);
                sTxtFileName.Format(_T("%sR_Nt%d_EM%d_%d.txt"), sSavePrefix.c_str(), _HC_Lt, uiEM, uiN);
            }
            //appGeneral(_T("checking %s ..."), sFileName);            
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, eLoadType);
            TArray<CFieldGauge*> gauge;
            gauge.AddItem(appGetLattice()->m_pGaugeField);

            switch (eJob)
            {
                case EDJKSEM_Polyakov:
                {
                    pPL->OnConfigurationAccepted(1, 0, gauge.GetData(), NULL, NULL);
                }
                break;
                case EDJKSEM_Chiral:
                {
                    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField);
                    for (UINT i = 0; i < iFieldCount; ++i)
                    {
                        if (0 == (1 & uiLoadFermion))
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
                                sFermionFile.Format(_T("%s_Light_Nt%d_EM%d_%d_F%d"), sFermionHead.c_str(), _HC_Lt, uiEM, uiN, uiSaveFermionStart + i);
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
                                sFileContent = sFileContent + _T("Beta: ") + appAnyToString(CCommonData::m_fBeta) + _T("\n");
                                sFileContent = sFileContent + _T("Electric: ") + appAnyToString(CCommonData::m_fEz) + _T("\n");
                                sFileContent = sFileContent + _T("Magnetic: ") + appAnyToString(CCommonData::m_fBz) + _T("\n");
                                sFileContent = sFileContent + _T("Mass: ") + appAnyToString(pF1Light->m_f2am) + _T("\n");
                                sFileContent = sFileContent + _T("Chage: ") + appAnyToString(dynamic_cast<CFieldFermionKSSU3EM*>(pF1Light)->GetQ()) + _T("\n");
                                sFileContent = sFileContent + _T("ShiftCenter: ") + (pF1Light->m_bEachSiteEta ? _T("TRUE") : _T("FALSE")) + _T("\n");
                                appGetFileSystem()->WriteAllText(sFermionFile + _T(".txt"), sFileContent);
                            }
                        }
                        else
                        {
                            CCString sF1FileName = "";
                            CCString sF2FileName = "";
                            sF1FileName.Format(_T("%s/EM%d/Light/%s_Light_Nt%d_EM%d_%d_F%d_F1.con"),
                                sLoadFermionHead.c_str(), uiEM, sLoadFermionFile.c_str(), _HC_Lt, uiEM, uiN, i + 1);
                            sF2FileName.Format(_T("%s/EM%d/Light/%s_Light_Nt%d_EM%d_%d_F%d_F2.con"),
                                sLoadFermionHead.c_str(), uiEM, sLoadFermionFile.c_str(), _HC_Lt, uiEM, uiN, i + 1);
                            pF1Light->InitialFieldWithFile(sF1FileName, EFFT_CLGBin);
                            pF2Light->InitialFieldWithFile(sF2FileName, EFFT_CLGBin);
                        }

                        pCCLight->OnConfigurationAcceptedZ4(
                            1, 0, gauge.GetData(), NULL,
                            NULL,
                            pF2Light,
                            pF1Light,
                            0 == i,
                            iFieldCount == i + 1);

                        pFALight->OnConfigurationAcceptedZ4(
                            1, 0, gauge.GetData(), NULL,
                            NULL,
                            pF2Light,
                            pF1Light,
                            0 == i,
                            iFieldCount == i + 1);

                        if (0 == (2 & uiLoadFermion))
                        {
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
                                sFermionFile.Format(_T("%s_Heavy_Nt%d_EM%d_%d_F%d"), sFermionHead.c_str(), _HC_Lt, uiEM, uiN, uiSaveFermionStart + i);
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
                                sFileContent = sFileContent + _T("Beta: ") + appAnyToString(CCommonData::m_fBeta) + _T("\n");
                                sFileContent = sFileContent + _T("Electric: ") + appAnyToString(CCommonData::m_fEz) + _T("\n");
                                sFileContent = sFileContent + _T("Magnetic: ") + appAnyToString(CCommonData::m_fBz) + _T("\n");
                                sFileContent = sFileContent + _T("Mass: ") + appAnyToString(pF1Heavy->m_f2am) + _T("\n");
                                sFileContent = sFileContent + _T("Chage: ") + appAnyToString(dynamic_cast<CFieldFermionKSSU3EM*>(pF1Heavy)->GetQ()) + _T("\n");
                                sFileContent = sFileContent + _T("ShiftCenter: ") + (pF1Heavy->m_bEachSiteEta ? _T("TRUE") : _T("FALSE")) + _T("\n");
                                appGetFileSystem()->WriteAllText(sFermionFile + _T(".txt"), sFileContent);
                            }
                        }
                        else
                        {
                            CCString sF1FileName = "";
                            CCString sF2FileName = "";
                            sF1FileName.Format(_T("%s/EM%d/Heavy/%s_Heavy_Nt%d_EM%d_%d_F%d_F1.con"),
                                sLoadFermionHead.c_str(), uiEM, sLoadFermionFile.c_str(), _HC_Lt, uiEM, uiN, i + 1);
                            sF2FileName.Format(_T("%s/EM%d/Heavy/%s_Heavy_Nt%d_EM%d_%d_F%d_F2.con"),
                                sLoadFermionHead.c_str(), uiEM, sLoadFermionFile.c_str(), _HC_Lt, uiEM, uiN, i + 1);
                            pF1Heavy->InitialFieldWithFile(sF1FileName, EFFT_CLGBin);
                            pF2Heavy->InitialFieldWithFile(sF2FileName, EFFT_CLGBin);
                        }

                        pCCHeavy->OnConfigurationAcceptedZ4(
                            1, 0, gauge.GetData(), NULL,
                            NULL,
                            pF2Heavy,
                            pF1Heavy,
                            0 == i,
                            iFieldCount == i + 1);

                        pFAHeavy->OnConfigurationAcceptedZ4(
                            1, 0, gauge.GetData(), NULL,
                            NULL,
                            pF2Heavy,
                            pF1Heavy,
                            0 == i,
                            iFieldCount == i + 1);
                    }

                }
                break;
                case EDJKSEM_BerryPhase:
                {
                    pBPu->OnConfigurationAccepted(1, 0, gauge.GetData(), NULL, NULL);
                    pBPd->OnConfigurationAccepted(1, 0, gauge.GetData(), NULL, NULL);
                }
                break;
                case EDJKSEM_Meson:
                {
                    pMeson->OnConfigurationAccepted(1, 0, gauge.GetData(), NULL, NULL);
                }
                break;
#if NotYet
                case EDJKS_AngularMomentum:
                {
                    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField);
                    pJG->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                }
                break;
                case EDJKS_ChiralAndFermionMomentum:
                {
                    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField);
                    pJG->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
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

                        pCCLight->OnConfigurationAcceptedZ4(
                            appGetLattice()->m_pGaugeField,
                            NULL,
                            pF2Light,
                            pF1Light,
                            0 == i,
                            iFieldCount == i + 1);

                        pFALight->OnConfigurationAcceptedZ4(
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

                        pCCHeavy->OnConfigurationAcceptedZ4(
                            appGetLattice()->m_pGaugeField,
                            NULL,
                            pF2Heavy,
                            pF1Heavy,
                            0 == i,
                            iFieldCount == i + 1);

                        pFAHeavy->OnConfigurationAcceptedZ4(
                            appGetLattice()->m_pGaugeField,
                            NULL,
                            pF2Heavy,
                            pF1Heavy,
                            0 == i,
                            iFieldCount == i + 1);
                    }
                }
                break;
                case EDJKS_PlaqutteEnergy:
                {
                    //pPE->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
                }
                break;
#endif
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

#pragma endregion

        switch (eJob)
        {
            case EDJKSEM_Polyakov:
            {
                CCString sFileNameWrite1;
                CCString sFileNameWrite2;
                sFileNameWrite1.Format(_T("%s_polyakov_Nt%d_R.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
                sFileNameWrite2.Format(_T("%s_polyakov_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                
                //extract result
                assert(static_cast<INT>(iEndN - iStartN + 1) * pPL->m_lstR.Num() == pPL->m_lstP.Num());
                
                if (uiEM == iStartEM)
                {
                    for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                    {
                        lstR.AddItem(F(0.5)* _hostsqrt(static_cast<Real>(pPL->m_lstR[i])));
                    }
                    WriteStringFileRealArray(sFileNameWrite1, lstR);
                }

                TArray<TArray<CLGComplex>> polyakovOmgR;
                TArray<TArray<CLGComplex>> polyakovOmgZSlice;
                TArray<CLGComplex> polyIn;
                TArray<CLGComplex> polyOut;

                for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                {
                    TArray<CLGComplex> thisConfiguration;
                    for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                    {
                        thisConfiguration.AddItem(pPL->m_lstP[j * pPL->m_lstR.Num() + i]);
                    }
                    if (pPL->m_bMeasureZSlice)
                    {
                        TArray<CLGComplex> thisConfigurationZSlice;
                        for (UINT i = 0; i < _HC_Lz; ++i)
                        {
                            thisConfigurationZSlice.AddItem(pPL->m_lstPZSlice[j * _HC_Lz + i]);
                        }
                        polyakovOmgZSlice.AddItem(thisConfigurationZSlice);
                    }
                    polyakovOmgR.AddItem(thisConfiguration);
                    polyIn.AddItem(pPL->m_lstLoopInner[j]);
                    polyOut.AddItem(pPL->m_lstLoop[j]);
                }
                lstPolyIn.AddItem(polyIn);
                lstPolyOut.AddItem(polyOut);
                WriteStringFileComplexArray2(sFileNameWrite2, polyakovOmgR);
                if (pPL->m_bMeasureZSlice)
                {
                    CCString sFileNameWrite3;
                    sFileNameWrite3.Format(_T("%s_polyakovZSlice_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                    WriteStringFileComplexArray2(sFileNameWrite3, polyakovOmgZSlice);
                }

                if (pPL->m_bMeasureLoopZ)
                {
                    CCString sFileNameWrite3;
                    sFileNameWrite3.Format(_T("%s_polyakovZ_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                    polyakovOmgR.RemoveAll();
                    polyIn.RemoveAll();
                    polyOut.RemoveAll();

                    for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                    {
                        TArray<CLGComplex> thisConfiguration;
                        for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                        {
                            thisConfiguration.AddItem(pPL->m_lstPZ[j * pPL->m_lstR.Num() + i]);
                        }
                        polyakovOmgR.AddItem(thisConfiguration);
                        polyIn.AddItem(pPL->m_lstLoopZInner[j]);
                        polyOut.AddItem(pPL->m_lstLoopZ[j]);
                    }
                    lstPolyInZ.AddItem(polyIn);
                    lstPolyOutZ.AddItem(polyOut);
                    WriteStringFileComplexArray2(sFileNameWrite3, polyakovOmgR);
                }
            }
            break;
            case EDJKSEM_Chiral:
            {
                _CLG_EXPORT_CHIRAL(pCCLight, ChiralKS, uiEM, EM);
                if (pCCLight->m_bMeasureConnect)
                {
                    _CLG_EXPORT_CHIRAL(pCCLight, ConnectSusp, uiEM, EM);
                }
                
                _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma3, uiEM, EM);
                _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma4, uiEM, EM);
                _CLG_EXPORT_CHIRAL(pCCHeavy, ChiralKS, uiEM, EM);
                if (pCCHeavy->m_bMeasureConnect)
                {
                    _CLG_EXPORT_CHIRAL(pCCHeavy, ConnectSusp, uiEM, EM);
                }
                
                _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma3, uiEM, EM);
                _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma4, uiEM, EM);

                if (uiEM == iStartEM)
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

                if (pCCLight->m_bMeasureZSlice)
                {
                    _CLG_EXPORT_CHIRALZSLICE(pCCLight, ChiralKS, uiEM, EM);
                    _CLG_EXPORT_CHIRALZSLICE(pCCLight, CMTKSGamma3, uiEM, EM);
                    _CLG_EXPORT_CHIRALZSLICE(pCCLight, CMTKSGamma4, uiEM, EM);
                }

                if (pCCHeavy->m_bMeasureZSlice)
                {
                    _CLG_EXPORT_CHIRALZSLICE(pCCHeavy, ChiralKS, uiEM, EM);
                    _CLG_EXPORT_CHIRALZSLICE(pCCHeavy, CMTKSGamma3, uiEM, EM);
                    _CLG_EXPORT_CHIRALZSLICE(pCCHeavy, CMTKSGamma4, uiEM, EM);
                }
            }
            break;
            case EDJKSEM_BerryPhase:
            {
                CCString sFileNameWriteBP;
                sFileNameWriteBP.Format(_T("%s_BerryPhaseU_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPu->m_lstData);
                sFileNameWriteBP.Format(_T("%s_BerryPhaseD_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPd->m_lstData);

                sFileNameWriteBP.Format(_T("%s_BerryPhaseUXY_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPu->m_lstDataXY);
                sFileNameWriteBP.Format(_T("%s_BerryPhaseDXY_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPd->m_lstDataXY);

                sFileNameWriteBP.Format(_T("%s_BerryPhaseUXZ_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPu->m_lstDataXZ);
                sFileNameWriteBP.Format(_T("%s_BerryPhaseDXZ_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPd->m_lstDataXZ);

                sFileNameWriteBP.Format(_T("%s_BerryPhaseUXT_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPu->m_lstDataXT);
                sFileNameWriteBP.Format(_T("%s_BerryPhaseDXT_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPd->m_lstDataXT);

                sFileNameWriteBP.Format(_T("%s_BerryPhaseUYZ_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPu->m_lstDataYZ);
                sFileNameWriteBP.Format(_T("%s_BerryPhaseDYZ_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPd->m_lstDataYZ);

                sFileNameWriteBP.Format(_T("%s_BerryPhaseUYT_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPu->m_lstDataYT);
                sFileNameWriteBP.Format(_T("%s_BerryPhaseDYT_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPd->m_lstDataYT);

                sFileNameWriteBP.Format(_T("%s_BerryPhaseUZT_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPu->m_lstDataZT);
                sFileNameWriteBP.Format(_T("%s_BerryPhaseDZT_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), _HC_Lt, uiEM);
                WriteStringFileRealArray2(sFileNameWriteBP, pBPd->m_lstDataZT);
            }
            break;
            case EDJKSEM_Meson:
                {
                    static const TCHAR* heads[16] =
                    {
                        "PSuu",
                        "PSud",
                        "PSdu",
                        "PSdd",
                        "VTuu",
                        "VTud",
                        "VTdu",
                        "VTdd",
                        "PVuu",
                        "PVud",
                        "PVdu",
                        "PVdd",
                        "Suu",
                        "Sud",
                        "Sdu",
                        "Sdd"
                    };
                    for (INT i = 0; i < 16; ++i)
                    {
                        TArray<TArray<DOUBLE>> onemesonconfig;
                        for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
                        {
                            TArray<DOUBLE> onemeson_oneconfig;
                            for (INT uiT = 0; uiT < _HC_Lti - 1; ++uiT)
                            {
                                onemeson_oneconfig.AddItem(pMeson->m_lstResults[j][i][uiT]);
                            }
                            onemesonconfig.AddItem(onemeson_oneconfig);
                        }

                        CCString sFileNameMeson;
                        sFileNameMeson.Format(_T("%s_%s_Nt%d_EM%d.csv"), sCSVSavePrefix.c_str(), heads[i], _HC_Lt, uiEM);
                        WriteStringFileRealArray2(sFileNameMeson, onemesonconfig);
                    }
                }
                break;
#if NotYet
            case EDJKS_AngularMomentum:
            {
                _CLG_EXPORT_ANGULAR(pJG, JG);
                _CLG_EXPORT_ANGULAR(pJG, JGS);
                _CLG_EXPORT_ANGULAR(pJG, JGChen);
                _CLG_EXPORT_ANGULAR(pJG, JGSurf);
                _CLG_EXPORT_ANGULAR(pJG, JGPot);
            }
            break;
            case EDJKS_ChiralAndFermionMomentum:
            {
                _CLG_EXPORT_ANGULAR(pJG, JG);
                _CLG_EXPORT_ANGULAR(pJG, JGS);
                _CLG_EXPORT_ANGULAR(pJG, JGChen);
                _CLG_EXPORT_ANGULAR(pJG, JGSurf);
                _CLG_EXPORT_ANGULAR(pJG, JGPot);

                _CLG_EXPORT_CHIRAL(pCCLight, ChiralKS);
                _CLG_EXPORT_CHIRAL(pCCLight, ConnectSusp);
                _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma3);
                _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma4);
                _CLG_EXPORT_CHIRAL(pCCHeavy, ChiralKS);
                _CLG_EXPORT_CHIRAL(pCCHeavy, ConnectSusp);
                _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma3);
                _CLG_EXPORT_CHIRAL(pCCHeavy, CMTKSGamma4);

                _CLG_EXPORT_CHIRAL(pFALight, OrbitalKS);
                _CLG_EXPORT_CHIRAL(pFALight, SpinKS);
                _CLG_EXPORT_CHIRAL(pFALight, PotentialKS);
                _CLG_EXPORT_CHIRAL(pFAHeavy, OrbitalKS);
                _CLG_EXPORT_CHIRAL(pFAHeavy, SpinKS);
                _CLG_EXPORT_CHIRAL(pFAHeavy, PotentialKS);

                if (uiOmega == iStartOmega)
                {
                    TArray<Real> lstRadius;
                    for (INT i = 0; i < pCCLight->m_lstR.Num(); ++i)
                    {
                        lstRadius.AddItem(F(0.5) * _hostsqrt(static_cast<Real>(pCCLight->m_lstR[i])));
                    }
                    CCString sRadiousFile;
                    sRadiousFile.Format(_T("%s_condensateR.csv"), sCSVSavePrefix.c_str());
                    WriteStringFileRealArray(sRadiousFile, lstRadius);
                }
            }
            break;
            case EDJKS_PlaqutteEnergy:
            {
                //CCString sFileName;
                //sFileName.Format(_T("%s_plaqutte.csv"), sCSVSavePrefix.c_str());
                //WriteStringFileRealArray(sFileName, pPE->m_lstData);
            }
            break;
#endif
            default:
                break;
        }

        appGeneral(_T("\n"));
    }

    switch (eJob)
    {
        case EDJKSEM_Polyakov:
        {
            CCString sFileNameWrite1;
            CCString sFileNameWrite2;
            sFileNameWrite1.Format(_T("%s_polyakov_Nt%d_In.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            sFileNameWrite2.Format(_T("%s_polyakov_Nt%d_Out.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            WriteStringFileComplexArray2(sFileNameWrite1, lstPolyIn);
            WriteStringFileComplexArray2(sFileNameWrite2, lstPolyOut);

            sFileNameWrite1.Format(_T("%s_polyakovZ_Nt%d_In.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            sFileNameWrite2.Format(_T("%s_polyakovZ_Nt%d_Out.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            WriteStringFileComplexArray2(sFileNameWrite1, lstPolyInZ);
            WriteStringFileComplexArray2(sFileNameWrite2, lstPolyOutZ);
        }
        break;
        case EDJKSEM_Chiral:
        {
            //nothing to do
        }
        break;
        case EDJKSEM_BerryPhase:
        {
            //nothing to do
        }
        break;
#if NotYet
        case EDJKS_AngularMomentum:
        {
            //nothing to do
        }
        break;
        case EDJKS_ChiralAndFermionMomentum:
        {
            //nothing to do
        }
        break;
#endif
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


