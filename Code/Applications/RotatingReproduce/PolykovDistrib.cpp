//=============================================================================
// FILENAME : PolykovDistrib.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [07/09/2019 nbale]
//=============================================================================

#include "RotatingReproduce.h"

__DEFINE_ENUM(EDistributionJob,
    EDJ_Polyakov,
    EDJ_Chiral,
    EDJ_AngularMomentum,
    EDJ_ChiralAndFermionMomentum,
    EDJ_PlaqutteEnergy,
    EDJ_RotatingPlaqutteEnergy,
    EDJ_CheckMD5,
    )


enum { kExportDigital = 20, };


INT MeasurePolyakovDist(CParameters& params)
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

    iVaule = 1;
    params.FetchValueINT(_T("FermionMomentum"), iVaule);
    UBOOL bJF = 0 != iVaule;

    iVaule = 200;
    params.FetchValueINT(_T("EndN"), iVaule);
    UINT iEndN = static_cast<UINT>(iVaule);

    iVaule = 10;
    params.FetchValueINT(_T("StochasticFieldCount"), iVaule);
    UINT iFieldCount = static_cast<UINT>(iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("CheckGaugeFixing"), iVaule);
    UBOOL bCheckGaugeFixing = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("UseZ4"), iVaule);
    UBOOL bZ4 = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("SubFolder"), iVaule);
    UBOOL bSubFolder = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("FreeFermion"), iVaule);
    const UBOOL bFreeFermion = 0 != iVaule;

    CCString sValue = _T("EDJ_Polyakov");
    params.FetchStringValue(_T("DistributionJob"), sValue);
    EDistributionJob eJob = __STRING_TO_ENUM(EDistributionJob, sValue);

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
    TArray<TArray<CLGComplex>> lstPolyInZ;
    TArray<TArray<CLGComplex>> lstPolyOutZ;

    //CCommonData::m_fBeta = fBeta;
    UINT uiNewLine = (iEndN - iStartN + 1) / 5;
    if (uiNewLine < 1)
    {
        uiNewLine = 1;
    }
    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureChiralCondensate* pCC = dynamic_cast<CMeasureChiralCondensate*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CMeasureAMomentumJG* pJG = dynamic_cast<CMeasureAMomentumJG*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
    CMeasureAMomentumStochastic* pJF = dynamic_cast<CMeasureAMomentumStochastic*>(appGetLattice()->m_pMeasurements->GetMeasureById(4));
    CMeasureAction * pPE = dynamic_cast<CMeasureAction*>(appGetLattice()->m_pMeasurements->GetMeasureById(5));
    CMeasureRotatingAction * pRA = dynamic_cast<CMeasureRotatingAction*>(appGetLattice()->m_pMeasurements->GetMeasureById(6));
    CActionGaugePlaquetteRotating* pAG = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->m_pActionList.Num() > 0 ? appGetLattice()->m_pActionList[0] : NULL);
    CFieldFermionWilsonSquareSU3DR* pFermion = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(appGetLattice()->GetFieldById(2));
    appPushLogDate(FALSE);

    pJG->m_fBetaOverN = fBeta / 3.0;

    CFieldFermionWilsonSquareSU3* pF1 = NULL;
    CFieldFermionWilsonSquareSU3* pF2 = NULL;

    if (EDJ_ChiralAndFermionMomentum == eJob
     || (EDJ_AngularMomentum == eJob && bJF)
     || EDJ_Chiral == eJob)
    {
        pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));
        pF2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));
    }

    for (UINT uiOmega = iStartOmega; uiOmega <= iEndOmega; ++uiOmega)
    {
        if (NULL != pAG)
        {
            pAG->SetGaugeOmega(fOmega * uiOmega);
        }
        if (NULL != pFermion)
        {
            pFermion->SetFermionOmega(fOmega * uiOmega);
        }
        appGeneral(_T("(* ==== Omega(%f) ========= *)\n"), fOmega * uiOmega);
        pPL->Reset();
        pCC->Reset();
        pJG->Reset();
        pJF->Reset();
        pCC->SetFieldCount(iFieldCount);
        pJF->SetFieldCount(iFieldCount);
        pRA->Reset();

#pragma region Measure

        appGeneral(_T("(*"));
        for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
        {
            CCString sFileName;
            if (bSubFolder)
            {
                sFileName.Format(_T("%s/O%d/%sRotate_Nt%d_O%d_%d.con"), sSubFolderPrefix.c_str(), uiOmega, sSavePrefix.c_str(), _HC_Lt, uiOmega, uiN);
            }
            else
            {
                sFileName.Format(_T("%sRotate_Nt%d_O%d_%d.con"), sSavePrefix.c_str(), _HC_Lt, uiOmega, uiN);
            }

            if (EDJ_CheckMD5 == eJob)
            {

                break;
            }

            if (bFreeFermion)
            {
                appGetLattice()->m_pGaugeField[0]->InitialField(EFIT_Identity);
            }
            else
            {
                appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName, EFFT_CLGBin);
            }
            
            switch (eJob)
            {
                case EDJ_Polyakov:
                {
                    pPL->OnConfigurationAccepted(_FIELDS, NULL);
                }
                break;
                case EDJ_Chiral:
                {
                    for (UINT i = 0; i < iFieldCount; ++i)
                    {
                        if (bZ4)
                        {
                            pF1->InitialField(EFIT_RandomZ4);
                        }
                        else
                        {
                            pF1->InitialField(EFIT_RandomGaussian);
                        }
                        pF1->FixBoundary(EFB_Field);
                        pF1->CopyTo(pF2);
                        pF1->InverseD(_FIELDS);
                        pF1->FixBoundary(EFB_Field);

                        pCC->OnConfigurationAcceptedZ4(
                            _FIELDS,
                            NULL,
                            pF2,
                            pF1,
                            0 == i,
                            iFieldCount == i + 1);
                    }
                }
                break;
                case EDJ_AngularMomentum:
                {
                    if (bCheckGaugeFixing && NULL != appGetLattice()->m_pGaugeFixing)
                    {
#if !_CLG_DOUBLEFLOAT
                        DOUBLE fError = appGetLattice()->m_pGaugeFixing->CheckRes(appGetLattice()->m_pGaugeField[0]);
                        if (appAbs(fError) > F(0.000001))
#else
                        Real fError = appGetLattice()->m_pGaugeFixing->CheckRes(appGetLattice()->m_pGaugeField[0]);
                        if (appAbs(fError) > F(0.000000000001))
#endif
                        {
                            appGeneral(_T("Bad Gauge Fixing\n"));
                        }
                    }
                    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField[0]);
                    pJG->OnConfigurationAccepted(_FIELDS, NULL);
                    if (bJF)
                    {
                        for (UINT i = 0; i < iFieldCount; ++i)
                        {
                            if (bZ4)
                            {
                                pF1->InitialField(EFIT_RandomZ4);
                            }
                            else
                            {
                                pF1->InitialField(EFIT_RandomGaussian);
                            }
                            pF1->FixBoundary(EFB_Field);
                            pF1->CopyTo(pF2);
                            pF1->InverseD(_FIELDS);

                            pJF->OnConfigurationAcceptedZ4(
                                _FIELDS,
                                NULL, 
                                pF2, 
                                pF1, 
                                0 == i, 
                                iFieldCount == i + 1);
                        }
                    }
                }
                break;
                case EDJ_ChiralAndFermionMomentum:
                {
                    appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField[0]);
                    pJG->OnConfigurationAccepted(_FIELDS, NULL);
                    if (NULL != pJF && NULL != pCC)
                    {
                        for (UINT i = 0; i < iFieldCount; ++i)
                        {
                            if (bZ4)
                            {
                                pF1->InitialField(EFIT_RandomZ4);
                            }
                            else
                            {
                                pF1->InitialField(EFIT_RandomGaussian);
                            }
                            pF1->FixBoundary(EFB_Field);
                            pF1->CopyTo(pF2);
                            pF1->InverseD(_FIELDS);


                            pJF->OnConfigurationAcceptedZ4(
                                _FIELDS,
                                NULL,
                                pF2,
                                pF1,
                                0 == i,
                                iFieldCount == i + 1);

                            pCC->OnConfigurationAcceptedZ4(
                                _FIELDS,
                                NULL,
                                pF2,
                                pF1,
                                0 == i,
                                iFieldCount == i + 1);
                        }
                    }
                }
                break;
                case EDJ_PlaqutteEnergy:
                {
                    pPE->OnConfigurationAccepted(_FIELDS, NULL);
                }
                break;
                case EDJ_RotatingPlaqutteEnergy:
                {
                    pRA->OnConfigurationAccepted(_FIELDS, NULL);
                }
                break;
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
            case EDJ_Polyakov:
            {
                pPL->Export(sCSVSavePrefix, iStartN, iEndN, uiOmega, iStartOmega);
            }
            break;
            case EDJ_Chiral:
            {
                _CLG_EXPORT_CHIRAL_WD(pCC, Chiral);
                _CLG_EXPORT_CHIRAL_WD(pCC, Gamma1);
                _CLG_EXPORT_CHIRAL_WD(pCC, Gamma2);
                _CLG_EXPORT_CHIRAL_WD(pCC, Gamma3);
                _CLG_EXPORT_CHIRAL_WD(pCC, Gamma4);
                _CLG_EXPORT_CHIRAL_WD(pCC, Gamma5);
                _CLG_EXPORT_CHIRAL_WD(pCC, Gamma45);
                _CLG_EXPORT_CHIRAL_WD(pCC, GammaX);
                _CLG_EXPORT_CHIRAL_WD(pCC, GammaY);

                if (uiOmega == iStartOmega)
                {
                    TArray<Real> lstRadius;
                    for (INT i = 0; i < pCC->m_lstR.Num(); ++i)
                    {
                        lstRadius.AddItem(_hostsqrt(static_cast<Real>(pCC->m_lstR[i])));
                    }
                    CCString sRadiousFile;
                    sRadiousFile.Format(_T("%s_condensateR.csv"), sCSVSavePrefix.c_str());
                    WriteRealArray(sRadiousFile, lstRadius);
                }
            }
            break;
            case EDJ_AngularMomentum:
            {
                _CLG_EXPORT_ANGULAR(pJG, JG);
                _CLG_EXPORT_ANGULAR(pJG, JGS2);
                _CLG_EXPORT_ANGULAR(pJG, JGS);
                _CLG_EXPORT_ANGULAR(pJG, JGChen);
                _CLG_EXPORT_ANGULAR(pJG, JGSurf);
                _CLG_EXPORT_ANGULAR(pJG, JGPot);
                //_CLG_EXPORT_ANGULAR(pJG, JGChenApprox);

                if (bJF && NULL != pJF)
                {
                    _CLG_EXPORT_ANGULAR(pJF, JL);
                    _CLG_EXPORT_ANGULAR(pJF, JS);
                    //_CLG_EXPORT_ANGULAR(pJF, JLPure);
                    //_CLG_EXPORT_ANGULAR(pJF, JLJM);
                    _CLG_EXPORT_ANGULAR(pJF, JPot);
                }

                if (uiOmega == iStartOmega)
                {
                    TArray<Real> lstRadius;
                    for (INT i = 0; i < pJG->m_lstR.Num(); ++i)
                    {
                        lstRadius.AddItem(_hostsqrt(static_cast<Real>(pJG->m_lstR[i])));
                    }
                    CCString sRadiousFile;
                    sRadiousFile.Format(_T("%s_angularR.csv"), sCSVSavePrefix.c_str());
                    WriteRealArray(sRadiousFile, lstRadius);
                }
            }
            break;
            case EDJ_ChiralAndFermionMomentum:
            {
                _CLG_EXPORT_ANGULAR(pJG, JG);
                _CLG_EXPORT_ANGULAR(pJG, JGS2);
                _CLG_EXPORT_ANGULAR(pJG, JGS);
                _CLG_EXPORT_ANGULAR(pJG, JGChen);
                _CLG_EXPORT_ANGULAR(pJG, JGSurf);
                _CLG_EXPORT_ANGULAR(pJG, JGPot);
                //_CLG_EXPORT_ANGULAR(pJG, JGChenApprox);

                if (NULL != pJF)
                {
                    _CLG_EXPORT_ANGULAR(pJF, JL);
                    _CLG_EXPORT_ANGULAR(pJF, JS);
                    //_CLG_EXPORT_ANGULAR(pJF, JLPure);
                    //_CLG_EXPORT_ANGULAR(pJF, JLJM);
                    _CLG_EXPORT_ANGULAR(pJF, JPot);
                }

                if (NULL != pCC)
                {
                    _CLG_EXPORT_CHIRAL_WD(pCC, Chiral);
                    _CLG_EXPORT_CHIRAL_WD(pCC, Gamma1);
                    _CLG_EXPORT_CHIRAL_WD(pCC, Gamma2);
                    _CLG_EXPORT_CHIRAL_WD(pCC, Gamma3);
                    _CLG_EXPORT_CHIRAL_WD(pCC, Gamma4);
                    _CLG_EXPORT_CHIRAL_WD(pCC, Gamma5);
                    _CLG_EXPORT_CHIRAL_WD(pCC, Gamma45);
                    _CLG_EXPORT_CHIRAL_WD(pCC, GammaX);
                    _CLG_EXPORT_CHIRAL_WD(pCC, GammaY);
                }

                if (uiOmega == iStartOmega)
                {
                    TArray<Real> lstRadius;
                    for (INT i = 0; i < pJG->m_lstR.Num(); ++i)
                    {
                        lstRadius.AddItem(_hostsqrt(static_cast<Real>(pJG->m_lstR[i])));
                    }
                    CCString sRadiousFile;
                    sRadiousFile.Format(_T("%s_angularR.csv"), sCSVSavePrefix.c_str());
                    WriteRealArray(sRadiousFile, lstRadius);
                }
            }
            break;
            case EDJ_PlaqutteEnergy:
            {
                CCString sFileName;
                sFileName.Format(_T("%s_plaqutte.csv"), sCSVSavePrefix.c_str());
                pPE->WriteRealListToFile(sFileName);
            }
            break;
            case EDJ_RotatingPlaqutteEnergy:
            {
                CCString sFileNameTotal;
                sFileNameTotal.Format(_T("%s_rotating_plaqutte_O%d.csv"), sCSVSavePrefix.c_str(), uiOmega);
                pRA->WriteRealListToFile(sFileNameTotal);

                CCString sFileNameS0;
                sFileNameS0.Format(_T("%s_rotating_plaqutte_S0_O%d.csv"), sCSVSavePrefix.c_str(), uiOmega);
                WriteRealArray(sFileNameS0, pRA->m_lstS0);

                CCString sFileNameS1;
                sFileNameS1.Format(_T("%s_rotating_plaqutte_S1_O%d.csv"), sCSVSavePrefix.c_str(), uiOmega);
                WriteRealArray(sFileNameS1, pRA->m_lstS1);

                CCString sFileNameS2;
                sFileNameS2.Format(_T("%s_rotating_plaqutte_S2_O%d.csv"), sCSVSavePrefix.c_str(), uiOmega);
                WriteRealArray(sFileNameS2, pRA->m_lstS2);
            }
            break;
            default:
                break;
        }

        appGeneral(_T("\n"));
    }

    switch (eJob)
    {
        case EDJ_Polyakov:
        {
            CCString sFileNameWrite1;
            CCString sFileNameWrite2;
            sFileNameWrite1.Format(_T("%s_polyakov_Nt%d_In.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            sFileNameWrite2.Format(_T("%s_polyakov_Nt%d_Out.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
            WriteComplexArray2(sFileNameWrite1, lstPolyIn);
            WriteComplexArray2(sFileNameWrite2, lstPolyOut);

            if (NULL != pPL && pPL->m_bMeasureLoopZ)
            {
                CCString sFileNameWrite3;
                CCString sFileNameWrite4;
                sFileNameWrite3.Format(_T("%s_polyakovZ_Nt%d_In.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
                sFileNameWrite4.Format(_T("%s_polyakovZ_Nt%d_Out.csv"), sCSVSavePrefix.c_str(), _HC_Lt);
                WriteComplexArray2(sFileNameWrite3, lstPolyInZ);
                WriteComplexArray2(sFileNameWrite4, lstPolyOutZ);
            }
        }
        break;
        case EDJ_Chiral:
        {
            //nothing to do
        }
        break;
        case EDJ_AngularMomentum:
        {
            //nothing to do
        }
        break;
        case EDJ_ChiralAndFermionMomentum:
        {
            //nothing to do
        }
        break;
        default:
            break;
    }

    appGeneral(_T("\n(*"));
    appPopLogDate();

    appGeneral(_T("\n=====================================\n========= finished! ==========\n*)"));
    
    if (NULL != pF1)
    {
        pF1->Return();
        pF2->Return();
    }

    appQuitCLG();

    return 0;
}


