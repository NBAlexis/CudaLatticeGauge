//=============================================================================
// FILENAME : Measure.cpp
// 
// DESCRIPTION:
//   Measurement for RotationImprovedFermion.
//   Two measure branches share this file:
//     - hisq measure      : ERIF_MeasureHISQ       -> MeasureHISQ()       (plain HISQ fermion, CFieldFermionHISQSU3)
//     - hisq rot measure  : ERIF_MeasureHISQRotation -> MeasureHISQRotation() (rotating HISQ fermion, CFieldFermionHISQSU3R)
//
//   The chiral condensate measurement follows HISQ/Measure.cpp :
//     - Gaussian/Z4 stochastic source, InverseD with the YAML configured
//       (non-CG) solver, then CMeasureChiralCondensateKS::OnConfigurationAcceptedZ4.
//   For the rotating fermion branch, the fermion field must be created with
//   "CachedGauge : 1" and the StapleCache must cache the rotation link
//   ("CacheRotationLink : 1"), which is refreshed before each measurement.
//   Do NOT use even pseudo-fermion for the measurement fermion (Even : 0).
//
// REVISION:
//  [08/06/2026 nbale]
//=============================================================================

#include "RotationImprovedFermion.h"

/**
 * Shared measurement driver.
 * bRotation selects the pooled fermion type (HISQSU3 vs HISQSU3R) used for the
 * chiral condensate; gauge part (polyakov) is identical for both branches.
 */
static INT MeasureInternal(CParameters& params, UBOOL bRotation)
{
    appSetupLog(params);

#pragma region read parameters

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

    CCString sValue = _T("ERIFMJ_Polyakov");
    params.FetchStringValue(_T("DistributionJob"), sValue);
    ERIFMeasureJob eJob = __STRING_TO_ENUM(ERIFMeasureJob, sValue);

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

    const UINT uiNewLine = (iEndN - iStartN + 1) / 5;
    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureChiralCondensateKS* pCC = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CFieldGaugeSU3* pStaple = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField[0]->GetCopy());

    // The chiral condensate needs two pooled copies of the valence fermion.
    // bRotation == FALSE : CFieldFermionHISQSU3 (plain HISQ, like HISQ/Measure.cpp)
    // bRotation == TRUE  : CFieldFermionHISQSU3R (rotating HISQ, CachedGauge must be 1)
    CFieldFermionKSSU3* pF1 = NULL;
    CFieldFermionKSSU3* pF2 = NULL;
    CFieldFermionKSSU3R* pF1R = NULL;
    CFieldFermionKSSU3R* pF2R = NULL;

    if (ERIFMJ_Chiral == eJob)
    {
        if (bRotation)
        {
            pF1R = dynamic_cast<CFieldFermionKSSU3R*>(appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));
            pF2R = dynamic_cast<CFieldFermionKSSU3R*>(appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));
        }
        else
        {
            pF1 = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));
            pF2 = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));
        }
    }

    appPushLogDate(FALSE);

    for (UINT uiOmega = iListStart; uiOmega < static_cast<UINT>(BetaList.Num()) && uiOmega < static_cast<UINT>(iListEnd); ++uiOmega)
    {
        appGeneral(_T("(* ==== Beta(%f) ========= *)\n"), BetaList[uiOmega]);
        pPL->Reset();
        if (NULL != pCC)
        {
            pCC->Reset();
            pCC->SetFieldCount(iFieldCount);
        }

#pragma region Measure

        appGeneral(_T("(*"));
        for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
        {
            CCString sFileName;
            CCString sTxtFileName;
            if (bSubFolder)
            {
                sFileName.Format(_T("%s/%s/%s%s_%d.con"), sSubFolderName.c_str(), PrefixList[uiOmega].c_str(), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
                sTxtFileName.Format(_T("%s/%s/%s%s_%d.txt"), sSubFolderName.c_str(), PrefixList[uiOmega].c_str(), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
            }
            else
            {
                sFileName.Format(_T("%s%s_%d.con"), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
                sTxtFileName.Format(_T("%s%s_%d.txt"), sSavePrefix.c_str(), PrefixList[uiOmega].c_str(), uiN);
            }

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
                case ERIFMJ_Polyakov:
                {
                    pPL->OnConfigurationAccepted(_FIELDS, NULL);
                }
                break;
                case ERIFMJ_Chiral:
                {
                    CGaugeSmearing* pSmearing = appGetGaugeSmearing(appGetLattice()->m_pGaugeField[0]->m_byFieldId);
                    if (NULL == pSmearing)
                    {
                        appCrucial(_T("No gauge smearing for field %d!\n"), appGetLattice()->m_pGaugeField[0]->m_byFieldId);
                        break;
                    }
                    pSmearing->GaugeSmearingC(appGetLattice()->m_pGaugeField[0]);
                    TArray<const CFieldGauge*> effectivegauge;
                    effectivegauge.AddItem(pSmearing->GetEffectiveGauge());

                    // The rotating fermion reads the rotation link from the staple cache
                    // (CachedGauge : 1), so the cache must be refreshed after the
                    // (HISQ) smearing, exactly like in the HMC integrator.
                    if (bRotation)
                    {
                        CStapleCache* pCache = appGetStapleCache(appGetLattice()->m_pGaugeField[0]->m_byFieldId);
                        if (NULL != pCache)
                        {
                            pCache->Cache(appGetLattice()->m_pGaugeField[0], ECC_BeforeGaugeUpdate);
                            pCache->Cache(appGetLattice()->m_pGaugeField[0], ECC_BeforeAllUpdateBeforeSmearing);
                            pCache->Cache(appGetLattice()->m_pGaugeField[0], ECC_BeforeAllUpdateAfterSmearing);
                        }
                        else
                        {
                            appCrucial(_T("Rotating fermion needs a StapleCache with CacheRotationLink : 1!\n"));
                        }
                    }

                    for (UINT i = 0; i < iFieldCount; ++i)
                    {
                        if (bRotation)
                        {
                            if (bZ4)
                            {
                                pF1R->InitialField(EFIT_RandomZ4);
                            }
                            else
                            {
                                pF1R->InitialField(EFIT_RandomGaussian);
                            }
                            pF1R->FixBoundary(EFB_Field);
                            pF1R->CopyTo(pF2R);
                            pF1R->InverseD(1, 0, 0, effectivegauge.GetData(), NULL, NULL);
                            pF1R->FixBoundary(EFB_Field);

                            pCC->OnConfigurationAcceptedZ4(
                                1, 0, 0, effectivegauge.GetData(), NULL,
                                NULL,
                                NULL,
                                pF2R,
                                pF1R,
                                0 == i,
                                iFieldCount == i + 1);
                        }
                        else
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
                            pF1->InverseD(1, 0, 0, effectivegauge.GetData(), NULL, NULL);
                            pF1->FixBoundary(EFB_Field);

                            pCC->OnConfigurationAcceptedZ4(
                                1, 0, 0, effectivegauge.GetData(), NULL,
                                NULL,
                                NULL,
                                pF2,
                                pF1,
                                0 == i,
                                iFieldCount == i + 1);
                        }
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
                appPushLogDate(FALSE);
                appGeneral(_T("="));
                appPopLogDate();
            }
        }
        appGeneral(_T("\n*)\n"));

#pragma endregion

        switch (eJob)
        {
            case ERIFMJ_Polyakov:
            {
                pPL->Export(sCSVSavePrefix, iStartN, iEndN, PrefixList[uiOmega], uiOmega, iListStart);
            }
            break;
            case ERIFMJ_Chiral:
            {
                _CLG_EXPORT_CHIRAL(pCC, ChiralKS, uiOmega);
                if (pCC->m_bMeasureConnect)
                {
                    _CLG_EXPORT_CHIRAL(pCC, ConnectSusp, uiOmega);
                }
                _CLG_EXPORT_CHIRAL(pCC, CMTKSGamma3, uiOmega);
                _CLG_EXPORT_CHIRAL(pCC, CMTKSGamma4, uiOmega);
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

    if (NULL != pF1)
    {
        pF1->Return();
        pF2->Return();
    }
    if (NULL != pF1R)
    {
        pF1R->Return();
        pF2R->Return();
    }

    appSafeDelete(pStaple);

    appQuitCLG();

    return 0;
}

INT MeasureHISQ(CParameters& params)
{
    return MeasureInternal(params, FALSE);
}

INT MeasureHISQRotation(CParameters& params)
{
    return MeasureInternal(params, TRUE);
}

//=============================================================================
// END OF FILE
//=============================================================================
