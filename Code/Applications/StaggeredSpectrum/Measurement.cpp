//=============================================================================
// FILENAME : Measurement.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [09/25/2020 nbale]
//=============================================================================

#include "StaggeredSpectrum.h"

static void SaveComplexMatrixNpyAndCsv(const CCString& sBaseFile, const TArray<TArray<cuDoubleComplex>>& data)
{
    TArray<cuDoubleComplex> flat;
    for (INT i = 0; i < data.Num(); ++i)
    {
        for (INT j = 0; j < data[i].Num(); ++j)
        {
            flat.AddItem(data[i][j]);
        }
    }

    TArray<INT> shape;
    shape.AddItem(data.Num());
    shape.AddItem(data.Num() > 0 ? data[0].Num() : 0);
    SaveAsNumpyFile<cuDoubleComplex>(sBaseFile + _T(".npy"), flat.GetData(), shape);
    WriteComplexArray2(sBaseFile + _T(".csv"), data);
}

static void SaveRealMatrixNpyAndCsv(const CCString& sBaseFile, const TArray<TArray<cuDoubleComplex>>& data)
{
    TArray<DOUBLE> flat;
    TArray<TArray<DOUBLE>> realData;
    for (INT i = 0; i < data.Num(); ++i)
    {
        TArray<DOUBLE> oneConf;
        for (INT j = 0; j < data[i].Num(); ++j)
        {
            flat.AddItem(data[i][j].x);
            oneConf.AddItem(data[i][j].x);
        }
        realData.AddItem(oneConf);
    }

    TArray<INT> shape;
    shape.AddItem(data.Num());
    shape.AddItem(data.Num() > 0 ? data[0].Num() : 0);
    SaveAsNumpyFile<DOUBLE>(sBaseFile + _T(".npy"), flat.GetData(), shape);
    WriteRealArray2(sBaseFile + _T(".csv"), realData);
}

static void SaveMesonPArrays(
    const CMeasureMesonCorrelatorStaggered* pMC,
    const CCString& sCSVSavePrefix,
    UINT uiN)
{
    TArray<INT> shape;
    shape.AddItem(_HC_Lti);
    shape.AddItem(8);
    shape.AddItem(8);
    shape.AddItem(8);

    CCString sFile;
    sFile.Format(_T("%s_p_%d.npy"), sCSVSavePrefix.c_str(), uiN);
    SaveAsNumpyFile<cuDoubleComplex>(sFile, pMC->m_pP2PPArray, shape);
    sFile.Format(_T("%s_w2w_p_%d.npy"), sCSVSavePrefix.c_str(), uiN);
    SaveAsNumpyFile<cuDoubleComplex>(sFile, pMC->m_pW2WPArray, shape);
}

static void SaveMesonCorrelators(
    const CMeasureMesonCorrelatorStaggered* pMC,
    const CCString& sCSVSavePrefix)
{
    const INT nConf = pMC->m_lstW2WCombinedCorrelator.Num();
    const INT nt = _HC_Lti;

    for (INT ty = 0; ty < CMeasureMesonCorrelatorStaggered::_kMesonCorrelatorType; ++ty)
    {
        for (INT sub = 0; sub < pMC->m_nSubChannels[ty]; ++sub)
        {
            TArray<TArray<cuDoubleComplex>> p2pData;
            TArray<TArray<cuDoubleComplex>> w2wData;
            for (INT conf = 0; conf < nConf; ++conf)
            {
                TArray<cuDoubleComplex> p2pOneConf;
                TArray<cuDoubleComplex> w2wOneConf;
                for (INT t = 0; t < nt; ++t)
                {
                    p2pOneConf.AddItem(pMC->m_lstP2PCorrelator[conf][ty][sub][t]);
                    w2wOneConf.AddItem(pMC->m_lstW2WCorrelator[conf][ty][sub][t]);
                }
                p2pData.AddItem(p2pOneConf);
                w2wData.AddItem(w2wOneConf);
            }

            CCString sFile;
            sFile.Format(_T("%s_correlationp2p_%d_%d"), sCSVSavePrefix.c_str(), ty, sub);
            SaveComplexMatrixNpyAndCsv(sFile, p2pData);
            sFile.Format(_T("%s_correlationw2w_%d_%d"), sCSVSavePrefix.c_str(), ty, sub);
            SaveComplexMatrixNpyAndCsv(sFile, w2wData);
        }

        TArray<TArray<cuDoubleComplex>> p2pCombined;
        TArray<TArray<cuDoubleComplex>> w2wCombined;
        for (INT conf = 0; conf < nConf; ++conf)
        {
            TArray<cuDoubleComplex> p2pOneConf;
            TArray<cuDoubleComplex> w2wOneConf;
            for (INT t = 0; t < nt; ++t)
            {
                p2pOneConf.AddItem(pMC->m_lstP2PCombinedCorrelator[conf][ty][t]);
                w2wOneConf.AddItem(pMC->m_lstW2WCombinedCorrelator[conf][ty][t]);
            }
            p2pCombined.AddItem(p2pOneConf);
            w2wCombined.AddItem(w2wOneConf);
        }

        CCString sFile;
        sFile.Format(_T("%s_%d"), sCSVSavePrefix.c_str(), ty);
        SaveRealMatrixNpyAndCsv(sFile, p2pCombined);
        sFile.Format(_T("%s_w2w_%d"), sCSVSavePrefix.c_str(), ty);
        SaveRealMatrixNpyAndCsv(sFile, w2wCombined);
    }
}

void AppendStringFile(const CCString& sFileName, const CCString& sContent)
{
    appGetFileSystem()->AppendAllText(sFileName, sContent);
}

#define _CLG_EXPORT_CHIRAL_SINGLE(measureName, lstName) \
CCString sFileNameWrite##measureName##lstName = _T("%s_condensate"); \
CCString sFileNameWrite##measureName##lstName##All = _T("%s_condensate"); \
CCString sFileNameWrite##measureName##lstName##In = _T("%s_condensate"); \
sFileNameWrite##measureName##lstName = sFileNameWrite##measureName##lstName + _T(#measureName) + _T(#lstName) + _T(".csv"); \
sFileNameWrite##measureName##lstName##All = sFileNameWrite##measureName##lstName##All + _T(#measureName) + _T(#lstName) + _T("_All.csv"); \
sFileNameWrite##measureName##lstName##In = sFileNameWrite##measureName##lstName##In + _T(#measureName) + _T(#lstName) + _T("_In.csv"); \
sFileNameWrite##measureName##lstName.Format(sFileNameWrite##measureName##lstName, sCSVSavePrefix.c_str()); \
sFileNameWrite##measureName##lstName##All.Format(sFileNameWrite##measureName##lstName##All, sCSVSavePrefix.c_str()); \
sFileNameWrite##measureName##lstName##In.Format(sFileNameWrite##measureName##lstName##In, sCSVSavePrefix.c_str()); \
TArray<TArray<CLGComplex>> lstName##measureName##OverR; \
TArray<CLGComplex> lstName##measureName##All; \
TArray<CLGComplex> lstName##measureName##In; \
for (UINT j = 0; j < (iEndN - iStartN + 1); ++j) \
{ \
    TArray<CLGComplex> thisConfiguration; \
    for (INT i = 0; i < measureName->m_lstR.Num(); ++i) \
    { \
        thisConfiguration.AddItem(measureName->m_lstCond[lstName][j * measureName->m_lstR.Num() + i]); \
    } \
    lstName##measureName##OverR.AddItem(thisConfiguration); \
    lstName##measureName##All.AddItem(measureName->m_lstCondAll[lstName][j]); \
    lstName##measureName##In.AddItem(measureName->m_lstCondIn[lstName][j]); \
} \
WriteComplexArray2(sFileNameWrite##measureName##lstName, lstName##measureName##OverR); \
WriteComplexArray(sFileNameWrite##measureName##lstName##All, lstName##measureName##All); \
WriteComplexArray(sFileNameWrite##measureName##lstName##In, lstName##measureName##In); 

INT Measurement(CParameters& params)
{

#pragma region read parameters

    appSetupLog(params);

    INT iVaule = 1;
    params.FetchValueINT(_T("StartN"), iVaule);
    UINT iStartN = static_cast<UINT>(iVaule);

    iVaule = 200;
    params.FetchValueINT(_T("EndN"), iVaule);
    UINT iEndN = static_cast<UINT>(iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("DoSmearing"), iVaule);
    UBOOL bDoSmearing = (0 != iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("UseZ4"), iVaule);
    UBOOL bZ4 = 0 != iVaule;

    iVaule = 10;
    params.FetchValueINT(_T("StochasticFieldCount"), iVaule);
    UINT iFieldCount = static_cast<UINT>(iVaule);

    //iVaule = 1;
    //params.FetchValueINT(_T("CheckGaugeFixing"), iVaule);
    //UBOOL bCheckGaugeFixing = 0 != iVaule;

    //iVaule = 0;
    //params.FetchValueINT(_T("UseZ4"), iVaule);
    //UBOOL bZ4 = 0 != iVaule;

    CCString sValue = _T("ESSM_Polyakov");
    params.FetchStringValue(_T("MeasureType"), sValue);
    EStaggeredSpectrumMeasure eJob = __STRING_TO_ENUM(EStaggeredSpectrumMeasure, sValue);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    CCString sCSVSavePrefix;
    params.FetchStringValue(_T("CSVSavePrefix"), sCSVSavePrefix);
    appGeneral(_T("csv save prefix: %s\n"), sCSVSavePrefix.c_str());

    CCString sSubFolderPrefix;
    params.FetchStringValue(_T("SubFolderPrefix"), sSubFolderPrefix);
    appGeneral(_T("sub folder prefix: %s\n"), sSubFolderPrefix.c_str());

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
    CFieldGaugeSU3* pStaple = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField[0]->GetCopy());
    CMeasureWilsonLoop* pPL = dynamic_cast<CMeasureWilsonLoop*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasurePolyakovXY* pPXY = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(6));
    CMeasureMesonCorrelatorStaggered* pMC = dynamic_cast<CMeasureMesonCorrelatorStaggered*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    CMeasureMesonCorrelatorStaggeredSimple2* pMCSimple = dynamic_cast<CMeasureMesonCorrelatorStaggeredSimple2*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
    CMeasureChiralCondensateKS* pCCLight = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(4));
    CMeasureChiralCondensateKS* pCCHeavy = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(5));
    pPL->Reset();
    pPXY->Reset();
    pMC->Reset();
    pMCSimple->Reset();
    pCCLight->Reset();
    pCCHeavy->Reset();
    pCCLight->SetFieldCount(iFieldCount);
    pCCHeavy->SetFieldCount(iFieldCount);

#pragma region Measure

    appGeneral(_T("(*\n"));
    appPushLogDate(FALSE);

    CFieldFermionKSSU3* pF1Light = NULL;
    CFieldFermionKSSU3* pF2Light = NULL;
    CFieldFermionKSSU3* pF1Heavy = NULL;
    CFieldFermionKSSU3* pF2Heavy = NULL;
    if (ESSM_All == eJob || ESSM_Chiral == eJob)
    {
        pF1Light = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(5, _T(__FILE__), __LINE__));
        pF2Light = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(5, _T(__FILE__), __LINE__));
        pF1Heavy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));
        pF2Heavy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));
    }

    for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
    {
        CCString sFileName;
        sFileName.Format(_T("%s_%d.con"), sSavePrefix.c_str(), uiN);
        appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName, eLoadType);

        switch (eJob)
        {
        case ESSM_Polyakov:
        {
            pPXY->OnConfigurationAccepted(_FIELDS, NULL);
        }
        break;
        case ESSM_Wilson:
        {
            if (bDoSmearing)
            {
                appGetLattice()->m_pGaugeField[0]->CalculateOnlyStaple(pStaple);
                appGetLattice()->m_pGaugeSmearing[appGetLattice()->m_pGaugeField[0]->m_byFieldId]->GaugeSmearing(appGetLattice()->m_pGaugeField[0], NULL, pStaple);
            }

            pPL->OnConfigurationAccepted(_FIELDS, NULL);
            if (uiN == iStartN)
            {
                TArray<Real> lstRadius;
                for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                {
                    lstRadius.AddItem(_hostsqrt(static_cast<Real>(pPL->m_lstR[i])));
                }
                CCString sRadiousFile;
                sRadiousFile.Format(_T("%s_VR_R.csv"), sCSVSavePrefix.c_str());
                WriteRealArray(sRadiousFile, lstRadius);
                
            }
        }
        break;
        case ESSM_Correlator:
        {
            pMC->OnConfigurationAccepted(_FIELDS, NULL);
            SaveMesonPArrays(pMC, sCSVSavePrefix, uiN);
        }
        break;
        case ESSM_CorrelatorSimple:
        {
            pMCSimple->OnConfigurationAccepted(_FIELDS, NULL);
        }
        break;
        case ESSM_Chiral:
        {
            for (UINT i = 0; i < iFieldCount; ++i)
            {
                if (NULL != pF1Light)
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

                    pCCLight->OnConfigurationAcceptedZ4(
                        _FIELDS,
                        NULL,
                        pF2Light,
                        pF1Light,
                        0 == i,
                        iFieldCount == i + 1);
                }

                if (NULL != pF1Heavy)
                {
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

                    pCCHeavy->OnConfigurationAcceptedZ4(
                        _FIELDS,
                        NULL,
                        pF2Heavy,
                        pF1Heavy,
                        0 == i,
                        iFieldCount == i + 1);
                }
            }
        }
        break;
        case ESSM_All:
        {
            pMC->OnConfigurationAccepted(_FIELDS, NULL);
            SaveMesonPArrays(pMC, sCSVSavePrefix, uiN);
            pMCSimple->OnConfigurationAccepted(_FIELDS, NULL);

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

                pCCHeavy->OnConfigurationAcceptedZ4(
                    _FIELDS,
                    NULL,
                    pF2Heavy,
                    pF1Heavy,
                    0 == i,
                    iFieldCount == i + 1);
            }

            if (bDoSmearing)
            {
                appGetLattice()->m_pGaugeField[0]->CalculateOnlyStaple(pStaple);
                appGetLattice()->m_pGaugeSmearing[appGetLattice()->m_pGaugeField[0]->m_byFieldId]->GaugeSmearing(appGetLattice()->m_pGaugeField[0], NULL, pStaple);
            }

            pPL->OnConfigurationAccepted(_FIELDS, NULL);
            if (uiN == iStartN)
            {
                TArray<Real> lstRadius;
                for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
                {
                    lstRadius.AddItem(_hostsqrt(static_cast<Real>(pPL->m_lstR[i])));
                }
                CCString sRadiousFile;
                sRadiousFile.Format(_T("%s_VR_R.csv"), sCSVSavePrefix.c_str());
                WriteRealArray(sRadiousFile, lstRadius);
            }
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
            appPushLogDate(FALSE);
            appGeneral(_T("="));
            appPopLogDate();
        }

    }
    appGeneral(_T("\n*)\n"));

    if (ESSM_All == eJob || ESSM_Chiral == eJob)
    {
        if (NULL != pF1Light)
        {
            pF1Light->Return();
            pF2Light->Return();
        }
        if (NULL != pF1Heavy)
        {
            pF1Heavy->Return();
            pF2Heavy->Return();
        }
    }

#pragma endregion

    switch (eJob)
    {
    case ESSM_Polyakov:
    {
        //Write result to file
        //CCString sCSVFile;
        //sCSVFile.Format(_T("%s_polya.csv"), sCSVSavePrefix.c_str());
        //TArray<CLGComplex> polyas;
        //for (INT j = 0; j < pPXY->m_lstLoop.Num(); ++j)
        //{
        //    polyas.AddItem(pPXY->m_lstLoop[j]);
        //}
        //WriteStringFileComplexArray(sCSVFile, polyas);
        pPXY->Export(sCSVSavePrefix, iStartN, iEndN, _T(""), 0, 0);
    }
    break;
    case ESSM_Wilson:
    {
        //Write result to file
        CCString sCSVFile;
        sCSVFile.Format(_T("%s_VR.csv"), sCSVSavePrefix.c_str());
        TArray<TArray<CLGComplex>> vrs;
        for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
        {
            TArray<CLGComplex> thisConfiguration;
            for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
            {
                for (UINT t = 0; t < _HC_Lt / 2; ++t)
                {
                    thisConfiguration.AddItem(pPL->m_lstC[j][i][t]);
                }
            }
            vrs.AddItem(thisConfiguration);
        }
        WriteComplexArray2(sCSVFile, vrs);
    }
    break;
    case ESSM_Correlator:
    {
        pMC->Report();
        SaveMesonCorrelators(pMC, sCSVSavePrefix);
    }
    break;
    case ESSM_CorrelatorSimple:
    {
        for (INT ty = 0; ty < CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2; ++ty)
        {
            CCString sCSVFile;
            sCSVFile.Format(_T("%s_mesonsimple%d.csv"), sCSVSavePrefix.c_str(), ty);
            TArray<TArray<DOUBLE>> res;
            for (INT conf = 0; conf < pMCSimple->m_lstResults.Num(); ++conf)
            {
                TArray<DOUBLE> oneConf;
                for (INT t = 0; t < _HC_Lti; ++t)
                {
                    oneConf.AddItem(pMCSimple->m_lstResults[conf][ty][t]);
                }
                res.AddItem(oneConf);
            }
            WriteRealArray2(sCSVFile, res);
        }
    }
    break;
    case ESSM_Chiral:
    {
		if (NULL != pF1Light)
		{
            _CLG_EXPORT_CHIRAL_SINGLE(pCCLight, ChiralKS);
            _CLG_EXPORT_CHIRAL_SINGLE(pCCLight, ConnectSusp);
            _CLG_EXPORT_CHIRAL_SINGLE(pCCLight, CMTKSGamma3);
            _CLG_EXPORT_CHIRAL_SINGLE(pCCLight, CMTKSGamma4);
		}
        if (NULL != pF1Heavy)
        {
            _CLG_EXPORT_CHIRAL_SINGLE(pCCHeavy, ChiralKS);
            _CLG_EXPORT_CHIRAL_SINGLE(pCCHeavy, ConnectSusp);
            _CLG_EXPORT_CHIRAL_SINGLE(pCCHeavy, CMTKSGamma3);
            _CLG_EXPORT_CHIRAL_SINGLE(pCCHeavy, CMTKSGamma4);
        }
    }
    break;
    case ESSM_All:
    {
        pMC->Report();
        SaveMesonCorrelators(pMC, sCSVSavePrefix);

        for (INT ty = 0; ty < CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2; ++ty)
        {
            CCString sCSVFile;
            sCSVFile.Format(_T("%s_mesonsimple%d.csv"), sCSVSavePrefix.c_str(), ty);
            TArray<TArray<DOUBLE>> res;
            for (INT conf = 0; conf < pMCSimple->m_lstResults.Num(); ++conf)
            {
                TArray<DOUBLE> oneConf;
                for (INT t = 0; t < _HC_Lti; ++t)
                {
                    oneConf.AddItem(pMCSimple->m_lstResults[conf][ty][t]);
                }
                res.AddItem(oneConf);
            }
            WriteRealArray2(sCSVFile, res);
        }

        _CLG_EXPORT_CHIRAL_SINGLE(pCCLight, ChiralKS);
        _CLG_EXPORT_CHIRAL_SINGLE(pCCLight, ConnectSusp);
        _CLG_EXPORT_CHIRAL_SINGLE(pCCLight, CMTKSGamma3);
        _CLG_EXPORT_CHIRAL_SINGLE(pCCLight, CMTKSGamma4);
        _CLG_EXPORT_CHIRAL_SINGLE(pCCHeavy, ChiralKS);
        _CLG_EXPORT_CHIRAL_SINGLE(pCCHeavy, ConnectSusp);
        _CLG_EXPORT_CHIRAL_SINGLE(pCCHeavy, CMTKSGamma3);
        _CLG_EXPORT_CHIRAL_SINGLE(pCCHeavy, CMTKSGamma4);

        //Write result to file
        CCString sCSVFile;
        sCSVFile.Format(_T("%s_VR.csv"), sCSVSavePrefix.c_str());
        TArray<TArray<CLGComplex>> vrs;
        for (UINT j = 0; j < (iEndN - iStartN + 1); ++j)
        {
            TArray<CLGComplex> thisConfiguration;
            for (INT i = 0; i < pPL->m_lstR.Num(); ++i)
            {
                for (UINT t = 0; t < _HC_Lt / 2; ++t)
                {
                    thisConfiguration.AddItem(pPL->m_lstC[j][i][t]);
                }
            }
            vrs.AddItem(thisConfiguration);
        }
        WriteComplexArray2(sCSVFile, vrs);
    }
    break;
    default:
        break;
    }

    appGeneral(_T("\n"));

    appGeneral(_T("\n(*"));
    appPopLogDate();

    appGeneral(_T("\n=====================================\n========= finished! ==========\n*)"));

    appSafeDelete(pStaple);

    appQuitCLG();

    return 0;
}

