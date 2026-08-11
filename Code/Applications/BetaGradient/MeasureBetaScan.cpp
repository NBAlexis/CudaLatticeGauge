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
    EBSMJ_Angular,
    EBSMJ_Meson,
    EBSMJ_MesonSimple,
    EBSMJ_DoubleToFloat,
    EBSMJ_WilsonPath,
    EBSMJ_WilsonPathTwoPath,
    )


/**
* read text file line-by-line
* ignore # started lines
* for each line, change string 1, 2, 3, 4 to TArray<INT>
*/
TArray<TArray<INT>> ParseWilsonPath(const CCString& sFileName)
{
    TArray<TArray<INT>> ret;
    CCString sFileContent = appGetFileSystem()->ReadAllText(sFileName);
    TArray<CCString> sLines = appGetStringList(sFileContent, _T('\n'), 15);
    for (INT i = 0; i < sLines.Num(); ++i)
    {
        if (!sLines[i].IsEmpty() && sLines[i].GetAt(0) != _T('#'))
        {
            TArray<CCString> sNumbers = appGetStringList(sLines[i], _T(','), 15);
            TArray<INT> iNumbers;
            for (INT j = 0; j < sNumbers.Num(); ++j)
            {
                iNumbers.AddItem(appStrToINT(sNumbers[j]));
            }
            ret.AddItem(iNumbers);
        }
    }
    return ret;
}

/**
* 1 - up
* 2 - down
* 3 - right
* 4 - left
*/
TArray<SCHAR> GetOnePath(const TArray<INT>& dirs, SCHAR mu, SCHAR nu)
{
    TArray<SCHAR> ret;
    for (INT i = 0; i < dirs.Num(); ++i)
    {
        if (1 == dirs[i])
        {
            ret.AddItem(-nu);
        }
        else if (2 == dirs[i])
        {
            ret.AddItem(nu);
        }
        else if (3 == dirs[i])
        {
            ret.AddItem(mu);
        }
        else if (4 == dirs[i])
        {
            ret.AddItem(-mu);
        }
    }
    return ret;
}

/**
* 1 - up
* 2 - down
* 3 - right
* 4 - left
* 5 - in
* 6 - out
*/
TArray<SCHAR> GetOnePath3D(const TArray<INT>& dirs, SCHAR mu, SCHAR nu, SCHAR rho)
{
    TArray<SCHAR> ret;
    for (INT i = 0; i < dirs.Num(); ++i)
    {
        if (1 == dirs[i])
        {
            ret.AddItem(-nu);
        }
        else if (2 == dirs[i])
        {
            ret.AddItem(nu);
        }
        else if (3 == dirs[i])
        {
            ret.AddItem(mu);
        }
        else if (4 == dirs[i])
        {
            ret.AddItem(-mu);
        }
        else if (5 == dirs[i])
        {
            ret.AddItem(rho);
        }
        else if (6 == dirs[i])
        {
            ret.AddItem(-rho);
        }
    }
    return ret;
}

void GetTwoPath(INT maxK, INT maxL, TArray<TArray<SCHAR>> &path, TArray<SSmallInt4>& shift)
{
    SCHAR mulst[3] = {1, 1, 2};
    SCHAR nulst[3] = {2, 3, 3};
    SCHAR shiftlst[3] = {3, 2, 1};
    for (INT k = 1; k <= maxK; ++k)
    {
        for (INT mu = 0; mu < 3; ++mu)
        {
            TArray<SCHAR> onePath;
            for (INT p = 0; p < k; ++p)
            {
                onePath.AddItem(mulst[mu]);
            }
            for (INT p = 0; p < k; ++p)
            {
                onePath.AddItem(nulst[mu]);
            }
            for (INT p = 0; p < k; ++p)
            {
                onePath.AddItem(-mulst[mu]);
            }
            for (INT p = 0; p < k; ++p)
            {
                onePath.AddItem(-nulst[mu]);
            }

            //2D shift
            for (INT shiftx = 0; shiftx <= maxL; ++shiftx)
            {
                for (INT shifty = 0; shifty <= maxL; ++shifty)
                {
                    if (shiftx >= k || shifty >= k)
                    {
                        path.AddItem(onePath);
                        SSmallInt4 sShift;
                        sShift.m_byData4[0] = 0;
                        sShift.m_byData4[1] = 0;
                        sShift.m_byData4[2] = 0;
                        sShift.m_byData4[3] = 0;
                        sShift.m_byData4[mulst[mu] - 1] = static_cast<SCHAR>(shiftx);
                        sShift.m_byData4[nulst[mu] - 1] = static_cast<SCHAR>(shifty);
                        shift.AddItem(sShift);
                    }
                }
            }

            //3D shift
            for (INT shiftz = 1; shiftz <= maxL; ++shiftz)
            {
                path.AddItem(onePath);
                SSmallInt4 sShift;
                sShift.m_byData4[0] = 0;
                sShift.m_byData4[1] = 0;
                sShift.m_byData4[2] = 0;
                sShift.m_byData4[3] = 0;
                sShift.m_byData4[shiftlst[mu] - 1] = static_cast<SCHAR>(shiftz);
                shift.AddItem(sShift);
            }
        }
    }
}

static CCString MesonOutputBase(const CCString& sCSVSavePrefix, const CCString& sBetaPrefix)
{
    CCString sBase;
    sBase.Format(_T("%s_%s"), sCSVSavePrefix.c_str(), sBetaPrefix.c_str());
    return sBase;
}

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
    const CCString& sBetaPrefix,
    UINT uiN)
{
    TArray<INT> shape;
    shape.AddItem(_HC_Lti);
    shape.AddItem(8);
    shape.AddItem(8);
    shape.AddItem(8);

    const CCString sBase = MesonOutputBase(sCSVSavePrefix, sBetaPrefix);
    CCString sFile;
    sFile.Format(_T("%s_p_%d.npy"), sBase.c_str(), uiN);
    SaveAsNumpyFile<cuDoubleComplex>(sFile, pMC->m_pP2PPArray, shape);
    sFile.Format(_T("%s_w2w_p_%d.npy"), sBase.c_str(), uiN);
    SaveAsNumpyFile<cuDoubleComplex>(sFile, pMC->m_pW2WPArray, shape);
}

static void SaveMesonCorrelators(
    const CMeasureMesonCorrelatorStaggered* pMC,
    const CCString& sCSVSavePrefix,
    const CCString& sBetaPrefix)
{
    const INT nConf = pMC->m_lstW2WCombinedCorrelator.Num();
    const INT nt = _HC_Lti;
    const CCString sBase = MesonOutputBase(sCSVSavePrefix, sBetaPrefix);

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
            sFile.Format(_T("%s_correlationp2p_%d_%d"), sBase.c_str(), ty, sub);
            SaveComplexMatrixNpyAndCsv(sFile, p2pData);
            sFile.Format(_T("%s_correlationw2w_%d_%d"), sBase.c_str(), ty, sub);
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
        sFile.Format(_T("%s_%d"), sBase.c_str(), ty);
        SaveRealMatrixNpyAndCsv(sFile, p2pCombined);
        sFile.Format(_T("%s_w2w_%d"), sBase.c_str(), ty);
        SaveRealMatrixNpyAndCsv(sFile, w2wCombined);
    }
}
    
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

    CCString sWilsonPathFile;
    params.FetchStringValue(_T("WilsonPathFile"), sWilsonPathFile);
    appGeneral(_T("Wilson Path File: %s\n"), sWilsonPathFile.c_str());

    iVaule = 3;
    params.FetchValueINT(_T("WilsonPathType"), iVaule);
    const UINT iWilsonPathDim = static_cast<UINT>(iVaule);

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
    CMeasureAMomentumJG* pAMJG = dynamic_cast<CMeasureAMomentumJG*>(appGetLattice()->m_pMeasurements->GetMeasureById(4));
    CFieldGaugeSU3* pStaple = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField[0]->GetCopy());

    CMeasureChiralCondensateKS* pCCLight = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    //CMeasureChiralCondensateKS* pCCHeavy = dynamic_cast<CMeasureChiralCondensateKS*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));

    CActionGaugePlaquette* pAG = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->m_pActionList.Num() > 0 ? appGetLattice()->m_pActionList[0] : NULL);

    CMeasureMesonCorrelatorStaggered* pMC = dynamic_cast<CMeasureMesonCorrelatorStaggered*>(appGetLattice()->m_pMeasurements->GetMeasureById(5));
    CMeasureMesonCorrelatorStaggeredSimple2* pMCSimple = dynamic_cast<CMeasureMesonCorrelatorStaggeredSimple2*>(appGetLattice()->m_pMeasurements->GetMeasureById(6));

    CMeasureWilsonLoopWithPath* pWilsonPath = dynamic_cast<CMeasureWilsonLoopWithPath*>(appGetLattice()->m_pMeasurements->GetMeasureById(7));

    CFieldFermionKSSU3* pF1Light = NULL;
    CFieldFermionKSSU3* pF2Light = NULL;
    //CFieldFermionKSSU3* pF1Heavy = NULL;
    //CFieldFermionKSSU3* pF2Heavy = NULL;


    if (EBSMJ_Chiral == eJob)
    {
        pF1Light = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));
        pF2Light = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));
        //pF1Heavy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(3));
        //pF2Heavy = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(3));
    }

    if (EBSMJ_WilsonPath == eJob)
    {
        TArray<TArray<SCHAR>> wilsonPaths;
        TArray<TArray<INT>> measurePathes = ParseWilsonPath(sWilsonPathFile);
        for (INT i = 0; i < measurePathes.Num(); ++i)
        {
            if (2 == iWilsonPathDim)
            {
                for (SCHAR mu = 1; mu <= 4; ++mu)
                {
                    for (SCHAR nu = mu + 1; nu <= 4; ++nu)
                    {
                        wilsonPaths.AddItem(GetOnePath(measurePathes[i], mu, nu));
                    }
                }
            }
            else if (3 == iWilsonPathDim)
            {
                for (SCHAR skip = 1; skip <= 4; ++skip)
                {
                    for (SCHAR mu = 1; mu <= 4; ++mu)
                    {
                        if (mu == skip)
                        {
                            continue;
                        }
                        for (SCHAR nu = 1; nu <= 4; ++nu)
                        {
                            if (nu == skip || nu == mu)
                            {
                                continue;
                            }
                            for (SCHAR rho = 1; rho <= 4; ++rho)
                            {
                                if (rho == skip || rho == mu || rho == nu)
                                {
                                    continue;
                                }
                                wilsonPaths.AddItem(GetOnePath3D(measurePathes[i], mu, nu, rho));
                            }
                        }
                    }
                }
            }
            else if (4 == iWilsonPathDim)
            {
                wilsonPaths.AddItem(GetOnePath(measurePathes[i], 1, 2));
                wilsonPaths.AddItem(GetOnePath(measurePathes[i], 2, 3));
                wilsonPaths.AddItem(GetOnePath(measurePathes[i], 1, 3));
            }
            else if (5 == iWilsonPathDim)
            {
                wilsonPaths.AddItem(GetOnePath3D(measurePathes[i], 1, 2, 3));
                wilsonPaths.AddItem(GetOnePath3D(measurePathes[i], 1, 3, 2));
                //wilsonPaths.AddItem(GetOnePath3D(measurePathes[i], 2, 1, 3));
                wilsonPaths.AddItem(GetOnePath3D(measurePathes[i], 2, 3, 1));
                //wilsonPaths.AddItem(GetOnePath3D(measurePathes[i], 3, 1, 2));
                //wilsonPaths.AddItem(GetOnePath3D(measurePathes[i], 3, 2, 1));
            }
            else if (6 == iWilsonPathDim)
            {
                wilsonPaths.AddItem(GetOnePath(measurePathes[i], 1, 2));
                wilsonPaths.AddItem(GetOnePath(measurePathes[i], 2, 1));
                wilsonPaths.AddItem(GetOnePath(measurePathes[i], 2, 3));
                wilsonPaths.AddItem(GetOnePath(measurePathes[i], 3, 2));
                wilsonPaths.AddItem(GetOnePath(measurePathes[i], 1, 3));
                wilsonPaths.AddItem(GetOnePath(measurePathes[i], 3, 1));
            }
        }
        pWilsonPath->SetPath(wilsonPaths);
    }
    else if (EBSMJ_WilsonPathTwoPath == eJob)
    {
        TArray<TArray<SCHAR>> wilsonPaths;
        TArray<SSmallInt4> shifts;
        GetTwoPath(_HC_Lx / 4, _HC_Lx / 2, wilsonPaths, shifts);
        pWilsonPath->SetAsTwoPathBackForward(TRUE, wilsonPaths, shifts);
    }

    appPushLogDate(FALSE);

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
        pAMJG->Reset();
        if (NULL != pAMJG)
        {
            pAMJG->m_fBetaOverN = BetaList[uiOmega] / 3.0;
        }

        pCCLight->SetFieldCount(iFieldCount);
        //pCCHeavy->SetFieldCount(iFieldCount);

        pMC->Reset();
        pMCSimple->Reset();
        pWilsonPath->Reset();

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
            
            appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName, eLoadType);

            switch (eJob)
            {
                case EBSMJ_Polyakov:
                {
                    pPL->OnConfigurationAccepted(_FIELDS, NULL);
                }
                break;
                case EBSMJ_Wilson:
                {
                    appGetLattice()->m_pGaugeField[0]->CalculateOnlyStaple(pStaple);
                    appGetLattice()->m_pGaugeSmearing[appGetLattice()->m_pGaugeField[0]->m_byFieldId]->GaugeSmearing(appGetLattice()->m_pGaugeField[0], NULL, pStaple);
                    pWL->OnConfigurationAccepted(_FIELDS, NULL);
                    if (uiN == iStartN)
                    {
                        TArray<Real> lstRadius;
                        for (INT i = 0; i < pWL->m_lstR.Num(); ++i)
                        {
                            lstRadius.AddItem(_hostsqrt(static_cast<Real>(pWL->m_lstR[i])));
                        }
                        CCString sRadiousFile;
                        sRadiousFile.Format(_T("%s_VR_R.csv"), sCSVSavePrefix.c_str());
                        WriteRealArray(sRadiousFile, lstRadius);
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
                        pF1Light->FixBoundary(EFB_Field);
                        pF1Light->CopyTo(pF2Light);
                        pF1Light->InverseD(_FIELDS);
                        pF1Light->FixBoundary(EFB_Field);
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
                case EBSMJ_Angular:
                    {
                        appGetLattice()->SetAPhys(appGetLattice()->m_pGaugeField[0]);
                        pAMJG->OnConfigurationAccepted(_FIELDS, NULL);
                    }
                    break;

                case EBSMJ_Meson:
                    {
                        pMC->OnConfigurationAccepted(_FIELDS, NULL);
                        SaveMesonPArrays(pMC, sCSVSavePrefix, PrefixList[uiOmega], uiN);
                    }
                    break;
                case EBSMJ_MesonSimple:
                    {
                        pMCSimple->OnConfigurationAccepted(_FIELDS, NULL);
                    }
                    break;
                case EBSMJ_WilsonPath:
                    {
                        pWilsonPath->OnConfigurationAccepted(_FIELDS, NULL);
                    }
                    break;
                case EBSMJ_WilsonPathTwoPath:
                    {
                        pWilsonPath->OnConfigurationAccepted(_FIELDS, NULL);
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

        // Call Report() for meson correlator text log output
        if (EBSMJ_Meson == eJob)
        {
            pMC->Report();
            SaveMesonCorrelators(pMC, sCSVSavePrefix, PrefixList[uiOmega]);
        }

        switch (eJob)
        {
            case EBSMJ_Polyakov:
            {
                pPL->Export(sCSVSavePrefix, iStartN, iEndN, PrefixList[uiOmega], uiOmega, iListStart);
            }
            break;
            case EBSMJ_Chiral:
            {
                _CLG_EXPORT_CHIRAL(pCCLight, ChiralKS, uiOmega);
                if (pCCLight->m_bMeasureConnect)
                {
                    _CLG_EXPORT_CHIRAL(pCCLight, ConnectSusp, uiOmega);
                }
                
                _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma3, uiOmega);
                _CLG_EXPORT_CHIRAL(pCCLight, CMTKSGamma4, uiOmega);

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
                sCSVFile.Format(_T("%s_VR_Nt%d_%s.csv"), sCSVSavePrefix.c_str(), _HC_Lt, PrefixList[uiOmega].c_str());
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
                WriteComplexArray2(sCSVFile, vrs);
            }
            break;
            case EBSMJ_Angular:
                {
                    _CLG_EXPORT_ANGULAR(pAMJG, JG, uiOmega, O);
                    _CLG_EXPORT_ANGULAR(pAMJG, JGS2, uiOmega, O);
                    _CLG_EXPORT_ANGULAR(pAMJG, JGS, uiOmega, O);
                    _CLG_EXPORT_ANGULAR(pAMJG, JGChen, uiOmega, O);
                    _CLG_EXPORT_ANGULAR(pAMJG, JGSurf, uiOmega, O);
                    _CLG_EXPORT_ANGULAR(pAMJG, JGPot, uiOmega, O);
                }
                break;
            case EBSMJ_MesonSimple:
            {
                for (INT ty = 0; ty < CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2; ++ty)
                {
                    CCString sCSVFile;
                    sCSVFile.Format(_T("%s_mesonsimple%d_%s.csv"), sCSVSavePrefix.c_str(), ty, PrefixList[uiOmega].c_str());
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
            case EBSMJ_WilsonPath:
                {
                    CCString sCSVFile;
                    sCSVFile.Format(_T("%s_%s_wilsonloops.csv"), sCSVSavePrefix.c_str(), PrefixList[uiOmega].c_str());
                    WriteComplexArray2(sCSVFile, pWilsonPath->m_lstV);
                }
                break;
            case EBSMJ_WilsonPathTwoPath:
                {
                    CCString sCSVFile;
                    sCSVFile.Format(_T("%s_%s_wilsonloops_c.csv"), sCSVSavePrefix.c_str(), PrefixList[uiOmega].c_str());
                    WriteComplexArray2(sCSVFile, pWilsonPath->m_lstV);

                    sCSVFile.Format(_T("%s_%s_wilsonloops_d.csv"), sCSVSavePrefix.c_str(), PrefixList[uiOmega].c_str());
                    WriteRealArray2(sCSVFile, pWilsonPath->m_lstDV);
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

    appSafeDelete(pStaple);

    appQuitCLG();

    return 0;
}
