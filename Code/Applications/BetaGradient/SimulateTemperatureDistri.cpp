//=============================================================================
// FILENAME : SimulateTemperatureDistri.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [mm/dd/yy]
//  [11/06/2024 nbale]
//=============================================================================
#include "BetaGradient.h"

INT SimulateTempDist(CParameters& params)
{
#pragma region Parameters

    appSetupLog(params);

    INT iVaule = 99;
    params.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    UINT iBeforeEquib = static_cast<UINT>(iVaule);

    iVaule = 6;
    params.FetchValueINT(_T("EquvibStep"), iVaule);
    UINT iEquib = static_cast<UINT>(iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("EquvibSkip"), iVaule);
    UINT iEquibSkip = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("SaveStartIndex"), iVaule);
    UINT iSaveStartIndex = static_cast<UINT>(iVaule);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix : %s\n"), sSavePrefix.c_str());

    TArray<CCString> sPrefixList;
    params.FetchStringVectorValue(_T("PrefixList"), sPrefixList);

    iVaule = 0;
    params.FetchValueINT(_T("ListStart"), iVaule);
    const INT iListStart = iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("ListEnd"), iVaule);
    const INT iListEnd = iVaule;

    TArray<DOUBLE> fBeta;
    params.FetchValueArrayDOUBLE(_T("BetaList"), fBeta);
    if (fBeta.Num() != sPrefixList.Num() || sPrefixList.Num() < 1)
    {
        appCrucial(_T("sSavePrefix and fMiddleBeta not corrected!\n"));
        return 0;
    }

    iVaule = 0;
    params.FetchValueINT(_T("Additive"), iVaule);
    UBOOL bAdditive = 0 != iVaule;

    TArray<Real> old_polyakovs;
    params.FetchValueArrayReal(_T("Polyakovs"), old_polyakovs);

    CCString sOldGaugeFileName;
    CCString sOldBosonFileName;
    const UBOOL bHasOldFile = params.FetchStringValue(_T("OldGaugeFileName"), sOldGaugeFileName);

    Real fOldPolyakov = F(0.0);
    Real fOldBosonV = F(0.0);
    if (bHasOldFile)
    {
        params.FetchValueReal(_T("OldPolyakov"), fOldPolyakov);
        params.FetchStringValue(_T("OldBosonFileName"), sOldBosonFileName);
        params.FetchValueReal(_T("OldBosonValue"), fOldBosonV);
    }

    CCString sFileName;
    TCHAR buff1[256];
    TCHAR buff2[256];
    CCString sInfo;

#pragma endregion

    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureBosonValueReal* pBV = dynamic_cast<CMeasureBosonValueReal*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    //TArray<TArray<CLGComplex>> polykovX_nx;
    TArray<CLGComplex> polykov;
    TArray<Real> polykovphase;

    CActionTemperatureDistribution* pGaugeAction = dynamic_cast<CActionTemperatureDistribution*>(appGetLattice()->GetActionById(1));

    pGaugeAction->SetBeta(fBeta[iListStart]);

    UBOOL bNeedBake = TRUE;
    if (!bAdditive && bHasOldFile)
    {
        appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sOldGaugeFileName, EFFT_CLGBin);
        appGetLattice()->m_pBosonField[0]->InitialFieldWithFile(sOldBosonFileName, EFFT_CLGBin);
        pPL->OnConfigurationAccepted(_FIELDS, NULL);
        pBV->OnConfigurationAccepted(_FIELDS, NULL);
        Real fError1 = static_cast<Real>(appAbs(cuCabs(pPL->m_lstLoop[0]) - fOldPolyakov));
#if _CLG_DOUBLEFLOAT
        Real fError2 = appAbs(_cuCabsf(pBV->m_lstEveryConfigurationC[0]) - fOldBosonV);
        if (fError1 < F(1E-07) && fError2 < F(1E-07))
#else
        Real fError2 = appAbs(_cuCabsf(_cToFloat(pBV->m_lstEveryConfigurationC[0])) - fOldBosonV);
        if (fError1 < F(1E-05) && fError2 < F(1E-05))
#endif
        {
            appGeneral(_T("\n ================ Bake using old file\nGauge: %s\nBoson: %s\n=================\n"), sOldGaugeFileName.c_str(), sOldBosonFileName.c_str());
        }
        else
        {
            appGeneral(_T("\n ================ have the initial file, but not matching....\nPolyakov: %2.12f, %2.12f, diff=%f\nBosonValue: %2.12f, %2.12f, diff=%f\n===========\n"),
                cuCabs(pPL->m_lstLoop[0]), fOldPolyakov, fError1,
#if _CLG_DOUBLEFLOAT
                _cuCabsf(pBV->m_lstEveryConfigurationC[0]), fOldBosonV, fError2);
#else
                _cuCabsf(_cToFloat(pBV->m_lstEveryConfigurationC[0])), fOldBosonV, fError2);
#endif
        }
        bNeedBake = FALSE;
    }

    if (bAdditive)
    {
        bNeedBake = FALSE;
    }

    if (bNeedBake && iBeforeEquib > 0)
    {
        if (bHasOldFile)
        {
            appGeneral(_T("!!!! Has old file but still bake !!!!\n"));
        }
        appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));
        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        appGetLattice()->m_pMeasurements->Reset();
        UINT uiAccepCountBeforeE = 0;
        while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
        {
            UINT uiAccepCountBeforeE2 = appGetLattice()->m_pUpdator->Update(1, FALSE);
            if (uiAccepCountBeforeE != uiAccepCountBeforeE2)
            {
                uiAccepCountBeforeE = uiAccepCountBeforeE2;
                pPL->OnConfigurationAccepted(_FIELDS, NULL);
            }
        }
        assert(pPL->m_lstLoop.Num() == static_cast<INT>(iBeforeEquib));
        appPushLogDate(FALSE);
        appGeneral(_T("\n|<P>|,arg<P>={\n"));
        for (INT i = 0; i < pPL->m_lstLoop.Num(); ++i)
        {
            appGeneral(_T("{%f, %f},\n"), cuCabs(pPL->m_lstLoop[i]), cuCarg(pPL->m_lstLoop[i]));
        }
        appGeneral(_T("}\n"));
        appPopLogDate();
    }
    else
    {
        appGeneral(_T("Not Baked\n"));
    }

    for (INT uiOmega = iListStart; uiOmega < sPrefixList.Num() && uiOmega < iListEnd; ++uiOmega)
    {
        CCString sHeader;
        sHeader.Format(_T("%s"), sPrefixList[uiOmega].c_str());
        appSetLogHeader(sHeader);
        appGeneral(_T("\n========= Beta=%f  ==========\n"), fBeta[uiOmega]);
        pGaugeAction->SetBeta(fBeta[uiOmega]);

        if (bAdditive)
        {
            Real fPolyaOld = F(0.0);
            if (old_polyakovs.Num() <= static_cast<INT>(uiOmega))
            {
                appGeneral(_T("\n ================ not have the initial value===========\n"));
                appFailQuitCLG();
                return 1;
            }
            fPolyaOld = old_polyakovs[uiOmega];

            appGetLattice()->m_pMeasurements->Reset();
            sFileName.Format(_T("%sGradient_%s_%d.con"), sSavePrefix.c_str(), sPrefixList[uiOmega].c_str(), iSaveStartIndex);
            appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName, EFFT_CLGBin);
            pPL->OnConfigurationAccepted(_FIELDS, NULL);
            Real fError = static_cast<Real>(appAbs(cuCabs(pPL->m_lstLoop[0]) - fPolyaOld));
#if _CLG_DOUBLEFLOAT
            if (fError < F(1E-07))
#else
            if (fError < F(1E-05))
#endif
            {
                appGeneral(_T("\n ================ using old file start from %d =================\n"), iSaveStartIndex);
            }
            else
            {
                appGeneral(_T("\n ================ have the initial file, but not matching.... %2.12f, %2.12f, diff=%f ===========\n"),
                    cuCabs(pPL->m_lstLoop[0]), fPolyaOld, fError);
                appFailQuitCLG();
                return 1;
            }
        }

        UINT iConfigNumberNow = 0;
        appGetLattice()->m_pMeasurements->Reset();
        appGetLattice()->m_pUpdator->SetConfigurationCount(0);

        while (iConfigNumberNow < iEquib)
        {
            appGetLattice()->m_pUpdator->Update(1, (iConfigNumberNow < iEquibSkip) ? FALSE : TRUE);
            const UINT uiAcce = appGetLattice()->m_pUpdator->GetConfigurationCount();
            if (uiAcce != iConfigNumberNow)
            {
                sFileName.Format(_T("%sGradient_%s_%d"), sSavePrefix.c_str(), sPrefixList[uiOmega].c_str(), uiAcce + iSaveStartIndex);

                //=================================
                //Save config
                const CCString MD5Gauge = appGetLattice()->m_pGaugeField[0]->SaveToFile(sFileName + _T("G.con"));
                const CCString MD5Boson = appGetLattice()->m_pBosonField[0]->SaveToFile(sFileName + _T("B.con"));

                //=================================
                //Save info
                appGetTimeNow(buff1, 256);
                appGetTimeUtc(buff2, 256);
                sInfo.Format(_T("TimeStamp : %d\nTime : %s\nTimeUTC : %s\nGauge MD5 : %s\nBoson MD5 : %s"),
                    appGetTimeStamp(),
                    buff1,
                    buff2,
                    MD5Gauge.c_str(),
                    MD5Boson.c_str());
                sInfo = sInfo + appGetLattice()->GetInfos(_T(""));
                appGetFileSystem()->WriteAllText(sFileName + _T(".txt"), sInfo);

                iConfigNumberNow = uiAcce;
            }
        }

#pragma region gather measurements

        appGetLattice()->m_pMeasurements->Report();

        //===================== Polyakov loop =====================
        assert(pPL->m_lstLoop.Num() == static_cast<INT>(iEquib - iEquibSkip));
        //assert(pPL->m_lstAverageLoopDensity.Num()
        //    == static_cast<INT>(CCommonData::m_sCenter.x));

        //============= polyakov gather =============
        polykov.AddItem(pPL->m_cAverageLoop);
        polykovphase.AddItem(__cuCargf(pPL->m_cAverageLoop));
        //for (UINT iX = 0; iX < static_cast<UINT>(CCommonData::m_sCenter.x); ++iX)
        //{
        //    polykovX_nx[iX].AddItem(pPL->m_lstAverageLoopDensity[iX]);
        //}

#pragma endregion

        appGetLattice()->m_pMeasurements->Reset();
        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    }

    appGeneral(_T("\n========= Finished! ==========\n\n"));

    appQuitCLG(); 

    return 0;
}
