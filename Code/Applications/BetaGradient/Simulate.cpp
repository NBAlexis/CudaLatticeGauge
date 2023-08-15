//=============================================================================
// FILENAME : Simulate.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [08/17/2022 nbale]
//=============================================================================
#include "BetaGradient.h"

INT Simulate(CParameters& params)
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

    iVaule = 0;
    params.FetchValueINT(_T("BetaStride"), iVaule);
    const INT iBetaStride = iVaule;

#if _CLG_DOUBLEFLOAT

    TArray<Real> fMiddleBeta;
    params.FetchValueArrayReal(_T("MiddleBetaList"), fMiddleBeta);
    if (fMiddleBeta.Num() != sPrefixList.Num() || sPrefixList.Num() < 1)
    {
        appCrucial(_T("sSavePrefix and fMiddleBeta not corrected!\n"));
        return 0;
    }

    Real fDeltaBeta = 3.0;
    params.FetchValueReal(_T("DeltaBeta"), fDeltaBeta);

#else

    TArray<DOUBLE> fMiddleBeta;
    params.FetchValueArrayDOUBLE(_T("MiddleBetaList"), fMiddleBeta);
    if (fMiddleBeta.Num() != sPrefixList.Num() || sPrefixList.Num() < 1)
    {
        appCrucial(_T("sSavePrefix and fMiddleBeta not corrected!\n"));
        return 0;
    }

    DOUBLE fDeltaBeta = 3.0;
    params.FetchValueDOUBLE(_T("DeltaBeta"), fDeltaBeta);

#endif

    iVaule = 0;
    params.FetchValueINT(_T("Additive"), iVaule);
    UBOOL bAdditive = 0 != iVaule;

    TArray<Real> old_polyakovs;
    params.FetchValueArrayReal(_T("Polyakovs"), old_polyakovs);

    CCString sOldFileName;
    const UBOOL bHasOldFile = params.FetchStringValue(_T("OldFileName"), sOldFileName);

    Real fOldPolyakov = F(0.0);
    if (bHasOldFile)
    {
        params.FetchValueReal(_T("OldPolyakov"), fOldPolyakov);
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
    TArray<TArray<CLGComplex>> polykovX_nx;
    TArray<CLGComplex> polykov;
    TArray<Real> polykovphase;
    for (UINT uiX = 0; uiX < static_cast<UINT>(CCommonData::m_sCenter.x); ++uiX)
    {
        TArray<CLGComplex> a;
        TArray<Real> b;
        polykovX_nx.AddItem(a);
    }

    CActionGaugePlaquetteGradient* pGaugeGradient = dynamic_cast<CActionGaugePlaquetteGradient*>(appGetLattice()->GetActionById(1));
    TArray<DOUBLE> betaArray;
    if (0 == iBetaStride)
    {
        for (INT i = 0; i < _HC_Lzi; ++i)
        {
            betaArray.AddItem(sin((-1.0 + 2.0 * i / static_cast<DOUBLE>(_HC_Lzi)) * PI) * fDeltaBeta + fMiddleBeta[0]);
        }
    }
    else
    {
        if (15 != _HC_Lzi)
        {
            appCrucial(_T("Lz must be 11"));
        }
        DOUBLE upper = fMiddleBeta[0] + fDeltaBeta;
        DOUBLE lower = fMiddleBeta[0] - fDeltaBeta;

        if (1 == iBetaStride)
        {
            // ----- . ++++++
            for (INT i = 0; i < 7; ++i)
            {
                betaArray.AddItem(lower);
            }

            betaArray.AddItem(fMiddleBeta[0]);

            for (INT i = 0; i < 7; ++i)
            {
                betaArray.AddItem(upper);
            }
        }
        else if (2 == iBetaStride)
        {
            // ---- ... +++++
            for (INT i = 0; i < 6; ++i)
            {
                betaArray.AddItem(lower);
            }

            for (INT i = -1; i <= 1; ++i)
            {
                betaArray.AddItem(fMiddleBeta[0] + i * 0.5 * fDeltaBeta);
            }

            for (INT i = 0; i < 6; ++i)
            {
                betaArray.AddItem(upper);
            }
        }
        else if (3 == iBetaStride)
        {
            // --- ..... +++
            for (INT i = 0; i < 5; ++i)
            {
                betaArray.AddItem(lower);
            }

            for (INT i = -2; i <= 2; ++i)
            {
                betaArray.AddItem(fMiddleBeta[0] + i * 0.333333333333333 * fDeltaBeta);
            }

            for (INT i = 0; i < 5; ++i)
            {
                betaArray.AddItem(upper);
            }
        }
        else if (4 == iBetaStride)
        {
            // --- ..... +++
            for (INT i = 0; i < 4; ++i)
            {
                betaArray.AddItem(lower);
            }

            for (INT i = -3; i <= 3; ++i)
            {
                betaArray.AddItem(fMiddleBeta[0] + i * 0.25 * fDeltaBeta);
            }

            for (INT i = 0; i < 4; ++i)
            {
                betaArray.AddItem(upper);
            }
        }
    }
    pGaugeGradient->SetBeta(betaArray);

    UBOOL bNeedBake = TRUE;
    if (!bAdditive && bHasOldFile)
    {
        appGetLattice()->m_pGaugeField->InitialFieldWithFile(sOldFileName, EFFT_CLGBin);
        pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
        Real fError = appAbs(_cuCabsf(pPL->m_lstLoop[0]) - fOldPolyakov);
#if _CLG_DOUBLEFLOAT
        if (fError < F(1E-07))
#else
        if (fError < F(1E-05))
#endif
        {
            appGeneral(_T("\n ================ Bake using old file %s =================\n"), sOldFileName.c_str());
        }
        else
        {
            appGeneral(_T("\n ================ have the initial file, but not matching.... %2.12f, %2.12f, diff=%f ===========\n"),
                _cuCabsf(pPL->m_lstLoop[0]), fOldPolyakov, fError);
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
                pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
            }
        }
        assert(pPL->m_lstLoop.Num() == static_cast<INT>(iBeforeEquib));
        appSetLogDate(FALSE);
        appGeneral(_T("\n|<P>|,arg<P>={\n"));
        for (INT i = 0; i < pPL->m_lstLoop.Num(); ++i)
        {
            appGeneral(_T("{%f, %f},\n"), _cuCabsf(pPL->m_lstLoop[i]), __cuCargf(pPL->m_lstLoop[i]));
        }
        appGeneral(_T("}\n"));
        appSetLogDate(TRUE);
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
        appGeneral(_T("\n========= Middle Beta=%f  ==========\n"), fMiddleBeta[uiOmega]);

        betaArray.RemoveAll();
        if (0 == iBetaStride)
        {
            for (INT i = 0; i < _HC_Lzi; ++i)
            {
                betaArray.AddItem(sin((-1.0 + 2.0 * i / static_cast<DOUBLE>(_HC_Lzi)) * PI) * fDeltaBeta + fMiddleBeta[uiOmega]);
            }
        }
        else
        {
            if (15 != _HC_Lzi)
            {
                appCrucial(_T("Lz must be 11"));
            }
            DOUBLE upper = fMiddleBeta[uiOmega] + fDeltaBeta;
            DOUBLE lower = fMiddleBeta[uiOmega] - fDeltaBeta;

            if (1 == iBetaStride)
            {
                // ----- . ++++++
                for (INT i = 0; i < 7; ++i)
                {
                    betaArray.AddItem(lower);
                }

                betaArray.AddItem(fMiddleBeta[uiOmega]);

                for (INT i = 0; i < 7; ++i)
                {
                    betaArray.AddItem(upper);
                }
            }
            else if (2 == iBetaStride)
            {
                // ---- ... +++++
                for (INT i = 0; i < 6; ++i)
                {
                    betaArray.AddItem(lower);
                }

                for (INT i = -1; i <= 1; ++i)
                {
                    betaArray.AddItem(fMiddleBeta[uiOmega] + i * 0.5 * fDeltaBeta);
                }

                for (INT i = 0; i < 6; ++i)
                {
                    betaArray.AddItem(upper);
                }
            }
            else if (3 == iBetaStride)
            {
                // --- ..... +++
                for (INT i = 0; i < 5; ++i)
                {
                    betaArray.AddItem(lower);
                }

                for (INT i = -2; i <= 2; ++i)
                {
                    betaArray.AddItem(fMiddleBeta[uiOmega] + i * 0.333333333333333 * fDeltaBeta);
                }

                for (INT i = 0; i < 5; ++i)
                {
                    betaArray.AddItem(upper);
                }
            }
            else if (4 == iBetaStride)
            {
                // --- ..... +++
                for (INT i = 0; i < 4; ++i)
                {
                    betaArray.AddItem(lower);
                }

                for (INT i = -3; i <= 3; ++i)
                {
                    betaArray.AddItem(fMiddleBeta[uiOmega] + i * 0.25 * fDeltaBeta);
                }

                for (INT i = 0; i < 4; ++i)
                {
                    betaArray.AddItem(upper);
                }
            }
        }
        pGaugeGradient->SetBeta(betaArray);

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
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sFileName, EFFT_CLGBin);
            pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
            Real fError = appAbs(_cuCabsf(pPL->m_lstLoop[0]) - fPolyaOld);
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
                    _cuCabsf(pPL->m_lstLoop[0]), fPolyaOld, fError);
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
                const CCString MD5 = appGetLattice()->m_pGaugeField->SaveToFile(sFileName + _T(".con"));

                //=================================
                //Save info
                appGetTimeNow(buff1, 256);
                appGetTimeUtc(buff2, 256);
                sInfo.Format(_T("TimeStamp : %d\nTime : %s\nTimeUTC : %s\nMD5 : %s\n"),
                    appGetTimeStamp(),
                    buff1,
                    buff2,
                    MD5.c_str());
                sInfo = sInfo + appGetLattice()->GetInfos(_T(""));
                appGetFileSystem()->WriteAllText(sFileName + _T(".txt"), sInfo);

                iConfigNumberNow = uiAcce;
            }
        }

#pragma region gather measurements

        appGetLattice()->m_pMeasurements->Report();

        //===================== Polyakov loop =====================
        assert(pPL->m_lstLoop.Num() == static_cast<INT>(iEquib - iEquibSkip));
        assert(pPL->m_lstAverageLoopDensity.Num()
            == static_cast<INT>(CCommonData::m_sCenter.x));

        //============= polyakov gather =============
        polykov.AddItem(pPL->m_cAverageLoop);
        polykovphase.AddItem(__cuCargf(pPL->m_cAverageLoop));
        for (UINT iX = 0; iX < static_cast<UINT>(CCommonData::m_sCenter.x); ++iX)
        {
            polykovX_nx[iX].AddItem(pPL->m_lstAverageLoopDensity[iX]);
        }

#pragma endregion

        appGetLattice()->m_pMeasurements->Reset();
        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    }

    appGeneral(_T("\n========= Finished! ==========\n\n"));

    appQuitCLG(); 

    return 0;
}
