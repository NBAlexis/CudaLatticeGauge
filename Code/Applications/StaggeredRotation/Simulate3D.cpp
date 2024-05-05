//=============================================================================
// FILENAME : Simulate.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [09/24/2020 nbale]
//=============================================================================

#include "StaggeredRotation.h"

INT SimulateStaggeredRotation3D(CParameters& params)
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

    iVaule = 250;
    params.FetchValueINT(_T("OmegaSep"), iVaule);
    UINT iAfterEquib = static_cast<UINT>(iVaule);

    iVaule = 2;
    params.FetchValueINT(_T("MinNt"), iVaule);
    UINT iMinNt = static_cast<UINT>(iVaule);

    iVaule = 6;
    params.FetchValueINT(_T("MaxNt"), iVaule);
    UINT iMaxNt = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("SaveStartIndex"), iVaule);
    UINT iSaveStartIndex = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("OmegaStart"), iVaule);
    UINT iOmegaStart = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("Additive"), iVaule);
    UBOOL bAdditive = 0 != iVaule;

    TArray<Real> old_polyakovs;
    params.FetchValueArrayReal(_T("Polyakovs"), old_polyakovs);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    TArray<CCString> sOldFileNames;
    TArray<Real> fOldFilePolyakov;

    for (UINT i = iMinNt; i <= iMaxNt; ++i)
    {
        CCString sFileName;
        Real fPolyakov = F(0.0);
        CCString sFileKeyName;
        CCString sPolyaKeyName;
        sFileKeyName.Format(_T("Nt%dFileName"), i);
        sPolyaKeyName.Format(_T("Nt%dPolyakov"), i);
        params.FetchStringValue(sFileKeyName, sFileName);
        params.FetchValueReal(sPolyaKeyName, fPolyakov);
        sOldFileNames.AddItem(sFileName);
        fOldFilePolyakov.AddItem(fPolyakov);

        appGeneral(_T("file: %s, |p| : %f\n"), sFileName.c_str(), fPolyakov);
    }

    Real fMaxOmega = F(0.1);
    params.FetchValueReal(_T("MaxOmega"), fMaxOmega);

    CCString sFileName;
    TCHAR buff1[256];
    TCHAR buff2[256];
    CCString sInfo;
#pragma endregion

    for (UINT uiNt = iMinNt; uiNt <= iMaxNt; ++uiNt)
    {
        TArray<INT> latticeDecomp;
        params.FetchValueArrayINT(_T("LatticeLength"), latticeDecomp);
        latticeDecomp[3] = uiNt;
        TArray<CCString> sLatticeDecomp;
        sLatticeDecomp.AddItem(appToString(latticeDecomp[0]));
        sLatticeDecomp.AddItem(appToString(latticeDecomp[1]));
        sLatticeDecomp.AddItem(appToString(latticeDecomp[2]));
        sLatticeDecomp.AddItem(appToString(latticeDecomp[3]));
        params.SetStringVectorVaule(_T("LatticeLength"), sLatticeDecomp);

        if (!appInitialCLG(params))
        {
            appCrucial(_T("Initial Failed!\n"));
            return 1;
        }

        CMeasurePolyakovXY3D* pPL = dynamic_cast<CMeasurePolyakovXY3D*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
        //TArray<TArray<CLGComplex>> polykovX_nx;
        TArray<CLGComplex> polykov;
        TArray<Real> polykovphase;
        //for (UINT uiX = 0; uiX < static_cast<UINT>(CCommonData::m_sCenter.x); ++uiX)
        //{
        //    TArray<CLGComplex> a;
        //    TArray<Real> b;
        //    polykovX_nx.AddItem(a);
        //}

        CActionGaugePlaquetteRotating* pGaugeRotation = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->GetActionById(1));
        //CActionGaugePlaquette* pGaugeNoRotation = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
        CCString sHeader;
        sHeader.Format(_T("Nt%d"), uiNt);
        appSetLogHeader(sHeader);
        appGeneral(_T("Run for Nt = %d, start baking.\n"), uiNt);

        //=============== Check oldfiles ==================
        UBOOL bNeedBake = TRUE;
        if (!bAdditive && !sOldFileNames[uiNt - iMinNt].IsEmpty())
        {
            appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sOldFileNames[uiNt - iMinNt], EFFT_CLGBin);
            pPL->OnConfigurationAccepted(_FIELDS, NULL);
            Real fError = appAbs(_cuCabsf(pPL->m_lstLoop[0]) - fOldFilePolyakov[uiNt - iMinNt]);
#if _CLG_DOUBLEFLOAT
            if (fError < F(1E-07))
#else
            if (fError < F(1E-05))
#endif
            {
                appGeneral(_T("\n ================ Bake using old file =================\n"));
                bNeedBake = FALSE;
            }
            else
            {
                appGeneral(_T("\n ================ have the initial file, but not matching.... %2.12f, %2.12f, diff=%f ===========\n"),
                    _cuCabsf(pPL->m_lstLoop[0]), fOldFilePolyakov[uiNt - iMinNt], fError);
            }
        }

        if (bAdditive)
        {
            bNeedBake = FALSE;
        }

        if (bNeedBake && iBeforeEquib > 0)
        {
            appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));
            if (NULL != pGaugeRotation)
            {
                pGaugeRotation->SetOmega(F(0.0));
            }
            else
            {
                appCrucial(_T("!!! Note you are using a no rotating action!\n"));
            }

            appGetLattice()->m_pGaugeField[0]->InitialField(EFIT_Random);

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
                appGeneral(_T("{%f, %f},\n"), _cuCabsf(pPL->m_lstLoop[i]), __cuCargf(pPL->m_lstLoop[i]));
            }
            appGeneral(_T("}\n"));
            appPopLogDate();
        }
        else
        {
            appGeneral(_T("Not Baked\n"));
        }

        UINT uiOmega = iOmegaStart;
        Real fSep = fMaxOmega / iAfterEquib;
        while (uiOmega <= iAfterEquib)
        {
            sHeader.Format(_T("Nt%dO%d"), uiNt, uiOmega);
            appSetLogHeader(sHeader);
            appGeneral(_T("\n========= Omega=%f  ==========\n"), fSep * uiOmega);

            if (NULL != pGaugeRotation)
            {
                pGaugeRotation->SetOmega(fSep * uiOmega);
            }
            else
            {
                appCrucial(_T("!!! NOTE: you are using non-rotating action!!!\n"));
            }

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
                sFileName.Format(_T("%sR_Nt%d_O%d_%d.con"), sSavePrefix.c_str(), uiNt, uiOmega, iSaveStartIndex);
                appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName, EFFT_CLGBin);
                pPL->OnConfigurationAccepted(_FIELDS, NULL);
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
                UINT uiAcce = appGetLattice()->m_pUpdator->GetConfigurationCount();
                if (uiAcce != iConfigNumberNow)
                {
                    sFileName.Format(_T("R_Nt%d_O%d_%d"), uiNt, uiOmega, uiAcce + iSaveStartIndex);
                    sFileName = sSavePrefix + sFileName;

                    //=================================
                    //Save config
                    const CCString MD5 = appGetLattice()->m_pGaugeField[0]->SaveToFile(sFileName + _T(".con"));

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

            ++uiOmega;

            appGetLattice()->m_pMeasurements->Reset();
            appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        }

        appGeneral(_T("\n========= Nt=%d finished! ==========\n\n"), uiNt);
        appPushLogDate(FALSE);
        assert(polykov.Num() == static_cast<INT>(iAfterEquib + 1));

        appGeneral(_T("|Polyakov|={\n"));
        for (UINT i = 0; i <= iAfterEquib; ++i)
        {
            appGeneral(i == iAfterEquib ? _T("%2.10f\n ") : _T("%2.10f,\n"),
                _cuCabsf(polykov[i]));
        }
        appGeneral(_T("}\n\narg(Polyakov)={\n"));

        for (UINT i = 0; i <= iAfterEquib; ++i)
        {
            appGeneral(i == iAfterEquib ? _T("%2.10f\n ") : _T("%2.10f,\n"),
                polykovphase[i]);
        }

        //for (UINT x = 0; x < static_cast<UINT>(CCommonData::m_sCenter.x); ++x)
        //{
        //    appGeneral(_T("Polyakov[x=%d]={\n"), x);
        //    for (UINT i = 0; i <= iAfterEquib; ++i)
        //    {
        //        appGeneral(i == iAfterEquib ? _T("%2.10f %s %2.10f I\n") : _T("%2.10f %s %2.10f I,\n"),
        //            polykovX_nx[x][i].x,
        //            polykovX_nx[x][i].y < F(0.0) ? _T("-") : _T("+"),
        //            appAbs(polykovX_nx[x][i].y)
        //        );
        //    }
        //    appGeneral(_T("}\n\n"));
        //}

        appPopLogDate();

        appGeneral(_T("\n=====================================\n========= Nt=%d finished! ==========\n"), uiNt);
        appQuitCLG();

    }

    return 0;
}
