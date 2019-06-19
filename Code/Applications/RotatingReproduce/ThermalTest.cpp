//=============================================================================
// FILENAME : Thermal.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [06/12/2019 nbale]
//=============================================================================

#include "RotatingReproduce.h"

INT TestThermal(CParameters& params)
{
    appSetupLog(params);

    INT iVaule = 99;
    params.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    UINT iBeforeEquib = static_cast<UINT>(iVaule);

    iVaule = 6;
    params.FetchValueINT(_T("EquvibStep"), iVaule);
    UINT iEquib = static_cast<UINT>(iVaule);

    iVaule = 250;
    params.FetchValueINT(_T("OmegaSep"), iVaule);
    UINT iAfterEquib = static_cast<UINT>(iVaule);

    iVaule = 2;
    params.FetchValueINT(_T("MinNt"), iVaule);
    UINT iMinNt = static_cast<UINT>(iVaule);

    iVaule = 6;
    params.FetchValueINT(_T("MaxNt"), iVaule);
    UINT iMaxNt = static_cast<UINT>(iVaule);

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

        appGeneral(_T("f: %s, p : %f\n"), sFileName, fPolyakov);
    }
   
    Real fMaxOmega = F(0.1);
    params.FetchValueReal(_T("MaxOmega"), fMaxOmega);

    CCString sFileName;
    TCHAR buff1[256];
    TCHAR buff2[256];
    CCString sInfo;

    for (UINT uiNt = iMinNt; uiNt <= iMaxNt; ++uiNt)
    {
        TArray<INT> latticeDecomp;
        params.FetchValueArrayINT(_T("LatticeLength"), latticeDecomp);
        latticeDecomp[3] = uiNt;
        TArray<CCString> sLatticeDecomp;
        sLatticeDecomp.AddItem(appIntToString(latticeDecomp[0]));
        sLatticeDecomp.AddItem(appIntToString(latticeDecomp[1]));
        sLatticeDecomp.AddItem(appIntToString(latticeDecomp[2]));
        sLatticeDecomp.AddItem(appIntToString(latticeDecomp[3]));
        params.SetStringVectorVaule(_T("LatticeLength"), sLatticeDecomp);

        if (!appInitialCLG(params))
        {
            appCrucial(_T("Initial Failed!\n"));
            return 1;
        }

        CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
        CMeasureChiralCondensate* pCC = dynamic_cast<CMeasureChiralCondensate*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
        TArray<TArray<CLGComplex>> polykovX_nx;
        TArray<Real> polykov;
        TArray<Real> polykovphase;
        TArray<TArray<Real>> chiral_nx;
        TArray<Real> chiral;
        for (UINT uiX = 0; uiX < static_cast<UINT>(CCommonData::m_sCenter.x); ++uiX)
        {
            TArray<CLGComplex> a;
            TArray<Real> b;
            polykovX_nx.AddItem(a);
            chiral_nx.AddItem(b);
        }

        CActionGaugePlaquetteRotating * pGauageAction = dynamic_cast<CActionGaugePlaquetteRotating *>(appGetLattice()->GetActionById(1));
        CCString sHeader;
        sHeader.Format(_T("Nt%d"), uiNt);
        appSetLogHeader(sHeader);
        appGeneral(_T("Run for Nt = %d, start baking.\n"), uiNt);

        //=============== Check oldfiles ==================
        UBOOL bNeedBake = TRUE;
        if (!sOldFileNames[uiNt - iMinNt].IsEmpty())
        {
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sOldFileNames[uiNt - iMinNt], EFFT_CLGBin);
            pPL->OnConfigurationAccepted(appGetLattice()->m_pGaugeField, NULL);
            Real fError = appAbs(_cuCabsf(pPL->m_lstLoop[0]) - fOldFilePolyakov[uiNt - iMinNt]);
            if (fError < F(1E-07))
            {
                appGeneral(_T("\n ================ Bake using old file =================\n"));
                bNeedBake = FALSE;
            }
            else
            {
                appGeneral(_T("\n ================ have the initial file, but not matching.... %f ===========\n"), fError);
            }
        }
        
        if (bNeedBake)
        {
            appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));

            pGauageAction->SetOmega(F(0.0));
            appGetLattice()->m_pGaugeField->InitialField(EFIT_Random);

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
            appGeneral(_T("\n|<P>|={\n"));
            for (INT i = 0; i < pPL->m_lstLoop.Num(); ++i)
            {
                appGeneral(_T("%f,\n"), _cuCabsf(pPL->m_lstLoop[i]));
            }
            appGeneral(_T("}\n"));
            appSetLogDate(TRUE);
        }

        UINT uiOmega = 0;
        Real fSep = fMaxOmega / iAfterEquib;
        while (uiOmega <= iAfterEquib)
        {
            appGeneral(_T("\n========= Omega=%f  ==========\n"), fSep * uiOmega);
            pGauageAction->SetOmega(fSep * uiOmega);
            UINT iConfigNumberNow = 0;
            appGetLattice()->m_pMeasurements->Reset();
            appGetLattice()->m_pUpdator->SetConfigurationCount(0);

            while (iConfigNumberNow < iEquib)
            {
                appGetLattice()->m_pUpdator->Update(1, (0 == iConfigNumberNow) ? FALSE : TRUE);
                UINT uiAcce = appGetLattice()->m_pUpdator->GetConfigurationCount();
                if (uiAcce != iConfigNumberNow)
                {
                    if (0 != iConfigNumberNow)
                    {
                        sFileName.Format(_T("Rotate_Nt%d_O%d_%d"), uiNt, uiOmega, uiAcce);

                        //=================================
                        //Save info
                        appGetTimeNow(buff1, 256);
                        appGetTimeUtc(buff2, 256);
                        sInfo.Format(_T("TimeStamp : %d\nTime : %s\nTimeUTC : %s\n"),
                            appGetTimeStamp(),
                            buff1,
                            buff2);
                        sInfo = sInfo + appGetLattice()->GetInfos(_T(""));
                        appGetFileSystem()->WriteAllText(sFileName + _T(".txt"), sInfo);

                        //=================================
                        //Save config
                        appGetLattice()->m_pGaugeField->SaveToFile(sFileName + _T(".con"));
                    }

                    iConfigNumberNow = uiAcce;
                }
            }
            appGetLattice()->m_pMeasurements->Report();
            
            //===================== Polyakov loop =====================
            assert(pPL->m_lstLoop.Num() == iEquib - 1);
            assert(pPL->m_lstAverageLoopDensity.Num()
                == static_cast<INT>(CCommonData::m_sCenter.x));

            //============= polyakov gather =============
            polykov.AddItem(_cuCabsf(pPL->m_cAverageLoop));
            polykovphase.AddItem(__cuCargf(pPL->m_cAverageLoop));
            for (UINT iX = 0; iX < static_cast<UINT>(CCommonData::m_sCenter.x); ++iX)
            {
                polykovX_nx[iX].AddItem(pPL->m_lstAverageLoopDensity[iX]);
            }

            //===================== Chiral condensate =====================
            assert(pCC->m_lstCondensate.Num() == iEquib - 1);
            assert(pCC->m_lstAverageCondensateDensity.Num() == static_cast<INT>(CCommonData::m_sCenter.x));
            chiral.AddItem(_cuCabsf(pCC->m_cAverageCondensate));
            for (UINT iX = 0; iX < static_cast<UINT>(CCommonData::m_sCenter.x); ++iX)
            {
                chiral_nx[iX].AddItem(_cuCabsf(pCC->m_lstAverageCondensateDensity[iX]));
            }

            ++uiOmega;

            appGetLattice()->m_pMeasurements->Reset();
            appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        }

        appGeneral(_T("\n========= Nt=%d finished! ==========\n\n"), uiNt);
        appSetLogDate(FALSE);
        assert(polykov.Num() == static_cast<INT>(iAfterEquib + 1));

        appGeneral(_T("Polyakov={\n"));
        for (UINT i = 0; i <= iAfterEquib; ++i)
        {
            appGeneral(i == iAfterEquib ? _T("%2.10f\n ") : _T("%2.10f,\n"),
                polykov[i]);
        }
        appGeneral(_T("}\n\narg(Polyakov)={\n"));

        for (UINT i = 0; i <= iAfterEquib; ++i)
        {
            appGeneral(i == iAfterEquib ? _T("%2.10f\n ") : _T("%2.10f,\n"),
                polykovphase[i]);
        }
        appGeneral(_T("}\n\nChiralCondensate={\n"));

        for (UINT i = 0; i <= iAfterEquib; ++i)
        {
            appGeneral(i == iAfterEquib ? _T("%2.10f\n") : _T("%2.10f,\n"),
                chiral[i]);
        }

        appGeneral(_T("}\n\n"));

        for (UINT x = 0; x < static_cast<UINT>(CCommonData::m_sCenter.x); ++x)
        {
            appGeneral(_T("Polyakov[x=%d]={\n"), x);
            for (UINT i = 0; i <= iAfterEquib; ++i)
            {
                appGeneral(i == iAfterEquib ? _T("%2.10f %s %2.10f I\n") : _T("%2.10f %s %2.10f I,\n"),
                    polykovX_nx[x][i].x,
                    polykovX_nx[x][i].y < F(0.0) ? _T("-") : _T("+"),
                    appAbs(polykovX_nx[x][i].y)
                );
            }
            appGeneral(_T("}\n\n"));
        }

        for (UINT x = 0; x < static_cast<UINT>(CCommonData::m_sCenter.x); ++x)
        {
            appGeneral(_T("ChiralCondensate[x=%d]={\n"), x);
            for (UINT i = 0; i <= iAfterEquib; ++i)
            {
                appGeneral(i == iAfterEquib ? _T("%2.10f\n") : _T("%2.10f,\n"),
                    chiral_nx[x][i]);
            }
            appGeneral(_T("}\n\n"));
        }
        appSetLogDate(TRUE);

        appGeneral(_T("\n=====================================\n========= Nt=%d finished! ==========\n"), uiNt);
        appQuitCLG();

    }

    return 0;
}
