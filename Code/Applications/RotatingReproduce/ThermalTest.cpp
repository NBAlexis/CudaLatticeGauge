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
    INT iVaule = 50;
    params.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    UINT iBeforeEquib = static_cast<UINT>(iVaule);
    iVaule = 200;
    params.FetchValueINT(_T("OmegaSep"), iVaule);
    UINT iAfterEquib = static_cast<UINT>(iVaule);

    iVaule = 2;
    params.FetchValueINT(_T("MinNt"), iVaule);
    UINT iMinNt = static_cast<UINT>(iVaule);

    iVaule = 8;
    params.FetchValueINT(_T("MaxNt"), iVaule);
    UINT iMaxNt = static_cast<UINT>(iVaule);

    iVaule = 5;
    params.FetchValueINT(_T("ConfigNumberEachOmega"), iVaule);
    UINT iConfigNumber = static_cast<UINT>(iVaule);
    
    Real fMaxOmega = F(0.1);
    params.FetchValueReal(_T("MaxOmega"), fMaxOmega);

    iVaule = 5;
    params.FetchValueINT(_T("MeasureListLength"), iVaule);
    params.SetStringVaule(_T("MeasureListLength"), appIntToString(iVaule + 1));

    appSetupLog(params);

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

        TArray<TArray<CLGComplex>> polykovX_nx;
        TArray<Real> polykov;
        TArray<TArray<CLGComplex>> chiral_nx;
        TArray<Real> chiral;
        for (UINT uiX = 0; uiX < _HC_Lx - 1; ++uiX)
        {
            TArray<CLGComplex> a;
            polykovX_nx.AddItem(a);
            chiral_nx.AddItem(a);
        }

        CActionGaugePlaquetteRotating * pGauageAction = dynamic_cast<CActionGaugePlaquetteRotating *>(appGetLattice()->GetActionById(1));
        appGeneral(_T("Run for Nt = %d, start baking.\n"), uiNt);
        appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));

        UINT uiLastAccept = 0;
        pGauageAction->SetOmega(F(0.0));
        appGetLattice()->m_pGaugeField->InitialField(EFIT_Random);

        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
        {
            appGetLattice()->m_pUpdator->Update(1, FALSE);
        }

        UINT uiOmega = 1;
        Real fSep = fMaxOmega / iAfterEquib;
        while (uiOmega <= iAfterEquib)
        {
            appGeneral(_T("\n========= Omega=%f  ==========\n"), fSep * uiOmega);
            pGauageAction->SetOmega(fSep * uiOmega);
            UINT iConfigNumberNow = 0;
            appGetLattice()->m_pMeasurements->Reset();
            appGetLattice()->m_pUpdator->SetConfigurationCount(0);

            while (iConfigNumberNow < iConfigNumber)
            {
                appGetLattice()->m_pUpdator->Update(1, TRUE);
                UINT uiAcce = appGetLattice()->m_pUpdator->GetConfigurationCount();
                if (uiAcce != iConfigNumberNow)
                {
                    CCString sFileName;
                    sFileName.Format(_T("Rotate_Nt%d_O%d_%d"), uiNt, uiOmega, uiAcce);

                    //=================================
                    //Save info
                    TCHAR buff1[256];
                    TCHAR buff2[256];
                    appGetTimeNow(buff1, 256);
                    appGetTimeUtc(buff2, 256);
                    CCString sInfo;
                    sInfo.Format(_T("TimeStamp : %d\nTime : %s\nTimeUTC : %s\n"),
                        appGetTimeStamp(),
                        buff1,
                        buff2);
                    sInfo = sInfo + appGetLattice()->GetInfos(_T(""));
                    appGetFileSystem()->WriteAllText(sFileName + _T(".txt"), sInfo);

                    //=================================
                    //Save config
                    appGetLattice()->m_pGaugeField->SaveToFile(sFileName + _T(".con_"));

                    iConfigNumberNow = uiAcce;
                }
            }
            appGetLattice()->m_pMeasurements->Report();

            CMeasureChargeAndCurrents* pCaC = dynamic_cast<CMeasureChargeAndCurrents*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
            CMeasureChiralCondensate* pCC = dynamic_cast<CMeasureChiralCondensate*>(appGetLattice()->m_pMeasurements->GetMeasureById(6));
            CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(5));
            //===================== Polyakov loop =====================
            assert(pPL->m_lstLoop.Num() == iConfigNumber);
            assert(pPL->m_lstLoopDensity.Num() == iConfigNumber * (_HC_Lx - 1) * (_HC_Ly - 1));
            for (UINT iconf = 0; iconf < iConfigNumber; ++iconf)
            {
                polykov.AddItem(_cuCabsf(pPL->m_lstLoop[iconf]));
                for (UINT iX = 0; iX < _HC_Lx - 1; ++iX)
                {
                    polykovX_nx[iX].AddItem(pPL->m_lstLoopDensity[
                        iconf * (_HC_Lx - 1) * (_HC_Ly - 1)
                            + CCommonData::m_sCenter.y * (_HC_Lx - 1)
                            + iX
                    ]);
                }
            }

            //===================== Chiral condensate =====================
            assert(pCC->m_lstCondensate.Num() == iConfigNumber);
            assert(pCaC->m_lstAllRes.Num() == iConfigNumber * (_HC_Lx - 1) * CMeasureChargeAndCurrents::_kGammaInInterests);
            for (UINT iconf = 0; iconf < iConfigNumber; ++iconf)
            {
                chiral.AddItem(_cuCabsf(pCC->m_lstCondensate[iconf]));
                for (UINT iX = 0; iX < _HC_Lx - 1; ++iX)
                {
                    chiral_nx[iX].AddItem(pCaC->m_lstAllRes[
                        iconf * (_HC_Lx - 1) * CMeasureChargeAndCurrents::_kGammaInInterests
                            + iX
                    ]);
                }
            }

            ++uiOmega;

            appGetLattice()->m_pMeasurements->Reset();
            appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        }

        appGeneral(_T("\n========= Nt=%d finished! ==========\n\n"), uiNt);

        assert(polykov.Num() == iConfigNumber * iAfterEquib);

        appGeneral(_T("Polyakov={\n"));
        for (UINT i = 0; i < iAfterEquib; ++i)
        {
            appGeneral(_T("{"));
            for (UINT j = 0; j < iConfigNumber; ++j)
            {
                appGeneral(j == (iConfigNumber - 1) ? _T("%2.10f ") : _T("%2.10f, "), polykov[i * iAfterEquib + j]);
            }
            appGeneral(i == (iAfterEquib - 1) ? _T("}\n") : _T("},\n"));
        }
        appGeneral(_T("}\n\nChiralCondensate={\n"));

        for (UINT i = 0; i < iAfterEquib; ++i)
        {
            appGeneral(_T("{"));
            for (UINT j = 0; j < iConfigNumber; ++j)
            {
                appGeneral(j == (iConfigNumber - 1) ? _T("%2.10f ") : _T("%2.10f, "), chiral[i * iAfterEquib + j]);
            }
            appGeneral(i == (iAfterEquib - 1) ? _T("}\n") : _T("},\n"));
        }

        appGeneral(_T("}\n\n"));

        for (UINT x = 0; x < _HC_Lx - 1; ++x)
        {
            appGeneral(_T("Polyakov[x=%d]={\n"), x);
            for (UINT i = 0; i < iAfterEquib; ++i)
            {
                appGeneral(_T("{"));
                for (UINT j = 0; j < iConfigNumber; ++j)
                {
                    appGeneral(j == (iConfigNumber - 1) ? _T("%2.10f %s %2.10f I") : _T("%2.10f %s %2.10f I, "), 
                        polykovX_nx[x][i * iAfterEquib + j].x,
                        polykovX_nx[x][i * iAfterEquib + j].y < F(0.0) ? _T("-") : _T("+"),
                        appAbs(polykovX_nx[x][i * iAfterEquib + j].y)
                    );
                }
                appGeneral(i == (iAfterEquib - 1) ? _T("}\n") : _T("},\n"));
            }
            appGeneral(_T("}\n\n"));
        }

        for (UINT x = 0; x < _HC_Lx - 1; ++x)
        {
            appGeneral(_T("ChiralCondensate[x=%d]={\n"), x);
            for (UINT i = 0; i < iAfterEquib; ++i)
            {
                appGeneral(_T("{"));
                for (UINT j = 0; j < iConfigNumber; ++j)
                {
                    appGeneral(j == (iConfigNumber - 1) ? _T("%2.10f %s %2.10f I") : _T("%2.10f %s %2.10f I, "),
                        chiral_nx[x][i * iAfterEquib + j].x,
                        chiral_nx[x][i * iAfterEquib + j].y < F(0.0) ? _T("-") : _T("+"),
                        appAbs(chiral_nx[x][i * iAfterEquib + j].y)
                    );
                }
                appGeneral(i == (iAfterEquib - 1) ? _T("}\n") : _T("},\n"));
            }
            appGeneral(_T("}\n\n"));
        }

        appGeneral(_T("\n=====================================\n========= Nt=%d finished! ==========\n"), uiNt);
        appQuitCLG();
    }
}
