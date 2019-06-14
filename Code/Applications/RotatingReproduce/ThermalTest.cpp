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
    INT iVaule = 49;
    params.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    UINT iBeforeEquib = static_cast<UINT>(iVaule);
    iVaule = 250;
    params.FetchValueINT(_T("OmegaSep"), iVaule);
    UINT iAfterEquib = static_cast<UINT>(iVaule);
    iVaule = 2;
    params.FetchValueINT(_T("MinNt"), iVaule);
    UINT iMinNt = static_cast<UINT>(iVaule);
    iVaule = 6;
    params.FetchValueINT(_T("MaxNt"), iVaule);
    UINT iMaxNt = static_cast<UINT>(iVaule);
   
    Real fMaxOmega = F(0.1);
    params.FetchValueReal(_T("MaxOmega"), fMaxOmega);

    appSetupLog(params);

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
        CCString sHeader;
        sHeader.Format(_T("Nt%d"), uiNt);
        appSetLogHeader(sHeader);
        appGeneral(_T("Run for Nt = %d, start baking.\n"), uiNt);
        
        appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));

        //UINT uiLastAccept = 0;
        pGauageAction->SetOmega(F(0.0));
        appGetLattice()->m_pGaugeField->InitialField(EFIT_Random);

        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
        {
            appGetLattice()->m_pUpdator->Update(1, FALSE);
        }

        //=================================
        //Save config
        appGetLattice()->m_pGaugeField->SaveToFile(sFileName + _T(".con"));

        UINT uiOmega = 0;
        Real fSep = fMaxOmega / iAfterEquib;
        while (uiOmega <= iAfterEquib)
        {
            appGeneral(_T("\n========= Omega=%f  ==========\n"), fSep * uiOmega);
            pGauageAction->SetOmega(fSep * uiOmega);
            UINT iConfigNumberNow = 0;
            appGetLattice()->m_pMeasurements->Reset();
            appGetLattice()->m_pUpdator->SetConfigurationCount(0);

            while (iConfigNumberNow < 1)
            {
                appGetLattice()->m_pUpdator->Update(1, TRUE);
                UINT uiAcce = appGetLattice()->m_pUpdator->GetConfigurationCount();
                if (uiAcce != iConfigNumberNow)
                {
                    sFileName.Format(_T("Rotate_Nt%d_%d"), uiNt, uiOmega);

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

                    iConfigNumberNow = uiAcce;
                }
            }
            appGetLattice()->m_pMeasurements->Report();

            CMeasureChargeAndCurrents* pCaC = dynamic_cast<CMeasureChargeAndCurrents*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
            CMeasureChiralCondensate* pCC = dynamic_cast<CMeasureChiralCondensate*>(appGetLattice()->m_pMeasurements->GetMeasureById(3));
            CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
            //===================== Polyakov loop =====================
            assert(pPL->m_lstLoop.Num() == 1);
            assert(pPL->m_lstLoopDensity.Num() == static_cast<INT>((_HC_Lx - 1) * (_HC_Ly - 1)));

            //============= polyakov gather =============
            polykov.AddItem(_cuCabsf(pPL->m_lstLoop[0]));
            for (UINT iX = 0; iX < _HC_Lx - 1; ++iX)
            {
                polykovX_nx[iX].AddItem(pPL->m_lstLoopDensity[
                    //iconf * (_HC_Lx - 1) * (_HC_Ly - 1) +
                        CCommonData::m_sCenter.y * (_HC_Lx - 1)
                        + iX
                ]);
            }

            //===================== Chiral condensate =====================
            assert(pCC->m_lstCondensate.Num() == 1);
            assert(pCaC->m_lstAllRes.Num() == static_cast<INT>((_HC_Lx - 1) * CMeasureChargeAndCurrents::_kGammaInInterests));
            chiral.AddItem(_cuCabsf(pCC->m_lstCondensate[0]));
            for (UINT iX = 0; iX < _HC_Lx - 1; ++iX)
            {
                chiral_nx[iX].AddItem(pCaC->m_lstAllRes[
                    //iconf * (_HC_Lx - 1) * CMeasureChargeAndCurrents::_kGammaInInterests
                        + iX
                ]);
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
        appGeneral(_T("}\n\nChiralCondensate={\n"));

        for (UINT i = 0; i <= iAfterEquib; ++i)
        {
            appGeneral(i == iAfterEquib ? _T("%2.10f\n") : _T("%2.10f,\n"),
                chiral[i]);
        }

        appGeneral(_T("}\n\n"));

        for (UINT x = 0; x < _HC_Lx - 1; ++x)
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

        for (UINT x = 0; x < _HC_Lx - 1; ++x)
        {
            appGeneral(_T("ChiralCondensate[x=%d]={\n"), x);
            for (UINT i = 0; i <= iAfterEquib; ++i)
            {
                appGeneral(i == iAfterEquib ? _T("%2.10f %s %2.10f I\n") : _T("%2.10f %s %2.10f I,\n"),
                    chiral_nx[x][i].x,
                    chiral_nx[x][i].y < F(0.0) ? _T("-") : _T("+"),
                    appAbs(chiral_nx[x][i].y)
                );
            }
            appGeneral(_T("}\n\n"));
        }
        appSetLogDate(TRUE);

        appGeneral(_T("\n=====================================\n========= Nt=%d finished! ==========\n"), uiNt);
        appQuitCLG();

    }

    return 0;
}
