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

    iVaule = 8;
    params.FetchValueINT(_T("MaxNt"), iVaule);
    UINT iMaxNt = static_cast<UINT>(iVaule);

    iVaule = 5;
    params.FetchValueINT(_T("ConfigNumberEachOmega"), iVaule);
    UINT iConfigNumber = static_cast<UINT>(iVaule);
    
    Real fMaxOmega = F(0.1);
    params.FetchValueReal(_T("MaxOmega"), fMaxOmega);


    appSetupLog(params);

    for (UINT uiNt = 2; uiNt <= iMaxNt; ++uiNt)
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
            ++uiOmega;

            appGetLattice()->m_pMeasurements->Reset();
            appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        }

        appGeneral(_T("\n========= Nt=%d finished! ==========\n"), uiNt);
        appQuitCLG();
    }
}
