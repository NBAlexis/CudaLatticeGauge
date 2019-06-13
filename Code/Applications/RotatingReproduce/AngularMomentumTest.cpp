//=============================================================================
// FILENAME : AngularMomentumTest.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [06/12/2019 nbale]
//=============================================================================

#include "RotatingReproduce.h"

INT TestAngularMomentum(CParameters& params)
{
    INT iVaule = 5;
    params.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    UINT iBeforeEquib = static_cast<UINT>(iVaule);
    iVaule = 200;
    params.FetchValueINT(_T("EquvibStep"), iVaule);
    UINT iAfterEquib = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("LoopStart"), iVaule);
    UINT iLoopStart = static_cast<UINT>(iVaule);

    iVaule = 6;
    params.FetchValueINT(_T("LoopEnd"), iVaule);
    UINT iLoopEnd = static_cast<UINT>(iVaule);

    appSetupLog(params);
    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    CActionGaugePlaquetteRotating * pGauageAction = dynamic_cast<CActionGaugePlaquetteRotating *>(appGetLattice()->GetActionById(1));
    appGeneral(_T("Start up.\n"));
    appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));

    for (UINT i = iLoopStart; i < iLoopEnd; ++i)
    {
        UINT uiLastAccept = 0;
        Real omega = F(0.02) * i;
        appGeneral(_T("%d run: Omega = %f.\n"), i + 1, omega);
        pGauageAction->SetOmega(omega);
        appGetLattice()->m_pGaugeField->InitialField(EFIT_Random);

        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
        {
            appGetLattice()->m_pUpdator->Update(1, FALSE);
        }

        appGetLattice()->m_pMeasurements->Reset();
        appGetLattice()->m_pUpdator->SetConfigurationCount(0);

        while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iAfterEquib)
        {
            appGetLattice()->m_pUpdator->Update(1, TRUE);
            UINT uiAcce = appGetLattice()->m_pUpdator->GetConfigurationCount();
            if (uiAcce != uiLastAccept)
            {
                CCString sFileName;
                sFileName.Format(_T("Rotate_00%d_%d"), i * 2, uiAcce);

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
                appGetLattice()->m_pGaugeField->SaveToFile(sFileName + _T(".con"));


                uiLastAccept = uiAcce;
            }
        }
        appGetLattice()->m_pMeasurements->Report();
    }

    appGeneral(_T("\n========= all finished! ==========\n"));
    appQuitCLG();
}