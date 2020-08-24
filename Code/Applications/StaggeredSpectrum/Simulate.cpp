//=============================================================================
// FILENAME : Simulate.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [08/21/2020 nbale]
//=============================================================================

#include "StaggeredSpectrum.h"

INT SimulateStaggered(CParameters& params)
{
    appSetupLog(params);

    INT iVaule = 200;
    params.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    const UINT iBeforeEquib = static_cast<UINT>(iVaule);

    iVaule = 2500;
    params.FetchValueINT(_T("EquvibStep"), iVaule);
    const UINT iEquib = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("CompressedFile"), iVaule);
    const UBOOL bCompressedFile = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("SaveStartIndex"), iVaule);
    const UINT iSaveIndexStart = static_cast<UINT>(iVaule);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    //=========================================================
    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    //================================================================
    //Themalization
    appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));
    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
    appGetLattice()->m_pMeasurements->Reset();
    UINT uiAccepCountBeforeE = 0;
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
    {
        const UINT uiAccepCountBeforeE2 = appGetLattice()->m_pUpdator->Update(1, FALSE);
        if (uiAccepCountBeforeE != uiAccepCountBeforeE2)
        {
            uiAccepCountBeforeE = uiAccepCountBeforeE2;
        }
    }

    //================================================================
    //Start working
    UINT uiAccepCountAfterE = 0;
    CCString sFileName;
    TCHAR buff1[256];
    TCHAR buff2[256];
    CCString sInfo;

    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    appGetLattice()->m_pMeasurements->Reset();
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iEquib)
    {
        const UINT uiAccepCountBeforeE2 = appGetLattice()->m_pUpdator->Update(1, TRUE);
        if (uiAccepCountAfterE != uiAccepCountBeforeE2)
        {
            uiAccepCountAfterE = uiAccepCountBeforeE2;

            //save configurations
            sFileName.Format(_T("Matching_%d"), uiAccepCountAfterE + iSaveIndexStart);
            sFileName = sSavePrefix + sFileName;
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
            if (bCompressedFile)
            {
                appGetLattice()->m_pGaugeField->SaveToCompressedFile(sFileName + _T(".cco"));
            }
            else
            {
                appGetLattice()->m_pGaugeField->SaveToFile(sFileName + _T(".con"));
            }

            appGeneral(_T("%s saved\n"), (sFileName + _T(".con")).c_str());
        }
    }
    appGetLattice()->m_pMeasurements->Report();
    appQuitCLG();
    return 0;
}
