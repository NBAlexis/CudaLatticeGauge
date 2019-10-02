//=============================================================================
// FILENAME : Thermal.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [06/12/2019 nbale]
//=============================================================================

#include "RotatingReproduce.h"

INT GaugeFixing(CParameters& params)
{
    appSetupLog(params);
    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    INT iVaule = 0;
    params.FetchValueINT(_T("OmegaStart"), iVaule);
    UINT iOmegaStart = static_cast<UINT>(iVaule);

    iVaule = 20;
    params.FetchValueINT(_T("OmegaEnd"), iVaule);
    UINT iOmegaEnd = static_cast<UINT>(iVaule);

    iVaule = 3;
    params.FetchValueINT(_T("Nt"), iVaule);
    UINT uiNt = static_cast<UINT>(iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("IndexStart"), iVaule);
    UINT iIndexStart = static_cast<UINT>(iVaule);

    iVaule = 200;
    params.FetchValueINT(_T("IndexEnd"), iVaule);
    UINT iIndexEnd = static_cast<UINT>(iVaule);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    CCString sLoadPrefix;
    params.FetchStringValue(_T("LoadPrefix"), sLoadPrefix);
    appGeneral(_T("load prefix: %s\n"), sLoadPrefix.c_str());

    for (UINT uiOmega = iOmegaStart; uiOmega <= iOmegaEnd; ++uiOmega)
    {
        for (UINT uiIndex = iIndexStart; uiIndex <= iIndexEnd; ++uiIndex)
        {
            CCString sLoadFile;
            CCString sSaveFile;
            sLoadFile.Format(_T("%sRotate_Nt%d_O%d_%d.con"), sLoadPrefix.c_str(), uiNt, uiOmega, uiIndex);
            sSaveFile.Format(_T("%sRotate_Nt%d_O%d_%d.con"), sSavePrefix.c_str(), uiNt, uiOmega, uiIndex);
            appGeneral(_T("Fixing O%d : %d \n"), uiOmega, uiIndex);
            appGetLattice()->m_pGaugeField->InitialFieldWithFile(sLoadFile, EFFT_CLGBin);
            appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField);
            appGetLattice()->m_pGaugeField->SaveToFile(sSaveFile);
        }
    }

    appGeneral(_T("\n=====================================\n========= finished! ==========\n"), uiNt);
    appQuitCLG();

    return 0;
}
