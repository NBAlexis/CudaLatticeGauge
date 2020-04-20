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
    const UINT iOmegaStart = static_cast<UINT>(iVaule);

    iVaule = 20;
    params.FetchValueINT(_T("OmegaEnd"), iVaule);
    const UINT iOmegaEnd = static_cast<UINT>(iVaule);

    iVaule = 3;
    params.FetchValueINT(_T("Nt"), iVaule);
    const UINT uiNt = static_cast<UINT>(iVaule);

    iVaule = 1;
    params.FetchValueINT(_T("IndexStart"), iVaule);
    const UINT iIndexStart = static_cast<UINT>(iVaule);

    iVaule = 200;
    params.FetchValueINT(_T("IndexEnd"), iVaule);
    const UINT iIndexEnd = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("OnlyCheck"), iVaule);
    const UBOOL bOnlyCheck = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("CheckAndFix"), iVaule);
    const UBOOL bAlsoFix = 0 != iVaule;

    iVaule = 0;
    params.FetchValueINT(_T("SubFolder"), iVaule);
    UBOOL bSubFolder = 0 != iVaule;

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    CCString sLoadPrefix;
    params.FetchStringValue(_T("LoadPrefix"), sLoadPrefix);
    appGeneral(_T("load prefix: %s\n"), sLoadPrefix.c_str());

    CCString sSubFolderPrefix;
    params.FetchStringValue(_T("SubFolderPrefix"), sSubFolderPrefix);
    appGeneral(_T("sub folder prefix: %s\n"), sSubFolderPrefix.c_str());

    if (bOnlyCheck)
    {
        for (UINT uiOmega = iOmegaStart; uiOmega <= iOmegaEnd; ++uiOmega)
        {
            appGeneral(_T("====== Start O%d : %d to %d ======\n"), uiOmega, iIndexStart, iIndexEnd);
            appSetLogDate(FALSE);
            for (UINT uiIndex = iIndexStart; uiIndex <= iIndexEnd; ++uiIndex)
            {
                CCString sSaveFile;
                sSaveFile.Format(_T("%sRotate_Nt%d_O%d_%d.con"), sSavePrefix.c_str(), uiNt, uiOmega, uiIndex);
                appGetLattice()->m_pGaugeField->InitialFieldWithFile(sSaveFile, EFFT_CLGBin);

                const Real fRes = appGetLattice()->m_pGaugeFixing->CheckRes(appGetLattice()->m_pGaugeField);
                if (fRes >= F(0.0) && fRes < appGetLattice()->m_pGaugeFixing->m_fAccuracy)
                {
                    if (0 == (uiIndex % 20))
                    {
                        appGeneral(_T("%1.2f,\n"), fRes);
                    }
                    else
                    {
                        appGeneral(_T("%1.2f,"), fRes);
                    }

                    continue;
                }

                if (bAlsoFix)
                {
                    if (fRes >= F(0.0) && fRes < F(0.01))
                    {
                        appGeneral(_T("\nNot good enough O%d : %d \n"), uiOmega, uiIndex);
                        appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField);
                        appGetLattice()->m_pGaugeField->SaveToFile(sSaveFile);
                    }
                    else
                    {
                        appGeneral(_T("\nBad O%d : %d \n"), uiOmega, uiIndex);
                        CCString sLoadFile;
                        if (bSubFolder)
                        {
                            sLoadFile.Format(_T("%s/O%d/%sRotate_Nt%d_O%d_%d.con"), sSubFolderPrefix.c_str(), uiOmega, sLoadPrefix.c_str(), _HC_Lt, uiOmega, uiIndex);
                        }
                        else
                        {
                            sLoadFile.Format(_T("%sRotate_Nt%d_O%d_%d.con"), sLoadPrefix.c_str(), uiNt, uiOmega, uiIndex);
                        }
                        appGetLattice()->m_pGaugeField->InitialFieldWithFile(sLoadFile, EFFT_CLGBin);
                        appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField);
                        appGetLattice()->m_pGaugeField->SaveToFile(sSaveFile);
                    }
                }
                else
                {
                    appGeneral(_T("\nBad O%d : %d, is %2.20f\n"), uiOmega, uiIndex, fRes);
                }
            }
            appSetLogDate(TRUE);
        }
    }
    else
    {
        for (UINT uiOmega = iOmegaStart; uiOmega <= iOmegaEnd; ++uiOmega)
        {
            for (UINT uiIndex = iIndexStart; uiIndex <= iIndexEnd; ++uiIndex)
            {
                CCString sLoadFile;
                CCString sSaveFile;
                if (bSubFolder)
                {
                    sLoadFile.Format(_T("%s/O%d/%sRotate_Nt%d_O%d_%d.con"), sSubFolderPrefix.c_str(), uiOmega, sLoadPrefix.c_str(), _HC_Lt, uiOmega, uiIndex);
                }
                else
                {
                    sLoadFile.Format(_T("%sRotate_Nt%d_O%d_%d.con"), sLoadPrefix.c_str(), uiNt, uiOmega, uiIndex);
                }
                sSaveFile.Format(_T("%sRotate_Nt%d_O%d_%d.con"), sSavePrefix.c_str(), uiNt, uiOmega, uiIndex);
                appGeneral(_T("Fixing O%d : %d \n"), uiOmega, uiIndex);
                appGetLattice()->m_pGaugeField->InitialFieldWithFile(sLoadFile, EFFT_CLGBin);
                appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField);
                appGetLattice()->m_pGaugeField->SaveToFile(sSaveFile);
            }
        }
    }

    appGeneral(_T("\n=====================================\n========= finished! ==========\n"), uiNt);
    appQuitCLG();

    return 0;
}
