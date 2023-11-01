//=============================================================================
// FILENAME : GaugeFixing.cpp
// 
// DESCRIPTION:
//  
//
// REVISION:
//  [01/17/2021 nbale]
//=============================================================================

#include "ElectricChemical.h"

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

    iVaule = 0;
    params.FetchValueINT(_T("CheckSubFolder"), iVaule);
    UBOOL bCheckSubFolder = 0 != iVaule;

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);
    appGeneral(_T("save prefix: %s\n"), sSavePrefix.c_str());

    CCString sLoadPrefix;
    params.FetchStringValue(_T("LoadPrefix"), sLoadPrefix);
    appGeneral(_T("load prefix: %s\n"), sLoadPrefix.c_str());

    CCString sSubFolderPrefix;
    params.FetchStringValue(_T("SubFolderPrefix"), sSubFolderPrefix);
    appGeneral(_T("sub folder prefix: %s\n"), sSubFolderPrefix.c_str());

    CCString sCheckSubFolderPrefix;
    params.FetchStringValue(_T("CheckSubFolderPrefix"), sCheckSubFolderPrefix);
    appGeneral(_T("checkgauge sub folder prefix: %s\n"), sCheckSubFolderPrefix.c_str());

    CCString sSaveType = _T("EFFT_CLGBin");
    EFieldFileType eSaveType = EFFT_CLGBin;
    if (params.FetchStringValue(_T("SaveType"), sSaveType))
    {
        eSaveType = __STRING_TO_ENUM(EFieldFileType, sSaveType);
    }
    appGeneral(_T("save type: %s\n"), __ENUM_TO_STRING(EFieldFileType, eSaveType).c_str());
    
    CCString sLoadType = _T("EFFT_CLGBin");
    EFieldFileType eLoadType = EFFT_CLGBin;
    if (params.FetchStringValue(_T("LoadType"), sLoadType))
    {
        eLoadType = __STRING_TO_ENUM(EFieldFileType, sLoadType);
    }
    appGeneral(_T("load type: %s\n"), __ENUM_TO_STRING(EFieldFileType, eLoadType).c_str());

    if (bOnlyCheck)
    {
        for (UINT uiOmega = iOmegaStart; uiOmega <= iOmegaEnd; ++uiOmega)
        {
            appGeneral(_T("====== Start %d : %d to %d ======\n"), uiOmega, iIndexStart, iIndexEnd);
            appPushLogDate(FALSE);
            for (UINT uiIndex = iIndexStart; uiIndex <= iIndexEnd; ++uiIndex)
            {
                CCString sSaveFile;
                if (!bCheckSubFolder)
                {
                    sSaveFile.Format(_T("%sEC_%d_%d.con"), sSavePrefix.c_str(), uiOmega, uiIndex);
                }
                else
                {
                    sSaveFile.Format(_T("%s%d/%sEC_%d_%d.con"), sCheckSubFolderPrefix.c_str(), uiOmega, sSavePrefix.c_str(), uiOmega, uiIndex);
                }
                appGetLattice()->m_pGaugeField->InitialFieldWithFile(sSaveFile, eLoadType);
#if !_CLG_DOUBLEFLOAT
                const DOUBLE fRes = appGetLattice()->m_pGaugeFixing->CheckRes(appGetLattice()->m_pGaugeField);
                if (fRes >= 0.0 && fRes < appGetLattice()->m_pGaugeFixing->m_fAccuracy)
#else
                const Real fRes = appGetLattice()->m_pGaugeFixing->CheckRes(appGetLattice()->m_pGaugeField);
                if (fRes >= F(0.0) && fRes < appGetLattice()->m_pGaugeFixing->m_fAccuracy)
#endif
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
                        appGeneral(_T("\nNot good enough %d : %d \n"), uiOmega, uiIndex);
                        appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField);
                        appGetLattice()->m_pGaugeField->SaveToFile(sSaveFile, eSaveType);
                    }
                    else
                    {
                        appGeneral(_T("\nBad %d : %d \n"), uiOmega, uiIndex);
                        CCString sLoadFile;
                        if (bSubFolder)
                        {
                            sLoadFile.Format(_T("%s%d/%sEC_%d_%d.con"), sSubFolderPrefix.c_str(), uiOmega, sLoadPrefix.c_str(), uiOmega, uiIndex);
                        }
                        else
                        {
                            sLoadFile.Format(_T("%sEC_%d_%d.con"), sLoadPrefix.c_str(), uiOmega, uiIndex);
                        }
                        appGetLattice()->m_pGaugeField->InitialFieldWithFile(sLoadFile, eLoadType);
                        appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField);
                        appGetLattice()->m_pGaugeField->SaveToFile(sSaveFile, eSaveType);
                    }
                }
                else
                {
                    appGeneral(_T("\nBad %d : %d, is %2.20f\n"), uiOmega, uiIndex, fRes);
                }
            }
            appPopLogDate();
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
                    sLoadFile.Format(_T("%s%d/%sEC_%d_%d.con"), sSubFolderPrefix.c_str(), uiOmega, sLoadPrefix.c_str(), uiOmega, uiIndex);
                }
                else
                {
                    sLoadFile.Format(_T("%sEC_%d_%d.con"), sLoadPrefix.c_str(), uiOmega, uiIndex);
                }
                sSaveFile.Format(_T("%sEC_%d_%d.con"), sSavePrefix.c_str(), uiOmega, uiIndex);
                appGeneral(_T("Fixing %d : %d \n"), uiOmega, uiIndex);
                appGetLattice()->m_pGaugeField->InitialFieldWithFile(sLoadFile, eLoadType);
                appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField);
                appGetLattice()->m_pGaugeField->SaveToFile(sSaveFile, eSaveType);
            }
        }
    }

    appGeneral(_T("\n=====================================\n========= %d finished! ==========\n"), uiNt);
    appQuitCLG();

    return 0;
}
