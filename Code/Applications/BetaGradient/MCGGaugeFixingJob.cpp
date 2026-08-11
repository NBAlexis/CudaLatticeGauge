//=============================================================================
// FILENAME : MCGGaugeFixingJob.cpp
//
// DESCRIPTION:
// MCG gauge fixing job: load configurations, apply MCG gauge fixing,
// save MCG-fixed configurations, project to Z3 center elements,
// and save Z3-projected configurations.
//
// REVISION:
//  [05/29/2026 nbale]
//=============================================================================

#include "BetaGradient.h"

INT MCGGaugeFixing(CParameters& params)
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

    iVaule = 1;
    params.FetchValueINT(_T("IndexStart"), iVaule);
    const UINT iIndexStart = static_cast<UINT>(iVaule);

    iVaule = 200;
    params.FetchValueINT(_T("IndexEnd"), iVaule);
    const UINT iIndexEnd = static_cast<UINT>(iVaule);

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

    TArray<CCString> FileList;
    params.FetchStringVectorValue(_T("FileList"), FileList);

    TArray<CCString> FolderList;
    params.FetchStringVectorValue(_T("FolderList"), FolderList);

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

    for (UINT uiOmega = iOmegaStart; uiOmega <= iOmegaEnd; ++uiOmega)
    {
        for (UINT uiIndex = iIndexStart; uiIndex <= iIndexEnd; ++uiIndex)
        {
            CCString sLoadFile;
            if (bSubFolder)
            {
                sLoadFile.Format(_T("%s/%s/%s_%s_%d.con"), sSubFolderPrefix.c_str(), FolderList[uiOmega].c_str(), sLoadPrefix.c_str(), FileList[uiOmega].c_str(), uiIndex);
            }
            else
            {
                sLoadFile.Format(_T("%s_%s_%d.con"), sLoadPrefix.c_str(), FileList[uiOmega].c_str(), uiIndex);
            }

            CCString sMCGSaveFile;
            sMCGSaveFile.Format(_T("%s_MCG_%s_%d.con"), sSavePrefix.c_str(), FileList[uiOmega].c_str(), uiIndex);

            CCString sZ3SaveFile;
            sZ3SaveFile.Format(_T("%s_Z3_%s_%d.con"), sSavePrefix.c_str(), FileList[uiOmega].c_str(), uiIndex);

            appGeneral(_T("Fixing %s : %d \n"), FileList[uiOmega].c_str(), uiIndex);
            appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sLoadFile, eLoadType);
            appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField[0]);
            appGetLattice()->m_pGaugeField[0]->SaveToFile(sMCGSaveFile, eSaveType);

            CFieldGauge* pZ3Field = dynamic_cast<CFieldGauge*>(
                appGetLattice()->GetPooledFieldById(
                    appGetLattice()->m_pGaugeField[0]->GetFieldId(),
                    _T(__FILE__), __LINE__));
            if (NULL != pZ3Field)
            {
                appGetLattice()->m_pGaugeField[0]->CopyTo(pZ3Field);

                CFieldGaugeSU3* pZ3SU3 = dynamic_cast<CFieldGaugeSU3*>(pZ3Field);
                if (NULL != pZ3SU3)
                {
                    CGaugeFixingMCGDirect::CenterProjection(pZ3SU3);
                    pZ3SU3->SaveToFile(sZ3SaveFile, eSaveType);
                }
                else
                {
                    appCrucial(_T("MCGGaugeFixing: Gauge field is not SU3, cannot perform Z3 projection\n"));
                }

                pZ3Field->Return();
            }
            else
            {
                appCrucial(_T("MCGGaugeFixing: Failed to get pooled field for Z3 projection\n"));
            }
        }
    }

    appGeneral(_T("\n=====================================\n========= MCG Gauge Fixing finished! ==========\n"));
    appQuitCLG();

    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
