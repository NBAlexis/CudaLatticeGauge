//=============================================================================
// FILENAME : Measure.cpp
// 
// DESCRIPTION:
// Agent-driven measurement: load configurations, run all registered measures,
// and export selected keys using CMeasurementManager::GetMeasureData.
//
// REVISION:
//  [06/14/2026]
//=============================================================================

#include "AgentDriver.h"
#include "AgentDriver.h"

__BEGIN_NAMESPACE

static void __ExportKey(const CCString& sKey, const CCString& sFile)
{
    appGeneral(_T("[Measure] Export key '%s' -> %s\n"), sKey.c_str(), sFile.c_str());
    CMeasureData* pData = appGetLattice()->m_pMeasurements->GetMeasureData(sKey);
    if (NULL != pData)
    {
        pData->Export(sFile);
        appGeneral(_T("[Measure] Exported '%s' OK.\n"), sKey.c_str());
    }
    else
    {
        appCrucial(_T("[Measure] Export failed: key '%s' not found.\n"), sKey.c_str());
    }
}

INT RunAgentMeasure(CParameters& params)
{
    appSetupLog(params);

    INT iValue = 0;
    params.FetchValueINT(_T("StartN"), iValue);
    const UINT iStartN = static_cast<UINT>(iValue);

    iValue = 0;
    params.FetchValueINT(_T("EndN"), iValue);
    const UINT iEndN = static_cast<UINT>(iValue);

    CCString sLoadPrefix;
    params.FetchStringValue(_T("LoadPrefix"), sLoadPrefix);

    CCString sLoadType = _T("EFFT_CLGBin");
    EFieldFileType eLoadType = EFFT_CLGBin;
    if (params.FetchStringValue(_T("LoadType"), sLoadType))
    {
        eLoadType = __STRING_TO_ENUM(EFieldFileType, sLoadType);
    }

    CParameters outputs = params.GetParameter(_T("Outputs"));
    TArray<CCString> names;
    outputs.FetchStringVectorValue(_T("Names"), names);
    TArray<CCString> files;
    outputs.FetchStringVectorValue(_T("Files"), files);

    appGeneral(_T("[Measure] StartN      = %d\n"), iStartN);
    appGeneral(_T("[Measure] EndN        = %d\n"), iEndN);
    appGeneral(_T("[Measure] LoadPrefix  = %s\n"), sLoadPrefix.c_str());
    appGeneral(_T("[Measure] LoadType    = %s\n"), sLoadType.c_str());
    appGeneral(_T("[Measure] Output keys = %d\n"), names.Num());
    for (INT i = 0; i < names.Num(); ++i)
    {
        appGeneral(_T("[Measure]   [%d] %s -> %s\n"), i, names[i].c_str(), files[i].c_str());
    }

    if (names.Num() != files.Num())
    {
        appCrucial(_T("Outputs.Names and Outputs.Files must have the same length.\n"));
        return 1;
    }

    appGeneral(_T("[Measure] Initializing CLGLib...\n"));
    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }
    appGeneral(_T("[Measure] CLGLib initialized.\n"));

    appGeneral(_T("[Measure] Resetting measurement manager.\n"));
    appGetLattice()->m_pMeasurements->Reset();

    appGeneral(_T("\n[Measure] ========== Measure: %d to %d ==========\n"), iStartN, iEndN);
    appPushLogDate(FALSE);

    for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
    {
        CCString sLoadFile;
        sLoadFile.Format(_T("%s_%d.con"), sLoadPrefix.c_str(), uiN);
        appGeneral(_T("[Measure] %d: loading %s\n"), uiN, sLoadFile.c_str());
        appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sLoadFile, eLoadType);

        appGeneral(_T("[Measure] %d: fixing field boundary\n"), uiN);
        appGetLattice()->FixAllFieldBoundary();

        appGeneral(_T("[Measure] %d: running measurements\n"), uiN);
        appGetLattice()->OnUpdatorConfigurationAccepted(
            appGetLattice()->m_pGaugeField.Num(),
            appGetLattice()->m_pBosonField.Num(),
            appGetLattice()->m_pTensor2Field.Num(),
            appGetLattice()->m_pGaugeField.GetData(),
            appGetLattice()->m_pBosonField.GetData(),
            appGetLattice()->m_pTensor2Field.GetData(),
            NULL);

        appGeneral(_T("[Measure] %d: done\n"), uiN);
    }

    appPopLogDate();

    appGeneral(_T("\n[Measure] Reporting measurements...\n"));
    appGetLattice()->m_pMeasurements->Report();

    appGeneral(_T("\n[Measure] Exporting %d output keys...\n"), names.Num());
    for (INT i = 0; i < names.Num(); ++i)
    {
        __ExportKey(names[i], files[i]);
    }

    appGeneral(_T("[Measure] Quitting CLGLib.\n"));
    appQuitCLG();

    appGeneral(_T("\n[Measure] ========== Finished! ==========\n\n"));
    return 0;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
