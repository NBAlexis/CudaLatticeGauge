//=============================================================================
// FILENAME : Simulate.cpp
// 
// DESCRIPTION:
// Agent-driven simulation: read YAML, thermalize, generate N configurations.
// No parameter scans, single run.
//
// REVISION:
//  [06/14/2026]
//=============================================================================

#include "AgentDriver.h"
#include "AgentDriver.h"

__BEGIN_NAMESPACE

INT RunAgentSimulate(CParameters& params)
{
    appSetupLog(params);

    INT iValue = 0;
    params.FetchValueINT(_T("WarmUp"), iValue);
    const UINT iWarmUp = static_cast<UINT>(iValue);

    iValue = 0;
    params.FetchValueINT(_T("SaveStartIndex"), iValue);
    const UINT iSaveStartIndex = static_cast<UINT>(iValue);

    iValue = 1;
    params.FetchValueINT(_T("ConfigurationNumber"), iValue);
    const UINT iConfigurationNumber = static_cast<UINT>(iValue);

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);

    iValue = 1;
    params.FetchValueINT(_T("SaveConfiguration"), iValue);
    const UINT iSaveConfiguration = static_cast<UINT>(iValue);

    CParameters outputs = params.GetParameter(_T("Outputs"));
    TArray<CCString> outputNames;
    TArray<CCString> outputFiles;
    outputs.FetchStringVectorValue(_T("Names"), outputNames);
    outputs.FetchStringVectorValue(_T("Files"), outputFiles);

    appGeneral(_T("[Simulate] WarmUp            = %d\n"), iWarmUp);
    appGeneral(_T("[Simulate] SaveStartIndex     = %d\n"), iSaveStartIndex);
    appGeneral(_T("[Simulate] ConfigNumber       = %d\n"), iConfigurationNumber);
    appGeneral(_T("[Simulate] SavePrefix         = %s\n"), sSavePrefix.c_str());
    appGeneral(_T("[Simulate] SaveConfiguration  = %d (0=only-last, 1=all)\n"), iSaveConfiguration);
    appGeneral(_T("[Simulate] Outputs            = %d\n"), outputNames.Num());

    appGeneral(_T("[Simulate] Initializing CLGLib...\n"));
    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }
    appGeneral(_T("[Simulate] CLGLib initialized.\n"));

    const EUpdatorType eUpdatorType = appGetLattice()->m_pUpdator->GetUpdatorType();
    const UBOOL bIsHeatbath = (EUT_Heatbath == eUpdatorType);
    appGeneral(_T("[Simulate] Updator type   = %d (Heatbath=%d)\n"),
        static_cast<INT>(eUpdatorType), static_cast<INT>(bIsHeatbath));

    // Warm-up: follow RotationImprovedFermion style.
    // Thermalize without auto-correction, then enable it for production.
    if (iWarmUp > 0)
    {
        appGeneral(_T("\n[Simulate] ========== Warm-up %d steps ==========\n"), iWarmUp);
        appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));
        appGetLattice()->m_pMeasurements->Reset();
        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
        appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
        appGetLattice()->m_pUpdator->Update(iWarmUp, TRUE);
        appGetLattice()->m_pUpdator->SetTestHdiff(FALSE);
        appGeneral(_T("[Simulate] Warm-up done.\n"));
    }
    else
    {
        appGeneral(_T("[Simulate] Skipping warm-up.\n"));
    }

    // Production: generate requested number of accepted configurations.
    appGeneral(_T("\n[Simulate] ========== Production %d configs ==========\n"), iConfigurationNumber);
    appGetLattice()->m_pMeasurements->Reset();
    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);

    if (bIsHeatbath)
    {
        // Heatbath updates every sweep, but does not use m_iAcceptedConfigurationCount.
        // We run one sweep at a time and save the field manually.
        appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));
        appGeneral(_T("[Simulate] Heatbath mode: %d sweeps will be saved manually.\n"), iConfigurationNumber);

        for (UINT i = 0; i < iConfigurationNumber; ++i)
        {
            appGetLattice()->m_pUpdator->Update(1, TRUE);
            if (1 == iSaveConfiguration || i == iConfigurationNumber - 1)
            {
                CCString sFile;
                sFile.Format(_T("%s_%d.con"), sSavePrefix.c_str(), i + iSaveStartIndex);
                appGetLattice()->m_pGaugeField[0]->SaveToFile(sFile);
                appGeneral(_T("[Simulate] Saved heatbath config %d/%d -> %s\n"),
                    i + 1, iConfigurationNumber, sFile.c_str());
            }
            else
            {
                appGeneral(_T("[Simulate] Heatbath config %d/%d accepted (not saved).\n"),
                    i + 1, iConfigurationNumber);
            }
        }
    }
    else
    {
        if (1 == iSaveConfiguration)
        {
            appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, sSavePrefix, iSaveStartIndex);
            appGeneral(_T("[Simulate] Auto-save enabled: prefix=%s start=%d\n"),
                sSavePrefix.c_str(), iSaveStartIndex);
        }
        else
        {
            appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));
            appGeneral(_T("[Simulate] Auto-save disabled: only the last config will be saved.\n"));
        }

        appGeneral(_T("[Simulate] Running until %d configurations are accepted...\n"), iConfigurationNumber);
        appGetLattice()->m_pUpdator->UpdateUntileAccept(iConfigurationNumber, TRUE);

        if (0 == iSaveConfiguration)
        {
            CCString sFile;
            sFile.Format(_T("%s_%d.con"), sSavePrefix.c_str(), iSaveStartIndex + iConfigurationNumber);
            appGetLattice()->m_pGaugeField[0]->SaveToFile(sFile);
            appGeneral(_T("[Simulate] Saved last config -> %s\n"), sFile.c_str());
        }

        appGeneral(_T("[Simulate] Production done: %d configs.\n"), iConfigurationNumber);
    }

    appGeneral(_T("\n[Simulate] Reporting measurements...\n"));
    appGetLattice()->m_pMeasurements->Report();

    if (outputNames.Num() > 0)
    {
        appGeneral(_T("\n[Simulate] Exporting %d outputs...\n"), outputNames.Num());
        for (INT i = 0; i < outputNames.Num(); ++i)
        {
            CMeasureData* pData = appGetLattice()->m_pMeasurements->GetMeasureData(outputNames[i]);
            if (NULL != pData)
            {
                appGeneral(_T("[Simulate] Export '%s' -> %s\n"), outputNames[i].c_str(), outputFiles[i].c_str());
                pData->Export(outputFiles[i]);
            }
            else
            {
                appCrucial(_T("[Simulate] Export failed: key '%s' not found.\n"), outputNames[i].c_str());
            }
        }
    }

    appGeneral(_T("[Simulate] Quitting CLGLib.\n"));
    appQuitCLG();

    appGeneral(_T("\n[Simulate] ========== Finished! Generated %d configs ==========\n\n"), iConfigurationNumber);
    return 0;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
