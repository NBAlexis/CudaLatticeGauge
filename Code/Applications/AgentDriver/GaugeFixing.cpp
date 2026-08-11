//=============================================================================
// FILENAME : GaugeFixing.cpp
// 
// DESCRIPTION:
// Agent-driven gauge fixing: load, fix, save over a range of configurations.
// No parameter scans, single run.
//
// REVISION:
//  [06/14/2026]
//=============================================================================

#include "AgentDriver.h"
#include "AgentDriver.h"

__BEGIN_NAMESPACE

INT RunAgentGaugeFixing(CParameters& params)
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

    CCString sSavePrefix;
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);

    CCString sLoadType = _T("EFFT_CLGBin");
    EFieldFileType eLoadType = EFFT_CLGBin;
    if (params.FetchStringValue(_T("LoadType"), sLoadType))
    {
        eLoadType = __STRING_TO_ENUM(EFieldFileType, sLoadType);
    }

    CCString sSaveType = _T("EFFT_CLGBin");
    EFieldFileType eSaveType = EFFT_CLGBin;
    if (params.FetchStringValue(_T("SaveType"), sSaveType))
    {
        eSaveType = __STRING_TO_ENUM(EFieldFileType, sSaveType);
    }

    iValue = 0;
    params.FetchValueINT(_T("OnlyCheck"), iValue);
    const UBOOL bOnlyCheck = (0 != iValue);

    iValue = 0;
    params.FetchValueINT(_T("CheckAndFix"), iValue);
    const UBOOL bCheckAndFix = (0 != iValue);

    appGeneral(_T("[GaugeFixing] StartN      = %d\n"), iStartN);
    appGeneral(_T("[GaugeFixing] EndN        = %d\n"), iEndN);
    appGeneral(_T("[GaugeFixing] LoadPrefix  = %s\n"), sLoadPrefix.c_str());
    appGeneral(_T("[GaugeFixing] SavePrefix  = %s\n"), sSavePrefix.c_str());
    appGeneral(_T("[GaugeFixing] LoadType    = %s\n"), sLoadType.c_str());
    appGeneral(_T("[GaugeFixing] SaveType    = %s\n"), sSaveType.c_str());
    appGeneral(_T("[GaugeFixing] OnlyCheck   = %d\n"), static_cast<INT>(bOnlyCheck));
    appGeneral(_T("[GaugeFixing] CheckAndFix = %d\n"), static_cast<INT>(bCheckAndFix));

    appGeneral(_T("[GaugeFixing] Initializing CLGLib...\n"));
    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }
    appGeneral(_T("[GaugeFixing] CLGLib initialized.\n"));

    appGeneral(_T("\n[GaugeFixing] ========== Gauge Fixing: %d to %d ==========\n"), iStartN, iEndN);
    appPushLogDate(FALSE);

    UINT uiFixed = 0;
    UINT uiGood = 0;
    UINT uiBad = 0;
    for (UINT uiN = iStartN; uiN <= iEndN; ++uiN)
    {
        CCString sLoadFile;
        sLoadFile.Format(_T("%s_%d.con"), sLoadPrefix.c_str(), uiN);
        appGeneral(_T("[GaugeFixing] %d: loading %s\n"), uiN, sLoadFile.c_str());
        appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sLoadFile, eLoadType);

        if (bOnlyCheck)
        {
#if !_CLG_DOUBLEFLOAT
            const DOUBLE fRes = appGetLattice()->m_pGaugeFixing->CheckRes(appGetLattice()->m_pGaugeField[0]);
            const UBOOL bGood = (fRes >= 0.0 && fRes < appGetLattice()->m_pGaugeFixing->m_fAccuracy);
#else
            const Real fRes = appGetLattice()->m_pGaugeFixing->CheckRes(appGetLattice()->m_pGaugeField[0]);
            const UBOOL bGood = (fRes >= F(0.0) && fRes < appGetLattice()->m_pGaugeFixing->m_fAccuracy);
#endif
            if (bGood)
            {
                appGeneral(_T("[GaugeFixing] %d: good (%f)\n"), uiN, static_cast<Real>(fRes));
                ++uiGood;
                continue;
            }

            if (bCheckAndFix)
            {
                appGeneral(_T("[GaugeFixing] %d: fixing (%f)\n"), uiN, static_cast<Real>(fRes));
                appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField[0]);
                CCString sSaveFile;
                sSaveFile.Format(_T("%s_%d.con"), sSavePrefix.c_str(), uiN);
                appGeneral(_T("[GaugeFixing] %d: saving %s\n"), uiN, sSaveFile.c_str());
                appGetLattice()->m_pGaugeField[0]->SaveToFile(sSaveFile, eSaveType);
                ++uiFixed;
            }
            else
            {
                appGeneral(_T("[GaugeFixing] %d: bad (%f)\n"), uiN, static_cast<Real>(fRes));
                ++uiBad;
            }
        }
        else
        {
            appGeneral(_T("[GaugeFixing] %d: running gauge fixing\n"), uiN);
            appGetLattice()->m_pGaugeFixing->GaugeFixing(appGetLattice()->m_pGaugeField[0]);
            CCString sSaveFile;
            sSaveFile.Format(_T("%s_%d.con"), sSavePrefix.c_str(), uiN);
            appGeneral(_T("[GaugeFixing] %d: saving %s\n"), uiN, sSaveFile.c_str());
            appGetLattice()->m_pGaugeField[0]->SaveToFile(sSaveFile, eSaveType);
            ++uiFixed;
        }
    }

    appPopLogDate();
    appGeneral(_T("\n[GaugeFixing] Summary: fixed=%d good=%d bad=%d\n"), uiFixed, uiGood, uiBad);
    appGeneral(_T("[GaugeFixing] Quitting CLGLib.\n"));
    appQuitCLG();

    appGeneral(_T("\n[GaugeFixing] ========== Finished! ==========\n\n"));
    return 0;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
