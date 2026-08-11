//=============================================================================
// FILENAME : Simulate.cpp
//
// DESCRIPTION:
// HMC simulation driver for the PSU(3) fundamental-lift action with a Z3
// 2-form boundary field B.
//
// The updator saves every accepted configuration:
//   <prefix>_<N>.con        gauge field
//   <prefix>_<N>_t<fieldId>.con   dynamic tensor2 (B) companion file
// (CUpdator::SaveConfiguration handles both.)
//
// REVISION:
//  [08/11/26]
//=============================================================================
#include "PSU3WithBoundary.h"

INT Simulate(CParameters& params)
{
    appSetupLog(params);

#pragma region read parameters
    INT iVaule = 5;
    params.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    UINT iBeforeEquib = static_cast<UINT>(iVaule);

    iVaule = 10;
    params.FetchValueINT(_T("EquvibStep"), iVaule);
    UINT iEquib = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("SaveStartIndex"), iVaule);
    UINT iSaveStartIndex = static_cast<UINT>(iVaule);

    CCString sSavePrefix = _T("PSU3");
    params.FetchStringValue(_T("SavePrefix"), sSavePrefix);

    CCString sSaveType = _T("EFFT_CLGBin");
    EFieldFileType eSaveType = EFFT_CLGBin;
    if (params.FetchStringValue(_T("SaveFileType"), sSaveType))
    {
        eSaveType = __STRING_TO_ENUM(EFieldFileType, sSaveType);
    }
#pragma endregion

    if (!appInitialCLG(params))
    {
        appCrucial(_T("PSU3WithBoundary Simulate: Initial Failed!\n"));
        return 1;
    }

    appGeneral(_T("PSU3WithBoundary Simulate: equilibration %u, production %u, save prefix %s\n"),
        iBeforeEquib, iEquib, sSavePrefix.c_str());

    // equilibration (no measurements, no saving during warm-up of the chain)
    appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));
    for (UINT i = 0; i < iBeforeEquib; ++i)
    {
        appGetLattice()->m_pUpdator->Update(1, FALSE);
    }

    // production with saving: gauge .con + dynamic tensor2 companion .con
    appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, sSavePrefix, iSaveStartIndex, eSaveType);
    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    for (UINT i = 0; i < iEquib; ++i)
    {
        appGetLattice()->m_pUpdator->Update(1, FALSE);
    }

    appQuitCLG();
    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
