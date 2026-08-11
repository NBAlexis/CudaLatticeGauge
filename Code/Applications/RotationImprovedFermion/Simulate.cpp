//=============================================================================
// FILENAME : RotationImprovedFermion.cpp
// 
// DESCRIPTION:
//
// REVISION:
// [mm/dd/yy]
// [08/22/2025 nbale]
//=============================================================================

#include "RotationImprovedFermion.h"

extern INT Simulate(CParameters& param)
{
    appSetupLog(param);

    param.Dump();

    appInitialCLG(param);

    //set omega for all field and action if exist
    DOUBLE fOmega = 0.0;
    if (param.FetchValueDOUBLE(_T("Omega"), fOmega))
    {
        //at most 5 fermion fields and 1 gauge action.
        CActionGaugePlaquetteRotating* pGR = dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->GetActionById(1));
        if (NULL != pGR)
        {
            appGeneral(_T("Set Omega for gauge action to %f\n"), fOmega);
            pGR->SetGaugeOmega(fOmega);
        }
        CFieldFermionHISQSU3R* pF1 = dynamic_cast<CFieldFermionHISQSU3R*>(appGetLattice()->GetFieldById(2));
        if (NULL != pF1)
        {
            appGeneral(_T("Set Omega for fermion1 to %f\n"), fOmega);
            pF1->SetFermionOmega(fOmega);
        }
        CFieldFermionHISQSU3R* pF2 = dynamic_cast<CFieldFermionHISQSU3R*>(appGetLattice()->GetFieldById(3));
        if (NULL != pF2)
        {
            appGeneral(_T("Set Omega for fermion2 to %f\n"), fOmega);
            pF2->SetFermionOmega(fOmega);
        }
        CFieldFermionHISQSU3R* pF3 = dynamic_cast<CFieldFermionHISQSU3R*>(appGetLattice()->GetFieldById(4));
        if (NULL != pF3)
        {
            appGeneral(_T("Set Omega for fermion3 to %f\n"), fOmega);
            pF3->SetFermionOmega(fOmega);
        }
        CFieldFermionHISQSU3R* pF4 = dynamic_cast<CFieldFermionHISQSU3R*>(appGetLattice()->GetFieldById(5));
        if (NULL != pF4)
        {
            appGeneral(_T("Set Omega for fermion4 to %f\n"), fOmega);
            pF4->SetFermionOmega(fOmega);
        }
        CFieldFermionHISQSU3R* pF5 = dynamic_cast<CFieldFermionHISQSU3R*>(appGetLattice()->GetFieldById(6));
        if (NULL != pF5)
        {
            appGeneral(_T("Set Omega for fermion5 to %f\n"), fOmega);
            pF5->SetFermionOmega(fOmega);
        }
    }

    INT iWarmUp = 5;
    param.FetchValueINT(_T("Warmup"), iWarmUp);
    INT iConfigurationNumber = 20;
    param.FetchValueINT(_T("ConfigurationNumber"), iConfigurationNumber);
    CCString sFileName = _T("Configuration");
    param.FetchStringValue(_T("FileName"), sFileName);
    INT iSaveStartIndex = 0;
    param.FetchValueINT(_T("SaveStartIndex"), iSaveStartIndex);
    CCString sSaveType = _T("EFFT_CLGBin");
    EFieldFileType eSaveType = EFFT_CLGBin;
    if (param.FetchStringValue(_T("SaveFileType"), sSaveType))
    {
        eSaveType = __STRING_TO_ENUM(EFieldFileType, sSaveType);
    }

    //warm up for 5 trajectories
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
    appGetLattice()->m_pUpdator->Update(iWarmUp, TRUE);

    //update for 10 trajectories and save
    //appClearProfiler();
    appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, sFileName, static_cast<UINT>(iSaveStartIndex), eSaveType);
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
    appGetLattice()->m_pUpdator->UpdateUntileAccept(iConfigurationNumber, TRUE);

    //appDumpProfiler();
    appQuitCLG();

    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
