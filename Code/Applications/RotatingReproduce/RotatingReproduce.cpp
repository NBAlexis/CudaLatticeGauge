//=============================================================================
// FILENAME : RotatingReproduce.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [02/26/2019 nbale]
//=============================================================================

#include "RotatingReproduce.h"

int main(int argc, char * argv[])
{
    CParameters params;
#if _CLG_DEBUG
    CYAMLParser::ParseFile(_T("RotatingReproduce.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/RotatingReproduce.yaml"), params);
#endif

    INT iVaule = 5;
    params.FetchValueINT(_T("BeforeEquvibStep"), iVaule);
    UINT iBeforeEquib = static_cast<UINT>(iVaule);
    iVaule = 200;
    params.FetchValueINT(_T("EquvibStep"), iVaule);
    UINT iAfterEquib = static_cast<UINT>(iVaule);

    iVaule = 0;
    params.FetchValueINT(_T("LoopStart"), iVaule);
    UINT iLoopStart = static_cast<UINT>(iVaule);

    iVaule = 6;
    params.FetchValueINT(_T("LoopEnd"), iVaule);
    UINT iLoopEnd = static_cast<UINT>(iVaule);

    appSetupLog(params);
    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    CActionGaugePlaquetteRotating * pGauageAction = dynamic_cast<CActionGaugePlaquetteRotating *>(appGetLattice()->GetActionById(1));    
    appGeneral(_T("Start up.\n"));
    

    for (UINT i = iLoopStart; i < iLoopEnd; ++i)
    {
        Real omega = F(0.02) * i;
        appGeneral(_T("%d run: Omega = %f.\n"), i + 1, omega);
        pGauageAction->SetOmega(omega);
        appGetLattice()->m_pGaugeField->InitialField(EFIT_Random);

        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        appGetLattice()->m_pUpdator->SetSaveConfiguration(FALSE, _T("notsave"));
        while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
        {
            appGetLattice()->m_pUpdator->Update(1, FALSE);
        }

        appGetLattice()->m_pMeasurements->Reset();
        appGetLattice()->m_pUpdator->SetConfigurationCount(0);
        CCString sFileName;
        sFileName.Format(_T("Rotate_00%d"), i);
        appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, sFileName);
        while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iAfterEquib)
        {
            appGetLattice()->m_pUpdator->Update(1, TRUE);
        }
        appGetLattice()->m_pMeasurements->Report();
    }

    appGeneral(_T("\n========= all finished! ==========\n"));
    appQuitCLG();
    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
