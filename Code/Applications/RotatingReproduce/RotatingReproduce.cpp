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

    INT iBeforeEquib = 5;
    params.FetchValueINT(_T("BeforeEquvibStep"), iBeforeEquib);
    INT iAfterEquib = 200;
    params.FetchValueINT(_T("EquvibStep"), iAfterEquib);

    appSetupLog(params);
    if (!appInitialCLG(params))
    {
        return 1;
    }

    CActionGaugePlaquetteRotating * pGauageAction = dynamic_cast<CActionGaugePlaquetteRotating *>(appGetLattice()->GetActionById(1));
    CMeasureAMomentumJG * pJG = dynamic_cast<CMeasureAMomentumJG *>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    CMeasureAMomentumJF * pJF = dynamic_cast<CMeasureAMomentumJF *>(appGetLattice()->m_pMeasurements->GetMeasureById(2));
    if (NULL == pGauageAction || NULL == pJG || NULL == pJF)
    {
        appCrucial(_T("Rotating gauge action or measurement not found!\n"));
        return 1;
    }
    
    appGeneral(_T("Start up.\n"));
    appGeneral(_T("1 run: Omega = 0.00.\n"));

    pGauageAction->SetOmega(F(0.0));
    appGetLattice()->m_pGaugeField->InitialField(EFIT_Random);

    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
    {
        appGetLattice()->m_pUpdator->Update(1, FALSE);
    }
    
    pJG->Reset();
    pJF->Reset();
    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iAfterEquib)
    {
        appGetLattice()->m_pUpdator->Update(1, TRUE);
    }
    pJG->Report();
    pJF->Report();


    appGeneral(_T("2 run: Omega = 0.02.\n"));

    pGauageAction->SetOmega(F(0.02));
    appGetLattice()->m_pGaugeField->InitialField(EFIT_Random);

    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
    {
        appGetLattice()->m_pUpdator->Update(1, FALSE);
    }

    pJG->Reset();
    pJF->Reset();
    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iAfterEquib)
    {
        appGetLattice()->m_pUpdator->Update(1, TRUE);
    }
    pJG->Report();
    pJF->Report();

    appGeneral(_T("3 run: Omega = 0.04.\n"));

    pGauageAction->SetOmega(F(0.04));
    appGetLattice()->m_pGaugeField->InitialField(EFIT_Random);

    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
    {
        appGetLattice()->m_pUpdator->Update(1, FALSE);
    }

    pJG->Reset();
    pJF->Reset();
    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iAfterEquib)
    {
        appGetLattice()->m_pUpdator->Update(1, TRUE);
    }
    pJG->Report();
    pJF->Report();

    appGeneral(_T("4 run: Omega = 0.06.\n"));

    pGauageAction->SetOmega(F(0.06));
    appGetLattice()->m_pGaugeField->InitialField(EFIT_Random);

    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
    {
        appGetLattice()->m_pUpdator->Update(1, FALSE);
    }

    pJG->Reset();
    pJF->Reset();
    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iAfterEquib)
    {
        appGetLattice()->m_pUpdator->Update(1, TRUE);
    }
    pJG->Report();
    pJF->Report();

    appGeneral(_T("5 run: Omega = 0.08.\n"));

    pGauageAction->SetOmega(F(0.08));
    appGetLattice()->m_pGaugeField->InitialField(EFIT_Random);

    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
    {
        appGetLattice()->m_pUpdator->Update(1, FALSE);
    }

    pJG->Reset();
    pJF->Reset();
    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iAfterEquib)
    {
        appGetLattice()->m_pUpdator->Update(1, TRUE);
    }
    pJG->Report();
    pJF->Report();

    appGeneral(_T("6 run: Omega = 0.10.\n"));

    pGauageAction->SetOmega(F(0.1));
    appGetLattice()->m_pGaugeField->InitialField(EFIT_Random);

    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iBeforeEquib)
    {
        appGetLattice()->m_pUpdator->Update(1, FALSE);
    }

    pJG->Reset();
    pJF->Reset();
    appGetLattice()->m_pUpdator->SetConfigurationCount(0);
    while (appGetLattice()->m_pUpdator->GetConfigurationCount() < iAfterEquib)
    {
        appGetLattice()->m_pUpdator->Update(1, TRUE);
    }
    pJG->Report();
    pJF->Report();

    appQuitCLG();
    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
