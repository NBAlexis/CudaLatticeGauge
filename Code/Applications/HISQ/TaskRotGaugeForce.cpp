// Offline golden-data task for the complete rotating tree-improved gauge force.
#include "HISQ.h"

INT TaskRotGaugeForceJob(CParameters& params)
{
    appSetupLog(params);

    CCString sCfgFile, sOutFile;
    params.FetchStringValue(_T("CfgFile"), sCfgFile);
    params.FetchStringValue(_T("OutputFile"), sOutFile);

    if (!appInitialCLG(params))
    {
        appCrucial(_T("[TaskRotGaugeForce] Initial Failed!\n"));
        return 1;
    }

    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1));
    CActionGaugePlaquetteRotating* pAction =
        dynamic_cast<CActionGaugePlaquetteRotating*>(appGetLattice()->GetActionById(1));
    if (NULL == pGauge || NULL == pAction)
    {
        appCrucial(_T("[TaskRotGaugeForce] required gauge field or rotating gauge action is missing\n"));
        appQuitCLG();
        return 2;
    }

    pGauge->InitialFieldWithFile(sCfgFile, EFFT_CLGBin);
    CFieldGauge* pForce = dynamic_cast<CFieldGauge*>(
        appGetLattice()->GetPooledFieldById(1, _T(__FILE__), __LINE__));
    if (NULL == pForce)
    {
        appCrucial(_T("[TaskRotGaugeForce] cannot allocate force field\n"));
        appQuitCLG();
        return 3;
    }
    pForce->Zero();

    TArray<const CFieldGauge*> gauges;
    TArray<CFieldGauge*> forces;
    gauges.AddItem(pGauge);
    forces.AddItem(pForce);
    if (!pAction->CalculateForce(
        1, 0, gauges.GetData(), NULL, forces.GetData(), NULL, NULL, ESP_Once))
    {
        appCrucial(_T("[TaskRotGaugeForce] production CalculateForce failed\n"));
        pForce->Return();
        appQuitCLG();
        return 4;
    }

    // Match CIntegrator::CalcForceOfActions and QUDA LatticeMom.
    pForce->LeftMul(pGauge, FALSE, TRUE);
    pForce->TA();

    const CCString md5 = pForce->SaveToFile(sOutFile, EFFT_CLGBinDouble);
    appGeneral(_T("[TaskRotGaugeForce] Omega=%.17g saved closed force md5=%s\n"),
        static_cast<DOUBLE>(pAction->GetOmega()), md5.c_str());
    pForce->Return();
    appQuitCLG();
    return 0;
}
