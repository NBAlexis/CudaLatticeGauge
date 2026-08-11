//=============================================================================
// FILENAME : TaskRotFinalForce.cpp
//
// DESCRIPTION:
//   Offline golden-data task for the complete rotating HISQ fermion force.
//   The loaded field is an even pseudofermion.  Force construction is
//   delegated to CActionFermionHISQCombined so the ordinary, Naik, epsilon,
//   rotating, and HISQ-smearing-backward terms all follow production code.
//=============================================================================

#include "HISQ.h"

INT TaskRotFinalForceJob(CParameters& params)
{
    appSetupLog(params);

    CCString sCfgFile, sVFile, sOutFile;
    params.FetchStringValue(_T("CfgFile"), sCfgFile);
    params.FetchStringValue(_T("VFile"), sVFile);
    params.FetchStringValue(_T("OutputFile"), sOutFile);

    appGeneral(_T("[TaskRotFinalForce] CfgFile=%s\n"), sCfgFile.c_str());
    appGeneral(_T("[TaskRotFinalForce] VFile=%s\n"), sVFile.c_str());
    appGeneral(_T("[TaskRotFinalForce] OutputFile=%s\n"), sOutFile.c_str());

    if (!appInitialCLG(params))
    {
        appCrucial(_T("[TaskRotFinalForce] Initial Failed!\n"));
        return 1;
    }

    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1));
    CFieldFermionHISQSU3R* pFerm = dynamic_cast<CFieldFermionHISQSU3R*>(appGetLattice()->GetFieldById(2));
    CGaugeSmearingHISQSU3* pSmear = dynamic_cast<CGaugeSmearingHISQSU3*>(appGetLattice()->m_pGaugeSmearing[1]);
    CActionFermionHISQCombined* pAction =
        dynamic_cast<CActionFermionHISQCombined*>(appGetLattice()->GetActionById(1));
    if (NULL == pGauge || NULL == pFerm || NULL == pSmear || NULL == pAction)
    {
        appCrucial(_T("[TaskRotFinalForce] required gauge, rotating HISQ field, smearing, or action is missing\n"));
        appQuitCLG();
        return 2;
    }

    pGauge->InitialFieldWithFile(sCfgFile, EFFT_CLGBin);
    pFerm->InitialFieldWithFile(sVFile, EFFT_CLGBin);
    pFerm->m_bEvenPseudofermion = TRUE;
    pFerm->SetMass(F(0.0));
    pFerm->ZeroOnEvenOdd(FALSE);

    appGeneral(
        _T("[TaskRotFinalForce] Omega=%.17g Epsilon=%.17g ShiftHalfCoord=%d Even=%d\n"),
        static_cast<DOUBLE>(pFerm->GetOmega()),
        static_cast<DOUBLE>(pFerm->m_fEpsilon),
        static_cast<INT>(pFerm->m_bShiftHalfCoord),
        static_cast<INT>(pFerm->m_bEvenPseudofermion));

    // CalculateForce consumes the effective gauge cached by the HISQ
    // smearing object but deposits the derivative in thin-link coordinates.
    pSmear->GaugeSmearingC(pGauge);

    CFieldGauge* pForce = dynamic_cast<CFieldGauge*>(
        appGetLattice()->GetPooledFieldById(1, _T(__FILE__), __LINE__));
    if (NULL == pForce)
    {
        appCrucial(_T("[TaskRotFinalForce] cannot allocate force field\n"));
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
        appCrucial(_T("[TaskRotFinalForce] production CalculateForce failed\n"));
        pForce->Return();
        appQuitCLG();
        return 4;
    }

    // Match CIntegrator::CalcForceOfActions and QUDA LatticeMom:
    // closed force = TA(U_thin * F_thin^dagger).
    pForce->LeftMul(pGauge, FALSE, TRUE);
    pForce->TA();

    const CCString md5 = pForce->SaveToFile(sOutFile, EFFT_CLGBinDouble);
    appGeneral(_T("[TaskRotFinalForce] saved closed force md5=%s\n"), md5.c_str());
    pForce->Return();

    appGeneral(_T("[TaskRotFinalForce] done.\n"));
    appQuitCLG();
    return 0;
}
