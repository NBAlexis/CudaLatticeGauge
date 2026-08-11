#include "HISQ.h"


INT TaskRotDJob(CParameters& params)
{
    appSetupLog(params);

    CCString cfgFile, sourceFile, outputFile;
    params.FetchStringValue(_T("CfgFile"), cfgFile);
    params.FetchStringValue(_T("VFile"), sourceFile);
    params.FetchStringValue(_T("OutputFile"), outputFile);

    appGeneral(_T("[TaskRotInverse] CfgFile   = %s\n"), cfgFile.c_str());
    appGeneral(_T("[TaskRotInverse] VFile     = %s\n"), sourceFile.c_str());
    appGeneral(_T("[TaskRotInverse] OutputFile= %s\n"), outputFile.c_str());

    if (!appInitialCLG(params))
    {
        appCrucial(_T("[TaskRotInverse] Initial Failed!\n"));
        return 1;
    }

    CFieldGaugeSU3* gauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1));
    CFieldFermionHISQSU3R* fermion =
        dynamic_cast<CFieldFermionHISQSU3R*>(appGetLattice()->GetFieldById(2));
    CGaugeSmearingHISQSU3* smearing =
        dynamic_cast<CGaugeSmearingHISQSU3*>(appGetLattice()->m_pGaugeSmearing[1]);
    if (NULL == gauge || NULL == fermion || NULL == smearing)
    {
        appCrucial(_T("[TaskRotInverse] required rotating HISQ fields not found\n"));
        appQuitCLG();
        return 2;
    }

    gauge->InitialFieldWithFile(cfgFile, EFFT_CLGBin);
    fermion->InitialFieldWithFile(sourceFile, EFFT_CLGBin);
    smearing->GaugeSmearingC(gauge);
    TArray<const CFieldGauge*> effectiveGauges;
    effectiveGauges.AddItem(smearing->GetEffectiveGauge());

    fermion->m_bEvenPseudofermion = FALSE;
    // PyQUDA's rotating matrix convention is M = 2m - D_CLG, while CLGLib's
    // InverseD solves (2am + D_CLG)x=b.  Negating both CLGLib's mass and the
    // loaded right-hand side gives the identical linear system.
    fermion->SetMass(-fermion->GetMass());
    fermion->ScalarMultply(F(-1.0));
    fermion->InverseD(1, 0, 0, effectiveGauges.GetData(), NULL, NULL);
    const CCString md5 = fermion->SaveToFile(outputFile, EFFT_CLGBinDouble);
    appGeneral(_T("[TaskRotInverse] saved %s md5=%s\n"), outputFile.c_str(), md5.c_str());

    appQuitCLG();
    return 0;
}
