//=============================================================================
// FILENAME : Task6FermionAction.cpp
//
// DESCRIPTION:
//   Task 6: Combined 3-field fermion action energy.
//
// REVISION:
//=============================================================================

#include "HISQ.h"

INT Task6FermionActionJob(CParameters& params)
{
    appSetupLog(params);

    CCString sCfgFile, sV1File, sV2File, sV3File, sOutFile;
    params.FetchStringValue(_T("CfgFile"), sCfgFile);
    params.FetchStringValue(_T("V1File"), sV1File);
    params.FetchStringValue(_T("V2File"), sV2File);
    params.FetchStringValue(_T("V3File"), sV3File);
    params.FetchStringValue(_T("OutputFile"), sOutFile);

    appGeneral(_T("[Task6FermionAction] CfgFile=%s\n"), sCfgFile.c_str());
    appGeneral(_T("[Task6FermionAction] V1File=%s V2File=%s V3File=%s\n"),
        sV1File.c_str(), sV2File.c_str(), sV3File.c_str());
    appGeneral(_T("[Task6FermionAction] OutputFile=%s\n"), sOutFile.c_str());

    if (!appInitialCLG(params))
    {
        appCrucial(_T("[Task6FermionAction] Initial Failed!\n"));
        return 1;
    }

    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1));
    if (NULL == pGauge)
    {
        appCrucial(_T("[Task6FermionAction] gauge field id=1 not found\n"));
        appQuitCLG();
        return 2;
    }

    CFieldFermionHISQSU3* pFerm1 = dynamic_cast<CFieldFermionHISQSU3*>(appGetLattice()->GetFieldById(2));
    CFieldFermionHISQSU3* pFerm2 = dynamic_cast<CFieldFermionHISQSU3*>(appGetLattice()->GetFieldById(3));
    CFieldFermionHISQSU3* pFerm3 = dynamic_cast<CFieldFermionHISQSU3*>(appGetLattice()->GetFieldById(4));
    if (NULL == pFerm1 || NULL == pFerm2 || NULL == pFerm3)
    {
        appCrucial(_T("[Task6FermionAction] one or more fermion fields not found\n"));
        appQuitCLG();
        return 3;
    }

    CGaugeSmearingHISQSU3* pSmear = dynamic_cast<CGaugeSmearingHISQSU3*>(appGetLattice()->m_pGaugeSmearing[1]);
    if (NULL == pSmear)
    {
        appCrucial(_T("[Task6FermionAction] HISQ gauge smearing not found\n"));
        appQuitCLG();
        return 4;
    }

    pGauge->InitialFieldWithFile(sCfgFile, EFFT_CLGBin);
    pFerm1->InitialFieldWithFile(sV1File, EFFT_CLGBin);
    pFerm2->InitialFieldWithFile(sV2File, EFFT_CLGBin);
    pFerm3->InitialFieldWithFile(sV3File, EFFT_CLGBin);

    pFerm1->m_bEvenPseudofermion = TRUE;
    pFerm1->SetMass(F(0.0));
    pFerm1->ZeroOnEvenOdd(FALSE);
    pFerm2->m_bEvenPseudofermion = TRUE;
    pFerm2->SetMass(F(0.0));
    pFerm2->ZeroOnEvenOdd(FALSE);
    pFerm3->m_bEvenPseudofermion = TRUE;
    pFerm3->SetMass(F(0.0));
    pFerm3->ZeroOnEvenOdd(FALSE);

    pSmear->GaugeSmearingC(pGauge);

    TArray<const CFieldGauge*> rawGauges;
    rawGauges.AddItem(pGauge);

    CAction* pAction = appGetLattice()->GetActionById(2);
    if (NULL == pAction)
    {
        appCrucial(_T("[Task6FermionAction] fermion action id=2 not found\n"));
        appQuitCLG();
        return 5;
    }

    appGeneral(_T("[Task6FermionAction] computing combined fermion action energy (3 fields)\n"));
    const DOUBLE S = pAction->Energy(FALSE, 1, 0, 0, rawGauges.GetData(), NULL, NULL, NULL);
    appGeneral(_T("[Task6FermionAction] S = %.17g\n"), S);

    CCString sData;
    sData.Format(_T("%.17g\n"), S);
    appGetFileSystem()->WriteAllText(sOutFile.c_str(), sData);
    appGeneral(_T("[Task6FermionAction] saved %s\n"), sOutFile.c_str());

    appGeneral(_T("[Task6FermionAction] done.\n"));
    appQuitCLG();
    return 0;
}
