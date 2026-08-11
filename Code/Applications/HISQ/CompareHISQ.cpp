//=============================================================================
// FILENAME : CompareHISQ.cpp
//
// DESCRIPTION:
//   Shared helpers for the CLGLib vs PyQUDA HISQ operator comparison tasks.
//   Each task has its own implementation file (Task1Full.cpp, Task1EO.cpp, etc.).
//
// REVISION:
//=============================================================================

#include "HISQ.h"

INT SetupCompare(CParameters& params, SCompareCtx& ctx, const CCString& sTag)
{
    appSetupLog(params);

    params.FetchStringValue(_T("CfgFile"), ctx.sCfgFile);
    params.FetchStringValue(_T("VFile"), ctx.sVFile);
    params.FetchStringValue(_T("OutputFile"), ctx.sOutFile);

    appGeneral(_T("[%s] CfgFile   = %s\n"), sTag.c_str(), ctx.sCfgFile.c_str());
    appGeneral(_T("[%s] VFile     = %s\n"), sTag.c_str(), ctx.sVFile.c_str());
    appGeneral(_T("[%s] OutputFile= %s\n"), sTag.c_str(), ctx.sOutFile.c_str());

    if (!appInitialCLG(params))
    {
        appCrucial(_T("[%s] Initial Failed!\n"), sTag.c_str());
        return 1;
    }

    ctx.pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1));
    if (NULL == ctx.pGauge)
    {
        appCrucial(_T("[%s] gauge field id=1 not found or wrong type\n"), sTag.c_str());
        return 2;
    }

    ctx.pFerm = dynamic_cast<CFieldFermionHISQSU3*>(appGetLattice()->GetFieldById(2));
    if (NULL == ctx.pFerm)
    {
        appCrucial(_T("[%s] fermion field id=2 not found or wrong type\n"), sTag.c_str());
        return 3;
    }

    ctx.pSmear = dynamic_cast<CGaugeSmearingHISQSU3*>(appGetLattice()->m_pGaugeSmearing[1]);
    if (NULL == ctx.pSmear)
    {
        appCrucial(_T("[%s] HISQ gauge smearing for fieldId=1 not found\n"), sTag.c_str());
        return 4;
    }

    ctx.pGauge->InitialFieldWithFile(ctx.sCfgFile, EFFT_CLGBin);
    ctx.pFerm->InitialFieldWithFile(ctx.sVFile, EFFT_CLGBin);

    ctx.pSmear->GaugeSmearingC(ctx.pGauge);
    ctx.effGauges.AddItem(ctx.pSmear->GetEffectiveGauge());
    ctx.rawGauges.AddItem(ctx.pGauge);

    appGeneral(_T("[%s] before run: m_f2am=%f m_bEvenPseudofermion=%d\n"),
        sTag.c_str(),
        static_cast<DOUBLE>(ctx.pFerm->GetMass()),
        static_cast<INT>(ctx.pFerm->m_bEvenPseudofermion));

    return 0;
}

INT SetupGaugeOnly(CParameters& params, SCompareCtx& ctx, const CCString& sTag)
{
    appSetupLog(params);

    params.FetchStringValue(_T("CfgFile"), ctx.sCfgFile);
    params.FetchStringValue(_T("OutputFile"), ctx.sOutFile);

    appGeneral(_T("[%s] CfgFile   = %s\n"), sTag.c_str(), ctx.sCfgFile.c_str());
    appGeneral(_T("[%s] OutputFile= %s\n"), sTag.c_str(), ctx.sOutFile.c_str());

    if (!appInitialCLG(params))
    {
        appCrucial(_T("[%s] Initial Failed!\n"), sTag.c_str());
        return 1;
    }

    ctx.pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->GetFieldById(1));
    if (NULL == ctx.pGauge)
    {
        appCrucial(_T("[%s] gauge field id=1 not found or wrong type\n"), sTag.c_str());
        return 2;
    }

    ctx.pFerm = dynamic_cast<CFieldFermionHISQSU3*>(appGetLattice()->GetFieldById(2));
    if (NULL == ctx.pFerm)
    {
        appCrucial(_T("[%s] fermion field id=2 not found or wrong type\n"), sTag.c_str());
        return 3;
    }

    ctx.pSmear = dynamic_cast<CGaugeSmearingHISQSU3*>(appGetLattice()->m_pGaugeSmearing[1]);
    if (NULL == ctx.pSmear)
    {
        appCrucial(_T("[%s] HISQ gauge smearing for fieldId=1 not found\n"), sTag.c_str());
        return 4;
    }

    ctx.pGauge->InitialFieldWithFile(ctx.sCfgFile, EFFT_CLGBin);
    ctx.rawGauges.AddItem(ctx.pGauge);

    CFieldGaugeSU3TreeImproved* pTree = dynamic_cast<CFieldGaugeSU3TreeImproved*>(ctx.pGauge);
    if (NULL != pTree)
    {
        appGeneral(_T("[%s] gauge is tree-improved, RectOverPlaq=%f\n"),
            sTag.c_str(), static_cast<DOUBLE>(pTree->m_fRectOverPlaq));
    }

    return 0;
}
