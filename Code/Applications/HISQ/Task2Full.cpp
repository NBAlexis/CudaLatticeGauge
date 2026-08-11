//=============================================================================
// FILENAME : Task2Full.cpp
//
// DESCRIPTION:
//   Task 2 Full: pFerm <- (M^dagger M)^{1/4} * v on the full lattice via D_MC.
//
// REVISION:
//=============================================================================

#include "HISQ.h"

INT Task2FullJob(CParameters& params)
{
    SCompareCtx ctx;
    const CCString sTag = _T("Task2Full");
    INT err = SetupCompare(params, ctx, sTag);
    if (0 != err)
    {
        appQuitCLG();
        return err;
    }

    ctx.pFerm->m_bEvenPseudofermion = FALSE;
    ctx.pFerm->SetMass(F(0.0));
    appGeneral(_T("[Task2Full] applying (M^dM)^{1/4}*v via D_MC with 2am=%f\n"),
        static_cast<DOUBLE>(ctx.pFerm->GetMass()));
    ctx.pFerm->D_MC(1, 0, 0, ctx.effGauges.GetData(), NULL, NULL);
    const CCString md5 = ctx.pFerm->SaveToFile(ctx.sOutFile, EFFT_CLGBinDouble);
    appGeneral(_T("[Task2Full] saved %s md5=%s\n"), ctx.sOutFile.c_str(), md5.c_str());

    appGeneral(_T("[Task2Full] done.\n"));
    appQuitCLG();
    return 0;
}
