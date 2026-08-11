//=============================================================================
// FILENAME : Task1Full.cpp
//
// DESCRIPTION:
//   Task 1 Full: pFerm <- (2am + D)*v on the full lattice.
//
// REVISION:
//=============================================================================

#include "HISQ.h"

INT Task1FullJob(CParameters& params)
{
    SCompareCtx ctx;
    const CCString sTag = _T("Task1Full");
    INT err = SetupCompare(params, ctx, sTag);
    if (0 != err)
    {
        appQuitCLG();
        return err;
    }

    ctx.pFerm->m_bEvenPseudofermion = FALSE;
    ctx.pFerm->SetMass(F(0.22));
    appGeneral(_T("[Task1Full] applying (2am + D)*v with 2am=%f\n"),
        static_cast<DOUBLE>(ctx.pFerm->GetMass()));
    ctx.pFerm->DWithMass(1, 0, 0, ctx.effGauges.GetData(), NULL, NULL, ctx.pFerm->GetMass());
    const CCString md5 = ctx.pFerm->SaveToFile(ctx.sOutFile, EFFT_CLGBinDouble);
    appGeneral(_T("[Task1Full] saved %s md5=%s\n"), ctx.sOutFile.c_str(), md5.c_str());

    appGeneral(_T("[Task1Full] done.\n"));
    appQuitCLG();
    return 0;
}
