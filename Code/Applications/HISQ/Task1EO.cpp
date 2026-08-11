//=============================================================================
// FILENAME : Task1EO.cpp
//
// DESCRIPTION:
//   Task 1 EO: pFerm <- (M^dagger M)_ee * v_e (even-even preconditioned).
//
// REVISION:
//=============================================================================

#include "HISQ.h"

INT Task1EOJob(CParameters& params)
{
    SCompareCtx ctx;
    const CCString sTag = _T("Task1EO");
    INT err = SetupCompare(params, ctx, sTag);
    if (0 != err)
    {
        appQuitCLG();
        return err;
    }

    ctx.pFerm->m_bEvenPseudofermion = TRUE;
    ctx.pFerm->SetMass(F(0.22));
    appGeneral(_T("[Task1EO] applying (M^dM)_ee*v_e with 2am=%f\n"),
        static_cast<DOUBLE>(ctx.pFerm->GetMass()));
    ctx.pFerm->DDdagger(1, 0, 0, ctx.effGauges.GetData(), NULL, NULL);
    const CCString md5 = ctx.pFerm->SaveToFile(ctx.sOutFile, EFFT_CLGBinDouble);
    appGeneral(_T("[Task1EO] saved %s md5=%s\n"), ctx.sOutFile.c_str(), md5.c_str());

    appGeneral(_T("[Task1EO] done.\n"));
    appQuitCLG();
    return 0;
}
