//=============================================================================
// FILENAME : Task2EO.cpp
//
// DESCRIPTION:
//   Task 2 EO: pFerm <- (M^dagger M)_ee^{1/4} * v_e via D_MC with Even=1.
//
// REVISION:
//=============================================================================

#include "HISQ.h"

INT Task2EOJob(CParameters& params)
{
    SCompareCtx ctx;
    const CCString sTag = _T("Task2EO");
    INT err = SetupCompare(params, ctx, sTag);
    if (0 != err)
    {
        appQuitCLG();
        return err;
    }

    ctx.pFerm->m_bEvenPseudofermion = TRUE;
    ctx.pFerm->SetMass(F(0.0));
    // Zero odd part of v so the EE-mode rational does not pass odd through "norm * v_odd".
    ctx.pFerm->ZeroOnEvenOdd(FALSE);
    appGeneral(_T("[Task2EO] applying (M^dM)_ee^{1/4}*v_e via D_MC with 2am=%f\n"),
        static_cast<DOUBLE>(ctx.pFerm->GetMass()));
    ctx.pFerm->D_MC(1, 0, 0, ctx.effGauges.GetData(), NULL, NULL);
    // Defensive: zero odd again on the result so file matches PyQUDA's odd-zeroed convention.
    // ctx.pFerm->ZeroOnEvenOdd(FALSE);
    const CCString md5 = ctx.pFerm->SaveToFile(ctx.sOutFile, EFFT_CLGBinDouble);
    appGeneral(_T("[Task2EO] saved %s md5=%s\n"), ctx.sOutFile.c_str(), md5.c_str());

    appGeneral(_T("[Task2EO] done.\n"));
    appQuitCLG();
    return 0;
}
