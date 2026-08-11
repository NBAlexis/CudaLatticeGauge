//=============================================================================
// FILENAME : Task5Plaq.cpp
//
// DESCRIPTION:
//   Task 5: Average plaquette.
//
// REVISION:
//=============================================================================

#include "HISQ.h"

INT Task5PlaqJob(CParameters& params)
{
    SCompareCtx ctx;
    const CCString sTag = _T("Task5Plaq");
    INT err = SetupGaugeOnly(params, ctx, sTag);
    if (0 != err)
    {
        appQuitCLG();
        return err;
    }

    appGeneral(_T("[Task5Plaq] computing average plaquette\n"));
    // Use CalculatePlaqutteEnergyOriginal to get the pure plaquette term
    // (unaffected by any tree-improved / rectangular overrides).
    const DOUBLE energy = ctx.pGauge->CalculatePlaqutteEnergyOriginal(1.0 / 3.0);
    const DOUBLE avgPlaq = 1.0 - energy / (6.0 * _HC_Volume);
    appGeneral(_T("[Task5Plaq] energy=%.17g avgPlaq=%.17g\n"), energy, avgPlaq);

    CCString sData5;
    sData5.Format(_T("%.17g\n"), avgPlaq);
    appGetFileSystem()->WriteAllText(ctx.sOutFile.c_str(), sData5);
    appGeneral(_T("[Task5Plaq] saved %s\n"), ctx.sOutFile.c_str());

    appGeneral(_T("[Task5Plaq] done.\n"));
    appQuitCLG();
    return 0;
}
