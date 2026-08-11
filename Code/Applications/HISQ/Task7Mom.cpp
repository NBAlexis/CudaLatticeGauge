//=============================================================================
// FILENAME : Task7Mom.cpp
//
// DESCRIPTION:
//   Task 7: Verify gauge momentum conventions with GaugeMomentumFactor.
//   Two modes:
//     1. Load momentum from binary and compute kinematic energy.
//     2. Generate random Gaussian momentum and compute kinematic energy.
//
//   Used to verify CLGLib vs QUDA momentum conventions.
//=============================================================================

#include "HISQ.h"

INT Task7MomJob(CParameters& params)
{
    SCompareCtx ctx;
    const CCString sTag = _T("Task7Mom");
    INT err = SetupGaugeOnly(params, ctx, sTag);
    if (0 != err)
    {
        appQuitCLG();
        return err;
    }

    // Check if we should generate random momentum
    INT generateRandom = 0;
    params.FetchValueINT(_T("GenerateRandom"), generateRandom);

    if (generateRandom > 0)
    {
        // Generate random Gaussian momentum
        appGeneral(_T("[%s] Generating random Gaussian momentum (GaugeMomentumFactor applied)\n"), sTag.c_str());
        ctx.pGauge->InitialField(EFIT_RandomGenerator);
        appGeneral(_T("[%s] random momentum generated\n"), sTag.c_str());

        // Save generated momentum for external verification
        CCString sOutMomFile;
        if (params.FetchStringValue(_T("OutMomFile"), sOutMomFile))
        {
            ctx.pGauge->SaveToFile(sOutMomFile);
            appGeneral(_T("[%s] saved generated momentum to %s\n"), sTag.c_str(), sOutMomFile.c_str());
        }
    }
    else
    {
        // Load momentum from binary file into the gauge field
        CCString sMomFile;
        params.FetchStringValue(_T("MomFile"), sMomFile);
        appGeneral(_T("[%s] MomFile = %s\n"), sTag.c_str(), sMomFile.c_str());

        ctx.pGauge->InitialFieldWithFile(sMomFile, EFFT_CLGBinDouble);
        appGeneral(_T("[%s] momentum loaded\n"), sTag.c_str());
    }

    // Compute kinematic energy
    const DOUBLE kin = ctx.pGauge->CalculateKinematicEnergy();
    appGeneral(_T("[%s] CalculateKinematicEnergy = %.17g\n"), sTag.c_str(), kin);

    // Write result to output file
    CCString sData;
    sData.Format(_T("%.17g\n"), kin);
    appGetFileSystem()->WriteAllText(ctx.sOutFile.c_str(), sData);
    appGeneral(_T("[%s] saved %s\n"), sTag.c_str(), ctx.sOutFile.c_str());

    appGeneral(_T("[%s] done.\n"), sTag.c_str());
    appQuitCLG();
    return 0;
}
