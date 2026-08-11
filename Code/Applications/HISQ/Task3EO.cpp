//=============================================================================
// FILENAME : Task3EO.cpp
//
// DESCRIPTION:
//   Task 3 EO: S = v_e^dagger * (M^dagger M)_ee^{-1/4} * v_e via Action::Energy.
//
// REVISION:
//=============================================================================

#include "HISQ.h"

INT Task3EOJob(CParameters& params)
{
    SCompareCtx ctx;
    const CCString sTag = _T("Task3EO");
    INT err = SetupCompare(params, ctx, sTag);
    if (0 != err)
    {
        appQuitCLG();
        return err;
    }

    ctx.pFerm->m_bEvenPseudofermion = TRUE;
    ctx.pFerm->SetMass(F(0.0));
    // Zero odd part of v so Energy()'s Dot only counts even contribution.
    ctx.pFerm->ZeroOnEvenOdd(FALSE);

    // Action2 in YAML => byId=2 is the fermion action.  byId=1 is plaquette.
    CAction* pAction = appGetLattice()->GetActionById(2);
    if (NULL == pAction)
    {
        appCrucial(_T("[Task3EO] fermion action id=2 not found\n"));
        appQuitCLG();
        return 5;
    }
    appGeneral(_T("[Task3EO] computing v_e^d (M^dM)_ee^{-1/4} v_e via Action::Energy\n"));
    const DOUBLE S = pAction->Energy(FALSE, 1, 0, 0, ctx.rawGauges.GetData(), NULL, NULL, NULL);
    appGeneral(_T("[Task3EO] S = %.17g\n"), S);

    CCString sData3;
    sData3.Format(_T("%.17g\n"), S);
    appGetFileSystem()->WriteAllText(ctx.sOutFile.c_str(), sData3);
    appGeneral(_T("[Task3EO] saved %s\n"), ctx.sOutFile.c_str());

    appGeneral(_T("[Task3EO] done.\n"));
    appQuitCLG();
    return 0;
}
