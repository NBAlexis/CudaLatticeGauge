//=============================================================================
// FILENAME : Task4Action.cpp
//
// DESCRIPTION:
//   Task 4: Tree-improved Symanzik gauge action energy.
//
// REVISION:
//=============================================================================

#include "HISQ.h"

INT Task4ActionJob(CParameters& params)
{
    SCompareCtx ctx;
    const CCString sTag = _T("Task4Action");
    INT err = SetupGaugeOnly(params, ctx, sTag);
    if (0 != err)
    {
        appQuitCLG();
        return err;
    }

    CFieldGaugeSU3TreeImproved* pTree = dynamic_cast<CFieldGaugeSU3TreeImproved*>(ctx.pGauge);
    if (NULL == pTree)
    {
        appCrucial(_T("[Task4Action] gauge field is not CFieldGaugeSU3TreeImproved\n"));
        appQuitCLG();
        return 5;
    }

    CAction* pAction = appGetLattice()->GetActionById(1);
    if (NULL == pAction)
    {
        appCrucial(_T("[Task4Action] gauge action id=1 not found\n"));
        appQuitCLG();
        return 6;
    }

    appGeneral(_T("[Task4Action] computing tree-improved gauge action energy\n"));
    const DOUBLE S = pAction->Energy(FALSE, 1, 0, 0, ctx.rawGauges.GetData(), NULL, NULL, NULL);
    appGeneral(_T("[Task4Action] S = %.17g\n"), S);

    CCString sData4;
    sData4.Format(_T("%.17g\n"), S);
    appGetFileSystem()->WriteAllText(ctx.sOutFile.c_str(), sData4);
    appGeneral(_T("[Task4Action] saved %s\n"), ctx.sOutFile.c_str());

    appGeneral(_T("[Task4Action] done.\n"));
    appQuitCLG();
    return 0;
}
