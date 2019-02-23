//=============================================================================
// FILENAME : TestUpdator.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [01/28/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestUpdator(CParameters& sParam)
{
    Real fExpected = F(0.2064);
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);

    //we calculate staple energy from beta = 1 - 6
    CActionGaugePlaquette * pAction = dynamic_cast<CActionGaugePlaquette*>(appGetLattice()->GetActionById(1));
    if (NULL == pAction)
    {
        return 1;
    }
    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    //pAction->SetBeta(F(3.0));

    //Equilibration
    appGetLattice()->m_pUpdator->Update(5, FALSE);

    //Measure
    pMeasure->Reset();
    appGetLattice()->m_pUpdator->Update(20, TRUE);

    Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f"), fExpected, fRes);

    if (appAbs(fRes - fExpected) > F(0.005))
    {
        return 1;
    }

    return 0;
}

__REGIST_TEST(TestUpdator, Updator, TestUpdatorLeapFrog);

__REGIST_TEST(TestUpdator, Updator, TestUpdatorOmelyan);


//=============================================================================
// END OF FILE
//=============================================================================
