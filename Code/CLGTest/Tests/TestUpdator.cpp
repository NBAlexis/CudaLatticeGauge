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

    pAction->SetBeta(10.0f);
    appGetLattice()->m_pUpdator->Update(10, TRUE);
    //TArray<Real> allRes;
    //for (INT i = 0; i < 14; ++i)
    //{
    //    //update 10 steps for equilibrate
    //    Real fBeta = F(1.0) + i * F(0.5);
    //    pAction->SetBeta(fBeta);

    //    appGetLattice()->m_pUpdator->Update(20, FALSE);
    //    pMeasure->Reset();
    //    appGetLattice()->m_pUpdator->Update(50, TRUE);
    //    pMeasure->Report();
    //    allRes.AddItem(fBeta);
    //    allRes.AddItem(pMeasure->m_fLastRealResult);
    //}
    //
    //for (UINT i = 0; i < allRes.Num(); i+=2)
    //{
    //    appGeneral(_T("Final Res: (beta, <S>) = (%f, %f)\n"), allRes[i], allRes[i + 1]);
    //}

    return 0;
}

__REGIST_TEST(TestUpdator, Updator, TestUpdatorLeapFrog);

//__REGIST_TEST(TestUpdator, Updator, TestRandomXORWOW);

//=============================================================================
// END OF FILE
//=============================================================================
