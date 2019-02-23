//=============================================================================
// FILENAME : TestFermionUpdator.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [02/07/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestFermionUpdator(CParameters& sParam)
{
    Real fExpected = F(0.17);
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
    CMeasurePlaqutteEnergy* pMeasure = dynamic_cast<CMeasurePlaqutteEnergy*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }

    appGetLattice()->m_pUpdator->Update(20, FALSE);
    appGetLattice()->m_pUpdator->Update(20, TRUE);
    Real fRes = pMeasure->m_fLastRealResult;
    appGeneral(_T("res : expected=%f res=%f"), fExpected, fRes);

    if (appAbs(fRes - fExpected) > F(0.01))
    {
        return 1;
    }
    return 0;
}

__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdator);

UINT TestFermionUpdatorWithMesonCorrelator(CParameters& sParam)
{
    CMeasureMesonCorrelator* pMeasure = dynamic_cast<CMeasureMesonCorrelator*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    if (NULL == pMeasure)
    {
        return 1;
    }
    Real fExpected = F(0.625);
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);
#if _CLG_DEBUG
    appGetLattice()->m_pUpdator->Update(10, FALSE);
    appGetLattice()->m_pUpdator->Update(20, TRUE);

    Real fRes = pMeasure->m_lstResults[0][0];
    appGeneral(_T("res : expected=%f res=%f"), fExpected, fRes);
    if (appAbs(fRes - fExpected) > F(0.01))
    {
        return 1;
    }
#else
    appGetLattice()->m_pUpdator->Update(20, FALSE);
    appGetLattice()->m_pUpdator->Update(50, TRUE);

    Real fRes = pMeasure->m_lstResults[0][0];
    appGeneral(_T("res : expected=%f res=%f"), fExpected, fRes);
    if (appAbs(fRes - fExpected) > F(0.005))
    {
        return 1;
    }
#endif
    return 0;
}

#if _CLG_DEBUG
__REGIST_TEST(TestFermionUpdatorWithMesonCorrelator, Updator, TestFermionUpdatorWithMesonCorrelator);
#else
__REGIST_TEST(TestFermionUpdatorWithMesonCorrelator, Updator, TestFermionUpdatorWithMesonCorrelatorRelease);
#endif

UINT TestFermionUpdatorL(CParameters& sParam)
{
    appGetLattice()->m_pUpdator->Update(1, TRUE);
    return 0;
}


__REGIST_TEST(TestFermionUpdatorL, Updator, TestFermionUpdatorLargeScale);


//=============================================================================
// END OF FILE
//=============================================================================
