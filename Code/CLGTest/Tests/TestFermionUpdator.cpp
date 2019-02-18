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
    appGetLattice()->m_pUpdator->Update(20, FALSE);
    appGetLattice()->m_pUpdator->Update(20, TRUE);
    return 0;
}

__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdator);

UINT TestFermionUpdatorL(CParameters& sParam)
{
    appGetLattice()->m_pUpdator->Update(1, TRUE);
    return 0;
}


__REGIST_TEST(TestFermionUpdatorL, Updator, TestFermionUpdatorLargeScale);


//=============================================================================
// END OF FILE
//=============================================================================
