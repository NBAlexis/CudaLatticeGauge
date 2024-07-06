//=============================================================================
// FILENAME : TestBoson.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [06/21/2024 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestUpdateSU4(CParameters& sParam)
{
    appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
    appGetLattice()->m_pUpdator->Update(5, FALSE);
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->Update(5, TRUE);

    return 0;
}

__REGIST_TEST(TestUpdateSU4, SUN, TestSU4, SUN);
__REGIST_TEST(TestUpdateSU4, SUN, TestSU4D, SUND);


//=============================================================================
// END OF FILE
//=============================================================================
