//=============================================================================
// FILENAME : TestBoson.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [06/21/2024 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestUpdateBoson(CParameters& sParam)
{
    //appGetLattice()->m_pUpdator->Update(1, FALSE);
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->Update(100, TRUE);

    return 0;
}

__REGIST_TEST(TestUpdateBoson, Updator, TestBosonU1NoGauge, BosonU1NoGauge);

//=============================================================================
// END OF FILE
//=============================================================================
