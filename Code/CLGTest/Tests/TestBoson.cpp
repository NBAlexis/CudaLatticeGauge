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
    appGetLattice()->m_pUpdator->Update(20, TRUE);

    return 0;
}

__REGIST_TEST(TestUpdateBoson, UpdatorBoson, TestBosonU1NoGauge, U1NoGauge);
__REGIST_TEST(TestUpdateBoson, UpdatorBoson, TestBosonU1, U1);
__REGIST_TEST(TestUpdateBoson, UpdatorBoson, TestBosonU1ExternalGaugeField, ExternalU1);

//=============================================================================
// END OF FILE
//=============================================================================
