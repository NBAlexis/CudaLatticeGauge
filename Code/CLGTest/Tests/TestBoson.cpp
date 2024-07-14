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
    appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
    appGetLattice()->m_pUpdator->Update(3, FALSE);
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->Update(10, TRUE);

    return 0;
}

__REGIST_TEST(TestUpdateBoson, Boson, TestBosonU1NoGauge, U1NoGauge);
__REGIST_TEST(TestUpdateBoson, Boson, TestBosonU1, U1);
__REGIST_TEST(TestUpdateBoson, Boson, TestBosonSU4, SU4);
__REGIST_TEST(TestUpdateBoson, Boson, TestBosonU1ExternalGaugeField, ExternalU1);
__REGIST_TEST(TestUpdateBoson, Boson, TestBosonU1D, U1D);
__REGIST_TEST(TestUpdateBoson, Boson, TestBosonU1P, U1P);
__REGIST_TEST(TestUpdateBoson, Boson, TestBosonRotation, Rotation);
__REGIST_TEST(TestUpdateBoson, Boson, TestBosonRotationD, RotationD);
__REGIST_TEST(TestUpdateBoson, Boson, TestBosonRotationP, RotationP);

//=============================================================================
// END OF FILE
//=============================================================================
