//=============================================================================
// FILENAME : TestZ3Symmetry.cpp
// 
// DESCRIPTION:
//
//     Test the Z3 Symmetry
//
// REVISION:
//  [06/23/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestZ3Symmetry(CParameters&)
{
    appGetLattice()->m_pUpdator->Update(10, FALSE);
    appGetLattice()->m_pUpdator->Update(20, TRUE);
    appGetLattice()->m_pMeasurements->Report();

    return 0;
}


__REGIST_TEST(TestZ3Symmetry, Misc, TestZ3Symmetry);



//=============================================================================
// END OF FILE
//=============================================================================
