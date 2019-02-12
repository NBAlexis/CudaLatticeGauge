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
    appGetLattice()->m_pUpdator->Update(5, TRUE);
    //CCudaHelper::DebugFunction();
    return 0;
}

__REGIST_TEST(TestFermionUpdator, Updator, TestFermionUpdator);


//=============================================================================
// END OF FILE
//=============================================================================
