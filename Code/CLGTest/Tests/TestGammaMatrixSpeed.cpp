//=============================================================================
// FILENAME : TestOperators.cpp
// 
// DESCRIPTION:
//
//     Test the operations on fields
//
// REVISION:
//  [02/10/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestGammaMatrix(CParameters& )
{
    CFieldFermionWilsonSquareSU3* pFermion = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));

    for (UINT i = 0; i < 1000; ++i)
    {
        for (UINT j = 0; j < EGM_MAX; ++j)
        {
            pFermion->ApplyGamma((EGammaMatrix)j);
        }
    }

    return 0;
}

__REGIST_TEST(TestGammaMatrix, Misc, TestGammaMatrixSpeed);


//=============================================================================
// END OF FILE
//=============================================================================
