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

UINT TestGamma5Hermiticity(CParameters& param)
{
    //test Ddagger
    INT iG5 = 0;
    param.FetchValueINT(_T("GAMM5Test"), iG5);

    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField);
    CFieldFermionWilsonSquareSU3* pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    UINT uiErrors = pF1->TestGamma5Hermitian(pGauge, 0 != iG5);
    //CCudaHelper::DebugFunction();
    pF1->Return();

    return uiErrors;
}

UINT TestAnitiHermiticity(CParameters& )
{
    //test Ddagger
    appGeneral(_T("omega?:%f\n"), CCommonData::m_fOmega);
    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField);
    CFieldFermionKSSU3* pF1 = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(2));
    UINT uiErrors = pF1->TestAntiHermitian(pGauge);
    //CCudaHelper::DebugFunction();
    pF1->Return();

    return uiErrors;
}

//930 - 990 ms
//__REGIST_TEST(TestGammaMatrix, Misc, TestGammaMatrixSpeed);

UINT TestDebugFunction(CParameters&)
{
    CCudaHelper::DebugFunction();
    return 0;
}

__REGIST_TEST(TestGamma5Hermiticity, Misc, TestGamm5Hermiticity);

__REGIST_TEST(TestAnitiHermiticity, Misc, TestAnitiHermiticity);

__REGIST_TEST(TestDebugFunction, Misc, TestDebug);

//=============================================================================
// END OF FILE
//=============================================================================
