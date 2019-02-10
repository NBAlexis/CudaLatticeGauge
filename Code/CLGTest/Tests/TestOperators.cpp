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

UINT TestOperators(CParameters& sParam)
{
    //appGetLattice()->m_pUpdator->Update(2, TRUE);
    //CCudaHelper::DebugFunction();

    //test Ddagger
    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField);
    CFieldFermionWilsonSquareSU3* pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    CFieldFermionWilsonSquareSU3* pF2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    CFieldFermionWilsonSquareSU3* pF3 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    CFieldFermionWilsonSquareSU3* pF4 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    pF1->InitialField(EFIT_RandomGaussian);
    pF2->InitialField(EFIT_RandomGaussian);

    _Complex dot1 = pF1->Dot(pF1);
    _Complex dot2 = pF2->Dot(pF2);
    _Complex dot3 = pF1->Dot(pF2);
    _Complex dot4 = pF2->Dot(pF1);
    appGeneral(_T("Estimator f1.f1 = (%f %f); f2.f2 = (%f %f); f1.f2 = (%f %f); f2.f1 = (%f %f)\n"),
        dot1.x / (12 * _HC_Volumn), dot1.y, dot2.x / (12 * _HC_Volumn), dot2.y, dot3.x / (12 * _HC_Volumn), dot3.y / (12 * _HC_Volumn), dot4.x / (12 * _HC_Volumn), dot4.y / (12 * _HC_Volumn));

    pF1->CopyTo(pF3);
    pF3->D(pGauge);
    pF2->CopyTo(pF4);
    pF4->Ddagger(pGauge);
    _Complex dot5 = pF2->Dot(pF3);
    _Complex dot6 = pF4->Dot(pF1);
    appGeneral(_T("DDagger f2.(D.f1) = (%f %f); (D+.f2).f1 = (%f %f);\n"),
        dot5.x, dot5.y, dot6.x, dot6.y);

    return 0;
}

__REGIST_TEST(TestOperators, Field, TestOperators);


//=============================================================================
// END OF FILE
//=============================================================================
