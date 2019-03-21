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

UINT TestOperators(CParameters& )
{
    UINT uiErrors = 0;

    //test Ddagger
    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField);
    CFieldFermionWilsonSquareSU3* pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    CFieldFermionWilsonSquareSU3* pF2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    CFieldFermionWilsonSquareSU3* pF3 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    CFieldFermionWilsonSquareSU3* pF4 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    pF1->InitialField(EFIT_RandomGaussian);
    pF2->InitialField(EFIT_RandomGaussian);

    CLGComplex dot1 = pF1->Dot(pF1);
    CLGComplex dot2 = pF2->Dot(pF2);
    CLGComplex dot3 = pF1->Dot(pF2);
    CLGComplex dot4 = pF2->Dot(pF1);
    appGeneral(_T("Estimator f1.f1 = (%f %f); f2.f2 = (%f %f); f1.f2 = (%f %f); f2.f1 = (%f %f)\n"),
        dot1.x / (12 * _HC_Volumn), dot1.y, dot2.x / (12 * _HC_Volumn), dot2.y, dot3.x / (12 * _HC_Volumn), dot3.y / (12 * _HC_Volumn), dot4.x / (12 * _HC_Volumn), dot4.y / (12 * _HC_Volumn));

    if (appAbs(1.0f - dot1.x / (12 * _HC_Volumn)) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(dot1.y) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(1.0f - dot2.x / (12 * _HC_Volumn)) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(dot2.y) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(dot3.x / (12 * _HC_Volumn)) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(dot3.y / (12 * _HC_Volumn)) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(dot4.x / (12 * _HC_Volumn)) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(dot4.y / (12 * _HC_Volumn)) > 0.01f)
    {
        ++uiErrors;
    }

    pF1->CopyTo(pF3);
    pF3->D(pGauge);
    pF2->CopyTo(pF4);
    pF4->Ddagger(pGauge);
    CLGComplex dot5 = pF2->Dot(pF3);
    CLGComplex dot6 = pF4->Dot(pF1);
    appGeneral(_T("DDagger f2.(D.f1) = (%f %f); (D+.f2).f1 = (%f %f);\n"),
        dot5.x, dot5.y, dot6.x, dot6.y);

    if (appAbs(dot5.x - dot6.x) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(dot5.y - dot6.y)  > 0.01f)
    {
        ++uiErrors;
    }
    //CCudaHelper::DebugFunction();

    return uiErrors;
}


UINT TestSmallMatrix(CParameters&)
{
    CLinearAlgebraHelper::TestSmallMatrix();
    return 0;
}

UINT TestQuickAxpy(CParameters&)
{
    CFieldFermionWilsonSquareSU3* pF0 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    CFieldFermionWilsonSquareSU3* pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pF0->GetCopy());
    CFieldFermionWilsonSquareSU3* pF2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pF0->GetCopy());
    CFieldFermionWilsonSquareSU3* pF3 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pF0->GetCopy());
    CFieldFermionWilsonSquareSU3* pF4 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pF0->GetCopy());

    CTimer timer1;
    timer1.Start();
    for (UINT i = 0; i < 20; ++i)
    {
        pF2->ScalarMultply(F(2.0)); //f2 = 2f0
        pF1->AxpyPlus(pF2); //f1 = f1+f2 = 3f0
        pF2->ScalarMultply(F(0.8)); //f2 = 0.8f2 = 1.6f0
        pF1->AxpyMinus(pF2); //f1 = f1-f2 = 3f0-1.6f0=1.4f0
        pF1->Axpy(F(-0.5), pF2); //f1 = f1 - 0.5f2 = 1.4f0 - 0.8f0 = 0.6f0
        pF1->Axpy(_make_cuComplex(F(0.25), F(0.0)), pF2); //f1 = f1+0.25f2=0.6f0+0.4f0 = f0
        pF1->CopyTo(pF2);
    }
    pF1->AxpyMinus(pF0);
    CLGComplex res = pF1->Dot(pF1);
    timer1.Stop();
    appGeneral(_T("res = %f %f  t=%f (ms)\n"), res.x, res.y, timer1.Elapsed());

    CTimer timer2;
    timer2.Start();
    for (UINT i = 0; i < 20; ++i)
    {
        pF4->ScalarMultply1(F(2.0)); //f2 = 2f0
        pF3->AxpyPlus1(pF4); //f1 = f1+f2 = 3f0
        pF4->ScalarMultply1(F(0.8)); //f2 = 0.8f2 = 1.6f0
        pF3->AxpyMinus1(pF4); //f1 = f1-f2 = 3f0-1.6f0=1.4f0
        pF3->Axpy1(F(-0.5), pF4); //f1 = f1 - 0.5f2 = 1.4f0 - 0.8f0 = 0.6f0
        pF3->Axpy1(_make_cuComplex(F(0.25), F(0.0)), pF4); //f1 = f1+0.25f2=0.6f0+0.4f0 = f0
        pF3->CopyTo(pF4);
    }
    pF3->AxpyMinus1(pF0);
    res = pF3->Dot1(pF3);
    timer2.Stop();
    appGeneral(_T("res = %f %f  t=%f (ms)\n"), res.x, res.y, timer2.Elapsed());
    return 0;
}

__REGIST_TEST(TestSmallMatrix, Misc, TestSmallMatrix);

__REGIST_TEST(TestOperators, Misc, TestOperators);

__REGIST_TEST(TestQuickAxpy, Misc, TestQuickAxpy);


//=============================================================================
// END OF FILE
//=============================================================================
