//=============================================================================
// FILENAME : TestCommon.cpp
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
    //CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField[0]);
    CFieldFermionWilsonSquareSU3* pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    CFieldFermionWilsonSquareSU3* pF2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    CFieldFermionWilsonSquareSU3* pF3 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    CFieldFermionWilsonSquareSU3* pF4 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    pF1->InitialField(EFIT_RandomGaussian);
    pF2->InitialField(EFIT_RandomGaussian);

    const CLGComplex dot1 = pF1->DotReal(pF1);
    const CLGComplex dot2 = pF2->DotReal(pF2);
    const CLGComplex dot3 = pF1->DotReal(pF2);
    const CLGComplex dot4 = pF2->DotReal(pF1);
    appGeneral(_T("Estimator f1.f1 = (%f %f); f2.f2 = (%f %f); f1.f2 = (%f %f); f2.f1 = (%f %f)\n"),
        dot1.x / (12 * _HC_Volume), dot1.y, dot2.x / (12 * _HC_Volume), dot2.y, dot3.x / (12 * _HC_Volume), dot3.y / (12 * _HC_Volume), dot4.x / (12 * _HC_Volume), dot4.y / (12 * _HC_Volume));

    if (appAbs(1.0f - dot1.x / (12 * _HC_Volume)) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(dot1.y) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(1.0f - dot2.x / (12 * _HC_Volume)) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(dot2.y) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(dot3.x / (12 * _HC_Volume)) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(dot3.y / (12 * _HC_Volume)) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(dot4.x / (12 * _HC_Volume)) > 0.01f)
    {
        ++uiErrors;
    }
    if (appAbs(dot4.y / (12 * _HC_Volume)) > 0.01f)
    {
        ++uiErrors;
    }

    pF1->CopyTo(pF3);
    pF3->D(_FIELDS);
    pF2->CopyTo(pF4);
    pF4->Ddagger(_FIELDS);
    const CLGComplex dot5 = pF2->DotReal(pF3);
    const CLGComplex dot6 = pF4->DotReal(pF1);
    appGeneral(_T("DDagger f2.(D.f1) = (%f %f); (D+.f2).f1 = (%f %f);\n"),
        dot5.x, dot5.y, dot6.x, dot6.y);

    if (appAbs(dot5.x - dot6.x) > F(0.02))
    {
        ++uiErrors;
    }
    if (appAbs(dot5.y - dot6.y)  > F(0.01))
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

UINT TestALogDefinition(CParameters&)
{
    //create a random field
    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField[0]);

    CFieldGaugeSU3* pTestGauge = dynamic_cast<CFieldGaugeSU3*>(pGauge->GetCopy());
    pGauge->TransformToIA();
    pGauge->TransformToU();
    pTestGauge->AxpyMinus(pGauge);

    const CLGComplex testDot = pTestGauge->DotReal(pTestGauge);
    appGeneral(_T("Dot result = %f + %f I\n"), testDot.x, testDot.y);
    if (_cuCabsf(testDot) > F(0.000000001))
    {
        return 1;
    }
    return 0;
}

__REGIST_TEST(TestSmallMatrix, Verify, TestSmallMatrix, SmallMatrix);

__REGIST_TEST(TestOperators, Misc, TestOperators, Operators);

//__REGIST_TEST(TestQuickAxpy, Misc, TestQuickAxpy);

__REGIST_TEST(TestALogDefinition, Misc, TestALogDefinition, ALog);

//__REGIST_TEST(TestDirichletDOperator, Misc, TestRotationOperator);

UINT TestGammaMatrix(CParameters&)
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

    CFieldGaugeSU3* pGauge = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pGaugeField[0]);
    CFieldFermionWilsonSquareSU3* pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(2));
    UINT uiErrors = pF1->TestGamma5Hermitian(pGauge, 0 != iG5);
    //CCudaHelper::DebugFunction();
    pF1->Return();

    return uiErrors;
}

UINT TestAnitiHermiticity(CParameters&)
{
    //test Ddagger
    //CFieldGauge* pGauge = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField);
    CFieldFermionKS* pF1 = dynamic_cast<CFieldFermionKS*>(appGetLattice()->GetPooledFieldById(2));
    UINT uiErrors = pF1->TestAntiHermitian(_FIELDS);
    appGeneral(_T("=== Tested Fermion: \n %s \n"), pF1->GetInfos(_T("     ")).c_str());
    pF1->Return();

    return uiErrors;
}

UINT TestBosonHermiticity(CParameters&)
{
    //test Ddagger
    //CFieldGauge* pGauge = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField);
    CFieldBoson* pF1 = dynamic_cast<CFieldBoson*>(appGetLattice()->GetPooledFieldById(2));
    UINT uiErrors = pF1->CheckHermitian(_FIELDS);
    appGeneral(_T("=== Tested Boson: \n %s \n"), pF1->GetInfos(_T("     ")).c_str());
    pF1->Return();

    return uiErrors;
}

//930 - 990 ms
//__REGIST_TEST(TestGammaMatrix, Misc, TestGammaMatrixSpeed);

UINT TestDebugFunction(CParameters&)
{
    CCudaHelper::DebugFunction();
    //CIndexData::DebugLinkDirichletOrDagger(1);
    return 0;
}

___REGIST_TEST(TestGamma5Hermiticity, Tools, TestGamm5Hermiticity, Gamm5Hermiticity, _TEST_NOCHECK);

___REGIST_TEST(TestAnitiHermiticity, Tools, TestAnitiHermiticity, AnitiHermiticity, _TEST_NOCHECK);

___REGIST_TEST(TestBosonHermiticity, Tools, TestBosonHermiticity, BosonHermiticity, _TEST_NOCHECK);

___REGIST_TEST(TestDebugFunction, Tools, TestDebug, Debug, _TEST_NOCHECK);

UINT TestGaugeInvarience(CParameters&)
{
    UINT uiError = 0;
    TArray<Real> beforeGaugeTransform;
    //TArray<CFieldGauge*> gaugefields;
    //gaugefields.AddItem(appGetLattice()->m_pGaugeField);

    for (INT i = 0; i < appGetLattice()->m_pActionList.Num(); ++i)
    {
        Real fEnergy = static_cast<Real>(appGetLattice()->GetActionById(static_cast<BYTE>(i + 1))->Energy(FALSE, _FIELDS, NULL));
        beforeGaugeTransform.AddItem(fEnergy);
    }

    CGaugeFixingRandom* pRandom = dynamic_cast<CGaugeFixingRandom*>(appGetLattice()->m_pGaugeFixing);

    //make sure the gauge field is really changed
    CFieldGauge* pGaugeCopy = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField[0]->GetCopy());
    pRandom->GaugeFixing(appGetLattice()->m_pGaugeField[0]);
    pGaugeCopy->AxpyMinus(appGetLattice()->m_pGaugeField[0]);
    CMeasure::LogGeneralComplex(pGaugeCopy->Dot(pGaugeCopy));

    appGetLattice()->m_pGaugeField[0]->CopyTo(pGaugeCopy);
    pGaugeCopy->AxpyMinus(appGetLattice()->m_pGaugeField[0]);
    CMeasure::LogGeneralComplex(pGaugeCopy->Dot(pGaugeCopy));
    
    for (INT j = 0; j < appGetLattice()->m_pFermionField.Num(); ++j)
    {
        pRandom->AlsoFixingFermion(dynamic_cast<CFieldFermion*>(appGetLattice()->GetFieldById(static_cast<BYTE>(j + 2))));
    }

    for (INT i = 0; i < appGetLattice()->m_pActionList.Num(); ++i)
    {
        Real fEnergy = static_cast<Real>(appGetLattice()->GetActionById(static_cast<BYTE>(i + 1))->Energy(FALSE, _FIELDS, NULL));
        appGeneral(_T("Action%d, Before:%2.20f, After:%2.20f\n"), i, beforeGaugeTransform[i], fEnergy);
        if (appAbs(beforeGaugeTransform[i] - fEnergy) > F(0.0000001))
        {
            ++uiError;
        }
    }

    return 0;
}

__REGIST_TEST(TestGaugeInvarience, Misc, TestGaugeInvarience, GaugeInvarience);

UINT TestBackgroundField(CParameters&)
{
    UINT uiError = 0;
    CFieldGaugeU1Real* pU1 = dynamic_cast<CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(2));

    TArray<BYTE> xyz;
    xyz.AddItem(0);
    xyz.AddItem(1);
    xyz.AddItem(2);

    TArray<BYTE> t;
    t.AddItem(3);

    TArray<BYTE> xy;
    xy.AddItem(0);
    xy.AddItem(1);

    TArray<BYTE> zt;
    zt.AddItem(2);
    zt.AddItem(3);

    TArray<BYTE> x;
    x.AddItem(0);

    TArray<BYTE> y;
    y.AddItem(1);

    TArray<BYTE> xzt;
    xzt.AddItem(0);
    xzt.AddItem(2);
    xzt.AddItem(3);

    TArray<BYTE> yzt;
    yzt.AddItem(1);
    yzt.AddItem(2);
    yzt.AddItem(3);

    appGeneral(_T("============ chemical potential ===========\n"));

    Real fCheck = pU1->CheckSliceSame(0, 1);
    appGeneral(_T("xy slice same? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    fCheck = pU1->CheckSliceSame(2, 3);
    appGeneral(_T("zt slice same? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    fCheck = pU1->CheckZero(0, 1, xyz);
    appGeneral(_T("xy slice[xyz] zero? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    fCheck = pU1->CheckZero(2, 3, xyz);
    appGeneral(_T("zt slice[xyz] zero? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    appGeneral(_T("\n\n"));
    pU1->DebugPrintSlice(0, 1, t);
    appGeneral(_T("\n\n"));

    appGeneral(_T("============ EZ 0 ===========\n"));
    pU1->InitialU1Real(EURT_None, EURT_E_t, EURT_None, F(0.0), F(0.1), F(0.0), TRUE);
    
    fCheck = pU1->CheckSliceSame(2, 3);
    appGeneral(_T("zt slice same? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    fCheck = pU1->CheckZero(2, 3, xy);
    appGeneral(_T("zt slice[xy] zero? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    appGeneral(_T("\n\n"));
    pU1->DebugPrintSlice(2, 3, zt);
    appGeneral(_T("\n\n"));

    appGeneral(_T("============ EZ 1 ===========\n"));
    pU1->InitialU1Real(EURT_None, EURT_E_z, EURT_None, F(0.0), F(0.1), F(0.0), TRUE);

    fCheck = pU1->CheckSliceSame(2, 3);
    appGeneral(_T("zt slice same? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    fCheck = pU1->CheckZero(2, 3, xy);
    appGeneral(_T("zt slice[xy] zero? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    appGeneral(_T("\n\n"));
    pU1->DebugPrintSlice(2, 3, zt);
    appGeneral(_T("\n\n"));

    appGeneral(_T("============ BZ 0 ===========\n"));
    pU1->InitialU1Real(EURT_None, EURT_None, EURT_Bp_y, F(0.0), F(0.0), F(0.1), TRUE);

    fCheck = pU1->CheckSliceSame(0, 1);
    appGeneral(_T("xy slice same? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    fCheck = pU1->CheckZero(0, 1, zt);
    appGeneral(_T("xy slice[zt] zero? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    appGeneral(_T("\n\n"));
    pU1->DebugPrintSlice(0, 1, xy);
    appGeneral(_T("\n\n"));

    appGeneral(_T("============ BZ 0 no twist ===========\n"));
    pU1->InitialU1Real(EURT_None, EURT_None, EURT_Bp_y_notwist, F(0.0), F(0.0), F(0.1), TRUE);

    fCheck = pU1->CheckSliceSame(0, 1);
    appGeneral(_T("xy slice same? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    fCheck = pU1->CheckZero(0, 1, yzt);
    appGeneral(_T("xy slice[yzt] zero? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    appGeneral(_T("\n\n"));
    pU1->DebugPrintSlice(0, 1, x);
    appGeneral(_T("\n\n"));

    appGeneral(_T("============ BZ 1 ===========\n"));
    pU1->InitialU1Real(EURT_None, EURT_None, EURT_Bp_x, F(0.0), F(0.0), F(0.1), TRUE);

    fCheck = pU1->CheckSliceSame(0, 1);
    appGeneral(_T("xy slice same? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    fCheck = pU1->CheckZero(0, 1, zt);
    appGeneral(_T("xy slice[zt] zero? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    appGeneral(_T("\n\n"));
    pU1->DebugPrintSlice(0, 1, xy);
    appGeneral(_T("\n\n"));

    appGeneral(_T("============ BZ 1 no twist ===========\n"));
    pU1->InitialU1Real(EURT_None, EURT_None, EURT_Bp_x_notwist, F(0.0), F(0.0), F(0.1), TRUE);

    fCheck = pU1->CheckSliceSame(0, 1);
    appGeneral(_T("xy slice same? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    fCheck = pU1->CheckZero(0, 1, xzt);
    appGeneral(_T("xy slice[xzt] zero? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    appGeneral(_T("\n\n"));
    pU1->DebugPrintSlice(0, 1, y);
    appGeneral(_T("\n\n"));

    appGeneral(_T("============ BZ 2 ===========\n"));
    pU1->InitialU1Real(EURT_None, EURT_None, EURT_Bp_xy, F(0.0), F(0.0), F(0.1), TRUE);

    fCheck = pU1->CheckSliceSame(0, 1);
    appGeneral(_T("xy slice same? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    fCheck = pU1->CheckZero(0, 1, zt);
    appGeneral(_T("xy slice[zt] zero? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    appGeneral(_T("\n\n"));
    pU1->DebugPrintSlice(0, 1, xy);
    appGeneral(_T("\n\n"));

    appGeneral(_T("============ BZ 2 no twist ===========\n"));
    pU1->InitialU1Real(EURT_None, EURT_None, EURT_Bp_xy_notwist, F(0.0), F(0.0), F(0.1), TRUE);

    fCheck = pU1->CheckSliceSame(0, 1);
    appGeneral(_T("xy slice same? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    fCheck = pU1->CheckZero(0, 1, zt);
    appGeneral(_T("xy slice[zt] zero? %f\n"), fCheck);
    if (fCheck > F(0.001))
    {
        ++uiError;
    }

    appGeneral(_T("\n\n"));
    pU1->DebugPrintSlice(0, 1, xy);
    appGeneral(_T("\n\n"));

    return 0;
}

___REGIST_TEST(TestBackgroundField, Tools, TestBackgroundField, PrintBackgroundField, _TEST_NOCHECK);

UINT TestPlaqutteTable(CParameters&)
{
    //appGetLattice()->m_pIndexCache->DebugEdgeGlue(1, SSmallInt4(-1, -1, -1, -1));
    appGetLattice()->m_pIndexCache->DebugStapleTable();
    return 0;
}

___REGIST_TEST(TestPlaqutteTable, Tools, TestPlaqutteTable, PlaqutteTable, _TEST_NOCHECK);



//=============================================================================
// END OF FILE
//=============================================================================
