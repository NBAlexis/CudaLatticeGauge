//=============================================================================
// FILENAME : TestConfigurationFileIO.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [01/31/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestFileIO(CParameters& sParam)
{
    //appGeneral(_T("size of sindex = %d\n"), sizeof(SIndex));
    //CCudaHelper::DebugFunction();

    Real fExpected = F(0.625129946974942);
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);

    Real fPlaqutteEneregy = appGetLattice()->m_pGaugeField->CalculatePlaqutteEnergy(F(1.0) / F(3.0)) / (6 * _HC_Volumn);

    CFieldGaugeSU3* pStable = dynamic_cast<CFieldGaugeSU3*>(appCreate(_T("CFieldGaugeSU3")));
    CFieldGaugeSU3* pForce = dynamic_cast<CFieldGaugeSU3*>(appCreate(_T("CFieldGaugeSU3")));

    appGetLattice()->m_pGaugeField->CalculateForceAndStaple(pForce, pStable, F(1.0) / F(3.0));
    Real fPlaqutteEneregy2 = appGetLattice()->m_pGaugeField->CalculatePlaqutteEnergyUsingStable(F(1.0) / F(3.0), pStable) / (6 * _HC_Volumn);

    appGeneral(_T("Plaqutte Energy (expected:0.625129946974942)= %1.10f and %1.10f\n"), F(1.0) - fPlaqutteEneregy, F(1.0) - fPlaqutteEneregy2);

    UINT uiError = 0;
    if (appAbs(F(1.0) - fPlaqutteEneregy - F(0.625129946974942)) > F(0.000001))
    {
        ++uiError;
    }
    if (appAbs(F(1.0) - fPlaqutteEneregy2 - F(0.625129946974942)) > F(0.000001))
    {
        ++uiError;
    }
    return uiError;
}

__REGIST_TEST(TestFileIO, FileIO, TestFileIOBridgePPText);

__REGIST_TEST(TestFileIO, FileIO, TestFileIOBridgePPBin);

UINT TestFileIOWithUpdate(CParameters& sParam)
{
    appGetLattice()->m_pUpdator->Update(150, TRUE);
    return 0;
}

__REGIST_TEST(TestFileIOWithUpdate, FileIO, TestSaveConfigurationLowMode);

UINT TestFileIOCLG(CParameters& sParam)
{
    UINT uiError = 0;
    appGetLattice()->m_pGaugeField->SaveToFile(_T("testGauge.con"));
    appGetLattice()->GetFieldById(2)->SaveToFile(_T("testFermion.con"));

    CFieldGaugeSU3* pNewGauge = dynamic_cast<CFieldGaugeSU3*>(appCreate(_T("CFieldGaugeSU3")));
    pNewGauge->InitialFieldWithFile(_T("testGauge.con"), EFFT_CLGBin);
    CFieldFermionWilsonSquareSU3* pNewFermion = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appCreate(_T("CFieldFermionWilsonSquareSU3")));
    pNewFermion->InitialFieldWithFile(_T("testFermion.con"), EFFT_CLGBin);

    _Complex res1 = cuCmulf_cr(pNewGauge->Dot(appGetLattice()->m_pGaugeField), __div(F(1.0), _HC_Volumn * _HC_Dir));
    pNewFermion->AxpyMinus(appGetLattice()->GetFieldById(2));
    _Complex res2 = pNewFermion->Dot(pNewFermion);

    appGeneral(_T("Gauge file test: expeted 3.0 + 0.0i, res = %f%s%f"), res1.x, res1.y > 0 ? _T("+") : _T(""), res1.y);
    appGeneral(_T("Fermion file test: expeted 0.0 + 0.0i, res = %f%s%f"), res2.x, res2.y > 0 ? _T("+") : _T(""), res2.y);

    if (appAbs(res1.x - F(3.0)) > F(0.000001)
     || appAbs(res1.y) > F(0.000001))
    {
        ++uiError;
    }

    if (appAbs(res2.x) > F(0.000001)
     || appAbs(res2.y) > F(0.000001))
    {
        ++uiError;
    }

    return uiError;
}

__REGIST_TEST(TestFileIOCLG, FileIO, TestSaveConfiguration);

//=============================================================================
// END OF FILE
//=============================================================================
