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

    Real fExpected = F(0.625129946974942);
    sParam.FetchValueReal(_T("ExpectedRes"), fExpected);

    const Real fPlaqutteEneregy = static_cast<Real>(appGetLattice()->m_pGaugeField->CalculatePlaqutteEnergy(F(1.0) / F(3.0)) / (6 * _HC_Volume));

    CFieldGaugeSU3* pStable = dynamic_cast<CFieldGaugeSU3*>(appCreate(_T("CFieldGaugeSU3")));
    CFieldGaugeSU3* pForce = dynamic_cast<CFieldGaugeSU3*>(appCreate(_T("CFieldGaugeSU3")));

    appGetLattice()->m_pGaugeField->CalculateForceAndStaple(pForce, pStable, F(1.0) / F(3.0));
    const Real fPlaqutteEneregy2 = static_cast<Real>(appGetLattice()->m_pGaugeField->CalculatePlaqutteEnergyUsingStable(F(1.0) / F(3.0), pStable) / (6 * _HC_Volume));

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

//__REGIST_TEST(TestFileIOWithUpdate, FileIO, TestSaveConfigurationLowMode);

UINT TestFileIOCLG(CParameters& sParam)
{
    UINT uiError = 0;
    appGetLattice()->m_pGaugeField->SaveToFile(_T("testGauge.con"));
    appGetLattice()->GetFieldById(2)->SaveToFile(_T("testFermion.con"));

    CFieldGaugeSU3* pNewGauge = dynamic_cast<CFieldGaugeSU3*>(appCreate(_T("CFieldGaugeSU3")));
    pNewGauge->InitialFieldWithFile(_T("testGauge.con"), EFFT_CLGBin);
    CFieldFermionWilsonSquareSU3* pNewFermion = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appCreate(_T("CFieldFermionWilsonSquareSU3")));
    pNewFermion->InitialFieldWithFile(_T("testFermion.con"), EFFT_CLGBin);

    const CLGComplex res1 = cuCmulf_cr(pNewGauge->DotReal(appGetLattice()->m_pGaugeField), __div(F(1.0), _HC_Volume * _HC_Dir));
    pNewFermion->AxpyMinus(appGetLattice()->GetFieldById(2));
    const CLGComplex res2 = pNewFermion->DotReal(pNewFermion);

    appGeneral(_T("Gauge file test: expeted 3.0 + 0.0i, res = %f %s %f\n"), res1.x, res1.y > 0 ? _T("+") : _T(""), res1.y);
    appGeneral(_T("Fermion file test: expeted 0.0 + 0.0i, res = %f %s %f\n"), res2.x, res2.y > 0 ? _T("+") : _T(""), res2.y);

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

    appSafeDelete(pNewGauge);
    appSafeDelete(pNewFermion);

    return uiError;
}

UINT TestFileIOCLGCompressed(CParameters& sParam)
{
    UINT uiError = 0;

    appGetLattice()->m_pGaugeField->SaveToCompressedFile(_T("testGaugeCompressed.con"));
    CFieldGaugeSU3* pNewGauge = dynamic_cast<CFieldGaugeSU3*>(appCreate(_T("CFieldGaugeSU3")));
    pNewGauge->InitialFieldWithFile(_T("testGaugeCompressed.con"), EFFT_CLGBinCompressed);

    //appGeneral(_T("=====================\n"));
    //appGetLattice()->m_pGaugeField->DebugPrintMe();

    pNewGauge->AxpyMinus(appGetLattice()->m_pGaugeField);

    const CLGComplex res1 = pNewGauge->DotReal(pNewGauge);

    appGeneral(_T("Gauge file test: expeted 0.0 + 0.0i, res = %2.16f %s %2.16f\n"), res1.x, res1.y > 0 ? _T("+") : _T(""), res1.y);

    if (_cuCabsf(res1) > F(0.00000001))
    {
        ++uiError;
    }
    appSafeDelete(pNewGauge);
    return uiError;
}

__REGIST_TEST(TestFileIOCLG, FileIO, TestSaveConfiguration);
__REGIST_TEST(TestFileIOCLGCompressed, FileIO, TestFileIOCLGCompressed);

//=============================================================================
// END OF FILE
//=============================================================================
