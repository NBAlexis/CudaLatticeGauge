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

    const Real fPlaqutteEneregy = static_cast<Real>(appGetLattice()->m_pGaugeField[0]->CalculatePlaqutteEnergy(F(1.0) / F(3.0)) / (6 * _HC_Volume));

    CFieldGaugeSU3* pStable = dynamic_cast<CFieldGaugeSU3*>(appCreate(_T("CFieldGaugeSU3")));
    CFieldGaugeSU3* pForce = dynamic_cast<CFieldGaugeSU3*>(appCreate(_T("CFieldGaugeSU3")));

    appGetLattice()->m_pGaugeField[0]->CalculateForceAndStaple(pForce, pStable, F(1.0) / F(3.0));
    const Real fPlaqutteEneregy2 = static_cast<Real>(appGetLattice()->m_pGaugeField[0]->CalculatePlaqutteEnergyUsingStable(F(1.0) / F(3.0), pStable) / (6 * _HC_Volume));

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

__REGIST_TEST(TestFileIO, FileIO, TestFileIOBridgePPText, BridgePPText);

__REGIST_TEST(TestFileIO, FileIO, TestFileIOBridgePPBin, BridgePPBin);

UINT TestFileIOWithUpdate(CParameters& sParam)
{
    appGetLattice()->m_pUpdator->Update(150, TRUE);
    return 0;
}

//__REGIST_TEST(TestFileIOWithUpdate, FileIO, TestSaveConfigurationLowMode);

UINT TestFileIOCLG(CParameters& sParam)
{
    UINT uiError = 0;
    const CCString MD51 = appGetLattice()->m_pGaugeField[0]->SaveToFile(_T("testGauge.con"));
    const CCString MD52 = appGetLattice()->GetFieldById(2)->SaveToFile(_T("testFermion.con"));

    CFieldGaugeSU3* pNewGauge = dynamic_cast<CFieldGaugeSU3*>(appCreate(_T("CFieldGaugeSU3")));
    pNewGauge->InitialFieldWithFile(_T("testGauge.con"), EFFT_CLGBin);
    CFieldFermionWilsonSquareSU3* pNewFermion = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appCreate(_T("CFieldFermionWilsonSquareSU3")));
    pNewFermion->InitialFieldWithFile(_T("testFermion.con"), EFFT_CLGBin);

    const CLGComplex res1 = cuCmulf_cr(pNewGauge->DotReal(appGetLattice()->m_pGaugeField[0]), __div(F(1.0), _HC_Volume * _HC_Dir));
    //const CLGComplex res3 = pNewGauge->DotReal(pNewGauge);
    //appGeneral(_T("dot res:%2.20f\n"), res3.x);
    //appGeneral(_T("dot res:%2.20f\n"), res3.y);
    pNewFermion->AxpyMinus(appGetLattice()->GetFieldById(2));
    const CLGComplex res2 = pNewFermion->DotReal(pNewFermion);

    appGeneral(_T("Gauge file test: expeted 3.0 + 0.0i, res = %f %s %f\n"), res1.x, res1.y > 0 ? _T("+") : _T(""), res1.y);
    appGeneral(_T("Fermion file test: expeted 0.0 + 0.0i, res = %f %s %f\n"), res2.x, res2.y > 0 ? _T("+") : _T(""), res2.y);

    appGeneral(_T("Gauge file MD5: %s\n"), MD51.c_str());
    appGeneral(_T("Fermion file MD5: %s\n"), MD52.c_str());

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

UINT TestFileDOUBLE(CParameters&)
{
    const DOUBLE expectedres = 49152.00000000000000000000;
    CCString sFile = _T("testGaugeDouble.con_");
#if !_CLG_DEBUG
    sFile = _T("../Debug/") + sFile;
#endif
    appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFile.c_str(), EFFT_CLGBinDouble);
    const DOUBLE res = appGetLattice()->m_pGaugeField[0]->Dot(appGetLattice()->m_pGaugeField[0]).x;

    appGeneral(_T("load:%2.10f, expect:%2.10f, delta:%2.10f\n"), res, expectedres, appAbs(res - expectedres));
    if (appAbs(res - expectedres) > F(0.000001))
    {
        return 1;
    }
    return 0;
}

UINT TestFileFloat(CParameters&)
{
    const DOUBLE expectedres = 49152.00000000000000000000;
    CCString sFile = _T("testGaugeSingle.con_");
#if !_CLG_DEBUG
    sFile = _T("../Debug/") + sFile;
#endif
    appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFile.c_str(), EFFT_CLGBinFloat);
    const DOUBLE res = appGetLattice()->m_pGaugeField[0]->Dot(appGetLattice()->m_pGaugeField[0]).x;

    appGeneral(_T("load:%2.10f, expect:%2.10f, delta:%2.10f\n"), res, expectedres, appAbs(res - expectedres));
    if (appAbs(res - expectedres) > F(0.000001))
    {
        return 1;
    }
    return 0;
}

UINT TestFileIOCLGCompressed(CParameters& sParam)
{
    UINT uiError = 0;

    appGetLattice()->m_pGaugeField[0]->SaveToCompressedFile(_T("testGaugeCompressed.con"));
    CFieldGaugeSU3* pNewGauge = dynamic_cast<CFieldGaugeSU3*>(appCreate(_T("CFieldGaugeSU3")));
    pNewGauge->InitialFieldWithFile(_T("testGaugeCompressed.con"), EFFT_CLGBinCompressed);

    //appGeneral(_T("=====================\n"));
    //appGetLattice()->m_pGaugeField->DebugPrintMe();
    //appGeneral(_T("=====================\n"));
    //pNewGauge->DebugPrintMe();

    pNewGauge->AxpyMinus(appGetLattice()->m_pGaugeField[0]);

    const CLGComplex res1 = pNewGauge->DotReal(pNewGauge);

    appGeneral(_T("Gauge file test: expeted 0.0 + 0.0i, res = %2.16f %s %2.16f\n"), res1.x, res1.y > 0 ? _T("+") : _T(""), res1.y);

#if _CLG_DOUBLEFLOAT
    if (_cuCabsf(res1) > F(0.00000001))
#else
    if (_cuCabsf(res1) > F(0.02))
#endif
    {
        ++uiError;
    }
    appSafeDelete(pNewGauge);
    return uiError;
}

__REGIST_TEST(TestFileIOCLG, FileIO, TestSaveConfiguration, Save);
#if _CLG_DEBUG
___REGIST_TEST(TestFileIOCLGCompressed, FileIO, TestFileIOCLGCompressedDebug, CLGCompressed, _TEST_DOUBLE);
#else
__REGIST_TEST(TestFileIOCLGCompressed, FileIO, TestFileIOCLGCompressed, CLGCompressed);
#endif

__REGIST_TEST(TestFileDOUBLE, FileIO, TestSaveConfigurationDouble, Double);
__REGIST_TEST(TestFileFloat, FileIO, TestSaveConfigurationFloat, Single);

//=============================================================================
// END OF FILE
//=============================================================================
