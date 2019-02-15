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


//=============================================================================
// END OF FILE
//=============================================================================
