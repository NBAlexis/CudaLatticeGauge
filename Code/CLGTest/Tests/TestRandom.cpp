//=============================================================================
// FILENAME : TestRandom.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestRandom(CParameters& sParam)
{
    //Test Host Random
    Real hostRandom = F(0.0);
    INT hostSampleCount = 1000;
    sParam.FetchValueINT(_T("HostSample"), hostSampleCount);
    for (INT i = 0; i < hostSampleCount; ++i)
    {
        hostRandom = hostRandom + GetRandomReal();
    }
    appGeneral(_T("------- Host random result:%f\n"), hostRandom / hostSampleCount);

    UINT uiError = 0;
    Real accuracy = F(0.001);
    TArray<UINT> decompPi;
    TArray<UINT> decompGaussian; 

    sParam.FetchValueArrayUINT(_T("PiDecomp"), decompPi);
    sParam.FetchValueArrayUINT(_T("GaussianDecomp"), decompGaussian);
    sParam.FetchValueReal(_T("TestAccuracy"), accuracy);

    const Real piv = CalculatePi(decompPi);
    appGeneral(_T("------- PI result:%f\n"), piv);

    const Real ev = CalculateE(decompGaussian);
    appGeneral(_T("------- 1/_sqrt(2) (should be 0.707) result:%f\n"), ev);

    if (appAbs(hostRandom / hostSampleCount - F(0.5)) > accuracy * F(50.0))
    {
        ++uiError;
    }
    if (appAbs(piv - PI) > accuracy)
    {
        ++uiError;
    }
    if (appAbs(ev - InvSqrt2) > accuracy)
    {
        ++uiError;
    }
    return uiError;
}

__REGIST_TEST(TestRandom, Random, TestRandomSchrage);

__REGIST_TEST(TestRandom, Random, TestRandomXORWOW);

__REGIST_TEST(TestRandom, Random, TestRandomMRG32K3A);

__REGIST_TEST(TestRandom, Random, TestRandomPhilox);

__REGIST_TEST(TestRandom, Random, TestRandomSOBOL32);

__REGIST_TEST(TestRandom, Random, TestRandomScrambledSOBOL32);

//=============================================================================
// END OF FILE
//=============================================================================
