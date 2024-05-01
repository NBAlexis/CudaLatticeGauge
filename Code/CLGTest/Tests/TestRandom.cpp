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
    Real accuracy = F(0.02);
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

//=====================================================================================
// speed: TestRandomSchrage > TestRandomXORWOW ~ TestRandomScrambledSOBOL32
// debug mode: TestRandomSchrage > TestRandomScrambledSOBOL32 > TestRandomXORWOW
//=====================================================================================

__REGIST_TEST(TestRandom, Random, TestRandomSchrage, Schrage);

__REGIST_TEST(TestRandom, Random, TestRandomXORWOW, XORWOW);

//__REGIST_TEST(TestRandom, Random, TestRandomMRG32K3A);

//__REGIST_TEST(TestRandom, Random, TestRandomPhilox);

//__REGIST_TEST(TestRandom, Random, TestRandomSOBOL32);

__REGIST_TEST(TestRandom, Random, TestRandomScrambledSOBOL32, ScrambledSOBOL32);



UINT TestSU3Generator(CParameters& sParam)
{
    //Test Host Random
    UINT uiError = 0;
    INT fieldSampleCount = 1000;
    Real accuracy = 0;
    sParam.FetchValueINT(_T("FieldSampleCount"), fieldSampleCount);
    sParam.FetchValueReal(_T("TestAccuracy"), accuracy);
    CBase* pGaugeField = appCreate(_T("CFieldGaugeSU3"));
    CFieldGauge* pGauge = (NULL != pGaugeField) ? (dynamic_cast<CFieldGauge*>(pGaugeField)) : NULL;
    CBase* pGaugeField2 = appCreate(_T("CFieldGaugeSU3"));
    CFieldGauge* pGauge2 = (NULL != pGaugeField2) ? (dynamic_cast<CFieldGauge*>(pGaugeField2)) : NULL;

    if (NULL == pGauge || NULL == pGauge2)
    {
        return 1;
    }

    CLGComplex average = _make_cuComplex(0, 0);
    Real sigma = 0;
    //size_t freeMemory;
    //size_t totalMemory;
    appPushLogDate(FALSE);
    for (INT i = 0; i < fieldSampleCount; ++i)
    {
        pGauge->InitialField(EFIT_RandomGenerator);
        pGauge2->InitialField(EFIT_SumGenerator);
        average = _cuCaddf(average, pGauge2->DotReal(pGauge));
        sigma += static_cast<Real>(pGauge->CalculateKinematicEnergy());

        if (0 == (i % 50))
        {
            appGeneral(_T("\n"));
        }
        else
        {
            appGeneral(_T("="));
        }
    }
    appPopLogDate();
    average.x = average.x / (_HC_LinkCount * 8);
    average.y = average.y / (_HC_LinkCount * 8);
    sigma = sigma / (_HC_LinkCount * 8);

    appGeneral(_T("\n expected:res = 0+0i:%f+%fi, 0.5:%f\n"), average.x / fieldSampleCount, average.y / fieldSampleCount, sigma / fieldSampleCount);

    if (appAbs(average.x) > accuracy)
    {
        ++uiError;
    }
    if (appAbs(average.y / fieldSampleCount) > accuracy)
    {
        ++uiError;
    }
    if (appAbs(sigma / fieldSampleCount - 0.5) > accuracy)
    {
        ++uiError;
    }

    appSafeDelete(pGauge);
    appSafeDelete(pGauge2);

    return uiError;
}

__REGIST_TEST(TestSU3Generator, Random, TestSU3GeneratorSchrage, GeneratorSchrage);

__REGIST_TEST(TestSU3Generator, Random, TestSU3GeneratorXORWOW, GeneratorXORWOW);

//__REGIST_TEST(TestSU3Generator, Random, TestSU3GeneratorPhilox);

//__REGIST_TEST(TestSU3Generator, Random, TestSU3GeneratorMRG32K3A);

//__REGIST_TEST(TestSU3Generator, Random, TestSU3GeneratorSOBOL32);

__REGIST_TEST(TestSU3Generator, Random, TestSU3GeneratorScrambledSOBOL32, GeneratorScrambledSOBOL32);


//=============================================================================
// END OF FILE
//=============================================================================
