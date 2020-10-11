//=============================================================================
// FILENAME : TestSU3Generator.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [01/30/2019 nbale]
//=============================================================================

#include "CLGTest.h"

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

    CLGComplex average = _make_cuComplex(0,0);
    Real sigma = 0;
    //size_t freeMemory;
    //size_t totalMemory;
    appSetLogDate(FALSE);
    for (INT i = 0; i < fieldSampleCount; ++i)
    {
        pGauge->InitialField(EFIT_RandomGenerator);
        pGauge2->InitialField(EFIT_SumGenerator);
        average = _cuCaddf(average, pGauge2->DotReal(pGauge));
        sigma += pGauge->CalculateKinematicEnergy();

        if (0 == (i % 50))
        {
            appGeneral(_T("\n"));
        }
        else
        {
            appGeneral(_T("="));
        }
    }
    appSetLogDate(TRUE);
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

__REGIST_TEST(TestSU3Generator, Random, TestSU3GeneratorSchrage);

__REGIST_TEST(TestSU3Generator, Random, TestSU3GeneratorXORWOW);

__REGIST_TEST(TestSU3Generator, Random, TestSU3GeneratorPhilox);

__REGIST_TEST(TestSU3Generator, Random, TestSU3GeneratorMRG32K3A);

__REGIST_TEST(TestSU3Generator, Random, TestSU3GeneratorSOBOL32);

__REGIST_TEST(TestSU3Generator, Random, TestSU3GeneratorScrambledSOBOL32);


//=============================================================================
// END OF FILE
//=============================================================================
