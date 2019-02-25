//=============================================================================
// FILENAME : TestGaugeSmearing.cpp
// 
// DESCRIPTION:
//
//     
//
// REVISION:
//  [02/24/2019 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestGaugeSmearing(CParameters& sParam)
{
    //CMeasureMesonCorrelator* pMeasure = dynamic_cast<CMeasureMesonCorrelator*>(appGetLattice()->m_pMeasurements->GetMeasureById(1));
    //if (NULL == pMeasure)
    //{
    //    return 1;
    //}
    //Real fExpected = F(0.625);
    //sParam.FetchValueReal(_T("ExpectedRes"), fExpected);

    //appGetLattice()->m_pUpdator->Update(10, FALSE);
    //appGetLattice()->m_pUpdator->Update(20, TRUE);

    //Real fRes = pMeasure->m_lstResults[0][0];
    //appGeneral(_T("res : expected=%f res=%f"), fExpected, fRes);
    //if (appAbs(fRes - fExpected) > F(0.01))
    //{
    //    return 1;
    //}

    return 0;
}

//930 - 990 ms
//__REGIST_TEST(TestGaugeSmearing, Updator, TestGaugeSmearingAPEProj);

//__REGIST_TEST(TestGaugeSmearing, Updator, TestGaugeSmearingAPEStout);


//=============================================================================
// END OF FILE
//=============================================================================
