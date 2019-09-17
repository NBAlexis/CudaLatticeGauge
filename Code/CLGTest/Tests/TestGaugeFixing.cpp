//=============================================================================
// FILENAME : TestZ3Symmetry.cpp
// 
// DESCRIPTION:
//
//     Test the Z3 Symmetry
//
// REVISION:
//  [06/23/2019 nbale]
//=============================================================================

#include "CLGTest.h"

#define _tfftMX 9
#define _tfftMY 10
#define _tfftMZ 11
#define _tfftMT 12

UINT TestFFT(CParameters&)
{
    CCLGFFTHelper::TestFFT();
    return 0;
}


__REGIST_TEST(TestFFT, Misc, TestFFT);



//=============================================================================
// END OF FILE
//=============================================================================
