//=============================================================================
// FILENAME : CLGTest.cpp
// 
// DESCRIPTION:
//
// REVISION:
//  [12/2/2018 nbale]
//=============================================================================

#include "CLGTest.h"

UINT TestRandom()
{
    //In Bridge++, it use 20 seeds, each for 10 000 000 samples.
    //We are not easy to change the seeds, so we use only 1 seed for 200 000 000 samples
    //1000 blocks, 1024 threads each block, and 200 loops each thread 204 800 000 samples
    TArray<UINT> decomp;
    decomp.AddItem(10);
    decomp.AddItem(10);
    decomp.AddItem(10);

    //fat index is 0 - 20479
    //so we can use a 1024 (max)
    decomp.AddItem(16);
    decomp.AddItem(16);
    decomp.AddItem(4);

    decomp.AddItem(200);

    Real piv = CalculatePi(decomp);
    appGeneral(_T("------- PI result:%f\n"), piv);

    //In bridge++, 1 000 000 samples
    //We use       1 228 800 samples
    decomp.RemoveAll();
    decomp.AddItem(3);
    decomp.AddItem(2);
    decomp.AddItem(2);

    //fat index is 0 - 20479
    //so we can use a 1024 (max)
    decomp.AddItem(16);
    decomp.AddItem(16);
    decomp.AddItem(4);

    decomp.AddItem(100);

    Real ev = CalculateE(decomp);
    appGeneral(_T("------- 1/_sqrt(2) (should be 0.707) result:%f\n"), ev);

    return 0;
}

//=============================================================================
// END OF FILE
//=============================================================================
