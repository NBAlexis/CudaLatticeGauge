//=============================================================================
// FILENAME : CLGLibMananger.cpp
// 
// DESCRIPTION:
// This is the class for global start-up, control, shut-down
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

extern "C" int GetCudaGPUCount()
{
    int GPU_N;
    checkCudaErrors(cudaGetDeviceCount(&GPU_N));
    return GPU_N;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================