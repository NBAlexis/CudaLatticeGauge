//=============================================================================
// FILENAME : Random.h
// 
// DESCRIPTION:
//
//
// REVISION:
//  [12/6/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__global__ 
void _kernalAllocateSeedTable(UINT* pDevicePtr)
{
    intokernal;

    UINT uiSeed = _DC_Seed;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir + 1; ++idir)
        {
            UINT fatIndex = _deviceGetFatIndex(coord, idir);

            CRandomSchrage::_deviceAsignSeeds(pDevicePtr, uiSeed, fatIndex);
        }
    }
}

extern "C" {
    void _callKernelInitialRandomTable(UINT* devicePtr)
    {
        preparethread;
        _kernalAllocateSeedTable << <block, threads >> > (devicePtr);
    }
}

CRandomSchrage::CRandomSchrage(UINT uiSeed)
{
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceSeedTable, sizeof(UINT) * _HC_Volumn * (_HC_Dir + 1)));
    _callKernelInitialRandomTable(m_pDeviceSeedTable);
}

CRandomSchrage::~CRandomSchrage()
{
    checkCudaErrors(cudaFree(m_pDeviceSeedTable));
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
