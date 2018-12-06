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
void _kernalAllocateSeedTable(UINT uiSeed, CDeviceLattice* pLattice, CRandomSchrage* pRandom)
{
    intokernal;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT fatIndex = _deviceGetFatIndex(pLattice, coord, idir + 1);

            pRandom->_deviceAsignSeeds(uiSeed, fatIndex);
        }
    }
}

CRandomSchrage::CRandomSchrage(UINT uiSeed, CDeviceLattice* pDeviceLattice)
    : m_pOwner(pDeviceLattice)
{
    preparethread;
    checkCudaErrors(cudaMalloc((void **)&m_pDeviceSeedTable, sizeof(UINT) * pLattice->m_uiVolumn * (pLattice->m_uiDir + 1)));

    _kernalAllocateSeedTable << <block, threads>> > (uiSeed, pDeviceLattice, this);
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
