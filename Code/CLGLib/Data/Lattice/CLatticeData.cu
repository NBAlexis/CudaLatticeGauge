//=============================================================================
// FILENAME : CLatticeData.cpp
// 
// DESCRIPTION:
// This is the class for the lattce data
// NOTE:: We only have 4D case, 3D = 1xLxLxL, and 2D= 1x1xLxL
// REVISION:
//  [12/4/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE


__global__ void _kernelDeletePtrs(CIndex * pdeviceIndex)
{
    if (NULL != pdeviceIndex)
    {
        delete pdeviceIndex;
    }
}

extern "C" {
    void _cDeletePtrs(CIndex * pdeviceIndex)
    {
        _kernelDeletePtrs << <1, 1 >> > (pdeviceIndex);
    }
}

/**
* m_uiLatticeDecompose[0,1,2] is the blocks
* m_uiLatticeDecompose[3,4,5] is the threads in blocks
*/
CLatticeData::CLatticeData()
    : m_pRandom(NULL)
    , m_pRandomSchrage(NULL)
    , m_pGaugeField(NULL)
    , m_pGaugeFieldStaple(NULL)
    , m_pUpdator(NULL)

    , m_pDeviceRandom(NULL)
    , m_pDeviceRandomSchrage(NULL)
    , m_pDeviceGaugeField(NULL)
    , m_pDeviceGaugeFieldStaple(NULL)
    , m_pDeviceIndex(NULL)
    //, m_pDeviceActionList(NULL)
{
    
}

CLatticeData::~CLatticeData()
{
    if (NULL != m_pUpdator)
    {
        appSafeDelete(m_pUpdator);
    }

    for (INT i = 0; i < m_pActionList.Num(); ++i)
    {
        appSafeDelete(m_pActionList[i]);
    }

    if (NULL != m_pDeviceGaugeField)
    {
        checkCudaErrors(cudaFree(m_pDeviceGaugeField));
        m_pDeviceGaugeField = NULL;
    }
    if (NULL != m_pDeviceIndex)
    {
        _cDeletePtrs(m_pDeviceIndex);
        m_pDeviceIndex = NULL;
    }
    if (NULL != m_pRandom)
    {
        checkCudaErrors(cudaFree(m_pRandom));
        m_pRandom = NULL;
    }
    if (NULL != m_pDeviceRandomSchrage)
    {
        checkCudaErrors(cudaFree(m_pDeviceRandomSchrage));
        m_pDeviceRandomSchrage = NULL;
    }

    appSafeDelete(m_pGaugeField);
    appSafeDelete(m_pRandom);
    appSafeDelete(m_pRandomSchrage);
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================

