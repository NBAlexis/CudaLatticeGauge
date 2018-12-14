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


//__global__ 
//void _kernelPrepareDeviceData(CDeviceLattice *& deviceLattice, CLatticeData* pLattice)
//{
    /*deviceLattice = new CDeviceLattice();
    
    memcpy(deviceLattice->m_uiLatticeLength, pLattice->m_uiLatticeLength, sizeof(UINT) * CCommonData::kMaxDim);
    memcpy(deviceLattice->m_uiLatticeDecompose, pLattice->m_uiLatticeDecompose, sizeof(UINT) * (CCommonData::kLatticeDecompose * 2));
    memcpy(deviceLattice->m_uiLatticeMultipy, pLattice->m_uiLatticeMultipy, sizeof(UINT) * (CCommonData::kMaxDim - 1));

    deviceLattice->m_uiVolumn = pLattice->m_uiVolumn;
    deviceLattice->m_uiDim = pLattice->m_uiDim;
    deviceLattice->m_uiDir = pLattice->m_uiDir;
    deviceLattice->m_uiTLength = pLattice->m_uiTLength;
    deviceLattice->m_fBeta = pLattice->m_fBeta;
    
    //create index
    deviceLattice->m_pIndex = new CIndexSquare(deviceLattice);*/
//}

//__global__ 
//void _kernelReleaseDeviceData(CDeviceLattice *& deviceLattice, CRandomSchrage *& random)
//{
    //delete random;
    //random = NULL;

    //delete deviceLattice;
    //deviceLattice = NULL;
//}

/**
* m_uiLatticeDecompose[0,1,2] is the blocks
* m_uiLatticeDecompose[3,4,5] is the threads in blocks
*/
CLatticeData::CLatticeData()
    : m_pRandom(NULL)
    , m_pRandomSchrage(NULL)
    , m_pGaugeField(NULL)
    , m_pGaugeFieldStaple(NULL)
    , m_pIndex(NULL)

    , m_pDeviceRandom(NULL)
    , m_pDeviceRandomSchrage(NULL)
    , m_pDeviceGaugeField(NULL)
    , m_pDeviceGaugeFieldStaple(NULL)
    , m_pDeviceIndex(NULL)
{
    //m_uiDim = CCommonData::m_uiDim;
    //m_uiDir = CCommonData::m_uiDir;
    //m_fBeta = CCommonData::m_fBeta;
    //memcpy(m_uiLatticeLength, CCommonData::m_uiLatticeLength, sizeof(UINT) * CCommonData::kMaxDim);
    //m_uiVolumn = m_uiLatticeLength[0] * m_uiLatticeLength[1] * m_uiLatticeLength[2] * m_uiLatticeLength[3];

    //m_uiLatticeMultipy[0] = m_uiLatticeLength[1] * m_uiLatticeLength[2] * m_uiLatticeLength[3];
    //m_uiLatticeMultipy[1] = m_uiLatticeLength[2] * m_uiLatticeLength[3];
    //m_uiLatticeMultipy[2] = m_uiLatticeLength[3];

    //if (m_uiLatticeLength[0] * m_uiLatticeLength[1] * m_uiLatticeLength[2] <= CCommonData::m_uiMaxThread)
    //{
    //    m_uiLatticeDecompose[0] = 1;
    //    m_uiLatticeDecompose[1] = 1;
    //    m_uiLatticeDecompose[2] = 1;

    //    m_uiLatticeDecompose[3] = m_uiLatticeLength[0];
    //    m_uiLatticeDecompose[4] = m_uiLatticeLength[1];
    //    m_uiLatticeDecompose[5] = m_uiLatticeLength[2];
    //}
    //else if (m_uiLatticeLength[1] * m_uiLatticeLength[2] <= CCommonData::m_uiMaxThread)
    //{
    //    m_uiLatticeDecompose[0] = m_uiLatticeLength[0];
    //    m_uiLatticeDecompose[1] = 1;
    //    m_uiLatticeDecompose[2] = 1;

    //    m_uiLatticeDecompose[3] = 1;
    //    m_uiLatticeDecompose[4] = m_uiLatticeLength[1];
    //    m_uiLatticeDecompose[5] = m_uiLatticeLength[2];
    //}
    //else if (m_uiLatticeLength[2] <= CCommonData::m_uiMaxThread)
    //{
    //    m_uiLatticeDecompose[0] = m_uiLatticeLength[0];
    //    m_uiLatticeDecompose[1] = m_uiLatticeLength[1];
    //    m_uiLatticeDecompose[2] = 1;

    //    m_uiLatticeDecompose[3] = 1;
    //    m_uiLatticeDecompose[4] = 1;
    //    m_uiLatticeDecompose[5] = m_uiLatticeLength[2];
    //}
    //else
    //{
    //    appCrucial("CLatticeData:: Fail to divide the blocks!");
    //    m_uiLatticeDecompose[0] = m_uiLatticeLength[0];
    //    m_uiLatticeDecompose[1] = m_uiLatticeLength[1];
    //    m_uiLatticeDecompose[2] = m_uiLatticeLength[2];

    //    m_uiLatticeDecompose[3] = 1;
    //    m_uiLatticeDecompose[4] = 1;
    //    m_uiLatticeDecompose[5] = 1;
    //}

    #pragma region Create Index and Boundary

    //for now we only support square
    //_kernelPrepareDeviceData << <1, 1 >> > (m_pDeviceInstance, this);

    //checkCudaErrors(cudaGetLastError());
    //checkCudaErrors(cudaDeviceSynchronize());

    #pragma endregion Create Index and Boundary

    //m_pDeviceInstance->m_pDeviceRandom = new CRandomSchrage(CCommonData::m_uiSeed, m_pDeviceInstance);

#pragma region Create Fields

    //Always create a gauge field
    //m_pDeviceInstance->m_pGaugeField = new CFieldGaugeSU3(this);

#pragma endregion Create Fields
}

CLatticeData::~CLatticeData()
{
    if (NULL != m_pDeviceGaugeField)
    {
        checkCudaErrors(cudaFree(m_pDeviceGaugeField));
        m_pDeviceGaugeField = NULL;
    }
    if (NULL != m_pDeviceIndex)
    {
        checkCudaErrors(cudaFree(m_pDeviceIndex));
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
    appSafeDelete(m_pIndex);
    appSafeDelete(m_pRandom);
    appSafeDelete(m_pRandomSchrage);
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================

