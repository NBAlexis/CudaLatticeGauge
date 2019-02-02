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

Real CLGAPI CCommonData::m_fBeta = 0;
Real CLGAPI CCommonData::m_fKai = 0;

__global__ void _kernelDeletePtrs(CIndex * pdeviceIndex)
{
    if (NULL != pdeviceIndex)
    {
        delete pdeviceIndex;
    }
}

/**
* m_uiLatticeDecompose[0,1,2] is the blocks
* m_uiLatticeDecompose[3,4,5] is the threads in blocks
*/
CLatticeData::CLatticeData()
    : m_pRandom(NULL)
    , m_pGaugeField(NULL)
    , m_pGaugeFieldStaple(NULL)
    , m_pUpdator(NULL)

    , m_pDeviceRandom(NULL)
    , m_pDeviceGaugeField(NULL)
    , m_pDeviceGaugeFieldStaple(NULL)
    , m_pDeviceIndex(NULL)

    , m_pFermionSolver(NULL)
    , m_pMeasurements(NULL)
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
        _kernelDeletePtrs << <1, 1 >> > (m_pDeviceIndex);
        m_pDeviceIndex = NULL;
    }
    if (NULL != m_pDeviceRandom)
    {
        checkCudaErrors(cudaFree(m_pDeviceRandom));
        m_pDeviceRandom = NULL;
    }

    appSafeDelete(m_pGaugeField);
    appSafeDelete(m_pRandom);
    appSafeDelete(m_pMeasurements);
    appSafeDelete(m_pFermionSolver);
}

void CLatticeData::CreateFermionSolver(const CCString& sSolver, const CParameters& param, const CField* pFermionField)
{
    CBase* pSolver = appCreate(sSolver);
    m_pFermionSolver = dynamic_cast<CSLASolver*>(pSolver);
    if (NULL == m_pFermionSolver)
    {
        appCrucial(_T("Create Fermion Solver %s failed!\n"), sSolver.c_str());
    }
    m_pFermionSolver->Configurate(param);
    m_pFermionSolver->AllocateBuffers(pFermionField);
}

void CLatticeData::OnUpdatorConfigurationAccepted()
{
    if (NULL != m_pMeasurements)
    {
        m_pMeasurements->OnConfigurationAccepted();
    }
}

void CLatticeData::OnUpdatorFinished()
{
    if (NULL != m_pMeasurements)
    {
        m_pMeasurements->OnUpdateFinished(FALSE);
    }
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================

