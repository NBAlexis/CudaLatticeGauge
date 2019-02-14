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

__global__ void _kernelGetPlaqLengthCount(UINT* deviceData)
{
    UINT length, countPersite, countPerLink;
    __idx->_deviceGetPlaqutteCountLength(length, countPersite, countPerLink);

    deviceData[0] = length;
    deviceData[1] = countPersite;
    deviceData[2] = countPerLink;
}

/**
* m_uiLatticeDecompose[0,1,2] is the blocks
* m_uiLatticeDecompose[3,4,5] is the threads in blocks
*/
CLatticeData::CLatticeData()
    : m_pRandom(NULL)
    , m_pGaugeField(NULL)
    , m_pUpdator(NULL)

    , m_pDeviceRandom(NULL)
    , m_pDeviceIndex(NULL)

    , m_pFermionSolver(NULL)
    , m_pMeasurements(NULL)
    , m_pFieldCache(NULL)
{
    m_pFieldCache = new CFieldCache();
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

    for (INT i = 0; i < m_pOtherFields.Num(); ++i)
    {
        appSafeDelete(m_pOtherFields[i]);
    }

    for (INT i = 0; i < m_pFieldPools.Num(); ++i)
    {
        appSafeDelete(m_pFieldPools[i]);
    }

    appSafeDelete(m_pFieldCache);
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
        exit(EXIT_FAILURE);
    }
    m_pFermionSolver->Configurate(param);
    m_pFermionSolver->AllocateBuffers(pFermionField);

    appGeneral(_T("Create sparse linear algebra solver: %s \n"), sSolver.c_str());
}

void CLatticeData::CreateFieldPool(BYTE byFieldId, UINT uiCount)
{
    if (m_pFieldPoolMap.Exist(byFieldId))
    {
        appGeneral(_T("Create field pool, but field id already exist!\n"));
        return;
    }
    CField* pField = GetFieldById(byFieldId);
    if (NULL == pField)
    {
        appCrucial(_T("Create field pool, but field cannot be found!\n"));
        return;
    }

    CFieldPool * pFieldPool = new CFieldPool(pField, uiCount);
    m_pFieldPools.AddItem(pFieldPool);
    m_pFieldPoolMap.SetAt(byFieldId, pFieldPool);
}

CField* CLatticeData::GetPooledFieldById(BYTE byId)
{
    if (!m_pFieldPoolMap.Exist(byId))
    {
        appCrucial(_T("Get Pooled field failed!\n"));
        return NULL;
    }

    return m_pFieldPoolMap[byId]->GetOne();
}

void CLatticeData::OnUpdatorConfigurationAccepted()
{
    if (NULL != m_pMeasurements)
    {
        m_pMeasurements->OnConfigurationAccepted();
    }
}

void CLatticeData::OnUpdatorFinished(UBOOL bMeasured)
{
    if (NULL != m_pMeasurements && bMeasured)
    {
        m_pMeasurements->OnUpdateFinished(FALSE);
    }
}

void CLatticeData::GetPlaquetteLengthCount(UINT& plaqLength, UINT& countPerSite, UINT& countPerLink)
{
    UINT * deviceData;
    checkCudaErrors(cudaMalloc((void**)&deviceData, sizeof(UINT) * 3));

    _kernelGetPlaqLengthCount << <1, 1 >> > (deviceData);

    UINT hostData[3];

    checkCudaErrors(cudaMemcpy(hostData, deviceData, sizeof(UINT) * 3, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(deviceData));

    plaqLength = hostData[0];
    countPerSite = hostData[1];
    countPerLink = hostData[2];
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
