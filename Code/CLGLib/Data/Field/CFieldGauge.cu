//=============================================================================
// FILENAME : CFieldGauge.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [02/03/2018 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region Kernels

/**
* Initial SU3 Field with a value
*/
__global__ void
_CLG_LAUNCH_BOUND
_kernelCachePlaqIndex(SIndex *pDevicePlaqPerSite, SIndex *pDevicePlaqPerLink)
{
    intokernaldir;

    SIndex tmpPlaq[kMaxPlaqutteCache];
    UINT uiPlaqutteCount = 0;
    UINT uiPlaqutteLength = 0;
    __idx->_deviceGetPlaquttesAtSite(tmpPlaq, uiPlaqutteCount, uiPlaqutteLength, uiSiteIndex);

    for (UINT iplaqIndex = 0; iplaqIndex < uiPlaqutteCount * uiPlaqutteLength; ++iplaqIndex)
    {
        pDevicePlaqPerSite[uiSiteIndex * (uiPlaqutteCount * uiPlaqutteLength) + iplaqIndex] = tmpPlaq[iplaqIndex];
    }

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        __idx->_deviceGetPlaquttesAtLink(tmpPlaq, uiPlaqutteCount, uiPlaqutteLength, linkIndex);

        for (UINT iplaqIndex = 0; iplaqIndex < uiPlaqutteCount * (uiPlaqutteLength - 1); ++iplaqIndex)
        {
            pDevicePlaqPerLink[linkIndex * (uiPlaqutteCount * (uiPlaqutteLength - 1)) + iplaqIndex] = tmpPlaq[iplaqIndex];
        }
    }
}

#pragma endregion

CFieldGauge::CFieldGauge()
    : CField()
    , m_pPlaquttesPerSite(NULL)
    , m_pPlaquttesPerLink(NULL)
    , m_bPlaqutteIndexCached(FALSE)
    , m_uiPlaqutteLength(0)
    , m_uiPlaqutteCountPerSite(0)
    , m_uiPlaqutteCountPerLink(0)
    , m_uiLinkeCount(_HC_Volumn * _HC_Dir)
{

}

CFieldGauge::~CFieldGauge()
{
    if (NULL != m_pPlaquttesPerSite)
    {
        checkCudaErrors(cudaFree(m_pPlaquttesPerSite));
    }

    if (NULL != m_pPlaquttesPerLink)
    {
        checkCudaErrors(cudaFree(m_pPlaquttesPerLink));
    }
}

void CFieldGauge::CachePlaqutteIndexes()
{
    appGetLattice()->GetPlaquetteLengthCount(m_uiPlaqutteLength, m_uiPlaqutteCountPerSite, m_uiPlaqutteCountPerLink);

    checkCudaErrors(cudaMalloc((void**)&m_pPlaquttesPerSite, _HC_Volumn * sizeof(SIndex) * m_uiPlaqutteLength * m_uiPlaqutteCountPerSite));
    checkCudaErrors(cudaMalloc((void**)&m_pPlaquttesPerLink, _HC_Volumn * _HC_Dir * sizeof(SIndex) * (m_uiPlaqutteLength - 1) * m_uiPlaqutteCountPerLink));

    preparethread;
    _kernelCachePlaqIndex << <block, threads >> > (m_pPlaquttesPerSite, m_pPlaquttesPerLink);

    
    m_bPlaqutteIndexCached = TRUE;

    appParanoiac(_T(" Get Plaq Length Count %d, %d, %d \n"), m_uiPlaqutteLength, m_uiPlaqutteCountPerSite, m_uiPlaqutteCountPerLink);
}

void CFieldGauge::CopyTo(CField* U) const
{
    CFieldGauge* pFieldGauge = dynamic_cast<CFieldGauge*>(U);
    if (NULL == pFieldGauge)
    {
        return;
    }

    CField::CopyTo(U);

    pFieldGauge->m_uiLinkeCount = m_uiLinkeCount;

    if (m_bPlaqutteIndexCached)
    {
        pFieldGauge->m_bPlaqutteIndexCached = TRUE;
        pFieldGauge->m_uiPlaqutteLength = m_uiPlaqutteLength;
        pFieldGauge->m_uiPlaqutteCountPerSite = m_uiPlaqutteCountPerSite;
        pFieldGauge->m_uiPlaqutteCountPerLink = m_uiPlaqutteCountPerLink;
        if (NULL == pFieldGauge->m_pPlaquttesPerSite)
        {
            checkCudaErrors(cudaMalloc((void**)&pFieldGauge->m_pPlaquttesPerSite, _HC_Volumn * sizeof(SIndex) * m_uiPlaqutteLength * m_uiPlaqutteCountPerSite));
        }
        if (NULL == pFieldGauge->m_pPlaquttesPerLink)
        {
            checkCudaErrors(cudaMalloc((void**)&pFieldGauge->m_pPlaquttesPerLink, _HC_Volumn * _HC_Dir * sizeof(SIndex) * (m_uiPlaqutteLength - 1) * m_uiPlaqutteCountPerLink));
        }
        checkCudaErrors(cudaMemcpy(pFieldGauge->m_pPlaquttesPerSite, m_pPlaquttesPerSite, _HC_Volumn * sizeof(SIndex) * m_uiPlaqutteLength * m_uiPlaqutteCountPerSite, cudaMemcpyDeviceToDevice));
        checkCudaErrors(cudaMemcpy(pFieldGauge->m_pPlaquttesPerLink, m_pPlaquttesPerLink, _HC_Volumn * _HC_Dir * sizeof(SIndex) * (m_uiPlaqutteLength - 1) * m_uiPlaqutteCountPerSite, cudaMemcpyDeviceToDevice));
    }
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================