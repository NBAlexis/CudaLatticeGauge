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
__global__
void _kernelCachePlaqIndex(SIndex *pDevicePlaqPerSite, SIndex *pDevicePlaqPerLink)
{
    intokernal;

    SIndex plaquttes[kMaxPlaqutteCache];

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        UINT siteIndex = _deviceGetSiteIndex(coord);
        UINT uiPlaqutteCount = 0;
        UINT uiPlaqutteLength = 0;
        __idx->_deviceGetPlaquttesAtSite(plaquttes, uiPlaqutteCount, uiPlaqutteLength, siteIndex);

        for (UINT iplaqIndex = 0; iplaqIndex < uiPlaqutteCount * uiPlaqutteLength; ++iplaqIndex)
        {
            pDevicePlaqPerSite[siteIndex * (uiPlaqutteCount * uiPlaqutteLength) + iplaqIndex] = plaquttes[iplaqIndex];
        }

        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);
            __idx->_deviceGetPlaquttesAtLink(plaquttes, uiPlaqutteCount, uiPlaqutteLength, linkIndex);

            for (UINT iplaqIndex = 0; iplaqIndex < uiPlaqutteCount * (uiPlaqutteLength - 1); ++iplaqIndex)
            {
                pDevicePlaqPerLink[linkIndex * (uiPlaqutteCount * (uiPlaqutteLength - 1)) + iplaqIndex] = plaquttes[iplaqIndex];
            }
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
{
    m_uiLinkeCount = _HC_Volumn * _HC_Dir;
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

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================