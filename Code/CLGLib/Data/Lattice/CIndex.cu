//=============================================================================
// FILENAME : CIndex.cu
// 
// DESCRIPTION:
// This is the class for boundary conditions
// Note that, the boundary conditions should only make sense together with lattice!!
//
// REVISION:
//  [12/16/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernels

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernelCreateIndex(void** devicePtr, deviceBoundaryCondition ** pBC, UINT* size, EIndexType eIT)
{
    if (eIT == EIndexType_Square)
    {
        (*devicePtr) = (void*)(new CIndexSquare(*pBC));
        size[0] =(UINT)sizeof(CIndexSquare);
        return;
    }
    size[0] = 0;
}

/**
* Initial SU3 Field with a value
*/
__global__ void
_CLG_LAUNCH_BOUND
_kernelCachePlaqIndex(SIndex *pDevicePlaqPerSite, SIndex *pDevicePlaqPerLink, BYTE byPlaqPerLink, BYTE byPlaqPerSite, BYTE byPlaqLength)
{
    intokernaldir;

    __idx->_deviceGetPlaquttesAtSiteAll(pDevicePlaqPerSite + uiSiteIndex * (byPlaqPerSite * byPlaqLength), uiSiteIndex);

    for (BYTE idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        __idx->_deviceGetPlaquttesAtLinkAll(pDevicePlaqPerLink + linkIndex * (byPlaqPerLink * (byPlaqLength - 1)), linkIndex);

        //for (BYTE iplaqIndex = 0; iplaqIndex < uiPlaqutteCount * (uiPlaqutteLength - 1); ++iplaqIndex)
        //{
        //    pDevicePlaqPerLink[linkIndex * (uiPlaqutteCount * (uiPlaqutteLength - 1)) + iplaqIndex] = tmpPlaq[iplaqIndex];
        //}
    }
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelCacheMoveIndex(BYTE byFieldId, SIndex *pGaugeMove, SIndex *pFermionMove)
{
    intokernaldir;

    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        pGaugeMove[linkIndex] = __idx->_deviceGaugeIndexWalk(uiSiteIndex, -(idir + 1));
        pFermionMove[linkIndex * 2] = __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, (idir + 1));
        pFermionMove[linkIndex * 2 + 1] = __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, -(idir + 1));
    }
}

#pragma endregion

extern "C" 
{ 
    void _cCreateIndex(void** devicePtr, deviceBoundaryCondition ** pBC, UINT* size, EIndexType eIT)
    {
        _kernelCreateIndex << <1, 1 >> > (devicePtr, pBC, size, eIT);
    }
}

void CIndexCache::CachePlaquttes()
{
    appGetLattice()->GetPlaquetteLengthCount(m_uiPlaqutteLength, m_uiPlaqutteCountPerSite, m_uiPlaqutteCountPerLink);

    checkCudaErrors(cudaMalloc((void**)&m_pPlaqutteCache, _HC_Volumn * sizeof(SIndex) * m_uiPlaqutteLength * m_uiPlaqutteCountPerSite));
    checkCudaErrors(cudaMalloc((void**)&m_pStappleCache, _HC_Volumn * _HC_Dir * sizeof(SIndex) * (m_uiPlaqutteLength - 1) * m_uiPlaqutteCountPerLink));

    preparethread;
    _kernelCachePlaqIndex << <block, threads >> > (m_pPlaqutteCache, m_pStappleCache, m_uiPlaqutteCountPerLink, m_uiPlaqutteCountPerSite, m_uiPlaqutteLength);

    appParanoiac(_T(" Get Plaq Length Count %d, %d, %d \n"), m_uiPlaqutteLength, m_uiPlaqutteCountPerSite, m_uiPlaqutteCountPerLink);
}


void CIndexCache::CacheFermion(BYTE byFieldId)
{
    assert(byFieldId < kMaxFieldCount);
    checkCudaErrors(cudaMalloc((void**)&m_pGaugeMoveCache[byFieldId], _HC_Volumn * _HC_Dir * sizeof(SIndex)));
    checkCudaErrors(cudaMalloc((void**)&m_pFermionMoveCache[byFieldId], _HC_Volumn * _HC_Dir * 2 * sizeof(SIndex)));

    preparethread;
    _kernelCacheMoveIndex << <block, threads >> > (byFieldId, m_pGaugeMoveCache[byFieldId], m_pFermionMoveCache[byFieldId]);
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================