//=============================================================================
// FILENAME : CIndexSquare.cu
// 
// DESCRIPTION:
// This is the class for index on square lattice
//
// REVISION:
//  [12/5/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CIndexSquare)

#pragma region device Functions

static __device__ __inline__ SSmallInt4 _deviceCoordMoving(const SSmallInt4& sFrom, BYTE i)
{
    SBYTE offset = i < _DC_Dir ? -1 : 1;
    SSmallInt4 ret = sFrom;
    ret.m_byData4[i < _DC_Dir ? (4 - _DC_Dir + i) : (4 - _DC_Dir * 2 + i)] += offset;
    return ret;
}

static __device__ __inline__ UINT _deviceGetBigIndex(const SSmallInt4& sSite, const UINT* __restrict__ pSmallData)
{
    return (sSite.x + CIndexData::kCacheIndexEdge) * pSmallData[CIndexData::kMultX]
         + (sSite.y + CIndexData::kCacheIndexEdge) * pSmallData[CIndexData::kMultY]
         + (sSite.z + CIndexData::kCacheIndexEdge) * pSmallData[CIndexData::kMultZ]
         + (sSite.w + CIndexData::kCacheIndexEdge);
}

#pragma endregion

#pragma region Kernels

__global__ void _CLG_LAUNCH_BOUND_SINGLE
_kernalBakeSmallData(UINT* pDeviceData)
{
    pDeviceData[CIndexData::kMultX] = (_DC_Ly + 2 * CIndexData::kCacheIndexEdge)
        * (_DC_Lz + 2 * CIndexData::kCacheIndexEdge) 
        * (_DC_Lt + 2 * CIndexData::kCacheIndexEdge);
    pDeviceData[CIndexData::kMultY] = (_DC_Lz + 2 * CIndexData::kCacheIndexEdge)
        * (_DC_Lt + 2 * CIndexData::kCacheIndexEdge);
    pDeviceData[CIndexData::kMultZ] = _DC_Lt + 2 * CIndexData::kCacheIndexEdge;
    pDeviceData[CIndexData::kPlaqLengthIdx] = 4;
    pDeviceData[CIndexData::kPlaqPerSiteIdx] = _DC_Dim * (_DC_Dim - 1) / 2;
    pDeviceData[CIndexData::kPlaqPerLinkIdx] = 2 * (_DC_Dim - 1);
}

__global__ void _CLG_LAUNCH_BOUND
_kernalBakeWalkingTable(UINT* pDeviceData, uint3 mods)
{
    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
    SSmallInt4 coord;
    coord.x = static_cast<SBYTE>(idxAll / mods.x);
    coord.y = static_cast<SBYTE>((idxAll % mods.x) / mods.y);
    coord.z = static_cast<SBYTE>((idxAll % mods.y) / mods.z);
    coord.w = static_cast<SBYTE>(idxAll % mods.z);

    for (BYTE i = 0; i < 2 * _DC_Dir; ++i)
    {
        SSmallInt4 movingCoord = _deviceCoordMoving(coord, i);
        if (movingCoord.x < 0) movingCoord.x = 0;
        if (movingCoord.y < 0) movingCoord.y = 0;
        if (movingCoord.z < 0) movingCoord.z = 0;
        if (movingCoord.w < 0) movingCoord.w = 0;
        if (movingCoord.x >= _DC_Lx + 2 * CIndexData::kCacheIndexEdge) movingCoord.x = _DC_Lx + 2 * CIndexData::kCacheIndexEdge - 1;
        if (movingCoord.y >= _DC_Ly + 2 * CIndexData::kCacheIndexEdge) movingCoord.y = _DC_Ly + 2 * CIndexData::kCacheIndexEdge - 1;
        if (movingCoord.z >= _DC_Lz + 2 * CIndexData::kCacheIndexEdge) movingCoord.z = _DC_Lz + 2 * CIndexData::kCacheIndexEdge - 1;
        if (movingCoord.w >= _DC_Lt + 2 * CIndexData::kCacheIndexEdge) movingCoord.w = _DC_Lt + 2 * CIndexData::kCacheIndexEdge - 1;

        pDeviceData[idxAll * 2 * _DC_Dir + i] = movingCoord.x * mods.x + movingCoord.y * mods.y + movingCoord.z * mods.z + movingCoord.w;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernalBakeMappingTable(SSmallInt4* pDeviceData, uint3 mods)
{
    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
    SSmallInt4 coord;
    coord.x = static_cast<SBYTE>(idxAll / mods.x) - CIndexData::kCacheIndexEdge;
    coord.y = static_cast<SBYTE>((idxAll % mods.x) / mods.y) - CIndexData::kCacheIndexEdge;
    coord.z = static_cast<SBYTE>((idxAll % mods.y) / mods.z) - CIndexData::kCacheIndexEdge;
    coord.w = static_cast<SBYTE>(idxAll % mods.z) - CIndexData::kCacheIndexEdge;

    pDeviceData[idxAll] = coord;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelBakePlaqIndexAtSite(SIndex* pResult, 
    const UINT* __restrict__ pWalkingTable, 
    const SIndex* __restrict__ pMappingTable, 
    const UINT* __restrict__ pSmallDataTable)
{
    intokernalInt4;

    UINT uiDim = _DC_Dim;
    UINT uiBigSiteIndex = _deviceGetBigIndex(sSite4, pSmallDataTable);

    //24
    UINT iListIndex = uiSiteIndex 
        * (pSmallDataTable[CIndexData::kPlaqPerSiteIdx] 
         * pSmallDataTable[CIndexData::kPlaqLengthIdx]);

    for (BYTE uiLink = 0; uiLink < uiDim; ++uiLink)
    {
        for (BYTE uiPlaq = uiLink + 1; uiPlaq < uiDim; ++uiPlaq)
        {
            pResult[iListIndex] = SIndex(uiSiteIndex, uiLink);
            ++iListIndex;

            //start ----> uiLink
            UINT movedSite = pWalkingTable[uiBigSiteIndex * 2 * uiDim + (uiDim + uiLink)];
            pResult[iListIndex] = pMappingTable[movedSite];
            pResult[iListIndex].m_byDir = uiPlaq;
            ++iListIndex;

            //start ----> uiPlaq
            movedSite = pWalkingTable[uiBigSiteIndex * 2 * uiDim + (uiDim + uiPlaq)];
            pResult[iListIndex] = pMappingTable[movedSite];
            pResult[iListIndex].m_byDir = uiLink;
            pResult[iListIndex].m_byTag |= _kDaggerOrOpposite;
            ++iListIndex;

            pResult[iListIndex] = SIndex(uiSiteIndex, uiPlaq);
            pResult[iListIndex].m_byTag |= _kDaggerOrOpposite;
            ++iListIndex;
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelBakePlaqIndexAtLink(SIndex* pResult, 
    const UINT* __restrict__ pWalkingTable,
    const SIndex* __restrict__ pMappingTable,
    const UINT* __restrict__ pSmallDataTable)
{
    intokernalInt4;

    UINT uiDim = _DC_Dim;
    UINT uiBigSiteIndex = _deviceGetBigIndex(sSite4, pSmallDataTable);
    
    for (UINT uiLinkDir = 0; uiLinkDir < uiDim; ++uiLinkDir)
    {
        UINT uiLinkIndex = uiSiteIndex * uiDim + uiLinkDir;

        UINT iListIndex = uiLinkIndex
            * (pSmallDataTable[CIndexData::kPlaqPerLinkIdx]
            * (pSmallDataTable[CIndexData::kPlaqLengthIdx] - 1));

        for (BYTE i = 0; i < uiDim; ++i)
        {
            if (i != uiLinkDir)
            {
                //=============================================
                //add forward
                //[site][p_dir], [site+p_dir][b_dir], [site+b_dir][p_dir]^1
                pResult[iListIndex] = SIndex(uiSiteIndex, i);
                ++iListIndex;

                //start ---> i
                UINT movedSite = pWalkingTable[uiBigSiteIndex * 2 * uiDim + (uiDim + i)];
                pResult[iListIndex] = pMappingTable[movedSite];
                pResult[iListIndex].m_byDir = uiLinkDir;
                ++iListIndex;

                //start ---> uiLinkDir
                movedSite = pWalkingTable[uiBigSiteIndex * 2 * uiDim + (uiDim + uiLinkDir)];
                pResult[iListIndex] = pMappingTable[movedSite];
                pResult[iListIndex].m_byDir = i;
                pResult[iListIndex].m_byTag |= _kDaggerOrOpposite;
                ++iListIndex;

                //=============================================
                //add backward
                //[site-p_dir][p_dir]^-1, [site-p_dir][b_dir], [site-p_dir+b_dir][p_dir]
                //i <---- start
                movedSite = pWalkingTable[uiBigSiteIndex * 2 * uiDim + i];
                pResult[iListIndex] = pMappingTable[movedSite];
                pResult[iListIndex].m_byDir = i;
                pResult[iListIndex].m_byTag |= _kDaggerOrOpposite;
                ++iListIndex;

                //
                pResult[iListIndex] = SIndex(pResult[iListIndex - 1].m_uiSiteIndex, uiLinkDir);
                if (pResult[iListIndex - 1].IsDirichlet())
                {
                    pResult[iListIndex].m_byTag |= _kDirichlet;
                }
                ++iListIndex;

                //last ----> uiLinkDir
                movedSite = pWalkingTable[movedSite * 2 * uiDim + (uiDim + uiLinkDir)];
                pResult[iListIndex] = pMappingTable[movedSite];
                pResult[iListIndex].m_byDir = i;
                ++iListIndex;
            }
        }
    }
}

/**
* gaugemove[linkIndex] = gauge[uiSite - linkIndex]_{linkIndex}
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCacheGaugeMove(SIndex* pCached, 
    const UINT* __restrict__ pWalkingTable,
    const SIndex* __restrict__ pMappingTable,
    const UINT* __restrict__ pSmallDataTable)
{
    intokernalInt4;

    UINT uiDir = _DC_Dir;
    UINT uiBigSiteIndex = _deviceGetBigIndex(sSite4, pSmallDataTable);
    for (UINT i = 0; i < uiDir; ++i)
    {
        pCached[uiSiteIndex * uiDir + i] = pMappingTable[pWalkingTable[uiBigSiteIndex * 2 * uiDir + i]];
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCacheFermionMove(SIndex* pCached, 
    const UINT* __restrict__ pWalkingTable,
    const SIndex* __restrict__ pMappingTable,
    const UINT* __restrict__ pSmallDataTable)
{
    intokernalInt4;

    UINT uiDir = _DC_Dir;
    UINT uiBigSiteIndex = _deviceGetBigIndex(sSite4, pSmallDataTable);
    for (UINT i = 0; i < uiDir; ++i)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, i);
        //first element is right, second element is left.
        pCached[linkIndex * 2]
            = pMappingTable[pWalkingTable[uiBigSiteIndex * 2 * uiDir + i + uiDir]];
        pCached[linkIndex * 2 + 1]
            = pMappingTable[pWalkingTable[uiBigSiteIndex * 2 * uiDir + i]];
    }
}

#pragma endregion

UINT CIndexSquare::GetDecompose(UINT volumn)
{
    TArray<UINT> factors = _getFactors(volumn);
    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock();
    UINT maxThreadPerBlock = deviceConstraints[1]; //we only use 1 dimension, so it is the constraint of blockDim.x

    UINT uiMax = 1;
    for (INT i = 0; i < factors.Num(); ++i)
    {
        if (factors[i] <= maxThreadPerBlock && factors[i] > uiMax)
        {
            uiMax = factors[i];
        }
    }
    
    return uiMax;
}

void CIndexSquare::BakeAllIndexBuffer(CIndexData* pData)
{
    uint4 biggerLattice;
    biggerLattice.x = _HC_Lx + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.y = _HC_Ly + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.z = _HC_Lz + 2 * CIndexData::kCacheIndexEdge;
    biggerLattice.w = _HC_Lt + 2 * CIndexData::kCacheIndexEdge;
    uint3 biggerLatticeMod;

    UINT uiVolumn = biggerLattice.x * biggerLattice.y * biggerLattice.z * biggerLattice.w;
    UINT threadPerSite = GetDecompose(uiVolumn);
    dim3 threads(threadPerSite, 1, 1);
    dim3 blocks(uiVolumn / threadPerSite, 1, 1);
    biggerLatticeMod.x = biggerLattice.y * biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.y = biggerLattice.z * biggerLattice.w;
    biggerLatticeMod.z = biggerLattice.w;

    pData->m_uiPlaqutteLength = 4;
    pData->m_uiPlaqutteCountPerSite = static_cast<BYTE>(_HC_Dim * (_HC_Dim - 1) / 2);
    pData->m_uiPlaqutteCountPerLink = static_cast<BYTE>(2 * (_HC_Dim - 1));

    //bake small data
    _kernalBakeSmallData << <1, 1 >> > (pData->m_pSmallData);
    //bake walking index
    _kernalBakeWalkingTable << <blocks, threads >> > (pData->m_pWalkingTable, biggerLatticeMod);

    //bake index mappings
    for (BYTE i = 1; i < kMaxFieldCount; ++i)
    {
        if (NULL != appGetLattice()->GetFieldById(i))
        {
            checkCudaErrors(cudaMalloc((void**)&pData->m_pIndexPositionToSIndex[i], sizeof(SIndex)
                * (_HC_Lx + 2 * CIndexData::kCacheIndexEdge) * (_HC_Ly + 2 * CIndexData::kCacheIndexEdge)
                * (_HC_Lz + 2 * CIndexData::kCacheIndexEdge) * (_HC_Lt + 2 * CIndexData::kCacheIndexEdge)
                ));

            //bake map index
            _kernalBakeMappingTable << <blocks, threads >> > (pData->m_pMappingTable, biggerLatticeMod);

            //bake boundary condition
            m_pBoundaryCondition->BakeEdgePoints(i, pData->m_pMappingTable, pData->m_pIndexPositionToSIndex[i]);
        }
    }

    //bake region id table
    m_pBoundaryCondition->BakeRegionTable(pData->m_byRegionTable);

    //copy the mapping table to device
    checkCudaErrors(cudaMalloc((void**)&pData->m_pDeviceIndexPositionToSIndex, sizeof(SIndex*) * kMaxFieldCount));
    checkCudaErrors(cudaMemcpy(pData->m_pDeviceIndexPositionToSIndex, pData->m_pIndexPositionToSIndex, sizeof(SIndex*) * kMaxFieldCount, cudaMemcpyHostToDevice));
}

void CIndexSquare::BakePlaquttes(CIndexData* pData, BYTE byFieldId)
{
    //If Has Gauge field
    if (NULL != appGetLattice()->GetFieldById(byFieldId))
    {
        preparethread;

        checkCudaErrors(cudaMalloc((void**)&pData->m_pPlaqutteCache, sizeof(SIndex) * _HC_Volume * (_HC_Dim * (_HC_Dim - 1) / 2) * 4));
        checkCudaErrors(cudaMalloc((void**)&pData->m_pStappleCache, sizeof(SIndex) * _HC_Volume * _HC_Dim * (2 * (_HC_Dim - 1)) * 3));

        //bake plaqutte per site
        _kernelBakePlaqIndexAtSite << <block, threads >> > (pData->m_pPlaqutteCache, pData->m_pWalkingTable, pData->m_pIndexPositionToSIndex[byFieldId], pData->m_pSmallData);

        //bake plaqutte per link
        _kernelBakePlaqIndexAtLink << <block, threads >> > (pData->m_pStappleCache, pData->m_pWalkingTable, pData->m_pIndexPositionToSIndex[byFieldId], pData->m_pSmallData);
    }
}

void CIndexSquare::BakeMoveIndex(CIndexData* pData, BYTE byFieldId)
{
    preparethread;

    checkCudaErrors(cudaMalloc((void**)&pData->m_pGaugeMoveCache[byFieldId], sizeof(SIndex) * _HC_Volume * _HC_Dir));
    checkCudaErrors(cudaMalloc((void**)&pData->m_pFermionMoveCache[byFieldId], sizeof(SIndex) * _HC_Volume * _HC_Dir * 2));

    _kernelCacheGaugeMove << <block, threads >> > (pData->m_pGaugeMoveCache[byFieldId], pData->m_pWalkingTable, pData->m_pIndexPositionToSIndex[1], pData->m_pSmallData);
    _kernelCacheFermionMove << <block, threads >> > (pData->m_pFermionMoveCache[byFieldId], pData->m_pWalkingTable, pData->m_pIndexPositionToSIndex[byFieldId], pData->m_pSmallData);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================