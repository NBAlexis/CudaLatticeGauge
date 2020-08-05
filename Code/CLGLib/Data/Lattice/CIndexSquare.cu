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
    const UINT* __restrict__ pSmallDataTable,
    const BYTE* __restrict__ pBondInfoTable)
{
    intokernalInt4;

    const UINT uiDim = _DC_Dim;
    const UINT uiBigSiteIndex = _deviceGetBigIndex(sSite4, pSmallDataTable);

    //24
    UINT iListIndex = uiSiteIndex 
        * (pSmallDataTable[CIndexData::kPlaqPerSiteIdx] 
         * pSmallDataTable[CIndexData::kPlaqLengthIdx]);

    for (BYTE uiLink = 0; uiLink < uiDim; ++uiLink)
    {
        for (BYTE uiPlaq = uiLink + 1; uiPlaq < uiDim; ++uiPlaq)
        {
            pResult[iListIndex] = SIndex(uiSiteIndex, uiLink);
            pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, uiBigSiteIndex, uiLink) ? _kDirichlet : 0;
            ++iListIndex;

            //start ----> uiLink
            UINT movedSite = pWalkingTable[uiBigSiteIndex * 2 * uiDim + (uiDim + uiLink)];
            pResult[iListIndex] = pMappingTable[movedSite];
            pResult[iListIndex].m_byDir = uiPlaq;
            pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, movedSite, uiPlaq) ? _kDirichlet : 0;
            ++iListIndex;

            //start ----> uiPlaq
            movedSite = pWalkingTable[uiBigSiteIndex * 2 * uiDim + (uiDim + uiPlaq)];
            pResult[iListIndex] = pMappingTable[movedSite];
            pResult[iListIndex].m_byDir = uiLink;
            pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, movedSite, uiLink) ? _kDirichlet : 0;
            pResult[iListIndex].m_byTag |= _kDaggerOrOpposite;
            ++iListIndex;

            pResult[iListIndex] = SIndex(uiSiteIndex, uiPlaq);
            pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, uiBigSiteIndex, uiPlaq) ? _kDirichlet : 0;
            pResult[iListIndex].m_byTag |= _kDaggerOrOpposite;
            ++iListIndex;
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelBakePlaqIndexAtLink(SIndex* pResult, 
    const UINT* __restrict__ pWalkingTable,
    const SIndex* __restrict__ pMappingTable,
    const UINT* __restrict__ pSmallDataTable,
    const BYTE* __restrict__ pBondInfoTable)
{
    intokernalInt4;

    UINT uiDim = _DC_Dim;
    const UINT uiBigSiteIndex = _deviceGetBigIndex(sSite4, pSmallDataTable);
    
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
                pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, uiBigSiteIndex, i) ? _kDirichlet : 0;
                ++iListIndex;

                //start ---> i
                UINT movedSite = pWalkingTable[uiBigSiteIndex * 2 * uiDim + (uiDim + i)];
                pResult[iListIndex] = pMappingTable[movedSite];
                pResult[iListIndex].m_byDir = uiLinkDir;
                pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, movedSite, uiLinkDir) ? _kDirichlet : 0;
                ++iListIndex;

                //start ---> uiLinkDir
                movedSite = pWalkingTable[uiBigSiteIndex * 2 * uiDim + (uiDim + uiLinkDir)];
                pResult[iListIndex] = pMappingTable[movedSite];
                pResult[iListIndex].m_byDir = i;
                pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, movedSite, i) ? _kDirichlet : 0;
                pResult[iListIndex].m_byTag |= _kDaggerOrOpposite;
                ++iListIndex;

                //=============================================
                //add backward
                //[site-p_dir][p_dir]^-1, [site-p_dir][b_dir], [site-p_dir+b_dir][p_dir]
                //i <---- start
                movedSite = pWalkingTable[uiBigSiteIndex * 2 * uiDim + i];
                pResult[iListIndex] = pMappingTable[movedSite];
                pResult[iListIndex].m_byDir = i;
                pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, movedSite, i) ? _kDirichlet : 0;
                pResult[iListIndex].m_byTag |= _kDaggerOrOpposite;
                ++iListIndex;

                //
                pResult[iListIndex] = SIndex(pResult[iListIndex - 1].m_uiSiteIndex, uiLinkDir);
                pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, movedSite, uiLinkDir) ? _kDirichlet : 0;
                ++iListIndex;

                //last ----> uiLinkDir
                movedSite = pWalkingTable[movedSite * 2 * uiDim + (uiDim + uiLinkDir)];
                pResult[iListIndex] = pMappingTable[movedSite];
                pResult[iListIndex].m_byDir = i;
                pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, movedSite, i) ? _kDirichlet : 0;
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
    const UINT* __restrict__ pSmallDataTable,
    const BYTE* __restrict__ pBondInfoTable)
{
    intokernalInt4;

    UINT uiDir = _DC_Dir;
    const UINT uiBigSiteIndex = _deviceGetBigIndex(sSite4, pSmallDataTable);
    for (UINT i = 0; i < uiDir; ++i)
    {
        const UINT uiBigIdx = pWalkingTable[uiBigSiteIndex * 2 * uiDir + i];
        pCached[uiSiteIndex * uiDir + i] = pMappingTable[uiBigIdx];
        pCached[uiSiteIndex * uiDir + i].m_byDir = i;
        pCached[uiSiteIndex * uiDir + i].m_byTag = 
            _deviceIsBondDirichlet(pBondInfoTable, uiBigIdx, i) ? _kDirichlet : 0;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCacheFermionMove(SIndex* pCached, 
    const UINT* __restrict__ pWalkingTable,
    const SIndex* __restrict__ pMappingTable,
    const UINT* __restrict__ pSmallDataTable)
{
    intokernalInt4;

    const UINT uiDir = _DC_Dir;
    const UINT uiBigSiteIndex = _deviceGetBigIndex(sSite4, pSmallDataTable);
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

__global__ void _CLG_LAUNCH_BOUND
_kernelCacheEtaMu(BYTE* pCached)
{
    intokernalInt4;
    pCached[uiSiteIndex] =
          ((sSite4.EtaOdd(4) ? 1 : 0) << 4)
        & ((sSite4.EtaOdd(3) ? 1 : 0) << 3)
        & ((sSite4.EtaOdd(2) ? 1 : 0) << 2)
        & ((sSite4.EtaOdd(1) ? 1 : 0) << 1)
        &  (sSite4.EtaOdd(0) ? 1 : 0);
}


__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteCount(UINT* atomic)
{
    intokernal;

    const BYTE plaqCount = static_cast<BYTE>(__idx->m_pSmallData[CIndexData::kPlaqPerSiteIdx]);
    const BYTE plaqLength = static_cast<BYTE>(__idx->m_pSmallData[CIndexData::kPlaqLengthIdx]);

    UINT plaqCountAll = plaqCount * plaqLength;
    UINT toAdd = 0;

    for (BYTE i = 0; i < plaqCount; ++i)
    {
        BYTE byDirichletCount = 0;
        for (BYTE j = 0; j < plaqLength; ++j)
        {
            if (__idx->m_pPlaqutteCache[i * plaqLength + j + uiSiteIndex * plaqCountAll].IsDirichlet())
            {
                ++byDirichletCount;
            }
        }
        if (byDirichletCount < plaqLength)
        {
            ++toAdd;
        }
    }

    atomicAdd(atomic, toAdd);
}

#pragma endregion

UINT CIndexSquare::GetDecompose(UINT volumn)
{
    TArray<UINT> factors = _getFactors(volumn);
    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock();
    UINT maxThreadPerBlock = deviceConstraints[1]; //we only use 1 dimension, so it is the constraint of blockDim.x
#if _CLG_USE_LAUNCH_BOUND
    maxThreadPerBlock = (maxThreadPerBlock > _CLG_LAUNCH_MAX_THREAD) ? _CLG_LAUNCH_MAX_THREAD : maxThreadPerBlock;
#endif
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

    const UINT uiVolumn = biggerLattice.x * biggerLattice.y * biggerLattice.z * biggerLattice.w;
    const UINT threadPerSite = GetDecompose(uiVolumn);
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

    //bake bond infos
    m_pBoundaryCondition->BakeBondInfo(pData->m_pMappingTable, pData->m_pBondInfoTable);

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
        _kernelBakePlaqIndexAtSite << <block, threads >> > (pData->m_pPlaqutteCache, pData->m_pWalkingTable, 
            pData->m_pIndexPositionToSIndex[byFieldId], pData->m_pSmallData, pData->m_pBondInfoTable);

        //bake plaqutte per link
        _kernelBakePlaqIndexAtLink << <block, threads >> > (pData->m_pStappleCache, pData->m_pWalkingTable, 
            pData->m_pIndexPositionToSIndex[byFieldId], pData->m_pSmallData, pData->m_pBondInfoTable);
    }
}

void CIndexSquare::BakeMoveIndex(CIndexData* pData, BYTE byFieldId)
{
    appParanoiac(_T("CIndexSquare::BakeMoveIndex for field ID:%d\n"), byFieldId);

    preparethread;

    checkCudaErrors(cudaMalloc((void**)&pData->m_pGaugeMoveCache[byFieldId], sizeof(SIndex) * _HC_Volume * _HC_Dir));
    checkCudaErrors(cudaMalloc((void**)&pData->m_pFermionMoveCache[byFieldId], sizeof(SIndex) * _HC_Volume * _HC_Dir * 2));

    _kernelCacheGaugeMove << <block, threads >> > (
        pData->m_pGaugeMoveCache[byFieldId], 
        pData->m_pWalkingTable, 
        pData->m_pIndexPositionToSIndex[1], 
        pData->m_pSmallData,
        pData->m_pBondInfoTable);
    _kernelCacheFermionMove << <block, threads >> > (pData->m_pFermionMoveCache[byFieldId], pData->m_pWalkingTable, pData->m_pIndexPositionToSIndex[byFieldId], pData->m_pSmallData);
}

void CIndexSquare::BakeEtaMuTable(class CIndexData* pData)
{
    appParanoiac(_T("CIndexSquare::BakeEtaMuTable\n"));
    checkCudaErrors(cudaMalloc((void**)&pData->m_pEtaMu, sizeof(BYTE) * _HC_Volume));
    preparethread;
    _kernelCacheEtaMu << <block, threads >> > (pData->m_pEtaMu);
}

UINT CIndexSquare::GetPlaqutteCount() const
{
    UINT res[1] = { 0 };
    UINT* deviceRes;
    checkCudaErrors(cudaMalloc((void**)&deviceRes, sizeof(UINT)));
    checkCudaErrors(cudaMemcpy(deviceRes, res, sizeof(UINT), cudaMemcpyHostToDevice));

    preparethread;
    _kernelPlaqutteCount << <block, threads >> > (deviceRes);

    checkCudaErrors(cudaMemcpy(res, deviceRes, sizeof(UINT), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(deviceRes));

    return res[0];
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================