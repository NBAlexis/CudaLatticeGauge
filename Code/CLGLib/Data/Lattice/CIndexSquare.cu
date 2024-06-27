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

/**
 * This record multiply factors to recover site index to pDeviceData
 */
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

    //printf("kPlaqPerSiteIdx=%d kPlaqPerLinkIdx=%d\n", pDeviceData[CIndexData::kPlaqPerSiteIdx], pDeviceData[CIndexData::kPlaqPerLinkIdx]);
}

/**
 * It calculate SSmallInt4 for every big-site
 */
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

/**
 * This bake the plaquttes for one site
 * for each site, there are 6 plaquttes for 4D
 * xy xz xt yz yt zt
 * so the index count is 4x6 = 24
 *
 * for 3D there are 3 plaquttes
 * xy, xz yz
 * so the index count is 4x3 = 12
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelBakePlaqIndexAtSite(SIndex* pResult,
    const SIndex* __restrict__ pLinkTable,
    const UINT* __restrict__ pSmallDataTable,
    const BYTE* __restrict__ pBondInfoTable)
{
    intokernalInt4;

    const UINT uiDim = _DC_Dim;
    const UINT uiBigSiteIndex = _deviceGetBigIndex(sSite4, pSmallDataTable);

    //24
    UINT iListIndex = uiSiteIndex
        * (pSmallDataTable[CIndexData::kPlaqPerSiteIdx] * pSmallDataTable[CIndexData::kPlaqLengthIdx]);

    //Only save plus mu and plus nu
    for (BYTE uiLink = 0; uiLink < uiDim; ++uiLink)
    {
        for (BYTE uiPlaq = uiLink + 1; uiPlaq < uiDim; ++uiPlaq)
        {
            pResult[iListIndex] = SIndex(uiSiteIndex, uiLink);
            pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, uiBigSiteIndex, uiLink) ? _kDirichlet : 0;
            ++iListIndex;

            //start ----> uiLink
            SSmallInt4 sWalking = _deviceSmallInt4OffsetC(sSite4, uiLink + 1);
            UINT uiBigIdx = _deviceGetBigIndex(sWalking, pSmallDataTable);
            //const SIndex& n_p_link__plaq = pLinkTable[uiBigIdx * _DC_Dir + uiPlaq];
            pResult[iListIndex] = pLinkTable[uiBigIdx * _DC_Dir + uiPlaq];
            ++iListIndex;

            //start ----> uiPlaq
            sWalking = _deviceSmallInt4OffsetC(sSite4, uiPlaq + 1);
            uiBigIdx = _deviceGetBigIndex(sWalking, pSmallDataTable);
            const SIndex& n_p_plaq__link = pLinkTable[uiBigIdx * _DC_Dir + uiLink];
            pResult[iListIndex] = n_p_plaq__link;
            pResult[iListIndex].m_byTag = pResult[iListIndex].m_byTag ^ _kDaggerOrOpposite;
            ++iListIndex;

            pResult[iListIndex] = SIndex(uiSiteIndex, uiPlaq);
            pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, uiBigSiteIndex, uiPlaq) ? _kDirichlet : 0;
            pResult[iListIndex].m_byTag = pResult[iListIndex].m_byTag ^ _kDaggerOrOpposite; //Since this is inside lattice, | or ^ are both OK
            ++iListIndex;
        }
    }
}

/**
 * This bake the staple
 * for each link, there are 6 plaquttes for 4D
 * for example, x-link
 * the staple is xy, xz, xt for forward and backward
 *
 * for 3D there are 4 plaquttes
 * xy, xz for forward and backward
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelBakePlaqIndexAtLink(SIndex* pResult,
    const SIndex* __restrict__ pLinkTable,
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
                SSmallInt4 sWalking = _deviceSmallInt4OffsetC(sSite4, i + 1);
                UINT uiBigIdx = _deviceGetBigIndex(sWalking, pSmallDataTable);
                //UINT movedSite = pWalkingTable[uiBigSiteIndex * 2 * uiDim + (uiDim + i)];
                //pResult[iListIndex] = pMappingTable[movedSite];
                //pResult[iListIndex].m_byDir = uiLinkDir;
                //pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, movedSite, uiLinkDir) ? _kDirichlet : 0;
                pResult[iListIndex] = pLinkTable[uiBigIdx * _DC_Dir + uiLinkDir];
                ++iListIndex;

                //start ---> uiLinkDir
                //movedSite = pWalkingTable[uiBigSiteIndex * 2 * uiDim + (uiDim + uiLinkDir)];
                //pResult[iListIndex] = pMappingTable[movedSite];
                //pResult[iListIndex].m_byDir = i;
                //pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, movedSite, i) ? _kDirichlet : 0;
                //pResult[iListIndex].m_byTag |= _kDaggerOrOpposite;
                sWalking = _deviceSmallInt4OffsetC(sSite4, uiLinkDir + 1);
                uiBigIdx = _deviceGetBigIndex(sWalking, pSmallDataTable);
                pResult[iListIndex] = pLinkTable[uiBigIdx * _DC_Dir + i];
                pResult[iListIndex].m_byTag = pResult[iListIndex].m_byTag ^ _kDaggerOrOpposite;
                ++iListIndex;

                //=============================================
                //add backward
                //[site-p_dir][p_dir]^-1, [site-p_dir][b_dir], [site-p_dir+b_dir][p_dir]
                //i <---- start
                //movedSite = pWalkingTable[uiBigSiteIndex * 2 * uiDim + i];
                //pResult[iListIndex] = pMappingTable[movedSite];
                //pResult[iListIndex].m_byDir = i;
                //pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, movedSite, i) ? _kDirichlet : 0;
                //pResult[iListIndex].m_byTag |= _kDaggerOrOpposite;
                sWalking = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(i) - 1);
                uiBigIdx = _deviceGetBigIndex(sWalking, pSmallDataTable);
                pResult[iListIndex] = pLinkTable[uiBigIdx * _DC_Dir + i];
                pResult[iListIndex].m_byTag = pResult[iListIndex].m_byTag ^ _kDaggerOrOpposite;
                ++iListIndex;

                //
                //pResult[iListIndex] = SIndex(pResult[iListIndex - 1].m_uiSiteIndex, uiLinkDir);
                //pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, movedSite, uiLinkDir) ? _kDirichlet : 0;
                pResult[iListIndex] = pLinkTable[uiBigIdx * _DC_Dir + uiLinkDir];
                ++iListIndex;

                //last ----> uiLinkDir
                //movedSite = pWalkingTable[movedSite * 2 * uiDim + (uiDim + uiLinkDir)];
                //pResult[iListIndex] = pMappingTable[movedSite];
                //pResult[iListIndex].m_byDir = i;
                //pResult[iListIndex].m_byTag = _deviceIsBondDirichlet(pBondInfoTable, movedSite, i) ? _kDirichlet : 0;
                sWalking = _deviceSmallInt4OffsetC(sWalking, uiLinkDir + 1);
                uiBigIdx = _deviceGetBigIndex(sWalking, pSmallDataTable);
                pResult[iListIndex] = pLinkTable[uiBigIdx * _DC_Dir + i];
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
    //const UINT* __restrict__ pWalkingTable,
    //const SIndex* __restrict__ pMappingTable,
    const SIndex* __restrict__ pLinkTable,
    const UINT* __restrict__ pSmallDataTable
    //const BYTE* __restrict__ pBondInfoTable
)
{
    intokernalInt4;

    UINT uiDir = _DC_Dir;
    //const UINT uiBigSiteIndex = _deviceGetBigIndex(sSite4, pSmallDataTable);
    for (UINT i = 0; i < uiDir; ++i)
    {
        const SSmallInt4 sWalking = _deviceSmallInt4OffsetC(sSite4, -static_cast<INT>(i) - 1);
        const UINT uiBigIdx = _deviceGetBigIndex(sWalking, pSmallDataTable);
        //const UINT uiBigIdx = pWalkingTable[uiBigSiteIndex * 2 * uiDir + i];
        //pCached[uiSiteIndex * uiDir + i] = pMappingTable[uiBigIdx];
        //pCached[uiSiteIndex * uiDir + i].m_byDir = i;
        //pCached[uiSiteIndex * uiDir + i].m_byTag = 
        //    _deviceIsBondDirichlet(pBondInfoTable, uiBigIdx, i) ? _kDirichlet : 0;
        pCached[uiSiteIndex * uiDir + i] = pLinkTable[uiBigIdx * uiDir + i];
        pCached[uiSiteIndex * uiDir + i].m_byTag = pCached[uiSiteIndex * uiDir + i].m_byTag ^ _kDaggerOrOpposite;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCacheFermionMove(SIndex* pCached, 
    const SIndex* __restrict__ pMappingTable,
    const UINT* __restrict__ pSmallDataTable)
{
    intokernalInt4;

    for (UINT i = 0; i < _DC_Dir; ++i)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, i);
        //first element is right, second element is left.
        pCached[linkIndex * 2]
            = pMappingTable[_deviceGetBigIndex(_deviceSmallInt4OffsetC(sSite4, __fwd(i)), pSmallDataTable)];
        pCached[linkIndex * 2 + 1]
            = pMappingTable[_deviceGetBigIndex(_deviceSmallInt4OffsetC(sSite4, __bck(i)), pSmallDataTable)];
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCacheEtaMu(BYTE* pCached)
{
    intokernalInt4;
    pCached[uiSiteIndex] =
          ((sSite4.EtaOdd(4) ? 1 : 0) << 4)
        | ((sSite4.EtaOdd(3) ? 1 : 0) << 3)
        | ((sSite4.EtaOdd(2) ? 1 : 0) << 2)
        | ((sSite4.EtaOdd(1) ? 1 : 0) << 1)
        |  (sSite4.EtaOdd(0) ? 1 : 0);
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
    //_kernalBakeWalkingTable << <blocks, threads >> > (pData->m_pWalkingTable, biggerLatticeMod);

    //bake index mappings
    for (BYTE i = 1; i < kMaxFieldCount; ++i)
    {
        const CField* pField = appGetLattice()->GetFieldById(i);
        if (NULL != pField)
        {
            checkCudaErrors(cudaMalloc((void**)&pData->m_pIndexPositionToSIndex[i], sizeof(SIndex)
                * (_HC_Lx + 2 * CIndexData::kCacheIndexEdge) * (_HC_Ly + 2 * CIndexData::kCacheIndexEdge)
                * (_HC_Lz + 2 * CIndexData::kCacheIndexEdge) * (_HC_Lt + 2 * CIndexData::kCacheIndexEdge)
                ));

            //bake map index
            _kernalBakeMappingTable << <blocks, threads >> > (pData->m_pMappingTable, biggerLatticeMod);

            //bake boundary condition
            m_pBoundaryCondition->BakeEdgePoints(i, pData->m_pMappingTable, pData->m_pIndexPositionToSIndex[i]);

            if (pField->IsGaugeField())
            {
                checkCudaErrors(cudaMalloc((void**)&pData->m_pIndexLinkToSIndex[i], sizeof(SIndex)
                    * (_HC_Lx + 2 * CIndexData::kCacheIndexEdge) * (_HC_Ly + 2 * CIndexData::kCacheIndexEdge)
                    * (_HC_Lz + 2 * CIndexData::kCacheIndexEdge) * (_HC_Lt + 2 * CIndexData::kCacheIndexEdge)
                    * _HC_Dir
                ));

                m_pBoundaryCondition->BakeBondGlue(i, pData->m_pMappingTable, pData->m_pIndexLinkToSIndex[i]);
            }
        }
    }

    //bake bond infos
    m_pBoundaryCondition->BakeBondInfo(pData->m_pMappingTable, pData->m_pBondInfoTable);

    //bake region id table
    m_pBoundaryCondition->BakeRegionTable(pData->m_byRegionTable);

    //copy the mapping table to device
    checkCudaErrors(cudaMemcpy(pData->m_pDeviceIndexPositionToSIndex, pData->m_pIndexPositionToSIndex, sizeof(SIndex*) * kMaxFieldCount, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(pData->m_pDeviceIndexLinkToSIndex, pData->m_pIndexLinkToSIndex, sizeof(SIndex*) * kMaxFieldCount, cudaMemcpyHostToDevice));
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
        _kernelBakePlaqIndexAtSite << <block, threads >> > (
            pData->m_pPlaqutteCache, pData->m_pIndexLinkToSIndex[byFieldId], 
            pData->m_pSmallData, pData->m_pBondInfoTable);
        checkCudaErrors(cudaDeviceSynchronize());
        //bake plaqutte per link
        _kernelBakePlaqIndexAtLink << <block, threads >> > (
            pData->m_pStappleCache, pData->m_pIndexLinkToSIndex[byFieldId], 
            pData->m_pSmallData, pData->m_pBondInfoTable);
        checkCudaErrors(cudaDeviceSynchronize());
    }
}

void CIndexSquare::BakeMoveIndex(CIndexData* pData, BYTE byFieldId)
{
    appParanoiac(_T("CIndexSquare::BakeMoveIndex for field ID:%d\n"), byFieldId);

    preparethread;

    checkCudaErrors(cudaMalloc((void**)&pData->m_pGaugeMoveCache[byFieldId], sizeof(SIndex) * _HC_Volume * _HC_Dir));
    checkCudaErrors(cudaMalloc((void**)&pData->m_pMoveCache[byFieldId], sizeof(SIndex) * _HC_Volume * _HC_Dir * 2));

    const CField* defaultGauge = appGetLattice()->GetFieldById(1);
    if (NULL != defaultGauge && defaultGauge->IsGaugeField())
    {
        _kernelCacheGaugeMove << <block, threads >> > (
            pData->m_pGaugeMoveCache[byFieldId],
            pData->m_pIndexLinkToSIndex[1],
            pData->m_pSmallData);
        checkCudaErrors(cudaDeviceSynchronize());
    }

    _kernelCacheFermionMove << <block, threads >> > (pData->m_pMoveCache[byFieldId], pData->m_pIndexPositionToSIndex[byFieldId], pData->m_pSmallData);
    checkCudaErrors(cudaDeviceSynchronize());
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