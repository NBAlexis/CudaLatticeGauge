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

#undef preparethread
#undef intokernal
#undef intokernaldir
#undef intokernalInt4


#define preparethreadIndex \
const dim3 block(_HC_DecompX, _HC_DecompY, _HC_DecompZ); \
const dim3 threads(_HC_DecompLx, _HC_DecompLy, _HC_DecompLz);

#define intokernalInt4Index \
SSmallInt4 sSite4; \
const UINT _ixy = (threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.x = static_cast<SCHAR> (_ixy / _DC_Ly); \
sSite4.y = static_cast<SCHAR> (_ixy % _DC_Ly); \
sSite4.z = static_cast<SCHAR>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.w = static_cast<SCHAR>(threadIdx.z + blockIdx.z * blockDim.z); \
const UINT uiSiteIndex = _ixy * _DC_GridDimZT + sSite4.z * _DC_Lt + sSite4.w; 

#define intokernalInt4dirIndex \
SSmallInt4 sSite4; \
const UINT _ixy = (threadIdx.x + blockIdx.x * blockDim.x); \
sSite4.x = static_cast<SCHAR> (_ixy / _DC_Ly); \
sSite4.y = static_cast<SCHAR> (_ixy % _DC_Ly); \
sSite4.z = static_cast<SCHAR>(threadIdx.y + blockIdx.y * blockDim.y); \
sSite4.w = static_cast<SCHAR>(threadIdx.z + blockIdx.z * blockDim.z); \
const UINT uiSiteIndex = _ixy * _DC_GridDimZT + sSite4.z * _DC_Lt + sSite4.w; \
const BYTE uiDir = static_cast<BYTE>(_DC_Dir);

#define intokernalIndex \
const UINT uiSiteIndex = ((threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z)); 



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
_kernalBakeMappingTable(SSmallInt4* pDeviceData, SSmallInt4* pSiteDeviceData, uint3 mods, UINT maxuisiteindex)
{
    UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
    SSmallInt4 coord;
    coord.x = static_cast<SCHAR>(idxAll / mods.x) - CIndexData::kCacheIndexEdge;
    coord.y = static_cast<SCHAR>((idxAll % mods.x) / mods.y) - CIndexData::kCacheIndexEdge;
    coord.z = static_cast<SCHAR>((idxAll % mods.y) / mods.z) - CIndexData::kCacheIndexEdge;
    coord.w = static_cast<SCHAR>(idxAll % mods.z) - CIndexData::kCacheIndexEdge;

    pDeviceData[idxAll] = coord;
    if (idxAll < maxuisiteindex)
    {
        pSiteDeviceData[idxAll] = __deviceSiteIndexToInt4Baking(idxAll);
    }
}

#if _CLG_MULTI_GPU
/**
 * Multi-GPU (Phase 1 / P4-5): build the halo gather map. One thread per
 * big-index cell; a cell that redirects to a halo slot (face/edge/corner per
 * P4-5.2's _deviceHaloRedirectSite) records the LOCAL source site at that slot.
 * This reuses the exact same _deviceHaloRedirectSite decision the edge/glue
 * bake uses, so slot numbering can never diverge. Cells that do not redirect
 * write nothing. The source-coordinate override below is direction-generic:
 * every split direction that is out of range is moved to its near-boundary
 * plane, so face (1 dir), edge (2 dirs) and corner (3 dirs) cells all resolve
 * to the near-boundary plane/edge/corner the neighbour rank needs to receive.
 */
__global__ void _CLG_LAUNCH_BOUND
_kernalBakeHaloGatherSite(const SSmallInt4* __restrict__ pMapping, SIndex* pHaloGatherIndex)
{
    const UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
    const SSmallInt4 raw(pMapping[idxAll]);
    const INT iRaw[4] = { raw.x, raw.y, raw.z, raw.w };

    INT iWrapped[4] = { raw.x, raw.y, raw.z, raw.w };
    for (UINT d = 0; d < 4; ++d)
    {
        const INT iL = _constIntegers[ECI_Lx + d];
        while (iWrapped[d] < 0) { iWrapped[d] += iL; }
        while (iWrapped[d] >= iL) { iWrapped[d] -= iL; }
    }

    const UINT uiGrid[4] = { _DC_GpuGridX, _DC_GpuGridY, _DC_GpuGridZ, _DC_GpuGridT };
    const UINT uiLocalL[4] = { static_cast<UINT>(_constIntegers[ECI_Lx]),
        static_cast<UINT>(_constIntegers[ECI_Ly]), static_cast<UINT>(_constIntegers[ECI_Lz]),
        static_cast<UINT>(_constIntegers[ECI_Lt]) };
    UINT uiHaloSlot = 0;
    if (EHR_Halo == _deviceHaloRedirectSite(uiLocalL, uiGrid, _DC_Volume, _DC_HaloWidth,
        iRaw, iWrapped, uiHaloSlot))
    {
        //A halo slot is used both as this rank's SEND buffer for a face and (after
        //MPI) as the RECV storage the kernel reads. The data a neighbour needs is
        //this rank's NEAR-boundary plane in the split direction, NOT the periodic
        //wrap (which is the far plane). For the split dir the periodic wrap mirrors
        //the plane to the wrong side, so override the split-direction source
        //coordinate to the near-boundary plane. Transverse axes keep their in-range
        //(wrapped == identity) coordinate so the slot's transverse cell matches.
        INT iSrc[4] = { iWrapped[0], iWrapped[1], iWrapped[2], iWrapped[3] };
        for (UINT d = 0; d < 4; ++d)
        {
            if (uiGrid[d] <= 1)
            {
                continue; //Not split -> no halo face in this dir; keep wrap.
            }
            const INT iL = static_cast<INT>(uiLocalL[d]);
            if (iRaw[d] < 0)
            {
                //neg face, layer = -1 - raw: near-neg plane is local coord = layer.
                iSrc[d] = -1 - iRaw[d];
            }
            else if (iRaw[d] >= iL)
            {
                //pos face, layer = raw - iL: near-pos plane is local coord = iL-1-layer.
                iSrc[d] = iL - 1 - (iRaw[d] - iL);
            }
        }
        SSmallInt4 src;
        src.x = static_cast<SCHAR>(iSrc[0]);
        src.y = static_cast<SCHAR>(iSrc[1]);
        src.z = static_cast<SCHAR>(iSrc[2]);
        src.w = static_cast<SCHAR>(iSrc[3]);
        const UINT uiSource = _deviceGetSiteIndex(src);
        //Improve-1 (3.7): the gather source must always resolve to a legal
        //interior site of this rank -- never a halo slot, never invalid.
        assert(uiSource < _DC_Volume);
        pHaloGatherIndex[uiHaloSlot] = SIndex(uiSource);
    }
}
#endif

/**
 * Multi-GPU (Phase 1, Improve-1 I8/3.7): bake the 32-bit global coordinate
 * table, one thread per big-index cell. Cells fully inside the local lattice
 * write their local slot [0, _DC_Volume); cells redirected to a halo slot
 * (the exact same _deviceHaloRedirectSite decision as the gather-map bake, so
 * slot numbering can never diverge) write slot _DC_Volume + haloSlot. The
 * global coordinate is this rank's offset plus the raw (possibly
 * out-of-local-range) coordinate, wrapped to the canonical [0, GlobalL) per
 * axis. Baked on BOTH single- and multi-GPU (the halo part is empty on
 * single-GPU), so _deviceSIndexToGlobalInt4 is a plain lookup that never
 * narrows through SCHAR.
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelBakeGlobalCoordinateTable(const SSmallInt4* __restrict__ pMapping, SInt4* pTable)
{
    const UINT idxAll = threadIdx.x + blockDim.x * blockIdx.x;
    const SSmallInt4 raw(pMapping[idxAll]);
    const INT iRaw[4] = { raw.x, raw.y, raw.z, raw.w };
    const UINT uiLocalL[4] = { _DC_Lx, _DC_Ly, _DC_Lz, _DC_Lt };
    const UINT uiOffset[4] = { _DC_OffsetX, _DC_OffsetY, _DC_OffsetZ, _DC_OffsetT };

    UBOOL bLocal = TRUE;
    for (UINT d = 0; d < 4; ++d)
    {
        if (iRaw[d] < 0 || iRaw[d] >= static_cast<INT>(uiLocalL[d]))
        {
            bLocal = FALSE;
            break;
        }
    }
    if (bLocal)
    {
        pTable[_deviceGetSiteIndex(raw)] = SInt4(
            static_cast<INT>(uiOffset[0]) + iRaw[0],
            static_cast<INT>(uiOffset[1]) + iRaw[1],
            static_cast<INT>(uiOffset[2]) + iRaw[2],
            static_cast<INT>(uiOffset[3]) + iRaw[3]);
        return;
    }

#if _CLG_MULTI_GPU
    INT iWrapped[4] = { iRaw[0], iRaw[1], iRaw[2], iRaw[3] };
    for (UINT d = 0; d < 4; ++d)
    {
        const INT iL = static_cast<INT>(uiLocalL[d]);
        while (iWrapped[d] < 0) { iWrapped[d] += iL; }
        while (iWrapped[d] >= iL) { iWrapped[d] -= iL; }
    }

    const UINT uiGrid[4] = { _DC_GpuGridX, _DC_GpuGridY, _DC_GpuGridZ, _DC_GpuGridT };
    const UINT uiGlobalL[4] = { _DC_GlobalLx, _DC_GlobalLy, _DC_GlobalLz, _DC_GlobalLt };
    UINT uiHaloSlot = 0;
    if (EHR_Halo == _deviceHaloRedirectSite(uiLocalL, uiGrid, _DC_Volume, _DC_HaloWidth,
        iRaw, iWrapped, uiHaloSlot))
    {
        INT iGlobal[4];
        for (UINT d = 0; d < 4; ++d)
        {
            INT iCoord = static_cast<INT>(uiOffset[d]) + iRaw[d];
            const INT iG = static_cast<INT>(uiGlobalL[d]);
            while (iCoord < 0) { iCoord += iG; }
            while (iCoord >= iG) { iCoord -= iG; }
            iGlobal[d] = iCoord;
        }
        pTable[_DC_Volume + uiHaloSlot] = SInt4(iGlobal[0], iGlobal[1], iGlobal[2], iGlobal[3]);
    }
#endif
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
    const UINT* __restrict__ pSmallDataTable
    //, const BYTE* __restrict__ pBondInfoTable
)
{
    intokernalInt4Index;

    const UINT uiDim = _DC_Dim;
    const UINT uiBigSiteIndex = _deviceGetBigIndex(sSite4, pSmallDataTable);

    //24
    UINT iListIndex = uiSiteIndex
        * (pSmallDataTable[CIndexData::kPlaqPerSiteIdx] * pSmallDataTable[CIndexData::kPlaqLengthIdx]);

    //Only save plus mu and plus nu
    #pragma unroll
    for (BYTE uiLink = 0; uiLink < uiDim; ++uiLink)
    {
        #pragma unroll
        for (BYTE uiPlaq = uiLink + 1; uiPlaq < uiDim; ++uiPlaq)
        {
            pResult[iListIndex] = SIndex(uiSiteIndex, uiLink);
            pResult[iListIndex].m_byTag = pLinkTable[uiBigSiteIndex * _DC_Dir + uiLink].IsDirichlet() ? _kDirichlet : 0;
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
            pResult[iListIndex].m_byTag = pLinkTable[uiBigSiteIndex * _DC_Dir + uiPlaq].IsDirichlet() ? _kDirichlet : 0; //_deviceIsBondDirichlet(pBondInfoTable, uiBigSiteIndex, uiPlaq) ? _kDirichlet : 0;
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
    const UINT* __restrict__ pSmallDataTable
    //, const BYTE* __restrict__ pBondInfoTable
)
{
    intokernalInt4Index;

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
                pResult[iListIndex].m_byTag = pLinkTable[uiBigSiteIndex * _DC_Dir + i].IsDirichlet()  ? _kDirichlet : 0;  // _deviceIsBondDirichlet(pBondInfoTable, uiBigSiteIndex, i) ? _kDirichlet : 0;
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
    intokernalInt4Index;

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
    intokernalInt4Index;

    for (UINT i = 0; i < _DC_Dir; ++i)
    {
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, i);
        //first element is right, second element is left.
        pCached[linkIndex * 2]
            = pMappingTable[_deviceGetBigIndex(_deviceSmallInt4OffsetC(sSite4, __fwd(i)), pSmallDataTable)];
        pCached[linkIndex * 2].m_byDir = static_cast<BYTE>(i);
        pCached[linkIndex * 2 + 1]
            = pMappingTable[_deviceGetBigIndex(_deviceSmallInt4OffsetC(sSite4, __bck(i)), pSmallDataTable)];
        pCached[linkIndex * 2 + 1].m_byDir = static_cast<BYTE>(i);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCacheEtaMu(BYTE* pCached)
{
    intokernalInt4Index;
    pCached[uiSiteIndex] =
          ((sSite4.EtaOdd(4) ? 1 : 0) << 4)
        | ((sSite4.EtaOdd(3) ? 1 : 0) << 3)
        | ((sSite4.EtaOdd(2) ? 1 : 0) << 2)
        | ((sSite4.EtaOdd(1) ? 1 : 0) << 1)
        |  (sSite4.EtaOdd(0) ? 1 : 0);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelCacheNaikIndex(SIndex* pCached, 
    const SIndex* __restrict__ pMappingTable,
    const UINT* __restrict__ pSmallDataTable)
{
    intokernalInt4dirIndex;
    SCHAR move[3];

    for (BYTE dir = 0; dir < uiDir; ++dir)
    {
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);

        SSmallInt4 moveForward = sSite4;
        move[0] = dir + 1;
        move[1] = dir + 1;
        move[2] = dir + 1;
        _deviceSmallInt4Offset(moveForward, move, 3);
        pCached[linkIndex * 2] = pMappingTable[_deviceGetBigIndex(moveForward, pSmallDataTable)];
        pCached[linkIndex * 2].m_byDir = dir;

        moveForward = sSite4;
        move[0] = -static_cast<SCHAR>(dir) - 1;
        move[1] = -static_cast<SCHAR>(dir) - 1;
        move[2] = -static_cast<SCHAR>(dir) - 1;
        _deviceSmallInt4Offset(moveForward, move, 3);
        pCached[linkIndex * 2 + 1] = pMappingTable[_deviceGetBigIndex(moveForward, pSmallDataTable)];
        pCached[linkIndex * 2 + 1].m_byDir = dir;
    }
}


__global__ void _CLG_LAUNCH_BOUND
_kernelPlaqutteCount(UINT* atomic, BYTE byFieldId)
{
    intokernalIndex;

#if !_CLG_ASSUME_SQUARE_LATTICE
    const BYTE plaqCountPerSite = static_cast<BYTE>(__idx->m_pSmallData[CIndexData::kPlaqPerSiteIdx]);
    const BYTE plaqLength = static_cast<BYTE>(__idx->m_pSmallData[CIndexData::kPlaqLengthIdx]);
#endif

    UINT plaqCountAll = plaqCountPerSite * plaqLength;
    UINT toAdd = 0;

    for (BYTE i = 0; i < plaqCountPerSite; ++i)
    {
        BYTE byDirichletCount = 0;
        for (BYTE j = 0; j < plaqLength; ++j)
        {
            if (__idx->m_pPlaqutteCache[byFieldId][i * plaqLength + j + uiSiteIndex * plaqCountAll].IsDirichlet())
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

//__global__ void _CLG_LAUNCH_BOUND
//_kernelBakeEvenOdd(UINT* res, UBOOL bEven, UINT uiHalfVolumn)
//{
//    const UINT uiHalfSiteIndex = threadIdx.x + blockIdx.x * blockDim.x;
//    if (uiHalfSiteIndex >= uiHalfVolumn)
//    {
//        return;
//    }
//
//    if (bEven)
//    {
//        UINT uiResSiteIndex = 2 * uiHalfSiteIndex;
//        SSmallInt4 site4 = __deviceSiteIndexToInt4Baking(uiResSiteIndex);
//        if (site4.IsOdd())
//        {
//            ++uiResSiteIndex;
//        }
//        res[uiHalfSiteIndex] = uiResSiteIndex;
//    }
//    else
//    {
//        UINT uiResSiteIndex = 2 * uiHalfSiteIndex + 1;
//        SSmallInt4 site4 = __deviceSiteIndexToInt4Baking(uiResSiteIndex);
//        if (!site4.IsOdd())
//        {
//            --uiResSiteIndex;
//        }
//        res[uiHalfSiteIndex + uiHalfVolumn] = uiResSiteIndex;
//    }
//}

#pragma endregion

UINT CIndexSquare::GetDecompose(UINT volumn)
{
    TArray<UINT> factors = _getFactors(volumn);
    TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock(appGetDeviceId());
    UINT maxThreadPerBlock = min(deviceConstraints[0], deviceConstraints[1]); //we only use 1 dimension, so it is the constraint of blockDim.x
    maxThreadPerBlock = (maxThreadPerBlock > CCommonData::m_uiMaxThreadPerBlock) ? CCommonData::m_uiMaxThreadPerBlock : maxThreadPerBlock;
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

#if !_CLG_ASSUME_SQUARE_LATTICE
    pData->m_uiPlaqutteLength = 4;
    pData->m_uiPlaqutteCountPerSite = static_cast<BYTE>(_HC_Dim * (_HC_Dim - 1) / 2);
    pData->m_uiPlaqutteCountPerLink = static_cast<BYTE>(2 * (_HC_Dim - 1));
#endif

    //bake small data
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    _LAUNCH_KERNEL(_kernalBakeSmallData, 1, 1, pData->m_pSmallData);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
    //bake walking index
    //_LAUNCH_KERNEL(_kernalBakeWalkingTable, blocks, threads, pData->m_pWalkingTable, biggerLatticeMod);

    UBOOL indextableEverBaked = FALSE;

    //bake index mappings
    for (BYTE i = 1; i < kMaxFieldCount; ++i)
    {
        const CField* pField = appGetLattice()->GetFieldById(i);
        if (NULL != pField)
        {
            checkCudaErrors(__cudaMalloc((void**)&pData->m_pIndexPositionToSIndex[i], sizeof(SIndex)
                * (_HC_Lx + 2 * CIndexData::kCacheIndexEdge) * (_HC_Ly + 2 * CIndexData::kCacheIndexEdge)
                * (_HC_Lz + 2 * CIndexData::kCacheIndexEdge) * (_HC_Lt + 2 * CIndexData::kCacheIndexEdge)
                ));

            if (!indextableEverBaked)
            {
                indextableEverBaked = TRUE;
                //bake map index
                _LAUNCH_KERNEL(_kernalBakeMappingTable, blocks, threads, pData->m_pMappingTable, pData->m_pSiteMappingTable, biggerLatticeMod,
                    _HC_Lx * _HC_Ly * _HC_Lz * _HC_Lt);
            }

            //bake boundary condition
            m_pBoundaryCondition->BakeEdgePoints(i, pData->m_pMappingTable, pData->m_pIndexPositionToSIndex[i]);
            checkCudaErrors(cudaDeviceSynchronize());
            checkCudaErrors(cudaGetLastError());

            if (pField->IsGaugeField())
            {
                checkCudaErrors(__cudaMalloc((void**)&pData->m_pIndexLinkToSIndex[i], sizeof(SIndex)
                    * (_HC_Lx + 2 * CIndexData::kCacheIndexEdge) * (_HC_Ly + 2 * CIndexData::kCacheIndexEdge)
                    * (_HC_Lz + 2 * CIndexData::kCacheIndexEdge) * (_HC_Lt + 2 * CIndexData::kCacheIndexEdge)
                    * _HC_Dir
                ));

                //checkCudaErrors(cudaMalloc((void**)&pData->m_pBondInfoTable[i], sizeof(BYTE)
                //    * (_HC_Lx + 2 * CIndexData::kCacheIndexEdge) * (_HC_Ly + 2 * CIndexData::kCacheIndexEdge)
                //    * (_HC_Lz + 2 * CIndexData::kCacheIndexEdge) * (_HC_Lt + 2 * CIndexData::kCacheIndexEdge)
                //    * _HC_Dir
                //));

                m_pBoundaryCondition->BakeBondGlue(i, pData->m_pMappingTable, pData->m_pIndexLinkToSIndex[i]);
                checkCudaErrors(cudaDeviceSynchronize());
                checkCudaErrors(cudaGetLastError());
                //m_pBoundaryCondition->BakeBondInfo(pData->m_pMappingTable, pData->m_pBondInfoTable[i], i);
                //checkCudaErrors(cudaDeviceSynchronize());
            }
        }
    }

    //bake bond infos
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    //bake region id table
    m_pBoundaryCondition->BakeRegionTable(pData->m_byRegionTable);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

#if _CLG_MULTI_GPU
    //Multi-GPU (Phase 1, Improve-1 I8): bake the halo gather map (halo slot ->
    //SIndex of the local interior source site). Sized to the layout's halo slots;
    //0 (and skipped) when nothing is split, keeping single-GPU
    //allocation/behaviour unchanged.
    pData->m_uiHaloSiteCount = _HC_HaloSiteCount();
    if (pData->m_uiHaloSiteCount > 0)
    {
        checkCudaErrors(__cudaMalloc((void**)&pData->m_pHaloGatherIndex,
            sizeof(SIndex) * pData->m_uiHaloSiteCount));
        //Default every slot to site 0 so an unfilled slot is in-range (defensive;
        //every valid slot is overwritten by the kernel below).
        checkCudaErrors(cudaMemset(pData->m_pHaloGatherIndex, 0,
            sizeof(SIndex) * pData->m_uiHaloSiteCount));
        _LAUNCH_KERNEL(_kernalBakeHaloGatherSite, blocks, threads,
            pData->m_pMappingTable, pData->m_pHaloGatherIndex);
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaGetLastError());
    }
#endif

    //Multi-GPU (Phase 1, Improve-1 I8/3.7): bake the 32-bit global coordinate
    //table -- local slots [0, V) then halo slots [V, V + haloSiteCount). Baked
    //on BOTH single- and multi-GPU (m_uiHaloSiteCount stays 0 on single-GPU),
    //so _deviceSIndexToGlobalInt4 is a plain lookup.
    checkCudaErrors(__cudaMalloc((void**)&pData->m_pGlobalCoordinateTable,
        sizeof(SInt4) * (_HC_Volume + pData->m_uiHaloSiteCount)));
    checkCudaErrors(cudaMemset(pData->m_pGlobalCoordinateTable, 0,
        sizeof(SInt4) * (_HC_Volume + pData->m_uiHaloSiteCount)));
    _LAUNCH_KERNEL(_kernelBakeGlobalCoordinateTable, blocks, threads,
        pData->m_pMappingTable, pData->m_pGlobalCoordinateTable);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());

    //copy the mapping table to device
    checkCudaErrors(cudaMemcpy(pData->m_pDeviceIndexPositionToSIndex, pData->m_pIndexPositionToSIndex, sizeof(SIndex*) * kMaxFieldCount, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(pData->m_pDeviceIndexLinkToSIndex, pData->m_pIndexLinkToSIndex, sizeof(SIndex*) * kMaxFieldCount, cudaMemcpyHostToDevice));
}

void CIndexSquare::BakePlaquttes(CIndexData* pData, BYTE byFieldId)
{
    //If Has Gauge field
    if (NULL != appGetLattice()->GetFieldById(byFieldId))
    {
        preparethreadIndex;

        checkCudaErrors(__cudaMalloc((void**)&pData->m_pPlaqutteCache[byFieldId], sizeof(SIndex) * _HC_Volume * (_HC_Dim * (_HC_Dim - 1) / 2) * 4));
        checkCudaErrors(__cudaMalloc((void**)&pData->m_pStappleCache[byFieldId], sizeof(SIndex) * _HC_Volume * _HC_Dim * (2 * (_HC_Dim - 1)) * 3));
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaGetLastError());

        //bake plaqutte per site       
        _LAUNCH_KERNEL(_kernelBakePlaqIndexAtSite, block, threads, 
            pData->m_pPlaqutteCache[byFieldId], pData->m_pIndexLinkToSIndex[byFieldId],
            pData->m_pSmallData); // , pData->m_pBondInfoTable[byFieldId]);
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaGetLastError());

        //The register for this function is 80, so use half of the thread but double of the block
        TArray<UINT> deviceConstraints = CCudaHelper::GetMaxThreadCountAndThreadPerblock(GCLGManager.m_iDeviceId);
        deviceConstraints[0] = deviceConstraints[0] / 2;
        TArray<UINT> latticeDim;
        latticeDim.AddItem(_HC_Lx * _HC_Ly);
        latticeDim.AddItem(_HC_Lz);
        latticeDim.AddItem(_HC_Lt);
        TArray <UINT> decomp = _getDecompose(deviceConstraints, latticeDim);
        const dim3 blockhalf(decomp[0], decomp[1], decomp[2]);
        const dim3 threadshalf(decomp[3], decomp[4], decomp[5]);
        //bake plaqutte per link
        _LAUNCH_KERNEL(_kernelBakePlaqIndexAtLink, blockhalf, threadshalf, 
            pData->m_pStappleCache[byFieldId], pData->m_pIndexLinkToSIndex[byFieldId],
            pData->m_pSmallData); // , pData->m_pBondInfoTable[byFieldId]);
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaGetLastError());
    }
}

void CIndexSquare::BakeMoveIndex(CIndexData* pData, BYTE byFieldId)
{
    appParanoiac(_T("CIndexSquare::BakeMoveIndex for field ID:%d\n"), byFieldId);

    preparethreadIndex;

    checkCudaErrors(__cudaMalloc((void**)&pData->m_pGaugeMoveCache[byFieldId], sizeof(SIndex) * _HC_Volume * _HC_Dir));
    checkCudaErrors(__cudaMalloc((void**)&pData->m_pMoveCache[byFieldId], sizeof(SIndex) * _HC_Volume * _HC_Dir * 2));

    const CField* defaultGauge = appGetLattice()->GetFieldById(1);
    if (NULL != defaultGauge && defaultGauge->IsGaugeField())
    {
        _LAUNCH_KERNEL(_kernelCacheGaugeMove, block, threads, 
            pData->m_pGaugeMoveCache[byFieldId],
            pData->m_pIndexLinkToSIndex[1],
            pData->m_pSmallData);
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaGetLastError());
    }

    _LAUNCH_KERNEL(_kernelCacheFermionMove, block, threads, pData->m_pMoveCache[byFieldId], pData->m_pIndexPositionToSIndex[byFieldId], pData->m_pSmallData);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

void CIndexSquare::BakeEtaMuTable(class CIndexData* pData)
{
    appParanoiac(_T("CIndexSquare::BakeEtaMuTable\n"));
    checkCudaErrors(__cudaMalloc((void**)&pData->m_pEtaMu, sizeof(BYTE) * _HC_Volume));
    preparethreadIndex;
    _LAUNCH_KERNEL(_kernelCacheEtaMu, block, threads, pData->m_pEtaMu);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

//void CIndexSquare::BakeEvenOddTable(class CIndexData* pData)
//{
//    //since only when the number of sites is even, when we use even-odd so we only bake this table when sites are even
//    if (_HC_Volume & 1)
//    {
//        appGeneral(_T("The number of sites is odd, not support Even-Odd\n"));
//        return;
//    }
//    const UINT uiVolumn = _HC_VolumeHalf;
//    _LAUNCH_KERNEL(_kernelBakeEvenOdd, _HC_DecompBlockHalf, _HC_DecompThreadHalf, pData->m_uiEvenOddTable, TRUE, uiVolumn);
//    _LAUNCH_KERNEL(_kernelBakeEvenOdd, _HC_DecompBlockHalf, _HC_DecompThreadHalf, pData->m_uiEvenOddTable, FALSE, uiVolumn);
//    checkCudaErrors(cudaDeviceSynchronize());
//    checkCudaErrors(cudaGetLastError());
//}

void CIndexSquare::BakeNaikTable(class CIndexData* pData, BYTE byFieldId)
{
    appParanoiac(_T("CIndexSquare::BakeNaikTable for field ID:%d\n"), byFieldId);
#if _CLG_MULTI_GPU
    //Improve-1 (multi-GPU-improve1.md I7): the Naik cache walks 3 links in one
    //step -- the widest stencil reach in the code. On a decomposed lattice every
    //such 3-hop must resolve inside the halo, so fail fast when the configured
    //width cannot cover it; otherwise the coordinate walk silently wraps to a
    //wrong local site (see I8: out-of-range must become Invalid, not wrap).
    if (NULL != appGetComm() && appGetComm()->Size() > 1 && _HC_HaloWidth < 3)
    {
        appCrucial(_T("CIndexSquare::BakeNaikTable: the Naik 3-hop path requires HaloWidth >= 3 on a decomposed lattice (configured HaloWidth=%d). Set HaloWidth: 3 (or wider) in the parameter file, or run unsplit.\n"),
            static_cast<INT>(_HC_HaloWidth));
        //Improve-1 (multi-GPU-improve1.md 3.8): an over-wide access must fail,
        //never continue into a silently wrong Naik table.
        _FAIL_EXIT;
    }
#endif
    checkCudaErrors(__cudaMalloc((void**)&pData->m_pNaikCache[byFieldId], sizeof(SIndex) * _HC_Volume * _HC_Dir * 2));
    preparethreadIndex;
    _LAUNCH_KERNEL(_kernelCacheNaikIndex, block, threads, pData->m_pNaikCache[byFieldId], pData->m_pIndexPositionToSIndex[byFieldId], pData->m_pSmallData);
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaGetLastError());
}

UINT CIndexSquare::GetPlaqutteCount(BYTE byFieldId) const
{
    UINT res[1] = { 0 };
    UINT* deviceRes;
    checkCudaErrors(__cudaMalloc((void**)&deviceRes, sizeof(UINT)));
    checkCudaErrors(cudaMemcpy(deviceRes, res, sizeof(UINT), cudaMemcpyHostToDevice));

    preparethreadIndex;
    _LAUNCH_KERNEL(_kernelPlaqutteCount, block, threads, deviceRes, byFieldId);

    checkCudaErrors(cudaMemcpy(res, deviceRes, sizeof(UINT), cudaMemcpyDeviceToHost));
    checkCudaErrors(__cudaFree(deviceRes));

    return res[0];
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================