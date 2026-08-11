//=============================================================================
// FILENAME : CIndexData.h
// 
// DESCRIPTION:
// To get rid of all the virtual device functions, we use a data-oriented index class
//
// Concepts:
// The walk of index is cached as data (index array), 
// So one do not care the types of index and boundary condition. (They are only used in cacheing)
// 
//
// REVISION:
//  [mm/dd/yy]
//  [04/22/2019 nbale]
//=============================================================================

#ifndef _CINDEXDATA_H_
#define _CINDEXDATA_H_

__BEGIN_NAMESPACE

class CLGAPI CIndexData
{
public:
    enum 
    { 
        kCacheIndexEdge = 4, 
        kCacheIndexSmallDataCount = 8, 

        kMultX = 0,
        kMultY = 1,
        kMultZ = 2,
        kPlaqLengthIdx = 3,
        kPlaqPerSiteIdx = 4,
        kPlaqPerLinkIdx = 5,
    };

    CIndexData()
        : m_pSmallData(NULL)
        , m_byRegionTable(NULL)
        , m_pMappingTable(NULL)
        , m_pSiteMappingTable(NULL)
        , m_pEtaMu(NULL)
        , m_pGlobalCoordinateTable(NULL)
        , m_pHaloGatherIndex(NULL)
        , m_uiHaloSiteCount(0)
        , m_uiSiteXYZT(1)
        , m_uiSiteXYZ(1)
        , m_uiLinkNumber(1)
        //, m_uiEvenOddTable(NULL)
    {
        checkCudaErrors(__cudaMalloc((void**)&m_pSmallData, sizeof(UINT) * kCacheIndexSmallDataCount));
        checkCudaErrors(__cudaMalloc((void**)&m_pMappingTable, sizeof(SSmallInt4)
            * (_HC_Lx + 2 * kCacheIndexEdge) * (_HC_Ly + 2 * kCacheIndexEdge)
            * (_HC_Lz + 2 * kCacheIndexEdge) * (_HC_Lt + 2 * kCacheIndexEdge) ));
        checkCudaErrors(__cudaMalloc((void**)&m_pSiteMappingTable, sizeof(SSmallInt4) * _HC_Lx * _HC_Ly * _HC_Lz * _HC_Lt));

        //checkCudaErrors(cudaMalloc((void**)&m_pBondInfoTable, sizeof(BYTE)
        //    * (_HC_Lx + 2 * kCacheIndexEdge) * (_HC_Ly + 2 * kCacheIndexEdge)
        //    * (_HC_Lz + 2 * kCacheIndexEdge) * (_HC_Lt + 2 * kCacheIndexEdge)
        //    * _HC_Dir));
        //memset(m_pBondInfoTable, 0, sizeof(BYTE*) * kMaxFieldCount);

        //region id is a byte, so max is 256
        checkCudaErrors(__cudaMalloc((void**)&m_byRegionTable, sizeof(UINT) * 256));

        checkCudaErrors(__cudaMalloc((void**)&m_pDeviceIndexPositionToSIndex, sizeof(SIndex*) * kMaxFieldCount));
        checkCudaErrors(__cudaMalloc((void**)&m_pDeviceIndexLinkToSIndex, sizeof(SIndex*) * kMaxFieldCount));
        memset(m_pIndexPositionToSIndex, 0, sizeof(SIndex*) * kMaxFieldCount);
        memset(m_pIndexLinkToSIndex, 0, sizeof(SIndex*) * kMaxFieldCount);

        memset(m_pPlaqutteCache, 0, sizeof(SIndex*) * kMaxFieldCount);
        memset(m_pStappleCache, 0, sizeof(SIndex*) * kMaxFieldCount);

        memset(m_pGaugeMoveCache, 0, sizeof(SIndex*) * kMaxFieldCount);
        memset(m_pMoveCache, 0, sizeof(SIndex*) * kMaxFieldCount);
        memset(m_pNaikCache, 0, sizeof(SIndex*) * kMaxFieldCount);

        memset(m_uiSiteNumber, 0, sizeof(UINT) * kMaxFieldCount);

        //checkCudaErrors(cudaMalloc((void**)&m_uiEvenOddTable, sizeof(UINT) * _HC_Lx * _HC_Ly * _HC_Lz * _HC_Lt));

        checkCudaErrors(cudaDeviceSynchronize());

    }

    ~CIndexData()
    {
        checkCudaErrors(__cudaFree(m_pSmallData));
        checkCudaErrors(__cudaFree(m_pMappingTable));
        checkCudaErrors(__cudaFree(m_pSiteMappingTable));
        checkCudaErrors(__cudaFree(m_byRegionTable));

        cudaSafeFree(m_pEtaMu);

        for (BYTE i = 0; i < kMaxFieldCount; ++i)
        {
            //if (NULL != m_pBondInfoTable[i])
            //{
            //    checkCudaErrors(cudaFree(m_pBondInfoTable[i]));
            //    m_pBondInfoTable[i] = NULL;
            //}

            if (NULL != m_pPlaqutteCache[i])
            {
                checkCudaErrors(__cudaFree(m_pPlaqutteCache[i]));
                m_pPlaqutteCache[i] = NULL;
            }

            if (NULL != m_pStappleCache[i])
            {
                checkCudaErrors(__cudaFree(m_pStappleCache[i]));
                m_pStappleCache[i] = NULL;
            }

            if (NULL != m_pIndexPositionToSIndex[i])
            {
                checkCudaErrors(__cudaFree(m_pIndexPositionToSIndex[i]));
                m_pIndexPositionToSIndex[i] = NULL;
            }
            if (NULL != m_pIndexLinkToSIndex[i])
            {
                checkCudaErrors(__cudaFree(m_pIndexLinkToSIndex[i]));
                m_pIndexLinkToSIndex[i] = NULL;
            }
            if (NULL != m_pGaugeMoveCache[i])
            {
                checkCudaErrors(__cudaFree(m_pGaugeMoveCache[i]));
                m_pGaugeMoveCache[i] = NULL;
            }
            if (NULL != m_pMoveCache[i])
            {
                checkCudaErrors(__cudaFree(m_pMoveCache[i]));
                m_pMoveCache[i] = NULL;
            }
            if (NULL != m_pNaikCache[i])
            {
                checkCudaErrors(__cudaFree(m_pNaikCache[i]));
                m_pNaikCache[i] = NULL;
            }
        }

        checkCudaErrors(__cudaFree(m_pDeviceIndexPositionToSIndex));
        checkCudaErrors(__cudaFree(m_pDeviceIndexLinkToSIndex));

        if (NULL != m_pGlobalCoordinateTable)
        {
            checkCudaErrors(__cudaFree(m_pGlobalCoordinateTable));
            m_pGlobalCoordinateTable = NULL;
        }
        if (NULL != m_pHaloGatherIndex)
        {
            checkCudaErrors(__cudaFree(m_pHaloGatherIndex));
            m_pHaloGatherIndex = NULL;
        }
        //checkCudaErrors(cudaFree(m_uiEvenOddTable));
    }

    __device__ __inline__ SSmallInt4 _deviceBigIndexToInt4(UINT uiBigIdx) const
    {
        return m_pMappingTable[uiBigIdx];
    }

    __device__ __inline__ UINT _deviceGetBigIndex(const SSmallInt4& inSite) const
    {
        return (inSite.x + CIndexData::kCacheIndexEdge) * m_pSmallData[CIndexData::kMultX]
             + (inSite.y + CIndexData::kCacheIndexEdge) * m_pSmallData[CIndexData::kMultY]
             + (inSite.z + CIndexData::kCacheIndexEdge) * m_pSmallData[CIndexData::kMultZ]
             + (inSite.w + CIndexData::kCacheIndexEdge);
    } 

    __device__ __inline__ SIndex _deviceGetMappingIndex(const SSmallInt4& inSite, BYTE byFieldId) const
    {
        return m_pDeviceIndexPositionToSIndex[byFieldId][_deviceGetBigIndex(inSite)];
    }

    __device__ __inline__ SIndex _deviceGetMappingLink(const SSmallInt4& inSite, BYTE byDir, BYTE byFieldId) const
    {
        return m_pDeviceIndexLinkToSIndex[byFieldId][_deviceGetBigIndex(inSite) * 4 + byDir];
    }

    __device__ __inline__ UBOOL _deviceIsBondOnSurface(UINT uiBigIdx, BYTE byFieldId, BYTE byDir) const
    {
        return (m_pDeviceIndexLinkToSIndex[byFieldId][uiBigIdx * 4 + byDir].m_byTag & _kDirichlet) != 0;
    }

    //====================================================
    // Directly using m_pDeviceIndexPositionToSIndex
    //====================================================
    //__device__ __inline__ SIndex _deviceIndexWalk(
    //    BYTE byFieldId, const SSmallInt4& inSite, SCHAR uiWalkDir) const
    //{
    //    //walking
    //    return m_pDeviceIndexPositionToSIndex[byFieldId]
    //        [_deviceIndexWalkBI(byFieldId, _deviceGetBigIndex(inSite), uiWalkDir)];
    //}


    //    return m_pDeviceIndexPositionToSIndex[byFieldId]
    //        [_deviceIndexWalkDoubleBI(
    //            byFieldId, 
    //            _deviceGetBigIndex(inSite), 
    //            uiWalkDir1, 
    //            uiWalkDir2
    //        )];
    //}

    __device__ __inline__ UINT _devcieExchangeBoundaryFieldSiteIndexBI(BYTE byField, UINT bigIdx) const
    {
        return NULL == m_byRegionTable ? 0 : m_byRegionTable[m_pDeviceIndexPositionToSIndex[byField][bigIdx].m_byReginId];
    }

    __device__ __inline__ UINT _devcieExchangeBoundaryFieldSiteIndex(const SIndex &site) const
    {
        return NULL == m_byRegionTable ? 0 : m_byRegionTable[site.m_byReginId];
    }

    static void DebugPlaqutteTable(const SSmallInt4& sSite, BYTE byFieldId);

    static void DebugPlaqutteTable(BYTE byFieldId);

    static void DebugEdgeMapping(BYTE byFieldId, const SSmallInt4& xyzt);

    static void DebugEdgeGlue(BYTE byFieldId, const SSmallInt4& xyzt);

    static void DebugStapleTable(BYTE byFieldId);
    static void DebugStapleTable(BYTE byFieldId, const SSmallInt4& xyzt);
    static void DebugStapleTable(BYTE byFieldId, const SSmallInt4& xyzt, UINT uiDir);
    static void DebugStapleTable(BYTE byFieldId, UINT uiIndex);

    static void DebugLinkDirichletOrDagger(BYTE byFieldId);

    //static void DebugPositionTable(BYTE byFieldId);
    //static void DebugLinkTable(BYTE byFieldId);

    //=============================================================
    //Small Data
    UINT* m_pSmallData;
    UINT* m_byRegionTable;

    //BigIndex mappingtable
    SSmallInt4* m_pMappingTable;
    //SiteIndex mapping table
    SSmallInt4* m_pSiteMappingTable;
    //BYTE* m_pBondInfoTable[kMaxFieldCount];

    //extend site position to SIndex mapping (i.e. m_pIndexPositionToSIndex[index])
    SIndex* m_pIndexPositionToSIndex[kMaxFieldCount];
    SIndex* m_pIndexLinkToSIndex[kMaxFieldCount];

    //used for device function
    //map big-index to sindex, the site outside lattice can have a big-index out of lattice,
    //however, after this mapping, it maps to a site inside lattice based on boundary condition
    SIndex** m_pDeviceIndexPositionToSIndex;
    //similar as m_pDeviceIndexPositionToSIndex, but with link
    SIndex** m_pDeviceIndexLinkToSIndex;

    //16*site
    SIndex* m_pPlaqutteCache[kMaxFieldCount];
    //18*links
    SIndex* m_pStappleCache[kMaxFieldCount];

    SIndex* m_pGaugeMoveCache[kMaxFieldCount];
    SIndex* m_pMoveCache[kMaxFieldCount];
    //SIndex* m_pBosonMoveCache[kMaxFieldCount];

    //eta mu table
    BYTE* m_pEtaMu;

    //Multi-GPU (Phase 1, Improve-1 I8/3.7): 32-bit global coordinate for every
    //slot of the field buffer, local slots [0, _DC_Volume) then halo slots
    //[_DC_Volume, _DC_Volume + m_uiHaloSiteCount). Baked in BakeAllIndexBuffer
    //for BOTH single- and multi-GPU (halo part empty on single-GPU), so
    //_deviceSIndexToGlobalInt4 is a plain table lookup and never narrows through
    //SCHAR. NULL until baked.
    SInt4* m_pGlobalCoordinateTable;

    //Multi-GPU (Phase 1, Improve-1 I8/3.7): halo gather map as SIndex. Index by
    //halo SITE slot (0..haloSites-1, i.e. relative to _DC_Volume); the entry's
    //m_uiSiteIndex is this rank's legal INTERIOR source site whose data feeds
    //that halo slot -- exactly the periodic wrap of the out-of-lattice neighbour.
    //CHaloManager packs from these sources and (for self-exchange) writes straight
    //into the field buffer's halo tail. NULL and 0-sized on single-GPU / unsplit.
    SIndex* m_pHaloGatherIndex;
    UINT m_uiHaloSiteCount;

    //[2*dir + 0] is the  3 site offset fermion
    //[2*dir + 1] is the -3 site offset fermion, and with link to the cached 3-link gauge field
    SIndex* m_pNaikCache[kMaxFieldCount];

#if !_CLG_ASSUME_SQUARE_LATTICE
    //For square lattice this is 4
    BYTE m_uiPlaqutteLength;
    //For 4D square lattice this is 6 ( C(d,2) = 6 )
    BYTE m_uiPlaqutteCountPerSite;
    //For 4D square lattice this is 6 ( 2(d-1) = 6 )
    BYTE m_uiPlaqutteCountPerLink;
#endif

    //Real size
    UINT m_uiSiteNumber[kMaxFieldCount];
    UINT m_uiSiteXYZT;
    UINT m_uiSiteXYZ;
    UINT m_uiLinkNumber;

    //this is slower
    //UINT* m_uiEvenOddTable;

};

#pragma region device functions

static __device__ __inline__ UINT _deviceGetBigIndex(const SSmallInt4& sSite, const UINT* __restrict__ pSmallData)
{
    return (sSite.x + CIndexData::kCacheIndexEdge) * pSmallData[CIndexData::kMultX]
        + (sSite.y + CIndexData::kCacheIndexEdge) * pSmallData[CIndexData::kMultY]
        + (sSite.z + CIndexData::kCacheIndexEdge) * pSmallData[CIndexData::kMultZ]
        + (sSite.w + CIndexData::kCacheIndexEdge);
}

static __device__ __inline__ SSmallInt4 _deviceBigIndexToInt4(UINT uiBigIdx)
{
    //SSmallInt4 coord;
    //coord.x = static_cast<SCHAR>(uiBigIdx / pSmallData[CIndexData::kMultX]) - CIndexData::kCacheIndexEdge;
    //coord.y = static_cast<SCHAR>((uiBigIdx % pSmallData[CIndexData::kMultX]) / pSmallData[CIndexData::kMultY]) - CIndexData::kCacheIndexEdge;
    //coord.z = static_cast<SCHAR>((uiBigIdx % pSmallData[CIndexData::kMultY]) / pSmallData[CIndexData::kMultZ]) - CIndexData::kCacheIndexEdge;
    //coord.w = static_cast<SCHAR>(uiBigIdx % pSmallData[CIndexData::kMultZ]) - CIndexData::kCacheIndexEdge;
    //return coord;
    return __idx->m_pMappingTable[uiBigIdx];
}

static __device__ __inline__ SSmallInt4 __deviceSiteIndexToInt4(UINT siteIndex)
{
    return __idx->m_pSiteMappingTable[siteIndex];
}

/**
 * Improve-1 (3.7): global coordinate of any SIndex -- local site or halo
 * slot -- as a plain lookup into the baked 32-bit table (local slots
 * [0, _DC_Volume), then halo slots). Never narrows through SCHAR, so global
 * extents beyond the SCHAR range are exact. An invalid SIndex has no legal
 * global coordinate and fails the assert in debug builds.
 */
static __device__ __inline__ SInt4 _deviceSIndexToGlobalInt4(const SIndex& sIndex)
{
    assert(!sIndex.IsInvalid());
    return __idx->m_pGlobalCoordinateTable[sIndex.m_uiSiteIndex];
}

#define __fwd(dir) static_cast<SCHAR>(dir + 1)
#define __bck(dir) static_cast<SCHAR>(-static_cast<SCHAR>(dir) - 1)
#define __bi(site) __idx->_deviceGetBigIndex(site)
#if _CLG_ASSUME_SQUARE_LATTICE
#define __bi4(site) (__idx->_deviceGetBigIndex(site) << 2U)
#else
#define __bi4(site) (__idx->_deviceGetBigIndex(site) * _DC_Dir)
#endif

///**
//* When both index and offset are known
//*/
//static __device__ __inline__ void _deviceSmallInt4Offset(SSmallInt4& sStart, INT idx, INT offset)
//{
//    sStart.m_byData4[idx] += offset;
//}

/**
* assume dir != 0, otherwise let it crush
*/
static __device__ __inline__ void _deviceSmallInt4Offset(SSmallInt4& sStart, SCHAR dir)
{
    const SCHAR sign_mask = (dir >> (sizeof(SCHAR) * 8 - 1)); //-1 (0xFFFFFFFF) of <0, otherwise 0
    // -1, -2, -3, -4, is 11111111, 11111110, 11111101, 11111100
    //so (-1) ^ (-1) = 0
    //   (-1) ^ (-2) = 1
    //   (-1) ^ (-3) = 2
    //   (-1) ^ (-4) = 3
    //so when dir < 0, -dir-1 = sign_mask ^ dir
    //when dir > 0, dir - 1   = sign_mask ^ (dir - 1)
    //so, dir - sign_mask - 1, if dir < 0, it is dir, if dir >= 0, it is dir - 1
    //idx = sign_mask ^ (dir - sign_mask - 1)
    sStart.m_byData4[sign_mask ^ (dir - sign_mask - 1)] += 1 - (sign_mask & 2);
}

/**
 * dir = 1,2,3,4 for +x,+y,+z,+t
 * dir = -1,-2,-3,-4 for -x,-y,-z,-t
 */
static __device__ __inline__ SSmallInt4 _deviceSmallInt4OffsetC(const SSmallInt4& sStart, SCHAR dir)
{
    SSmallInt4 ret = sStart;
    _deviceSmallInt4Offset(ret, dir);
    return ret;
}

static __device__ __inline__ void _deviceSmallInt4Offset(SSmallInt4& sStart, SCHAR* path, BYTE byLength)
{
    for (BYTE i = 0U; i < byLength; ++i)
    {
        _deviceSmallInt4Offset(sStart, path[i]);
    }
}

static __device__ __inline__ SSmallInt4 _deviceSmallInt4OffsetC(
    const SSmallInt4& sStart, const SCHAR* __restrict__ path, BYTE byLength)
{
    SSmallInt4 ret = sStart;
    for (BYTE i = 0U; i < byLength; ++i)
    {
        _deviceSmallInt4Offset(ret, path[i]);
    }
    return ret;
}

#pragma endregion

#pragma region Host Functions

inline static SSmallInt4 _hostBigIndexToInt4(UINT uiBigIdx)
{
    const UINT uiMX = (_HC_Ly + 2 * CIndexData::kCacheIndexEdge)
        * (_HC_Lz + 2 * CIndexData::kCacheIndexEdge)
        * (_HC_Lt + 2 * CIndexData::kCacheIndexEdge);
    const UINT uiMY = (_HC_Lz + 2 * CIndexData::kCacheIndexEdge)
        * (_HC_Lt + 2 * CIndexData::kCacheIndexEdge);
    const UINT uiMZ = _HC_Lt + 2 * CIndexData::kCacheIndexEdge;

    SSmallInt4 coord;
    coord.x = static_cast<SCHAR>(uiBigIdx / uiMX) - CIndexData::kCacheIndexEdge;
    coord.y = static_cast<SCHAR>((uiBigIdx % uiMX) / uiMY) - CIndexData::kCacheIndexEdge;
    coord.z = static_cast<SCHAR>((uiBigIdx % uiMY) / uiMZ) - CIndexData::kCacheIndexEdge;
    coord.w = static_cast<SCHAR>(uiBigIdx % uiMZ) - CIndexData::kCacheIndexEdge;
    return coord;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CINDEXDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================