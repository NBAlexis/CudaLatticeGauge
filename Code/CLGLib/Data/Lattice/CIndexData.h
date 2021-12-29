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
        kCacheIndexEdge = 2, 
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
        , m_pBondInfoTable(NULL)
        , m_pPlaqutteCache(NULL)
        , m_pStappleCache(NULL)
        , m_pEtaMu(NULL)
        , m_uiSiteXYZT(1)
        , m_uiSiteXYZ(1)
        , m_uiLinkNumber(1)
    {
        checkCudaErrors(cudaMalloc((void**)&m_pSmallData, sizeof(UINT) * kCacheIndexSmallDataCount));
        checkCudaErrors(cudaMalloc((void**)&m_pMappingTable, sizeof(SSmallInt4)
            * (_HC_Lx + 2 * kCacheIndexEdge) * (_HC_Ly + 2 * kCacheIndexEdge)
            * (_HC_Lz + 2 * kCacheIndexEdge) * (_HC_Lt + 2 * kCacheIndexEdge) ));

        checkCudaErrors(cudaMalloc((void**)&m_pBondInfoTable, sizeof(BYTE)
            * (_HC_Lx + 2 * kCacheIndexEdge) * (_HC_Ly + 2 * kCacheIndexEdge)
            * (_HC_Lz + 2 * kCacheIndexEdge) * (_HC_Lt + 2 * kCacheIndexEdge)
            * _HC_Dir));

        //region id is a byte, so max is 256
        checkCudaErrors(cudaMalloc((void**)&m_byRegionTable, sizeof(UINT) * 256));

        checkCudaErrors(cudaMalloc((void**)&m_pDeviceIndexPositionToSIndex, sizeof(SIndex*) * kMaxFieldCount));
        checkCudaErrors(cudaMalloc((void**)&m_pDeviceIndexLinkToSIndex, sizeof(SIndex*) * kMaxFieldCount));
        memset(m_pIndexPositionToSIndex, 0, sizeof(SIndex*) * kMaxFieldCount);
        memset(m_pIndexLinkToSIndex, 0, sizeof(SIndex*) * kMaxFieldCount);

        memset(m_pGaugeMoveCache, 0, sizeof(SIndex*) * kMaxFieldCount);
        memset(m_pFermionMoveCache, 0, sizeof(SIndex*) * kMaxFieldCount);

        memset(m_uiSiteNumber, 0, sizeof(UINT) * kMaxFieldCount);

    }

    ~CIndexData()
    {
        checkCudaErrors(cudaFree(m_pSmallData));
        checkCudaErrors(cudaFree(m_pMappingTable));
        checkCudaErrors(cudaFree(m_byRegionTable));
        checkCudaErrors(cudaFree(m_pBondInfoTable));

        if (NULL != m_pPlaqutteCache)
        {
            checkCudaErrors(cudaFree(m_pPlaqutteCache));
        }

        if (NULL != m_pStappleCache)
        {
            checkCudaErrors(cudaFree(m_pStappleCache));
        }
        cudaSafeFree(m_pEtaMu);

        for (BYTE i = 0; i < kMaxFieldCount; ++i)
        {
            if (NULL != m_pIndexPositionToSIndex[i])
            {
                checkCudaErrors(cudaFree(m_pIndexPositionToSIndex[i]));
                m_pIndexPositionToSIndex[i] = NULL;
            }
            if (NULL != m_pIndexPositionToSIndex[i])
            {
                checkCudaErrors(cudaFree(m_pIndexLinkToSIndex[i]));
                m_pIndexLinkToSIndex[i] = NULL;
            }
            if (NULL != m_pGaugeMoveCache[i])
            {
                checkCudaErrors(cudaFree(m_pGaugeMoveCache[i]));
                m_pGaugeMoveCache[i] = NULL;
            }
            if (NULL != m_pFermionMoveCache[i])
            {
                checkCudaErrors(cudaFree(m_pFermionMoveCache[i]));
                m_pFermionMoveCache[i] = NULL;
            }
        }

        checkCudaErrors(cudaFree(m_pDeviceIndexPositionToSIndex));
        checkCudaErrors(cudaFree(m_pDeviceIndexLinkToSIndex));
    }

    __device__ __inline__ SSmallInt4 _deviceBigIndexToInt4(UINT uiBigIdx) const
    {
        SSmallInt4 coord;
        coord.x = static_cast<SBYTE>(uiBigIdx / m_pSmallData[kMultX]) - CIndexData::kCacheIndexEdge;
        coord.y = static_cast<SBYTE>((uiBigIdx % m_pSmallData[kMultX]) / m_pSmallData[kMultY]) - CIndexData::kCacheIndexEdge;
        coord.z = static_cast<SBYTE>((uiBigIdx % m_pSmallData[kMultY]) / m_pSmallData[kMultZ]) - CIndexData::kCacheIndexEdge;
        coord.w = static_cast<SBYTE>(uiBigIdx % m_pSmallData[kMultZ]) - CIndexData::kCacheIndexEdge;
        return coord;
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

    __device__ __inline__ UBOOL _deviceIsBondOnSurface(UINT uiBigIdx, BYTE byDir) const
    {
        return (m_pBondInfoTable[uiBigIdx * _DC_Dir + byDir] & _kDirichlet) != 0;
    }

    //====================================================
    // Directly using m_pDeviceIndexPositionToSIndex
    //====================================================
    //__device__ __inline__ SIndex _deviceIndexWalk(
    //    BYTE byFieldId, const SSmallInt4& inSite, SBYTE uiWalkDir) const
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

    static void DebugPlaqutteTable(const SSmallInt4& sSite);

    static void DebugPlaqutteTable();

    static void DebugEdgeMapping(BYTE byFieldId, const SSmallInt4& xyzt);

    static void DebugEdgeGlue(BYTE byFieldId, const SSmallInt4& xyzt);

    //=============================================================
    //Small Data
    UINT* m_pSmallData;
    UINT* m_byRegionTable;

    SSmallInt4* m_pMappingTable;
    BYTE* m_pBondInfoTable;

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
    SIndex* m_pPlaqutteCache;

    //24*links
    SIndex* m_pStappleCache;

    SIndex* m_pGaugeMoveCache[kMaxFieldCount];
    SIndex* m_pFermionMoveCache[kMaxFieldCount];

    //eta mu table
    BYTE* m_pEtaMu;

    BYTE m_uiPlaqutteLength;
    BYTE m_uiPlaqutteCountPerSite;
    BYTE m_uiPlaqutteCountPerLink;

    //Real size
    UINT m_uiSiteNumber[kMaxFieldCount];
    UINT m_uiSiteXYZT;
    UINT m_uiSiteXYZ;
    UINT m_uiLinkNumber;

};

#pragma region device functions

static __device__ __inline__ UINT _deviceGetBigIndex(const SSmallInt4& sSite, const UINT* __restrict__ pSmallData)
{
    return (sSite.x + CIndexData::kCacheIndexEdge) * pSmallData[CIndexData::kMultX]
        + (sSite.y + CIndexData::kCacheIndexEdge) * pSmallData[CIndexData::kMultY]
        + (sSite.z + CIndexData::kCacheIndexEdge) * pSmallData[CIndexData::kMultZ]
        + (sSite.w + CIndexData::kCacheIndexEdge);
}

static __device__ __inline__ SSmallInt4 _deviceBigIndexToInt4(UINT uiBigIdx, const UINT* __restrict__ pSmallData)
{
    SSmallInt4 coord;
    coord.x = static_cast<SBYTE>(uiBigIdx / pSmallData[CIndexData::kMultX]) - CIndexData::kCacheIndexEdge;
    coord.y = static_cast<SBYTE>((uiBigIdx % pSmallData[CIndexData::kMultX]) / pSmallData[CIndexData::kMultY]) - CIndexData::kCacheIndexEdge;
    coord.z = static_cast<SBYTE>((uiBigIdx % pSmallData[CIndexData::kMultY]) / pSmallData[CIndexData::kMultZ]) - CIndexData::kCacheIndexEdge;
    coord.w = static_cast<SBYTE>(uiBigIdx % pSmallData[CIndexData::kMultZ]) - CIndexData::kCacheIndexEdge;
    return coord;
}

#define __fwd(dir) (dir + 1)
#define __bck(dir) (-static_cast<INT>(dir) - 1)
#define __bi(site) __idx->_deviceGetBigIndex(site)
#define __bi4(site) __idx->_deviceGetBigIndex(site) * _DC_Dir

/**
 * dir = 1,2,3,4 for +x,+y,+z,+t
 * dir = -1,-2,-3,-4 for -x,-y,-z,-t
 */
static __device__ __inline__ SSmallInt4 _deviceSmallInt4OffsetC(
    const SSmallInt4& sStart, INT dir)
{
    SSmallInt4 ret = sStart;
    if (0 == dir)
    {
        return ret;
    }
    const INT idx = dir < 0 ? (-dir - 1) : (dir - 1);
    ret.m_byData4[idx] = ret.m_byData4[idx] + (dir > 0 ? 1 : (-1));
    return ret;
}

static __device__ __inline__ void _deviceSmallInt4Offset(SSmallInt4& sStart, INT dir)
{
    if (0 != dir)
    {
        const INT idx = dir < 0 ? (-dir - 1) : (dir - 1);
        sStart.m_byData4[idx] = sStart.m_byData4[idx] + (dir > 0 ? 1 : (-1));
    }
}

static __device__ __inline__ void _deviceSmallInt4Offset(
    SSmallInt4& sStart, INT* path, BYTE byLength)
{
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == path[i])
        {
            continue;
        }
        const INT idx = path[i] < 0 ? (-path[i] - 1) : (path[i] - 1);
        sStart.m_byData4[idx] = sStart.m_byData4[idx] + (path[i] > 0 ? 1 : (-1));
    }
}

static __device__ __inline__ SSmallInt4 _deviceSmallInt4OffsetC(
    const SSmallInt4& sStart, const INT* __restrict__ path, BYTE byLength)
{
    SSmallInt4 ret = sStart;
    for (BYTE i = 0; i < byLength; ++i)
    {
        if (0 == path[i])
        {
            continue;
        }
        const INT idx = path[i] < 0 ? (-path[i] - 1) : (path[i] - 1);
        ret.m_byData4[idx] = ret.m_byData4[idx] + (path[i] > 0 ? 1 : (-1));
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
    coord.x = static_cast<SBYTE>(uiBigIdx / uiMX) - CIndexData::kCacheIndexEdge;
    coord.y = static_cast<SBYTE>((uiBigIdx % uiMX) / uiMY) - CIndexData::kCacheIndexEdge;
    coord.z = static_cast<SBYTE>((uiBigIdx % uiMY) / uiMZ) - CIndexData::kCacheIndexEdge;
    coord.w = static_cast<SBYTE>(uiBigIdx % uiMZ) - CIndexData::kCacheIndexEdge;
    return coord;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CINDEXDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================