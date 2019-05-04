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
        , m_pWalkingTable(NULL)
        , m_pPlaqutteCache(NULL)
        , m_pStappleCache(NULL)
        , m_byRegionTable(NULL)
    {
        checkCudaErrors(cudaMalloc((void**)&m_pSmallData, sizeof(UINT) * kCacheIndexSmallDataCount));
        checkCudaErrors(cudaMalloc((void**)&m_pWalkingTable, sizeof(UINT)
            * (_HC_Lx + 2 * kCacheIndexEdge) * (_HC_Ly + 2 * kCacheIndexEdge) 
            * (_HC_Lz + 2 * kCacheIndexEdge) * (_HC_Lt + 2 * kCacheIndexEdge)
            * _HC_Dir * 2));

        //region id is a byte, so max is 256
        checkCudaErrors(cudaMalloc((void**)&m_byRegionTable, sizeof(UINT) * 256));
        checkCudaErrors(cudaMalloc((void**)&m_pIndexPositionToSIndex, sizeof(SIndex*) * kMaxFieldCount));
        
        memset(m_pIndexPositionToSIndex, 0, sizeof(SIndex*) * kMaxFieldCount);
        memset(m_pGaugeMoveCache, 0, sizeof(SIndex*) * kMaxFieldCount);
        memset(m_pFermionMoveCache, 0, sizeof(SIndex*) * kMaxFieldCount);

    }

    ~CIndexData()
    {
        checkCudaErrors(cudaFree(m_pSmallData));
        checkCudaErrors(cudaFree(m_pWalkingTable));
        checkCudaErrors(cudaFree(m_byRegionTable));

        if (NULL != m_pPlaqutteCache)
        {
            checkCudaErrors(cudaFree(m_pPlaqutteCache));
        }

        if (NULL != m_pStappleCache)
        {
            checkCudaErrors(cudaFree(m_pStappleCache));
        }

        for (BYTE i = 0; i < kMaxFieldCount; ++i)
        {
            if (NULL != m_pIndexPositionToSIndex[i])
            {
                checkCudaErrors(cudaFree(m_pIndexPositionToSIndex[i]));
                m_pIndexPositionToSIndex[i] = NULL;
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
    }

    __device__ __inline__ UINT _deviceGetBigIndex(const SSmallInt4& inSite) const
    {
        return inSite.x * m_pSmallData[CIndexData::kMultX]
             + inSite.y * m_pSmallData[CIndexData::kMultY]
             + inSite.z * m_pSmallData[CIndexData::kMultZ]
             + inSite.w;
    } 

    __device__ __inline__ SIndex _deviceIndexWalk(
        BYTE byFieldId, 
        const SSmallInt4& inSite, 
        SBYTE uiWalkDir) const
    {
        //walking
        return m_pDeviceIndexPositionToSIndex[byFieldId]
            [m_pWalkingTable
                [
                _deviceGetBigIndex(inSite) * 2 * _DC_Dir 
                + (uiWalkDir > 0) ? (uiWalkDir - 1) : (_DC_Dir - uiWalkDir - 1)
                ]
            ];
    }

    __device__ void _deviceIndexWalkChain(
        BYTE byFieldId, 
        const SSmallInt4& inSite,
        const SBYTE* uiWalkDir, 
        SIndex* res, 
        BYTE byCount) const
    {
        //map to Index Position
        UINT idxPos = _deviceGetBigIndex(inSite);
        for (BYTE i = 0; i < byCount; ++i)
        {
            BYTE dirNow = uiWalkDir[i];
            //walking
            dirNow = (dirNow > 0) ? (dirNow - 1) : static_cast<BYTE>(_DC_Dir - dirNow - 1);
            idxPos = m_pWalkingTable[idxPos * 2 * _DC_Dir + dirNow];
            //return a SIndex
            res[i] = m_pDeviceIndexPositionToSIndex[byFieldId][idxPos];
        }
    }

    __device__ __inline__ UINT _devcieExchangeBoundaryFieldSiteIndex(const SIndex &site) const
    {
        return NULL == m_byRegionTable ? 0 : m_byRegionTable[site.m_byReginId];
    }

    __device__ void DebugPrintWalkingTable()
    {

    }

    //=============================================================
    //Small Data
    UINT* m_pSmallData;
    UINT* m_byRegionTable;

    //extend site * dir * 2
    //cached the neighbours of a site, cached as a index of m_pWalkingTable[index]
    UINT* m_pWalkingTable;

    //extend site position to SIndex mapping (i.e. m_pIndexPositionToSIndex[index])
    SIndex* m_pIndexPositionToSIndex[kMaxFieldCount];
    //used for device function
    SIndex** m_pDeviceIndexPositionToSIndex;

    //16*site
    SIndex* m_pPlaqutteCache;

    //24*links
    SIndex* m_pStappleCache;

    SIndex* m_pGaugeMoveCache[kMaxFieldCount];
    SIndex* m_pFermionMoveCache[kMaxFieldCount];

    BYTE m_uiPlaqutteLength;
    BYTE m_uiPlaqutteCountPerSite;
    BYTE m_uiPlaqutteCountPerLink;
};

__END_NAMESPACE

#endif //#ifndef _CINDEXDATA_H_

//=============================================================================
// END OF FILE
//=============================================================================