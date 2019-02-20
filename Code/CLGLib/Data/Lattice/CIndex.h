//=============================================================================
// FILENAME : CIndex.h
// 
// DESCRIPTION:
// This is the class for index on lattice
//
// Concepts:
//  site index: UINT, unique for every site, 
//      siteIndex = x * lengthY * lengthZ * lengthT + y * lengthZ * lengthT + z * lengthT + t
//  link index: UINT, unique for every link,
//      linkIndex = siteIndex * dir + dir
//  fat index: UINT, unique for both site and link
//      fatIndex = siteIndex * (dir + 1) + bSite ? 0 : (dir + 1)
//  data index: UINT, unique for every element
//      for gauge field, dataIndex = linkIndex * elementCount + n (for SU3, linkIndex * 9 + n, for SU2 linkIndex * 4 + n, etc)
//  mult array
//      mult array is use to calculate dataIndex
//      for gauge field:
//          dataIndex = x * mult[0] + y * mult[1] + z * mult[2] + t * mult[3] + dir * mult[4] + n
//          mult[0] = lengthY * lengthZ * lengthT * dir * elementCount
//          mult[1] = lengthZ * lengthT * dir * elementCount
//          mult[2] = lengthT * dir * elementCount
//          mult[3] = dir * elementCount
//          mult[4] = elementCount
//
// REVISION:
//  [12/5/2018 nbale]
//=============================================================================

#ifndef _CINDEX_H_
#define _CINDEX_H_

__BEGIN_NAMESPACE

//For 4D cubic, this is 24
#define kMaxPlaqutteCache (24)

//EIT been used with Enum Integrator Type

__DEFINE_ENUM(EIndexType,

    EIndexType_Square,
    EIndexType_Max,
    EIndexType_ForceDWORD = 0x7fffffff,
    )


extern "C" 
{ 
    extern void _cCreateIndex(void** devicePtr, deviceBoundaryCondition ** pBC, UINT* size, EIndexType eIT); 
}

class CLGAPI CIndex
{
public:
    enum { kSpace = 0x01, kTime = 0x02, kSpaceTime = 0x03, };

    __device__ CIndex(class deviceBoundaryCondition * devicePtr) : m_pBoundaryCondition(devicePtr) { ; }
    __device__ ~CIndex()
    {
        appSafeDelete(m_pBoundaryCondition);
    }

    /**
    * This function to get the plaquttes indexes.
    * one can specify to get only spatial plaquttes or spacetime
    */
    __device__ virtual void _deviceGetPlaquttesAtLink(SIndex* retV, BYTE& count, BYTE& plaqutteLength, UINT uiLinkIndex, BYTE st = kSpaceTime) const = 0;

    /**
    * Different from _deviceGetPlaquttesAtLink which return all plaquttes related to the link (For D dimension, square lattice, there are 2(D-1))
    * This function get all plaquttes related to the site (For D dimension, square lattice, there are D(D-1)/2 for each site, thus is (D-1)/2 per link which is 1/4 of above because each plaqutte has 4 edges)
    */
    __device__ virtual void _deviceGetPlaquttesAtSite(SIndex* retV, BYTE& count, BYTE& plaqutteLength, UINT uiSiteIndex, BYTE st = kSpaceTime) const = 0;

    /**
    * Use for cache
    */
    __device__ virtual void _deviceGetPlaqutteCountLength(BYTE& plaqLength, BYTE& countPerSite, BYTE& countPerLink) = 0;
    __device__ virtual void _deviceGetPlaquttesAtLinkAll(SIndex* retV, UINT uiLinkIndex) const = 0;
    __device__ virtual void _deviceGetPlaquttesAtSiteAll(SIndex* retV, UINT uiSiteIndex) const = 0;

#pragma region Index Walking

    __device__ virtual SIndex _deviceFermionIndexWalk(BYTE uiFieldId, UINT uiSiteIndex, SBYTE uiWalkDir) const = 0;
    __device__ virtual SIndex _deviceGaugeIndexWalk(UINT uiSiteIndex, SBYTE uiWalkDir) const = 0;

#pragma endregion

    class deviceBoundaryCondition * m_pBoundaryCondition;
};

class CIndexCache
{
public:
    CIndexCache() 
        : m_pPlaqutteCache(NULL)
        , m_pStappleCache(NULL)
        , m_uiPlaqutteLength(0)
        , m_uiPlaqutteCountPerSite(0)
        , m_uiPlaqutteCountPerLink(0)
    {
        for (UINT i = 0; i < kMaxFieldCount; ++i)
        {
            m_pGaugeMoveCache[i] = NULL;
            m_pFermionMoveCache[i] = NULL;
        }
    }
    ~CIndexCache()
    {
        if (NULL != m_pPlaqutteCache)
        {
            checkCudaErrors(__cudaFree(m_pPlaqutteCache));
        }
        if (NULL != m_pStappleCache)
        {
            checkCudaErrors(__cudaFree(m_pStappleCache));
        }

        for (UINT i = 0; i < kMaxFieldCount; ++i)
        {
            if (NULL != m_pGaugeMoveCache[i])
            {
                checkCudaErrors(__cudaFree(m_pGaugeMoveCache[i]));
            }
            if (NULL != m_pFermionMoveCache[i])
            {
                checkCudaErrors(__cudaFree(m_pFermionMoveCache[i]));
            }
        }
    }

    void CachePlaquttes();
    void CacheFermion(BYTE byFieldId);

    SIndex* m_pPlaqutteCache;
    SIndex* m_pStappleCache;
    SIndex* m_pGaugeMoveCache[kMaxFieldCount];
    SIndex* m_pFermionMoveCache[kMaxFieldCount];

    BYTE m_uiPlaqutteLength;
    BYTE m_uiPlaqutteCountPerSite;
    BYTE m_uiPlaqutteCountPerLink;
};

__END_NAMESPACE

#endif //#ifndef _CINDEX_H_

//=============================================================================
// END OF FILE
//=============================================================================