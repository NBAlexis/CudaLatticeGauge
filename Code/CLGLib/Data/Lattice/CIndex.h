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

//For 4D cubic, this is 18
#define kMaxPlaqutteCache (32)

//EIT been used with Enum Integrator Type

__DEFINE_ENUM(EIndexType,

    EIndexType_Square,
    EIndexType_Max,
    EIndexType_ForceDWORD = 0x7fffffff,
    )


extern "C" { extern void _cCreateIndex(void** devicePtr, deviceBoundaryCondition ** pBC, UINT* size, EIndexType eIT); }

class CLGAPI CIndex
{
public:
    enum { kSpace = 0x01, kTime = 0x10, kSpaceTime = 0x11, };

    __device__ CIndex(class deviceBoundaryCondition * devicePtr) : m_pBoundaryCondition(devicePtr) { ; }
    __device__ ~CIndex()
    {
        appSafeDelete(m_pBoundaryCondition);
    }

    /**
    * This function to get the plaquttes indexes.
    * one can specify to get only spatial plaquttes or spacetime
    */
    __device__ virtual void _deviceGetPlaquttesAtLink(SIndex* retV, UINT& count, UINT& plaqutteLength, UINT uiLinkIndex, UINT st = kSpaceTime) const = 0;

    /**
    * Different from _deviceGetPlaquttesAtLink which return all plaquttes related to the link (For D dimension, square lattice, there are 2(D-1))
    * This function get all plaquttes related to the site (For D dimension, square lattice, there are D(D-1)/2 for each site, thus is (D-1)/2 per link which is 1/4 of above because each plaqutte has 4 edges)
    */
    __device__ virtual void _deviceGetPlaquttesAtSite(SIndex* retV, UINT& count, UINT& plaqutteLength, UINT uiSiteIndex, UINT st = kSpaceTime) const = 0;

    /**
    * Use for cache
    */
    __device__ virtual void _deviceGetPlaqutteCountLength(UINT& plaqLength, UINT& countPerSite, UINT& countPerLink) = 0;

#pragma region Index Walking

    __device__ virtual SIndex _deviceFermionIndexWalk(BYTE uiFieldId, UINT uiSiteIndex, INT uiWalkDir) const = 0;
    __device__ virtual SIndex _deviceGaugeIndexWalk(UINT uiSiteIndex, INT uiWalkDir) const = 0;

#pragma endregion

    class deviceBoundaryCondition * m_pBoundaryCondition;
};

__END_NAMESPACE

#endif //#ifndef _CINDEX_H_

//=============================================================================
// END OF FILE
//=============================================================================