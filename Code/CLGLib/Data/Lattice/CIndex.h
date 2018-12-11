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

#pragma region index functions

//========================================
// implement after CLatticeData is known

__device__ __inline__ static UINT _deviceGetSiteIndex(const class CDeviceLattice * pLattice, const UINT* coord);
__device__ __inline__ static UINT _deviceGetLinkIndex(const class CDeviceLattice * pLattice, UINT siteIndex, UINT dir);
__device__ __inline__ static UINT _deviceGetLinkIndex(const class CDeviceLattice * pLattice, const UINT* coord, UINT dir);

/**
* for site, dir_plus_one = 0
* for link, dir_plus_one = dir + 1
*/
__device__ __inline__ static UINT _deviceGetFatIndex(const class CDeviceLattice * pLattice, const UINT* coord, UINT dir_plus_one);

/**
* int4.xyzw = x, y, z, t
*/
__device__ static __inline__ int4 __deviceSiteIndexToInt4(const class CDeviceLattice * pLattice, UINT siteIndex);
__device__ static __inline__ int4 __deviceLinkIndexToInt4(const class CDeviceLattice * pLattice, UINT linkIndex);
__device__ static __inline__ int4 __deviceFatIndexToInt4(const class CDeviceLattice * pLattice, UINT linkIndex);

#pragma endregion index functions

class CLGAPI CIndex
{

public:
    enum { kSpace = 0x01, kTime = 0x10, kSpaceTime = 0x11, };

    __device__ CIndex(class CDeviceLattice * pOwner) : m_pOwner(pOwner), m_pBoundaryCondition(NULL) { ; }

    __device__ class CDeviceLattice * GetOwner() const { return m_pOwner; }

    __device__ class CBoundaryCondition * GetBoundaryCondition() const { return m_pBoundaryCondition; }
    __device__ void SetBoundaryCondition(class CBoundaryCondition * pBoundaryCondition) { m_pBoundaryCondition = pBoundaryCondition; }

    /**
    * WOW, we can use virtual device functions!
    * This function to get the plaquttes indexes.
    * one can specify to get only spatial plaquttes or spacetime
    */
    __device__ virtual int2* _deviceGetPlaquttesAtLink(UINT& count, UINT& plaqutteLength, UINT uiLinkIndex, UINT st = kSpaceTime) const = 0;

protected:

    class CDeviceLattice * m_pOwner;
    class CBoundaryCondition * m_pBoundaryCondition;
};

__END_NAMESPACE

#endif //#ifndef _CINDEX_H_

//=============================================================================
// END OF FILE
//=============================================================================