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

__DEFINE_ENUM(EIndexType,

    EIndexType_Square,
    EIndexType_Max,
    EIndexType_ForceDWORD = 0x7fffffff,
    )



__DEFINE_ENUM(EIndexEvenOdd,

    EIE_All,
    EIE_EvenToOdd,
    EIE_OddToEven,
    EIE_EvenToOddOnEven,
    EIE_OddToEvenOnOdd,
    EIE_ForceDWORD = 0x7fffffff,
    )


class CLGAPI CIndex : public CBase
{
public:

    CIndex() : m_pBoundaryCondition(NULL), m_bNeedToFixBoundary(FALSE) {  }
    ~CIndex()
    {
        appSafeDelete(m_pBoundaryCondition);
    }

    void SetBoundaryCondition(class CBoundaryCondition * pBc) 
    { 
        m_pBoundaryCondition = pBc; 
        m_bNeedToFixBoundary = pBc->NeedToFixBoundary();
    }

    /**
    * Now the functional of the Index is data-based.
    * All index walking thing is cached and BakeMappingTable
    */
    virtual void BakeAllIndexBuffer(class CIndexData* pData) = 0;

    /**
    * Cache the plaquttes index so that we do not have to walk every time
    */
    virtual void BakePlaquttes(class CIndexData* pData, BYTE byFieldId) = 0;

    /**
    * Cache the neighbour index used for fermion fiels
    */
    virtual void BakeMoveIndex(class CIndexData* pData, BYTE byFieldId) = 0;

    /**
     * pre-calculate eta mu
     */
    virtual void BakeEtaMuTable(class CIndexData* pData) = 0;

    /**
    * Call this after BakePlaquttes
    */
    virtual UINT GetPlaqutteCount(BYTE byFieldId) const = 0;

    /**
    * Because of Dirchlet boundary condition, 
    * when calculate average, we need the real volume
    */
    virtual void CalculateSiteCount(class CIndexData* pData) const;

    UBOOL NeedToFixBoundary() const { return m_bNeedToFixBoundary; }

    class CBoundaryCondition* GetBoudanryCondition() const { return m_pBoundaryCondition; }

    class CBoundaryCondition * m_pBoundaryCondition;
    UBOOL m_bNeedToFixBoundary;
};


__END_NAMESPACE

#endif //#ifndef _CINDEX_H_

//=============================================================================
// END OF FILE
//=============================================================================