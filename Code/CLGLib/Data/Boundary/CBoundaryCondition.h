//=============================================================================
// FILENAME : CBoundaryCondition.h
// 
// DESCRIPTION:
// This is the class for boundary conditions
// Note that, the boundary conditions should only make sense together with lattice!!
//
// REVISION:
//  [12/5/2018 nbale]
//=============================================================================

#ifndef _CBOUNDARYCONDITION_H_
#define _CBOUNDARYCONDITION_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(EBoundaryCondition,
    EBC_TorusSquare,
    EBC_TorusAndDirichlet,
    EBC_Max,
    EBC_ForceDWORD = 0x7fffffff,
    )


struct SBoundCondition
{
    SBoundCondition() {}

    union
    {
        SSmallInt4 m_sPeriodic;
        INT m_iPeriodic;
    };

    void* m_pBCFieldDevicePtr[8];
};


class CLGAPI CBoundaryCondition : public CBase
{
public:

    CBoundaryCondition()
    {
        
    }

    virtual void BakeEdgePoints(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceBuffer) const = 0;

    /**
    * For example, set field Id of BC or set anti-periodic condition
    */
    virtual void SetFieldSpecificBc(BYTE byFieldId, const SBoundCondition& bc) = 0;

    virtual void BakeRegionTable(UINT* deviceTable) const {}

    /**
    * The bondary condition of the bonds(links) sometimes should be different from the sites.
    */
    virtual void BakeBondInfo(const SSmallInt4* deviceMappingTable, BYTE* deviceTable) const = 0;

    virtual UBOOL NeedToFixBoundary() const { return FALSE; }

protected:

    /**
    * 0 For Dirichlet, 1 and -1 for periodic and anti-periodic
    * Note that the gauge field id is 1
    */
    SSmallInt4 m_FieldBC[kMaxFieldCount];

};

__END_NAMESPACE

#endif //#ifndef _CBOUNDARYCONDITION_H_

//=============================================================================
// END OF FILE
//=============================================================================