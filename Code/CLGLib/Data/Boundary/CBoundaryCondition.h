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
    virtual void SetFieldSpecificBc(BYTE byFieldId, const SBoundCondition& bc);

    virtual void BakeRegionTable(UINT* deviceTable) const {}

    /**
    * The bondary condition of the bonds(links) sometimes should be different from the sites.
    * 
    * Discarded: bond info are all put into link table
    */
    //virtual void BakeBondInfo(const SSmallInt4* deviceMappingTable, BYTE* deviceTable, BYTE byFieldId) const = 0;

    /**
     * For the glue of links, we must bake the target.
     * For example, 0,0,0,5 _1 can be map to 0,1,0,0 _2^+
     */
    virtual void BakeBondGlue(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceTable) const = 0;

    virtual UBOOL NeedToFixBoundary() const { return FALSE; }

    SSmallInt4 GetFieldBC(BYTE byFieldId) const { return m_FieldBC[byFieldId]; }

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