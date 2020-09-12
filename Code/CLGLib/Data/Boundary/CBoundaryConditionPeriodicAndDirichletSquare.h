//=============================================================================
// FILENAME : CBoundaryConditionPeriodicAndDirichletSquare.h
// 
// DESCRIPTION:
// This is the periodic boundary and Dirichlet condition
// For Dirichlet boundary, we set the boundary to fixed values before evaluation
// That is because, if we only get element values from fixed values, the Hamitonian will changed (not adiabat)
// Therefor, the Metropolis step will fail.
//
// REVISION:
//  [04/20/2019 nbale]
//=============================================================================

#ifndef _CBOUNDARYCONDITIONPERIODICANDDIRICHLETSQUARE_H_
#define _CBOUNDARYCONDITIONPERIODICANDDIRICHLETSQUARE_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CBoundaryConditionPeriodicAndDirichletSquare)

class CLGAPI CBoundaryConditionPeriodicAndDirichletSquare : public CBoundaryCondition
{
    __CLGDECLARE_CLASS(CBoundaryConditionPeriodicAndDirichletSquare)

public:


    enum
    {
        byXLeft = (1 << 0),
        byYLeft = (1 << 1),
        byZLeft = (1 << 2),
        byTLeft = (1 << 3),
        byXRight = (1 << 4),
        byYRight = (1 << 5),
        byZRight = (1 << 6),
        byTRight = (1 << 7),
    };

    CBoundaryConditionPeriodicAndDirichletSquare();

    void BakeEdgePoints(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceBuffer) const override;

    void SetFieldSpecificBc(BYTE byFieldId, const SBoundCondition& bc) override;

    void BakeRegionTable(UINT* deviceTable) const override;

    //It is only neccessary when simulating with holes or inpuries
    UBOOL NeedToFixBoundary() const override { return TRUE; }

    void BakeBondInfo(const SSmallInt4* deviceMappingTable, BYTE* deviceTable) const override;

    void BakeBondGlue(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceTable) const override;
};

__END_NAMESPACE

#endif //#ifndef _CBOUNDARYCONDITIONPERIODICANDDIRICHLETSQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================