//=============================================================================
// FILENAME : CBoundaryConditionPeriodicAndDirichletSquare.h
// 
// DESCRIPTION:
// This is the periodic boundary and Dirichlet condition
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

    virtual void BakeEdgePoints(BYTE byFieldId, SIndex* deviceBuffer) const;

    virtual void SetFieldSpecificBc(BYTE byFieldId, const SBoundCondition& bc);

    virtual void BakeRegionTable(UINT* deviceTable) const;

};

__END_NAMESPACE

#endif //#ifndef _CBOUNDARYCONDITIONPERIODICANDDIRICHLETSQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================