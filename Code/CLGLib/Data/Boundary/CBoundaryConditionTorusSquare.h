//=============================================================================
// FILENAME : CBoundaryConditionTorusSquare.h
// 
// DESCRIPTION:
// This is the periodic boundary condition
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CBOUNDARYCONDITIONTORUSSQUARE_H_
#define _CBOUNDARYCONDITIONTORUSSQUARE_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CBoundaryConditionTorusSquare)

class CLGAPI CBoundaryConditionTorusSquare : public CBoundaryCondition
{
    __CLGDECLARE_CLASS(CBoundaryConditionTorusSquare)

public:

    CBoundaryConditionTorusSquare();

    void BakeEdgePoints(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceBuffer) const override;

    void BakeBondInfo(const SSmallInt4* deviceMappingTable, BYTE* deviceTable, BYTE byFieldId) const override;

    void BakeBondGlue(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceTable) const override;
};

__END_NAMESPACE

#endif //#ifndef _CBOUNDARYCONDITIONTORUSSQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================