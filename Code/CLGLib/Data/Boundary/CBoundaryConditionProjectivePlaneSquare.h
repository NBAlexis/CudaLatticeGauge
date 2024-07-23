//=============================================================================
// FILENAME : CBoundaryConditionProjectivePlaneSquare.h
// 
// DESCRIPTION:
// This is the periodic boundary condition
// The X-Y Plane is projective plane
//
// REVISION:
//  [09/10/2020 nbale]
//=============================================================================

#ifndef _CBOUNDARYCONDITIONPROJECTIVEPLANESQUARE_H_
#define _CBOUNDARYCONDITIONPROJECTIVEPLANESQUARE_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CBoundaryConditionProjectivePlaneSquare)

class CLGAPI CBoundaryConditionProjectivePlaneSquare : public CBoundaryCondition
{
    __CLGDECLARE_CLASS(CBoundaryConditionProjectivePlaneSquare)

public:

    CBoundaryConditionProjectivePlaneSquare();

    void BakeEdgePoints(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceBuffer) const override;

    void BakeBondInfo(const SSmallInt4* deviceMappingTable, BYTE* deviceTable, BYTE byFieldId) const override;

    void BakeBondGlue(BYTE byFieldId, const SSmallInt4* deviceMappingTable, SIndex* deviceTable) const override;
};

__END_NAMESPACE

#endif //#ifndef _CBOUNDARYCONDITIONPROJECTIVEPLANESQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================