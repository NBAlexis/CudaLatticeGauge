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

class CLGAPI CBoundaryConditionTorusSquare : public CBoundaryCondition
{
public:

    CBoundaryConditionTorusSquare(class CIndex * pOwner) : CBoundaryCondition(pOwner) { ; }

    virtual __device__ uint2 GetMappedIndex(const int4 &site, const int4 &fromsite, const UINT* length, const UINT* mult);
};

__END_NAMESPACE

#endif //#ifndef _CBOUNDARYCONDITIONTORUSSQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================