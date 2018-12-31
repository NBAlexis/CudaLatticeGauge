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

class CLGAPI CBoundaryConditionTorusSquare : public deviceBoundaryCondition
{

public:

    __device__ CBoundaryConditionTorusSquare() : deviceBoundaryCondition() { ; }

    __device__ virtual SIndex _devcieGetMappedIndex(const int4 &site, const int4 &fromsite) const;

    __device__ virtual SIndex _devcieGetFermionMappedIndex(BYTE byFieldId, const int4 &site, const int4 &fromsite) const;
};

__END_NAMESPACE

#endif //#ifndef _CBOUNDARYCONDITIONTORUSSQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================