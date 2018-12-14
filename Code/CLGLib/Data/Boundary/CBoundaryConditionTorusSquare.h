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
    __CLGDECLARE_CLASS(CBoundaryConditionTorusSquare)

public:

    CBoundaryConditionTorusSquare() : CBoundaryCondition() { ; }

    __device__ virtual uint2 _devcieGetMappedIndex(const int4 &site, const int4 &fromsite) const;
};

__END_NAMESPACE

#endif //#ifndef _CBOUNDARYCONDITIONTORUSSQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================