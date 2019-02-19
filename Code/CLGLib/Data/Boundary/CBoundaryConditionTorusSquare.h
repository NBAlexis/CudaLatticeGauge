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

    __device__ CBoundaryConditionTorusSquare();

    __device__ virtual SIndex _devcieGetMappedIndex(const SSmallInt4 &site, const SIndex &fromsite) const;

    __device__ virtual SIndex _devcieGetFermionMappedIndex(BYTE byFieldId, const SSmallInt4 &site, const SIndex &fromsite) const;

    __device__ virtual void SetFieldSpecificBc(BYTE byFieldId, const SBoundCondition& bc);

    SSmallInt4 m_FermionBC[_kMaxFieldCount];
};

__END_NAMESPACE

#endif //#ifndef _CBOUNDARYCONDITIONTORUSSQUARE_H_

//=============================================================================
// END OF FILE
//=============================================================================