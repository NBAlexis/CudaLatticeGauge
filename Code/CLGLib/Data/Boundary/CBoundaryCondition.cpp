//=============================================================================
// FILENAME : CBoundaryCondition.cpp
// 
// DESCRIPTION:
// This is the class for boundary conditions
// Note that, the boundary conditions should only make sense together with lattice!!
//
// REVISION:
//  [12/7/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__device__ 
CBoundaryCondition::CBoundaryCondition(class CIndex * pOwner)
    : m_pOwner(pOwner)
    , m_pLattice(pOwner->GetOwner())
{
    ;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================