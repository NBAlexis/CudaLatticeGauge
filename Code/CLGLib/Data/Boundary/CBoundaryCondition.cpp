//=============================================================================
// FILENAME : CBoundaryCondition.cpp
// 
// DESCRIPTION:
// This is the class for boundary conditions
// Note that, the boundary conditions should only make sense together with lattice!!
//
// REVISION:
//  [07/23/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

void CBoundaryCondition::SetFieldSpecificBc(BYTE byFieldId, const SBoundCondition& bc)
{
    assert(byFieldId < kMaxFieldCount);
    m_FieldBC[byFieldId] = bc.m_sPeriodic;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================