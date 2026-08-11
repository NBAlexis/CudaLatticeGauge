//=============================================================================
// FILENAME : CFieldBosonVNRotation.cpp
// 
// DESCRIPTION:
// This is the class for the spin fields
//
// REVISION:
//  [07/13/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "CFieldBosonVNRotation.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationU1)
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU2)
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU3)

#if _CLG_SU4_BOSON
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU4)
#endif
#if _CLG_SU5_BOSON
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU5)
#endif
#if _CLG_SU6_BOSON
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU6)
#endif
#if _CLG_SU7_BOSON
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU7)
#endif
#if _CLG_SU8_BOSON
__CLGIMPLEMENT_CLASS(CFieldBosonVNRotationSU8)
#endif

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================