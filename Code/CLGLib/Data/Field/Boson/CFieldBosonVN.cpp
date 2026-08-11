//=============================================================================
// FILENAME : CFieldBosonvn.cpp
// 
// DESCRIPTION:
// This is the class for the spin fields
//
// REVISION:
//  [mm/dd/yy]
//  [07/04/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "CFieldBosonVN.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldBosonU1)
__CLGIMPLEMENT_CLASS(CFieldBosonSU2)
__CLGIMPLEMENT_CLASS(CFieldBosonSU3)

#if _CLG_SU4_BOSON
__CLGIMPLEMENT_CLASS(CFieldBosonSU4)
#endif
#if _CLG_SU5_BOSON
__CLGIMPLEMENT_CLASS(CFieldBosonSU5)
#endif
#if _CLG_SU6_BOSON
__CLGIMPLEMENT_CLASS(CFieldBosonSU6)
#endif
#if _CLG_SU7_BOSON
__CLGIMPLEMENT_CLASS(CFieldBosonSU7)
#endif
#if _CLG_SU8_BOSON
__CLGIMPLEMENT_CLASS(CFieldBosonSU8)
#endif

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================