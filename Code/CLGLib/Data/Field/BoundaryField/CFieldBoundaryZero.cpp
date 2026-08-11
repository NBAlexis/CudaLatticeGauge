//=============================================================================
// FILENAME : CFieldBoundaryGaugeSUN.cu
// 
// DESCRIPTION:
// This is the class for index on square lattice
//
// REVISION:
//  [07/03/2024 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldBoundaryZero.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldBoundarySU3Vector)
__CLGIMPLEMENT_CLASS(CFieldBoundaryWilsonSquareSU3)

__CLGIMPLEMENT_CLASS(CFieldBoundaryBosonU1)
__CLGIMPLEMENT_CLASS(CFieldBoundaryBosonSU2)
__CLGIMPLEMENT_CLASS(CFieldBoundaryBosonSU3)

#if _CLG_SU4_BOSON
__CLGIMPLEMENT_CLASS(CFieldBoundaryBosonSU4)
#endif
#if _CLG_SU5_BOSON
__CLGIMPLEMENT_CLASS(CFieldBoundaryBosonSU5)
#endif
#if _CLG_SU6_BOSON
__CLGIMPLEMENT_CLASS(CFieldBoundaryBosonSU6)
#endif
#if _CLG_SU7_BOSON
__CLGIMPLEMENT_CLASS(CFieldBoundaryBosonSU7)
#endif
#if _CLG_SU8_BOSON
__CLGIMPLEMENT_CLASS(CFieldBoundaryBosonSU8)
#endif

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================