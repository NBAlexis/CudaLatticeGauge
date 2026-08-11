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
#include "CFieldBoundaryOne.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeU1)
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU2)
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU3)

#if _CLG_SU4_GAUGE
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU4)
#endif
#if _CLG_SU5_GAUGE
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU5)
#endif
#if _CLG_SU6_GAUGE
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU6)
#endif
#if _CLG_SU7_GAUGE
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU7)
#endif
#if _CLG_SU8_GAUGE
__CLGIMPLEMENT_CLASS(CFieldBoundaryGaugeSU8)
#endif

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================