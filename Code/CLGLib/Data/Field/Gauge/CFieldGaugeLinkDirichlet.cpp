//=============================================================================
// FILENAME : CFieldGaugeLinkDirichlet.cu
// 
// DESCRIPTION:
//
// REVISION:
//  [07/06/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldGaugeLinkDirichlet.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldGaugeU1D)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU2D)
__CLGIMPLEMENT_CLASS(CFieldGaugeSU3D)

#if _CLG_SU4_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeSU4D)
#endif
#if _CLG_SU5_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeSU5D)
#endif
#if _CLG_SU6_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeSU6D)
#endif
#if _CLG_SU7_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeSU7D)
#endif
#if _CLG_SU8_GAUGE
__CLGIMPLEMENT_CLASS(CFieldGaugeSU8D)
#endif

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================