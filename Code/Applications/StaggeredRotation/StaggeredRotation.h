//=============================================================================
// FILENAME : StaggeredRotation.h
// 
// DESCRIPTION:
//   This is to compress the configuration files to half the size
//
// REVISION:
//  [09/24/2020 nbale]
//=============================================================================

#include "CLGLib.h"

__DEFINE_ENUM(EStaggeredRotationJob,
    ESR_Simulate,
    ESR_Measure,
    ESR_GaugeFixing,
    )

extern INT SimulateStaggeredRotation(CParameters& params);

//=============================================================================
// END OF FILE
//=============================================================================
