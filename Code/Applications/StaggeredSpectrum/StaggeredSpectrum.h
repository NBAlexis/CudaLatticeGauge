//=============================================================================
// FILENAME : StaggeredSpectrum.h
// 
// DESCRIPTION:
//   This is to compress the configuration files to half the size
//
// REVISION:
//  [08/21/2020 nbale]
//=============================================================================

#include "CLGLib.h"

__DEFINE_ENUM(EStaggeredSpectrumJob,
    ESSJ_Simulate,
    ESSJ_GaugeFixing,
    ESSJ_Measure,
    )

__DEFINE_ENUM(EStaggeredSpectrumMeasure,
    ESSM_Mass,
    ESSM_Vr,
    )

extern INT SimulateStaggered(CParameters& params);

//=============================================================================
// END OF FILE
//=============================================================================
