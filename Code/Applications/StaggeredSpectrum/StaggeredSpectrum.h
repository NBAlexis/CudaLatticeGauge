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
    ESSJ_SimulateQuenched,
    )

__DEFINE_ENUM(EStaggeredSpectrumMeasure,
    ESSM_Polyakov,
    ESSM_Wilson,
    ESSM_Correlator,
    ESSM_CorrelatorSimple,
    ESSM_Chiral,
    ESSM_All,
)

extern INT SimulateStaggered(CParameters& params);
extern INT Measurement(CParameters& params);
extern INT StaggeredGaugeFixing(CParameters& params);

//=============================================================================
// END OF FILE
//=============================================================================
