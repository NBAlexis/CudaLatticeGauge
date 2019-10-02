//=============================================================================
// FILENAME : RotatingReproduce.h
// 
// DESCRIPTION:
// Reproduce the results of PRL
//
// REVISION:
//  [02/26/2019 nbale]
//=============================================================================

#include "CLGLib.h"

__DEFINE_ENUM(ERotatingJob,
    ERJ_AngularMomentum,
    ERJ_Thermal,
    ERJ_PolyakovDist,
    ERJ_GaugeFixing,
    )


extern INT TestAngularMomentum(CParameters& params);
extern INT TestThermal(CParameters& params);
extern INT MeasurePolyakovDist(CParameters& params);
extern INT GaugeFixing(CParameters& params);




//=============================================================================
// END OF FILE
//=============================================================================
