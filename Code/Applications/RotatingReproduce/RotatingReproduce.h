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
    )


extern INT TestAngularMomentum(CParameters& params);
extern INT TestThermal(CParameters& params);



//=============================================================================
// END OF FILE
//=============================================================================
