//=============================================================================
// FILENAME : ConstAcc.h
// 
// DESCRIPTION:
//   This is to compress the configuration files to half the size
//
// REVISION:
//  [04/10/2020 nbale]
//=============================================================================

#include "CLGLib.h"

__DEFINE_ENUM(EAccJob,
    EAJ_Simulate,
    EAJ_SimulateQ,
    EAJ_Measure,
    )

extern INT SimulateAcc(CParameters& params);
extern INT Measurement(CParameters& params);

//=============================================================================
// END OF FILE
//=============================================================================
