//=============================================================================
// FILENAME : PSU3WithBoundary.h
//
// DESCRIPTION:
// Application driver for the PSU(3) fundamental-lift action with a Z3 2-form
// boundary field B (CActionGaugePlaquettePSU3WithBoundary + CFieldTensor2Z3).
//
//   S(U,B) = -(Beta/3) * sum_p Re[B_p Tr(U_p)]
//
// Two work jobs:
//   EPSU3J_Simulate - HMC simulation; saves gauge configurations together with
//                     the dynamic tensor2 (B) companion files
//   EPSU3J_Measure  - loads saved configurations (gauge + B companion) and
//                     measures Polyakov loop, action energy and the Z3
//                     monopole count of B
//
// REVISION:
//  [08/11/26]
//=============================================================================
#pragma once

#include "CLGLib.h"

__DEFINE_ENUM(EPSU3Job,
    EPSU3J_Simulate,
    EPSU3J_Measure,
    )

extern INT Simulate(CParameters& params);
extern INT Measure(CParameters& params);

//=============================================================================
// END OF FILE
//=============================================================================
