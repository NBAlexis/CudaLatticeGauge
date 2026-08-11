//=============================================================================
// FILENAME : RotationImprovedFermion.h
// 
// DESCRIPTION:
//
// REVISION:
// [mm/dd/yy]
// [08/22/2025 nbale]
//=============================================================================
#pragma once

#include "CLGLib.h"

__DEFINE_ENUM(ERotationImproveFermionJob,
    ERIF_SimulateHISQ,
    ERIF_SimulateHISQRotation,
    ERIF_SimulateStoutlink,
    ERIF_SimulateStoutlinkRotation,
    ERIF_MeasureHISQ,
    ERIF_MeasureHISQRotation,
    )

/**
 * Measurement jobs shared by the two measure branches.
 * - ERIFMJ_Polyakov : measure Polyakov loop (CMeasurePolyakovXY)
 * - ERIFMJ_Chiral   : measure chiral condensate (CMeasureChiralCondensateKS, stochastic source)
 */
__DEFINE_ENUM(ERIFMeasureJob,
    ERIFMJ_Polyakov,
    ERIFMJ_Chiral,
    )

extern INT Simulate(CParameters& param);
extern INT SimulateRotation(CParameters& param);
extern INT MeasureHISQ(CParameters& param);
extern INT MeasureHISQRotation(CParameters& param);

//=============================================================================
// END OF FILE
//=============================================================================
