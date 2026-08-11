//=============================================================================
// FILENAME : AgentDriver.h
// 
// DESCRIPTION:
// Unified YAML-driven driver for agent use.
// Supports three tasks: Simulate, GaugeFixing, Measure.
//
// REVISION:
//  [06/14/2026]
//=============================================================================

#ifndef _AGENTDRIVER_H_
#define _AGENTDRIVER_H_

#include "CLGLib.h"

__BEGIN_NAMESPACE

INT RunAgentSimulate(CParameters& params);
INT RunAgentGaugeFixing(CParameters& params);
INT RunAgentMeasure(CParameters& params);

__END_NAMESPACE

#endif //#ifndef _AGENTDRIVER_H_

//=============================================================================
// END OF FILE
//=============================================================================
