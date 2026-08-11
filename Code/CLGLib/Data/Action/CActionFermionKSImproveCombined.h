//=============================================================================
// FILENAME : CActionFermionKSImproveCombined.h
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [mm/dd/yy]
//  [08/02/2025 nbale]
//=============================================================================
#pragma once

#include "CActionFermionKSCombined.h"

#ifndef _CACTIONFERMIONKSIMPROVECOMBINED_H_
#define _CACTIONFERMIONKSIMPROVECOMBINED_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionFermionKSImproveCombined)

class CLGAPI CActionFermionKSImproveCombined : public CActionFermionKSCombined
{
    __CLGDECLARE_CLASS(CActionFermionKSImproveCombined)
public:

    /**
    * Make sure this is called after lattice and fields are created.
    */
    CActionFermionKSImproveCombined();

    DOUBLE Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* stapleFields) override;
    UBOOL CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
        CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
        CFieldGauge* const* stapleFields, ESolverPhase ePhase) const override;
    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate) override;

};

__END_NAMESPACE

#endif //#ifndef _CACTIONFERMIONKSIMPROVECOMBINED_H_

//=============================================================================
// END OF FILE
//=============================================================================