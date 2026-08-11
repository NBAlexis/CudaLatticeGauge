//=============================================================================
// FILENAME : CActionTemperatureDistribution.h
// 
// DESCRIPTION:
// This is like At Gradient, but use a real-boson instead
// since z-boundary is strange, not support Dirichlet gauge field anymore, this simplifies significantly the calculation of force
//
// REVISION:
//  [mm/dd/yy]
//  [11/04/2024 nbale]
//=============================================================================
#pragma once

#ifndef _CACTIONTEMPERATUREDISTRIBUTION_H_
#define _CACTIONTEMPERATUREDISTRIBUTION_H_

#include "Data/Field/Boson/CFieldBosonReal.h"

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionTemperatureDistribution)

class CLGAPI CActionTemperatureDistribution : public CAction
{
    __CLGDECLARE_CLASS(CActionTemperatureDistribution)
public:
    CActionTemperatureDistribution();

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;

    DOUBLE Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldGauge* const* stapleFields) override;

    //virtual void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId);

    UBOOL CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
        CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
        CFieldGauge* const* stapleFields, ESolverPhase ePhase) const override;

    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate) override;

    CCString GetInfos(const CCString& tab) const override;

    DOUBLE m_fBosonKinetic;
    class CFieldBosonReal* m_pStaticBosonField;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONTEMPERATUREDISTRIBUTION_H_

//=============================================================================
// END OF FILE
//=============================================================================