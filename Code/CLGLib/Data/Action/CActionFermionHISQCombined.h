//=============================================================================
// FILENAME : CActionFermionHISQCombined.h
// 
// DESCRIPTION:
// This only supports Even-odd fermion, and EO fermion with phase
//
// REVISION:
//  [mm/dd/yy]
//  [01/26/2025 nbale]
//=============================================================================
#pragma once

#include "Data/Field/Staggered/CFieldFermionKSHISQ.h"

#ifndef _CACTIONFERMIONHISQCOMBINED_H_
#define _CACTIONFERMIONHISQCOMBINED_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionFermionHISQCombined)

class CLGAPI CActionFermionHISQCombined : public CAction
{
    __CLGDECLARE_CLASS(CActionFermionHISQCombined)

public:

    /**
    * Make sure this is called after lattice and fields are created.
    */
    CActionFermionHISQCombined();

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    UBOOL IsFermion() const override { return TRUE; }

    DOUBLE Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* stapleFields) override;
    UBOOL CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
        CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
        CFieldGauge* const* stapleFields, ESolverPhase ePhase) const override;
    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate) override;

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CAction::GetInfos(tab);
        for (INT i = 0; i < m_pFerimionField.Num(); ++i)
        {
            sRet = sRet + tab + _T("Fermion ") + appToString(i) + _T(": \n") + m_pFerimionField[i]->GetInfos(tab + _T("    ")) + _T("\n");
        }
        return sRet;
    }

protected:

    class TArray<CFieldFermion*> m_pFerimionField;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONFERMIONHISQCOMBINED_H_

//=============================================================================
// END OF FILE
//=============================================================================