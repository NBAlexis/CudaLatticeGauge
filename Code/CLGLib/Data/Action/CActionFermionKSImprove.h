//=============================================================================
// FILENAME : CActionFermionKSImprove.h
// 
// DESCRIPTION:
// This is for the preparation of implementing improved KS action
// It changes the force calculation machinism
// by firstly computing X^+(n+mu)Y(n) and X^+(n)Y(n+mu)
// only work with one gauge field now
// 
//
// REVISION:
//  [mm/dd/yy]
//  [11/10/2024 nbale]
//=============================================================================
#pragma once

#include "CActionFermionKS.h"

#ifndef _CACTIONFERMIONKSIMPROVE_H_
#define _CACTIONFERMIONKSIMPROVE_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionFermionKSImprove)

class CLGAPI CActionFermionKSImprove : public CActionFermionKS
{
    __CLGDECLARE_CLASS(CActionFermionKSImprove)
public:

    /**
    * Make sure this is called after lattice and fields are created.
    */
    CActionFermionKSImprove();
    ~CActionFermionKSImprove();

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    UBOOL IsFermion() const override { return TRUE; }

    DOUBLE Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* stapleFields) override;
    UBOOL CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
        CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
        CFieldGauge* const* stapleFields, ESolverPhase ePhase) const override;
    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate) override;

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CActionFermionKS::GetInfos(tab);
        sRet = sRet + tab + _T("HISQ : ") + appToString(m_bHISQ) + _T("\n");
        return sRet;
    }

protected:

    //class CFieldFermion* m_pFerimionField;
    //class CFieldGauge* m_pEffectiveGauge;
    //class CFieldGauge* m_pForce;

    //HISQ is special to have a Naik term
    UBOOL m_bHISQ;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONFERMIONKSIMPROVE_H_

//=============================================================================
// END OF FILE
//=============================================================================