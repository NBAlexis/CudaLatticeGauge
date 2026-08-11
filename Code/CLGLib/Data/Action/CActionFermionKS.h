//=============================================================================
// FILENAME : CActionFermionKS.h
// 
// DESCRIPTION:
// This is a naive KS fermion action, for Nf=1 or 2, or 1+1 etc,
// Fermions are independent, not optimized for Nf=2+1 with a heavy fermion
// 
// For historical reason, this was implemented for only staggered fermion, therefore, with a name 'CActionFermionKS'
// However, this is in fact the RHMC, which now supports Wilson-Dirac fermion despite the name 'KS'
//
// REVISION:
//  [mm/dd/yy]
//  [06/30/2020 nbale]
//=============================================================================
#pragma once

#ifndef _CACTIONFERMIONKS_H_
#define _CACTIONFERMIONKS_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionFermionKS)

class CLGAPI CActionFermionKS : public CAction
{
    __CLGDECLARE_CLASS(CActionFermionKS)

public:

    /**
    * Make sure this is called after lattice and fields are created.
    */
    CActionFermionKS();

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    UBOOL IsFermion() const override { return TRUE; }

    //Use in test only
    void SetFermionFieldTest(CFieldFermion* pField) { m_pFerimionField = pField; }

    DOUBLE Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, INT tensor2Num, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* stapleFields) override;
    UBOOL CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
        CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
        CFieldGauge* const* stapleFields, ESolverPhase ePhase) const override;
    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate) override;

    CCString GetInfos(const CCString& tab) const override
    {
        CCString sRet = CAction::GetInfos(tab);
        sRet = sRet + tab + _T("Fermion : \n") + m_pFerimionField->GetInfos(tab + _T("    ")) + _T("\n");
        return sRet;
    }

protected:

    class CFieldFermion* m_pFerimionField;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONFERMIONKS_H_

//=============================================================================
// END OF FILE
//=============================================================================