//=============================================================================
// FILENAME : CActionFermionKS.h
// 
// DESCRIPTION:
// This is a naive KS fermion action, for Nf=1 or 2, or 1+1 etc,
// Fermions are independent, not optimized for Nf=2+1 with a heavy fermion
//
// REVISION:
//  [06/30/2020 nbale]
//=============================================================================

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
    CCString GetInfos(const CCString &tab) const override;
    UBOOL IsFermion() const override { return TRUE; }

    //Use in test only
    void SetFermionFieldTest(CFieldFermionKS* pField) { m_pFerimionField = pField; }

    DOUBLE Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldGauge* const* stableFields) override;
    UBOOL CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
        CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
        CFieldGauge* const* stapleFields, ESolverPhase ePhase) const override;
    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate) override;

protected:

    class CFieldFermionKS* m_pFerimionField;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONFERMIONKS_H_

//=============================================================================
// END OF FILE
//=============================================================================