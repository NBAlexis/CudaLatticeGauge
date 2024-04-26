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

    DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable) override;

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;
    CCString GetInfos(const CCString &tab) const override;
    UBOOL IsFermion() const override { return TRUE; }
    class CFieldFermionKS* m_pFerimionField;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONFERMIONKS_H_

//=============================================================================
// END OF FILE
//=============================================================================