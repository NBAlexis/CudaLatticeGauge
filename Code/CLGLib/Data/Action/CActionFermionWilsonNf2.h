//=============================================================================
// FILENAME : CActionFermionWilsonNf2.h
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [02/06/2019 nbale]
//=============================================================================

#ifndef _CACTIONFERMIONWILSONNF2_H_
#define _CACTIONFERMIONWILSONNF2_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionFermionWilsonNf2)

class CLGAPI CActionFermionWilsonNf2 : public CAction
{
    __CLGDECLARE_CLASS(CActionFermionWilsonNf2)
public:
    /**
    * Make sure this is called after lattice and fields are created.
    */
    CActionFermionWilsonNf2();

    Real Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) override;
    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    UBOOL CalculateForceOnGauge(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const override;
    void PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate) override;
    CCString GetInfos(const CCString &tab) const override;

protected:

    class CFieldFermionWilsonSquareSU3* m_pFerimionField;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONFERMIONWILSONNF2_H_

//=============================================================================
// END OF FILE
//=============================================================================