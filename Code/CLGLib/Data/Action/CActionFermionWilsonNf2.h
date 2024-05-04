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

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;
    CCString GetInfos(const CCString &tab) const override;
    UBOOL IsFermion() const override { return TRUE; }

    //Use in test only
    void SetFermionFieldTest(CFieldFermionWilsonSquareSU3* pField) { m_pFerimionField = pField; }

    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate) override;
    DOUBLE Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldGauge* const* stableFields) override;
    UBOOL CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
        CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
        CFieldGauge* const* stapleFields, ESolverPhase ePhase) const override;

protected:

    class CFieldFermionWilsonSquareSU3* m_pFerimionField;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONFERMIONWILSONNF2_H_

//=============================================================================
// END OF FILE
//=============================================================================