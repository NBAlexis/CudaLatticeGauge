//=============================================================================
// FILENAME : CActionPhi4.h
// 
// DESCRIPTION:
// (nabla phi)*(nabla phi) + m phi* phi + lambda (phi*phi)^2
//
// REVISION:
//  [06/13/2024 nbale]
//=============================================================================

#ifndef _CACTIONPHI4_H_
#define _CACTIONPHI4_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CActionPhi4)

class CLGAPI CActionPhi4 : public CAction
{
    __CLGDECLARE_CLASS(CActionPhi4)
public:
    /**
    * Make sure this is called after lattice and fields are created.
    */
    CActionPhi4() : CAction()
        , m_fM(F(1.0))
        , m_fLambda(F(1.0))
    {

    }

    void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) override;

    DOUBLE Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldGauge* const* stableFields) override;

    UBOOL CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
        CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
        CFieldGauge* const* stapleFields, ESolverPhase ePhase) const override;

    void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate) override;

    Real m_fM;
    Real m_fLambda;

};

__END_NAMESPACE

#endif //#ifndef _CACTIONPHI4_H_

//=============================================================================
// END OF FILE
//=============================================================================