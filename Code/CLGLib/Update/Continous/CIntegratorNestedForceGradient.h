//=============================================================================
// FILENAME : CIntegratorOmelyan.h
// 
// DESCRIPTION:
// This is the Omelyan integrator for HMC
//
// REVISION:
//  [12/12/2018 nbale]
//=============================================================================

#ifndef _CINTEGRATORNESTEDFORCEGRADIENT_H_
#define _CINTEGRATORNESTEDFORCEGRADIENT_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CIntegratorNestedForceGradient)

class CLGAPI CIntegratorNestedForceGradient : public CNestedIntegrator
{
    __CLGDECLARE_CLASS(CIntegratorNestedForceGradient)

public:

    CIntegratorNestedForceGradient() : CNestedIntegrator() {}

    void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params) override;

    void Evaluate() override;

    void NestedEvaluate(UBOOL bLast);

    CCString GetInfos(const CCString& sTab) const override;

    void ChangeStepCount(UBOOL bGrow) override
    {
        CNestedIntegrator::ChangeStepCount(bGrow);
        m_fNestedStepLength = F(0.5) * m_fNestedStepLength;
    }

protected:

};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATORNESTEDLEAPFORG_H_

//=============================================================================
// END OF FILE
//=============================================================================