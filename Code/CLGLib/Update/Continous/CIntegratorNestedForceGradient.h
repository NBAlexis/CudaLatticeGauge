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

    CIntegratorNestedForceGradient() : CNestedIntegrator(), m_pUPrime(NULL) {}
    ~CIntegratorNestedForceGradient();

    virtual void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params);

    virtual void Evaluate();

    void NestedEvaluate(UBOOL bLast);

    virtual CCString GetInfos(const CCString& sTab) const;

protected:

    class CFieldGauge * m_pUPrime;
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATORNESTEDLEAPFORG_H_

//=============================================================================
// END OF FILE
//=============================================================================