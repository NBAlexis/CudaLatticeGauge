//=============================================================================
// FILENAME : CIntegratorOmelyan.h
// 
// DESCRIPTION:
// This is the Omelyan integrator for HMC
//
// REVISION:
//  [12/12/2018 nbale]
//=============================================================================

#ifndef _CINTEGRATORFORCEGRADIENT_H_
#define _CINTEGRATORFORCEGRADIENT_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CIntegratorForceGradient)

class CLGAPI CIntegratorForceGradient : public CIntegrator
{
    __CLGDECLARE_CLASS(CIntegratorForceGradient)

public:

    CIntegratorForceGradient() : CIntegrator(), m_pUPrime(NULL) {}
    ~CIntegratorForceGradient();

    virtual void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params);

    virtual void Evaluate();

    virtual CCString GetInfos(const CCString& sTab) const;

protected:

    class CFieldGauge * m_pUPrime;
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATORFORCEGRADIENT_H_

//=============================================================================
// END OF FILE
//=============================================================================