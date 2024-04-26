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

    CIntegratorForceGradient() : CIntegrator() {}

    void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params) override;

    void Evaluate() override;

    CCString GetInfos(const CCString& sTab) const override;

protected:


};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATORFORCEGRADIENT_H_

//=============================================================================
// END OF FILE
//=============================================================================