//=============================================================================
// FILENAME : CIntegratorOmelyan.h
// 
// DESCRIPTION:
// This is the Omelyan integrator for HMC
//
// REVISION:
//  [12/12/2018 nbale]
//=============================================================================

#ifndef _CINTEGRATORNESTEDOMELYAN_H_
#define _CINTEGRATORNESTEDOMELYAN_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CIntegratorNestedOmelyan)

class CLGAPI CIntegratorNestedOmelyan : public CNestedIntegrator
{
    __CLGDECLARE_CLASS(CIntegratorNestedOmelyan)

public:

    CIntegratorNestedOmelyan() : CNestedIntegrator(), m_f2Lambda(OmelyanLambda2) {}

    virtual void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params);

    virtual void Evaluate();

    void NestedEvaluate(UBOOL bLast);

    virtual CCString GetInfos(const CCString& sTab) const;

protected:

    Real m_f2Lambda;
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATORNESTEDOMELYAN_H_

//=============================================================================
// END OF FILE
//=============================================================================