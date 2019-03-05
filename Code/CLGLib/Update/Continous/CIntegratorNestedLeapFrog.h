//=============================================================================
// FILENAME : CIntegratorOmelyan.h
// 
// DESCRIPTION:
// This is the Omelyan integrator for HMC
//
// REVISION:
//  [12/12/2018 nbale]
//=============================================================================

#ifndef _CINTEGRATORNESTEDLEAPFORG_H_
#define _CINTEGRATORNESTEDLEAPFORG_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CIntegratorNestedLeapFrog)

class CLGAPI CIntegratorNestedLeapFrog : public CNestedIntegrator
{
    __CLGDECLARE_CLASS(CIntegratorNestedLeapFrog)

public:

    CIntegratorNestedLeapFrog() : CNestedIntegrator() {}

    virtual void Evaluate();

    void NestedEvaluate(UBOOL bLast);

    virtual CCString GetInfos(const CCString& sTab) const;

};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATORNESTEDLEAPFORG_H_

//=============================================================================
// END OF FILE
//=============================================================================