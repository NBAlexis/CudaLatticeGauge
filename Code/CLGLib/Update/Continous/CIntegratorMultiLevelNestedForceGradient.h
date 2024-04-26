//=============================================================================
// FILENAME : CIntegratorMultiLevelNestedForceGradient.h
// 
// DESCRIPTION:
// This is the Force Gradient integrator for HMC
//
// REVISION:
//  [08/17/2020 nbale]
//=============================================================================

#ifndef _CINTEGRATORMULTILEVELNESTEDFORCEGRADIENT_H_
#define _CINTEGRATORMULTILEVELNESTEDFORCEGRADIENT_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CIntegratorMultiLevelNestedForceGradient)

class CLGAPI CIntegratorMultiLevelNestedForceGradient : public CMultiLevelNestedIntegrator
{
    __CLGDECLARE_CLASS(CIntegratorMultiLevelNestedForceGradient)

public:

    CIntegratorMultiLevelNestedForceGradient() : CMultiLevelNestedIntegrator() {}

    void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params) override;

    void Evaluate() override;
    CCString GetInfos(const CCString& sTab) const override;

protected:

    void NestedEvaluate(INT iLevel, Real fNestedStep, UBOOL bFirst, UBOOL bLast);

};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATORMULTILEVELNESTEDFORCEGRADIENT_H_

//=============================================================================
// END OF FILE
//=============================================================================