//=============================================================================
// FILENAME : CIntegratorMultiLevelOmelyan.h
// 
// DESCRIPTION:
// This is the Force Gradient integrator for HMC
//
// REVISION:
//  [08/18/2020 nbale]
//=============================================================================

#ifndef _CINTEGRATORMULTILEVELNESTEDOMELYAN_H_
#define _CINTEGRATORMULTILEVELNESTEDOMELYAN_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CIntegratorMultiLevelOmelyan)

class CLGAPI CIntegratorMultiLevelOmelyan : public CMultiLevelNestedIntegrator
{
    __CLGDECLARE_CLASS(CIntegratorMultiLevelOmelyan)

public:

    CIntegratorMultiLevelOmelyan()
        : CMultiLevelNestedIntegrator()
        , m_f2Lambda(OmelyanLambda2)
    {
        
    }

    void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params) override;

    void Evaluate() override;
    CCString GetInfos(const CCString& sTab) const override;

protected:

    void NestedEvaluate(INT iLevel, Real fNestedStepLength, UBOOL bFirst, UBOOL bLast);

    Real m_f2Lambda;
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATORMULTILEVELNESTEDOMELYAN_H_

//=============================================================================
// END OF FILE
//=============================================================================