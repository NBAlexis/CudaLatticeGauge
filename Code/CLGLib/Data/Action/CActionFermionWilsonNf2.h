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

class CLGAPI CActionFermionWilsonNf2 : public CAction
{
    __CLGDECLARE_CLASS(CActionFermionWilsonNf2)
public:
    /**
    * Make sure this is called after lattice and fields are created.
    */
    CActionFermionWilsonNf2();

    virtual Real Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL);
    virtual void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId);
    virtual UBOOL CalculateForceOnGauge(UINT uiStep, const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple) const;
    virtual void PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate);
    virtual CCString GetInfos(const CCString &tab) const;

protected:

    class CFieldFermionWilsonSquareSU3* m_pFerimionField;
    UINT m_uiMutiStep;
};

__END_NAMESPACE

#endif //#ifndef _CACTIONFERMIONWILSONNF2_H_

//=============================================================================
// END OF FILE
//=============================================================================