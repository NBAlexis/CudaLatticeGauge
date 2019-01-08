//=============================================================================
// FILENAME : CAction.h
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [12/4/2018 nbale]
//=============================================================================

#ifndef _CACTION_H_
#define _CACTION_H_

__BEGIN_NAMESPACE

class CLGAPI CAction : public CBase
{
public:
    CAction() : m_pOwner(NULL) { ; }

    /**
    * This is called langevin in Bridge++
    * This is S for specific configuration using for exp(-S)/exp(-S0) update
    * Because of the presence of Fermions, we can no longer just calculate a local(ultral-local) change os S
    */
    virtual Real Energy(const class CFieldGauge* pGauge) const = 0;

    /**
    * Obtain the pointer of the fields
    */
    virtual void Initial(class CLatticeData* pOwner, const CParameters& param) = 0;

    /**
    * Calculate the force on gauge fields to update the fields
    * Note that, gauge field will be changed, the change will be accepte or rejected.
    * So the "pGauge" must be a copy of real gauge field, not the orignal one!
    */
    virtual void CalculateForceOnGauge(class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple) const = 0;

protected:

    class CLatticeData* m_pOwner;
};

__END_NAMESPACE

#endif //#ifndef _CACTION_H_

//=============================================================================
// END OF FILE
//=============================================================================