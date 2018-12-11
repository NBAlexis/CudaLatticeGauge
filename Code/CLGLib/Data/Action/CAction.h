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

class CLGAPI CAction
{
public:
    CAction()
        : m_pOwner(NULL)
        , m_pDeviceLattce(NULL){ ; }

    /**
    * This is called langevin in Bridge++
    * This is S for specific configuration using for exp(-S)/exp(-S0) update
    * Because of the presence of Fermions, we can no longer just calculate a local(ultral-local) change os S
    */
    virtual FLOAT Energy() = 0;

    /**
    * Obtain the pointer of the fields
    */
    virtual void Prepare() = 0;

    /**
    * Calculate the force on gauge fields to update the fields
    */
    virtual void CalculateForceOnGauge(class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple) const = 0;

protected:

    class CLatticeData* m_pOwner;
    class CDeviceLattice* m_pDeviceLattce;
};

__END_NAMESPACE

#endif //#ifndef _CACTION_H_

//=============================================================================
// END OF FILE
//=============================================================================