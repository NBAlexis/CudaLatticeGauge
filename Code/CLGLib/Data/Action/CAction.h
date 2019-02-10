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
    CAction() : m_pOwner(NULL), m_byActionId(0) { ; }

    /**
    * This is called langevin in Bridge++
    * This is S for specific configuration using for exp(-S)/exp(-S0) update
    * Because of the presence of Fermions, we can no longer just calculate a local(ultral-local) change os S
    */
    virtual Real Energy(const class CFieldGauge* pGauge) const = 0;

    /**
    * \brief Special version of Energy using pre-calculated stables (so this action must be gauge action)
    *
    * The default implementation is just call Energy. For those support calculate using stable, override this function
    */
    virtual Real Energy(const class CFieldGauge* pGauge, const class CFieldGauge* ) const
    {
        //The default implementation is just call Energy. For those support calculate using stable, override this function
        return Energy(pGauge);
    }

    /**
    * Obtain the pointer of the fields
    */
    virtual void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId) = 0;

    /**
    * Calculate the force on gauge fields to update the fields
    * Note that, gauge field will be changed, the change will be accepte or rejected.
    * So the "pGauge" must be a copy of real gauge field, not the orignal one!
    * Can fail due to solver
    */
    virtual UBOOL CalculateForceOnGauge(const class CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple) const = 0;

    /**
    * Generate randoms
    */
    virtual void PrepareForHMC(const CFieldGauge* ) { ; }

    BYTE GetActionId() const { return m_byActionId; }

protected:

    class CLatticeData* m_pOwner;
    BYTE m_byActionId;
};

__END_NAMESPACE

#endif //#ifndef _CACTION_H_

//=============================================================================
// END OF FILE
//=============================================================================