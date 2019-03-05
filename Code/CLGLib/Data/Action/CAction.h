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
    *
    * To set pStable, if energy can be calculate using pre-calculated stables (so this action must be gauge action)
    */
    virtual Real Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable = NULL) = 0;

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
    * For pure gauge action, only in the first update, we need to calculate the energy, 
    * in latter updates, we just use the calculated energy depend on whether it is accepted.
    * see also OnFinishOneTrajotory
    */
    virtual void PrepareForHMC(const CFieldGauge* , UINT /* iUpdateIterate */) { ; }

    virtual void OnFinishTrajectory(UBOOL /* bAccepted */) { ; }

    virtual CCString GetInfos(const CCString &tab) const = 0;

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