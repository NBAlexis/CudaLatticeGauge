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

#define m_fBetaOverNR static_cast<Real>(m_fBetaOverN)

__BEGIN_NAMESPACE

class CLGAPI CAction : public CBase
{
public:
    CAction() 
        : m_pOwner(NULL)
        , m_byActionId(0) 
        , m_fLastEnergy(0.0)
        , m_fNewEnergy(0.0)
        , m_fBetaOverN(0.1)
    {  
    }

    virtual ~CAction()
    {

    }

    /**
    * This is called langevin in Bridge++
    *
    * To set pStable, if energy can be calculate using pre-calculated stables (so this action must be gauge action)
    * bBeforeEvolution is set to be TRUE if call this just after an update (therefore the energy is already calculate for Metroplis of the last step),
    * and will return the last result.
    */
    virtual DOUBLE Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldGauge* const* stableFields);

    /**
    * Obtain the pointer of the fields
    */
    virtual void Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId);

    /**
    * Calculate the force on gauge fields to update the fields
    * Note that, gauge field will be changed, the change will be accepte or rejected.
    * So the "pGauge" must be a copy of real gauge field, not the orignal one!
    * Can fail due to solver
    */
    virtual UBOOL CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, 
        CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
        CFieldGauge* const* stapleFields, ESolverPhase ePhase) const;

    /**
    * Generate randoms
    * For pure gauge action, only in the first update, we need to calculate the energy, 
    * in latter updates, we just use the calculated energy depend on whether it is accepted.
    * see also OnFinishOneTrajotory
    */
    virtual void PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate);

    virtual void OnFinishTrajectory(UBOOL bAccepted)
    { 
        if (bAccepted && !IsFermion())
        {
            m_fLastEnergy = m_fNewEnergy;
        }
    }

    CCString GetInfos(const CCString &tab) const override;

    BYTE GetActionId() const { return m_byActionId; }

    virtual UBOOL IsFermion() const { return FALSE; }

    DOUBLE GetBetaOverN() const { return m_fBetaOverN; }

    void SetBeta(DOUBLE fBeta)
    {
        m_fBetaOverN = fBeta / static_cast<DOUBLE>(GetDefaultMatrixN());
    }

protected:

    virtual UINT GetDefaultMatrixN() const;

    virtual UBOOL CalculateForceOnGaugeSingleField(const class CFieldGauge* pGauge, class CFieldGauge* pForce, class CFieldGauge* pStaple, ESolverPhase ePhase) const
    {
        appCrucial(_T("CalculateForceOnGauge not implemented\n"));
        return FALSE;
    }

    virtual void PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
    {
        appCrucial(_T("PrepareForHMC not implemented\n"));
    }

    virtual DOUBLE EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
    {
        appCrucial(_T("Energy not implemented\n"));
        return 0.0;
    }

    class CLatticeData* m_pOwner;
    BYTE m_byActionId;

    DOUBLE m_fLastEnergy;
    DOUBLE m_fNewEnergy;
    DOUBLE m_fBetaOverN;

    TArray<BYTE> m_byGaugeFieldIds;
    TArray<BYTE> m_byBosonFieldIds;
};

__END_NAMESPACE

#endif //#ifndef _CACTION_H_

//=============================================================================
// END OF FILE
//=============================================================================