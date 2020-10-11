//=============================================================================
// FILENAME : CIntegrator.h
// 
// DESCRIPTION:
// This is the class for hibrid Monte Carlo
//
// REVISION:
//  [12/8/2018 nbale]
//=============================================================================

#ifndef _CINTEGRATOR_H_
#define _CINTEGRATOR_H_

__BEGIN_NAMESPACE

__DEFINE_ENUM(EIntegratorType,
    EIT_LeapFrog,
    EIT_Omelyan,
    EIT_Max,

    EIT_ForceDWORD = 0x7fffffff,
    )

class CLGAPI CIntegrator : public CBase
{
public:
    CIntegrator() 
        : m_byBindDir(0)
        , m_uiStepCount(1)
        , m_fEStep(0)
        , m_bDebugForce(FALSE)
        , m_bStapleCached(FALSE)
        , m_fUpdateResultEnery(F(0.0))
        , m_pOwner(NULL)
        , m_pLattice(NULL)
        , m_pGaugeField(NULL)
        , m_pForceField(NULL)
        , m_pMomentumField(NULL)
        , m_pStapleField(NULL)
        
    {  }

    ~CIntegrator();

    virtual void Evaluate() = 0;
    virtual void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params);
    void OnFinishTrajectory(UBOOL bAccepted);
    void Prepare(UBOOL bLastAccepted, UINT uiStep);
    void UpdateP(Real fStep, UBOOL bCacheStaple, ESolverPhase ePhase);
    void UpdateU(Real fStep) const;
#if !_CLG_DOUBLEFLOAT
    virtual DOUBLE GetEnergy(UBOOL bBeforeEvolution) const;
#else
    virtual Real GetEnergy(UBOOL bBeforeEvolution) const;
#endif
    virtual CCString GetInfos(const CCString& sTab) const = 0;
    void FinishEvaluate() const;
    virtual void ChangeStepCount(UBOOL bGrow) 
    {
        if (bGrow)
        {
            m_fEStep = m_fEStep * m_uiStepCount / (m_uiStepCount + F(1.0));
            m_uiStepCount++;
        }
        else if (m_uiStepCount > 2)
        {
            m_fEStep = m_fEStep * m_uiStepCount / (m_uiStepCount - F(1.0));
            m_uiStepCount--;
        }
    }
    UINT GetStepCount() const { return m_uiStepCount; }

protected:

    BYTE m_byBindDir;
    UINT m_uiStepCount;
    Real m_fEStep;
    UBOOL m_bDebugForce;
    UBOOL m_bStapleCached;

    Real m_fUpdateResultEnery;

    class CHMC* m_pOwner;
    CLatticeData* m_pLattice;
    TArray<class CAction*> m_lstActions;

public:

    CFieldGauge* m_pGaugeField;
    CFieldGauge* m_pForceField;
    CFieldGauge* m_pMomentumField;
    CFieldGauge* m_pStapleField;


};

class CLGAPI CNestedIntegrator : public CIntegrator
{
public:
    CNestedIntegrator()
        : CIntegrator()
        , m_uiNestedStep(1)
        , m_fNestedStepLength(F(0.0))
        , m_bInnerLeapFrog(FALSE)
    {
    }

    void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params) override;
    CCString GetNestedInfo(const CCString & sTab) const;

    void UpdatePF(Real fStep, ESolverPhase ePhase);

    //Gauge force is irrelevant from solver
    void UpdatePG(Real fStep, UBOOL bCacheStaple);

    void ChangeStepCount(UBOOL bGrow) override
    {
        CIntegrator::ChangeStepCount(bGrow);
        m_fNestedStepLength = m_fEStep / m_uiNestedStep;
    }

protected:

    void NestedEvaluateLeapfrog(UBOOL bLast);

    UINT m_uiNestedStep;
    Real m_fNestedStepLength;
    UBOOL m_bInnerLeapFrog;
};

class CLGAPI CMultiLevelNestedIntegrator : public CIntegrator
{
public:
    CMultiLevelNestedIntegrator()
        : CIntegrator()
        , m_fTotalStepLength(F(1.0))
        , m_bInnerLeapFrog(FALSE)
    {
    }

    void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params) override;
    CCString GetNestedInfo(const CCString& sTab) const;

    void ChangeStepCount(UBOOL bGrow) override
    {
        CIntegrator::ChangeStepCount(bGrow);
        m_fNestedStepLengths.RemoveAll();
        m_fNestedStepLengths.AddItem(m_fEStep);
        Real fStep = m_fEStep;
        for (INT i = 0; i < m_uiNestedStep.Num(); ++i)
        {
            fStep = fStep / m_uiNestedStep[i];
            m_fNestedStepLengths.AddItem(fStep);
        }
    }

    /**
     * In force gradiant, sometimes we only cauclate pForce, but not update Momentum
     * So there is a 'bUpdateP'
     */
    void UpdateP(Real fStep, TArray<UINT> actionList, ESolverPhase ePhase, UBOOL bCacheStaple, UBOOL bUpdateP);
    void UpdateP(Real fStep, INT iLevel, ESolverPhase ePhase, UBOOL bCacheStaple, UBOOL bUpdateP)
    {
        UpdateP(fStep, m_iNestedActionId[iLevel], ePhase, bCacheStaple, bUpdateP);
    }

protected:

    void NestedEvaluateLeapfrog(INT iLevel, Real fNestedStepLength, UBOOL bFirst, UBOOL bLast);

    TArray<UINT> m_uiNestedStep;
    TArray<Real> m_fNestedStepLengths;
    TArray<TArray<UINT>> m_iNestedActionId;
    Real m_fTotalStepLength;
    UBOOL m_bInnerLeapFrog;
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================