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
        
    {  
    }

    ~CIntegrator();

    virtual void Evaluate() = 0;
    virtual void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params);
    void OnFinishTrajectory(UBOOL bAccepted);
    void Prepare(UBOOL bLastAccepted, UINT uiStep);
    void UpdateP(Real fStep, UBOOL bCacheStaple, ESolverPhase ePhase);
    void UpdateU(Real fStep) const;
    virtual DOUBLE GetEnergy(UBOOL bBeforeEvolution) const;

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

    void ZeroForce() const
    {
        for (INT i = 0; i < m_pForceField.Num(); ++i)
        {
            m_pForceField[i]->Zero();
        }

        for (INT i = 0; i < m_pBosonForceFields.Num(); ++i)
        {
            m_pBosonForceFields[i]->Zero();
        }
    }

    void AddForce(Real fStep, UBOOL bBind) const
    {
        for (INT i = 0; i < m_pMomentumField.Num(); ++i)
        {
            m_pMomentumField[i]->Axpy(fStep, m_pForceField[i]);
            if (bBind)
            {
                m_pMomentumField[i]->SetOneDirectionZero(m_byBindDir);
            }
        }

        for (INT i = 0; i < m_pBosonMomentumFields.Num(); ++i)
        {
            m_pBosonMomentumFields[i]->Axpy(fStep, m_pBosonForceFields[i]);
        }
    }

    /**
    * The next step is to update momentum using U'
    * this function calculate m_pGaugeField and m_pBosonFields (they are U') using 'force'
    */
    void AddForceToFieldDirectly(Real fStep) const
    {
        for (INT i = 0; i < m_pForceField.Num(); ++i)
        {
            m_pForceField[i]->ExpMult(fStep, m_pGaugeField[i]);
        }

        for (INT i = 0; i < m_pBosonForceFields.Num(); ++i)
        {
            m_pBosonFields[i]->Axpy(fStep, m_pBosonForceFields[i]);
        }
    }

    static void SetOneDirZero(const TArray<CFieldGauge*>& fields, BYTE byDir)
    {
        for (INT i = 0; i < fields.Num(); ++i)
        {
            fields[i]->SetOneDirectionZero(byDir);
        }
    }

    static void SetOneDirOne(const TArray<CFieldGauge*>& fields, BYTE byDir)
    {
        for (INT i = 0; i < fields.Num(); ++i)
        {
            fields[i]->SetOneDirectionUnity(byDir);
        }
    }

    Real CalcForce() const
    {
        Real fRet = F(0.0);

        for (INT i = 0; i < m_pForceField.Num(); ++i)
        {
#if !_CLG_DOUBLEFLOAT
            const CLGComplex force = _cToFloat(m_pForceField[i]->Dot(m_pForceField[i]));
#else
            const CLGComplex force = m_pForceField->Dot(m_pForceField);
#endif
            fRet += force.x;
        }

        for (INT i = 0; i < m_pBosonFields.Num(); ++i)
        {
#if !_CLG_DOUBLEFLOAT
            const CLGComplex force = _cToFloat(m_pBosonFields[i]->Dot(m_pBosonFields[i]));
#else
            const CLGComplex force = m_pBosonFields->Dot(m_pBosonFields);
#endif
            fRet += force.x;
        }
        return fRet;
    }

    DOUBLE CalcMomentumEnery() const
    {
        DOUBLE ret = 0.0;
        for (INT i = 0; i < m_pMomentumField.Num(); ++i)
        {
            ret += m_pMomentumField[i]->CalculateKinematicEnergy();
        }

        for (INT i = 0; i < m_pBosonMomentumFields.Num(); ++i)
        {
            ret += m_pBosonMomentumFields[i]->Dot(m_pBosonMomentumFields[i]).x;
        }
        return ret;
    }

    void InitialMomentumNoise() const
    {
        for (INT i = 0; i < m_pMomentumField.Num(); ++i)
        {
            m_pMomentumField[i]->MakeRandomGenerator();
            m_pMomentumField[i]->SetOneDirectionZero(m_byBindDir);
        }

        for (INT i = 0; i < m_pBosonMomentumFields.Num(); ++i)
        {
            m_pBosonMomentumFields[i]->MakeRandomMomentum();
        }
    }

    void PreserveFields() const
    {
        for (INT i = 0; i < m_pGaugeField.Num(); ++i)
        {
            m_pGaugeField[i]->CopyTo(m_pUPrime[i]);
        }

        for (INT i = 0; i < m_pBosonFields.Num(); ++i)
        {
            m_pBosonFields[i]->CopyTo(m_pPhiPrime[i]);
        }
    }

    void RecoverFields() const
    {
        for (INT i = 0; i < m_pGaugeField.Num(); ++i)
        {
            m_pUPrime[i]->CopyTo(m_pGaugeField[i]);
        }

        for (INT i = 0; i < m_pBosonFields.Num(); ++i)
        {
            m_pPhiPrime[i]->CopyTo(m_pBosonFields[i]);
        }
    }

    void CreateBackupFields()
    {
        for (INT i = 0; i < m_pGaugeField.Num(); ++i)
        {
            m_pUPrime.AddItem(dynamic_cast<CFieldGauge*>(m_pGaugeField[i]->GetCopy()));
        }

        for (INT i = 0; i < m_pBosonFields.Num(); ++i)
        {
            m_pPhiPrime.AddItem(dynamic_cast<CFieldBoson*>(m_pBosonFields[i]->GetCopy()));
        }
    }

public:

    TArray <class CFieldGauge*> m_pGaugeField;
    TArray <class CFieldGauge*> m_pForceField;
    TArray <class CFieldGauge*> m_pMomentumField;
    TArray <class CFieldGauge*> m_pStapleField;

    TArray<class CFieldBoson*> m_pBosonFields;
    TArray<class CFieldBoson*> m_pBosonForceFields;
    TArray<class CFieldBoson*> m_pBosonMomentumFields;

protected:

    /**
    * For force gradient updaters
    */
    TArray<class CFieldGauge*> m_pUPrime;
    TArray<class CFieldBoson*> m_pPhiPrime;
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