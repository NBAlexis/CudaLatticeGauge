//=============================================================================
// FILENAME : CIntegrator.h
// 
// DESCRIPTION:
// This is the class for hibrid Monte Carlo
//
// REVISION:
//  [mm/dd/yy]
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

__DEFINE_ENUM(EForceCalc,
    EFC_All,
    EFC_Gauge,
    EFC_Fermion,
    )

__DEFINE_ENUM(EBackupFieldType,
    EBFT_Same,
    EBFT_SU3_12,
    )

inline class CStapleCache* appGetStapleCache(BYTE byFieldId);

class CLGAPI CIntegrator : public CBase
{
public:
    CIntegrator()
        : m_byBindDir(0)
        , m_uiStepCount(1)
        , m_uiStepCountMetropolis(1)
        , m_uiStepCountWarmup(0)
        , m_fEStep(0)
        , m_bDebugForce(FALSE)
        //, m_bStapleCached(FALSE)
        , m_fUpdateResultEnery(F(0.0))
        , m_pOwner(NULL)
        , m_pLattice(NULL)
        , m_eBackupFieldType(EBFT_Same)
    {
        m_bUDirty[0] = TRUE;
        m_bUDirty[1] = TRUE;
    }

    ~CIntegrator();

    virtual void Evaluate() = 0;
    virtual void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params);
    void OnFinishTrajectory(UBOOL bAccepted);
    void Prepare(UBOOL bLastAccepted, UINT uiStep);
    void UpdateP(Real fStep, ESolverPhase ePhase);
    void UpdateU(Real fStep);
    virtual DOUBLE GetEnergy(UBOOL bBeforeEvolution, TArray<DOUBLE>& actions);

    void OnGaugeChanged()
    {
        m_bUDirty[0] = TRUE;
        m_bUDirty[1] = TRUE;
    }

    /**
    * This function just call gauge smearing,
    * the action or fermion handle the replacement with effective gauge fields
    * 
    * It is called,
    * before prepare gaussian fields (need to inverse)
    * before calculate energy
    * before calculate force on fermion
    */
    virtual void RequireGaugeSmearing();

    /**
    * if updateMode & 1, it was gauge update
    * if updateMode & 2, it was fermion update
    * Cache is always after gauge smearing because the cache can be used for smeared gauge
    */
    virtual void OnCacheStaple(ECacheCall eCall);

    /**
    * if updateMode & 1, it was gauge update
    * if updateMode & 2, it was fermion update
    * 
    * If smearing need to cache orignal gauge field for staple calculation:
    * Cache staple first cache plaqutte or staple for orignal gauge field
    * 
    * Then, do the smearing
    * 
    * Then, cache Fmunu, etc of effective gauge
    */
    virtual void OnCacheAndSmearing(INT updateMode)
    {
        if (0 == (updateMode & 2))
        {
            if (!m_bUDirty[0])
            {
                return;
            }
            //only gauge update, no gauge smearing
            OnCacheStaple(ECC_BeforeGaugeUpdate);
            m_bUDirty[0] = FALSE;
        }
        else
        {
            if (updateMode & 1) //both gauge and fermion
            {
                if (m_bUDirty[0])
                {
                    OnCacheStaple(ECC_BeforeGaugeUpdate);
                }
                if (m_bUDirty[1])
                {
                    OnCacheStaple(ECC_BeforeAllUpdateBeforeSmearing);
                    RequireGaugeSmearing();
                    OnCacheStaple(ECC_BeforeAllUpdateAfterSmearing);
                }

                m_bUDirty[0] = FALSE;
                m_bUDirty[1] = FALSE;
            }
            else //only fermion
            {
                if (!m_bUDirty[1])
                {
                    return;
                }
                m_bUDirty[1] = FALSE;
                OnCacheStaple(ECC_BeforeFermionUpdateBeforeSmearing);
                RequireGaugeSmearing();
                OnCacheStaple(ECC_BeforeFermionUpdateAfterSmearing);
            }
        }
    }

    virtual CCString GetInfos(const CCString& sTab) const = 0;
    void FinishEvaluate();
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

    virtual void ChangeStepCountTo(UINT uiStep)
    {
        m_fEStep = m_fEStep * m_uiStepCount / uiStep;
        m_uiStepCount = uiStep;
    }

    UINT GetStepCount() const { return m_uiStepCount; }
    BYTE GetBindDir() const { return m_byBindDir; }

    void FixStep(UBOOL bWarmup)
    {
        if (bWarmup && m_uiStepCount != m_uiStepCountWarmup && m_uiStepCountWarmup > 0)
        {
            ChangeStepCountTo(m_uiStepCountWarmup);
        }

        if (!bWarmup && m_uiStepCount != m_uiStepCountMetropolis)
        {
            ChangeStepCountTo(m_uiStepCountMetropolis);
        }
    }

protected:

    BYTE m_byBindDir;
    UINT m_uiStepCount;
    UINT m_uiStepCountMetropolis;
    UINT m_uiStepCountWarmup;
    Real m_fEStep;
    UBOOL m_bDebugForce;

    //staple cache to be discarded in the furture, use CStapleCache instead (if nesscerary)
    //UBOOL m_bStapleCached;

    Real m_fUpdateResultEnery;

    class CHMC* m_pOwner;
    CLatticeData* m_pLattice;

    //updator should never change actions, so set to be constant
    TArray<class CAction*> m_lstActions;

    void ZeroForce() const
    {
        for (INT i = 0; i < m_pForceField.Num(); ++i)
        {
            if (NULL != m_pForceField[i])
            {
                m_pForceField[i]->Zero();
            }
        }

        for (INT i = 0; i < m_pBosonForceFields.Num(); ++i)
        {
            if (NULL != m_pBosonForceFields[i])
            {
                m_pBosonForceFields[i]->Zero();
                //m_pBosonForceFields[i]->DebugPrintMe();
            }
        }
    }

    void AddForce(Real fStep, UBOOL bBind) const
    {
        FixGaugeBondary(m_pForceField, EFB_Force);
        for (INT i = 0; i < m_pMomentumField.Num(); ++i)
        {
            if (NULL != m_pMomentumField[i])
            {
                m_pMomentumField[i]->Axpy(fStep, m_pForceField[i]);
                if (bBind)
                {
                    m_pMomentumField[i]->SetOneDirectionZero(m_byBindDir);
                }
            }
        }
        FixGaugeBondary(m_pMomentumField, EFB_Momentum);

        FixBosonBondary(m_pBosonForceFields, EFB_Force);
        for (INT i = 0; i < m_pBosonMomentumFields.Num(); ++i)
        {
            if (NULL != m_pBosonMomentumFields[i])
            {
                m_pBosonMomentumFields[i]->Axpy(fStep, m_pBosonForceFields[i]);
            }
        }
        FixBosonBondary(m_pBosonMomentumFields, EFB_Momentum);
    }

    /**
    * The next step is to update momentum using U'
    * this function calculate m_pGaugeField and m_pBosonFields (they are U') using 'force'
    */
    void AddForceToFieldDirectly(Real fStep)
    {
        FixGaugeBondary(m_pForceField, EFB_Force);
        for (INT i = 0; i < m_pForceField.Num(); ++i)
        {
            if (NULL != m_pForceField[i])
            {
                m_pForceField[i]->ExpMult(fStep / _HC_GaugeMomentumFactor, m_pGaugeField[i]);
            }
        }
        FixGaugeBondary(m_pGaugeField, EFB_Field);

        FixBosonBondary(m_pBosonForceFields, EFB_Force);
        for (INT i = 0; i < m_pBosonForceFields.Num(); ++i)
        {
            if (NULL != m_pBosonFields[i])
            {
                m_pBosonFields[i]->Axpy(fStep / _HC_GaugeMomentumFactor, m_pBosonForceFields[i]);
            }
        }
        FixBosonBondary(m_pBosonFields, EFB_Field);
        OnGaugeChanged();
    }

    static void SetOneDirZero(const TArray<CFieldGauge*>& fields, BYTE byDir)
    {
        for (INT i = 0; i < fields.Num(); ++i)
        {
            if (NULL != fields[i])
            {
                fields[i]->SetOneDirectionZero(byDir);
            }
        }
    }

    static void SetOneDirOne(const TArray<CFieldGauge*>& fields, BYTE byDir)
    {
        for (INT i = 0; i < fields.Num(); ++i)
        {
            if (NULL != fields[i])
            {
                fields[i]->SetOneDirectionUnity(byDir);
            }
        }
    }

    static void FixBosonBondary(const TArray<CFieldBoson*>& fields, EFixBoundary eType)
    {
        for (INT i = 0; i < fields.Num(); ++i)
        {
            if (NULL != fields[i])
            {
                fields[i]->FixBoundary(eType);
            }
        }
    }

    static void FixGaugeBondary(const TArray<CFieldGauge*>& fields, EFixBoundary eType)
    {
        for (INT i = 0; i < fields.Num(); ++i)
        {
            if (NULL != fields[i])
            {
                fields[i]->FixBoundary(eType);
            }
        }
    }

    /**
    * actionlst is accually const TArray<const CAction*>&
    */
    void CalcForceOfActions(const TArray<CAction*>& actionlst, EForceCalc eMode, ESolverPhase ePhase);

    Real CalcForce() const
    {
        Real fRet = F(0.0);

        for (INT i = 0; i < m_pForceField.Num(); ++i)
        {
            if (NULL == m_pForceField[i])
            {
                continue;
            }
#if !_CLG_DOUBLEFLOAT
            const CLGComplex force = _cToFloat(m_pForceField[i]->Dot(m_pForceField[i]));
#else
            const CLGComplex force = m_pForceField[i]->Dot(m_pForceField[i]);
#endif
            fRet += force.x;
        }

        for (INT i = 0; i < m_pBosonFields.Num(); ++i)
        {
            if (NULL == m_pBosonFields[i])
            {
                continue;
            }
#if !_CLG_DOUBLEFLOAT
            const CLGComplex force = _cToFloat(m_pBosonFields[i]->Dot(m_pBosonFields[i]));
#else
            const CLGComplex force = m_pBosonFields[i]->Dot(m_pBosonFields[i]);
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
            if (NULL != m_pMomentumField[i])
            {
                //m_pMomentumField[i]->DebugPrintMe();
                DOUBLE kin = m_pMomentumField[i]->CalculateKinematicEnergy();
                appDetailed(_T("kin %d: %f\n"), i, kin);
                ret += kin;
            }
        }

        for (INT i = 0; i < m_pBosonMomentumFields.Num(); ++i)
        {
            if (NULL != m_pBosonMomentumFields[i])
            {
                ret += m_pBosonMomentumFields[i]->Dot(m_pBosonMomentumFields[i]).x;
            }
        }
        return ret;
    }

    void InitialMomentumNoise() const
    {
        for (INT i = 0; i < m_pMomentumField.Num(); ++i)
        {
            if (NULL != m_pMomentumField[i])
            {
                m_pMomentumField[i]->MakeRandomGenerator();
                m_pMomentumField[i]->SetOneDirectionZero(m_byBindDir);
            }
        }
        FixGaugeBondary(m_pMomentumField, EFB_Momentum);

        for (INT i = 0; i < m_pBosonMomentumFields.Num(); ++i)
        {
            if (NULL != m_pBosonMomentumFields[i])
            {
                m_pBosonMomentumFields[i]->MakeRandomMomentum();
            }
        }
        FixBosonBondary(m_pBosonMomentumFields, EFB_Momentum);
        //m_pBosonMomentumFields[0]->DebugPrintMe();
    }

    void PreserveFields() const
    {
        for (INT i = 0; i < m_pGaugeField.Num(); ++i)
        {
            if (NULL != m_pGaugeField[i])
            {
                m_pGaugeField[i]->CopyTo(m_pUPrime[i]);
            }
        }

        for (INT i = 0; i < m_pBosonFields.Num(); ++i)
        {
            if (NULL != m_pBosonFields[i])
            {
                m_pBosonFields[i]->CopyTo(m_pPhiPrime[i]);
            }
        }
    }

    void RecoverFields()
    {
        for (INT i = 0; i < m_pGaugeField.Num(); ++i)
        {
            if (NULL != m_pGaugeField[i])
            {
                m_pUPrime[i]->CopyTo(m_pGaugeField[i]);
            }
        }

        for (INT i = 0; i < m_pBosonFields.Num(); ++i)
        {
            if (NULL != m_pBosonFields[i])
            {
                m_pPhiPrime[i]->CopyTo(m_pBosonFields[i]);
            }
        }
        OnGaugeChanged();
    }

    void CreateBackupFields()
    {
        for (INT i = 0; i < m_pGaugeField.Num(); ++i)
        {
            if (NULL != m_pGaugeField[i])
            {
                if (EBFT_SU3_12 == m_eBackupFieldType)
                {
                    CFieldGaugeSU3_12* pBackup = new CFieldGaugeSU3_12();
                    m_pGaugeField[i]->CopyTo(pBackup);
                    m_pUPrime.AddItem(pBackup);
                }
                else
                {
                    m_pUPrime.AddItem(dynamic_cast<CFieldGauge*>(m_pGaugeField[i]->GetCopy()));
                }
            }
            else
            {
                m_pUPrime.AddItem(NULL);
            }
        }

        for (INT i = 0; i < m_pBosonFields.Num(); ++i)
        {
            if (NULL != m_pBosonFields[i])
            {
                m_pPhiPrime.AddItem(dynamic_cast<CFieldBoson*>(m_pBosonFields[i]->GetCopy()));
            }
            else
            {
                m_pPhiPrime.AddItem(NULL);
            }
        }
    }

public:

    TArray <class CFieldGauge*> m_pGaugeField;
    TArray <class CFieldGauge*> m_pForceField;
    TArray <class CFieldGauge*> m_pMomentumField;
    //TArray <class CFieldGauge*> m_pStapleField;

    TArray<class CFieldBoson*> m_pBosonFields;
    TArray<class CFieldBoson*> m_pBosonForceFields;
    TArray<class CFieldBoson*> m_pBosonMomentumFields;

    TArray<class CFieldTensor2*> m_pTensor2Field;

protected:

    /**
    * For force gradient updaters
    */
    TArray<class CFieldGauge*> m_pUPrime;
    TArray<class CFieldBoson*> m_pPhiPrime;
    EBackupFieldType m_eBackupFieldType;

    //tag for gauge smearing and cache, if U not changed, do not cache or smearing
    UBOOL m_bUDirty[2];
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
    void UpdatePG(Real fStep);

    void ChangeStepCount(UBOOL bGrow) override
    {
        CIntegrator::ChangeStepCount(bGrow);
        m_fNestedStepLength = m_fEStep / m_uiNestedStep;
    }

    void ChangeStepCountTo(UINT uiStep) override
    {
        CIntegrator::ChangeStepCountTo(uiStep);
        m_fNestedStepLength = m_fEStep / m_uiNestedStep;
    }

protected:

    void NestedEvaluateLeapfrog(UBOOL bLast);

    UINT m_uiNestedStep;
    Real m_fNestedStepLength;
    UBOOL m_bInnerLeapFrog;
};

/**
* The multi-level nested integrator seems broken
* No time for testing and fixing it
* Just do not use it plz
*/
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

    void ChangeStepCountTo(UINT uiStep) override
    {
        CIntegrator::ChangeStepCountTo(uiStep);
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
    void UpdateP(Real fStep, TArray<UINT> actionList, ESolverPhase ePhase, UBOOL bUpdateP);
    void UpdateP(Real fStep, INT iLevel, ESolverPhase ePhase, UBOOL bUpdateP)
    {
        UpdateP(fStep, m_iNestedActionId[iLevel], ePhase, bUpdateP);
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