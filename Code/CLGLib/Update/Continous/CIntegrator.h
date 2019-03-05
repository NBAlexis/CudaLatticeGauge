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
        : m_uiStepCount(1)
        , m_fEStep(0)
        , m_pOwner(NULL)
        , m_pLattice(NULL)
        , m_pGaugeField(NULL)
        , m_pForceField(NULL)
        , m_pMomentumField(NULL)
        , m_pStapleField(NULL)
        , m_bStapleCached(FALSE)
    { ; }

    ~CIntegrator();

    virtual void Evaluate() = 0;
    virtual void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params);
    void OnFinishTrajectory(UBOOL bAccepted);
    void Prepare(UBOOL bLastAccepted, UINT uiStep);
    void UpdateP(Real fStep, UBOOL bCacheStaple);
    void UpdateU(Real fStep);
    virtual Real GetEnergy(UBOOL bBeforeEvolution) const;
    virtual CCString GetInfos(const CCString& sTab) const = 0;
    void FinishEvaluate();

protected:

    UBOOL m_bStapleCached;
    Real m_fEStep;
    UINT m_uiStepCount;

    Real m_fUpdateResultEnery;

public:

    CFieldGauge* m_pGaugeField;
    CFieldGauge* m_pForceField;
    CFieldGauge* m_pMomentumField;
    CFieldGauge* m_pStapleField;

protected:

    class CHMC* m_pOwner;
    CLatticeData* m_pLattice;
    TArray<class CAction*> m_lstActions;
};

class CLGAPI CNestedIntegrator : public CIntegrator
{
public:
    CNestedIntegrator()
        : CIntegrator()
        , m_uiNestedStep(1)
        , m_fNestedStepLength(F(0.0))
    {
    }

    virtual void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params);
    CCString GetNestedInfo(const CCString & sTab) const;

    void UpdatePF(Real fStep);
    void UpdatePG(Real fStep, UBOOL bCacheStaple);

protected:

    UINT m_uiNestedStep;
    Real m_fNestedStepLength;
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================