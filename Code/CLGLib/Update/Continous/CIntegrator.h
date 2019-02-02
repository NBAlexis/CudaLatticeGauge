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
        , m_bPCached(FALSE)
    { ; }

    ~CIntegrator();

    virtual void Evaluate() = 0;
    virtual void Initial(class CHMC* pOwner, class CLatticeData* pLattice, const CParameters& params);
    void Accept();
    void Prepare(UBOOL bLastAccepted);
    void UpdateP(Real fStep);
    void UpdateU(Real fStep);
    virtual Real GetEnergy() const;

protected:

    UBOOL m_bPCached;
    Real m_fEStep;
    UINT m_uiStepCount;

    Real m_fUpdateResultEnery;

    CFieldGauge* m_pGaugeField;
    CFieldGauge* m_pForceField;
    CFieldGauge* m_pMomentumField;
    CFieldGauge* m_pStapleField;

    class CHMC* m_pOwner;
    CLatticeData* m_pLattice;
    TArray<class CAction*> m_lstActions;
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================