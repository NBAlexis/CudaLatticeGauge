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

class CLGAPI CIntegrator
{
public:
    CIntegrator(class CHMC* pOwner);
    virtual void Evaluate() = 0;

    //std::vector<CAction*>

protected:

    UBOOL m_bPCached;
    UINT m_uiExpPrecision;
    FLOAT m_fEStep;
    UINT m_uiStepCount;

    CFieldGauge* m_pGaugeField;
    CFieldGauge* m_pForceField;
    CFieldGauge* m_pMomentumField;
    CFieldGauge* m_pStapleField;

    class CHMC* m_pOwner;
    CLatticeData* m_pLattice;
    CDeviceLattice* m_pDeviceLattice;
    TArray<class CAction*> m_lstActions;

    virtual void Initial(const TArray<class CAction*>& actionList, UINT uiStepCount, UINT uiExpPrecision);
    void Prepare();
    void UpdateP(FLOAT fStep);
    void UpdateU(FLOAT fStep);
};

__END_NAMESPACE

#endif //#ifndef _CINTEGRATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================