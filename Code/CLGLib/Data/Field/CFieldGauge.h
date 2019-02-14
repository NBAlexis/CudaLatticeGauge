//=============================================================================
// FILENAME : CFieldGauge.h
// 
// DESCRIPTION:
// This is the class for the gauge fields
// Gauge fields are defined on links, so the total number of elements are N X Dir
//
// REVISION:
//  [12/3/2018 nbale]
//=============================================================================

#ifndef _CFIELDGAUGE_H_
#define _CFIELDGAUGE_H_

__BEGIN_NAMESPACE

class CLGAPI CFieldGauge : public CField
{
public:
    CFieldGauge();
    ~CFieldGauge();

#pragma region HMC update

    /**
    * Before many other steps
    */
    virtual void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, Real betaOverN) const = 0;

    virtual void MakeRandomGenerator() = 0;

    virtual Real CalculatePlaqutteEnergy(Real betaOverN) const = 0;

    virtual Real CalculatePlaqutteEnergyUsingStable(Real betaOverN, const CFieldGauge *pStable) const = 0;

    virtual Real CalculateKinematicEnergy() const = 0;

    /**
    * U = exp(a this)U
    */
    virtual void ExpMult(const _Complex& a, CField* U) const = 0;
    virtual void ExpMult(Real a, CField* U) const = 0;

    /**
    * Due to numerical precision, sometimes, the Unitary matrix will diviate from Unitary a little bit.
    * This is to make the elements Unitary again.
    */
    virtual void ElementNormalize() = 0;

    /**
    * Optimize
    */
    void CachePlaqutteIndexes();

#pragma endregion

    virtual UBOOL ApplyOperator(EFieldOperator , const CField*, EOperatorCoefficientType , Real , Real )
    {
        appCrucial("CFieldGauge: Do Operator implimented yet");
        return FALSE;
    }

    virtual void CopyTo(CField* U) const;

protected:

    UINT m_uiLinkeCount;

    UBOOL m_bPlaqutteIndexCached;
    UINT m_uiPlaqutteLength;
    UINT m_uiPlaqutteCountPerSite;
    UINT m_uiPlaqutteCountPerLink;

    SIndex* m_pPlaquttesPerSite;
    SIndex* m_pPlaquttesPerLink;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_H_

//=============================================================================
// END OF FILE
//=============================================================================