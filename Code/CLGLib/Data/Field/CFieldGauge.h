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

#pragma region HMC update

    /**
    * Before many other steps
    */
    virtual void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStable, const _Complex& minusBetaOverN) const = 0;

    virtual void MakeRandomGenerator() = 0;

    virtual Real CalculatePlaqutteEnergy(const _Complex& minusBetaOverN) const = 0;

    virtual Real CalculateKinematicEnergy() const = 0;

#pragma endregion HMC update

    virtual _Complex Dot(const CField* other) const
    {
        UN_USE(other);
        appCrucial("CFieldGauge: Do NOT know how to DOT SU3");
        return _make_cuComplex(0, 0);
    }

    virtual void ApplyOperator(EFieldOperator op, const CField* otherfield)
    {
        UN_USE(op);
        UN_USE(otherfield);
        appCrucial("CFieldGauge: Do Operator implimented yet");
    }

protected:

    UINT m_uiLinkeCount;
};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_H_

//=============================================================================
// END OF FILE
//=============================================================================