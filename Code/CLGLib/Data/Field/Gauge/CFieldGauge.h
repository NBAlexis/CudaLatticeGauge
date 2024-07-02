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

    virtual void CalculateOnlyStaple(CFieldGauge* pStable) const = 0;

    virtual void MakeRandomGenerator() = 0;

    virtual DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const = 0;

    virtual DOUBLE CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const = 0;

    virtual DOUBLE CalculatePlaqutteEnergyUsingStable(DOUBLE betaOverN, const CFieldGauge* pStable) const = 0;

    virtual DOUBLE CalculateKinematicEnergy() const = 0;

    /**
    * U = exp(a this)U
    */
    virtual void ExpMult(Real a, CField* U) const = 0;

    /**
    * Due to numerical precision, sometimes, the Unitary matrix will diviate from Unitary a little bit.
    * This is to make the elements Unitary again.
    */
    virtual void ElementNormalize() = 0;

    /**
     * Use to simplify the implementation of gauge force
     */
    virtual void SetOneDirectionUnity(BYTE byDir) = 0;
    virtual void SetOneDirectionZero(BYTE byDir) = 0;

    virtual UINT MatrixN() const = 0;

#pragma endregion

#pragma region Test Functions to test gauge invarience of angular momentum

    /**
     * iA = U.TA() / 2
     */
    virtual void TransformToIA() = 0;

    /**
     * U=exp(iA)
     */
    virtual void TransformToU() = 0;

    /**
     * E_i(n) = U_{4i}(n) for i = 0,1,2
     */
    virtual void CalculateE_Using_U(CFieldGauge* pResoult) const = 0;

    /**
     * X_0(n) = nabla E
     */
    virtual void CalculateNablaE_Using_U(CFieldGauge* pResoult, UBOOL bNaive = FALSE) const = 0;

#pragma endregion


    UBOOL ApplyOperator(EFieldOperator , INT , INT , const CFieldGauge* const* , const CFieldBoson* const* , EOperatorCoefficientType , Real , Real , void* ) override
    {
        appCrucial("CFieldGauge: Do Operator implimented yet");
        return FALSE;
    }

    void CopyTo(CField* U) const override;

    UBOOL IsGaugeField() const override { return TRUE; }
    UBOOL IsFermionField() const override { return FALSE; }

protected:

    UINT m_uiLinkeCount;

};

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_H_

//=============================================================================
// END OF FILE
//=============================================================================