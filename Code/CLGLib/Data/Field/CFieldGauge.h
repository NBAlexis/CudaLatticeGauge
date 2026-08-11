//=============================================================================
// FILENAME : CFieldGauge.h
// 
// DESCRIPTION:
// This is the class for the gauge fields
// Gauge fields are defined on links, so the total number of elements are N X Dir
//
// REVISION:
//  [mm/dd/yy]
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
    virtual void CalculateForceAndStaple(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const = 0;

    virtual void CalculateForceAnisotropy(CFieldGauge* pForce, DOUBLE betaOverN, DOUBLE xi) const
    {
        appCrucial(_T("CalculateForceAnisotropy not implemented for this case!\n"));
    }

    /**
    * For Dirichlet boundary condition, we need a different interface
    * For other cases, just use CalculateForceAndStaple
    */
    virtual void CalculateForceAndStapleClover(CFieldGauge* pForce, CFieldGauge* pStaple, Real betaOverN) const
    {
        CalculateForceAndStaple(pForce, pStaple, betaOverN);
    }
    virtual void CalculateForceCloverAnisotropy(CFieldGauge* pForce, DOUBLE betaOverN, DOUBLE xi) const
    {
        CalculateForceAnisotropy(pForce, betaOverN, xi);
    }

    virtual void CalculateOnlyStaple(CFieldGauge* pStaple) const = 0;

    virtual void MakeRandomGenerator() = 0;

    virtual DOUBLE CalculatePlaqutteEnergy(DOUBLE betaOverN) const = 0;

    virtual DOUBLE CalculatePlaqutteEnergyAnisotropy(DOUBLE betaOverN, DOUBLE xi) const
    {
        appCrucial(_T("CalculatePlaqutteEnergyAnisotropy not implemented for this case!\n"));
        return 0.0;
    }

    /**
    * CalculatePlaqutteEnergy accutally include the improved action
    * To measure for example u0, use this CalculatePlaqutteEnergyOriginal
    */
    virtual DOUBLE CalculatePlaqutteEnergyOriginal(DOUBLE betaOverN) const = 0;

    virtual DOUBLE CalculatePlaqutteEnergyUseClover(DOUBLE betaOverN) const = 0;

    virtual DOUBLE CalculatePlaqutteEnergyUseCloverAnisotropy(DOUBLE betaOverN, DOUBLE xi) const
    {
        appCrucial(_T("CalculatePlaqutteEnergyUseCloverAnisotropy not implemented for this case!\n"));
        return 0.0;
    }

    virtual DOUBLE CalculatePlaqutteEnergyUsingStaple(DOUBLE betaOverN, const CFieldGauge* pStaple) const = 0;

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

    virtual void TA() = 0;

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


    UBOOL ApplyOperator(EFieldOperator , INT , INT , INT , const CFieldGauge* const* , const CFieldBoson* const* , const CFieldTensor2* const* , EOperatorCoefficientType , Real , Real , void* ) override
    {
        appCrucial("CFieldGauge: Do Operator implimented yet\n");
        return FALSE;
    }

    /**
    * add for test
    * type = 0: x,y,z,t standard
    * type = 1: MILC convention (t,x,y,z)
    */
    virtual void ApplyStaggeredPhase(UINT iType = 0)
    {
        appCrucial("CFieldGauge: ApplyStaggeredPhase not implimented yet\n");
    }

    virtual void PolyakovOnSpatialSite(cuDoubleComplex* buffer, BYTE byDir = 3) const = 0;

    void CopyTo(CField* U) const override;

    UBOOL IsGaugeField() const override { return TRUE; }

    virtual void SetAsConnection(const CField* n, const CField* n_p_m)
    {
        appCrucial("CFieldGauge: SetAsConnection not implimented yet");
    }

    virtual void SetAsConnection(const CField* n, Real fCoeff)
    {
        appCrucial("CFieldGauge: SetAsConnection not implimented yet");
    }

    virtual void AddConnection(const CField* n, Real fCoeff)
    {
        appCrucial("CFieldGauge: AddConnection not implimented yet");
    }

    virtual void AddConnection(const CField* n, const CField* n_p_m, Real fCoeff)
    {
        appCrucial("CFieldGauge: AddConnection not implimented yet");
    }

    virtual void AddLinkTo(CField* target, const SCHAR* devicePath, BYTE byPathLen, BYTE mu, Real fCoeff) const
    {
        appCrucial("CFieldGauge: AddLinkTo not implimented yet");
    }

    virtual void AddNaikForce(const CFieldGauge* naikf0)
    {
        appCrucial("CFieldGauge: AddNaikForce not implimented yet");
    }

    UINT GetLinkCount() const { return m_uiLinkeCount; }

    void CopyParamTo(CField* U) const override
    {
        CField::CopyParamTo(U);
        CFieldGauge* pField = dynamic_cast<CFieldGauge*>(U);
        pField->m_uiLinkeCount = m_uiLinkeCount;
        pField->m_uiHaloLinkCount = m_uiHaloLinkCount;
    }

    //Multi-GPU (Phase 1): number of halo link slots appended after the local
    //links in m_pDeviceData. 0 on single-GPU / unsplit builds. See
    //Docs/MultiGPU-Plan.md section 10 and Core/Distributed/CLGHaloLayout.h.
    UINT GetHaloLinkCount() const { return m_uiHaloLinkCount; }

protected:

    UINT m_uiLinkeCount;

    //Halo capacity appended to m_pDeviceData (links). Set at allocation time.
    UINT m_uiHaloLinkCount;

};

extern void CLGAPI appSetGaugeLink(CFieldGauge* pGauge, UINT uiLinkIndex, const CLGComplex* pMatrix);

__END_NAMESPACE

#endif //#ifndef _CFIELDGAUGE_H_

//=============================================================================
// END OF FILE
//=============================================================================