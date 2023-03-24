//=============================================================================
// FILENAME : CMeasureAngularMomentumKS.h
// 
// DESCRIPTION:
// NOTE: 
// 
// JL = qbar GA4 (y Dx - x Dy)q
// JS = qbar (i/2) GA4 S12E q
// 
//
// REVISION:
//  [01/17/2021 nbale]
//=============================================================================

#ifndef _CMEASUREANGULARMOMENTUMKS_H_
#define _CMEASUREANGULARMOMENTUMKS_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureAngularMomentumKS)

enum EAngularMeasureTypeKS
{
    OrbitalKS = 0,
    SpinKS = 1,
    PotentialKS = 2,

    EAngularMeasureMax,
};

class CLGAPI CMeasureAngularMomentumKS : public CMeasureStochastic
{
    __CLGDECLARE_CLASS(CMeasureAngularMomentumKS)
public:

    CMeasureAngularMomentumKS()
        : CMeasureStochastic()
        , m_uiConfigurationCount(0)
        , m_pHostXYBuffer(NULL)
        , m_pHostZBuffer(NULL)

        , m_pDistributionR(NULL)
        , m_pDistribution(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistribution(NULL)

        , m_uiMaxR(1)
        , m_uiEdge(1)
        , m_bShiftCenter(FALSE)
        , m_bShowResult(FALSE)
        , m_bMeasureZSlice(FALSE)
    {
        
    }

    ~CMeasureAngularMomentumKS();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedZ4(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd) override;
    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Average(UINT uiConfigurationCount) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return FALSE; }
    UBOOL IsZ4Source() const override { return TRUE; }

protected:

    virtual void ApplyOrbitalMatrix(deviceSU3Vector* pAppliedBuffer, const deviceSU3Vector* pInverseZ4, const deviceSU3* pGauge) const;
    virtual void ApplySpinMatrix(deviceSU3Vector* pAppliedBuffer, const deviceSU3Vector* pInverseZ4, const deviceSU3* pGauge) const;
    virtual void ApplyPotentialMatrix(deviceSU3Vector* pAppliedBuffer, const deviceSU3Vector* pInverseZ4, const deviceSU3* pGauge) const;
    
    UINT m_uiConfigurationCount;

    CLGComplex* m_pDeviceXYBuffer[EAngularMeasureMax];
    CLGComplex* m_pDeviceZBuffer[EAngularMeasureMax];
    CLGComplex* m_pHostXYBuffer;
    CLGComplex* m_pHostZBuffer;
    //CLGComplex m_cTmpSum[EAngularMeasureMax];

    UINT* m_pDistributionR;
    CLGComplex* m_pDistribution;
    UINT* m_pHostDistributionR;
    CLGComplex* m_pHostDistribution;
    UINT m_uiMaxR;
    UINT m_uiEdge;
    UBOOL m_bShiftCenter;
    UBOOL m_bShowResult;

public:

    UBOOL m_bMeasureZSlice;
    TArray<UINT> m_lstR;
    TArray<CLGComplex> m_lstCondAll[EAngularMeasureMax];
    TArray<CLGComplex> m_lstCondIn[EAngularMeasureMax];
    TArray<CLGComplex> m_lstCond[EAngularMeasureMax];
    TArray<CLGComplex> m_lstCondZSlice[EAngularMeasureMax];
};

__END_NAMESPACE

#endif //#ifndef _CMEASURECHIRALCONDENSATEKS_H_

//=============================================================================
// END OF FILE
//=============================================================================