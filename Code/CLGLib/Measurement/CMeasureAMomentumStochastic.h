//=============================================================================
// FILENAME : CMeasureAMomentumStochastic.h
// 
// DESCRIPTION:
// NOTE: 
//
// REVISION:
//  [07/10/2019 nbale]
//=============================================================================

#ifndef _CMEASUREAMOMENTUMSTOCHASTIC_H_
#define _CMEASUREAMOMENTUMSTOCHASTIC_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureAMomentumStochastic)

class CLGAPI CMeasureAMomentumStochastic : public CMeasureStochastic
{
    __CLGDECLARE_CLASS(CMeasureAMomentumStochastic)
public:
    CMeasureAMomentumStochastic()
        : CMeasureStochastic()
        , m_uiConfigurationCount(0)
        , m_pDeviceXYBufferJL(NULL)
        , m_pDeviceXYBufferJS(NULL)
        , m_pDeviceXYBufferJLPure(NULL)
        , m_pDeviceXYBufferJLJM(NULL)
        , m_pDeviceXYBufferJPot(NULL)
        , m_bExponential(TRUE)
        , m_bNaive(TRUE)
        , m_bMeasureJLPure(FALSE)

        , m_pDistributionR(NULL)
        , m_pDistributionJL(NULL)
        , m_pDistributionJS(NULL)
        , m_pDistributionJLPure(NULL)
        , m_pDistributionJLJM(NULL)
        , m_pDistributionJPot(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistributionJL(NULL)
        , m_pHostDistributionJS(NULL)
        , m_pHostDistributionJLPure(NULL)
        , m_pHostDistributionJLJM(NULL)
        , m_pHostDistributionJPot(NULL)

        , m_uiMaxR(1)
        , m_uiEdgeR(1)
        , m_bShowResult(FALSE)
    {
        
    }

    ~CMeasureAMomentumStochastic();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedZ4(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd) override;
    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Average(UINT uiConfigurationCount) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return FALSE; }
    UBOOL IsZ4Source() const override  { return TRUE; }

protected:
    
    UINT m_uiConfigurationCount;
    Real* m_pDeviceXYBufferJL;
    Real* m_pDeviceXYBufferJS;
    Real* m_pDeviceXYBufferJLPure;
    Real* m_pDeviceXYBufferJLJM;
    Real* m_pDeviceXYBufferJPot;

    UBOOL m_bExponential;
    UBOOL m_bNaive;
    UBOOL m_bMeasureJLPure;

    UINT* m_pDistributionR;
    Real* m_pDistributionJL;
    Real* m_pDistributionJS;
    Real* m_pDistributionJLPure;
    Real* m_pDistributionJLJM;
    Real* m_pDistributionJPot;
    UINT* m_pHostDistributionR;
    Real* m_pHostDistributionJL;
    Real* m_pHostDistributionJS;
    Real* m_pHostDistributionJLPure;
    Real* m_pHostDistributionJLJM;
    Real* m_pHostDistributionJPot;

    UINT m_uiMaxR;
    UINT m_uiEdgeR;
    UBOOL m_bShowResult;

public:

    TArray<UINT> m_lstR;
    TArray<Real> m_lstJL;
    TArray<Real> m_lstJS;
    TArray<Real> m_lstJLPure;
    TArray<Real> m_lstJLJM;
    TArray<Real> m_lstJPot;

    TArray<Real> m_lstJLAll;
    TArray<Real> m_lstJLInner;
    TArray<Real> m_lstJSAll;
    TArray<Real> m_lstJSInner;
    TArray<Real> m_lstJLPureAll;
    TArray<Real> m_lstJLPureInner;
    TArray<Real> m_lstJLJMAll;
    TArray<Real> m_lstJLJMInner;
    TArray<Real> m_lstJPotAll;
    TArray<Real> m_lstJPotInner;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREAMOMENTUMSTOCHASTIC_H_

//=============================================================================
// END OF FILE
//=============================================================================