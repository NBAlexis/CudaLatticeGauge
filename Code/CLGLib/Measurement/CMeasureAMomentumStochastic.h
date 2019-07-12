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

class CLGAPI CMeasureAMomentumStochastic : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureAMomentumStochastic)
public:
    CMeasureAMomentumStochastic()
        : CMeasure()
        , m_uiFieldCount(100)
        , m_uiConfigurationCount(0)
        , m_pDeviceXYBufferJL(NULL)
        , m_pDeviceXYBufferJS(NULL)
        , m_bExponential(TRUE)
        , m_bNaive(TRUE)

        , m_pDistributionR(NULL)
        , m_pDistributionJL(NULL)
        , m_pDistributionJS(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistributionJL(NULL)
        , m_pHostDistributionJS(NULL)

        , m_uiMaxR(1)
        , m_uiEdgeR(1)
        , m_bShowResult(FALSE)
    {
        
    }

    ~CMeasureAMomentumStochastic();

    virtual void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId);
    virtual void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple);
    virtual void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) {}
    virtual void Average(UINT uiConfigurationCount);
    virtual void Report();
    virtual void Reset();

    virtual UBOOL IsGaugeMeasurement() const { return TRUE; }
    virtual UBOOL IsSourceScanning() const { return FALSE; }

protected:
    
    UINT m_uiFieldCount;
    UINT m_uiConfigurationCount;
    Real* m_pDeviceXYBufferJL;
    Real* m_pDeviceXYBufferJS;

    UBOOL m_bExponential;
    UBOOL m_bNaive;

    UINT* m_pDistributionR;
    Real* m_pDistributionJL;
    Real* m_pDistributionJS;
    UINT* m_pHostDistributionR;
    Real* m_pHostDistributionJL;
    Real* m_pHostDistributionJS;

    UINT m_uiMaxR;
    UINT m_uiEdgeR;
    UBOOL m_bShowResult;

public:

    TArray<UINT> m_lstR;
    TArray<Real> m_lstJL;
    TArray<Real> m_lstJS;

    TArray<Real> m_lstJLAll;
    TArray<Real> m_lstJLInner;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREAMOMENTUMSTOCHASTIC_H_

//=============================================================================
// END OF FILE
//=============================================================================