//=============================================================================
// FILENAME : CMeasureChiralCondensate.h
// 
// DESCRIPTION:
// NOTE: 
//
// REVISION:
//  [05/28/2019 nbale]
//=============================================================================

#ifndef _CMEASURECHIRALCONDENSATE_H_
#define _CMEASURECHIRALCONDENSATE_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureChiralCondensate)

class CLGAPI CMeasureChiralCondensate : public CMeasureStochastic
{
    __CLGDECLARE_CLASS(CMeasureChiralCondensate)
public:
    CMeasureChiralCondensate()
        : CMeasureStochastic()
        , m_uiConfigurationCount(0)
        , m_pDeviceXYBuffer(NULL)
        , m_pHostXYBuffer(NULL)

        , m_pDistributionR(NULL)
        , m_pDistributionC(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistributionC(NULL)

        , m_uiMaxR(1)
        , m_bMeasureDistribution(FALSE)
        , m_bShowResult(FALSE)
    {
        
    }

    ~CMeasureChiralCondensate();

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
    
    UINT m_uiConfigurationCount;
    CLGComplex* m_pDeviceXYBuffer;
    CLGComplex* m_pHostXYBuffer;
    CLGComplex m_cTmpSum;

    UINT* m_pDistributionR;
    Real* m_pDistributionC;
    UINT* m_pHostDistributionR;
    Real* m_pHostDistributionC;

    UINT m_uiMaxR;
    UBOOL m_bMeasureDistribution;
    UBOOL m_bShowResult;

public:

    TArray<CLGComplex> m_lstCondensate;
    TArray<CLGComplex> m_lstCondensateDensity;
    CLGComplex m_cAverageCondensate;
    TArray<CLGComplex> m_lstAverageCondensateDensity;

    //c(R)
    TArray<UINT> m_lstR;
    TArray<Real> m_lstC;
};

__END_NAMESPACE

#endif //#ifndef _CMEASURECHIRALCONDENSATE_H_

//=============================================================================
// END OF FILE
//=============================================================================