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
        , m_pDeviceXYBufferChiral(NULL)
        , m_pDeviceXYBufferPion(NULL)
        , m_pDeviceXYBufferRhon(NULL)
        , m_pHostXYBuffer(NULL)

        , m_pDistributionR(NULL)
        , m_pDistributionChiral(NULL)
        , m_pDistributionPion(NULL)
        , m_pDistributionRhon(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistributionChiral(NULL)
        , m_pHostDistributionPion(NULL)
        , m_pHostDistributionRhon(NULL)

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
    CLGComplex* m_pDeviceXYBufferChiral;
    CLGComplex* m_pDeviceXYBufferPion;
    CLGComplex* m_pDeviceXYBufferRhon;
    CLGComplex* m_pHostXYBuffer;
    CLGComplex m_cTmpSumChiral;
    CLGComplex m_cTmpSumPion;
    CLGComplex m_cTmpSumRhon;

    UINT* m_pDistributionR;
    CLGComplex* m_pDistributionChiral;
    CLGComplex* m_pDistributionPion;
    CLGComplex* m_pDistributionRhon;
    UINT* m_pHostDistributionR;
    CLGComplex* m_pHostDistributionChiral;
    CLGComplex* m_pHostDistributionPion;
    CLGComplex* m_pHostDistributionRhon;

    UINT m_uiMaxR;
    UBOOL m_bMeasureDistribution;
    UBOOL m_bShowResult;

public:

    TArray<CLGComplex> m_lstChiralAll;
    TArray<CLGComplex> m_lstPionAll;
    TArray<CLGComplex> m_lstRhonAll;
    TArray<CLGComplex> m_lstCondensateDensity;

    CLGComplex m_cAverageCondensate;
    TArray<CLGComplex> m_lstAverageCondensateDensity;

    //c(R)
    TArray<UINT> m_lstR;
    TArray<CLGComplex> m_lstChiral;
    TArray<CLGComplex> m_lstPion;
    TArray<CLGComplex> m_lstRhon;
};

__END_NAMESPACE

#endif //#ifndef _CMEASURECHIRALCONDENSATE_H_

//=============================================================================
// END OF FILE
//=============================================================================