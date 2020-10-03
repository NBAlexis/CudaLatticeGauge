//=============================================================================
// FILENAME : CMeasureChiralCondensateKS.h
// 
// DESCRIPTION:
// NOTE: 
//
// REVISION:
//  [10/01/2020 nbale]
//=============================================================================

#ifndef _CMEASURECHIRALCONDENSATEKS_H_
#define _CMEASURECHIRALCONDENSATEKS_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureChiralCondensateKS)

class CLGAPI CMeasureChiralCondensateKS : public CMeasureStochastic
{
    __CLGDECLARE_CLASS(CMeasureChiralCondensateKS)
public:

    //Only Chiral condensate, others will be added later
    enum { _kCondMeasureCountKS = 1 };

    CMeasureChiralCondensateKS()
        : CMeasureStochastic()
        , m_uiConfigurationCount(0)
        , m_pHostXYBuffer(NULL)

        , m_pDistributionR(NULL)
        , m_pDistribution(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistribution(NULL)

        , m_uiMaxR(1)
        , m_bMeasureDistribution(FALSE)
        , m_bShowResult(FALSE)
    {
        
    }

    ~CMeasureChiralCondensateKS();

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

    CLGComplex* m_pDeviceXYBuffer[_kCondMeasureCountKS];
    CLGComplex* m_pHostXYBuffer;
    CLGComplex m_cTmpSum[_kCondMeasureCountKS];

    UINT* m_pDistributionR;
    CLGComplex* m_pDistribution;
    UINT* m_pHostDistributionR;
    CLGComplex* m_pHostDistribution;
    UINT m_uiMaxR;
    UBOOL m_bMeasureDistribution;
    UBOOL m_bShowResult;

public:

    TArray<UINT> m_lstR;
    TArray<CLGComplex> m_lstCondAll[_kCondMeasureCountKS];
    TArray<CLGComplex> m_lstCond[_kCondMeasureCountKS];
};

__END_NAMESPACE

#endif //#ifndef _CMEASURECHIRALCONDENSATEKS_H_

//=============================================================================
// END OF FILE
//=============================================================================