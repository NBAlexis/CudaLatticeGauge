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

    enum { _kCondMeasureCount = 9 };

    CMeasureChiralCondensate()
        : CMeasureStochastic()
        , m_pHostXYBuffer(NULL)

        , m_pDistributionR(NULL)
        , m_pDistribution(NULL)
        , m_pHostDistributionR(NULL)
        , m_pHostDistribution(NULL)

        , m_uiMaxR(1)
        , m_bMeasureDistribution(FALSE)
    {
        
    }

    ~CMeasureChiralCondensate();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedZ4SingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return FALSE; }
    UBOOL IsZ4Source() const override { return TRUE; }

protected:
    
    CLGComplex* m_pDeviceXYBuffer[_kCondMeasureCount];
    CLGComplex* m_pHostXYBuffer;
    CLGComplex m_cTmpSum[_kCondMeasureCount];

    UINT* m_pDistributionR;
    CLGComplex* m_pDistribution;
    UINT* m_pHostDistributionR;
    CLGComplex* m_pHostDistribution;
    UINT m_uiMaxR;
    UBOOL m_bMeasureDistribution;

public:

    TArray<UINT> m_lstR;
    TArray<CLGComplex> m_lstCondAll[_kCondMeasureCount];
    TArray<CLGComplex> m_lstCond[_kCondMeasureCount];  
};

__END_NAMESPACE

#endif //#ifndef _CMEASURECHIRALCONDENSATE_H_

//=============================================================================
// END OF FILE
//=============================================================================