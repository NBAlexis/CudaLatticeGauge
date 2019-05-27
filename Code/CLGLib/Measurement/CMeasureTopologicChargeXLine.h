//=============================================================================
// FILENAME : CMeasureTopologicChargeXLine.h
// 
// DESCRIPTION:
//
// REVISION:
//  [05/28/2019 nbale]
//=============================================================================

#ifndef _CMEASURETOPOLOGICALCHARGEXLINE_H_
#define _CMEASURETOPOLOGICALCHARGEXLINE_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureTopologicChargeXLine)

class CLGAPI CMeasureTopologicChargeXLine : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureTopologicChargeXLine)
public:
    CMeasureTopologicChargeXLine()
        : CMeasure()
        , m_pHostDataBuffer(NULL)
        , m_pDeviceDataBuffer(NULL) 
        , m_pOperatorData(NULL)
        , m_uiConfigurationCount(0)
        , m_bShowResult(TRUE)
    {
    }

    ~CMeasureTopologicChargeXLine();

    virtual void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId);
    virtual void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple);
    virtual void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) {}
    virtual void Average(UINT uiConfigurationCount);
    virtual void Report();
    virtual void Reset();

    virtual UBOOL IsGaugeMeasurement() const { return TRUE; }
    virtual UBOOL IsSourceScanning() const { return FALSE; }

protected:

    CLGComplex * m_pHostDataBuffer;
    CLGComplex * m_pDeviceDataBuffer;
    deviceWilsonVectorSU3* m_pOperatorData;
    
    UINT m_uiConfigurationCount;
    UBOOL m_bShowResult;
    TArray<Real> m_lstAllRes;
};

__END_NAMESPACE

#endif //#ifndef _CMEASURETOPOLOGICALCHARGEXLINE_H_

//=============================================================================
// END OF FILE
//=============================================================================