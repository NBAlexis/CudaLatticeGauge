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

#define OneOver32PI2 (F(1.0) / (F(32.0) * PISQ))

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureChiralCondensate)

class CLGAPI CMeasureChiralCondensate : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureChiralCondensate)
public:
    CMeasureChiralCondensate()
        : CMeasure()
        , m_uiFieldCount(100)
        , m_uiConfigurationCount(0)
        , m_pDeviceXYBuffer(NULL)
        , m_pHostXYBuffer(NULL)
    {
        
    }

    ~CMeasureChiralCondensate();

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
    CLGComplex* m_pDeviceXYBuffer;
    CLGComplex* m_pHostXYBuffer;

public:

    TArray<CLGComplex> m_lstCondensate;
    TArray<CLGComplex> m_lstCondensateDensity;
};

__END_NAMESPACE

#endif //#ifndef _CMEASURECHIRALCONDENSATE_H_

//=============================================================================
// END OF FILE
//=============================================================================