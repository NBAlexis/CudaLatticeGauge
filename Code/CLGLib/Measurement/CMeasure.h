//=============================================================================
// FILENAME : CMeasure.h
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [01/29/2019 nbale]
//=============================================================================

#ifndef _CMEASURE_H_
#define _CMEASURE_H_

__BEGIN_NAMESPACE

class CLGAPI CMeasure : public CBase
{
public:
    CMeasure() 
        : m_pOwner(NULL)
        , m_bNeedSmearing(FALSE)
        , m_byFieldId(0)
        , m_fLastRealResult(F(0.0))
    {
    }

    virtual void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
    {
        m_pOwner = pOwner;
        m_pLatticeData = pLatticeData;
        m_byId = byId;

        INT iNeedGaugeSmearing = 0;
        param.FetchValueINT(_T("GaugeSmearing"), iNeedGaugeSmearing);
        m_bNeedSmearing = 0 != iNeedGaugeSmearing;
    }

    /**
    * Accept gauge can be smoothed.
    * pCorrespondingStaple Might be NULL.
    */
    virtual void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) = 0;

    /**
    * NOTE: sources will be passed to multiple measures, do NOT change the content!
    */
    virtual void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) = 0;
    virtual void Average(UINT uiConfigurationCount) = 0;
    virtual void Report() = 0;
    virtual void Reset() = 0;

    virtual UBOOL IsGaugeMeasurement() const = 0;
    virtual UBOOL IsSourceScanning() const = 0;
    virtual UBOOL NeedGaugeSmearing() const { return m_bNeedSmearing; }

    BYTE GetFieldId() const { return m_byFieldId; }

protected:

    class CMeasurementManager* m_pOwner;
    class CLatticeData* m_pLatticeData;
    UBOOL m_bNeedSmearing;
    BYTE m_byId;
    BYTE m_byFieldId;

public:
    //============================================================
    //some simple measurement only produce real or complex results
    Real m_fLastRealResult;
    CLGComplex m_cLastComplexResult;
    
};

__END_NAMESPACE

#endif //#ifndef _CMEASURE_H_

//=============================================================================
// END OF FILE
//=============================================================================