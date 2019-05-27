//=============================================================================
// FILENAME : CMeasureChargeAndCurrents.h
// 
// DESCRIPTION:
// This is measurement for some charge and current densities
//
// REVISION:
//  [05/28/2019 nbale]
//=============================================================================

#ifndef _CMEASURECHARGEANDCURRENT_H_
#define _CMEASURECHARGEANDCURRENT_H_

__BEGIN_NAMESPACE

typedef void(*_deviceMeasureFunc)(const SSmallInt4& site, const SSmallInt4& sCenter, Real fKappa, Real fOmega, deviceWilsonVectorSU3& element);

__CLG_REGISTER_HELPER_HEADER(CMeasureChargeAndCurrents)

class CLGAPI CMeasureChargeAndCurrents : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureChargeAndCurrents)

public:

    /**
    * We want to know
    * gamma_1,2,3,4
    * gamma_4 sigma_12
    * gamma_1 -i y Omega gamma_4
    * gamma_2 +i x Omega gamma_4
    * gamma_4 gamma _5
    */
    enum { _kGammaInInterests = 8, };

    CMeasureChargeAndCurrents();
    ~CMeasureChargeAndCurrents();

    virtual void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId);
    virtual void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) {}
    virtual void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site);
    virtual void Average(UINT uiConfigurationCount);
    virtual void Report();
    virtual void Reset();

    virtual UBOOL IsGaugeMeasurement() const { return FALSE; }
    virtual UBOOL IsSourceScanning() const { return TRUE; }

protected:

    CLGComplex * m_pHostDataBuffer;
    CLGComplex * m_pDeviceDataBuffer;
    deviceWilsonVectorSU3* m_pOperatorData;
    
    UINT m_uiConfigurationCount;
    UBOOL m_bShowResult;
    //_deviceMeasureFunc * m_pMeasureFunctions;

    TArray<CLGComplex> m_lstAllRes;
};

__END_NAMESPACE

#endif //#ifndef _CMEASURECHARGEANDCURRENT_H_

//=============================================================================
// END OF FILE
//=============================================================================