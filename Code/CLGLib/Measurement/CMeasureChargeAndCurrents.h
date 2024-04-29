//=============================================================================
// FILENAME : CMeasureChargeAndCurrents.h
// 
// DESCRIPTION:
// This is measurement for some charge and current densities
//
//  -- Discarded
//
// REVISION:
//  [05/28/2019 nbale]
//=============================================================================

#ifndef _CMEASURECHARGEANDCURRENT_H_
#define _CMEASURECHARGEANDCURRENT_H_

__BEGIN_NAMESPACE

typedef void(*_deviceMeasureFunc)(const SSmallInt4& site, Real fKappa, Real fOmega, deviceWilsonVectorSU3& element);

__CLG_REGISTER_HELPER_HEADER(CMeasureChargeAndCurrents)

class CLGAPI CMeasureChargeAndCurrents : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureChargeAndCurrents)

public:

    /**
    * We want to know
    * just D^-1
    * gamma_1,2,3,4
    * gamma_4 sigma_12
    * gamma_1 -i y Omega gamma_4
    * gamma_2 +i x Omega gamma_4
    * gamma_4 gamma _5
    */
    enum { _kGammaInInterests = 9, };

    CMeasureChargeAndCurrents();
    ~CMeasureChargeAndCurrents();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override {}
    void SourceSanningSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return FALSE; }
    UBOOL IsSourceScanning() const override { return TRUE; }

protected:

    CLGComplex * m_pHostDataBuffer;
    CLGComplex * m_pDeviceDataBuffer;
    deviceWilsonVectorSU3* m_pOperatorData;
    
    //_deviceMeasureFunc * m_pMeasureFunctions;

public:

    TArray<CLGComplex> m_lstAllRes;
};

__END_NAMESPACE

#endif //#ifndef _CMEASURECHARGEANDCURRENT_H_

//=============================================================================
// END OF FILE
//=============================================================================