//=============================================================================
// FILENAME : CMeasureMesonCorrelator.h
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [02/22/2019 nbale]
//=============================================================================

#ifndef _CMEASUREMESONCORRELATOR_H_
#define _CMEASUREMESONCORRELATOR_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureMesonCorrelator)

class CLGAPI CMeasureMesonCorrelator : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureMesonCorrelator)
public:
    CMeasureMesonCorrelator() : CMeasure(), m_uiResoultCount(0) {}
    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;

    void OnConfigurationAccepted(const CFieldGauge* pGaugeField, const CFieldGauge* pStapleField) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Average(UINT uiConfigurationCount) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    void CalculateCorrelator(const CFieldGauge* pGaugeField, const CFieldGauge* pStapleField);

public:

    UINT m_uiLt;
    Real m_f2OverVolumnSqrt;
    TArray<EGammaMatrix> m_lstGammas;

    //When a configuration is generated, result[nt] = (result[nt] * N + newresult[nt])/(N+1), N=N+1,
    TArray<TArray<Real>> m_lstResults;
    TArray<TArray<Real>> m_lstResultsLastConf;
    UINT m_uiResoultCount;
    //This is a complex field at each site
    CLGComplex * m_pDeviceCorrelator;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREMESONCORRELATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================