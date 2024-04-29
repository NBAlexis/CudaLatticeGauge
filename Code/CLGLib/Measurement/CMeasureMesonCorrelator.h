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
    CMeasureMesonCorrelator() : CMeasure() {}
    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;

    void OnConfigurationAcceptedSingleField(const CFieldGauge* pGaugeField, const CFieldGauge* pStapleField) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    void CalculateCorrelator(const CFieldGauge* pGaugeField, const CFieldGauge* pStapleField);

public:

    UINT m_uiLt;
#if !_CLG_DOUBLEFLOAT
    DOUBLE m_f2OverVolumnSqrt;
#else
    Real m_f2OverVolumnSqrt;
#endif
    TArray<EGammaMatrix> m_lstGammas;

    //When a configuration is generated, result[nt] = (result[nt] * N + newresult[nt])/(N+1), N=N+1,
#if !_CLG_DOUBLEFLOAT
    TArray<TArray<DOUBLE>> m_lstResults;
    TArray<TArray<DOUBLE>> m_lstResultsLastConf;
#else
    TArray<TArray<Real>> m_lstResults;
    TArray<TArray<Real>> m_lstResultsLastConf;
#endif
    //This is a complex field at each site
    CLGComplex * m_pDeviceCorrelator;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREMESONCORRELATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================