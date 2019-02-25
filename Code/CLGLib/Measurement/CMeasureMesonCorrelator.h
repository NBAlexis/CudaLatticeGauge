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

class CLGAPI CMeasureMesonCorrelator : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureMesonCorrelator)
public:
    CMeasureMesonCorrelator() : CMeasure(), m_uiResoultCount(0) {}
    virtual void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId);

    virtual void OnConfigurationAccepted(const CFieldGauge* pGaugeField, const CFieldGauge* pStapleField);
    virtual void Average(UINT uiConfigurationCount);
    virtual void Report();
    virtual void Reset();

protected:

    void CalculateCorrelator(const CFieldGauge* pGaugeField, const CFieldGauge* pStapleField);

public:

    BYTE m_byFermionFieldId;
    UINT m_uiLt;
    Real m_f2OverVolumnSqrt;
    TArray<EGammaMatrix> m_lstGammas;

    //When a configuration is generated, result[nt] = (result[nt] * N + newresult[nt])/(N+1), N=N+1,
    TArray<TArray<Real>> m_lstResults;
    UINT m_uiResoultCount;
    //This is a complex field at each site
    _Complex * m_pDeviceCorrelator;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREMESONCORRELATOR_H_

//=============================================================================
// END OF FILE
//=============================================================================