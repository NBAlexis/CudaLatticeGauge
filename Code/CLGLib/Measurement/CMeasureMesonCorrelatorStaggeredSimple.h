//=============================================================================
// FILENAME : CMeasureMesonCorrelatorStaggeredSimple.h
// 
// DESCRIPTION:
// This is the 10.1016/0550-3213(88)90396-3
// It use only ONE propagator, and not gauge fixing things
//
// REVISION:
//  [10/04/2020 nbale]
//=============================================================================

#ifndef _CMEASUREMESONCORRELATORSTAGGEREDSIMPLE_H_
#define _CMEASUREMESONCORRELATORSTAGGEREDSIMPLE_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureMesonCorrelatorStaggeredSimple)

class CLGAPI CMeasureMesonCorrelatorStaggeredSimple : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureMesonCorrelatorStaggeredSimple)
public:

    enum { _kMesonCorrelatorTypeSimple = 4 };

    CMeasureMesonCorrelatorStaggeredSimple() : CMeasure()
        , m_pDevicePropogators(NULL)
        , m_pResPropogators(NULL)
        , m_uiConfigurationCount(0)
        , m_bShowResult(FALSE)
    {
        
    }
    ~CMeasureMesonCorrelatorStaggeredSimple();
    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;

    void OnConfigurationAccepted(const CFieldGauge* pGaugeField, const CFieldGauge* pStapleField) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Average(UINT uiConfigurationCount) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    Real* m_pDevicePropogators;

#if !_CLG_DOUBLEFLOAT
    DOUBLE* m_pResPropogators;
#else
    //This is 4 x (Lt - 1)
    Real* m_pResPropogators;
#endif

public:

#if !_CLG_DOUBLEFLOAT
    TArray<TArray<DOUBLE>> m_lstAverageResults;

    //m_lstResults[conf][type][t]
    TArray<TArray<TArray<DOUBLE>>> m_lstResults;
#else
    TArray<TArray<Real>> m_lstAverageResults;

    //m_lstResults[conf][type][t]
    TArray<TArray<TArray<Real>>> m_lstResults;
#endif
    UINT m_uiConfigurationCount;
    UBOOL m_bShowResult;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREMESONCORRELATORSTAGGEREDSIMPLE_H_

//=============================================================================
// END OF FILE
//=============================================================================