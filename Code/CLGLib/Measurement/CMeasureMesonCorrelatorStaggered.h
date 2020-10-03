//=============================================================================
// FILENAME : CMeasureMesonCorrelatorStaggered.h
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [02/22/2019 nbale]
//=============================================================================

#ifndef _CMEASUREMESONCORRELATORSTAGGERED_H_
#define _CMEASUREMESONCORRELATORSTAGGERED_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureMesonCorrelatorStaggered)

class CLGAPI CMeasureMesonCorrelatorStaggered : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureMesonCorrelatorStaggered)
public:

    enum { _kMesonCorrelatorType = 20 };

    CMeasureMesonCorrelatorStaggered() : CMeasure()
        , m_pDeviceW1(NULL)
        , m_pDeviceW2(NULL)
        , m_pGaugeFixing(NULL)
        , m_pDeviceSignTable(NULL)
        , m_pDeviceDeltaTable(NULL)
        , m_pDevicePropogators(NULL)
        , m_pDevicePropogatorsEveryTimeSlice(NULL)
        , m_pResPropogators(NULL)
        , m_uiConfigurationCount(0)
        , m_bGaugeFixing(FALSE)
        , m_bShowResult(FALSE)
    {
        
    }
    ~CMeasureMesonCorrelatorStaggered();
    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;

    void OnConfigurationAccepted(const CFieldGauge* pGaugeField, const CFieldGauge* pStapleField) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Average(UINT uiConfigurationCount) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    //This is a 20 x 8 table
    BYTE m_pSignTable[_kMesonCorrelatorType * 8];
    //This is a 20 x 1 table
    BYTE m_pDeltaTable[_kMesonCorrelatorType];
    void InitialSignTable();

    CFieldFermionKSSU3* m_pW1[24];
    CFieldFermionKSSU3* m_pW2[24];
    deviceSU3Vector** m_pDeviceW1;
    deviceSU3Vector** m_pDeviceW2;
    void CalculateSources(const CFieldGauge* pGauge);

    //uiSite
    CFieldGauge* m_pGaugeFixing;
    void CalculatePropogators();

    BYTE* m_pDeviceSignTable;
    BYTE* m_pDeviceDeltaTable;
    CLGComplex* m_pDevicePropogators;
    CLGComplex* m_pDevicePropogatorsEveryTimeSlice;

    //This is 20 x (Lt - 2)
    CLGComplex* m_pResPropogators;

    void InitialBuffers();

public:

    TArray<TArray<Real>> m_lstAverageResults;

    //m_lstResults[conf][type][t]
    TArray<TArray<TArray<CLGComplex>>> m_lstResults;
    UINT m_uiConfigurationCount;
    UBOOL m_bGaugeFixing;
    UBOOL m_bShowResult;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREMESONCORRELATORSTAGGERED_H_

//=============================================================================
// END OF FILE
//=============================================================================