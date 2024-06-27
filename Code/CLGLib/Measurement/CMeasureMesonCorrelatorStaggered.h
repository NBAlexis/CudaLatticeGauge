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
        , m_bSimpleVersion(FALSE)
        , m_pDeviceW1(NULL)
        , m_pDeviceW2(NULL)
        , m_pDeviceSignTable(NULL)
        , m_pDeviceDeltaTable(NULL)
        , m_pDevicePropogators(NULL)
        , m_pDevicePropogatorsEveryTimeSlice(NULL)
        , m_pResPropogators(NULL)
        , m_bGaugeFixing(FALSE)
    {
        
    }
    ~CMeasureMesonCorrelatorStaggered();
    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;

    void OnConfigurationAccepted(INT gn, INT bn, const CFieldGauge* const* gs, const CFieldBoson* const* bs, const CFieldGauge* const* stp) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    UBOOL m_bSimpleVersion;
    //This is a 20 x 8 table
    BYTE m_pSignTable[_kMesonCorrelatorType * 8];
    //This is a 20 x 1 table
    BYTE m_pDeltaTable[_kMesonCorrelatorType];
    void InitialSignTable();

    CFieldFermionKSSU3* m_pW1[24];
    CFieldFermionKSSU3* m_pW2[24];
    deviceSU3Vector** m_pDeviceW1;
    deviceSU3Vector** m_pDeviceW2;
    void CalculateSources(INT gn, INT bn, const CFieldGauge* const* gs, const CFieldBoson* const* bs);

    //uiSite
    void CalculatePropogators();
    void SimplerVersion();

    BYTE* m_pDeviceSignTable;
    BYTE* m_pDeviceDeltaTable;
    CLGComplex* m_pDevicePropogators;
#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex* m_pDevicePropogatorsEveryTimeSlice;
    cuDoubleComplex* m_pResPropogators;
#else
    CLGComplex* m_pDevicePropogatorsEveryTimeSlice;
    CLGComplex* m_pResPropogators;
#endif

    //This is 20 x (Lt - 2)
    

    void InitialBuffers();

public:

#if !_CLG_DOUBLEFLOAT
    TArray<TArray<DOUBLE>> m_lstAverageResults;

    //m_lstResults[conf][type][t]
    TArray<TArray<TArray<cuDoubleComplex>>> m_lstResults;
#else
    TArray<TArray<Real>> m_lstAverageResults;

    //m_lstResults[conf][type][t]
    TArray<TArray<TArray<CLGComplex>>> m_lstResults;
#endif

    UBOOL m_bGaugeFixing;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREMESONCORRELATORSTAGGERED_H_

//=============================================================================
// END OF FILE
//=============================================================================