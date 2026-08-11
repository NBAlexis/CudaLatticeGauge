//=============================================================================
// FILENAME : CMeasureMesonCorrelatorStaggeredSimple2.h
// 
// DESCRIPTION:
// This is the 10.1103/PhysRevD.38.2245
// It uses 3 propagators, if wall source is used, gauge fixing is needed
// It supports charged mesons
//
// REVISION: [dd-mm-yy]
//  [08/12/2022 nbale]
//=============================================================================

#ifndef _CMEASUREMESONCORRELATORSTAGGEREDSIMPLE2_H_
#define _CMEASUREMESONCORRELATORSTAGGEREDSIMPLE2_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureMesonCorrelatorStaggeredSimple2)

class CLGAPI CMeasureMesonCorrelatorStaggeredSimple2 : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureMesonCorrelatorStaggeredSimple2)
public:

    enum { _kMesonCorrelatorTypeSimple2 = 16 };

    CMeasureMesonCorrelatorStaggeredSimple2() : CMeasure()
        , m_pDevicePropogators(NULL)
        , m_pResPropogators(NULL)
        , m_pResPropogatorsX(NULL)
        , m_pResPropogatorsY(NULL)
        , m_pResPropogatorsZ(NULL)
        , m_eSource(EFS_Point)
        , m_byFieldID2(0)
    {
        
    }
    ~CMeasureMesonCorrelatorStaggeredSimple2();
    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;

    void OnConfigurationAccepted(INT gn, INT bn, INT tensor2Num, const CFieldGauge* const* gs, const CFieldBoson* const* bs, const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* stp) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }
    UBOOL HasOtherField() const
    {
        return m_byFieldID2 > 0 && m_byFieldID2 != GetFermionFieldId();
    }

protected:

    void BuildSource();
    void IniverseSource(INT gn, INT bn, const CFieldGauge* const* gs, const CFieldBoson* const* bs);
    void ReleaseSource();

    TArray<CFieldFermionKSSU3*> m_pSources;

    Real* m_pDevicePropogators;

    //This is _kMesonCorrelatorTypeSimple2 x (Lt - 1) x (flavour^2)
    //In the case of flavour=2, it is uu, ud, du, dd
    DOUBLE* m_pResPropogators;

    //This is _kMesonCorrelatorTypeSimple2 x (Lz - 1) x (flavour^2)
    DOUBLE* m_pResPropogatorsX;
    DOUBLE* m_pResPropogatorsY;
    DOUBLE* m_pResPropogatorsZ;

public:

    TArray<TArray<DOUBLE>> m_lstAverageResults;
    TArray<TArray<DOUBLE>> m_lstAverageResultsX;
    TArray<TArray<DOUBLE>> m_lstAverageResultsY;
    TArray<TArray<DOUBLE>> m_lstAverageResultsZ;

    //m_lstResults[conf][type][t]
    //for nf = 1 + 1, type is: type1uu, type1ud, type1dd, type1du, type2uu, type2ud, type2dd, type2du, ...
    TArray<TArray<TArray<DOUBLE>>> m_lstResults;
    TArray<TArray<TArray<DOUBLE>>> m_lstResultsX;
    TArray<TArray<TArray<DOUBLE>>> m_lstResultsY;
    TArray<TArray<TArray<DOUBLE>>> m_lstResultsZ;

    EFermionBosonSource m_eSource;
    BYTE m_byFieldID2;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREMESONCORRELATORSTAGGEREDSIMPLE2_H_

//=============================================================================
// END OF FILE
//=============================================================================