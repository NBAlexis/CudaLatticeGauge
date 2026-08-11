//=============================================================================
// FILENAME : CMeasureMesonCorrelatorStaggered.h
// 
// DESCRIPTION:
// This is the class for staggered meson correlator measurement.
// Supports 20 meson channels with up to 3 SO3-degenerate sub-channels each.
// Computes both wall-to-wall (W2W) and point-to-point (P2P) contractions.
//
// REVISION:
//  [02/22/2019 nbale]
//  [06/08/2026 nbale] Added SO3 sub-channels and P2P contraction
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
    enum { _kMaxSubChannels = 3 };

    CMeasureMesonCorrelatorStaggered() : CMeasure()
        , m_pDeviceW1(NULL)
        , m_pDeviceW2(NULL)
        , m_pDeviceKernelBuffer(NULL)
        , m_pDeviceKernelEveryTimeSlice(NULL)
        , m_pW2WCorrelator(NULL)
        , m_pP2PCorrelator(NULL)
        , m_pW2WPArray(NULL)
        , m_pP2PPArray(NULL)
        , m_bGaugeFixing(FALSE)
    {
        memset(m_pSignTable, 0, sizeof(m_pSignTable));
        memset(m_pDeltaTable, 0, sizeof(m_pDeltaTable));
        memset(m_nSubChannels, 0, sizeof(m_nSubChannels));
    }
    ~CMeasureMesonCorrelatorStaggered();
    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;

    void OnConfigurationAccepted(INT gn, INT bn, INT tensor2Num, const CFieldGauge* const* gs, const CFieldBoson* const* bs, const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* stp) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

    static void InitialSignTable(
        BYTE signtable[_kMesonCorrelatorType][_kMaxSubChannels][8],
        BYTE deltatable[_kMesonCorrelatorType][_kMaxSubChannels],
        BYTE nsubchannels[_kMesonCorrelatorType]);

protected:

    // Sign table: [type][sub-channel][A], values 0 or 1
    // Sign contribution = m_pSignTable[type][sub][A], even=+1, odd=-1
    BYTE m_pSignTable[_kMesonCorrelatorType][_kMaxSubChannels][8];
    // Delta table: [type][sub-channel], values 0-7 (spatial displacement index)
    BYTE m_pDeltaTable[_kMesonCorrelatorType][_kMaxSubChannels];

    // 24 wall-source fermion fields: 8 shifts x 3 colors
    CFieldFermionKSSU3* m_pW1[24];
    CFieldFermionKSSU3* m_pW2[24];
    // Device pointer arrays: 24 field pointers each (W2W/P2P shared)
    deviceSU3Vector** m_pDeviceW1;   // DSource (D^-1)
    deviceSU3Vector** m_pDeviceW2;   // DdSource (Dd^-1)
    void CalculateSources(INT gn, INT bn, const CFieldGauge* const* gs, const CFieldBoson* const* bs);

    // W2W contraction: two kernel calls per (A,B) + CPU combination
    void CalculateW2W();
    // P2P contraction: one kernel call per (A,B,delta)
    void CalculateP2P();
    // Shared projection: C_d(t) = sum_{A,B} sign(A)sign(B) p(t,A,B,d)
    void ProjectPArray(const cuDoubleComplex* pArray, cuDoubleComplex* correlator);

    // GPU buffers (per-site, shared by W2W and P2P)
    cuDoubleComplex* m_pDeviceKernelBuffer;         // Volume * 9
    cuDoubleComplex* m_pDeviceKernelEveryTimeSlice; // Volume_xyz

    // Host per-config correlator buffers
    cuDoubleComplex* m_pW2WCorrelator;   // 20 * 3 * Lt
    cuDoubleComplex* m_pP2PCorrelator;   // 20 * 3 * Lt

    void InitialBuffers();

public:

    // Per-config p arrays (overwritten each config)
    cuDoubleComplex* m_pW2WPArray;       // Lt * 8 * 8 * 8
    cuDoubleComplex* m_pP2PPArray;       // Lt * 8 * 8 * 8

    TArray<TArray<DOUBLE>> m_lstAverageResults;

    // W2W per-sub-channel correlator: [conf][type][sub][t]
    TArray<TArray<TArray<TArray<cuDoubleComplex>>>> m_lstW2WCorrelator;
    // W2W SO3 combined correlator: [conf][type][t]
    TArray<TArray<TArray<cuDoubleComplex>>> m_lstW2WCombinedCorrelator;
    // P2P per-sub-channel correlator: [conf][type][sub][t]
    TArray<TArray<TArray<TArray<cuDoubleComplex>>>> m_lstP2PCorrelator;
    // P2P SO3 combined correlator: [conf][type][t]
    TArray<TArray<TArray<cuDoubleComplex>>> m_lstP2PCombinedCorrelator;

    UBOOL m_bGaugeFixing;
    // Number of SO3-degenerate sub-channels per type (1 or 3)
    BYTE m_nSubChannels[_kMesonCorrelatorType];
};

/**
 * Even: 1, Odd: -1
 */
static inline INT __eta(INT x, INT y, INT z, INT i)
{
    switch (i)
    {
    case 1:
        return 0;
    case 2:
        return x;
    case 3:
        return x + y;
    case 4:
        return x + y + z;
    default:
        break;
    }
    return x + y + z;
}

static inline INT __xi(INT x, INT y, INT z, INT i)
{
    switch (i)
    {
    case 1:
        return y + z;
    case 2:
        return z;
    case 3:
        return 1;
    case 4:
        return 1;
    default:
        break;
    }
    return x + y + z;
}

__END_NAMESPACE

#endif //#ifndef _CMEASUREMESONCORRELATORSTAGGERED_H_

//=============================================================================
// END OF FILE
//=============================================================================
