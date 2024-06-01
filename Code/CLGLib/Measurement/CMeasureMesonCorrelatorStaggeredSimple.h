//=============================================================================
// FILENAME : CMeasureMesonCorrelatorStaggeredSimple.h
// 
// DESCRIPTION:
// This is the 10.1016/0550-3213(88)90396-3 (see the appendix)
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
    {
        
    }
    ~CMeasureMesonCorrelatorStaggeredSimple();
    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;

    void OnConfigurationAccepted(INT gn, INT bn, const CFieldGauge* const* gs, const CFieldBoson* const* bs, const CFieldGauge* const* pStapleField) override;
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
};

#pragma region device functions

static __device__ __inline__ SBYTE _deviceStaggeredFermionSimplePhase(const SSmallInt4& sSite, BYTE byType)
{
    SBYTE ret = 1;
    switch (byType)
    {
    case 1:
        ret = 3 - ((sSite.x & 1) << 1)
                - ((sSite.y & 1) << 1)
                - ((sSite.z & 1) << 1);
        //printf("shift check%d = %d\n", static_cast<INT>(ret), ((sSite.x & 1) ? -1 : 1) + ((sSite.y & 1) ? -1 : 1) + ((sSite.z & 1) ? -1 : 1));
        break;
    case 2:
        ret = 3 - (((sSite.x + sSite.y) & 1) << 1)
                - (((sSite.y + sSite.z) & 1) << 1)
                - (((sSite.x + sSite.z) & 1) << 1);
        break;
    case 3:
        ret = 1 - (((sSite.x + sSite.y + sSite.z) & 1) << 1);
        break;
    case 4:
        ret = 1 - ((sSite.x & 1) << 1);
        break;
    case 5:
        ret = 1 - ((sSite.y & 1) << 1);
        break;
    case 6:
        ret = 1 - ((sSite.z & 1) << 1);
        break;
    case 7:
        ret = 1 - (((sSite.x + sSite.y) & 1) << 1);
        break;
    case 8:
        ret = 1 - (((sSite.y + sSite.z) & 1) << 1);
        break;
    case 9:
        ret = 1 - (((sSite.x + sSite.z) & 1) << 1);
        break;
    default:
        break;
    }
    return ret;
}

#pragma endregion

__END_NAMESPACE

#endif //#ifndef _CMEASUREMESONCORRELATORSTAGGEREDSIMPLE_H_

//=============================================================================
// END OF FILE
//=============================================================================