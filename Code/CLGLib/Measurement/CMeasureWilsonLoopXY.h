//=============================================================================
// FILENAME : CMeasureWilsonLoopXY.h
// 
// DESCRIPTION:
// 
// 
//
// REVISION:
//  [09/24/2021 nbale]
//=============================================================================

#ifndef _CMEASUREWILSONLOOPXY_H_
#define _CMEASUREWILSONLOOPXY_H_

__BEGIN_NAMESPACE


__CLG_REGISTER_HELPER_HEADER(CMeasureWilsonLoopXY)

class CLGAPI CMeasureWilsonLoopXY : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureWilsonLoopXY)

public:

    CMeasureWilsonLoopXY()
        : CMeasure()
        , m_pTmpLoop(NULL)

        , m_pCorrelator(NULL)
        , m_pCorrelatorCounter(NULL)
        , m_pHostCorrelator(NULL)
        , m_pHostCorrelatorCounter(NULL)
        , m_uiMaxLengthSq(1)
    {
    }

    ~CMeasureWilsonLoopXY();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    deviceSU3* m_pTmpLoop;

    CLGComplex* m_pCorrelator;
    UINT* m_pCorrelatorCounter;
    CLGComplex* m_pHostCorrelator;
    UINT* m_pHostCorrelatorCounter;

    UINT m_uiMaxLengthSq;

public:

    TArray<UINT> m_lstR;
    //m_lstC[conf][r][t]
    TArray<TArray<TArray<CLGComplex>>> m_lstC;
    TArray<TArray<CLGComplex>> m_lstAverageC;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREWILSONLOOPXY_H_

//=============================================================================
// END OF FILE
//=============================================================================