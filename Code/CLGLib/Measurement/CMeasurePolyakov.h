//=============================================================================
// FILENAME : CMeasurePolyakov.h
// 
// DESCRIPTION:
// This is measurement for Polyakov loop and static quark potential
// Suitable for translational invarient case
//
// REVISION:
//  [07/07/2019 nbale]
//=============================================================================

#ifndef _CMEASUREPOLYAKOV_H_
#define _CMEASUREPOLYAKOV_H_

__BEGIN_NAMESPACE


__CLG_REGISTER_HELPER_HEADER(CMeasurePolyakov)

class CLGAPI CMeasurePolyakov : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasurePolyakov)

public:

    enum { _kGammaInInterests = 8, };

    CMeasurePolyakov()
        : CMeasure()
          , m_pTmpLoop(NULL)
          , m_pTraceRes(NULL)
          , m_pTmpDeviceSum(NULL)

          , m_pCorrelator(NULL)
          , m_pCorrelatorCounter(NULL)
          , m_pHostCorrelator(NULL)
          , m_pHostCorrelatorCounter(NULL)
          , m_uiMaxLengthSq(1)
    {
    }

    ~CMeasurePolyakov();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    deviceSU3* m_pTmpLoop;
    CLGComplex* m_pTraceRes;
    CLGComplex* m_pTmpDeviceSum;
    
    CLGComplex* m_pCorrelator;
    UINT* m_pCorrelatorCounter;
    CLGComplex* m_pHostCorrelator;
    UINT* m_pHostCorrelatorCounter;

    UINT m_uiMaxLengthSq;

public:

    TArray<UINT> m_lstR;
    TArray<CLGComplex> m_lstC;
    TArray<CLGComplex> m_lstAverageC;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPOLYAKOV_H_

//=============================================================================
// END OF FILE
//=============================================================================