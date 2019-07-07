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

        , m_uiConfigurationCount(0)
        , m_bShowResult(FALSE)
    {
    }

    ~CMeasurePolyakov();

    virtual void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId);
    virtual void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple);
    virtual void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) {}
    virtual void Average(UINT uiConfigurationCount);
    virtual void Report();
    virtual void Reset();

    virtual UBOOL IsGaugeMeasurement() const { return TRUE; }
    virtual UBOOL IsSourceScanning() const { return FALSE; }

protected:

    deviceSU3* m_pTmpLoop;
    CLGComplex* m_pTraceRes;
    CLGComplex* m_pTmpDeviceSum;
    
    CLGComplex* m_pCorrelator;
    UINT* m_pCorrelatorCounter;
    CLGComplex* m_pHostCorrelator;
    UINT* m_pHostCorrelatorCounter;

    UINT m_uiMaxLengthSq;

    UINT m_uiConfigurationCount;
    UBOOL m_bShowResult;

public:

    TArray<UINT> m_lstR;
    TArray<CLGComplex> m_lstC;
    TArray<CLGComplex> m_lstAverageLoop;

    CLGComplex m_cAverageLoop;
    TArray<CLGComplex> m_lstAverageC;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPOLYAKOV_H_

//=============================================================================
// END OF FILE
//=============================================================================