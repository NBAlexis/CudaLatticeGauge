//=============================================================================
// FILENAME : CMeasurePolyakovSUN.h
// 
// DESCRIPTION:
// This is measurement for Polyakov loop and static quark potential
// Suitable for translational invarient case
//
// REVISION:
//  [07/02/2024 nbale]
//=============================================================================

#ifndef _CMEASUREPOLYAKOVSUN_H_
#define _CMEASUREPOLYAKOVSUN_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasurePolyakov2)

class CLGAPI CMeasurePolyakov2 : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasurePolyakov2)
public:

    CMeasurePolyakov2()
        : CMeasure()
        , m_pTraceRes(NULL)

        , m_pCorrelator(NULL)
        , m_pCorrelatorCounter(NULL)
        , m_pHostCorrelator(NULL)
        , m_pHostCorrelatorCounter(NULL)
        , m_uiMaxLengthSq(1)
    {

    }

    ~CMeasurePolyakov2();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    cuDoubleComplex* m_pTraceRes;
    
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

#endif //#ifndef _CMEASUREPOLYAKOVSUN_H_

//=============================================================================
// END OF FILE
//=============================================================================