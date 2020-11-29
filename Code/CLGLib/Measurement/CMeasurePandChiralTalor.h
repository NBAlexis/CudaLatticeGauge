//=============================================================================
// FILENAME : CMeasurePandChiralTalor.h
// 
// DESCRIPTION:
// This is measurement for Talor expansion of rotation frame
// The angular momentum of each site of two lines is calculated.
// It is impractical to calculate the angular momentum for each site (a 16^4 lattice has 65536 sites)
//
// REVISION:
//  [11/25/2020 nbale]
//=============================================================================

#ifndef _CMEASUREPANDCHIRALTALOR_H_
#define _CMEASUREPANDCHIRALTALOR_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasurePandChiralTalor)

__DEFINE_ENUM(ECPCTTraceType,
    ECPCTTT_D,
    ECPCTTT_MD,
    ECPCTTT_DMD,
    ECPCTTT_MDMD,
    ECPCTTT_DMDMD,

    ECPCTTT_Max,
    );

class CLGAPI CMeasurePandChiralTalor : public CMeasureStochastic
{
    __CLGDECLARE_CLASS(CMeasurePandChiralTalor)
public:
    CMeasurePandChiralTalor()
        : CMeasureStochastic()
        , m_uiConfigurationCount(0)
        , m_bShowResult(TRUE)
    {
    }

    ~CMeasurePandChiralTalor();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAcceptedZ4(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const class CFieldFermion* pZ4, const class CFieldFermion* pInverseZ4, UBOOL bStart, UBOOL bEnd) override;
    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Average(UINT uiConfigurationCount) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return FALSE; }
    UBOOL IsZ4Source() const override { return TRUE; }

protected:
    
    UINT m_uiConfigurationCount;
    UBOOL m_bShowResult;

public:

#if _CLG_DOUBLEFLOAT
    CLGComplex m_cTmpSum[ECPCTTT_Max];
    TArray<CLGComplex> m_lstTraceRes[ECPCTTT_Max];
    TArray<CLGComplex> m_lstPolyakov;
    TArray<Real> m_lstPolyakovSOmega;
    TArray<Real> m_lstPolyakovSOmegaSq;
#else
    cuDoubleComplex m_cTmpSum[ECPCTTT_Max];
    TArray<cuDoubleComplex> m_lstTraceRes[ECPCTTT_Max];
    TArray<cuDoubleComplex> m_lstPolyakov;
    TArray<DOUBLE> m_lstPolyakovSOmega;
    TArray<DOUBLE> m_lstPolyakovSOmegaSq;
#endif
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPANDCHIRALTALOR_H_

//=============================================================================
// END OF FILE
//=============================================================================