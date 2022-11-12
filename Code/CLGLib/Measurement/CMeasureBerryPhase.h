//=============================================================================
// FILENAME : CMeasureBerryPhase.h
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [10/13/2020 nbale]
//=============================================================================

#ifndef _CMEASUREBERRYPHASE_H_
#define _CMEASUREBERRYPHASE_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureBerryPhase)

class CLGAPI CMeasureBerryPhase : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureBerryPhase)
public:
    CMeasureBerryPhase()
        : CMeasure()
        , m_bWilsonDirac(TRUE)
        , m_bGuageFixing(FALSE)
        , m_bShowRes(FALSE)
        , m_uiConfigurationCount(0)
        , m_pMomentumField(NULL)
        , m_pU1Field(NULL)
        , m_pGaugeFixing(NULL)
    {
        
    }

    ~CMeasureBerryPhase();

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple = NULL) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Average(UINT uiConfigurationCount) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

    void CalculateMomentumSpacePhiWilsonDiracForPoint(const SSmallInt4& xprime, const CFieldGaugeSU3* pGauge);
    void CalculateMomentumSpacePhiWilsonDirac(const CFieldGaugeSU3* pGauge);
    void CalculateU1FieldWilsonDirac();
    void CalculateMomentumSpacePhiKSForPoint(const SSmallInt4& xprime, const CFieldGaugeSU3* pGauge);
    void CalculateMomentumSpacePhiKS(const CFieldGaugeSU3* pGauge);
    void CalculateU1FieldKS();

    void CalculateBerryPhase(BYTE byGaugeFieldId);
    void AllocateBuffers();

#if !_CLG_DOUBLEFLOAT
    TArray<TArray<DOUBLE>> m_lstData;
    TArray<TArray<DOUBLE>> m_lstDataXY;
    TArray<TArray<DOUBLE>> m_lstDataXZ;
    TArray<TArray<DOUBLE>> m_lstDataXT;
    TArray<TArray<DOUBLE>> m_lstDataYZ;
    TArray<TArray<DOUBLE>> m_lstDataYT;
    TArray<TArray<DOUBLE>> m_lstDataZT;
#else
    TArray<TArray<Real>> m_lstData;
    TArray<TArray<Real>> m_lstDataXY;
    TArray<TArray<Real>> m_lstDataXZ;
    TArray<TArray<Real>> m_lstDataXT;
    TArray<TArray<Real>> m_lstDataYZ;
    TArray<TArray<Real>> m_lstDataYT;
    TArray<TArray<Real>> m_lstDataZT;
#endif

    UBOOL m_bWilsonDirac;
    UBOOL m_bGuageFixing;
    UBOOL m_bShowRes;
    UINT m_uiConfigurationCount;

    CFieldFermion* m_pMomentumField;
    CFieldGaugeU1* m_pU1Field;
    CFieldGaugeSU3* m_pGaugeFixing;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREBERRYPHASE_H_

//=============================================================================
// END OF FILE
//=============================================================================