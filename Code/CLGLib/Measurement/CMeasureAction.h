//=============================================================================
// FILENAME : CMeasureAction.h
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [10/06/2020 nbale]
//=============================================================================

#ifndef _CMEASUREACTION_H_
#define _CMEASUREACTION_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureAction)

class CLGAPI CMeasureAction : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureAction)
public:
    CMeasureAction() : CMeasure(), m_iFermionFieldCount(1) {}

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;
    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Average(UINT uiConfigurationCount) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

    TArray<Real> m_lstData;

    UINT m_iFermionFieldCount;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREACTION_H_

//=============================================================================
// END OF FILE
//=============================================================================