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
    void OnConfigurationAccepted(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple) override;
    void Report() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

    UINT m_iFermionFieldCount;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREACTION_H_

//=============================================================================
// END OF FILE
//=============================================================================