//=============================================================================
// FILENAME : CMeasurePlaqutteEnergy.h
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [01/29/2019 nbale]
//=============================================================================

#ifndef _CMEASUREPLAQUTTEENERGY_H_
#define _CMEASUREPLAQUTTEENERGY_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasurePlaqutteEnergy)

class CLGAPI CMeasurePlaqutteEnergy : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasurePlaqutteEnergy)
public:
    CMeasurePlaqutteEnergy() : CMeasure() {}

    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Report() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPLAQUTTEENERGY_H_

//=============================================================================
// END OF FILE
//=============================================================================