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

    virtual void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple);
    virtual void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) {}
    virtual void Average(UINT uiConfigurationCount);
    virtual void Report();
    virtual void Reset();

    virtual UBOOL IsGaugeMeasurement() const { return TRUE; }
    virtual UBOOL IsSourceScanning() const { return FALSE; }

protected:

    TArray<Real> m_lstData;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPLAQUTTEENERGY_H_

//=============================================================================
// END OF FILE
//=============================================================================