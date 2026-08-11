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
    CMeasurePlaqutteEnergy() 
        : CMeasure() 
        , m_bV0(FALSE)
    {
    }

    void OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple) override;
    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId) override;
    void Report() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

    UBOOL m_bV0;
    TArray<DOUBLE> m_lstV0;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREPLAQUTTEENERGY_H_

//=============================================================================
// END OF FILE
//=============================================================================