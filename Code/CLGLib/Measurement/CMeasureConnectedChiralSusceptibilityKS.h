//=============================================================================
// FILENAME : CMeasureConnectedSusceptibilityKS.h
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [02/22/2019 nbale]
//=============================================================================

#ifndef _CMEASURECONNECTEDCHIRALSUSCEPTIBILITYKS_H_
#define _CMEASURECONNECTEDCHIRALSUSCEPTIBILITYKS_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureConnectedSusceptibilityKS)

class CLGAPI CMeasureConnectedSusceptibilityKS : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureConnectedSusceptibilityKS)
public:

    CMeasureConnectedSusceptibilityKS() : CMeasure()
        , m_pSourceZero(NULL)
        , m_bShowResult(FALSE)
    {
        
    }
    ~CMeasureConnectedSusceptibilityKS();
    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters&, BYTE byId) override;

    void OnConfigurationAccepted(const CFieldGauge* pGaugeField, const CFieldGauge* pStapleField) override;
    void SourceSanning(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& site) override {}
    void Report() override;

    UBOOL IsGaugeMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

protected:

    CFieldFermion* m_pSourceZero;

public:

    UBOOL m_bShowResult;
};

__END_NAMESPACE

#endif //#ifndef _CMEASURECONNECTEDCHIRALSUSCEPTIBILITYKS_H_

//=============================================================================
// END OF FILE
//=============================================================================