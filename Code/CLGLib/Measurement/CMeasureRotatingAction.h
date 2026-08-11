//=============================================================================
// FILENAME : CMeasureRotatingAction.h
//
// DESCRIPTION:
// This is the class for rotating gauge action energy measurement
//
// REVISION:
//  [05/22/2026 nbale]
//=============================================================================

#ifndef _CMEASUREROTATINGACTION_H_
#define _CMEASUREROTATINGACTION_H_

__BEGIN_NAMESPACE

__CLG_REGISTER_HELPER_HEADER(CMeasureRotatingAction)

class CLGAPI CMeasureRotatingAction : public CMeasure
{
    __CLGDECLARE_CLASS(CMeasureRotatingAction)
public:
    CMeasureRotatingAction()
        : CMeasure()
        , m_iActionIndex(1)
    {
    }

    void Initial(class CMeasurementManager* pOwner, class CLatticeData* pLatticeData, const CParameters& param, BYTE byId) override;
    void OnConfigurationAccepted(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, const class CFieldGauge* const* pCorrespondingStaple) override;
    void Report() override;
    void Reset() override;

    UBOOL IsGaugeOrBosonMeasurement() const override { return TRUE; }
    UBOOL IsSourceScanning() const override { return FALSE; }

    INT m_iActionIndex;

    TArray<DOUBLE> m_lstS0;
    TArray<DOUBLE> m_lstS1;
    TArray<DOUBLE> m_lstS2;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREROTATINGACTION_H_

//=============================================================================
// END OF FILE
//=============================================================================
