//=============================================================================
// FILENAME : CMeasurementManager.h
// 
// DESCRIPTION:
// This is the class collecting all measurements
//
// REVISION:
//  [01/29/2019 nbale]
//=============================================================================

#ifndef _CMEASUREMENTMANAGER_H_
#define _CMEASUREMENTMANAGER_H_

__BEGIN_NAMESPACE

class CLGAPI CMeasurementManager
{
public:
    CMeasurementManager(class CLatticeData* pOwner) 
        : m_iAcceptedConfigurationCount(0)
        , m_pOwner(pOwner) 
        , m_bNeedGaugeSmearing(FALSE)
        , m_bEverResetted(FALSE)
    {
    }

    void OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple);
    void OnUpdateFinished(UBOOL bReport = TRUE);
    void Reset();
    void Report();

    CMeasure* GetMeasureById(BYTE byId) const;

    TArray<CMeasure*> m_lstAllMeasures;
    THashMap<BYTE, CMeasure*> m_mapMeasures;

protected:

    UINT m_iAcceptedConfigurationCount;
    class CLatticeData* m_pOwner;
    UBOOL m_bNeedGaugeSmearing;
    UBOOL m_bEverResetted;

    THashMap<BYTE, TArray<CMeasure*>> HasSourceScanning() const;
    UBOOL NeedSmearing() const;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREMENTMANAGER_H_

//=============================================================================
// END OF FILE
//=============================================================================