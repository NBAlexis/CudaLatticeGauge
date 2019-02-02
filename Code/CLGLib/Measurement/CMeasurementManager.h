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
    CMeasurementManager(class CLatticeData* pOwner) : m_iAcceptedConfigurationCount(0), m_pOwner(pOwner) {}

    void OnConfigurationAccepted();
    void OnUpdateFinished(UBOOL bReport = TRUE);
    void Reset();

    CMeasure* GetMeasureById(BYTE byId) const;

    TArray<CMeasure*> m_lstAllMeasures;
    THashMap<BYTE, CMeasure*> m_mapMeasures;

protected:

    UINT m_iAcceptedConfigurationCount;
    class CLatticeData* m_pOwner;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREMENTMANAGER_H_

//=============================================================================
// END OF FILE
//=============================================================================