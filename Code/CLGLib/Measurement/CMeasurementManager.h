//=============================================================================
// FILENAME : CMeasurementManager.h
// 
// DESCRIPTION:
// This is the class collecting all measurements
// 
// NOTE: m_bNeedGaugeSmearing will be infected
//
// REVISION:
//  [mm/dd/yy]
//  [01/29/2019 nbale]
//=============================================================================

#ifndef _CMEASUREMENTMANAGER_H_
#define _CMEASUREMENTMANAGER_H_

#include "CMeasureData.h"

__BEGIN_NAMESPACE

class CLGAPI CMeasurementManager
{
public:
    CMeasurementManager(class CLatticeData* pOwner);
    ~CMeasurementManager();

    void OnConfigurationAccepted(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, const class CFieldGauge* const* pCorrespondingStaple);
    void OnUpdateFinished(UBOOL bReport = TRUE);
    void Reset();
    void Report();
    void AverageAll();
    TArray<Real> AverageReals() const;
    TArray<Real> LastReals() const;

    CMeasure* GetMeasureById(BYTE byId) const;

    CMeasureData* GetMeasureData(const CCString& sKey) const;

    template<typename T>
    void AddOneConfigurationResult(CMeasure* pMeasure, const CCString& sName, T data)
    {
        appAssert(NULL != pMeasure);
        const CCString sKey = _T("measure") + appToString(pMeasure->GetId()) + _T(".") + sName;
        CMeasureDataT<T>* pData = GetOrCreateDataArray<T>(sKey);
        pData->Get().AddItem(data);
        appGeneral(_T("[CMeasurementManager] Add result: key='%s' count=%d\n"),
            sKey.c_str(), pData->Get().Num());
    }

    TArray<CMeasure*> m_lstAllMeasures;
    THashMap<BYTE, CMeasure*> m_mapMeasures;

protected:

    UINT m_iAcceptedConfigurationCount;
    class CLatticeData* m_pOwner;
    UBOOL m_bNeedGaugeSmearing;
    UBOOL m_bEverResetted;

    CMemStack m_MeasureDataMem;
    THashMap<CCString, CMeasureData*> m_mapMeasureData;

    template<typename T>
    CMeasureDataT<T>* GetOrCreateDataArray(const CCString& sKey)
    {
        CMeasureData* pBase = m_mapMeasureData.GetAt(sKey);
        if (NULL != pBase)
        {
            CMeasureDataT<T>* pTyped = dynamic_cast<CMeasureDataT<T>*>(pBase);
            appAssert(NULL != pTyped);
            return pTyped;
        }

        void* pMem = m_MeasureDataMem.PushBytes(static_cast<INT>(sizeof(CMeasureDataT<T>)));
        CMeasureDataT<T>* pNewData = new (pMem) CMeasureDataT<T>();
        m_mapMeasureData.SetAt(sKey, pNewData);
        return pNewData;
    }

    THashMap<BYTE, TArray<CMeasure*>> HasSourceScanning(UBOOL& bHasSourceScanning) const;
    THashMap<BYTE, TArray<CMeasureStochastic*>> HasZ4(UINT& uiFieldCount) const;
    UBOOL NeedSmearing() const;
};

__END_NAMESPACE

#endif //#ifndef _CMEASUREMENTMANAGER_H_

//=============================================================================
// END OF FILE
//=============================================================================