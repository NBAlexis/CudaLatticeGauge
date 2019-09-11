//=============================================================================
// FILENAME : CMeasurementManager.cpp
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [01/29/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

void CMeasurementManager::OnConfigurationAccepted(const CFieldGauge* pAcceptGauge, const CFieldGauge* pCorrespondingStaple)
{
    if (!m_bEverResetted && 0 == m_iAcceptedConfigurationCount)
    {
        m_bNeedGaugeSmearing = NeedSmearing();
    }

    ++m_iAcceptedConfigurationCount;

    CFieldGauge* pSmearing = NULL;
    CFieldGauge* pSmearingStaple = NULL;
    if (m_bNeedGaugeSmearing && NULL != appGetGaugeSmearing())
    {
        //for smearing, we have to use staple
        pSmearing = dynamic_cast<CFieldGauge*>(pAcceptGauge->GetCopy());
        if (NULL != pCorrespondingStaple)
        {
            pSmearingStaple = dynamic_cast<CFieldGauge*>(pCorrespondingStaple->GetCopy());
        }
        else
        {
            pSmearingStaple = dynamic_cast<CFieldGauge*>(pAcceptGauge->GetCopy());
            pAcceptGauge->CalculateOnlyStaple(pSmearingStaple);
        }
        appGetGaugeSmearing()->GaugeSmearing(pSmearing, pSmearingStaple);
    }

    //gauge measurement
    for (INT i = 0; i < m_lstAllMeasures.Num(); ++i)
    {
        if (NULL != m_lstAllMeasures[i] && m_lstAllMeasures[i]->IsGaugeMeasurement())
        {
            m_lstAllMeasures[i]->OnConfigurationAccepted(
                m_lstAllMeasures[i]->NeedGaugeSmearing() ? pSmearing : pAcceptGauge,
                m_lstAllMeasures[i]->NeedGaugeSmearing() ? pSmearingStaple : pCorrespondingStaple);
        }
    }

    //source scanning measurement
    const THashMap<BYTE, TArray<CMeasure*>> allScanningFields = HasSourceScanning();
    TArray<BYTE> allFieldIds = allScanningFields.GetAllKeys();
    for (INT i = 0; i < allFieldIds.Num(); ++i)
    {
        BYTE byFieldId = allFieldIds[i];
        TArray<CMeasure*> measures = allScanningFields.GetAt(byFieldId);
        for (UINT x = 1; x < _HC_Lx; ++x)
        {
            SSmallInt4 sourceSite;
            sourceSite.x = static_cast<SBYTE>(x);
            sourceSite.y = CCommonData::m_sCenter.y;
            sourceSite.z = CCommonData::m_sCenter.z;
            sourceSite.w = CCommonData::m_sCenter.w;

            CFieldFermion* pFermion = dynamic_cast<CFieldFermion*>(appGetLattice()->GetFieldById(byFieldId));
            TArray<CFieldFermion*> sources = pFermion->GetSourcesAtSiteFromPool(
                m_lstAllMeasures[i]->NeedGaugeSmearing() ? pSmearing : pAcceptGauge,
                sourceSite);

            for (INT j = 0; j < measures.Num(); ++j)
            {
                measures[j]->SourceSanning(
                    m_lstAllMeasures[i]->NeedGaugeSmearing() ? pSmearing : pAcceptGauge,
                    m_lstAllMeasures[i]->NeedGaugeSmearing() ? pSmearingStaple : pCorrespondingStaple,
                    sources, 
                    sourceSite);
            }

            for (INT j = 0; j < sources.Num(); ++j)
            {
                sources[j]->Return();
            }
        }
    }

    appSafeDelete(pSmearing);
    appSafeDelete(pSmearingStaple);
}

void CMeasurementManager::OnUpdateFinished(UBOOL bReport)
{
    for (INT i = 0; i < m_lstAllMeasures.Num(); ++i)
    {
        if (NULL != m_lstAllMeasures[i])
        {
            m_lstAllMeasures[i]->Average(m_iAcceptedConfigurationCount);
            if (bReport)
            {
                m_lstAllMeasures[i]->Report();
            }
        }
    }
    m_iAcceptedConfigurationCount = 0;
}

void CMeasurementManager::Reset()
{
    m_iAcceptedConfigurationCount = 0;
    for (INT i = 0; i < m_lstAllMeasures.Num(); ++i)
    {
        if (NULL != m_lstAllMeasures[i])
        {
            m_lstAllMeasures[i]->Reset();
        }
    }

    m_bNeedGaugeSmearing = NeedSmearing();
    m_bEverResetted = TRUE;
}

void CMeasurementManager::Report()
{
    for (INT i = 0; i < m_lstAllMeasures.Num(); ++i)
    {
        if (NULL != m_lstAllMeasures[i])
        {
            m_lstAllMeasures[i]->Report();
        }
    }
    
#if !_CLG_DEBUG
    appFlushLog();
#endif
}

CMeasure* CMeasurementManager::GetMeasureById(BYTE byId) const
{
    return m_mapMeasures.GetAt(byId);
}

THashMap<BYTE, TArray<CMeasure*>> CMeasurementManager::HasSourceScanning() const
{
    THashMap<BYTE, TArray<CMeasure*>> ret;
    for (INT i = 0; i < m_lstAllMeasures.Num(); ++i)
    {
        if (NULL != m_lstAllMeasures[i] && m_lstAllMeasures[i]->IsSourceScanning())
        {
            BYTE byFieldId = m_lstAllMeasures[i]->GetFieldId();
            if (ret.Exist(byFieldId))
            {
                TArray<CMeasure*> lst = ret.GetAt(byFieldId);
                lst.AddItem(m_lstAllMeasures[i]);
                ret.SetAt(byFieldId, lst);
            }
            else
            {
                TArray<CMeasure*> lst;
                lst.AddItem(m_lstAllMeasures[i]);
                ret.SetAt(byFieldId, lst);
            }
        }
    }
    return ret;
}

UBOOL CMeasurementManager::NeedSmearing() const
{
    for (INT i = 0; i < m_lstAllMeasures.Num(); ++i)
    {
        if (NULL != m_lstAllMeasures[i] && m_lstAllMeasures[i]->NeedGaugeSmearing())
        {
            return TRUE;
        }
    }
    return FALSE;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================