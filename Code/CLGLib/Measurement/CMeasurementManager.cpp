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

void CMeasurementManager::OnConfigurationAccepted(INT gaugeNum, INT bosonNum, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldGauge* const* pCorrespondingStaple)
{
    if (!m_bEverResetted && 0 == m_iAcceptedConfigurationCount)
    {
        m_bNeedGaugeSmearing = NeedSmearing();
    }

    ++m_iAcceptedConfigurationCount;

    TArray<const CFieldGauge*> allGauges;
    TArray<const CFieldGauge*> allStaples;
    if (m_bNeedGaugeSmearing && NULL != appGetGaugeSmearing())
    {
        CFieldGauge* pSmearing = NULL;
        CFieldGauge* pSmearingStaple = NULL;
        for (INT i = 0; i < gaugeNum; ++i)
        {
            if (NULL != pAcceptGauge[i])
            {
                //for smearing, we have to use staple
                pSmearing = dynamic_cast<CFieldGauge*>(pAcceptGauge[i]->GetCopy());
                if (NULL != pCorrespondingStaple)
                {
                    pSmearingStaple = dynamic_cast<CFieldGauge*>(pCorrespondingStaple[i]->GetCopy());
                }
                else
                {
                    pSmearingStaple = dynamic_cast<CFieldGauge*>(pAcceptGauge[i]->GetCopy());
                    pAcceptGauge[i]->CalculateOnlyStaple(pSmearingStaple);
                }
                appGetGaugeSmearing()->GaugeSmearing(pSmearing, pSmearingStaple);
                allGauges.AddItem(pSmearing);
                allStaples.AddItem(pSmearingStaple);
            }
            else
            {
                allGauges.AddItem(NULL);
                allStaples.AddItem(NULL);
            }
        }
    }

    //gauge measurement
    for (INT i = 0; i < m_lstAllMeasures.Num(); ++i)
    {
        if (NULL != m_lstAllMeasures[i] && m_lstAllMeasures[i]->IsGaugeMeasurement())
        {
            m_lstAllMeasures[i]->OnConfigurationAccepted(
                gaugeNum,
                bosonNum,
                m_lstAllMeasures[i]->NeedGaugeSmearing() ? allGauges.GetData() : pAcceptGauge,
                pAcceptBoson,
                m_lstAllMeasures[i]->NeedGaugeSmearing() ? allStaples.GetData() : pCorrespondingStaple);
        }
    }

    //z4 source
    UINT uiFieldCount = 0;
    const THashMap<BYTE, TArray<CMeasureStochastic*>> allZ4Fields = HasZ4(uiFieldCount);
    if (uiFieldCount > 0)
    {
        TArray<BYTE> allFieldIdsz4 = allZ4Fields.GetAllKeys();
        for (INT i = 0; i < allFieldIdsz4.Num(); ++i)
        {
            BYTE byFieldIdz4 = allFieldIdsz4[i];
            TArray<CMeasureStochastic*> measures = allZ4Fields.GetAt(byFieldIdz4);
            TArray<CFieldFermion*> leftField;
            TArray<CFieldFermion*> rightField;

            for (UINT j = 0; j < uiFieldCount; ++j)
            {
                CFieldFermion* pF1 = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(byFieldIdz4));
                CFieldFermion* pF2 = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(byFieldIdz4));
                if (CCommonData::m_bStochasticGaussian)
                {
                    pF1->InitialField(EFIT_RandomGaussian);
                }
                else
                {
                    pF1->InitialField(EFIT_RandomZ4);
                }
                pF1->FixBoundary();
                pF1->CopyTo(pF2);
                pF1->InverseD(gaugeNum, bosonNum, 
                    m_bNeedGaugeSmearing ? allGauges.GetData() : pAcceptGauge,
                    pAcceptBoson);

                for (INT k = 0; k < measures.GetCount(); ++k)
                {
                    measures[k]->OnConfigurationAcceptedZ4(
                        gaugeNum,
                        bosonNum,
                        measures[k]->NeedGaugeSmearing() ? allGauges.GetData() : pAcceptGauge,
                        pAcceptBoson,
                        measures[k]->NeedGaugeSmearing() ? allStaples.GetData() : pCorrespondingStaple,
                        pF2, pF1,
                        0 == j, uiFieldCount == j + 1);
                }

                pF1->Return();
                pF2->Return();
            }
        }
    }

    //source scanning measurement
    UBOOL bHasSourceScanning = FALSE;
    const THashMap<BYTE, TArray<CMeasure*>> allScanningFields = HasSourceScanning(bHasSourceScanning);
    if (bHasSourceScanning)
    {
        TArray<BYTE> allFieldIds = allScanningFields.GetAllKeys();
        for (INT i = 0; i < allFieldIds.Num(); ++i)
        {
            BYTE byFieldId = allFieldIds[i];
            TArray<CMeasure*> measures = allScanningFields.GetAt(byFieldId);
            for (UINT x = 1; x < _HC_Lx; ++x)
            {
                SSmallInt4 sourceSite;
                sourceSite.x = static_cast<SBYTE>(x);
                sourceSite.y = static_cast<SBYTE>(_HC_Centery);
                sourceSite.z = static_cast<SBYTE>(_HC_Centerz);
                sourceSite.w = static_cast<SBYTE>(_HC_Centert);

                CFieldFermion* pFermion = dynamic_cast<CFieldFermion*>(appGetLattice()->GetFieldById(byFieldId));
                TArray<CFieldFermion*> sources = pFermion->GetSourcesAtSiteFromPool(
                    gaugeNum, bosonNum,
                    m_lstAllMeasures[i]->NeedGaugeSmearing() ? allGauges.GetData() : pAcceptGauge,
                    pAcceptBoson,
                    sourceSite);

                for (INT j = 0; j < measures.Num(); ++j)
                {
                    measures[j]->SourceSanning(
                        gaugeNum,
                        bosonNum,
                        m_lstAllMeasures[i]->NeedGaugeSmearing() ? allGauges.GetData() : pAcceptGauge,
                        pAcceptBoson,
                        m_lstAllMeasures[i]->NeedGaugeSmearing() ? allStaples.GetData() : pCorrespondingStaple,
                        sources,
                        sourceSite);
                }

                for (INT j = 0; j < sources.Num(); ++j)
                {
                    sources[j]->Return();
                }
            }
        }
    }

    INT iCreated = allGauges.Num();
    for (INT i = 0; i < iCreated; ++i)
    {
        appSafeDelete(allGauges[i]);
        appSafeDelete(allStaples[i]);
    }
}

void CMeasurementManager::OnUpdateFinished(UBOOL bReport)
{
    
    for (INT i = 0; i < m_lstAllMeasures.Num(); ++i)
    {
        if (NULL != m_lstAllMeasures[i])
        {
            assert(m_iAcceptedConfigurationCount == m_lstAllMeasures[i]->GetConfigurationCount());
            m_lstAllMeasures[i]->Average();
            if (bReport)
            {
                m_lstAllMeasures[i]->Report();
            }
        }
    }
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

THashMap<BYTE, TArray<CMeasure*>> CMeasurementManager::HasSourceScanning(UBOOL& bHasSourceScanning) const
{
    bHasSourceScanning = FALSE;
    THashMap<BYTE, TArray<CMeasure*>> ret;
    for (INT i = 0; i < m_lstAllMeasures.Num(); ++i)
    {
        if (NULL != m_lstAllMeasures[i] && m_lstAllMeasures[i]->IsSourceScanning())
        {
            bHasSourceScanning = TRUE;
            BYTE byFieldId = m_lstAllMeasures[i]->GetFermionFieldId();
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

THashMap<BYTE, TArray<CMeasureStochastic*>> CMeasurementManager::HasZ4(UINT &uiFieldCount) const
{
    THashMap<BYTE, TArray<CMeasureStochastic*>> ret;
    uiFieldCount = 0;
    for (INT i = 0; i < m_lstAllMeasures.Num(); ++i)
    {
        CMeasureStochastic* pStoch = dynamic_cast<CMeasureStochastic*>(m_lstAllMeasures[i]);
        if (NULL != pStoch && pStoch->IsZ4Source())
        {
            BYTE byFieldId = m_lstAllMeasures[i]->GetFermionFieldId();
            if (ret.Exist(byFieldId))
            {
                TArray<CMeasureStochastic*> lst = ret.GetAt(byFieldId);
                lst.AddItem(pStoch);
                ret.SetAt(byFieldId, lst);
            }
            else
            {
                TArray<CMeasureStochastic*> lst;
                lst.AddItem(pStoch);
                ret.SetAt(byFieldId, lst);
            }
            if (pStoch->GetFieldCount() > uiFieldCount)
            {
                uiFieldCount = pStoch->GetFieldCount();
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