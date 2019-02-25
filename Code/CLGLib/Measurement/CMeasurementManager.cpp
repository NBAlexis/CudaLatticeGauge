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
    ++m_iAcceptedConfigurationCount;
    for (INT i = 0; i < m_lstAllMeasures.Num(); ++i)
    {
        if (NULL != m_lstAllMeasures[i])
        {
            m_lstAllMeasures[i]->OnConfigurationAccepted(pAcceptGauge, pCorrespondingStaple);
        }
    }
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
}

CMeasure* CMeasurementManager::GetMeasureById(BYTE byId) const
{
    return m_mapMeasures.GetAt(byId);
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================