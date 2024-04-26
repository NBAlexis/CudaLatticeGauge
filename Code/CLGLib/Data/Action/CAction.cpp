//=============================================================================
// FILENAME : CAction.cpp
// 
// DESCRIPTION:
// Some common implementations
//
// REVISION:
//  [04/25/2024 nbale]
//=============================================================================
#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

/**
* Default case is only one gauge field with fieldid = 1
*/
void CAction::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    param.FetchValueArrayBYTE(_T("GaugeFields"), m_byGaugeFieldIds);
    param.FetchValueArrayBYTE(_T("BosonFields"), m_byBosonFieldIds);

    if (0 == m_byGaugeFieldIds.Num() && 0 == m_byBosonFieldIds.Num())
    {
        m_byGaugeFieldIds.AddItem(1);
    }

    m_pOwner = pOwner;
    m_byActionId = byId;

    DOUBLE fBeta = 0.1;
    if (param.FetchValueDOUBLE(_T("Beta"), fBeta))
    {
        CCommonData::m_fBeta = fBeta;
    }
}

DOUBLE CAction::Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldGauge* const* stableFields)
{
    if (1 == m_byGaugeFieldIds.Num() && 0 == m_byBosonFieldIds.Num())
    {
        INT idx = GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
        return EnergySingleField(bBeforeEvolution, gaugeFields[idx], (NULL == stableFields) ? NULL : stableFields[idx]);
    }
    return 0.0;
}

UBOOL CAction::CalculateForceOnGauge(INT num, const CFieldGauge* const* gaugeFields, CFieldGauge* const* forceFields, CFieldGauge* const* stapleFields, ESolverPhase ePhase) const
{
    if (1 == m_byGaugeFieldIds.Num() && 0 == m_byBosonFieldIds.Num())
    {
        INT idx = GetGaugeFieldIndexById(num, gaugeFields, m_byGaugeFieldIds[0]);
        return CalculateForceOnGaugeSingleField(gaugeFields[idx], forceFields[idx], (NULL == stapleFields) ? NULL : stapleFields[idx], ePhase);
    }
    return FALSE;
}

void CAction::PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate)
{
    if (1 == m_byGaugeFieldIds.Num() && 0 == m_byBosonFieldIds.Num())
    {
        INT idx = GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
        PrepareForHMCSingleField(gaugeFields[idx], iUpdateIterate);
    }
}

CCString CAction::GetInfos(const CCString& tab) const
{
    CCString sRet = tab + _T("GaugeFields : ") + appAnyToString(m_byGaugeFieldIds) + _T("\n");
    sRet = sRet + tab + _T("BosonFields : ") + appAnyToString(m_byBosonFieldIds) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================