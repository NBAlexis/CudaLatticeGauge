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
        m_fBetaOverN = fBeta / static_cast<DOUBLE>(GetDefaultMatrixN());
    }
}

DOUBLE CAction::Energy(UBOOL bBeforeEvolution, INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, const CFieldGauge* const* stableFields)
{
    if (1 == m_byGaugeFieldIds.Num() && 0 == m_byBosonFieldIds.Num())
    {
        INT idx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
        return EnergySingleField(bBeforeEvolution, gaugeFields[idx], (NULL == stableFields) ? NULL : stableFields[idx]);
    }
    appCrucial(_T("Energy not implemented!\n"));
    return 0.0;
}

UBOOL CAction::CalculateForce(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields,
    CFieldGauge* const* gaugeForces, CFieldBoson* const* bosonForces,
    CFieldGauge* const* stapleFields, ESolverPhase ePhase) const
{
    if (1 == m_byGaugeFieldIds.Num() && 0 == m_byBosonFieldIds.Num())
    {
        INT idx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
        return CalculateForceOnGaugeSingleField(gaugeFields[idx], gaugeForces[idx], (NULL == stapleFields) ? NULL : stapleFields[idx], ePhase);
    }
    appCrucial(_T("CalculateForce not implemented!\n"));
    return FALSE;
}

void CAction::PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* bosonFields, UINT iUpdateIterate)
{
    if (1 == m_byGaugeFieldIds.Num() && 0 == m_byBosonFieldIds.Num())
    {
        INT idx = CLatticeData::GetGaugeFieldIndexById(gaugeNum, gaugeFields, m_byGaugeFieldIds[0]);
        PrepareForHMCSingleField(gaugeFields[idx], iUpdateIterate);
        return;
    }
    appCrucial(_T("CAction PrepareForHMC not implemented!\n"));
}

UINT CAction::GetDefaultMatrixN() const
{
    if (m_byGaugeFieldIds.Num() > 0)
    {
        const CFieldGauge* gauge = dynamic_cast<const CFieldGauge*>(appGetLattice()->GetFieldById(m_byGaugeFieldIds[0]));
        if (NULL != gauge)
        {
            return gauge->MatrixN();
        }
    }

    return appGetLattice()->GetDefaultSUN();
}

CCString CAction::GetInfos(const CCString& tab) const
{
    CCString sRet = CBase::GetInfos(tab);
    sRet = sRet + tab + _T("ActionId : ") + appToString(m_byActionId) + _T("\n");
    sRet = sRet + tab + _T("Nc : ") + appToString(GetDefaultMatrixN()) + _T("\n");
    sRet = sRet + tab + _T("Beta / Nc : ") + appToString(m_fBetaOverN) + _T("\n");
    sRet = sRet + tab + _T("GaugeFields : ") + appToString(m_byGaugeFieldIds) + _T("\n");
    sRet = sRet + tab + _T("BosonFields : ") + appToString(m_byBosonFieldIds) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================