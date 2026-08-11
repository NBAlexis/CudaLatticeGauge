//=============================================================================
// FILENAME : CGaugeSmearingHISQWithPhase.cpp
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [mm/dd/yy]
//  [12/30/2024 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "Data/Field/Gauge/CFieldGaugeU1Real.h"
#include "CGaugeSmearingASQTAD.h"
#include "CGaugeSmearingHISQ.h"
#include "CGaugeSmearingHISQWithPhase.h"

__BEGIN_NAMESPACE

#pragma endregion

template<typename gaugetype, INT matrixN>
void CGaugeSmearingHISQWithPhase<gaugetype, matrixN>::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    CGaugeSmearingHISQ<gaugetype, matrixN>::Initial(pOwner, params);

    INT iVaule = 0;
    params.FetchValueINT(_T("U1FieldId"), iVaule);
    m_byU1FieldId = static_cast<BYTE>(iVaule);

    params.FetchValueReal(_T("Charge"), m_fCharge);
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingHISQWithPhase<gaugetype, matrixN>::GaugeSmearing(CFieldGauge* pGauge, const class CFieldGauge* pOrignalGauge, CFieldGauge* pStaple, UBOOL bProject)
{
    const CFieldGaugeU1Real* phaseField = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byU1FieldId));
    if (NULL != phaseField)
    {
        pGauge->ApplyPhaseR(phaseField, m_fCharge);
    }
    else
    {
        appDetailed(_T("Phase not applied because u1 field not found\n"));
    }

    if (NULL == m_pOrignalGauge)
    {
        m_pOrignalGauge = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledCopy(pGauge, _T(__FILE__), __LINE__));
    }
    else
    {
        pGauge->CopyTo(m_pOrignalGauge);
    }

    //We need a 'orignal gauge' to preserve orignal gauge after applied phase
    CGaugeSmearingHISQ<gaugetype, matrixN>::GaugeSmearing(pGauge, m_pOrignalGauge, pStaple, bProject);
}

template<typename gaugetype, INT matrixN>
void CGaugeSmearingHISQWithPhase<gaugetype, matrixN>::DerivateOnU(const CFieldGauge* pEffectiveGauge, const class CFieldGauge* pOrignalGauge, const class CFieldGauge* pNaikForce, CFieldGauge* pf0) const
{
    //const CFieldGaugeU1Real* phaseField = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byU1FieldId));
    //if (NULL != phaseField)
    //{
    //    pf0->ApplyPhaseR(phaseField, -m_fCharge);
    //}
    CGaugeSmearingHISQ<gaugetype, matrixN>::DerivateOnU(pEffectiveGauge, m_pOrignalGauge, pNaikForce, pf0);
    //const CFieldGaugeU1Real* phaseField = dynamic_cast<const CFieldGaugeU1Real*>(appGetLattice()->GetFieldById(m_byU1FieldId));
    //if (NULL != phaseField)
    //{
    //    pf0->ApplyPhaseR(phaseField, -m_fCharge);
    //}
}

template<typename gaugetype, INT matrixN>
CCString CGaugeSmearingHISQWithPhase<gaugetype, matrixN>::GetInfos(const CCString &tab) const
{
    CCString sRet = CGaugeSmearingHISQ<gaugetype, matrixN>::GetInfos(tab);
    sRet = sRet + tab + _T("UIFieldID   : ") + appToString(m_byU1FieldId) + _T("\n");
    return sRet;
}

__CLGIMPLEMENT_CLASS(CGaugeSmearingHISQWithPhaseSU3)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================