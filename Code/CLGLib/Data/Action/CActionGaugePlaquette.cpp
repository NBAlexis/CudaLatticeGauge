//=============================================================================
// FILENAME : CActionGaugePlaquette.cu
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [12/4/2018 nbale]
//=============================================================================
#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquette)

CActionGaugePlaquette::CActionGaugePlaquette()
    : CAction()
    , m_uiPlaqutteCount(0)
{
}

void CActionGaugePlaquette::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    m_pOwner = pOwner;
    m_byActionId = byId;
    Real fBeta = 0.1f;
    param.FetchValueReal(_T("Beta"), fBeta);
    CCommonData::m_fBeta = fBeta;
    if (NULL != pOwner->m_pGaugeField && EFT_GaugeSU3 == pOwner->m_pGaugeField->GetFieldType())
    {
        fBeta = fBeta / F(3.0);
    }
    m_fBetaOverN = fBeta;
    m_uiPlaqutteCount = _HC_Volumn * (_HC_Dir - 1) * (_HC_Dir - 2);
}

void CActionGaugePlaquette::SetBeta(Real fBeta)
{
    CCommonData::m_fBeta = fBeta;
    if (NULL != m_pOwner->m_pGaugeField && EFT_GaugeSU3 == m_pOwner->m_pGaugeField->GetFieldType())
    {
        fBeta = fBeta / F(3.0);
    }
    m_fBetaOverN = fBeta;
}

void CActionGaugePlaquette::CalculateForceOnGauge(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple) const
{
    pGauge->CalculateForceAndStaple(pForce, pStaple, m_fBetaOverN);
    checkCudaErrors(cudaDeviceSynchronize());
}

/**
* The implementation depends on the type of gauge field
*/
Real CActionGaugePlaquette::Energy(const class CFieldGauge* pGauge) const
{
    return pGauge->CalculatePlaqutteEnergy(m_fBetaOverN);
}

Real CActionGaugePlaquette::Energy(const class CFieldGauge* pGauge, const class CFieldGauge* pStable) const
{
    return pGauge->CalculatePlaqutteEnergyUsingStable(m_fBetaOverN, pStable);
}

Real CActionGaugePlaquette::GetEnergyPerPlaqutte() const
{
    return m_pOwner->m_pGaugeField->CalculatePlaqutteEnergy(m_fBetaOverN) / m_uiPlaqutteCount;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================