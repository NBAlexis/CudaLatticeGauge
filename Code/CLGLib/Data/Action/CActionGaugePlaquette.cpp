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
    , m_fLastEnergy(F(0.0))
    , m_fNewEnergy(F(0.0))
    , m_fBetaOverN(F(0.1))
{
}

void CActionGaugePlaquette::PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        m_fLastEnergy = pGauge->CalculatePlaqutteEnergy(m_fBetaOverN);
    }
}

void CActionGaugePlaquette::OnFinishTrajectory(UBOOL bAccepted)
{
    if (bAccepted)
    {
        m_fLastEnergy = m_fNewEnergy;
    }
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

UBOOL CActionGaugePlaquette::CalculateForceOnGauge(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple) const
{
    pGauge->CalculateForceAndStaple(pForce, pStaple, m_fBetaOverN);
    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
Real CActionGaugePlaquette::Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }
    if (NULL == pStable)
    {
        m_fNewEnergy = pGauge->CalculatePlaqutteEnergy(m_fBetaOverN);
    }
    else
    {
        m_fNewEnergy = pGauge->CalculatePlaqutteEnergyUsingStable(m_fBetaOverN, pStable);
    }
    return m_fNewEnergy;
}

Real CActionGaugePlaquette::GetEnergyPerPlaqutte() const
{
    return m_pOwner->m_pGaugeField->CalculatePlaqutteEnergy(m_fBetaOverN) / m_uiPlaqutteCount;
}

CCString CActionGaugePlaquette::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CActionGaugePlaquette\n");
    sRet = sRet + tab + _T("Beta : ") + appFloatToString(CCommonData::m_fBeta) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================