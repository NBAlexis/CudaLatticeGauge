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
#include "CActionGaugePlaquette.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquette)

CActionGaugePlaquette::CActionGaugePlaquette()
    : CAction()
    , m_bCloverEnergy(FALSE)
    , m_uiPlaqutteCount(0)
{
}

void CActionGaugePlaquette::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        if (m_bCloverEnergy)
        {
            m_fLastEnergy = pGauge->CalculatePlaqutteEnergyUseClover(m_fBetaOverN);
        }
        else
        {
            m_fLastEnergy = pGauge->CalculatePlaqutteEnergy(m_fBetaOverN);
        }
    }
}

void CActionGaugePlaquette::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    m_fBetaOverN = CCommonData::m_fBeta / static_cast<DOUBLE>(GetDefaultMatrixN());
    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

    INT iUsing4Plaq = 0;
    if (param.FetchValueINT(_T("CloverEnergy"), iUsing4Plaq))
    {
        if (1 == iUsing4Plaq)
        {
            m_bCloverEnergy = TRUE;
        }
    }
}
void CActionGaugePlaquette::SetBeta(DOUBLE fBeta)
{
    CCommonData::m_fBeta = fBeta;
    m_fBetaOverN = fBeta / static_cast<DOUBLE>(GetDefaultMatrixN());
}

UBOOL CActionGaugePlaquette::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    pGauge->CalculateForceAndStaple(pForce, pStaple, static_cast<Real>(m_fBetaOverN));
    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
DOUBLE CActionGaugePlaquette::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

    if (NULL == pStable)
    {
        if (m_bCloverEnergy)
        {
            m_fNewEnergy = pGauge->CalculatePlaqutteEnergyUseClover(m_fBetaOverN);
        }
        else
        {
            m_fNewEnergy = pGauge->CalculatePlaqutteEnergy(m_fBetaOverN);
        }
    }
    else
    {
        m_fNewEnergy = pGauge->CalculatePlaqutteEnergyUsingStable(m_fBetaOverN, pStable);
    }
    return m_fNewEnergy;
}

//Real CActionGaugePlaquette::GetEnergyPerPlaqutte() const
//{
//    return m_pOwner->m_pGaugeField->CalculatePlaqutteEnergy(m_fBetaOverN) / m_uiPlaqutteCount;
//}

CCString CActionGaugePlaquette::GetInfos(const CCString &tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("Beta : ") + appToString(CCommonData::m_fBeta) + _T("\n");
    sRet = sRet + tab + _T("Clover : ") + appToString(m_bCloverEnergy) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================