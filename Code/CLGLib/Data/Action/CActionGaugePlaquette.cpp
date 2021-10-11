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
    , m_bCloverEnergy(FALSE)
    , m_fLastEnergy(F(0.0))
    , m_fNewEnergy(F(0.0))
    , m_fBetaOverN(F(0.1))
    , m_uiPlaqutteCount(0)
{
}

void CActionGaugePlaquette::PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate)
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
#if !_CLG_DOUBLEFLOAT
    DOUBLE fBeta = 0.1;
    param.FetchValueDOUBLE(_T("Beta"), fBeta);
#else
    Real fBeta = F(0.1);
    param.FetchValueReal(_T("Beta"), fBeta);
#endif
    CCommonData::m_fBeta = fBeta;
    m_fBetaOverN = fBeta / static_cast<DOUBLE>(_HC_SUN);
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
#if !_CLG_DOUBLEFLOAT
void CActionGaugePlaquette::SetBeta(DOUBLE fBeta)
#else
void CActionGaugePlaquette::SetBeta(Real fBeta)
#endif
{
    CCommonData::m_fBeta = fBeta;
    m_fBetaOverN = fBeta / static_cast<DOUBLE>(_HC_SUN);
}

UBOOL CActionGaugePlaquette::CalculateForceOnGauge(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
#if !_CLG_DOUBLEFLOAT
    pGauge->CalculateForceAndStaple(pForce, pStaple, static_cast<Real>(m_fBetaOverN));
#else
    pGauge->CalculateForceAndStaple(pForce, pStaple, m_fBetaOverN);
#endif
    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
#if !_CLG_DOUBLEFLOAT
DOUBLE CActionGaugePlaquette::Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
#else
Real CActionGaugePlaquette::Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
#endif
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
    CCString sRet;
    sRet = tab + _T("Name : CActionGaugePlaquette\n");
    sRet = sRet + tab + _T("Beta : ") + appFloatToString(CCommonData::m_fBeta) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================