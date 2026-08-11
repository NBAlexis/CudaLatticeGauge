//=============================================================================
// FILENAME : CActionGaugePlaquette.cu
// 
// DESCRIPTION:
// This is the class for all fields, gauge, fermion and spin fields are inherent from it
//
// REVISION:
//  [mm/dd/yy]
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
        //CRITICAL: this cached "before" energy is returned verbatim by
        //EnergySingleField(bBeforeEvolution=TRUE). The "after" energy there is
        //Allreduced to the global action, so this MUST be globalised the same way
        //or every trajectory sees a constant ΔH offset (= the other ranks' local
        //action) and rejects. CalculatePlaqutteEnergy* sums only this rank's local
        //sub-lattice; sum across ranks for the global action. No-op on one rank.
        appGlobalSum(m_fLastEnergy);
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaGetLastError());
    }
}

void CActionGaugePlaquette::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    //m_fBetaOverN = m_fBetaOverN / static_cast<DOUBLE>(GetDefaultMatrixN());
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

UBOOL CActionGaugePlaquette::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    //CFieldGauge* pForceOfThisAction = dynamic_cast<CFieldGauge*>(appGetLattice()->GetPooledFieldById(pGauge->m_byFieldId));
    if (m_bCloverEnergy)
    {
        pGauge->CalculateForceAndStapleClover(pForce, pStaple, static_cast<Real>(m_fBetaOverN));
    }
    else
    {
        pGauge->CalculateForceAndStaple(pForce, pStaple, static_cast<Real>(m_fBetaOverN));
    }
    
    //pForceOfThisAction->LeftMul(pGauge);
    //pForceOfThisAction->TA();
    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
DOUBLE CActionGaugePlaquette::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

    if (m_bCloverEnergy)
    {
        m_fNewEnergy = pGauge->CalculatePlaqutteEnergyUseClover(m_fBetaOverN);
    }
    else
    {
        if (NULL == pStaple)
        {
            m_fNewEnergy = pGauge->CalculatePlaqutteEnergy(m_fBetaOverN);
        }
        else
        {
            m_fNewEnergy = pGauge->CalculatePlaqutteEnergyUsingStaple(m_fBetaOverN, pStaple);
        }
    }
    //CalculatePlaqutteEnergy* sums only over this rank's local sub-lattice (each
    //plaquette is counted once, at its owned base site), so the global action
    //energy is the sum across ranks. No-op on a single rank. This is the HMC
    //energy used by the Metropolis test -- it must be identical on every rank.
    appGlobalSum(m_fNewEnergy);
    return m_fNewEnergy;
}

//Real CActionGaugePlaquette::GetEnergyPerPlaqutte() const
//{
//    return m_pOwner->m_pGaugeField->CalculatePlaqutteEnergy(m_fBetaOverN) / m_uiPlaqutteCount;
//}

CCString CActionGaugePlaquette::GetInfos(const CCString &tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("Clover : ") + appToString(m_bCloverEnergy) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================