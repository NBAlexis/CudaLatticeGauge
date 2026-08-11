//=============================================================================
// FILENAME : CActionGaugePlaquetteAnisotropic.cpp
//
// DESCRIPTION:
// Anisotropic tree-level Symanzik gauge action
//
// REVISION:
//  [mm/dd/yy]
//  [07/25/2026 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionGaugePlaquetteAnisotropic.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteAnisotropic)

CActionGaugePlaquetteAnisotropic::CActionGaugePlaquetteAnisotropic()
    : CActionGaugePlaquette()
    , m_fXi(1.0)
{
}

void CActionGaugePlaquetteAnisotropic::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CActionGaugePlaquette::Initial(pOwner, param, byId);

    param.FetchValueDOUBLE(_T("Xi"), m_fXi);
    if (m_fXi < _CLG_FLT_EPSILON)
    {
        appCrucial(_T("CActionGaugePlaquetteAnisotropic: Xi (%f) must be positive!\n"), m_fXi);
        _FAIL_EXIT;
    }
}

void CActionGaugePlaquetteAnisotropic::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        if (abs(m_fXi - 1.0) < _CLG_FLT_EPSILON)
        {
            CActionGaugePlaquette::PrepareForHMCSingleField(pGauge, uiUpdateIterate);
            return;
        }

        //the field level kernels use the inverse convention, so 1/Xi is passed
        const DOUBLE fKernelXi = 1.0 / m_fXi;
        if (m_bCloverEnergy)
        {
            m_fLastEnergy = pGauge->CalculatePlaqutteEnergyUseCloverAnisotropy(m_fBetaOverN, fKernelXi);
        }
        else
        {
            m_fLastEnergy = pGauge->CalculatePlaqutteEnergyAnisotropy(m_fBetaOverN, fKernelXi);
        }
        //I9: the anisotropic branch must globalise the cached "before" energy
        //exactly like the base-class xi=1 branch -- CalculatePlaqutteEnergy*Anisotropy
        //sums only this rank's local sub-lattice, so without this every trajectory
        //sees a constant ΔH offset (= the other ranks' local action). No-op on one rank.
        appGlobalSum(m_fLastEnergy);
        checkCudaErrors(cudaDeviceSynchronize());
        checkCudaErrors(cudaGetLastError());
    }
}

UBOOL CActionGaugePlaquetteAnisotropic::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    if (abs(m_fXi - 1.0) < _CLG_FLT_EPSILON)
    {
        return CActionGaugePlaquette::CalculateForceOnGaugeSingleField(pGauge, pForce, pStaple, ePhase);
    }

    const DOUBLE fKernelXi = 1.0 / m_fXi;
    if (m_bCloverEnergy)
    {
        pGauge->CalculateForceCloverAnisotropy(pForce, m_fBetaOverN, fKernelXi);
    }
    else
    {
        pGauge->CalculateForceAnisotropy(pForce, m_fBetaOverN, fKernelXi);
    }
    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

/**
* The implementation depends on the type of gauge field
*/
DOUBLE CActionGaugePlaquetteAnisotropic::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple)
{
    if (abs(m_fXi - 1.0) < _CLG_FLT_EPSILON)
    {
        return CActionGaugePlaquette::EnergySingleField(bBeforeEvolution, pGauge, pStaple);
    }

    if (NULL != pStaple)
    {
        //there is no anisotropic staple-energy contract
        appCrucial(_T("CActionGaugePlaquetteAnisotropic: energy using staple is not supported in the anisotropic branch!\n"));
        _FAIL_EXIT;
    }

    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

    const DOUBLE fKernelXi = 1.0 / m_fXi;
    if (m_bCloverEnergy)
    {
        m_fNewEnergy = pGauge->CalculatePlaqutteEnergyUseCloverAnisotropy(m_fBetaOverN, fKernelXi);
    }
    else
    {
        m_fNewEnergy = pGauge->CalculatePlaqutteEnergyAnisotropy(m_fBetaOverN, fKernelXi);
    }
    //I9: local partial sum -> global action (see PrepareForHMCSingleField).
    appGlobalSum(m_fNewEnergy);
    return m_fNewEnergy;
}

CCString CActionGaugePlaquetteAnisotropic::GetInfos(const CCString &tab) const
{
    CCString sRet = CActionGaugePlaquette::GetInfos(tab);
    sRet = sRet + tab + _T("Xi (bare gauge anisotropy) : ") + appToString(m_fXi) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
