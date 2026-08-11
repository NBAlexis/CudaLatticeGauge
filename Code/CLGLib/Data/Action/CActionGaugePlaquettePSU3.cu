//=============================================================================
// FILENAME : CActionGaugePlaquettePSU3.cu
// 
// DESCRIPTION:
// Implementation of the PSU(N) adjoint plaquette action for SU2 and SU3.
//
// Dynamic action:
//   S_dyn = -beta_tilde * sum_{n,mu<nu} |Tr(U_{mu nu}(n))|^2
//
// Force contribution from a single surrounding plaquette P:
//   F_P = beta_tilde * [Tr(U_P^dag) U_P - Tr(U_P) U_P^dag]
//
// The trace must be taken per plaquette before summation. The resulting matrix
// is analytically traceless anti-Hermitian; a numerical TA projection is applied
// at the end of the kernel.
//
// REVISION:
//  [07/02/26]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionGaugePlaquettePSU3.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquettePSU3)

#pragma region kernels

/**
* Energy kernel: S_dyn = -beta_tilde * sum_P |Tr(U_P)|^2
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelEnergy_PSU3(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerSite,
#endif
    Real beta_tilde,
    DOUBLE* results
)
{
    intokernalInt4;

    DOUBLE resThisThread = 0.0;

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqCountAllSite = plaqCountPerSite * plaqLength;
#endif

    for (BYTE i = 0; i < plaqCountPerSite; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLength + uiSiteIndex * plaqCountAllSite];
        deviceGauge toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

        if (first.NeedToDagger())
        {
            toAdd.Dagger();
        }

        for (BYTE j = 1; j < plaqLength; ++j)
        {
            first = pCachedIndex[i * plaqLength + j + uiSiteIndex * plaqCountAllSite];
            deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
            if (first.NeedToDagger())
            {
                toAdd.MulDagger(toMul);
            }
            else
            {
                toAdd.Mul(toMul);
            }
        }

        const CLGComplex tr = toAdd.Tr();
        const Real tr_sq = tr.x * tr.x + tr.y * tr.y;
        resThisThread -= static_cast<DOUBLE>(tr_sq);
    }

    results[uiSiteIndex] = resThisThread * static_cast<DOUBLE>(beta_tilde);
}

/**
* Force kernel.
*
* For each link, loop over the surrounding staples and accumulate
*   Y = sum_P Tr(U_P) staple
* where U_P = U_mu(n) * staple^dag is the plaquette.
* The integrator will later apply U * Y^dag + TA projection, which yields
* the final traceless anti-Hermitian force
*   F = -beta_tilde/2 * sum_P [Tr(U_P^dag) U_P - Tr(U_P) U_P^dag].
*/
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelForce_PSU3(
    BYTE byFieldId,
    const deviceGauge* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceGauge* pForceData,
    Real beta_tilde)
{
    intokernalDir_NoDir;

#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1;
#endif

    deviceGauge res = _makeZero<deviceGauge>();
    const deviceGauge u = pDeviceData[uiLinkIndex];

    for (BYTE i = 0; i < plaqCountPerLink; ++i)
    {
        SIndex first = pCachedIndex[i * plaqLengthm1 + uiLinkIndex * plaqCountAllLink];
        deviceGauge staple(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));

        if (first.NeedToDagger())
        {
            staple.Dagger();
        }

        for (BYTE j = 1; j < plaqLengthm1; ++j)
        {
            SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + uiLinkIndex * plaqCountAllLink];
            deviceGauge toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink));
            if (nextlink.NeedToDagger())
            {
                staple.MulDagger(toMul);
            }
            else
            {
                staple.Mul(toMul);
            }
        }

        // U_P = U * staple^dag
        const deviceGauge u_p = u.MulDaggerC(staple);
        const CLGComplex tr = u_p.Tr();

        // Accumulate Tr(U_P) * staple.  The integrator applies LeftMul+TA.
        const deviceGauge term = staple.MulCompC(tr);
        res.Add(term);
    }

    res.MulReal(-beta_tilde);

    // Add to the shared force field. The integrator zeroes it once and every
    // action accumulates its raw contribution; assignment would erase the
    // forces produced by earlier actions.
    _add(pForceData[uiLinkIndex], res);
}

#pragma endregion


CActionGaugePlaquettePSU3::CActionGaugePlaquettePSU3()
    : CAction()
    , m_uiPlaqutteCount(0)
{
}

void CActionGaugePlaquettePSU3::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        m_fLastEnergy = EnergySingleField(FALSE, pGauge, NULL);
    }
}

void CActionGaugePlaquettePSU3::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

    // CAction::Initial stores Beta / Nc in m_fBetaOverN.
    // The adjoint action uses beta_tilde = Beta / (N^2 - 1):
    //   SU(3) -> Beta / 8,  SU(2) -> Beta / 3.
    const DOUBLE matrixN = static_cast<DOUBLE>(GetDefaultMatrixN());
    if (matrixN < 1.5 || matrixN > 3.5)
    {
        appCrucial(_T("CActionGaugePlaquettePSU3 only supports SU2 or SU3 gauge fields.\n"));
    }
    m_fBetaOverN = m_fBetaOverN * matrixN / (matrixN * matrixN - 1.0);
}

UBOOL CActionGaugePlaquettePSU3::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    if (NULL != pGaugeSU3 && NULL != pForceSU3)
    {
        preparethreadDir;

#if !_CLG_ASSUME_SQUARE_LATTICE
        _LAUNCH_KERNEL(_kernelForce_PSU3<deviceSU3>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pStappleCache[pGaugeSU3->m_byFieldId],
            appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
            appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
            pForceSU3->m_pDeviceData,
            m_fBetaOverNR);
#else
        _LAUNCH_KERNEL(_kernelForce_PSU3<deviceSU3>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pStappleCache[pGaugeSU3->m_byFieldId],
            pForceSU3->m_pDeviceData,
            m_fBetaOverNR);
#endif
        checkCudaErrors(cudaDeviceSynchronize());
        return TRUE;
    }

    const CFieldGaugeSU2* pGaugeSU2 = dynamic_cast<const CFieldGaugeSU2*>(pGauge);
    CFieldGaugeSU2* pForceSU2 = dynamic_cast<CFieldGaugeSU2*>(pForce);
    if (NULL != pGaugeSU2 && NULL != pForceSU2)
    {
        preparethreadDir;

#if !_CLG_ASSUME_SQUARE_LATTICE
        _LAUNCH_KERNEL(_kernelForce_PSU3<deviceSU2>, block, threads,
            pGaugeSU2->m_byFieldId,
            pGaugeSU2->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pStappleCache[pGaugeSU2->m_byFieldId],
            appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
            appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
            pForceSU2->m_pDeviceData,
            m_fBetaOverNR);
#else
        _LAUNCH_KERNEL(_kernelForce_PSU3<deviceSU2>, block, threads,
            pGaugeSU2->m_byFieldId,
            pGaugeSU2->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pStappleCache[pGaugeSU2->m_byFieldId],
            pForceSU2->m_pDeviceData,
            m_fBetaOverNR);
#endif
        checkCudaErrors(cudaDeviceSynchronize());
        return TRUE;
    }

    appCrucial(_T("CActionGaugePlaquettePSU3 only works with SU2 or SU3 gauge fields.\n"));
    return FALSE;
}

DOUBLE CActionGaugePlaquettePSU3::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    if (NULL != pGaugeSU3)
    {
        preparethread;
        appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

#if !_CLG_ASSUME_SQUARE_LATTICE
        _LAUNCH_KERNEL(_kernelEnergy_PSU3<deviceSU3>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGaugeSU3->m_byFieldId],
            appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
            appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
            m_fBetaOverNR,
            _D_RealThreadBuffer);
#else
        _LAUNCH_KERNEL(_kernelEnergy_PSU3<deviceSU3>, block, threads,
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGaugeSU3->m_byFieldId],
            m_fBetaOverNR,
            _D_RealThreadBuffer);
#endif
        m_fNewEnergy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
        return m_fNewEnergy;
    }

    const CFieldGaugeSU2* pGaugeSU2 = dynamic_cast<const CFieldGaugeSU2*>(pGauge);
    if (NULL != pGaugeSU2)
    {
        preparethread;
        appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

#if !_CLG_ASSUME_SQUARE_LATTICE
        _LAUNCH_KERNEL(_kernelEnergy_PSU3<deviceSU2>, block, threads,
            pGaugeSU2->m_byFieldId,
            pGaugeSU2->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGaugeSU2->m_byFieldId],
            appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
            appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerSite,
            m_fBetaOverNR,
            _D_RealThreadBuffer);
#else
        _LAUNCH_KERNEL(_kernelEnergy_PSU3<deviceSU2>, block, threads,
            pGaugeSU2->m_byFieldId,
            pGaugeSU2->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGaugeSU2->m_byFieldId],
            m_fBetaOverNR,
            _D_RealThreadBuffer);
#endif
        m_fNewEnergy = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
        return m_fNewEnergy;
    }

    appCrucial(_T("CActionGaugePlaquettePSU3 only works with SU2 or SU3 gauge fields.\n"));
    return 0.0;
}

CCString CActionGaugePlaquettePSU3::GetInfos(const CCString &tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("PSU3 Adjoint Plaquette Action\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================
