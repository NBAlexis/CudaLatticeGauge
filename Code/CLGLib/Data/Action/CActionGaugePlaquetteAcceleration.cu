//=============================================================================
// FILENAME : CActionGaugePlaquetteAcceleration.cu
// 
// DESCRIPTION:
// This is the class for rotating su3
//
// REVISION:
//  [mm/dd/yy]
//  [07/27/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionGaugePlaquetteAcceleration.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteAcceleration)


#pragma region kernels


/**
* Using plaqutte and (f(n)+f(n+mu)+f(n+nu)+f(n+mu+nu))/4 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermSU3_Acc(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
    Real betaOverN, Real fGsq,
#if !_CLG_DOUBLEFLOAT
    DOUBLE* results
#else
    Real* results
#endif
)
{
    intokernalInt4;

#if !_CLG_ASSUME_SQUARE_LATTICE
    UINT plaqLength = __idx->m_pSmallData[CIndexData::kPlaqLengthIdx];
    UINT plaqCountAllSite = __idx->m_pSmallData[CIndexData::kPlaqPerSiteIdx] * plaqLength;
#endif

    //Multi-GPU: the acceleration weight g*t^2 needs the GLOBAL t
    //(P4-1.4/1.4-R1); identity on single-GPU. sSite4 itself stays local for
    //the index-table lookups below.
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));
    fGsq = fGsq * sSite4G.w * sSite4G.w;
    //i=0: 12
    //  1: 13
    //  2: 14
    //  3: 23
    //  4: 24
    //  5: 34

    BYTE idx = 1;

    //Real resThisThread = F(0.0);

    //========================================
    //find plaqutte 1-3, or 2-3
    SIndex first = pCachedPlaqutte[idx * plaqLength + uiSiteIndex * plaqCountAllSite];
    deviceSU3 toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
    if (first.NeedToDagger())
    {
        toAdd.Dagger();
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[idx * plaqLength + j + uiSiteIndex * plaqCountAllSite];
        deviceSU3 toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
        if (first.NeedToDagger())
        {
            toAdd.MulDagger(toMul);
        }
        else
        {
            toAdd.Mul(toMul);
        }
    }
    results[uiSiteIndex] = static_cast<DOUBLE>(betaOverN * (F(3.0) - toAdd.ReTr()) * fGsq);

    idx = 3;
    first = pCachedPlaqutte[idx * plaqLength + uiSiteIndex * plaqCountAllSite];
    toAdd = _deviceGetGaugeBCT(byFieldId, pDeviceData, first);
    if (first.NeedToDagger())
    {
        toAdd.Dagger();
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[idx * plaqLength + j + uiSiteIndex * plaqCountAllSite];
        deviceSU3 toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
        if (first.NeedToDagger())
        {
            toAdd.MulDagger(toMul);
        }
        else
        {
            toAdd.Mul(toMul);
        }
    }
#if !_CLG_DOUBLEFLOAT
    results[uiSiteIndex] += static_cast<DOUBLE>(betaOverN * (F(3.0) - toAdd.ReTr()) * fGsq);
#else
    results[uiSiteIndex] += betaOverN * (F(3.0) - toAdd.ReTr()) * fGsq;
#endif

}


/**
* 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Chair_Acc(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
    Real betaOverN, Real fG,
#if !_CLG_DOUBLEFLOAT
    DOUBLE* results
#else
    Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = F(0.125) * betaOverN;

    //===============
    //V413
    const Real fV413 = _deviceChairTermT(pDeviceData, byFieldId, sSite4, 3, 0, 2, uiN);

    //===============
    //V423
    const Real fV423 = _deviceChairTermT(pDeviceData, byFieldId, sSite4, 3, 1, 2, uiN);

    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));
    results[uiSiteIndex] = -F(2.0) * (fV413 + fV423) * betaOverN * fG * sSite4G.w;
}

/**
* 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermSU3_Acc(
    const deviceSU3 * __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedIndex,
#if !_CLG_ASSUME_SQUARE_LATTICE
    BYTE plaqLength, BYTE plaqCountPerLink,
#endif
    deviceSU3 *pForceData,
    Real betaOverN, 
    Real fGSq,
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    //const SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    //Multi-GPU: the force factor (1 + g^2 t^2) needs the GLOBAL t (same as the
    //energy kernels above); take it once here instead of twice per staple.
    const SInt4 sSite4G = _deviceSIndexToGlobalInt4(__deviceSiteIndexToSIndex(uiSiteIndex));

    //Real test_force = F(0.0);
    betaOverN = betaOverN * F(-0.5);
#if !_CLG_ASSUME_SQUARE_LATTICE
    const UINT plaqLengthm1 = plaqLength - 1;
    const UINT plaqCountAllLink = plaqCountPerLink * plaqLengthm1;
#endif

    #pragma unroll
    for (UINT idir = 0; idir < 4; ++idir)
    {
        //if (__idx->_deviceIsBondOnSurface(uiBigIdx, idir))
        //{
        //    continue;
        //}

        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        deviceSU3 res = deviceSU3::makeSU3Zero();
        const BYTE mu = idir;
        //there are 6 staples,
        //3 'other directions' each is sum of two plaquttes
        for (int i = 0; i < plaqCountPerLink; ++i)
        {
            SIndex first = pCachedIndex[i * plaqLengthm1 + linkIndex * plaqCountAllLink];
            const BYTE nu = first.m_byDir;

            deviceSU3 toAdd(_deviceGetGaugeBCT(byFieldId, pDeviceData, first));
            Real fFactorG = F(1.0);

            if (
                   (0 == mu && 2 == nu)
                || (1 == mu && 2 == nu)
                || (2 == mu && 0 == nu)
                || (2 == mu && 1 == nu)
                )
            {
                fFactorG = fFactorG + fGSq * sSite4G.w * sSite4G.w;
            }

            if (first.NeedToDagger())
            {
                toAdd.Dagger();
            }

            for (BYTE j = 1; j < plaqLengthm1; ++j)
            {
                SIndex nextlink = pCachedIndex[i * plaqLengthm1 + j + linkIndex * plaqCountAllLink];
                deviceSU3 toMul(_deviceGetGaugeBCT(byFieldId, pDeviceData, nextlink));

                if (nextlink.NeedToDagger())
                {
                    toAdd.MulDagger(toMul);
                }
                else
                {
                    toAdd.Mul(toMul);
                }
            }

            toAdd.MulReal(fFactorG);
            res.Add(toAdd);
        }

        //staple calculated
        res.MulReal(betaOverN);

        //force is additive
        pForceData[linkIndex].Add(res);
    }
}


/**
* Split to 6 functions
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term413_Acc(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fG * F(0.125);

    //===============
    //add force for rho=3
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);
    deviceSU3 staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 0, 3, _deviceHi_Acc);

    staple_term.MulReal(betaOverN);
    pForceData[uiLink3].Add(staple_term);

    //===============
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);
    staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 0, 2, _deviceHi_Acc);

    staple_term.MulReal(betaOverN);
    pForceData[uiLink4].Add(staple_term);

    //===============
    //add force for nu=1
    const UINT uiLink1 = _deviceGetLinkIndex(uiSiteIndex, 0);
    staple_term = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 0, 2, _deviceHi_Acc);

    staple_term.MulReal(betaOverN);
    pForceData[uiLink1].Add(staple_term);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term423_Acc(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fG * F(0.125);

    //===============
    //g t V423
    //add force for rho=3
    const UINT uiLink3 = _deviceGetLinkIndex(uiSiteIndex, 2);
    deviceSU3 staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        2, 1, 3, _deviceHi_Acc);

    staple_term.MulReal(betaOverN);
    pForceData[uiLink3].Add(staple_term);

    //===============
    //add force for mu=4
    const UINT uiLink4 = _deviceGetLinkIndex(uiSiteIndex, 3);
    staple_term = _deviceStapleChairTerm1T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 1, 2, _deviceHi_Acc);

    staple_term.MulReal(betaOverN);
    pForceData[uiLink4].Add(staple_term);

    //===============
    //add force for nu=2
    const UINT uiLink2 = _deviceGetLinkIndex(uiSiteIndex, 1);

    staple_term = _deviceStapleChairTerm2T(byFieldId, pDeviceData, sSite4, uiSiteIndex, uiBigIdx,
        3, 1, 2, _deviceHi_Acc);

    staple_term.MulReal(betaOverN);
    pForceData[uiLink2].Add(staple_term);
}

#pragma endregion


CActionGaugePlaquetteAcceleration::CActionGaugePlaquetteAcceleration()
    : CAction()
    , m_uiPlaqutteCount(0)
{
}

void CActionGaugePlaquetteAcceleration::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        m_fLastEnergy = EnergySingleField(FALSE, pGauge, NULL);
    }
}

void CActionGaugePlaquetteAcceleration::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

    Real fG = 0.1f;
    param.FetchValueReal(_T("AccG"), fG);
    CCommonData::m_fG = fG;

}

UBOOL CActionGaugePlaquetteAcceleration::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    //pGauge->CalculateForceAndStaple(pForce, pStaple, m_fBetaOverN);

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    if (NULL == pGaugeSU3 || NULL == pForceSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3 now.\n"));
        return TRUE;
    }

    const Real fBeta = m_fBetaOverNR;
    const Real fG = CCommonData::m_fG;
    const Real fGSq = fG * fG;
    preparethread;

#if !_CLG_ASSUME_SQUARE_LATTICE
    _LAUNCH_KERNEL(_kernelAddForce4PlaqutteTermSU3_Acc, block, threads, 
        pGaugeSU3->m_pDeviceData, 
        appGetLattice()->m_pIndexCache->m_pStappleCache[pGaugeSU3->m_byFieldId],
        appGetLattice()->m_pIndexCache->m_uiPlaqutteLength,
        appGetLattice()->m_pIndexCache->m_uiPlaqutteCountPerLink,
        pForceSU3->m_pDeviceData, 
        fBeta,
        fGSq,
        pGauge->m_byFieldId);
#else
    _LAUNCH_KERNEL(_kernelAddForce4PlaqutteTermSU3_Acc, block, threads,
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pStappleCache[pGaugeSU3->m_byFieldId],
        pForceSU3->m_pDeviceData,
        fBeta,
        fGSq,
        pGauge->m_byFieldId);
#endif

    _LAUNCH_KERNEL(_kernelAddForceChairTermSU3_Term413_Acc, block, threads, pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, fBeta, fG);

    _LAUNCH_KERNEL(_kernelAddForceChairTermSU3_Term423_Acc, block, threads, pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, fBeta, fG);


    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

DOUBLE CActionGaugePlaquetteAcceleration::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStaple)
{
    if (bBeforeEvolution)
    {
        return m_fLastEnergy;
    }
    m_fNewEnergy = pGauge->CalculatePlaqutteEnergy(m_fBetaOverN);

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    if (NULL == pGaugeSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3-Dirichlet now.\n"));
        return m_fNewEnergy;
    }

    const Real fBeta = m_fBetaOverNR;
    const Real fG = CCommonData::m_fG;
    const Real fGSq = fG * fG;
    preparethread;

    appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

    _LAUNCH_KERNEL(_kernelAdd4PlaqutteTermSU3_Acc, block, threads, 
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData, 
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache[pGaugeSU3->m_byFieldId],
            fBeta,
            fGSq,
            _D_RealThreadBuffer);

    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    _LAUNCH_KERNEL(_kernelAddChairTermSU3_Chair_Acc, block, threads, 
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        fBeta,
        fG,
        _D_RealThreadBuffer);

     m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    //m_fNewEnergy is assembled from local partial sums (plaquette part plus
    //rotating-term ThreadBufferSums); sum the total across ranks for the
    //global action (P4-1.4). No-op on a lone rank.
    appGlobalSum(m_fNewEnergy);

    return m_fNewEnergy;
}

void CActionGaugePlaquetteAcceleration::SetG(Real fG)
{ 
    CCommonData::m_fG = fG;
}

CCString CActionGaugePlaquetteAcceleration::GetInfos(const CCString &tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("fG : ") + appToString(CCommonData::m_fG) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================