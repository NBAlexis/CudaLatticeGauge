//=============================================================================
// FILENAME : CActionGaugePlaquetteBoost.cu
// 
// DESCRIPTION:
// This is the class for rotating su3
//
// REVISION:
//  [07/27/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteBoost)


#pragma region kernels


/**
* Using plaqutte and (f(n)+f(n+mu)+f(n+nu)+f(n+mu+nu))/4 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermSU3_Boost(
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

    UINT plaqLength = __idx->m_pSmallData[CIndexData::kPlaqLengthIdx];
    UINT plaqCountAll = __idx->m_pSmallData[CIndexData::kPlaqPerSiteIdx] * plaqLength;
    
    //i=0: 12
    //  1: 13
    //  2: 14
    //  3: 23
    //  4: 24
    //  5: 34
    //0->2, 1->4
    BYTE idx = 2;

    //Real resThisThread = F(0.0);

    //========================================
    //find plaqutte 1-4, or 2-4
    SIndex first = pCachedPlaqutte[idx * plaqLength + uiSiteIndex * plaqCountAll];
    deviceSU3 toAdd(_deviceGetGaugeBCSU3(byFieldId, pDeviceData, first));
    if (first.NeedToDagger())
    {
        toAdd.Dagger();
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[idx * plaqLength + j + uiSiteIndex * plaqCountAll];
        deviceSU3 toMul(_deviceGetGaugeBCSU3(byFieldId, pDeviceData, first));
        if (first.NeedToDagger())
        {
            toAdd.MulDagger(toMul);
        }
        else
        {
            toAdd.Mul(toMul);
        }
    }
    results[uiSiteIndex] = static_cast<DOUBLE>(betaOverN * (3.0 - toAdd.ReTr()) * fGsq);

    idx = 4;
    first = pCachedPlaqutte[idx * plaqLength + uiSiteIndex * plaqCountAll];
    toAdd = _deviceGetGaugeBCSU3(byFieldId, pDeviceData, first);
    if (first.NeedToDagger())
    {
        toAdd.Dagger();
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[idx * plaqLength + j + uiSiteIndex * plaqCountAll];
        deviceSU3 toMul(_deviceGetGaugeBCSU3(byFieldId, pDeviceData, first));
        if (first.NeedToDagger())
        {
            toAdd.MulDagger(toMul);
        }
        else
        {
            toAdd.Mul(toMul);
        }
    }

    results[uiSiteIndex] += static_cast<DOUBLE>(betaOverN * (3.0 - toAdd.ReTr()) * fGsq);

}


/**
* 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Chair_Boost(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
    Real betaOverN, Real fG,
    Real* results)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = F(0.125) * betaOverN;

    //===============
    //V413
    const Real fV413 = _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 0, 2, uiN);

    //===============
    //V423
    const Real fV423 = _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 1, 2, uiN);

    results[uiSiteIndex] = (fV413 + fV423) * betaOverN * fG;
}

/**
* 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermSU3_Boost(
    const deviceSU3 * __restrict__ pDeviceData,
    deviceSU3 *pForceData,
    Real betaOverN, 
    Real fGSq,
    BYTE byFieldId)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(-0.5);
    #pragma unroll
    for (UINT idir = 0; idir < 3; ++idir)
    {
        if (2 == idir)
        {
            //We do not have z, BUT have T
            idir = 3;
            UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

            //mu = idir, nu = i = sum _1-3
            deviceSU3 stap(_deviceStapleTerm_Boost(byFieldId, pDeviceData, sSite4, uiBigIdx, idir, 0));
            stap.Add(_deviceStapleTerm_Boost(byFieldId, pDeviceData, sSite4, uiBigIdx, idir, 1));
            deviceSU3 force(pDeviceData[linkIndex]);
            force.MulDagger(stap);
            force.Ta();
            force.MulReal(betaOverN * fGSq);
            pForceData[linkIndex].Add(force);
        }
        else
        {
            UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

            //mu = idir, nu = 4, i = mu
            deviceSU3 stap(_deviceStapleTerm_Boost(byFieldId, pDeviceData, sSite4, uiBigIdx, idir, 3));
            deviceSU3 force(pDeviceData[linkIndex]);
            force.MulDagger(stap);
            force.Ta();
            force.MulReal(betaOverN * fGSq);
            pForceData[linkIndex].Add(force);
        }
    }
}


/**
* Split to 6 functions
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term413_3_Boost(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fG * F(0.125);

    //===============
    //g t V413
    //add force for rho=3
    UINT uiLink = _deviceGetLinkIndex(uiSiteIndex, 2);
    const deviceSU3 staple_term = _deviceStapleChairTerm1_Boost(byFieldId, pDeviceData, 
        sSite4, uiBigIdx,
        2, 0, 3);
    deviceSU3 force(pDeviceData[uiLink]);
    force.MulDagger(staple_term);
    force.Ta();
    force.MulReal(betaOverN);
    pForceData[uiLink].Add(force);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term413_4_Boost(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
    deviceSU3 *pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fG * F(0.125);

    //===============
    //g t V413
    //add force for mu=4
    UINT uiLink = _deviceGetLinkIndex(uiSiteIndex, 3);
    const deviceSU3 staple_term = _deviceStapleChairTerm1_Boost(byFieldId, pDeviceData, 
        sSite4, uiBigIdx,
        3, 0, 2);
    deviceSU3 force(pDeviceData[uiLink]);
    force.MulDagger(staple_term);
    force.Ta();
    force.MulReal(betaOverN);
    pForceData[uiLink].Add(force);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term413_1_Boost(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
    deviceSU3 *pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fG * F(0.125);

    //===============
    //+x Omega V413
    //add force for nu=1
    UINT uiLink = _deviceGetLinkIndex(uiSiteIndex, 0);

    const deviceSU3 staple_term = _deviceStapleChairTerm2_Boost(byFieldId, pDeviceData, 
        sSite4, uiBigIdx,
        3, 0, 2);
    deviceSU3 force(pDeviceData[uiLink]);
    force.MulDagger(staple_term);
    force.Ta();
    force.MulReal(betaOverN);
    pForceData[uiLink].Add(force);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term423_3_Boost(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fG * F(0.125);

    //===============
    //g t V413
    //add force for rho=3
    UINT uiLink = _deviceGetLinkIndex(uiSiteIndex, 2);
    const deviceSU3 staple_term = _deviceStapleChairTerm1_Boost(byFieldId, pDeviceData, 
        sSite4, uiBigIdx,
        2, 1, 3);
    deviceSU3 force(pDeviceData[uiLink]);
    force.MulDagger(staple_term);
    force.Ta();
    force.MulReal(betaOverN);
    pForceData[uiLink].Add(force);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term423_4_Boost(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fG * F(0.125);

    //===============
    //g t V413
    //add force for mu=4
    UINT uiLink = _deviceGetLinkIndex(uiSiteIndex, 3);
    const deviceSU3 staple_term = _deviceStapleChairTerm1_Boost(byFieldId, pDeviceData, 
        sSite4, uiBigIdx,
        3, 1, 2);
    deviceSU3 force(pDeviceData[uiLink]);
    force.MulDagger(staple_term);
    force.Ta();
    force.MulReal(betaOverN);
    pForceData[uiLink].Add(force);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term423_2_Boost(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * F(0.5) * fG * F(0.125);

    //===============
    //+x Omega V413
    //add force for nu=1
    UINT uiLink = _deviceGetLinkIndex(uiSiteIndex, 1);

    const deviceSU3 staple_term = _deviceStapleChairTerm2_Boost(byFieldId, pDeviceData, 
        sSite4, uiBigIdx,
        3, 1, 2);
    deviceSU3 force(pDeviceData[uiLink]);
    force.MulDagger(staple_term);
    force.Ta();
    force.MulReal(betaOverN);
    pForceData[uiLink].Add(force);
}

#pragma endregion


CActionGaugePlaquetteBoost::CActionGaugePlaquetteBoost()
    : CAction()
    , m_uiPlaqutteCount(0)
{
}

void CActionGaugePlaquetteBoost::PrepareForHMCSingleField(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        m_fLastEnergy = EnergySingleField(FALSE, pGauge, NULL);
    }
}

void CActionGaugePlaquetteBoost::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
{
    CAction::Initial(pOwner, param, byId);

    m_fBetaOverN = CCommonData::m_fBeta / static_cast<Real>(GetDefaultMatrixN());
    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

    Real fG = 0.1f;
    param.FetchValueReal(_T("Boost"), fG);
    CCommonData::m_fG = fG;

    //TArray<INT> centerArray;
    //param.FetchValueArrayINT(_T("Center"), centerArray);
    //if (centerArray.Num() > 3)
    //{
    //    SSmallInt4 sCenter;
    //    sCenter.x = static_cast<SBYTE>(centerArray[0]);
    //    sCenter.y = static_cast<SBYTE>(centerArray[1]);
    //    sCenter.z = static_cast<SBYTE>(centerArray[2]);
    //    sCenter.w = static_cast<SBYTE>(centerArray[3]);
    //    CCommonData::m_sCenter = sCenter;
    //}
    //else
    //{
    //    CCommonData::m_sCenter.w = 0;
    //}
}

void CActionGaugePlaquetteBoost::SetBeta(Real fBeta)
{
    CCommonData::m_fBeta = static_cast<DOUBLE>(fBeta);
    m_fBetaOverN = fBeta / static_cast<Real>(GetDefaultMatrixN());
}

UBOOL CActionGaugePlaquetteBoost::CalculateForceOnGaugeSingleField(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    pGauge->CalculateForceAndStaple(pForce, pStaple, m_fBetaOverNR);

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    if (NULL == pGaugeSU3 || NULL == pForceSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3 now.\n"));
        return TRUE;
    }

    preparethread;

    _kernelAddForce4PlaqutteTermSU3_Boost << <block, threads >> >(
        pGaugeSU3->m_pDeviceData, 
        pForceSU3->m_pDeviceData, 
        m_fBetaOverNR, 
        CCommonData::m_fG * CCommonData::m_fG,
        pGauge->m_byFieldId);

    //_kernelAddForceChairTermSU3_Term413_3_Boost << <block, threads >> >(pGaugeSU3->m_pDeviceData,
    //    pForceSU3->m_pDeviceData, m_fBetaOverN, CCommonData::m_fG);

    //_kernelAddForceChairTermSU3_Term413_4_Boost << <block, threads >> >(pGaugeSU3->m_pDeviceData,
    //    pForceSU3->m_pDeviceData, m_fBetaOverN, CCommonData::m_fG);

    //_kernelAddForceChairTermSU3_Term413_1_Boost << <block, threads >> >(pGaugeSU3->m_pDeviceData,
    //    pForceSU3->m_pDeviceData, m_fBetaOverN, CCommonData::m_fG);

    //_kernelAddForceChairTermSU3_Term423_3_Boost << <block, threads >> >(pGaugeSU3->m_pDeviceData,
    //    pForceSU3->m_pDeviceData, m_fBetaOverN, CCommonData::m_fG);

    //_kernelAddForceChairTermSU3_Term423_4_Boost << <block, threads >> >(pGaugeSU3->m_pDeviceData,
    //    pForceSU3->m_pDeviceData, m_fBetaOverN, CCommonData::m_fG);

    //_kernelAddForceChairTermSU3_Term423_2_Boost << <block, threads >> >(pGaugeSU3->m_pDeviceData,
    //    pForceSU3->m_pDeviceData, m_fBetaOverN, CCommonData::m_fG);

    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

DOUBLE CActionGaugePlaquetteBoost::EnergySingleField(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
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

    preparethread;

    appGetCudaHelper()->ThreadBufferZero(_D_RealThreadBuffer);

    _kernelAdd4PlaqutteTermSU3_Boost << <block, threads >> > (
            pGaugeSU3->m_byFieldId,
            pGaugeSU3->m_pDeviceData, 
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
            m_fBetaOverNR,
            CCommonData::m_fG * CCommonData::m_fG,
            _D_RealThreadBuffer);

    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);


    //_kernelAddChairTermSU3_Chair_Boost << <block, threads >> > (
    //    pGaugeSU3->m_pDeviceData,
    //    m_fBetaOverN,
    //    CCommonData::m_fG,
    //    _D_RealThreadBuffer);

    // m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    return m_fNewEnergy;
}

void CActionGaugePlaquetteBoost::SetG(Real fG)
{ 
    CCommonData::m_fG = fG;
}

CCString CActionGaugePlaquetteBoost::GetInfos(const CCString &tab) const
{
    CCString sRet = CAction::GetInfos(tab);
    sRet = sRet + tab + _T("Beta : ") + appToString(CCommonData::m_fBeta) + _T("\n");
    sRet = sRet + tab + _T("Boost : ") + appToString(CCommonData::m_fG) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================