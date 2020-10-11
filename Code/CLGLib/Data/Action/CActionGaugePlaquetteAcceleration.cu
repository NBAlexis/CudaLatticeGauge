//=============================================================================
// FILENAME : CActionGaugePlaquetteAcceleration.cu
// 
// DESCRIPTION:
// This is the class for rotating su3
//
// REVISION:
//  [07/27/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionGaugePlaquetteAcceleration)


#pragma region kernels


/**
* Using plaqutte and (f(n)+f(n+mu)+f(n+nu)+f(n+mu+nu))/4 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermSU3_Acc(
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
    //intokernalInt4;
    SSmallInt4 sSite4;
    UINT _ixy = (threadIdx.x + blockIdx.x * blockDim.x);
    UINT _iz_idx = (threadIdx.y + blockIdx.y * blockDim.y);

    sSite4.x = static_cast<SBYTE> (_ixy / _DC_Lx);
    sSite4.y = static_cast<SBYTE> (_ixy % _DC_Lx);
    sSite4.z = static_cast<SBYTE>(_iz_idx >> 1);
    sSite4.w = static_cast<SBYTE>(threadIdx.z + blockIdx.z * blockDim.z);
    UINT uiSiteIndex = _ixy * _DC_GridDimZT + sSite4.z * _DC_Lt + sSite4.w;

    UINT plaqLength = __idx->m_pSmallData[CIndexData::kPlaqLengthIdx];
    UINT plaqCountAll = __idx->m_pSmallData[CIndexData::kPlaqPerSiteIdx] * plaqLength;
    
    //i=0: 12
    //  1: 13
    //  2: 14
    //  3: 23
    //  4: 24
    //  5: 34
    //0->2, 1->4
    BYTE idx = ((_iz_idx & 1) + 1) * 2;

    //Real resThisThread = F(0.0);

    //========================================
    //find plaqutte 1-4, or 2-4
    SIndex first = pCachedPlaqutte[idx * plaqLength + uiSiteIndex * plaqCountAll];
    deviceSU3 toAdd(_deviceGetGaugeBCSU3(pDeviceData, first));
    if (first.NeedToDagger())
    {
        toAdd.Dagger();
    }
    for (BYTE j = 1; j < plaqLength; ++j)
    {
        first = pCachedPlaqutte[idx * plaqLength + j + uiSiteIndex * plaqCountAll];
        deviceSU3 toMul(_deviceGetGaugeBCSU3(pDeviceData, first));
        if (first.NeedToDagger())
        {
            toAdd.MulDagger(toMul);
        }
        else
        {
            toAdd.Mul(toMul);
        }
    }

    //Note that Retr[U14] = Retr[U41], Retr[U24] = Retr[U42], so it is OK
#if !_CLG_DOUBLEFLOAT
    atomicAdd(&results[uiSiteIndex], static_cast<DOUBLE>(betaOverN * (F(3.0) - toAdd.ReTr()) * _deviceGnAcc(sSite4, fGsq)));
#else
    atomicAdd(&results[uiSiteIndex], betaOverN * (F(3.0) - toAdd.ReTr()) * _deviceGnAcc(sSite4, fGsq));
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
    const Real fV413 = _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 0, 2, uiN);

    //===============
    //V423
    const Real fV423 = _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 1, 2, uiN);

    results[uiSiteIndex] = (fV413 + fV423) * betaOverN * fG * sSite4.w;
}

/**
* 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForce4PlaqutteTermSU3_Acc(
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
            deviceSU3 stap(_deviceStapleTerm_Acc_T(byFieldId, pDeviceData, sSite4, fGSq, uiBigIdx, idir, 0));
            stap.Add(_deviceStapleTerm_Acc_T(byFieldId, pDeviceData, sSite4, fGSq, uiBigIdx, idir, 1));
            deviceSU3 force(pDeviceData[linkIndex]);
            force.MulDagger(stap);
            force.Ta();
            force.MulReal(betaOverN);
            pForceData[linkIndex].Add(force);
        }
        else
        {
            UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

            //mu = idir, nu = 4, i = mu
            deviceSU3 stap(_deviceStapleTerm_Acc_XY(pDeviceData, sSite4, fGSq, byFieldId, uiBigIdx, idir, 3));
            deviceSU3 force(pDeviceData[linkIndex]);
            force.MulDagger(stap);
            force.Ta();
            force.MulReal(betaOverN);
            pForceData[linkIndex].Add(force);
        }
    }
}


/**
* Split to 6 functions
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term413_3_Acc(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * fG * F(0.125);

    //===============
    //g t V413
    //add force for rho=3
    UINT uiLink = _deviceGetLinkIndex(uiSiteIndex, 2);
    const deviceSU3 staple_term = _deviceStapleChairTerm1_Acc(byFieldId, pDeviceData, sSite4, uiBigIdx,
        2, 0, 3);
    deviceSU3 force(pDeviceData[uiLink]);
    force.MulDagger(staple_term);
    force.Ta();
    force.MulReal(betaOverN);
    pForceData[uiLink].Add(force);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term413_4_Acc(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
    deviceSU3 *pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * fG * F(0.125);

    //===============
    //g t V413
    //add force for mu=4
    UINT uiLink = _deviceGetLinkIndex(uiSiteIndex, 3);
    const deviceSU3 staple_term = _deviceStapleChairTerm1_Acc(byFieldId, pDeviceData, sSite4, uiBigIdx,
        3, 0, 2);
    deviceSU3 force(pDeviceData[uiLink]);
    force.MulDagger(staple_term);
    force.Ta();
    force.MulReal(betaOverN);
    pForceData[uiLink].Add(force);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term413_1_Acc(
    BYTE byFieldId,
    const deviceSU3 * __restrict__ pDeviceData,
    deviceSU3 *pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * fG * F(0.125);

    //===============
    //+x Omega V413
    //add force for nu=1
    UINT uiLink = _deviceGetLinkIndex(uiSiteIndex, 0);

    const deviceSU3 staple_term = _deviceStapleChairTerm2_Acc(byFieldId, pDeviceData, sSite4, uiBigIdx,
        3, 0, 2);
    deviceSU3 force(pDeviceData[uiLink]);
    force.MulDagger(staple_term);
    force.Ta();
    force.MulReal(betaOverN);
    pForceData[uiLink].Add(force);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term423_3_Acc(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * fG * F(0.125);

    //===============
    //g t V413
    //add force for rho=3
    UINT uiLink = _deviceGetLinkIndex(uiSiteIndex, 2);
    const deviceSU3 staple_term = _deviceStapleChairTerm1_Acc(byFieldId, pDeviceData, sSite4, uiBigIdx,
        2, 1, 3);
    deviceSU3 force(pDeviceData[uiLink]);
    force.MulDagger(staple_term);
    force.Ta();
    force.MulReal(betaOverN);
    pForceData[uiLink].Add(force);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term423_4_Acc(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * fG * F(0.125);

    //===============
    //g t V413
    //add force for mu=4
    UINT uiLink = _deviceGetLinkIndex(uiSiteIndex, 3);
    const deviceSU3 staple_term = _deviceStapleChairTerm1_Acc(byFieldId, pDeviceData, sSite4, uiBigIdx,
        3, 1, 2);
    deviceSU3 force(pDeviceData[uiLink]);
    force.MulDagger(staple_term);
    force.Ta();
    force.MulReal(betaOverN);
    pForceData[uiLink].Add(force);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddForceChairTermSU3_Term423_2_Acc(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    deviceSU3* pForceData,
    Real betaOverN, Real fG)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

    betaOverN = betaOverN * fG * F(0.125);

    //===============
    //+x Omega V413
    //add force for nu=1
    UINT uiLink = _deviceGetLinkIndex(uiSiteIndex, 1);

    const deviceSU3 staple_term = _deviceStapleChairTerm2_Acc(byFieldId, pDeviceData, 
        sSite4, uiBigIdx,
        3, 1, 2);
    deviceSU3 force(pDeviceData[uiLink]);
    force.MulDagger(staple_term);
    force.Ta();
    force.MulReal(betaOverN);
    pForceData[uiLink].Add(force);
}

#pragma endregion


CActionGaugePlaquetteAcceleration::CActionGaugePlaquetteAcceleration()
    : CAction()
    , m_fLastEnergy(F(0.0))
    , m_fNewEnergy(F(0.0))
    , m_fBetaOverN(F(0.1))
    , m_uiPlaqutteCount(0)
{
}

void CActionGaugePlaquetteAcceleration::PrepareForHMC(const CFieldGauge* pGauge, UINT uiUpdateIterate)
{
    if (0 == uiUpdateIterate)
    {
        m_fLastEnergy = Energy(FALSE, pGauge, NULL);
    }
}

void CActionGaugePlaquetteAcceleration::OnFinishTrajectory(UBOOL bAccepted)
{
    if (bAccepted)
    {
        m_fLastEnergy = m_fNewEnergy;
    }
}

void CActionGaugePlaquetteAcceleration::Initial(class CLatticeData* pOwner, const CParameters& param, BYTE byId)
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
    m_uiPlaqutteCount = _HC_Volume * (_HC_Dir - 1) * (_HC_Dir - 2);

    Real fG = 0.1f;
    param.FetchValueReal(_T("AccG"), fG);
    CCommonData::m_fG = fG;

    TArray<INT> centerArray;
    param.FetchValueArrayINT(_T("Center"), centerArray);
    if (centerArray.Num() > 3)
    {
        SSmallInt4 sCenter;
        sCenter.x = static_cast<SBYTE>(centerArray[0]);
        sCenter.y = static_cast<SBYTE>(centerArray[1]);
        sCenter.z = static_cast<SBYTE>(centerArray[2]);
        sCenter.w = static_cast<SBYTE>(centerArray[3]);
        CCommonData::m_sCenter = sCenter;
    }
    else
    {
        CCommonData::m_sCenter.w = 0;
    }
}

void CActionGaugePlaquetteAcceleration::SetBeta(Real fBeta)
{
    CCommonData::m_fBeta = fBeta;
    if (NULL != m_pOwner->m_pGaugeField && EFT_GaugeSU3 == m_pOwner->m_pGaugeField->GetFieldType())
    {
        fBeta = fBeta / F(3.0);
    }
    m_fBetaOverN = fBeta;
}

UBOOL CActionGaugePlaquetteAcceleration::CalculateForceOnGauge(const CFieldGauge * pGauge, class CFieldGauge * pForce, class CFieldGauge * pStaple, ESolverPhase ePhase) const
{
    pGauge->CalculateForceAndStaple(pForce, pStaple, m_fBetaOverN);

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    if (NULL == pGaugeSU3 || NULL == pForceSU3)
    {
        appCrucial(_T("CActionGaugePlaquetteRotating only work with SU3 now.\n"));
        return TRUE;
    }

    preparethread;

    _kernelAddForce4PlaqutteTermSU3_Acc << <block, threads >> >(
        pGaugeSU3->m_pDeviceData, 
        pForceSU3->m_pDeviceData, 
        m_fBetaOverN, 
        CCommonData::m_fG * CCommonData::m_fG,
        pGauge->m_byFieldId);

    _kernelAddForceChairTermSU3_Term413_3_Acc << <block, threads >> >(pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, CCommonData::m_fG);

    _kernelAddForceChairTermSU3_Term413_4_Acc << <block, threads >> >(pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, CCommonData::m_fG);

    _kernelAddForceChairTermSU3_Term413_1_Acc << <block, threads >> >(pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, CCommonData::m_fG);

    _kernelAddForceChairTermSU3_Term423_3_Acc << <block, threads >> >(pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, CCommonData::m_fG);

    _kernelAddForceChairTermSU3_Term423_4_Acc << <block, threads >> >(pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, CCommonData::m_fG);

    _kernelAddForceChairTermSU3_Term423_2_Acc << <block, threads >> >(pGaugeSU3->m_byFieldId, pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData, m_fBetaOverN, CCommonData::m_fG);

    checkCudaErrors(cudaDeviceSynchronize());
    return TRUE;
}

#if !_CLG_DOUBLEFLOAT
DOUBLE CActionGaugePlaquetteAcceleration::Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
#else
Real CActionGaugePlaquetteAcceleration::Energy(UBOOL bBeforeEvolution, const class CFieldGauge* pGauge, const class CFieldGauge* pStable)
#endif
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

    dim3 block2 = block;
    block2.y = block.y * 2;
    _kernelAdd4PlaqutteTermSU3_Acc << <block2, threads >> > (
            pGaugeSU3->m_pDeviceData, 
            appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
            m_fBetaOverN,
            CCommonData::m_fG * CCommonData::m_fG,
            _D_RealThreadBuffer);

    m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    _kernelAddChairTermSU3_Chair_Acc << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        m_fBetaOverN,
        CCommonData::m_fG,
        _D_RealThreadBuffer);

     m_fNewEnergy += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    return m_fNewEnergy;
}

void CActionGaugePlaquetteAcceleration::SetG(Real fG)
{ 
    CCommonData::m_fG = fG;
}

CCString CActionGaugePlaquetteAcceleration::GetInfos(const CCString &tab) const
{
    CCString sRet = tab + _T("Name : CActionGaugePlaquetteAcceleration\n");
    sRet = sRet + tab + _T("Beta : ") + appFloatToString(CCommonData::m_fBeta) + _T("\n");
    sRet = sRet + tab + _T("Omega : ") + appFloatToString(CCommonData::m_fG) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================