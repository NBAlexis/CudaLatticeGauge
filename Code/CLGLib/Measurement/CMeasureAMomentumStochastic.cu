//=============================================================================
// FILENAME : CMeasureAMomentumStochastic.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [07/10/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureAMomentumStochastic)

#pragma region kernels

/**
 * Here JL is using y(T^+ - T^-) - x(T^+ - T^-)
 * Note that we are using fmY = -static_cast<Real>(sSite4.y - sCenter.y) = -y
 * By definition this is in fact -JL
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDotAndGatherXYAMomentumJL(
    const deviceWilsonVectorSU3 * __restrict__ pLeft,
    const deviceWilsonVectorSU3 * __restrict__ pRight,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    BYTE byFieldId, UBOOL bNaive, SSmallInt4 sCenter,
    Real * resultXYPlaneJL)
{
    intokernalInt4;

    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    //SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    Real fmY = -static_cast<Real>(sSite4.y - sCenter.y);
    Real fmX = -static_cast<Real>(sSite4.x - sCenter.x);

    deviceWilsonVectorSU3 jl = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

#pragma region JL

    //idir = mu
    for (UINT idir = 0; idir < 2; ++idir)
    {
        //=========================
        //get things

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

        SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        //deviceSU3 x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir);
        deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
        x_m_mu_Gauge_element.Dagger();

        deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pRight, x_p_mu_Fermion, byFieldId);
        deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pRight, x_m_mu_Fermion, byFieldId);

        //hopping terms

        //U(x,mu) phi(x+ mu)
        deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            if (0 == idir)
            {
                //- k y Omega x_p_m
                //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
                u_phi_x_p_m.MulReal(fmY);
                if (!bNaive)
                {
                    jl.Sub(u_phi_x_p_m);
                }
                jl.Add(gamma4.MulWilsonC(u_phi_x_p_m));
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_phi_x_p_m.MulReal(fmX);
                if (!bNaive)
                {
                    jl.Add(u_phi_x_p_m);
                }
                jl.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
            }
        }
        else
        {
            if (0 == idir)
            {
                //- k y Omega x_p_m
                //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
                u_phi_x_p_m.MulReal(fmY);
                if (!bNaive)
                {
                    jl.Add(u_phi_x_p_m);
                }
                jl.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_phi_x_p_m.MulReal(fmX);
                if (!bNaive)
                {
                    jl.Sub(u_phi_x_p_m);
                }
                jl.Add(gamma4.MulWilsonC(u_phi_x_p_m));
            }
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            if (0 == idir)
            {
                //- k y Omega x_p_m
                //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
                u_dagger_phi_x_m_m.MulReal(fmY);
                if (!bNaive)
                {
                    jl.Sub(u_dagger_phi_x_m_m);
                }
                jl.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_dagger_phi_x_m_m.MulReal(fmX);
                if (!bNaive)
                {
                    jl.Add(u_dagger_phi_x_m_m);
                }
                jl.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
            }
        }
        else
        {
            if (0 == idir)
            {
                //- k y Omega x_p_m
                //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
                u_dagger_phi_x_m_m.MulReal(fmY);
                if (!bNaive)
                {
                    jl.Add(u_dagger_phi_x_m_m);
                }
                jl.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_dagger_phi_x_m_m.MulReal(fmX);
                if (!bNaive)
                {
                    jl.Sub(u_dagger_phi_x_m_m);
                }
                jl.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
            }
        }
    }

#pragma endregion

    CLGComplex cDotRes = pLeft[uiSiteIndex].ConjugateDotC(jl);

    atomicAdd(&resultXYPlaneJL[uiXY], cDotRes.x);
}

//Simplify it because we only use the naive version

__global__ void _CLG_LAUNCH_BOUND
_kernelDotAndGatherXYAMomentumJS(
    const deviceWilsonVectorSU3 * __restrict__ pLeft,
    const deviceWilsonVectorSU3 * __restrict__ pRight,
    BYTE byFieldId, 
    Real * resultXYPlaneJS)
{
    intokernal;

    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;

#pragma region JS

    const deviceWilsonVectorSU3 right_element(pRight[uiSiteIndex]);
    deviceWilsonVectorSU3 js(__chiralGamma[SIGMA12].MulWilsonC(right_element));
    js = __chiralGamma[GAMMA4].MulWilsonC(js);
    js.MulComp(_make_cuComplex(F(0.0), F(-1.0)));

#pragma endregion

    CLGComplex cDotRes = pLeft[uiSiteIndex].ConjugateDotC(js);

    atomicAdd(&resultXYPlaneJS[uiXY], cDotRes.x);
}

/**
 * JS = (gamma4-1) T+  +  (1+gamma4) T-
 * Note that, the T+ term multiply 0.5, it is in fact -0.5, that is -0.5(1-gamma4)T+
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDotAndGatherXYAMomentumJS_Exp(
    const deviceWilsonVectorSU3 * __restrict__ pLeft,
    const deviceWilsonVectorSU3 * __restrict__ pRight,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    BYTE byFieldId,
    Real * resultXYPlaneJS)
{
    intokernalInt4;

    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    //SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    const gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    const gammaMatrix sigma12 = __chiralGamma[SIGMA12];

#pragma region JS

    //idir = mu
    //=========================
    //get things

    //x, mu
    const BYTE idir = 3;
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

    const SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

    const SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    const deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
    x_m_mu_Gauge_element.Dagger();

    const deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pRight, x_p_mu_Fermion, byFieldId);
    const deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pRight, x_m_mu_Fermion, byFieldId);

    //hopping terms

    //U(x,mu) phi(x+ mu)
    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 js = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
    js.Sub(gamma4.MulWilsonC(js));
    //only sin part
    js.MulComp(_make_cuComplex(F(0.0), x_m_mu_Fermion.NeedToOpposite() ? F(0.5) : -F(0.5)));

    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    u_dagger_phi_x_m_m.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    u_dagger_phi_x_m_m.MulComp(_make_cuComplex(F(0.0), x_m_mu_Fermion.NeedToOpposite() ? -F(0.5) : F(0.5)));

    js.Add(u_dagger_phi_x_m_m);
    js = sigma12.MulWilsonC(js);

#pragma endregion

    CLGComplex cDotRes = pLeft[uiSiteIndex].ConjugateDotC(js);

    atomicAdd(&resultXYPlaneJS[uiXY], cDotRes.x);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAMomentumInitialDist(
    UINT* pCount, 
    Real* pValueJL, 
    Real* pValueJS,
    Real* pValueJLPure)
{
    pCount[threadIdx.x] = 0;
    pValueJL[threadIdx.x] = F(0.0);
    pValueJS[threadIdx.x] = F(0.0);
    if (NULL != pValueJLPure)
    {
        pValueJLPure[threadIdx.x] = F(0.0);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDotAndGatherXYAMomentumJPot(
    const deviceWilsonVectorSU3* __restrict__ pLeft,
    const deviceWilsonVectorSU3* __restrict__ pRight,
    const deviceSU3* __restrict__ pAphys,
    SSmallInt4 sCenter,
    BYTE byFieldId,
    Real* resultXYPlaneJPot)
{
    intokernalInt4;

    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    const Real fY = static_cast<Real>(sSite4.y - sCenter.y);
    const Real fX = static_cast<Real>(sSite4.x - sCenter.x);

    //x ay - y ax
    deviceSU3 midY = _deviceGetGaugeBCSU3DirZero(pAphys, uiBigIdx, 1);
    deviceSU3 midX = _deviceGetGaugeBCSU3DirZero(pAphys, uiBigIdx, 0);
    midY.MulReal(fX);
    midX.MulReal(fY);
    midY.Sub(midX);

    CLGComplex cDotRes = pLeft[uiSiteIndex].ConjugateDotC(midY.MulWilsonVector(pRight[uiSiteIndex]));

    atomicAdd(&resultXYPlaneJPot[uiXY], cDotRes.x);
}

#pragma endregion

CMeasureAMomentumStochastic::~CMeasureAMomentumStochastic()
{
    if (NULL != m_pDeviceXYBufferJL)
    {
        checkCudaErrors(cudaFree(m_pDeviceXYBufferJL));
    }

    if (NULL != m_pDeviceXYBufferJS)
    {
        checkCudaErrors(cudaFree(m_pDeviceXYBufferJS));
    }

    if (NULL != m_pDistributionR)
    {
        checkCudaErrors(cudaFree(m_pDistributionR));
    }

    if (NULL != m_pDistributionJL)
    {
        checkCudaErrors(cudaFree(m_pDistributionJL));
    }

    if (NULL != m_pDistributionJS)
    {
        checkCudaErrors(cudaFree(m_pDistributionJS));
    }

    if (NULL != m_pHostDistributionR)
    {
        free(m_pHostDistributionR);
    }

    if (NULL != m_pHostDistributionJL)
    {
        free(m_pHostDistributionJL);
    }

    if (NULL != m_pHostDistributionJS)
    {
        free(m_pHostDistributionJS);
    }

    cudaSafeFree(m_pDeviceXYBufferJLPure);
    cudaSafeFree(m_pDistributionJLPure);
    appSafeFree(m_pHostDistributionJLPure);

    cudaSafeFree(m_pDeviceXYBufferJPot);
    cudaSafeFree(m_pDistributionJPot);
    appSafeFree(m_pHostDistributionJPot);
}

void CMeasureAMomentumStochastic::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasureStochastic::Initial(pOwner, pLatticeData, param, byId);

    INT iValue = 0;
    param.FetchValueINT(_T("MeasurePure"), iValue);
    m_bMeasureJLPure = iValue != 0;

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceXYBufferJL, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceXYBufferJS, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    if (m_bMeasureJLPure)
    {
        checkCudaErrors(cudaMalloc((void**)& m_pDeviceXYBufferJLPure, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
        checkCudaErrors(cudaMalloc((void**)& m_pDeviceXYBufferJPot, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    }

    Reset();

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;

    iValue = 1;
    param.FetchValueINT(_T("Exponential"), iValue);
    m_bExponential = iValue != 0;

    iValue = 1;
    param.FetchValueINT(_T("Naive"), iValue);
    m_bNaive = iValue != 0;

    //assuming the center is really at center
    m_uiMaxR = ((_HC_Lx + 1) / 2 ) * ((_HC_Lx + 1) / 2 )
        + ((_HC_Ly + 1) / 2 ) * ((_HC_Ly + 1) / 2 );

    m_uiEdgeR = ((_HC_Lx + 1) / 2 - 1) * ((_HC_Lx + 1) / 2 - 1);

    checkCudaErrors(cudaMalloc((void**)&m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pDistributionJL, sizeof(Real) * (m_uiMaxR + 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pDistributionJS, sizeof(Real) * (m_uiMaxR + 1)));

    m_pHostDistributionR = (UINT*)malloc(sizeof(UINT) * (m_uiMaxR + 1));
    m_pHostDistributionJL = (Real*)malloc(sizeof(Real) * (m_uiMaxR + 1));
    m_pHostDistributionJS = (Real*)malloc(sizeof(Real) * (m_uiMaxR + 1));

    if (m_bMeasureJLPure)
    {
        checkCudaErrors(cudaMalloc((void**)& m_pDistributionJLPure, sizeof(Real) * (m_uiMaxR + 1)));
        m_pHostDistributionJLPure = (Real*)malloc(sizeof(Real) * (m_uiMaxR + 1));
        checkCudaErrors(cudaMalloc((void**)& m_pDistributionJPot, sizeof(Real) * (m_uiMaxR + 1)));
        m_pHostDistributionJPot = (Real*)malloc(sizeof(Real) * (m_uiMaxR + 1));
    }
}

void CMeasureAMomentumStochastic::OnConfigurationAcceptedZ4(
    const class CFieldGauge* pAcceptGauge, 
    const class CFieldGauge* pCorrespondingStaple, 
    const class CFieldFermion* pZ4, 
    const class CFieldFermion* pInverseZ4,
    UBOOL bStart,
    UBOOL bEnd)
{

    if (bStart)
    {
        _ZeroXYPlane(m_pDeviceXYBufferJL);
        _ZeroXYPlane(m_pDeviceXYBufferJS);
        _ZeroXYPlane(m_pDeviceXYBufferJLPure);
        _ZeroXYPlane(m_pDeviceXYBufferJPot);
    }

    const CFieldFermionWilsonSquareSU3 * pF1W = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(pInverseZ4);
    const CFieldFermionWilsonSquareSU3 * pF2W = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(pZ4);
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pAcceptGauge);

#pragma region Dot

    //The reuslts are stored by AtomicAdd into m_pDeviceXYBufferJL and m_pDeviceXYBufferJS

    preparethread;
    _kernelDotAndGatherXYAMomentumJL << <block, threads >> >(
        pF2W->m_pDeviceData,
        pF1W->m_pDeviceData,
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        m_byFieldId,
        m_bNaive,
        CCommonData::m_sCenter,
        m_pDeviceXYBufferJL
        );

    if (m_bExponential)
    {
        _kernelDotAndGatherXYAMomentumJS_Exp << <block, threads >> >(
            pF2W->m_pDeviceData,
            pF1W->m_pDeviceData,
            pGaugeSU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            m_byFieldId,
            m_pDeviceXYBufferJS
            );
    }
    else
    {
        _kernelDotAndGatherXYAMomentumJS << <block, threads >> >(
            pF2W->m_pDeviceData,
            pF1W->m_pDeviceData,
            m_byFieldId,
            m_pDeviceXYBufferJS
            );
    }

    if (m_bMeasureJLPure)
    {
        CFieldGaugeSU3* pAphys = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pAphys);
        CFieldGaugeSU3* pUpure = dynamic_cast<CFieldGaugeSU3*>(appGetLattice()->m_pUpure);
        if (NULL == pAphys || NULL == pUpure)
        {
            appCrucial(_T("CMeasureAMomentumStochastic: A phys undefined.\n"));
        }
        else
        {

            _kernelDotAndGatherXYAMomentumJL << <block, threads >> > (
                pF2W->m_pDeviceData,
                pF1W->m_pDeviceData,
                pUpure->m_pDeviceData,
                appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
                appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
                m_byFieldId,
                m_bNaive,
                CCommonData::m_sCenter,
                m_pDeviceXYBufferJLPure
                );

            _kernelDotAndGatherXYAMomentumJPot << <block, threads >> > (
                pF2W->m_pDeviceData,
                pF1W->m_pDeviceData,
                pAphys->m_pDeviceData,
                CCommonData::m_sCenter,
                m_byFieldId,
                m_pDeviceXYBufferJPot
                );
        }
    }

#pragma endregion

    if (bEnd)
    {
        const Real fDivider = CCommonData::m_fKai / (m_uiFieldCount * _HC_Lz * _HC_Lt);

        XYDataToRdistri_R(m_pDeviceXYBufferJL, m_pDistributionR, m_pDistributionJL, m_uiMaxR, TRUE, m_byFieldId);
        XYDataToRdistri_R(m_pDeviceXYBufferJS, m_pDistributionR, m_pDistributionJS, m_uiMaxR, false, m_byFieldId);

        checkCudaErrors(cudaMemcpy(m_pHostDistributionR, m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(m_pHostDistributionJL, m_pDistributionJL, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));        
        checkCudaErrors(cudaMemcpy(m_pHostDistributionJS, m_pDistributionJS, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));

        if (m_bMeasureJLPure)
        {
            XYDataToRdistri_R(m_pDeviceXYBufferJLPure, m_pDistributionR, m_pDistributionJLPure, m_uiMaxR, false, m_byFieldId);
            XYDataToRdistri_R(m_pDeviceXYBufferJPot, m_pDistributionR, m_pDistributionJPot, m_uiMaxR, false, m_byFieldId);
            checkCudaErrors(cudaMemcpy(m_pHostDistributionJLPure, m_pDistributionJLPure, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
            checkCudaErrors(cudaMemcpy(m_pHostDistributionJPot, m_pDistributionJPot, sizeof(Real)* (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
        }

        Real fAverageJLInner = F(0.0);
        Real fAverageJLAll = F(0.0);
        Real fAverageJSInner = F(0.0);
        Real fAverageJSAll = F(0.0);
        Real fAverageJLPureInner = F(0.0);
        Real fAverageJLPureAll = F(0.0);
        Real fAverageJPotInner = F(0.0);
        Real fAverageJPotAll = F(0.0);
        UINT uiInnerPointsJLInner = 0;
        UINT uiInnerPointsJLAll = 0;
        if (0 == m_uiConfigurationCount)
        {
            assert(0 == m_lstR.Num());
            assert(0 == m_lstJL.Num());
            assert(0 == m_lstJS.Num());
            assert(0 == m_lstJLPure.Num());
            assert(0 == m_lstJPot.Num());

            for (UINT uiL = 0; uiL <= m_uiMaxR; ++uiL)
            {
                if (m_pHostDistributionR[uiL] > 0)
                {
                    m_lstR.AddItem(uiL);
                    m_lstJL.AddItem(m_pHostDistributionJL[uiL] * fDivider);
                    m_lstJS.AddItem(m_pHostDistributionJS[uiL] * fDivider);
                    if (m_bMeasureJLPure)
                    {
                        m_lstJLPure.AddItem(m_pHostDistributionJLPure[uiL] * fDivider);
                        m_lstJPot.AddItem(m_pHostDistributionJPot[uiL] * fDivider);
                    }

                    uiInnerPointsJLAll += m_pHostDistributionR[uiL];
                    fAverageJLAll += m_pHostDistributionR[uiL] * m_pHostDistributionJL[uiL];
                    fAverageJSAll += m_pHostDistributionR[uiL] * m_pHostDistributionJS[uiL];
                    if (m_bMeasureJLPure)
                    {
                        fAverageJLPureAll += m_pHostDistributionR[uiL] * m_pHostDistributionJLPure[uiL];
                        fAverageJPotAll += m_pHostDistributionR[uiL] * m_pHostDistributionJPot[uiL];
                    }
                    if (uiL < m_uiEdgeR)
                    {
                        uiInnerPointsJLInner += m_pHostDistributionR[uiL];
                        fAverageJLInner += m_pHostDistributionR[uiL] * m_pHostDistributionJL[uiL];
                        fAverageJSInner += m_pHostDistributionR[uiL] * m_pHostDistributionJS[uiL];
                        if (m_bMeasureJLPure)
                        {
                            fAverageJLPureInner += m_pHostDistributionR[uiL] * m_pHostDistributionJLPure[uiL];
                            fAverageJPotInner += m_pHostDistributionR[uiL] * m_pHostDistributionJPot[uiL];
                        }
                    }

                    if (m_bShowResult)
                    {
                        appDetailed(_T("(%f)JL=%f, JS=%f, JLPure=%f, JPot=%f\n"),
                            _hostsqrt(static_cast<Real>(uiL)),
                            m_pHostDistributionJL[uiL],
                            m_pHostDistributionJS[uiL],
                            m_bMeasureJLPure ? m_pHostDistributionJLPure[uiL] : 0.0f,
                            m_bMeasureJLPure ? m_pHostDistributionJPot[uiL] : 0.0f);

                    }
                }
            }
        }
        else
        {
            for (INT i = 0; i < m_lstR.Num(); ++i)
            {
                assert(m_pHostDistributionR[m_lstR[i]] > 0);
                m_lstJL.AddItem(m_pHostDistributionJL[m_lstR[i]] * fDivider);
                m_lstJS.AddItem(m_pHostDistributionJS[m_lstR[i]] * fDivider);
                if (m_bMeasureJLPure)
                {
                    m_lstJLPure.AddItem(m_pHostDistributionJLPure[m_lstR[i]] * fDivider);
                    m_lstJPot.AddItem(m_pHostDistributionJPot[m_lstR[i]] * fDivider);
                }

                uiInnerPointsJLAll += m_pHostDistributionR[m_lstR[i]];
                fAverageJLAll += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJL[m_lstR[i]];
                fAverageJSAll += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJS[m_lstR[i]];
                if (m_bMeasureJLPure)
                {
                    fAverageJLPureAll += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJLPure[m_lstR[i]];
                    fAverageJPotAll += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJPot[m_lstR[i]];
                }
                if (m_lstR[i] < m_uiEdgeR)
                {
                    uiInnerPointsJLInner += m_pHostDistributionR[m_lstR[i]];
                    fAverageJLInner += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJL[m_lstR[i]];
                    fAverageJSInner += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJS[m_lstR[i]];
                    if (m_bMeasureJLPure)
                    {
                        fAverageJLPureInner += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJLPure[m_lstR[i]];
                        fAverageJPotInner += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJPot[m_lstR[i]];
                    }
                }

                if (m_bShowResult)
                {
                    appDetailed(_T("(%f)JL=%f, JS=%f, JLPure=%f, JPot=%f\n"),
                        _hostsqrt(static_cast<Real>(m_lstR[i])),
                        m_pHostDistributionJL[m_lstR[i]],
                        m_pHostDistributionJS[m_lstR[i]],
                        m_bMeasureJLPure ? m_pHostDistributionJLPure[m_lstR[i]] : 0.0f,
                        m_bMeasureJLPure ? m_pHostDistributionJPot[m_lstR[i]] : 0.0f);
                }
            }
        }

        if (uiInnerPointsJLAll > 0)
        {
            fAverageJLAll = fAverageJLAll / uiInnerPointsJLAll;
            fAverageJSAll = fAverageJSAll / uiInnerPointsJLAll;
            fAverageJLPureAll = fAverageJLPureAll / uiInnerPointsJLAll;
            fAverageJPotAll = fAverageJPotAll / uiInnerPointsJLAll;
        }
        if (uiInnerPointsJLInner > 0)
        {
            fAverageJLInner = fAverageJLInner / uiInnerPointsJLInner;
            fAverageJSInner = fAverageJSInner / uiInnerPointsJLAll;
            fAverageJLPureInner = fAverageJLPureInner / uiInnerPointsJLAll;
            fAverageJPotInner = fAverageJPotInner / uiInnerPointsJLAll;
        }
        m_lstJLAll.AddItem(fAverageJLAll * fDivider);
        m_lstJLInner.AddItem(fAverageJLInner * fDivider);
        m_lstJSAll.AddItem(fAverageJSAll* fDivider);
        m_lstJSInner.AddItem(fAverageJSInner* fDivider);
        if (m_bMeasureJLPure)
        {
            m_lstJLPureAll.AddItem(fAverageJLPureAll* fDivider);
            m_lstJLPureInner.AddItem(fAverageJLPureInner* fDivider);
            m_lstJPotAll.AddItem(fAverageJPotAll* fDivider);
            m_lstJPotInner.AddItem(fAverageJPotInner* fDivider);
        }
        ++m_uiConfigurationCount;
    }
}

void CMeasureAMomentumStochastic::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{
    //do nothing, we use OnConfigurationAcceptedZ4 instead
}

void CMeasureAMomentumStochastic::Average(UINT )
{
    //nothing to do
}

void CMeasureAMomentumStochastic::Report()
{
    assert(static_cast<INT>(m_uiConfigurationCount) * m_lstR.Num() == m_lstJL.Num());
    assert(static_cast<INT>(m_uiConfigurationCount) * m_lstR.Num() == m_lstJS.Num());

    appSetLogDate(FALSE);

    appGeneral(_T("jgr={"));

    for (INT i = 0; i < m_lstR.Num(); ++i)
    {
        appGeneral(_T("%2.12f%s"), _hostsqrt(static_cast<Real>(m_lstR[i])), (i == m_lstR.Num() - 1) ? _T("") : _T(", "));
    }

    appGeneral(_T("};\n"));

    appGeneral(_T("================ JFL ===================\n"));

    ReportDistributeWithR_R(m_uiConfigurationCount, m_lstR.Num(), m_lstJL);

    appGeneral(_T("================ JFS ===================\n"));

    ReportDistributeWithR_R(m_uiConfigurationCount, m_lstR.Num(), m_lstJS);

    appGeneral(_T("================ JFL pure ===================\n"));

    ReportDistributeWithR_R(m_uiConfigurationCount, m_lstR.Num(), m_lstJLPure);

    appGeneral(_T("================ Jpot ===================\n"));

    ReportDistributeWithR_R(m_uiConfigurationCount, m_lstR.Num(), m_lstJPot);

    appGeneral(_T("==========================================================================\n"));
    appSetLogDate(TRUE);
}

void CMeasureAMomentumStochastic::Reset()
{
    m_uiConfigurationCount = 0;

    m_lstR.RemoveAll();
    m_lstJL.RemoveAll();
    m_lstJLPure.RemoveAll();
    m_lstJPot.RemoveAll();
    m_lstJS.RemoveAll();

    m_lstJLAll.RemoveAll();
    m_lstJLInner.RemoveAll();
    m_lstJSAll.RemoveAll();
    m_lstJSInner.RemoveAll();
    m_lstJLPureAll.RemoveAll();
    m_lstJLPureInner.RemoveAll();
    m_lstJPotAll.RemoveAll();
    m_lstJPotInner.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================