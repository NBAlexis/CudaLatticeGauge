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
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
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

    deviceWilsonVectorSU3 right_element(pRight[uiSiteIndex]);
    deviceWilsonVectorSU3 js(__chiralGamma[SIGMA12].MulWilsonC(right_element));
    js = __chiralGamma[GAMMA4].MulWilsonC(js);
    js.MulComp(_make_cuComplex(F(0.0), F(-1.0)));

#pragma endregion

    CLGComplex cDotRes = pLeft[uiSiteIndex].ConjugateDotC(js);

    atomicAdd(&resultXYPlaneJS[uiXY], cDotRes.x);
}

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
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix sigma12 = __chiralGamma[SIGMA12];

#pragma region JS

    //idir = mu
    //=========================
    //get things

    //x, mu
    BYTE idir = 3;
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

__global__ void
_CLG_LAUNCH_BOUND
_kernelAMomentumInitialDist(UINT* pCount, Real* pValueJL, Real* pValueJS)
{
    pCount[threadIdx.x] = 0;
    pValueJL[threadIdx.x] = F(0.0);
    pValueJS[threadIdx.x] = F(0.0);
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelAMomentumMeasureDist(
    const Real* __restrict__ jlXY,
    const Real* __restrict__ jsXY,
    SSmallInt4 sCenter, UINT uiMax, BYTE byFieldId,
    UINT* counter, Real* jl, Real* js)
{
    UINT uiXY = (threadIdx.x + blockIdx.x * blockDim.x);
    SBYTE uiX = static_cast<SBYTE>(uiXY / _DC_Ly);
    SBYTE uiY = static_cast<SBYTE>(uiXY % _DC_Ly);
    UINT uiC = (sCenter.x - uiX) * (sCenter.x - uiX)
        + (sCenter.y - uiY) * (sCenter.y - uiY);

    SSmallInt4 sSite4;
    sSite4.z = sCenter.z;
    sSite4.w = sCenter.w;
    sSite4.x = uiX;
    sSite4.y = uiY;
    if (uiC <= uiMax && !__idx->_deviceGetMappingIndex(sSite4, byFieldId).IsDirichlet())
    {
        atomicAdd(&counter[uiC], 1);
        atomicAdd(&jl[uiC], jlXY[uiXY]);
        atomicAdd(&js[uiC], jsXY[uiXY]);
    }
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelAMomentumAverageDist(UINT* pCount, Real* pValueJL, Real* pValueJS)
{
    UINT uiIdx = threadIdx.x;
    if (pCount[uiIdx] > 0)
    {
        pValueJL[uiIdx] = pValueJL[uiIdx] / static_cast<Real>(pCount[uiIdx]);
        pValueJS[uiIdx] = pValueJS[uiIdx] / static_cast<Real>(pCount[uiIdx]);
    }
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
}

void CMeasureAMomentumStochastic::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasureStochastic::Initial(pOwner, pLatticeData, param, byId);

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceXYBufferJL, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceXYBufferJS, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));

    Reset();

    INT iValue = 1;
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

#pragma endregion

    if (bEnd)
    {
        Real fDivider = CCommonData::m_fKai / (m_uiFieldCount * _HC_Lz * _HC_Lt);
        dim3 block2(_HC_DecompX, 1, 1);
        dim3 threads2(_HC_DecompLx, 1, 1);
        dim3 block3(m_uiMaxR + 1, 1, 1);
        dim3 threads3(m_uiMaxR + 1, 1, 1);

        _kernelAMomentumInitialDist << <block3, threads3 >> >(m_pDistributionR, m_pDistributionJL, m_pDistributionJS);

        _kernelAMomentumMeasureDist << <block2, threads2 >> >(
            m_pDeviceXYBufferJL,
            m_pDeviceXYBufferJS,
            CCommonData::m_sCenter,
            m_uiMaxR,
            m_byFieldId,
            m_pDistributionR,
            m_pDistributionJL,
            m_pDistributionJS
            );

        _kernelAMomentumAverageDist << <block3, threads3 >> >(m_pDistributionR, m_pDistributionJL, m_pDistributionJS);

        //extract res
        checkCudaErrors(cudaMemcpy(m_pHostDistributionR, m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(m_pHostDistributionJL, m_pDistributionJL, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(m_pHostDistributionJS, m_pDistributionJS, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));

        Real fAverageJLInner = F(0.0);
        UINT uiInnerPointsJLInner = 0;
        Real fAverageJLAll = F(0.0);
        UINT uiInnerPointsJLAll = 0;
        if (0 == m_uiConfigurationCount)
        {
            assert(0 == m_lstR.Num());
            assert(0 == m_lstJL.Num());
            assert(0 == m_lstJS.Num());

            for (UINT uiL = 0; uiL <= m_uiMaxR; ++uiL)
            {
                if (m_pHostDistributionR[uiL] > 0)
                {
                    m_lstR.AddItem(uiL);
                    m_lstJL.AddItem(m_pHostDistributionJL[uiL] * fDivider);
                    m_lstJS.AddItem(m_pHostDistributionJS[uiL] * fDivider);

                    uiInnerPointsJLAll += m_pHostDistributionR[uiL];
                    fAverageJLAll += m_pHostDistributionR[uiL] * m_pHostDistributionJL[uiL];
                    if (uiL < m_uiEdgeR)
                    {
                        uiInnerPointsJLInner += m_pHostDistributionR[uiL];
                        fAverageJLInner += m_pHostDistributionR[uiL] * m_pHostDistributionJL[uiL];
                    }

                    if (m_bShowResult)
                    {
                        appDetailed(_T("(%f)JL=%f, JS=%f\n"),
                            _hostsqrt(static_cast<Real>(uiL)),
                            m_pHostDistributionJL[uiL],
                            m_pHostDistributionJS[uiL]);

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

                uiInnerPointsJLAll += m_pHostDistributionR[m_lstR[i]];
                fAverageJLAll += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJL[m_lstR[i]];
                if (m_lstR[i] < m_uiEdgeR)
                {
                    uiInnerPointsJLInner += m_pHostDistributionR[m_lstR[i]];
                    fAverageJLInner += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJL[m_lstR[i]];
                }

                if (m_bShowResult)
                {
                    appDetailed(_T("(%f)JL=%f, JS=%f\n"),
                        _hostsqrt(static_cast<Real>(m_lstR[i])),
                        m_pHostDistributionJL[m_lstR[i]],
                        m_pHostDistributionJS[m_lstR[i]]);
                }
            }
        }

        if (uiInnerPointsJLAll > 0)
        {
            fAverageJLAll = fAverageJLAll / uiInnerPointsJLAll;
        }
        if (uiInnerPointsJLInner > 0)
        {
            fAverageJLInner = fAverageJLInner / uiInnerPointsJLInner;
        }
        m_lstJLAll.AddItem(fAverageJLAll * fDivider);
        m_lstJLInner.AddItem(fAverageJLInner * fDivider);

        ++m_uiConfigurationCount;
    }
}

void CMeasureAMomentumStochastic::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{
    CFieldFermionWilsonSquareSU3* pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    CFieldFermionWilsonSquareSU3* pF2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    _ZeroXYPlane(m_pDeviceXYBufferJL);
    _ZeroXYPlane(m_pDeviceXYBufferJS);

    for (UINT i = 0; i < m_uiFieldCount; ++i)
    {
        pF1->InitialField(EFIT_RandomZ4);
        pF1->FixBoundary();
        //pF1->DebugPrintMe();
        //pGauge->DebugPrintMe();
        pF1->CopyTo(pF2);
        pF1->InverseD(pGauge);       

#pragma region Dot

        preparethread;
        _kernelDotAndGatherXYAMomentumJL<<<block, threads>>>(
            pF2->m_pDeviceData,
            pF1->m_pDeviceData,
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
                pF2->m_pDeviceData,
                pF1->m_pDeviceData,
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
                pF2->m_pDeviceData,
                pF1->m_pDeviceData,
                m_byFieldId,
                m_pDeviceXYBufferJS
            );
        }

#pragma endregion

    }

    Real fDivider = CCommonData::m_fKai / (m_uiFieldCount * _HC_Lz * _HC_Lt);
    dim3 block2(_HC_DecompX, 1, 1);
    dim3 threads2(_HC_DecompLx, 1, 1);
    dim3 block3(m_uiMaxR + 1, 1, 1);
    dim3 threads3(m_uiMaxR + 1, 1, 1);

    _kernelAMomentumInitialDist << <block3, threads3 >> >(m_pDistributionR, m_pDistributionJL, m_pDistributionJS);

    _kernelAMomentumMeasureDist << <block2, threads2 >> >(
        m_pDeviceXYBufferJL,
        m_pDeviceXYBufferJS,
        CCommonData::m_sCenter,
        m_uiMaxR,
        m_byFieldId,
        m_pDistributionR,
        m_pDistributionJL,
        m_pDistributionJS
        );

    _kernelAMomentumAverageDist << <block3, threads3 >> >(m_pDistributionR, m_pDistributionJL, m_pDistributionJS);

    //extract res
    checkCudaErrors(cudaMemcpy(m_pHostDistributionR, m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostDistributionJL, m_pDistributionJL, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostDistributionJS, m_pDistributionJS, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));

    Real fAverageJLInner = F(0.0);
    UINT uiInnerPointsJLInner = 0;
    Real fAverageJLAll = F(0.0);
    UINT uiInnerPointsJLAll = 0;
    if (0 == m_uiConfigurationCount)
    {
        assert(0 == m_lstR.Num());
        assert(0 == m_lstJL.Num());
        assert(0 == m_lstJS.Num());

        for (UINT uiL = 0; uiL <= m_uiMaxR; ++uiL)
        {
            if (m_pHostDistributionR[uiL] > 0)
            {
                m_lstR.AddItem(uiL);
                m_lstJL.AddItem(m_pHostDistributionJL[uiL] * fDivider);
                m_lstJS.AddItem(m_pHostDistributionJS[uiL] * fDivider);

                uiInnerPointsJLAll += m_pHostDistributionR[uiL];
                fAverageJLAll += m_pHostDistributionR[uiL] * m_pHostDistributionJL[uiL];
                if (uiL < m_uiEdgeR)
                {
                    uiInnerPointsJLInner += m_pHostDistributionR[uiL];
                    fAverageJLInner += m_pHostDistributionR[uiL] * m_pHostDistributionJL[uiL];
                }

                if (m_bShowResult)
                {
                    appDetailed(_T("(%f)JL=%f, JS=%f\n"),
                        _hostsqrt(static_cast<Real>(uiL)),
                        m_pHostDistributionJL[uiL],
                        m_pHostDistributionJS[uiL]);

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

            uiInnerPointsJLAll += m_pHostDistributionR[m_lstR[i]];
            fAverageJLAll += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJL[m_lstR[i]];
            if (m_lstR[i] < m_uiEdgeR)
            {
                uiInnerPointsJLInner += m_pHostDistributionR[m_lstR[i]];
                fAverageJLInner += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJL[m_lstR[i]];
            }

            if (m_bShowResult)
            {
                appDetailed(_T("(%f)JL=%f, JS=%f\n"),
                    _hostsqrt(static_cast<Real>(m_lstR[i])),
                    m_pHostDistributionJL[m_lstR[i]],
                    m_pHostDistributionJS[m_lstR[i]]);
            }
        }
    }

    if (uiInnerPointsJLAll > 0)
    {
        fAverageJLAll = fAverageJLAll / uiInnerPointsJLAll;
    }
    if (uiInnerPointsJLInner > 0)
    {
        fAverageJLInner = fAverageJLInner / uiInnerPointsJLInner;
    }
    m_lstJLAll.AddItem(fAverageJLAll * fDivider);
    m_lstJLInner.AddItem(fAverageJLInner * fDivider);

    pF1->Return();
    pF2->Return();
    ++m_uiConfigurationCount;
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

    appGeneral(_T("};\nJL={\n"));

    for (UINT j = 0; j < m_uiConfigurationCount; ++j)
    {
        appGeneral(_T("{"));
        for (INT i = 0; i < m_lstR.Num(); ++i)
        {
            appGeneral(_T("%2.12f%s"), m_lstJL[j * m_lstR.Num() + i], (i == m_lstR.Num() - 1) ? _T("") : _T(", "));
        }
        appGeneral(_T("}%s\n"), (j == m_uiConfigurationCount - 1) ? _T("") : _T(","));
    }

    appGeneral(_T("\n};\nJS={\n"));

    for (UINT j = 0; j < m_uiConfigurationCount; ++j)
    {
        appGeneral(_T("{"));
        for (INT i = 0; i < m_lstR.Num(); ++i)
        {
            appGeneral(_T("%2.12f%s"), m_lstJS[j * m_lstR.Num() + i], (i == m_lstR.Num() - 1) ? _T("") : _T(", "));
        }
        appGeneral(_T("}%s\n"), (j == m_uiConfigurationCount - 1) ? _T("") : _T(","));
    }

    appGeneral(_T("\n};"));

    appGeneral(_T("==========================================================================\n"));
    appSetLogDate(TRUE);
}

void CMeasureAMomentumStochastic::Reset()
{
    m_uiConfigurationCount = 0;

    m_lstR.RemoveAll();
    m_lstJL.RemoveAll();
    m_lstJS.RemoveAll();

    m_lstJLAll.RemoveAll();
    m_lstJLInner.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================