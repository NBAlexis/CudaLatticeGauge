//=============================================================================
// FILENAME : CMeasureAMomemtumJG.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/21/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureAMomentumJF)

#pragma region kernles

/**
* We only need to calculate one site
* res is a deviceWilsonVectorSU3 x L array
* res[c,s] is Dy on source[c,s]
*/
__global__ void
_CLG_LAUNCH_(12, 1)
_kernel_XDy_yDx(
    const deviceSU3* __restrict__ pGauge,
    deviceWilsonVectorSU3** pSources,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    SSmallInt4 sSite4,
    SSmallInt4 sCenter,
    BYTE byArrayIdx,
    BYTE byFieldId,
    UBOOL bNaive,
    deviceWilsonVectorSU3* res)
{
    //s * 3 + c
    UINT uiC = threadIdx.x;
    UINT uiS = threadIdx.y;
    UINT uiCS = uiS * 3 + uiC;

    UINT uiSiteIndex = _deviceGetSiteIndex(sSite4);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        res[byArrayIdx] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    Real fmY = -static_cast<Real>(sSite4.y - sCenter.y);
    Real fmX = -static_cast<Real>(sSite4.x - sCenter.x);

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

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

        deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pSources[uiCS], x_p_mu_Fermion, byFieldId);
        deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pSources[uiCS], x_m_mu_Fermion, byFieldId);

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
                    result.Sub(u_phi_x_p_m);
                }
                result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_phi_x_p_m.MulReal(fmX);
                if (!bNaive)
                {
                    result.Add(u_phi_x_p_m);
                }
                result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
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
                    result.Add(u_phi_x_p_m);
                }
                result.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_phi_x_p_m.MulReal(fmX);
                if (!bNaive)
                {
                    result.Sub(u_phi_x_p_m);
                }
                result.Add(gamma4.MulWilsonC(u_phi_x_p_m));
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
                    result.Sub(u_dagger_phi_x_m_m);
                }
                result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_dagger_phi_x_m_m.MulReal(fmX);
                if (!bNaive)
                {
                    result.Add(u_dagger_phi_x_m_m);
                }
                result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
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
                    result.Add(u_dagger_phi_x_m_m);
                }
                result.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_dagger_phi_x_m_m.MulReal(fmX);
                if (!bNaive)
                {
                    result.Sub(u_dagger_phi_x_m_m);
                }
                result.Sub(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
            }
        }
    }

    //Note that, it is not res[byArrayIdx] = result
    //It is res[c,s] = delta_{cc}delta_ss result[c,s]
    res[byArrayIdx].m_d[uiS].m_ve[uiC] = result.m_d[uiS].m_ve[uiC];
}

__global__ void
_CLG_LAUNCH_(12, 1)
_kernel_FS(
    deviceWilsonVectorSU3** pSources,
    SSmallInt4 sSite4,
    BYTE byArrayIdx,
    BYTE byFieldId,
    deviceWilsonVectorSU3* res)
{
    //s * 3 + c
    UINT uiC = threadIdx.x;
    UINT uiS = threadIdx.y;
    UINT uiCS = uiS * 3 + uiC;

    UINT uiSiteIndex = _deviceGetSiteIndex(sSite4);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        res[byArrayIdx] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    deviceWilsonVectorSU3 right_element(pSources[uiCS][uiSiteIndex]);
    deviceWilsonVectorSU3 term3(__chiralGamma[SIGMA12].MulWilsonC(right_element));
    term3 = __chiralGamma[GAMMA4].MulWilsonC(term3);
    term3.MulComp(_make_cuComplex(F(0.0), F(-1.0)));

    //Note that, it is not res[byArrayIdx] = result
    //It is res[c,s] = delta_{cc}delta_ss result[c,s]
    res[byArrayIdx].m_d[uiS].m_ve[uiC] = term3.m_d[uiS].m_ve[uiC];
}

__global__ void
_CLG_LAUNCH_(12, 1)
_kernel_FS_Exponential(
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3** pSources,
    SSmallInt4 sSite4,
    BYTE byArrayIdx,
    BYTE byFieldId,
    deviceWilsonVectorSU3* res)
{
    //s * 3 + c
    UINT uiC = threadIdx.x;
    UINT uiS = threadIdx.y;
    UINT uiCS = uiS * 3 + uiC;

    UINT uiSiteIndex = _deviceGetSiteIndex(sSite4);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        res[byArrayIdx] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    gammaMatrix gamma4 = __chiralGamma[GAMMA4];
    gammaMatrix sigma12 = __chiralGamma[SIGMA12];

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

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

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pSources[uiCS], x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pSources[uiCS], x_m_mu_Fermion, byFieldId);

    //hopping terms

    //U(x,mu) phi(x+ mu)
    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
    u_phi_x_p_m.Sub(gamma4.MulWilsonC(u_phi_x_p_m));
    //only sin part
    u_phi_x_p_m.MulComp(_make_cuComplex(F(0.0), x_m_mu_Fermion.NeedToOpposite() ? F(0.5) : - F(0.5)));

    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    u_dagger_phi_x_m_m.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    u_dagger_phi_x_m_m.MulComp(_make_cuComplex(F(0.0), x_m_mu_Fermion.NeedToOpposite() ? -F(0.5) : F(0.5)));
    
    u_phi_x_p_m.Add(u_dagger_phi_x_m_m);
    u_phi_x_p_m = sigma12.MulWilsonC(u_phi_x_p_m);

    //Note that, it is not res[byArrayIdx] = result
    //It is res[c,s] = delta_{cc}delta_ss result[c,s]
    res[byArrayIdx].m_d[uiS].m_ve[uiC] = u_phi_x_p_m.m_d[uiS].m_ve[uiC];
}

/**
*
*/
__global__ void
_CLG_LAUNCH_(128, 1)
_kernel_Trace_JFLS(
    const deviceWilsonVectorSU3* __restrict__ pOperator,
    CLGComplex* pResLine,
    Real fKappa)
{
    UINT uiIdx = threadIdx.x;
    pResLine[uiIdx] = cuCmulf_cr(pOperator[uiIdx].Sum(), fKappa);
}

#pragma endregion

CMeasureAMomentumJF::~CMeasureAMomentumJF()
{
    if (NULL != m_pHostDataBuffer)
    {
        free(m_pHostDataBuffer);
    }
    if (NULL != m_pDeviceDataBufferS)
    {
        checkCudaErrors(cudaFree(m_pDeviceDataBufferS));
    }
    if (NULL != m_pDeviceDataBufferL)
    {
        checkCudaErrors(cudaFree(m_pDeviceDataBufferL));
    }
    if (NULL != m_pOperatorDataS)
    {
        checkCudaErrors(cudaFree(m_pOperatorDataS));
    }
    if (NULL != m_pOperatorDataL)
    {
        checkCudaErrors(cudaFree(m_pOperatorDataL));
    }
}

void CMeasureAMomentumJF::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    m_pHostDataBuffer = (CLGComplex*)malloc(sizeof(CLGComplex) * (_HC_Lx - 1));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDataBufferS, sizeof(CLGComplex) * (_HC_Lx - 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDataBufferL, sizeof(CLGComplex) * (_HC_Lx - 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pOperatorDataS, sizeof(deviceWilsonVectorSU3) * (_HC_Lx - 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pOperatorDataL, sizeof(deviceWilsonVectorSU3) * (_HC_Lx - 1)));

    Reset();

    INT iValue = 1;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;

    iValue = 1;
    param.FetchValueINT(_T("Naive"), iValue);
    m_bNaive = iValue != 0;

    iValue = 1;
    param.FetchValueINT(_T("Exponential"), iValue);
    m_bExponential = iValue != 0;
}

void CMeasureAMomentumJF::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{
    //if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    //{
    //    appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
    //    return;
    //}
    //const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    //dim3 _blocks(1, 1, 1);
    //dim3 _thread1(12, 1, 1);
    //for (UINT i = 1; i < _HC_Lx; ++i)
    //{
    //    CFieldFermionWilsonSquareSU3DR* pFermionSources[12];
    //    for (UINT j = 0; j < 12; ++j)
    //    {
    //        pFermionSources[j] = dynamic_cast<CFieldFermionWilsonSquareSU3DR*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    //        if (NULL == pFermionSources[j])
    //        {
    //            appCrucial(_T("Meson correlator only implemented with Wilson SU3 Dirichlet Rotating\n"));
    //            _FAIL_EXIT;
    //        }
    //    }

    //    deviceWilsonVectorSU3* pDevicePtr[12];
    //    SSmallInt4 sourceSite;
    //    sourceSite.x = static_cast<SBYTE>(i);
    //    sourceSite.y = CCommonData::m_sCenter.y;
    //    sourceSite.z = CCommonData::m_sCenter.z;
    //    sourceSite.w = CCommonData::m_sCenter.w;
    //    for (BYTE s = 0; s < 4; ++s)
    //    {
    //        for (BYTE c = 0; c < 3; ++c)
    //        {
    //            SFermionSource sourceData;
    //            sourceData.m_eSourceType = EFS_Point;
    //            sourceData.m_sSourcePoint = sourceSite;
    //            sourceData.m_byColorIndex = c;
    //            sourceData.m_bySpinIndex = s;

    //            pFermionSources[s * 3 + c]->InitialAsSource(sourceData);

    //            if (NULL != appGetFermionSolver() && !appGetFermionSolver()->IsAbsoluteAccuracy())
    //            {
    //                pFermionSources[s * 3 + c]->m_fLength = pFermionSources[s * 3 + c]->Dot(pFermionSources[s * 3 + c]).x;
    //            }
    //            pFermionSources[s * 3 + c]->InverseD(pGaugeSU3);
    //            pDevicePtr[s * 3 + c] = pFermionSources[s * 3 + c]->m_pDeviceData;
    //        }
    //    }

    //    deviceWilsonVectorSU3** ppDevicePtr;
    //    checkCudaErrors(cudaMalloc((void**)&ppDevicePtr, sizeof(deviceWilsonVectorSU3*) * 12));
    //    checkCudaErrors(cudaMemcpy(ppDevicePtr, pDevicePtr, sizeof(deviceWilsonVectorSU3*) * 12, cudaMemcpyHostToDevice));

    //    _kernel_XDy_yDx << <_blocks, _thread1 >> > (
    //            pGaugeSU3->m_pDeviceData,
    //            ppDevicePtr,
    //            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
    //            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
    //            sourceSite,
    //            CCommonData::m_sCenter,
    //            static_cast<BYTE>(i - 1),
    //            m_byFieldId,
    //            m_bNaive,
    //            m_pOperatorDataL
    //            );

    //    if (m_bExponential)
    //    {
    //        _kernel_FS_Exponential << <_blocks, _thread1 >> > (
    //            pGaugeSU3->m_pDeviceData,
    //            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
    //            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
    //            ppDevicePtr,
    //            sourceSite,
    //            static_cast<BYTE>(i - 1),
    //            m_byFieldId,
    //            m_pOperatorDataS
    //            );
    //    }
    //    else
    //    {
    //        _kernel_FS << <_blocks, _thread1 >> > (
    //            ppDevicePtr,
    //            sourceSite,
    //            static_cast<BYTE>(i - 1),
    //            m_byFieldId,
    //            m_pOperatorDataS
    //            );
    //    }

    //    checkCudaErrors(cudaFree(ppDevicePtr));
    //    for (UINT j = 0; j < 12; ++j)
    //    {
    //        pFermionSources[j]->Return();
    //    }
    //}

    //++m_uiConfigurationCount;
    //dim3 _thread2(_HC_Lx - 1, 1, 1);
    //Real fKappa = CCommonData::m_fKai;
    //_kernel_Trace_JFLS << <_blocks, _thread2 >> > (m_pOperatorDataL, m_pDeviceDataBufferL, fKappa);
    //_kernel_Trace_JFLS << <_blocks, _thread2 >> > (m_pOperatorDataS, m_pDeviceDataBufferS, fKappa);

    //if (m_bShowResult)
    //{
    //    appDetailed(_T("\n\n ==================== Angular Momentum (%d con)============================ \n\n"), m_uiConfigurationCount);
    //    appDetailed(_T(" ----------- Orbital ------------- \n"));
    //}
    //checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBufferL, sizeof(CLGComplex) * (_HC_Lx - 1), cudaMemcpyDeviceToHost));

    //for (UINT i = 0; i < _HC_Lx - 1; ++i)
    //{
    //    m_lstAllRes.AddItem(m_pHostDataBuffer[i].x);
    //    if (m_bShowResult)
    //    {
    //        appDetailed(_T("%d=(%1.6f,%1.6f)   "), i, m_pHostDataBuffer[i].x, m_pHostDataBuffer[i].y);
    //    }
    //}

    //if (m_bShowResult)
    //{
    //    appDetailed(_T("\n\n ----------- Spin ------------- \n"));
    //}
    //checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBufferS, sizeof(CLGComplex) * (_HC_Lx - 1), cudaMemcpyDeviceToHost));

    //for (UINT i = 0; i < _HC_Lx - 1; ++i)
    //{
    //    m_lstAllRes.AddItem(m_pHostDataBuffer[i].x);
    //    if (m_bShowResult)
    //    {
    //        appDetailed(_T("%d=(%1.6f,%1.6f)   "), i, m_pHostDataBuffer[i].x, m_pHostDataBuffer[i].y);
    //    }
    //}
    //if (m_bShowResult)
    //{
    //    appDetailed(_T("\n\n ================================================ \n\n"));
    //}
}

void CMeasureAMomentumJF::SourceSanning(const class CFieldGauge* pGauge, const class CFieldGauge* pCorrespondingStaple, const TArray<CFieldFermion*>& sources, const SSmallInt4& sourceSite)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureAMomentumJF only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    if (12 != sources.Num())
    {
        appCrucial(_T("Wilson Dirac SU3 need 12 sources!\n"));
        return;
    }
    deviceWilsonVectorSU3* pDevicePtr[12];
    for (INT i = 0; i < 12; ++i)
    {
        CFieldFermionWilsonSquareSU3* fermionfield = dynamic_cast<CFieldFermionWilsonSquareSU3*>(sources[i]);
        pDevicePtr[i] = fermionfield->m_pDeviceData;
    }

    dim3 _blocks(1, 1, 1);
    dim3 _thread1(12, 1, 1);

    deviceWilsonVectorSU3** ppDevicePtr;
    checkCudaErrors(cudaMalloc((void**)&ppDevicePtr, sizeof(deviceWilsonVectorSU3*) * 12));
    checkCudaErrors(cudaMemcpy(ppDevicePtr, pDevicePtr, sizeof(deviceWilsonVectorSU3*) * 12, cudaMemcpyHostToDevice));
    BYTE byArrayIdx = static_cast<BYTE>(sourceSite.x - 1);

    _kernel_XDy_yDx << <_blocks, _thread1 >> > (
        pGaugeSU3->m_pDeviceData,
        ppDevicePtr,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        sourceSite,
        CCommonData::m_sCenter,
        byArrayIdx,
        m_byFieldId,
        m_bNaive,
        m_pOperatorDataL
        );

    if (m_bExponential)
    {
        _kernel_FS_Exponential << <_blocks, _thread1 >> > (
            pGaugeSU3->m_pDeviceData,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
            ppDevicePtr,
            sourceSite,
            byArrayIdx,
            m_byFieldId,
            m_pOperatorDataS
            );
    }
    else
    {
        _kernel_FS << <_blocks, _thread1 >> > (
            ppDevicePtr,
            sourceSite,
            byArrayIdx,
            m_byFieldId,
            m_pOperatorDataS
            );
    }

    checkCudaErrors(cudaFree(ppDevicePtr));

    if (sourceSite.x == static_cast<SBYTE>(_HC_Lx) - 1)
    {
        //all sites calculated
        ++m_uiConfigurationCount;
        dim3 _thread2(_HC_Lx - 1, 1, 1);
        Real fKappa = CCommonData::m_fKai;
        _kernel_Trace_JFLS << <_blocks, _thread2 >> > (m_pOperatorDataL, m_pDeviceDataBufferL, fKappa);
        _kernel_Trace_JFLS << <_blocks, _thread2 >> > (m_pOperatorDataS, m_pDeviceDataBufferS, fKappa);

        if (m_bShowResult)
        {
            appDetailed(_T("\n\n ==================== Angular Momentum (%d con)============================ \n\n"), m_uiConfigurationCount);
            appDetailed(_T(" ----------- Orbital (with a -1) ------------- \n"));
        }
        checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBufferL, sizeof(CLGComplex) * (_HC_Lx - 1), cudaMemcpyDeviceToHost));

        for (UINT i = 0; i < _HC_Lx - 1; ++i)
        {
            m_lstAllRes.AddItem(m_pHostDataBuffer[i].x);
            if (m_bShowResult)
            {
                appDetailed(_T("%d=(%1.6f,%1.6f)   "), i, m_pHostDataBuffer[i].x, m_pHostDataBuffer[i].y);
            }
        }

        if (m_bShowResult)
        {
            appDetailed(_T("\n\n ----------- Spin ------------- \n"));
        }
        checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBufferS, sizeof(CLGComplex) * (_HC_Lx - 1), cudaMemcpyDeviceToHost));

        for (UINT i = 0; i < _HC_Lx - 1; ++i)
        {
            m_lstAllRes.AddItem(m_pHostDataBuffer[i].x);
            if (m_bShowResult)
            {
                appDetailed(_T("%d=(%1.6f,%1.6f)   "), i, m_pHostDataBuffer[i].x, m_pHostDataBuffer[i].y);
            }
        }
        if (m_bShowResult)
        {
            appDetailed(_T("\n\n ================================================ \n\n"));
        }
    }
}

void CMeasureAMomentumJF::Average(UINT )
{
    //nothing to do
}

void CMeasureAMomentumJF::Report()
{
    appSetLogDate(FALSE);

    assert(m_uiConfigurationCount * (_HC_Lx - 1) * 2 == static_cast<UINT>(m_lstAllRes.Num()));
    TArray<Real> tmpSum;

    appGeneral(_T("\n\n==========================================================================\n"));
    appGeneral(_T("==================== Angular Momentum (%d con)============================\n"), m_uiConfigurationCount);
    appGeneral(_T("==========================================================================\n"));
    appGeneral(_T("\n ----------- Orbital (with a -1) ------------- \n"));

    appGeneral(_T("{\n"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        appGeneral(_T("{"));
        for (UINT j = 0; j < _HC_Lx - 1; ++j)
        {
            appGeneral(_T("%2.12f, "), m_lstAllRes[(i * 2) * (_HC_Lx - 1) + j]);
            if (0 == i)
            {
                tmpSum.AddItem(m_lstAllRes[(i * 2) * (_HC_Lx - 1) + j]);
            }
            else
            {
                tmpSum[j] += m_lstAllRes[(i * 2) * (_HC_Lx - 1) + j];
            }
        }
        appGeneral(_T("},\n"));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n ----------- Orbital average (with a -1) ------------- \n"));
    for (UINT j = 0; j < _HC_Lx - 1; ++j)
    {
        appGeneral(_T("%2.12f, "), tmpSum[j] / m_uiConfigurationCount);
    }
    tmpSum.RemoveAll();

    appGeneral(_T("\n -------------------------------------------- \n"));
    appGeneral(_T("\n\n ----------- Spin ------------- \n"));

    appGeneral(_T("{\n"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        appGeneral(_T("{"));
        for (UINT j = 0; j < _HC_Lx - 1; ++j)
        {
            appGeneral(_T("%2.12f, "), m_lstAllRes[(i * 2 + 1) * (_HC_Lx - 1) + j]);
            if (0 == i)
            {
                tmpSum.AddItem(m_lstAllRes[(i * 2 + 1) * (_HC_Lx - 1) + j]);
            }
            else
            {
                tmpSum[j] += m_lstAllRes[(i * 2 + 1) * (_HC_Lx - 1) + j];
            }
        }
        appGeneral(_T("},\n"));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n ----------- Spin average ------------- \n"));
    for (UINT j = 0; j < _HC_Lx - 1; ++j)
    {
        appGeneral(_T("%2.12f, "), tmpSum[j] / m_uiConfigurationCount);
    }
    tmpSum.RemoveAll();

    appGeneral(_T("\n==========================================================================\n"));
    appGeneral(_T("==========================================================================\n\n"));

    appSetLogDate(TRUE);
}

void CMeasureAMomentumJF::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstAllRes.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================