//=============================================================================
// FILENAME : CMeasurePandChiralTalor.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [11/26/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasurePandChiralTalor)

#pragma region kernels

/**
 * This is kappa * g4.(-y Dx + x Dy - i sigma12)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelTraceApplyM(
    const deviceWilsonVectorSU3* __restrict__ pRight,
    deviceWilsonVectorSU3* pRes,
    SSmallInt4 sCenter,
    const DOUBLE fKappa,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    BYTE byFieldId)
{
    intokernalInt4;

    
    const gammaMatrix& sigma12 = __chiralGamma[SIGMA12E];
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    //const SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    const Real fmY = -static_cast<Real>((sSite4.y - sCenter.y) * fKappa);
    const Real fmX = -static_cast<Real>((sSite4.x - sCenter.x) * fKappa);

    deviceWilsonVectorSU3 jl = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //this is now + y - x
    #pragma unroll
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
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

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
                jl.Add(u_phi_x_p_m);
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_phi_x_p_m.MulReal(fmX);
                jl.Sub(u_phi_x_p_m);
            }
        }
        else
        {
            if (0 == idir)
            {
                //- k y Omega x_p_m
                //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
                u_phi_x_p_m.MulReal(fmY);
                jl.Sub(u_phi_x_p_m);
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_phi_x_p_m.MulReal(fmX);
                jl.Add(u_phi_x_p_m);
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
                jl.Sub(u_dagger_phi_x_m_m);
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_dagger_phi_x_m_m.MulReal(fmX);
                jl.Add(u_dagger_phi_x_m_m);
            }
        }
        else
        {
            if (0 == idir)
            {
                //- k y Omega x_p_m
                //+ k y Omega gamma4 U(x,mu) phi(x+ mu)
                u_dagger_phi_x_m_m.MulReal(fmY);
                jl.Add(u_dagger_phi_x_m_m);
            }
            else if (1 == idir)
            {
                //+ k x Omega x_p_m
                //- k x Omega gamma4 U(x,mu) phi(x+ mu)
                u_dagger_phi_x_m_m.MulReal(fmX);
                jl.Sub(u_dagger_phi_x_m_m);
            }
        }
    }

    //=========================================
    //-i * sigma_12
    pRes[uiSiteIndex] = sigma12.MulWilsonC(pRight[uiSiteIndex]);
    pRes[uiSiteIndex].MulComp(_make_cuComplex(F(0.0), static_cast<Real>(-fKappa)));

    // -y + x -i * sigma_12
    pRes[uiSiteIndex].Sub(jl);

    // gamma4(-y + x -i * sigma_12)
    const gammaMatrix& gamma4 = __chiralGamma[GAMMA4];
    pRes[uiSiteIndex] = gamma4.MulWilsonC(pRes[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovLoopOfSiteTalor(
    const deviceSU3* __restrict__ pDeviceBuffer,
#if _CLG_DOUBLEFLOAT
    CLGComplex* res, Real* sumCount
#else
    cuDoubleComplex* res, DOUBLE* sumCount
#endif
)
{
    UINT uiXYZ = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    deviceSU3 beforeTrace;
#if _CLG_DOUBLEFLOAT
    sumCount[uiXYZ] = F(0.0);
#else
    sumCount[uiXYZ] = 0.0;
#endif
    for (UINT t = 0; t < _DC_Lt; ++t)
    {
        const UINT uiSiteIndex = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + t;
        UINT uiLinkIdx = _deviceGetLinkIndex(uiSiteIndex, 3);
        const SSmallInt4 site4 = __deviceSiteIndexToInt4(uiSiteIndex);
        const UINT uiBigIdx = __idx->_deviceGetBigIndex(site4);

        if (0 == t)
        {
            if (__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
            {
                beforeTrace = deviceSU3::makeSU3Zero();
            }
            else
            {
#if _CLG_DOUBLEFLOAT
                sumCount[uiXYZ] = F(1.0);
#else
                sumCount[uiXYZ] = 1.0;
#endif
                beforeTrace = pDeviceBuffer[uiLinkIdx];
            }
        }
        else
        {
            if (!__idx->_deviceIsBondOnSurface(uiBigIdx, 3))
            {
                beforeTrace.Mul(pDeviceBuffer[uiLinkIdx]);
            }
        }
    }
#if _CLG_DOUBLEFLOAT
    res[uiXYZ] = beforeTrace.Tr();
#else
    res[uiXYZ] = _cToDouble(beforeTrace.Tr());
#endif
}

/**
 * x (V412 + V432) - y (V421 + V431)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelActionTalorOmega(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    SSmallInt4 sCenterSite,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN, DOUBLE* results
#else
    Real betaOverN, Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet())
    {
        results[uiSiteIndex] = F(0.0);
        return;
    }

    betaOverN = F(0.125) * betaOverN;
    const Real fXOmega = (sSite4.x - sCenterSite.x);

    //===============
    //+x Omega V412
    const Real fV412 = fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 0, 1, uiN);

    //===============
    //+x Omega V432
    const Real fV432 = fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 1, uiN);


    const Real fYOmega = -(sSite4.y - sCenterSite.y);

    //===============
    //-y Omega V421
    const Real fV421 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 1, 0, uiN);

    //===============
    //-y Omega V431
    const Real fV431 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 0, uiN);

    results[uiSiteIndex] = (fV412 + fV432 + fV421 + fV431) * betaOverN;
}

/**
 * -(x^2 U23 + y^2 U13 + r^2 U12 - xy V132)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelActionTalorOmegaSq(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
    const SIndex* __restrict__ pCachedPlaqutte,
    SSmallInt4 sCenterSite,
#if !_CLG_DOUBLEFLOAT
    const DOUBLE betaOverN, DOUBLE* results
#else
    const Real betaOverN, Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    const UINT plaqLength = __idx->m_pSmallData[CIndexData::kPlaqLengthIdx];
    const UINT plaqCountAll = __idx->m_pSmallData[CIndexData::kPlaqPerSiteIdx] * plaqLength;

#if !_CLG_DOUBLEFLOAT
    DOUBLE res = 0.0;
#else
    Real res = F(0.0);
#endif
    #pragma unroll
    for (BYTE idx0 = 0; idx0 < 3; ++idx0)
    {
        //i=0: 12
        //  1: 13
        //  2: 14
        //  3: 23
        //  4: 24
        //  5: 34
        //0->0, 1->1, 2->3
        //0-> r^2, 1->y^2, 2(or 3)-> x^2
        const BYTE idx = (2 == idx0) ? (idx0 + 1) : idx0;

        //Real resThisThread = F(0.0);

        //========================================
        //find plaqutte 1-4, or 2-4, or 3-4
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

        //0 -> xy, 1 -> xz, 2 -> yz
        //x x y
        const BYTE mushift = (idx0 / 2);
        //y z z
        const BYTE nushift = ((idx0 + 1) / 2) + 1;

        //the Fi function is: 0 : x^2 + y^2, 1 : y^2, 2 : x^2
        //r^2 (1-U_{xy}) + y^2 U_{xz} + x^2 U_{yz}

#if !_CLG_DOUBLEFLOAT
        res += static_cast<DOUBLE>(betaOverN * (3.0 - toAdd.ReTr()) * _deviceFi(byFieldId, sSite4, sCenterSite, uiN, idx0, mushift, nushift));
#else
        res += betaOverN * (F(3.0) - toAdd.ReTr()) * _deviceFi(byFieldId, sSite4, sCenterSite, uiN, idx0, mushift, nushift);
#endif
    }

    if (!__idx->m_pDeviceIndexPositionToSIndex[1][uiN].IsDirichlet())
    {
        const Real fXYOmega2 = -(sSite4.x - sCenterSite.x) * (sSite4.y - sCenterSite.y);

        //===============
        //+Omega^2 xy V132
        const Real fV132 = fXYOmega2 * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);

        res += fV132 * F(0.125) * betaOverN;
    }

    results[uiSiteIndex] = -res;
}

#pragma endregion

CMeasurePandChiralTalor::~CMeasurePandChiralTalor()
{

}

void CMeasurePandChiralTalor::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasureStochastic::Initial(pOwner, pLatticeData, param, byId);
    Reset();

    INT iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;
}

void CMeasurePandChiralTalor::OnConfigurationAcceptedZ4(
    const class CFieldGauge* pAcceptGauge,
    const class CFieldGauge* pCorrespondingStaple,
    const class CFieldFermion* pZ4,
    const class CFieldFermion* pInverseZ4,
    UBOOL bStart,
    UBOOL bEnd)
{
    if (bStart)
    {
        for (UINT i = 0; i < ECPCTTT_Max; ++i)
        {
#if _CLG_DOUBLEFLOAT
            m_cTmpSum[i] = _zeroc;
#else
            m_cTmpSum[i] = make_cuDoubleComplex(0.0, 0.0);
#endif
        }
    }

    const CFieldFermionWilsonSquareSU3* pF1W = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(pZ4);

    CFieldFermionWilsonSquareSU3* pF2W = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(pInverseZ4->m_byFieldId));
    pInverseZ4->CopyTo(pF2W);
    CFieldFermionWilsonSquareSU3* pTmp = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(pInverseZ4->m_byFieldId));
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pAcceptGauge);

    preparethread;

    //======= D =========
#if _CLG_DOUBLEFLOAT
    m_cTmpSum[ECPCTTT_D] = _cuCaddf(m_cTmpSum[ECPCTTT_D], pF1W->Dot(pF2W));
#else
    m_cTmpSum[ECPCTTT_D] = cuCadd(m_cTmpSum[ECPCTTT_D], pF1W->Dot(pF2W));
#endif

    //======= MD =========
    pF2W->CopyTo(pTmp);
    _kernelTraceApplyM << <block, threads >> > (
        pTmp->m_pDeviceData,
        pF2W->m_pDeviceData,
        CCommonData::m_sCenter,
        CCommonData::m_fKai,
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        m_byFieldId
        );


#if _CLG_DOUBLEFLOAT
    m_cTmpSum[ECPCTTT_MD] = _cuCaddf(m_cTmpSum[ECPCTTT_MD], pF1W->Dot(pF2W));
#else
    m_cTmpSum[ECPCTTT_MD] = cuCadd(m_cTmpSum[ECPCTTT_MD], pF1W->Dot(pF2W));
#endif

    //======= DMD =========

    pF2W->InverseD(pAcceptGauge);


#if _CLG_DOUBLEFLOAT
    m_cTmpSum[ECPCTTT_DMD] = _cuCaddf(m_cTmpSum[ECPCTTT_DMD], pF1W->Dot(pF2W));
#else
    m_cTmpSum[ECPCTTT_DMD] = cuCadd(m_cTmpSum[ECPCTTT_DMD], pF1W->Dot(pF2W));
#endif

    //======= MDMD =========

    pF2W->CopyTo(pTmp);
    _kernelTraceApplyM << <block, threads >> > (
        pTmp->m_pDeviceData,
        pF2W->m_pDeviceData,
        CCommonData::m_sCenter,
        CCommonData::m_fKai,
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        m_byFieldId
        );


#if _CLG_DOUBLEFLOAT
    m_cTmpSum[ECPCTTT_MDMD] = _cuCaddf(m_cTmpSum[ECPCTTT_MDMD], pF1W->Dot(pF2W));
#else
    m_cTmpSum[ECPCTTT_MDMD] = cuCadd(m_cTmpSum[ECPCTTT_MDMD], pF1W->Dot(pF2W));
#endif

    //======= DMDMD =========

    pF2W->InverseD(pAcceptGauge);


#if _CLG_DOUBLEFLOAT
    m_cTmpSum[ECPCTTT_DMDMD] = _cuCaddf(m_cTmpSum[ECPCTTT_DMDMD], pF1W->Dot(pF2W));
#else
    m_cTmpSum[ECPCTTT_DMDMD] = cuCadd(m_cTmpSum[ECPCTTT_DMDMD], pF1W->Dot(pF2W));
#endif

    pF2W->Return();
    pTmp->Return();

    if (bEnd)
    {
#if _CLG_DOUBLEFLOAT
        const Real fDiv2 = F(1.0) / m_uiFieldCount;
#else
        const DOUBLE fDiv2 = 1.0 / m_uiFieldCount;
#endif

        for (UINT i = 0; i < ECPCTTT_Max; ++i)
        {
#if _CLG_DOUBLEFLOAT
            m_cTmpSum[i] = cuCmulf_cr(m_cTmpSum[i], fDiv2);
#else
            m_cTmpSum[i] = cuCmulf_cd(m_cTmpSum[i], fDiv2);
#endif
            appDetailed(_T("\n Condensate %d = %2.12f + %2.12f\n"), i, m_cTmpSum[i].x, m_cTmpSum[i].y);
            m_lstTraceRes[i].AddItem(m_cTmpSum[i]);
        }

        ++m_uiConfigurationCount;
    }
}

void CMeasurePandChiralTalor::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    //const BYTE byGaugeFiledId = pGaugeSU3->m_byFieldId;
    preparethread;
    const dim3 blockxyz(_HC_DecompX, _HC_DecompY, 1); 
    const dim3 threadsxyz(_HC_DecompLx, _HC_DecompLy, 1);

    //=========== Calculate Polyakov loop ================
    _kernelPolyakovLoopOfSiteTalor << <blockxyz, threadsxyz >> > (
        pGaugeSU3->m_pDeviceData,
        _D_ComplexThreadBuffer,
        _D_RealThreadBuffer
        );

#if _CLG_DOUBLEFLOAT
    CLGComplex polyakovSum = appGetCudaHelper()->ReduceComplex(_D_ComplexThreadBuffer, _HC_Volume_xyz);
    Real polyakovSiteCount = appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz);
    polyakovSum = cuCdivf_cr_host(polyakovSum, polyakovSiteCount);
#else
    cuDoubleComplex polyakovSum = appGetCudaHelper()->ReduceComplex(_D_ComplexThreadBuffer, _HC_Volume_xyz);
    DOUBLE polyakovSiteCount = appGetCudaHelper()->ReduceReal(_D_RealThreadBuffer, _HC_Volume_xyz);
    polyakovSum = cuCdivf_cd_host(polyakovSum, polyakovSiteCount);
#endif

    m_lstPolyakov.AddItem(polyakovSum);

    //=========== Calculate Omega term ================
    _kernelActionTalorOmega << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        CCommonData::m_sCenter,
        CCommonData::m_fBeta,
        _D_RealThreadBuffer);

#if _CLG_DOUBLEFLOAT
    const Real omegaterm = appGetCudaHelper()->ReduceRealWithThreadCount(_D_RealThreadBuffer);
#else
    const DOUBLE omegaterm = appGetCudaHelper()->ReduceRealWithThreadCount(_D_RealThreadBuffer);
#endif

    m_lstPolyakovSOmega.AddItem(omegaterm);

    //=========== Calculate Omega Squire term ================
    _kernelActionTalorOmegaSq << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pPlaqutteCache,
        CCommonData::m_sCenter,
        CCommonData::m_fBeta,
        _D_RealThreadBuffer);

#if _CLG_DOUBLEFLOAT
    const Real omegasqterm = appGetCudaHelper()->ReduceRealWithThreadCount(_D_RealThreadBuffer);
#else
    const DOUBLE omegasqterm = appGetCudaHelper()->ReduceRealWithThreadCount(_D_RealThreadBuffer);
#endif

    m_lstPolyakovSOmegaSq.AddItem(omegasqterm);
}

void CMeasurePandChiralTalor::Average(UINT)
{
    //nothing to do
}

void CMeasurePandChiralTalor::Report()
{
    for (UINT i = 0; i < ECPCTTT_Max; ++i)
    {
        assert(m_uiConfigurationCount == static_cast<UINT>(m_lstTraceRes[i].Num()));

        appGeneral(_T("\n==========================================================================\n"));
        appGeneral(_T("==================== Traces No %d (%d con)============================\n"), i, m_uiConfigurationCount);

#if _CLG_DOUBLEFLOAT
        CLGComplex tmpChargeSum = _zeroc;
#else
        cuDoubleComplex tmpChargeSum = make_cuDoubleComplex(0.0, 0.0);
#endif
        if (m_uiConfigurationCount > 1)
        {
            appGeneral(_T("\n ----------- each configuration ------------- \n"));
            appGeneral(_T("{"));

            for (UINT j = 0; j < m_uiConfigurationCount; ++j)
            {
                tmpChargeSum.x += m_lstTraceRes[i][j].x;
                tmpChargeSum.y += m_lstTraceRes[i][j].y;

                LogGeneralComplex(m_lstTraceRes[i][j]);
            }
            appGeneral(_T("}\n"));

            tmpChargeSum.x = tmpChargeSum.x / m_uiConfigurationCount;
            tmpChargeSum.y = tmpChargeSum.y / m_uiConfigurationCount;
            appGeneral(_T("\n ----------- average condensate = %2.12f + %2.12f ------------- \n"),
                tmpChargeSum.x, tmpChargeSum.y);

            //m_cAverageCondensate = tmpChargeSum;
        }
        else
        {
            appGeneral(_T("\n ----------- average condensate = %2.12f + %2.12f ------------- \n"),
                m_lstTraceRes[i][0].x,
                m_lstTraceRes[i][0].y);

            //m_cAverageCondensate = m_lstCondAll[i][0];
        }
    }

    appGeneral(_T("==========================================================================\n"));
    appSetLogDate(TRUE);
}

void CMeasurePandChiralTalor::Reset()
{
    m_uiConfigurationCount = 0;
    for (UINT i = 0; i < ECPCTTT_Max; ++i)
    {
        m_lstTraceRes[i].RemoveAll();
    }

    m_lstPolyakov.RemoveAll();
    m_lstPolyakovSOmega.RemoveAll();
    m_lstPolyakovSOmegaSq.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================