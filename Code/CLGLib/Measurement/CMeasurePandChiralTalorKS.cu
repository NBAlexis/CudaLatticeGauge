//=============================================================================
// FILENAME : CMeasurePandChiralTalorKS.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [16/11/2022 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasurePandChiralTalorKS)

#pragma region kernels

/**
 * This is nothing but remove omega of staggered fermion rotation
 * Those codes are just copy from CFieldFermionKSSU3R
 */
/**
* When link n and n+mu, the coordinate is stick with n
* When link n and n-mu, the coordinate is stick with n-mu
* Irrelavent with tau
* Optimization: bXorY removed, block.x *= 2 
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTerm_NoOmega(
    const deviceSU3Vector * __restrict__ pDeviceData,
    const deviceSU3 * __restrict__ pGauge,
    const BYTE * __restrict__ pEtaTable,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    SSmallInt4 sCenter)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
    //const INT eta_tau = ((pEtaTable[uiSiteIndex] >> 3) & 1);
    const INT eta_tau = pEtaTable[uiSiteIndex] >> 3;

    #pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusMu  = idx & 2;
        const UBOOL bPlusTau = idx & 4;
        //x or y, and y or x is the derivate, not coefficient
        const UINT bXorY = idx & 1;
        const UINT bYorX = 1 - bXorY;
        SSmallInt4 sTargetSite = sSite4;
        SSmallInt4 sMidSite = sSite4;
        sTargetSite.m_byData4[bYorX] = sTargetSite.m_byData4[bYorX] + (bPlusMu ? 2 : -2);
        sMidSite.m_byData4[bYorX] = sMidSite.m_byData4[bYorX] + (bPlusMu ? 1 : -1);
        sTargetSite.w = sTargetSite.w + (bPlusTau ? 1 : -1);
        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sTargetSite)];
        const SIndex& sMiddleBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sMidSite)];
        sMidSite = __deviceSiteIndexToInt4(sMiddleBigIndex.m_uiSiteIndex);

        //note that bYorX = 1, it is x partial_y term, therefore is '-'
        //INT this_eta_tau = (bPlusTau ? eta_tau : ((pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3) & 1))
        INT this_eta_tau = (bPlusTau ? eta_tau : (pEtaTable[sTargetBigIndex.m_uiSiteIndex] >> 3))
                         + bYorX;

        if (sTargetBigIndex.NeedToOpposite())
        {            
            this_eta_tau = this_eta_tau + 1;
        }

        deviceSU3Vector right = _deviceVXXTauOptimized(pGauge, sSite4, byGaugeFieldId, bXorY, bPlusMu, bPlusTau).MulVector(
            pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        //when bXorY = 1, it is y partial _x, so is [1]
        //when bXorY = 0, it is x partial _y, so is [0]
        right.MulReal(sMidSite.m_byData4[bXorY] - sCenter.m_byData4[bXorY] + F(0.5));

        if (!bPlusMu)
        {
            //for -2x, -2y terms, there is another minus sign
            this_eta_tau = this_eta_tau + 1;
        }

        if (this_eta_tau & 1)
        {
            result.Sub(right);
        }
        else
        {
            result.Add(right);
        }
    }

    //if (bDDagger)
    //{
    //    result.MulReal(F(-0.25));
    //}
    //else
    //{
        result.MulReal(F(0.25));
    //}

    //switch (eCoeff)
    //{
    //case EOCT_Real:
    //    result.MulReal(fCoeff);
    //    break;
    //case EOCT_Complex:
    //    result.MulComp(cCoeff);
    //    break;
    //}

    //change to equal
    pResultData[uiSiteIndex] = result;
}


__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTau_Term_NoOmega(
    const deviceSU3Vector* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    deviceSU3Vector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId)
{
    intokernalInt4;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();

    #pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusX = (0 != (idx & 1));
        const UBOOL bPlusY = (0 != (idx & 2));
        const UBOOL bPlusT = (0 != (idx & 4));

        SSmallInt4 sOffset = sSite4;
        sOffset.x = sOffset.x + (bPlusX ? 1 : -1);
        sOffset.y = sOffset.y + (bPlusY ? 1 : -1);
        sOffset.w = sOffset.w + (bPlusT ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];
        
        const deviceSU3Vector right = _deviceVXYTOptimized(pGauge, sSite4, byGaugeFieldId, bPlusX, bPlusY, bPlusT)
        .MulVector(pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta124 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta124 = bPlusT ? (sSite4.y + sSite4.z) : (site_target.y + site_target.z + 1);

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta124 = eta124 + 1;
        }

        if (eta124 & 1)
        {
            result.Sub(right);
        }
        else
        {
            result.Add(right);
        }
    }

    //if (bDDagger)
    //{
    //    result.MulReal(-F(0.125));
    //}
    //else
    //{
        result.MulReal(F(0.125));
    //}

    //switch (eCoeff)
    //{
    //case EOCT_Real:
    //    result.MulReal(fCoeff);
    //    break;
    //case EOCT_Complex:
    //    result.MulComp(cCoeff);
    //    break;
    //}

    pResultData[uiSiteIndex].Add(result);
}

/**
 * just a copy
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelPolyakovLoopOfSiteTalor2(
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
 * just copy from CActionGaugePlaquetteRotating
 */

#pragma region Clover

__global__ void _CLG_LAUNCH_BOUND
_kernelAdd4PlaqutteTermSU3_Shifted_NoOmegaSq(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN,
    DOUBLE* results
#else
    Real betaOverN,
    Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);

#if !_CLG_DOUBLEFLOAT
    DOUBLE fXSq = (sSite4.x - _DC_Centerx + 0.5);
    fXSq = fXSq * fXSq;
    DOUBLE fYSq = (sSite4.y - _DC_Centery + 0.5);
    fYSq = fYSq * fYSq;

    //======================================================
    //4-plaqutte terms
    //Omega^2 x^2 Retr[1 - U_2,3]
    const DOUBLE fU23 = fXSq * _device4PlaqutteTerm(pDeviceData, 1, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 y^2 Retr[1 - U_1,3]
    const DOUBLE fU13 = fYSq * _device4PlaqutteTerm(pDeviceData, 0, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 (x^2 + y^2) Retr[1 - U_1,2]
    const DOUBLE fU12 = (fXSq + fYSq) * _device4PlaqutteTerm(pDeviceData, 0, 1, uiBigIdx, sSite4, byFieldId);
#else
    Real fXSq = (sSite4.x - _DC_Centerx + F(0.5));
    fXSq = fXSq * fXSq;
    Real fYSq = (sSite4.y - _DC_Centery + F(0.5));
    fYSq = fYSq * fYSq;

    //======================================================
    //4-plaqutte terms
    //Omega^2 x^2 Retr[1 - U_2,3]
    const Real fU23 = fXSq * _device4PlaqutteTerm(pDeviceData, 1, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 y^2 Retr[1 - U_1,3]
    const Real fU13 = fYSq * _device4PlaqutteTerm(pDeviceData, 0, 2, uiBigIdx, sSite4, byFieldId);

    //Omega^2 (x^2 + y^2) Retr[1 - U_1,2]
    const Real fU12 = (fXSq + fYSq) * _device4PlaqutteTerm(pDeviceData, 0, 1, uiBigIdx, sSite4, byFieldId);
#endif

    results[uiSiteIndex] = (fU23 + fU13 + fU12) * betaOverN;
}

#pragma endregion

#pragma region Chair energy

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term12_Shifted_NoOmega(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN,
    DOUBLE* results
#else
    Real betaOverN,
    Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

#if !_CLG_DOUBLEFLOAT
    betaOverN = -0.125 * betaOverN;
    const DOUBLE fXOmega = (sSite4.x - _DC_Centerx + 0.5);

    //===============
    //- x Omega V412
    const DOUBLE fV412 = fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 0, 1, uiN);

    //===============
    //- x Omega V432
    const DOUBLE fV432 = fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 1, uiN);

#else
    betaOverN = -F(0.125) * betaOverN;
    const Real fXOmega = (sSite4.x - _DC_Centerx + F(0.5));

    //===============
    //+x Omega V412
    const Real fV412 = fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 0, 1, uiN);

    //===============
    //+x Omega V432
    const Real fV432 = fXOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 1, uiN);
#endif

    results[uiSiteIndex] = (fV412 + fV432) * betaOverN;
}


__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term34_Shifted_NoOmega(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN,
    DOUBLE* results
#else
    Real betaOverN,
    Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

#if !_CLG_DOUBLEFLOAT
    betaOverN = 0.125 * betaOverN;
    const DOUBLE fYOmega = (sSite4.y - _DC_Centery + 0.5);

    //===============
    //+ y Omega V421
    const DOUBLE fV421 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 1, 0, uiN);

    //===============
    //+ y Omega V431
    const DOUBLE fV431 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 0, uiN);
#else
    betaOverN = F(0.125) * betaOverN;
    const Real fYOmega = (sSite4.y - _DC_Centery + F(0.5));

    //===============
    //-y Omega V421
    const Real fV421 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 1, 0, uiN);

    //===============
    //-y Omega V431
    const Real fV431 = fYOmega * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 3, 2, 0, uiN);
#endif

    results[uiSiteIndex] = (fV421 + fV431) * betaOverN;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAddChairTermSU3_Term5_Shifted_NoOmegaSq(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pDeviceData,
#if !_CLG_DOUBLEFLOAT
    DOUBLE betaOverN,
    DOUBLE* results
#else
    Real betaOverN,
    Real* results
#endif
)
{
    intokernalInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);

#if !_CLG_DOUBLEFLOAT
    betaOverN = -0.125 * betaOverN;
    const DOUBLE fXYOmega2 = (sSite4.x - _DC_Centerx + 0.5) * (sSite4.y - _DC_Centery + 0.5);

    //===============
    //-Omega^2 xy V142
    const DOUBLE fV132 = fXYOmega2 * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);
#else
    betaOverN = -F(0.125) * betaOverN;
    const Real fXYOmega2 = (sSite4.x - _DC_Centerx + F(0.5)) * (sSite4.y - _DC_Centery + F(0.5));

    //===============
    //+Omega^2 xy V142
    const Real fV132 = fXYOmega2 * _deviceChairTerm(pDeviceData, byFieldId, sSite4, 0, 2, 1, uiN);
#endif

    results[uiSiteIndex] = fV132 * betaOverN;
}

#pragma endregion

#pragma endregion


void CMeasurePandChiralTalorKS::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasureStochastic::Initial(pOwner, pLatticeData, param, byId);
    Reset();

    INT iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;
}

void CMeasurePandChiralTalorKS::OnConfigurationAcceptedZ4(
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

    const CFieldFermionKSSU3* pF1W = dynamic_cast<const CFieldFermionKSSU3*>(pZ4);

    CFieldFermionKSSU3* pF2W = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(pInverseZ4->m_byFieldId));
    pInverseZ4->CopyTo(pF2W);
    CFieldFermionKSSU3* pTmp = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(pInverseZ4->m_byFieldId));
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pAcceptGauge);

    preparethread;

    //======= D =========
#if _CLG_DOUBLEFLOAT
    m_cTmpSum[ECPCTTTKS_D] = _cuCaddf(m_cTmpSum[ECPCTTTKS_D], pF1W->Dot(pF2W));
#else
    m_cTmpSum[ECPCTTTKS_D] = cuCadd(m_cTmpSum[ECPCTTTKS_D], pF1W->Dot(pF2W));
#endif

    //======= MD =========
    pF2W->CopyTo(pTmp);
    ApplyM(pF2W, pTmp, pGaugeSU3);


#if _CLG_DOUBLEFLOAT
    m_cTmpSum[ECPCTTTKS_MD] = _cuCaddf(m_cTmpSum[ECPCTTTKS_MD], pF1W->Dot(pF2W));
#else
    m_cTmpSum[ECPCTTTKS_MD] = cuCadd(m_cTmpSum[ECPCTTTKS_MD], pF1W->Dot(pF2W));
#endif

    //======= DMD =========

    pF2W->InverseD(pAcceptGauge);


#if _CLG_DOUBLEFLOAT
    m_cTmpSum[ECPCTTTKS_DMD] = _cuCaddf(m_cTmpSum[ECPCTTTKS_DMD], pF1W->Dot(pF2W));
#else
    m_cTmpSum[ECPCTTTKS_DMD] = cuCadd(m_cTmpSum[ECPCTTTKS_DMD], pF1W->Dot(pF2W));
#endif

    //======= MDMD =========

    pF2W->CopyTo(pTmp);
    ApplyM(pF2W, pTmp, pGaugeSU3);


#if _CLG_DOUBLEFLOAT
    m_cTmpSum[ECPCTTTKS_MDMD] = _cuCaddf(m_cTmpSum[ECPCTTTKS_MDMD], pF1W->Dot(pF2W));
#else
    m_cTmpSum[ECPCTTTKS_MDMD] = cuCadd(m_cTmpSum[ECPCTTTKS_MDMD], pF1W->Dot(pF2W));
#endif

    //======= DMDMD =========

    pF2W->InverseD(pAcceptGauge);


#if _CLG_DOUBLEFLOAT
    m_cTmpSum[ECPCTTTKS_DMDMD] = _cuCaddf(m_cTmpSum[ECPCTTTKS_DMDMD], pF1W->Dot(pF2W));
#else
    m_cTmpSum[ECPCTTTKS_DMDMD] = cuCadd(m_cTmpSum[ECPCTTTKS_DMDMD], pF1W->Dot(pF2W));
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

        for (UINT i = 0; i < ECPCTTTKS_Max; ++i)
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

void CMeasurePandChiralTalorKS::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
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
    _kernelPolyakovLoopOfSiteTalor2 << <blockxyz, threadsxyz >> > (
        pGaugeSU3->m_pDeviceData,
        _D_ComplexThreadBuffer,
        _D_RealThreadBuffer
        );

#if _CLG_DOUBLEFLOAT
    CLGComplex polyakovSum = appGetCudaHelper()->ReduceComplex(_D_ComplexThreadBuffer, _HC_Volume_xyz);
#else
    cuDoubleComplex polyakovSum = appGetCudaHelper()->ReduceComplex(_D_ComplexThreadBuffer, _HC_Volume_xyz);
#endif

    m_lstPolyakov.AddItem(polyakovSum);

    //=========== Calculate Omega term ================
    _kernelAddChairTermSU3_Term12_Shifted_NoOmega << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        CCommonData::m_fBeta,
        _D_RealThreadBuffer);
    DOUBLE omegaterm = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    _kernelAddChairTermSU3_Term34_Shifted_NoOmega << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        CCommonData::m_fBeta,
        _D_RealThreadBuffer);
    omegaterm += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    m_lstPolyakovSOmega.AddItem(omegaterm);

    //=========== Calculate Omega Squire term ================
    _kernelAdd4PlaqutteTermSU3_Shifted_NoOmegaSq << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        CCommonData::m_fBeta,
        _D_RealThreadBuffer);

    DOUBLE omegasqterm = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);


    _kernelAddChairTermSU3_Term5_Shifted_NoOmegaSq << <block, threads >> > (
        pGaugeSU3->m_byFieldId,
        pGaugeSU3->m_pDeviceData,
        CCommonData::m_fBeta,
        _D_RealThreadBuffer);
    omegasqterm += appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);

    m_lstPolyakovSOmegaSq.AddItem(omegasqterm);
}


void CMeasurePandChiralTalorKS::ApplyM(CFieldFermionKSSU3* pTarget, const CFieldFermionKSSU3* pSource, const CFieldGaugeSU3* pGauge)
{
    preparethread;
    _kernelDFermionKS_PR_XYTerm_NoOmega << <block, threads >> > (
        pSource->m_pDeviceData,
        pGauge->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget->m_pDeviceData,
        pSource->m_byFieldId,
        pGauge->m_byFieldId,
        _HC_Center);


    _kernelDFermionKS_PR_XYTau_Term_NoOmega << <block, threads >> > (
        pSource->m_pDeviceData,
        pGauge->m_pDeviceData,
        pTarget->m_pDeviceData,
        pSource->m_byFieldId,
        pGauge->m_byFieldId);

}

void CMeasurePandChiralTalorKS::Report()
{
    appPushLogDate(FALSE);
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
    appPopLogDate();
}

void CMeasurePandChiralTalorKS::Reset()
{
    CMeasureStochastic::Reset();
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