//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3Gamma.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [08/27/2023 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSquareSU3Gamma)

#pragma region kernels

/**
* + kappa * Gamma
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelApplyGammaSU3Gamma(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    deviceWilsonVectorSU3* pResult, 
    UINT uiGamma, 
    Real fGammaCoeff, 
    Real fKappa,
    BYTE byFieldId, 
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    deviceWilsonVectorSU3 result = pDeviceData[uiSiteIndex];
    if (sSite.IsDirichlet())
    {
        return;
    }

    if (bDDagger)
    {
        result = __chiralGamma[GAMMA5].MulWilsonC(result);
    }
    result = __chiralGamma[uiGamma].MulWilsonC(result);
    result.MulReal(fGammaCoeff * fKappa);
    if (GAMMA1 == uiGamma || GAMMA2 == uiGamma || GAMMA3 == uiGamma || GAMMA4 == uiGamma)
    {
        result.MulComp(_make_cuComplex(F(0.0), F(1.0)));
    }

    if (bDDagger)
    {
        result = __chiralGamma[GAMMA5].MulWilsonC(result);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        result.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        result.MulComp(cCoeff);
        break;
    }

    pResult[uiSiteIndex].Add(result);
}

/**
* In the case of exp, it must be simulation, so we do not care imaginary or real
* We use Exp(i Gamma) when Gamma is Hermitian(gamma_mu, sigma), and use Exp(Gamma) when Gamma is anti-Hermitian(gamma5.gamma_mu)
* 
* Exp(i c Gamma) = Cos(c) + i Sin(c) Gamma
* Exp(c Gamma) = Cos(c) + Sin(c) Gamma
* 
* The offset is along the t-direction
* 
* So it is: -kappa phi(n) g4. [ (-exp(Gamma)) U_t phi(n+t) +  (exp(-Gamma)) Ud_t phi(n-t) ]
* 
* It was gamma4.Gamma
* therefore, we only consider thoes not gamma5-Hermitian, which are:
* g4,s14,s24,s34,g54
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelApplyGammaSU3GammaExp(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    Real fGammaCoeff,
    UINT uiGamma,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    const gammaMatrix& gamma5 = __chiralGamma[GAMMA5];
    const gammaMatrix& gamma4 = __chiralGamma[GAMMA4];
    const gammaMatrix& gamma = __chiralGamma[uiGamma];
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);

    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    const deviceSU3& x_Gauge_element = _deviceGetGaugeBCSU3Dir(byGaugeFieldId, pGauge, uiBigIdx, 3);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(byGaugeFieldId, pGauge, x_m_mu_Gauge);
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.Dagger();
    }

    deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
    deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);

    if (bDDagger)
    {
        x_p_mu_Fermion_element = gamma5.MulWilsonC(x_p_mu_Fermion_element);
        x_m_mu_Fermion_element = gamma5.MulWilsonC(x_m_mu_Fermion_element);
    }

    //hopping terms
    const Real ksin = kai * _sin(fGammaCoeff);
    const Real kcos = kai * _cos(fGammaCoeff);

    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
    // cospart = k cos(f) U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 cospart = u_phi_x_p_m.MulRealC(kcos);
    // u_phi_x_p_m =  k sin(f) U(x,mu) phi(x+ mu)
    if (uiGamma == GAMMA4 || uiGamma == GAMMA5 || uiGamma == SIGMA41 || uiGamma == SIGMA42 || uiGamma == SIGMA43 || uiGamma == GAMMA45 || uiGamma == GAMMA54)
    {
        u_phi_x_p_m.MulReal(ksin);
    }
    else
    {
        u_phi_x_p_m.MulComp(_make_cuComplex(F(0.0), ksin));
    }
    // cospart =  k cos(f) U(x,mu) phi(x+ mu) + k sin(f) gamma U(x,mu) phi(x+ mu)
    cospart.Add(gamma.MulWilsonC(u_phi_x_p_m));
    // cospart =  k g4 [cos(f) U(x,mu) phi(x+ mu) + k sin(f) gamma U(x,mu) phi(x+ mu)]
    cospart = gamma4.MulWilsonC(cospart);

    if (x_p_mu_Fermion.NeedToOpposite())
    {
        result.Sub(cospart);
    }
    else
    {
        result.Add(cospart);
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);

    //u_dagger_phi_x_m_m.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));
    // cospart2 = k Cos(f) U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 cospart2 = u_dagger_phi_x_m_m.MulRealC(kcos);
    // u_dagger_phi_x_m_m =  -k Sin(f) U^{dagger}(x-mu) phi(x-mu)
    if (uiGamma == GAMMA4 || uiGamma == GAMMA5 || uiGamma == SIGMA41 || uiGamma == SIGMA42 || uiGamma == SIGMA43 || uiGamma == GAMMA45 || uiGamma == GAMMA54)
    {
        u_dagger_phi_x_m_m.MulReal(-ksin);
    }
    else
    {
        u_dagger_phi_x_m_m.MulComp(_make_cuComplex(F(0.0), -ksin));
    }
    // cospart2 =  k [cos(f) - Sin(f) gamma ] U^{dagger}(x-mu) phi(x-mu)
    cospart2.Add(gamma.MulWilsonC(u_dagger_phi_x_m_m));
    // cospart =  k g4 [cos(f) U(x,mu) phi(x+ mu) - k sin(f) gamma U(x,mu) phi(x+ mu)]
    cospart2 = gamma4.MulWilsonC(cospart2);

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        result.Add(cospart2);
    }
    else
    {
        result.Sub(cospart2);
    }

    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }

    //==============================
    //res = [gamma5 (orig - term3 - result) gamma5] * coeff
    //    = [gamma5 orig gamma5] * coeff - [gamma5 (term3 + result) gamma5] * coeff
    //    = res0 - [gamma5 (term3 + result) gamma5] * coeff
    //term3.Add(result);
    switch (eCoeff)
    {
    case EOCT_Real:
        result.MulReal(fCoeff);
        break;
    case EOCT_Complex:
        result.MulComp(cCoeff);
        break;
    }
    // -kappa, so it is sub
    pResultData[uiSiteIndex].Sub(result);
}

/**
* -kappa phi(n) g4. [ (-exp(Gamma)) U_t phi(n+t) +  (exp(-Gamma)) Ud_t phi(n-t) ]
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelForceGammaExp(
    const deviceWilsonVectorSU3* __restrict__ pInverseD,
    const deviceWilsonVectorSU3* __restrict__ pInverseDDdagger,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    UINT uiGamma,
    deviceSU3* pForce,
    Real fKai,
    Real fCoeff,
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    const gammaMatrix& gamma4 = __chiralGamma[GAMMA4];
    const gammaMatrix& gamma = __chiralGamma[uiGamma];

    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);
    const SIndex& x_p_mu_Fermion = pFermionMove[linkIndex * 2];

    //If one of the sites is on surface, it has no contribution.
    //Note that, the bond on surface is equivelant to both sites on surface.

    if (x_p_mu_Fermion.IsDirichlet() || sSite.IsDirichlet())
    {
        return;
    }

    //====================
    //Get things
    const deviceWilsonVectorSU3& x_Left = pInverseDDdagger[uiSiteIndex];
    deviceWilsonVectorSU3 x_Right(pInverseD[uiSiteIndex]);
    //all not on surface
    deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
    const deviceWilsonVectorSU3& x_p_mu_Left = pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex];

    const deviceSU3& x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

    const Real fFac = (x_p_mu_Fermion.NeedToOpposite() ? F(-1.0) : F(1.0)) * fKai;
    const Real fCos = fFac * _cos(fCoeff);
    const Real fSin = fFac * _sin(fCoeff);

    //-kappa (-exp(Gamma)) U_t phi(n+t) = exp(Gamma) phi(n+t)
    deviceWilsonVectorSU3 x_p_mu_Right_real = x_p_mu_Right.MulRealC(fCos);
    if (uiGamma == GAMMA4 || uiGamma == GAMMA5 || uiGamma == SIGMA41 || uiGamma == SIGMA42 || uiGamma == SIGMA43 || uiGamma == GAMMA45 || uiGamma == GAMMA54)
    {
        x_p_mu_Right.MulReal(fSin);
    }
    else
    {
        x_p_mu_Right.MulComp(_make_cuComplex(F(0.0), fSin));
    }
    x_p_mu_Right = gamma.MulWilsonC(x_p_mu_Right);
    x_p_mu_Right_real.Add(x_p_mu_Right);
    x_p_mu_Right_real = gamma4.MulWilsonC(x_p_mu_Right_real);
    deviceSU3 mid = deviceSU3::makeSU3Contract(x_Left, x_p_mu_Right_real);

    //-kappa (exp(Gamma)) U_t phi(n+t) = -exp(Gamma) phi(n+t)
    deviceWilsonVectorSU3 x_Right_real = x_Right.MulRealC(fCos);
    if (uiGamma == GAMMA4 || uiGamma == GAMMA5 || uiGamma == SIGMA41 || uiGamma == SIGMA42 || uiGamma == SIGMA43 || uiGamma == GAMMA45 || uiGamma == GAMMA54)
    {
        x_Right.MulReal(-fSin);
    }
    else
    {
        x_Right.MulComp(_make_cuComplex(F(0.0), -fSin));
    }
    x_Right = gamma.MulWilsonC(x_Right);
    x_Right_real.Add(x_Right);
    x_Right_real = gamma4.MulWilsonC(x_Right_real);
    mid.Sub(deviceSU3::makeSU3Contract(x_Right_real, x_p_mu_Left));

    deviceSU3 forceOfThisLink = x_Gauge_element.MulC(mid);
    forceOfThisLink.Ta();
    pForce[linkIndex].Add(forceOfThisLink);
}

#pragma endregion

void CFieldFermionWilsonSquareSU3Gamma::appApplyGammaExp(
    void* pTargetBuffer,
    const void* pBuffer,
    const void* pGaugeBuffer,
    const SIndex* pGaugeMove,
    const SIndex* pFermionMove,
    EGammaMatrix eGamma,
    UBOOL bDagger,
    UBOOL bExp,
    Real fGammaCoeff,
    Real fKappa,
    EOperatorCoefficientType eOCT,
    Real fRealCoeff,
    CLGComplex cCmpCoeff,
    BYTE byFieldID,
    BYTE byGaugeFieldID)
{
    deviceWilsonVectorSU3* pTarget = (deviceWilsonVectorSU3*)pTargetBuffer;
    const deviceWilsonVectorSU3* pSource = (const deviceWilsonVectorSU3*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
    
    preparethread;

    if (bExp)
    {
        _kernelApplyGammaSU3GammaExp << <block, threads >> > (
            pSource,
            pGauge,
            pGaugeMove,
            pFermionMove,
            pTarget,
            fKappa,
            fGammaCoeff,
            static_cast<UINT>(eGamma),
            byFieldID,
            byGaugeFieldID,
            bDagger,
            eOCT,
            fRealCoeff,
            cCmpCoeff
            );
    }
    else
    {
        _kernelApplyGammaSU3Gamma << <block, threads >> > (
            pSource,
            pTarget,
            static_cast<UINT>(eGamma),
            fGammaCoeff,
            fKappa,
            byFieldID,
            bDagger,
            eOCT,
            fRealCoeff,
            cCmpCoeff
            );
    }
}

void CFieldFermionWilsonSquareSU3Gamma::GammaForceExp(
    void* pForce,
    const void* pGaugeBuffer,
    const void* InverseD,
    const void* InverseDDagger,
    const SIndex* pFermionMove,
    UBOOL bExp,
    Real fGammaCoeff,
    Real fKappa,
    EGammaMatrix eGamma,
    BYTE byFieldID,
    BYTE byGaugeFieldID)
{
    if (bExp)
    {
        preparethread;

        const deviceWilsonVectorSU3* pInverseD = (const deviceWilsonVectorSU3*)InverseD;
        const deviceWilsonVectorSU3* pInverseDDagger = (const deviceWilsonVectorSU3*)InverseDDagger;
        const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
        deviceSU3* pForceSU3 = (deviceSU3*)pForce;

        _kernelForceGammaExp << <block, threads >> > (
            pInverseD,
            pInverseDDagger,
            pGauge,
            pFermionMove,
            static_cast<UINT>(eGamma),
            pForceSU3,
            fKappa,
            fGammaCoeff,
            byFieldID
            );
    }
}

#pragma region DOperator

void CFieldFermionWilsonSquareSU3Gamma::DOperator(void* pTargetBuffer, const void* pBuffer, const void* pGaugeBuffer, BYTE byGaugeFieldId,
    UBOOL bDagger, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    CFieldFermionWilsonSquareSU3D::DOperator(pTargetBuffer, pBuffer, pGaugeBuffer, byGaugeFieldId, bDagger, eOCT, fRealCoeff, cCmpCoeff);

    if (abs(m_fCoeffGamma1) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer, 
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            GAMMA1, bDagger, m_bExpGamma, m_fCoeffGamma1, m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffGamma2) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            GAMMA2, bDagger, m_bExpGamma, m_fCoeffGamma2, m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffGamma3) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            GAMMA3, bDagger, m_bExpGamma, m_fCoeffGamma3, m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffGamma4) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            GAMMA4, bDagger, m_bExpGamma, m_fCoeffGamma4, m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }

    if (abs(m_fCoeffSigma12) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            SIGMA12, bDagger, m_bExpGamma, m_fCoeffSigma12, m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffSigma13) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            SIGMA31, bDagger, m_bExpGamma, m_fCoeffSigma13 * F(-1.0), m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffSigma14) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            SIGMA41, bDagger, m_bExpGamma, m_fCoeffSigma14 * F(-1.0), m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffSigma23) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            SIGMA23, bDagger, m_bExpGamma, m_fCoeffSigma23, m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffSigma24) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            SIGMA42, bDagger, m_bExpGamma, m_fCoeffSigma24 * F(-1.0), m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffSigma34) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            SIGMA43, bDagger, m_bExpGamma, m_fCoeffSigma34 * F(-1.0), m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }

    if (abs(m_fCoeffGamma51) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            GAMMA51, bDagger, m_bExpGamma, m_fCoeffGamma51, m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffGamma52) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            GAMMA52, bDagger, m_bExpGamma, m_fCoeffGamma52, m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffGamma53) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            GAMMA53, bDagger, m_bExpGamma, m_fCoeffGamma53, m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
    if (abs(m_fCoeffGamma54) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            GAMMA54, bDagger, m_bExpGamma, m_fCoeffGamma54, m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }

    if (abs(m_fCoeffGamma5) > _CLG_FLT_EPSILON)
    {
        appApplyGammaExp(pTargetBuffer, pBuffer, pGaugeBuffer,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            GAMMA5, bDagger, m_bExpGamma, m_fCoeffGamma5, m_fKai, eOCT, fRealCoeff, cCmpCoeff, m_byFieldId, byGaugeFieldId);
    }
}

void CFieldFermionWilsonSquareSU3Gamma::DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer, BYTE byGaugeFieldId) const
{
    preparethread;
    CFieldFermionWilsonSquareSU3D::DerivateDOperator(pForce, pDphi, pDDphi, pGaugeBuffer, byGaugeFieldId);

    if (m_bExpGamma)
    {
        if (abs(m_fCoeffGamma1) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffGamma1, m_fKai, GAMMA1, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffGamma2) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffGamma2, m_fKai, GAMMA2, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffGamma3) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffGamma3, m_fKai, GAMMA3, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffGamma4) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffGamma4, m_fKai, GAMMA4, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffSigma12) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffSigma12, m_fKai, SIGMA12, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffSigma13) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffSigma13 * F(-1.0), m_fKai, SIGMA31, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffSigma14) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffSigma14 * F(-1.0), m_fKai, SIGMA41, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffSigma23) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffSigma23, m_fKai, SIGMA23, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffSigma24) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffSigma24 * F(-1.0), m_fKai, SIGMA42, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffSigma34) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffSigma34 * F(-1.0), m_fKai, SIGMA43, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffGamma51) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffGamma51, m_fKai, GAMMA51, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffGamma52) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffGamma52, m_fKai, GAMMA52, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffGamma53) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffGamma53, m_fKai, GAMMA53, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffGamma54) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffGamma54, m_fKai, GAMMA54, m_byFieldId, byGaugeFieldId);
        }

        if (abs(m_fCoeffGamma5) > _CLG_FLT_EPSILON)
        {
            GammaForceExp(pForce, pGaugeBuffer, pDphi, pDDphi,
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                TRUE, m_fCoeffGamma5, m_fKai, GAMMA5, m_byFieldId, byGaugeFieldId);
        }
    }
}

#pragma endregion

CFieldFermionWilsonSquareSU3Gamma::CFieldFermionWilsonSquareSU3Gamma()
    : CFieldFermionWilsonSquareSU3D()
    , m_bExpGamma(TRUE)
    , m_fCoeffGamma1(F(0.0))
    , m_fCoeffGamma2(F(0.0))
    , m_fCoeffGamma3(F(0.0))
    , m_fCoeffGamma4(F(0.0))
    , m_fCoeffGamma5(F(0.0))
    , m_fCoeffGamma51(F(0.0))
    , m_fCoeffGamma52(F(0.0))
    , m_fCoeffGamma53(F(0.0))
    , m_fCoeffGamma54(F(0.0))
    , m_fCoeffSigma12(F(0.0))
    , m_fCoeffSigma13(F(0.0))
    , m_fCoeffSigma14(F(0.0))
    , m_fCoeffSigma23(F(0.0))
    , m_fCoeffSigma24(F(0.0))
    , m_fCoeffSigma34(F(0.0))
{
    
}

CFieldFermionWilsonSquareSU3Gamma::~CFieldFermionWilsonSquareSU3Gamma()
{
    
}

void CFieldFermionWilsonSquareSU3Gamma::InitialOtherParameters(CParameters & params)
{
    CFieldFermionWilsonSquareSU3D::InitialOtherParameters(params);

    INT iValue = 1;
    if (params.FetchValueINT(_T("ExpGamma"), iValue))
    {
        m_bExpGamma = (0 != iValue);
    }

    Real fValue = F(0.0);
    if (params.FetchValueReal(_T("Gamma1"), fValue))
    {
        m_fCoeffGamma1 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma2"), fValue))
    {
        m_fCoeffGamma2 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma3"), fValue))
    {
        m_fCoeffGamma3 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma4"), fValue))
    {
        m_fCoeffGamma4 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma5"), fValue))
    {
        m_fCoeffGamma5 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma51"), fValue))
    {
        m_fCoeffGamma51 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma52"), fValue))
    {
        m_fCoeffGamma52 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma53"), fValue))
    {
        m_fCoeffGamma53 = fValue;
    }
    if (params.FetchValueReal(_T("Gamma54"), fValue))
    {
        m_fCoeffGamma54 = fValue;
    }
    if (params.FetchValueReal(_T("Sigma12"), fValue))
    {
        m_fCoeffSigma12 = fValue;
    }
    if (params.FetchValueReal(_T("Sigma13"), fValue))
    {
        m_fCoeffSigma13 = fValue;
    }
    if (params.FetchValueReal(_T("Sigma14"), fValue))
    {
        m_fCoeffSigma14 = fValue;
    }
    if (params.FetchValueReal(_T("Sigma23"), fValue))
    {
        m_fCoeffSigma23 = fValue;
    }
    if (params.FetchValueReal(_T("Sigma24"), fValue))
    {
        m_fCoeffSigma24 = fValue;
    }
    if (params.FetchValueReal(_T("Sigma34"), fValue))
    {
        m_fCoeffSigma34 = fValue;
    }
}

void CFieldFermionWilsonSquareSU3Gamma::CopyTo(CField* U) const
{
    CFieldFermionWilsonSquareSU3D::CopyTo(U);
    CFieldFermionWilsonSquareSU3Gamma* pOther = dynamic_cast<CFieldFermionWilsonSquareSU3Gamma*>(U);

    pOther->m_bExpGamma = m_bExpGamma;
    pOther->m_fCoeffGamma1 = m_fCoeffGamma1;
    pOther->m_fCoeffGamma2 = m_fCoeffGamma2;
    pOther->m_fCoeffGamma3 = m_fCoeffGamma3;
    pOther->m_fCoeffGamma4 = m_fCoeffGamma4;
    pOther->m_fCoeffGamma5 = m_fCoeffGamma5;
    pOther->m_fCoeffGamma51 = m_fCoeffGamma51;
    pOther->m_fCoeffGamma52 = m_fCoeffGamma52;
    pOther->m_fCoeffGamma53 = m_fCoeffGamma53;
    pOther->m_fCoeffGamma54 = m_fCoeffGamma54;
    pOther->m_fCoeffSigma12 = m_fCoeffSigma12;
    pOther->m_fCoeffSigma13 = m_fCoeffSigma13;
    pOther->m_fCoeffSigma14 = m_fCoeffSigma14;
    pOther->m_fCoeffSigma23 = m_fCoeffSigma23;
    pOther->m_fCoeffSigma24 = m_fCoeffSigma24;
    pOther->m_fCoeffSigma34 = m_fCoeffSigma34;
}

CCString CFieldFermionWilsonSquareSU3Gamma::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldFermionWilsonSquareSU3::GetInfos(tab);

    sRet = sRet + tab + _T("ExpGamma : ") + (m_bExpGamma ? _T("1") : _T("0")) + _T("\n");

    sRet = sRet + tab + _T("Gamma1 : ") + appToString(m_fCoeffGamma1) + _T("\n");
    sRet = sRet + tab + _T("Gamma2 : ") + appToString(m_fCoeffGamma2) + _T("\n");
    sRet = sRet + tab + _T("Gamma3 : ") + appToString(m_fCoeffGamma3) + _T("\n");
    sRet = sRet + tab + _T("Gamma4 : ") + appToString(m_fCoeffGamma4) + _T("\n");
    sRet = sRet + tab + _T("Gamma5 : ") + appToString(m_fCoeffGamma5) + _T("\n");

    sRet = sRet + tab + _T("Gamma51 : ") + appToString(m_fCoeffGamma51) + _T("\n");
    sRet = sRet + tab + _T("Gamma52 : ") + appToString(m_fCoeffGamma52) + _T("\n");
    sRet = sRet + tab + _T("Gamma53 : ") + appToString(m_fCoeffGamma53) + _T("\n");
    sRet = sRet + tab + _T("Gamma54 : ") + appToString(m_fCoeffGamma54) + _T("\n");

    sRet = sRet + tab + _T("Sigma12 : ") + appToString(m_fCoeffSigma12) + _T("\n");
    sRet = sRet + tab + _T("Sigma13 : ") + appToString(m_fCoeffSigma13) + _T("\n");
    sRet = sRet + tab + _T("Sigma14 : ") + appToString(m_fCoeffSigma14) + _T("\n");
    sRet = sRet + tab + _T("Sigma23 : ") + appToString(m_fCoeffSigma23) + _T("\n");
    sRet = sRet + tab + _T("Sigma24 : ") + appToString(m_fCoeffSigma24) + _T("\n");
    sRet = sRet + tab + _T("Sigma34 : ") + appToString(m_fCoeffSigma34) + _T("\n");

    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================