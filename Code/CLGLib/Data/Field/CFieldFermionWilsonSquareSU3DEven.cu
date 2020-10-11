//=============================================================================
// FILENAME : CFieldFermionWilsonSU3DEven.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/18/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSU3DEven)

#pragma region DOperator

#pragma region kernel

/**
* Step1 is phi_odd = Doo^{-1} Doe phi_even,
* in phi_even = Dee phi_even - Deo (Doo^{-1} Doe phi_even)
* Note: for M = c D, it neglect all the coefficient in Doo^{-1} Doe.
* Also, because gamma5 ^2 = 1, for M = Ddagger = gamma5 D gamma5,
* then it is: gamma5 (Deo Doo^{-1} Doe) gamma5, so we keep only the gamma5 on the right
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_D_Even_Step1(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    intokernalInt4_odd;
    BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    //
    if (sIdx.IsDirichlet() || !sSite4.IsOdd())
    {
        //also set Dirichlet to zero because in Doe, the matrix element is zero
        pResultData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    //Here, the right vector are all even sites
    const gammaMatrix& gamma5 = __chiralGamma[GAMMA5];
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        gammaMatrix gammaMu = __chiralGamma[GAMMA1 + idir];

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

        SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        //deviceSU3 x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir);
        //deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

        //deviceWilsonVectorSU3 x_p_mu_Fermion_element = pDeviceData[x_p_mu_Fermion.m_uiSiteIndex];
        //deviceWilsonVectorSU3 x_m_mu_Fermion_element = pDeviceData[x_m_mu_Fermion.m_uiSiteIndex];
        deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
        deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);

        if (bDDagger)
        {
            x_p_mu_Fermion_element = gamma5.MulWilsonC(x_p_mu_Fermion_element);
            x_m_mu_Fermion_element = gamma5.MulWilsonC(x_m_mu_Fermion_element);
        }

        //hopping terms

        //U(x,mu) phi(x+ mu)
        deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            //printf("Opposite x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
            result.Sub(u_phi_x_p_m);

            //- gammamu U(x,mu) phi(x+ mu)
            result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        }
        else
        {
            result.Add(u_phi_x_p_m);

            //- gammamu U(x,mu) phi(x+ mu)
            result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            result.Sub(u_dagger_phi_x_m_m);

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
        else
        {
            result.Add(u_dagger_phi_x_m_m);

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
    }

    //result = phi(x) - kai sum _mu result
    //result.MulReal(kai);
    //pResultData[uiSiteIndex].Sub(result);
    result.MulReal(-kai);
    pResultData[uiSiteIndex] = result;
}

/**
 * phi_even = Dee phi_even - Deo phi_odd
 * If we have coefficient, then, it is c (Dee phi_even - Deo phi_odd)
 * If we have dagger, then, it is c (Dee phi_even - gamma5 Deo phi_odd)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_D_Even_Step2(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4_even;
    BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        gammaMatrix gammaMu = __chiralGamma[GAMMA1 + idir];

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

        SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        //deviceSU3 x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir);
        //deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

        //deviceWilsonVectorSU3 x_p_mu_Fermion_element = pResultData[x_p_mu_Fermion.m_uiSiteIndex];
        //deviceWilsonVectorSU3 x_m_mu_Fermion_element = pResultData[x_m_mu_Fermion.m_uiSiteIndex];
        deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pResultData, x_p_mu_Fermion, byFieldId);
        deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pResultData, x_m_mu_Fermion, byFieldId);

        //hopping terms

        //U(x,mu) phi(x+ mu)
        deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            //printf("Opposite x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
            result.Sub(u_phi_x_p_m);

            //- gammamu U(x,mu) phi(x+ mu)
            result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        }
        else
        {
            result.Add(u_phi_x_p_m);

            //- gammamu U(x,mu) phi(x+ mu)
            result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            result.Sub(u_dagger_phi_x_m_m);

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
        else
        {
            result.Add(u_dagger_phi_x_m_m);

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
    }

    //result = phi(x) - kai sum _mu result
    result.MulReal(kai);
    //pResultData[uiSiteIndex].Sub(result);
    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }
    pResultData[uiSiteIndex].Add(result);

    switch (eCoeff)
    {
    case EOCT_Real:
        pResultData[uiSiteIndex].MulReal(fCoeff);
        break;
    case EOCT_Complex:
        pResultData[uiSiteIndex].MulComp(cCoeff);
        break;
    }
}

#pragma endregion

void CFieldFermionWilsonSU3DEven::DOperator(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceWilsonVectorSU3* pTarget = (deviceWilsonVectorSU3*)pTargetBuffer;
    const deviceWilsonVectorSU3* pSource = (deviceWilsonVectorSU3*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    
    preparethread_even;

    _kernelDFermionWilsonSquareSU3_D_Even_Step1 << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pTarget,
        m_fKai,
        m_byParentId,
        bDagger);

    _kernelDFermionWilsonSquareSU3_D_Even_Step2 << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pTarget,
        m_fKai,
        m_byParentId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    SetOddZero(pTarget);
}

void CFieldFermionWilsonSU3DEven::DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const
{
    appCrucial(_T("You shall never enter this function!"));
    _FAIL_EXIT;
}

#pragma endregion

#pragma region Kernel

/**
 * - Deo Doo^{-1} psi_odd
 * For Wilson-Dirac, Doo=1, so it is just
 * - Deo psi_odd
 *
 * For Ddagger
 *  - gamma5 Deo Doo^{-1} gamma5 psi_odd
 *
 * D is "- kai sum _mu (1-gamma _mu) U(x,mu) phi(x+ mu) + (1+gamma _mu) U^{dagger}(x-mu) phi(x-mu)"
 * This works when left is e, and right is o
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelCopyWilsonSquareSU3DEven_Step1(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData, /*Parent, right vector*/
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    BYTE byParentFieldId,
    UBOOL bDdagger)
{
    intokernalInt4_even;

    //make sure they are all even

    BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byParentFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    //To here, the left vector only has even sites
    const gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        gammaMatrix gammaMu = __chiralGamma[GAMMA1 + idir];

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

        SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        //deviceSU3 x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir);
        //deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

        //deviceWilsonVectorSU3 x_p_mu_Fermion_element = pDeviceData[x_p_mu_Fermion.m_uiSiteIndex];
        //deviceWilsonVectorSU3 x_m_mu_Fermion_element = pDeviceData[x_m_mu_Fermion.m_uiSiteIndex];
        deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byParentFieldId);
        deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byParentFieldId);
        if (bDdagger)
        {
            x_p_mu_Fermion_element = gamma5.MulWilsonC(x_p_mu_Fermion_element);
            x_m_mu_Fermion_element = gamma5.MulWilsonC(x_m_mu_Fermion_element);
        }
        //hopping terms

        //U(x,mu) phi(x+ mu)
        deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            //printf("Opposite x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
            result.Sub(u_phi_x_p_m);

            //- gammamu U(x,mu) phi(x+ mu)
            result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        }
        else
        {
            result.Add(u_phi_x_p_m);

            //- gammamu U(x,mu) phi(x+ mu)
            result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            result.Sub(u_dagger_phi_x_m_m);

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
        else
        {
            result.Add(u_dagger_phi_x_m_m);

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
    }

    //result = phi(x) - kai sum _mu result
    result.MulReal(kai);
    // pResultData[uiSiteIndex].Sub(result); <-- this is the orignal one, now we only need -Deo
    if (bDdagger)
    {
        result = gamma5.MulWilsonC(result);
    }

    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex].AddC(result);
}


/**
 * z_even obtained, we calculate D_oo^{-1}(phi_odd - D_oe z_even)
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_EOFinalStep(
    const deviceWilsonVectorSU3* __restrict__ pZEven,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    BYTE byFieldId,
    UBOOL bDDagger)
{
    UINT uiEvenIndex = (((threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z)) << 1);
    UINT uiSiteIndex = uiEvenIndex + 1; //odd index
    SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);
    if (!sSite4.IsOdd())
    {
        uiSiteIndex = uiSiteIndex - 1; //odd index
        uiEvenIndex = uiEvenIndex + 1;
        sSite4 = __deviceSiteIndexToInt4(uiSiteIndex); //odd index
    }

    //write z_even
    pResultData[uiEvenIndex] = pZEven[uiEvenIndex];

    //====================
    //Now we only care the odd sites
    BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        pResultData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    gammaMatrix gamma5 = __chiralGamma[GAMMA5];
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        gammaMatrix gammaMu = __chiralGamma[GAMMA1 + idir];

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        SIndex x_m_mu_Gauge = pGaugeMove[linkIndex];

        SIndex x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        SIndex x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        //deviceSU3 x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir);
        //deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

        //deviceWilsonVectorSU3 x_p_mu_Fermion_element = pDeviceData[x_p_mu_Fermion.m_uiSiteIndex];
        //deviceWilsonVectorSU3 x_m_mu_Fermion_element = pDeviceData[x_m_mu_Fermion.m_uiSiteIndex];
        deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pZEven, x_p_mu_Fermion, byFieldId);
        deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pZEven, x_m_mu_Fermion, byFieldId);

        if (bDDagger)
        {
            x_p_mu_Fermion_element = gamma5.MulWilsonC(x_p_mu_Fermion_element);
            x_m_mu_Fermion_element = gamma5.MulWilsonC(x_m_mu_Fermion_element);
        }

        //hopping terms

        //U(x,mu) phi(x+ mu)
        deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            //printf("Opposite x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
            result.Sub(u_phi_x_p_m);

            //- gammamu U(x,mu) phi(x+ mu)
            result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        }
        else
        {
            result.Add(u_phi_x_p_m);

            //- gammamu U(x,mu) phi(x+ mu)
            result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            result.Sub(u_dagger_phi_x_m_m);

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
        else
        {
            result.Add(u_dagger_phi_x_m_m);

            //gammamu U^{dagger}(x-mu) phi(x-mu)
            result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
    }

    //result = phi(x) - kai sum _mu result
    result.MulReal(kai);
    //result = -Deo z_even
    //       = -Deo gamma5 z_even

    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }
    //result = -Deo z_even
    //       = -gamma5 Deo gamma5 z_even
    
    //D_oo^ { -1 }(phi_odd - D_oe z_even)
    pResultData[uiSiteIndex].Add(result);

}


__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusFermionWilsonSquareSU3_Even(
    deviceWilsonVectorSU3* pMe,
    const deviceWilsonVectorSU3* __restrict__ pOther)
{
    intokernal_even;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusFermionWilsonSquareSU3_Even(
    deviceWilsonVectorSU3* pMe,
    const deviceWilsonVectorSU3* __restrict__ pOther)
{
    intokernal_even;
    pMe[uiSiteIndex].Sub(pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyComplexFermionWilsonSquareSU3_Even(
    deviceWilsonVectorSU3* pMe,
    const deviceWilsonVectorSU3* __restrict__ pOther, CLGComplex a)
{
    intokernal_even;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex].MulCompC(a));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyRealFermionWilsonSquareSU3_Even(
    deviceWilsonVectorSU3* pMe,
    const deviceWilsonVectorSU3* __restrict__ pOther, Real a)
{
    intokernal_even;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex].MulRealC(a));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDotFermionWilsonSquareSU3_Even(
    const deviceWilsonVectorSU3* __restrict__ pMe,
    const deviceWilsonVectorSU3* __restrict__ pOther,
#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex* result
#else
    CLGComplex* result
#endif
)
{
    const UINT uiEvenIndex = ((threadIdx.x + blockIdx.x * blockDim.x) * _DC_GridDimZT + (threadIdx.y + blockIdx.y * blockDim.y) * _DC_Lt + (threadIdx.z + blockIdx.z * blockDim.z));
    UINT uiSiteIndex = uiEvenIndex << 1; 
    if (__deviceSiteIndexToInt4(uiSiteIndex).IsOdd())
    { 
        uiSiteIndex = uiSiteIndex + 1; 
    }
    //SSmallInt4 sSite4 = __deviceSiteIndexToInt4(uiSiteIndex);
    //printf("uiEvenIndex = %d, idx = (%d, %d, %d, %d)\n", uiEvenIndex, sSite4.x, sSite4.y, sSite4.z, sSite4.w);
#if !_CLG_DOUBLEFLOAT
    result[uiEvenIndex] = _cToDouble(pMe[uiSiteIndex].ConjugateDotC(pOther[uiSiteIndex]));
#else
    result[uiEvenIndex] = pMe[uiSiteIndex].ConjugateDotC(pOther[uiSiteIndex]);
#endif
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyComplex_Even(
    deviceWilsonVectorSU3* pMe,
    CLGComplex a)
{
    intokernal_even;
    pMe[uiSiteIndex].MulComp(a);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyReal_Even(
    deviceWilsonVectorSU3* pMe,
    Real a)
{
    intokernal_even;
    pMe[uiSiteIndex].MulReal(a);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFermionWilsonSquareSU3Conjugate_Even(deviceWilsonVectorSU3* pDevicePtr)
{
    intokernal_even;
    pDevicePtr[uiSiteIndex].Conjugate();
}


__global__ void _CLG_LAUNCH_BOUND
_kernelFermionWilsonMakeOddZero(deviceWilsonVectorSU3* pDevicePtr)
{
    intokernal_odd;
    pDevicePtr[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
}

#pragma endregion


void CFieldFermionWilsonSU3DEven::CopyTo(CField* U) const
{
    CFieldFermionWilsonSquareSU3D::CopyTo(U);
    CFieldFermionWilsonSU3DEven* pField = dynamic_cast<CFieldFermionWilsonSU3DEven*>(U);
    if (NULL != pField)
    {
        pField->m_byParentId = m_byParentId;
    }
}

/**
 * phi_even = psi_even - Deo Doo^{-1} psi_odd
 */
void CFieldFermionWilsonSU3DEven::WriteEvenSites(const CFieldFermion* x, const CFieldGauge* pGauge, UBOOL bDdagger)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3* pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    m_byParentId = pField->m_byFieldId;

    preparethread_even;
    _kernelCopyWilsonSquareSU3DEven_Step1 << <block, threads >> > (
        pField->m_pDeviceData,
        pFieldSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        m_pDeviceData,
        CCommonData::m_fKai,
        m_byParentId,
        bDdagger);

    SetOddZero(m_pDeviceData);

    if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
    {
        m_fLength = Dot(this).x;
    }

    if (0 != (_HC_Lx & 1) || 0 != (_HC_Ly & 1) || 0 != (_HC_Lt & 1) || 0 != (_HC_Lt & 1))
    {
        appGeneral(_T("Warning: Even-Odd work only for even extent!\n"));
    }

    //DebugPrintMe();
}

void CFieldFermionWilsonSU3DEven::WriteBackEvenSites(CFieldFermion* pParentField, const CFieldGauge* pGauge, UBOOL bDdagger) const
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    if (NULL == pParentField || EFT_FermionWilsonSquareSU3 != pParentField->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3DEven can only play with gauge CFieldFermionWilsonSquareSU3!"));
        return;
    }
    CFieldFermionWilsonSquareSU3* pParentFieldSU3 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pParentField);

    preparethread_even;
    _kernelDFermionWilsonSquareSU3_EOFinalStep << <block, threads >> > (
        m_pDeviceData,
        pFieldSU3->m_pDeviceData,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byParentId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byParentId],
        pParentFieldSU3->m_pDeviceData,
        CCommonData::m_fKai,
        m_byParentId,
        bDdagger);
}

void CFieldFermionWilsonSU3DEven::AxpyPlus(const CField* x)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3* pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethread_even;
    _kernelAxpyPlusFermionWilsonSquareSU3_Even << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldFermionWilsonSU3DEven::AxpyMinus(const CField* x)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3* pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethread_even;
    _kernelAxpyMinusFermionWilsonSquareSU3_Even << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldFermionWilsonSU3DEven::Axpy(Real a, const CField* x)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3* pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethread_even;
    _kernelAxpyRealFermionWilsonSquareSU3_Even << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

void CFieldFermionWilsonSU3DEven::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3* pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethread_even;
    _kernelAxpyComplexFermionWilsonSquareSU3_Even << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

#if !_CLG_DOUBLEFLOAT
cuDoubleComplex CFieldFermionWilsonSU3DEven::Dot(const CField* x) const
#else
CLGComplex CFieldFermionWilsonSU3DEven::Dot(const CField* x) const
#endif
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return make_cuDoubleComplex(0, 0);
    }
    const CFieldFermionWilsonSquareSU3* pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethread_even;
    _kernelDotFermionWilsonSquareSU3_Even << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, _D_ComplexThreadBuffer);
    return appGetCudaHelper()->ReduceComplex(_D_ComplexThreadBuffer, _HC_Volume >> 1);
}

void CFieldFermionWilsonSU3DEven::ScalarMultply(const CLGComplex& a)
{
    preparethread_even;
    _kernelScalarMultiplyComplex_Even << <block, threads >> > (m_pDeviceData, a);
}

void CFieldFermionWilsonSU3DEven::ScalarMultply(Real a)
{
    preparethread_even;
    _kernelScalarMultiplyReal_Even << <block, threads >> > (m_pDeviceData, a);
}


void CFieldFermionWilsonSU3DEven::Dagger()
{
    preparethread_even;
    _kernelFermionWilsonSquareSU3Conjugate_Even << <block, threads >> > (m_pDeviceData);
}


void CFieldFermionWilsonSU3DEven::SetOddZero(deviceWilsonVectorSU3* pDeviceData)
{
    preparethread_even;
    _kernelFermionWilsonMakeOddZero << <block, threads >> > (pDeviceData);
}

CCString CFieldFermionWilsonSU3DEven::GetInfos(const CCString& tab) const
{
    CCString sRet = tab + _T("Name : CFieldFermionWilsonSU3DEven\n");
    sRet = sRet + tab + _T("Hopping : ") + appFloatToString(CCommonData::m_fKai) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================