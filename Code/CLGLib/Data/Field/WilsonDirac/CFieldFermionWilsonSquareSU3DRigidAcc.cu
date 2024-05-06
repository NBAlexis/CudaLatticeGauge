//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3DRigidAcc.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [08/01/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSquareSU3DRigidAcc)

#pragma region DOperator

#pragma region kernel

#pragma region hamitonian

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DRigAcc_XYZ(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    Real fG,
    BYTE byFieldId,
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
        pResultData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
        return;
    }

    const gammaMatrix& gamma5 = __chiralGamma[GAMMA5];
    const Real fOnePlusGZ = F(1.0) + sSite4.z * fG;

    /* gamma5 ^2 = 1
    pResultData[uiSiteIndex] = bDDagger ? gamma5.MulWilsonC(pDeviceData[uiSiteIndex]) : pDeviceData[uiSiteIndex];
    pResultData[uiSiteIndex].MulReal(fOnePlusGZ - F(2.0) * fG * sSite4.z * kai);

    if (bDDagger)
    {
        pResultData[uiSiteIndex] = gamma5.MulWilsonC(pResultData[uiSiteIndex]);
    }
    */
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];
    pResultData[uiSiteIndex].MulReal(fOnePlusGZ - F(2.0) * fG * sSite4.z * kai);

    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();

    //idir = mu
    #pragma unroll
    for (UINT idir = 0; idir < 3; ++idir)
    {
        //=========================
        //get things
        const gammaMatrix& gammaMu = __chiralGamma[GAMMA1 + idir];
        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];
        SSmallInt4 site_m_mu = __deviceSiteIndexToInt4(x_m_mu_Fermion.m_uiSiteIndex);
        const Real fOnePlusGZm1 = F(1.0) + site_m_mu.z * fG;

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        //deviceSU3 x_Gauge_element = pGauge[linkIndex];
        const deviceSU3& x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir);
        deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
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

        //U(x,mu) phi(x+ mu)
        deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
        u_phi_x_p_m.MulReal(fOnePlusGZ);

        //(1 - gammamu) U(x,mu) phi(x+ mu)
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            result.Sub(u_phi_x_p_m);
            result.Add(gammaMu.MulWilsonC(u_phi_x_p_m));
        }
        else
        {
            result.Add(u_phi_x_p_m);
            result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
        u_dagger_phi_x_m_m.MulReal(fOnePlusGZm1);

        // (1 + gammamu) U^{dagger}(x-mu) phi(x-mu)
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            result.Sub(u_dagger_phi_x_m_m);
            result.Sub(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
        else
        {
            result.Add(u_dagger_phi_x_m_m);
            result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
        }
    }

    //result.MulReal(fOnePlusGZ * kai);
    result.MulReal(kai);
    //result = phi(x) - kai sum _mu result
    //result.MulReal(kai);
    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }

    //==============================
    //res = [gamma5 (orig - term3 - result) gamma5] * coeff
    //    = [gamma5 orig gamma5] * coeff - [gamma5 (term3 + result) gamma5] * coeff
    //    = res0 - [gamma5 (term3 + result) gamma5] * coeff
    pResultData[uiSiteIndex].Sub(result);
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

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_DRigAcc_T(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    Real fG2,
    BYTE byFieldId,
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
    const gammaMatrix& sigma43 = __chiralGamma[SIGMA43];

    deviceWilsonVectorSU3 result = pDeviceData[uiSiteIndex];
    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }
    result = sigma43.MulWilsonC(result);
    result = gamma4.MulWilsonC(result);
    result.MulReal(-fG2);

    //idir = mu, we have x, y and t term but no z term
    //=========================
    //get things

    //x, mu
    UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);

    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    //Assuming periodic
    //get U(x,mu), U^{dagger}(x-mu), 
    //deviceSU3 x_Gauge_element = pGauge[linkIndex];
    const deviceSU3& x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, 3);
    deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3(pGauge, x_m_mu_Gauge);
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

    //U(x,mu) phi(x+ mu)
    deviceWilsonVectorSU3 u_phi_x_p_m = x_Gauge_element.MulWilsonVector(x_p_mu_Fermion_element);
    u_phi_x_p_m.Sub(gamma4.MulWilsonC(u_phi_x_p_m));

    if (x_p_mu_Fermion.NeedToOpposite())
    {
        //printf("OppositeT x=%d y=%d z=%d t=%d\n", static_cast<INT>(sSite4.x), static_cast<INT>(sSite4.y), static_cast<INT>(sSite4.z), static_cast<INT>(sSite4.w));
        //k(cos_m_1 - i sin sigma34)(1-gamma4)
        result.Sub(u_phi_x_p_m);
    }
    else
    {
        //-k(cos_m_1 - i sin sigma34)(1-gamma4)
        result.Add(u_phi_x_p_m);
    }

    //U^{dagger}(x-mu) phi(x-mu)
    deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
    u_dagger_phi_x_m_m.Add(gamma4.MulWilsonC(u_dagger_phi_x_m_m));

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        //-k(cos_m_1+isin)(1+gamma4)
        result.Sub(u_dagger_phi_x_m_m);
    }
    else
    {
        result.Add(u_dagger_phi_x_m_m);
    }

    //result = phi(x) - kai sum _mu result
    //result.MulReal(kai);
    if (bDDagger)
    {
        result = gamma5.MulWilsonC(result);
    }
    result.MulReal(kai);
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
    pResultData[uiSiteIndex].Sub(result);
}

#pragma endregion

#pragma region force

__global__ void _CLG_LAUNCH_BOUND
_kernelDWilsonForceSU3_DRigAcc_XYZ(
    const deviceWilsonVectorSU3* __restrict__ pInverseD,
    const deviceWilsonVectorSU3* __restrict__ pInverseDDdagger,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    deviceSU3* pForce,
    Real fKai,
    Real fG,
    BYTE byFieldId)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex& sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (sSite.IsDirichlet())
    {
        return;
    }

    const deviceWilsonVectorSU3 x_Left(_deviceGetFermionBCWilsonSU3(pInverseDDdagger, sSite, byFieldId));
    const deviceWilsonVectorSU3 x_Right(_deviceGetFermionBCWilsonSU3(pInverseD, sSite, byFieldId));
    Real fOnePlusGZ = F(1.0) + sSite4.z * fG;

    #pragma unroll
    for (UINT idir = 0; idir < 4; ++idir)
    {
        if (3 == idir)
        {
            fOnePlusGZ = F(1.0);
        }

        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        const SIndex x_p_mu_Fermion = pFermionMove[linkIndex * 2];
        if (x_p_mu_Fermion.IsDirichlet())
        {
            continue;
        }
        //const SSmallInt4 site_p_mu = __deviceSiteIndexToInt4(x_p_mu_Fermion.m_uiSiteIndex);
        //const Real fOnePlusGZ1 = F(1.0);// +site_p_mu.z * fG;

        //Get Gamma mu
        const gammaMatrix& gammaMu = __chiralGamma[GAMMA1 + idir];

        //all not on surface
        const deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
        const deviceWilsonVectorSU3 x_p_mu_Left(pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex]);

        const deviceSU3& x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

        deviceWilsonVectorSU3 right1(x_p_mu_Right);
        right1.Sub(gammaMu.MulWilsonC(right1));
        deviceSU3 mid = deviceSU3::makeSU3Contract(x_Left, right1);

        deviceWilsonVectorSU3 right2(x_Right);
        right2.Add(gammaMu.MulWilsonC(right2));

        mid.Add(deviceSU3::makeSU3Contract(right2, x_p_mu_Left));

        deviceSU3 forceOfThisLink = x_Gauge_element.MulC(mid);
        forceOfThisLink.Ta();
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            forceOfThisLink.MulReal(-fKai * fOnePlusGZ);
        }
        else
        {
            forceOfThisLink.MulReal(fKai * fOnePlusGZ);
        }

        pForce[linkIndex].Add(forceOfThisLink);
    }
}

#pragma endregion

#pragma endregion

/**
* It will not produce correct result unless _CLG_USE_LAUNCH_BOUND is set to 1
* It seems max-registor count is reached without throw an error
* To avoid it, we use slower method, split it into two functions.
*/
void CFieldFermionWilsonSquareSU3DRigidAcc::DOperator(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, 
    UBOOL bDagger, EOperatorCoefficientType eOCT, 
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceWilsonVectorSU3* pTarget = (deviceWilsonVectorSU3*)pTargetBuffer;
    const deviceWilsonVectorSU3* pSource = (deviceWilsonVectorSU3*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
    preparethread;
    _kernelDFermionWilsonSquareSU3_DRigAcc_XYZ << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        pTarget,
        m_fKai,
        CCommonData::m_fG,
        m_byFieldId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    _kernelDFermionWilsonSquareSU3_DRigAcc_T << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        pTarget,
        m_fKai,
        CCommonData::m_fG,
        m_byFieldId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

}

void CFieldFermionWilsonSquareSU3DRigidAcc::DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const
{
    deviceSU3* pForceSU3 = (deviceSU3*)pForce;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
    const deviceWilsonVectorSU3* pDphiBuffer = (deviceWilsonVectorSU3*)pDphi;
    const deviceWilsonVectorSU3* pDDphiBuffer = (deviceWilsonVectorSU3*)pDDphi;

    preparethread;

    _kernelDWilsonForceSU3_DRigAcc_XYZ << <block, threads >> > (
        pDphiBuffer,
        pDDphiBuffer,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        pForceSU3,
        m_fKai,
        CCommonData::m_fG,
        m_byFieldId);
}

#pragma endregion

void CFieldFermionWilsonSquareSU3DRigidAcc::InitialOtherParameters(CParameters& params)
{
    CFieldFermionWilsonSquareSU3D::InitialOtherParameters(params);

    Real fG2 = F(0.2);
    params.FetchValueReal(_T("G2"), fG2);
    m_fG2 = fG2;
}


void CFieldFermionWilsonSquareSU3DRigidAcc::CopyTo(CField* U) const
{
    CFieldFermionWilsonSquareSU3D::CopyTo(U);
}

CCString CFieldFermionWilsonSquareSU3DRigidAcc::GetInfos(const CCString &tab) const
{
    CCString sRet = tab + _T("Name : CFieldFermionWilsonSquareSU3DRigidAcc\n");
    sRet = sRet + CFieldFermion::GetInfos(tab);
    sRet = sRet + tab + _T("Hopping : ") + appToString(CCommonData::m_fKai) + _T("\n");
    sRet = sRet + tab + _T("G2 : ") + appToString(m_fG2) + _T("\n");

    SSmallInt4 boundary = appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(m_byFieldId);
    sRet = sRet + tab + _T("boundary : [") + appToString(boundary.x) + _T(", ") + appToString(boundary.y) + _T(", ") + appToString(boundary.z) + _T(", ") + appToString(boundary.w) + _T("]\n");

    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================