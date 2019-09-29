//=============================================================================
// FILENAME : CFieldFermionWilsonSU3EvenOdd.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/18/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"


__BEGIN_NAMESPACE

#if FURTURE

__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSU3EvenOdd)

#pragma region DOperator

#pragma region kernel

/**
* Dw phi(x) = phi(x) - kai sum _mu (1-gamma _mu) U(x,mu) phi(x+ mu) + (1+gamma _mu) U^{dagger}(x-mu) phi(x-mu)
* U act on su3
* gamma act on spinor
*
* byEvenOdd -> 0 : all, 1 : even to odd, 2 : odd to even
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3_D_EO(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    BYTE byFieldId,
    EIndexEvenOdd eEvenOdd)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    if (sIdx.IsDirichlet())
    {
        return;
    }

    if (sSite4.IsOdd() && eEvenOdd == EIE_OddToEven)
    {
        return;
    }

    if (!sSite4.IsOdd() && eEvenOdd == EIE_EvenToOdd)
    {
        return;
    }

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
        x_m_mu_Gauge_element.Dagger();

        //deviceWilsonVectorSU3 x_p_mu_Fermion_element = pDeviceData[x_p_mu_Fermion.m_uiSiteIndex];
        //deviceWilsonVectorSU3 x_m_mu_Fermion_element = pDeviceData[x_m_mu_Fermion.m_uiSiteIndex];
        deviceWilsonVectorSU3 x_p_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_p_mu_Fermion, byFieldId);
        deviceWilsonVectorSU3 x_m_mu_Fermion_element = _deviceGetFermionBCWilsonSU3(pDeviceData, x_m_mu_Fermion, byFieldId);

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
    result.MulReal(-kai);
    pResultData[uiSiteIndex] = result;
}

/**
* The output is on a gauge field
* Therefor cannot make together with _kernelDWilson
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDWilsonForceSU3_D(
    const deviceWilsonVectorSU3* __restrict__ pInverseD,
    const deviceWilsonVectorSU3* __restrict__ pInverseDDdagger,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    deviceSU3* pForce,
    Real fKai,
    BYTE byFieldId)
{
    intokernalInt4;
    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex sSite = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];
    const deviceWilsonVectorSU3 x_Left(_deviceGetFermionBCWilsonSU3(pInverseDDdagger, sSite, byFieldId));
    const deviceWilsonVectorSU3 x_Right(_deviceGetFermionBCWilsonSU3(pInverseD, sSite, byFieldId));

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        SIndex x_p_mu_Fermion = pFermionMove[linkIndex * 2];

        //If one of the sites is on surface, it has no contribution.
        //Note that, the bond on surface is equivelant to both sites on surface.

        if (//__idx->_deviceIsBondOnSurface(uiBigIdx, idir) ||
            x_p_mu_Fermion.IsDirichlet()
            || sSite.IsDirichlet())
        {
            //continue;
        }
        else
        {
            //Get Gamma mu
            gammaMatrix gammaMu = __chiralGamma[GAMMA1 + idir];

            //SIndex x_m_mu_Gauge = __idx->_deviceGaugeIndexWalk(uiSiteIndex, -(idir + 1));
             // __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, (idir + 1));

            //all not on surface
            const deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
            deviceWilsonVectorSU3 x_p_mu_Left(pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex]);
            //deviceWilsonVectorSU3 x_p_mu_Right = _deviceGetFermionBCWilsonSU3(pInverseD, x_p_mu_Fermion, byFieldId);
            //deviceWilsonVectorSU3 x_p_mu_Left = _deviceGetFermionBCWilsonSU3(pInverseDDdagger, x_p_mu_Fermion, byFieldId);

            deviceSU3 x_Gauge_element = pGauge[linkIndex]; // _deviceGetGaugeBCSU3Dir(pGauge, uiBigIdx, idir); //pGauge[linkIndex];

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
                forceOfThisLink.MulReal(-fKai);
            }
            else
            {
                forceOfThisLink.MulReal(fKai);
            }

            pForce[linkIndex].Add(forceOfThisLink);
        }
    }
}

#pragma endregion

void CFieldFermionWilsonSU3EvenOdd::DOperator(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceWilsonVectorSU3* pTarget = (deviceWilsonVectorSU3*)pTargetBuffer;
    const deviceWilsonVectorSU3* pSource = (deviceWilsonVectorSU3*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    preparethread;
    _kernelDFermionWilsonSquareSU3_D << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        pTarget,
        m_fKai,
        m_byFieldId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

}

void CFieldFermionWilsonSU3EvenOdd::DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const
{
    deviceSU3* pForceSU3 = (deviceSU3*)pForce;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
    const deviceWilsonVectorSU3* pDphiBuffer = (deviceWilsonVectorSU3*)pDphi;
    const deviceWilsonVectorSU3* pDDphiBuffer = (deviceWilsonVectorSU3*)pDDphi;

    preparethread;
    _kernelDWilsonForceSU3_D << <block, threads >> > (
        pDphiBuffer,
        pDDphiBuffer,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        pForceSU3,
        m_fKai, m_byFieldId);
}

#pragma endregion

#pragma region Kernel



#pragma endregion


void CFieldFermionWilsonSU3EvenOdd::CopyTo(CField* U) const
{
    CFieldFermionWilsonSquareSU3D::CopyTo(U);
}


CCString CFieldFermionWilsonSU3EvenOdd::GetInfos(const CCString& tab) const
{
    CCString sRet = tab + _T("Name : CFieldFermionWilsonSU3EvenOdd\n");
    sRet = sRet + tab + _T("Hopping : ") + appFloatToString(CCommonData::m_fKai) + _T("\n");
    return sRet;
}

#endif

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================