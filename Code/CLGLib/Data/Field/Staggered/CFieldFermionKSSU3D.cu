//=============================================================================
// FILENAME : CFieldFermionKSSU3D.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [09/19/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionKSSU3D)

#pragma region DOperator

#pragma region kernel

/**
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_D(
    const deviceSU3Vector* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceSU3Vector* pResultData,
    Real f2am,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    const UINT uiBigIndex = __bi(sSite4);
    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();

    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
    {
        pResultData[uiSiteIndex] = result;
        return;
    }
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < _DC_Dir; ++idir)
    {
        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);
        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];
        deviceSU3Vector res = deviceSU3Vector::makeZeroSU3Vector();

        //UBOOL btestpf = FALSE;
        if (!x_p_mu_Fermion.IsDirichlet())
        {
            const deviceSU3& x_Gauge_element = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIndex, idir);
            //U(x,mu) phi(x+ mu)
            res = x_Gauge_element.MulVector(pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
            if (x_p_mu_Fermion.NeedToOpposite())
            {
                res.MulReal(F(-1.0));
            }
            //btestpf = TRUE;
        }

        //UBOOL btestmf = FALSE;
        if (!x_m_mu_Fermion.IsDirichlet())
        {
            const deviceSU3 x_m_mu_Gauge_element = _deviceGetGaugeBCSU3DirSIndex(pGauge, x_m_mu_Gauge, 1);
            const deviceSU3Vector u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulVector(pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
            if (x_m_mu_Fermion.NeedToOpposite())
            {
                res.Add(u_dagger_phi_x_m_m);
            }
            else
            {
                res.Sub(u_dagger_phi_x_m_m);
            }
            //btestmf = TRUE;
        }
        //if (!btestpf && !btestmf)
        //{
        //    printf("both p fermion and m fermion are dirichlet?\n");
        //}

        res.MulReal(eta_mu);
        result.Add(res);
    }

    pResultData[uiSiteIndex].MulReal(f2am);
    //if (__cuCabsSqf(result.m_ve[0]) + __cuCabsSqf(result.m_ve[1]) + __cuCabsSqf(result.m_ve[2]) < 1.0e-10)
    //{
    //    printf("zero result!\n");
    //}

    if (bDDagger)
    {
        pResultData[uiSiteIndex].Sub(result);
    }
    else
    {
        pResultData[uiSiteIndex].Add(result);
    }

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

/**
 * Calculate Force
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_D(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIndex = __bi(sSite4);
    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
    {
        return;
    }

    //idir = mu
    for (UINT idir = 0; idir < _DC_Dir; ++idir)
    {
        if (__idx->_deviceIsBondOnSurface(uiBigIndex, static_cast<BYTE>(idir)))
        {
            continue;
        }

        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        if (x_p_mu_Fermion.IsDirichlet())
        {
            continue;
        }

        for (UINT uiR = 0; uiR < uiRational; ++uiR)
        {
            const deviceSU3Vector* phi_i = pFermionPointers[uiR];
            const deviceSU3Vector* phi_id = pFermionPointers[uiR + uiRational];
            //const deviceSU3 gaugelink = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIndex, idir);
            deviceSU3Vector toContract = pGauge[linkIndex].MulVector(phi_i[x_p_mu_Fermion.m_uiSiteIndex]);
            deviceSU3 thisTerm = deviceSU3::makeSU3ContractV(phi_id[uiSiteIndex], toContract);

            toContract = pGauge[linkIndex].MulVector(phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
            thisTerm.Add(deviceSU3::makeSU3ContractV(toContract, phi_i[uiSiteIndex]));

            if (x_p_mu_Fermion.NeedToOpposite())
            {
                thisTerm.MulReal(eta_mu * pNumerators[uiR] * F(-1.0));
            }
            else
            {
                thisTerm.MulReal(eta_mu * pNumerators[uiR]);
            }

            thisTerm.Ta();
            pForce[linkIndex].Sub(thisTerm);
        }
    }

}


#pragma endregion


void CFieldFermionKSSU3D::DOperatorKS(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceSU3Vector* pTarget = (deviceSU3Vector*)pTargetBuffer;
    const deviceSU3Vector* pSource = (const deviceSU3Vector*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    preparethread;
    _kernelDFermionKS_D << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        f2am,
        m_byFieldId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

/**
 * partial D_{st0} / partial omega
 * Make sure m_pMDNumerator and m_pRationalFieldPointers are filled
 */
void CFieldFermionKSSU3D::DerivateD0(
    void* pForce,
    const void* pGaugeBuffer) const
{
    preparethread;
    _kernelDFermionKSForce_D << <block, threads >> > (
        (const deviceSU3*)pGaugeBuffer,
        (deviceSU3*)pForce,
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        m_pRationalFieldPointers,
        m_pMDNumerator,
        m_rMD.m_uiDegree,
        m_byFieldId);
}

#pragma endregion

#pragma region Other Kernel

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_FixBoundary(
    deviceSU3Vector* pDeviceData,
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIndex = __bi(sSite4);
    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
    {
        //printf("has dirichlet\n");
        pDeviceData[uiSiteIndex] = deviceSU3Vector::makeZeroSU3Vector();
    }
}

#pragma endregion

void CFieldFermionKSSU3D::FixBoundary()
{
    preparethread;
    _kernelDFermionKS_FixBoundary << <block, threads >> > (m_pDeviceData, m_byFieldId);
}

void CFieldFermionKSSU3D::PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson)
{
    InitialField(EFIT_RandomGaussian);
    FixBoundary();
    D_MC(gaugeNum, bosonNum, gaugeFields, pBoson);
    FixBoundary();

    if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
    {
        m_fLength = Dot(this).x;
    }

}

void CFieldFermionKSSU3D::CopyTo(CField* U) const
{
    CFieldFermionKSSU3::CopyTo(U);
}

CCString CFieldFermionKSSU3D::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldFermionKSSU3::GetInfos(tab);

    SSmallInt4 boundary = appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(m_byFieldId);
    sRet = sRet + tab +appToString(boundary) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================