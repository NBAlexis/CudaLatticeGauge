//=============================================================================
// FILENAME : CFieldFermionKSTD.cu
// 
// DESCRIPTION:
// 
//
// REVISION:
//  [09/19/2020 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "CFieldFermionKSTD.h"

__BEGIN_NAMESPACE

#pragma region DOperator

#pragma region kernel

/**
*
*/
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_DT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    Real f2am,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    const UINT uiBigIndex = __bi(sSite4);
    deviceVector result = _makeZero<deviceVector>();

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
        deviceVector res = _makeZero<deviceVector>();

        //UBOOL btestpf = FALSE;
        if (!x_p_mu_Fermion.IsDirichlet())
        {
            const deviceGauge& x_Gauge_element = _deviceGetGaugeBCDirT(byGaugeFieldId, pGauge, uiBigIndex, idir);
            //U(x,mu) phi(x+ mu)
            res = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
            if (x_p_mu_Fermion.NeedToOpposite())
            {
                _mul(res, F(-1.0));
            }
            //btestpf = TRUE;
        }

        //UBOOL btestmf = FALSE;
        if (!x_m_mu_Fermion.IsDirichlet())
        {
            const deviceGauge x_m_mu_Gauge_element = _deviceGetGaugeBCDirT(pGauge, x_m_mu_Gauge, 1);
            const deviceVector u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
            if (x_m_mu_Fermion.NeedToOpposite())
            {
                _add(res, u_dagger_phi_x_m_m);
            }
            else
            {
                _sub(res, u_dagger_phi_x_m_m);
            }
            //btestmf = TRUE;
        }
        //if (!btestpf && !btestmf)
        //{
        //    printf("both p fermion and m fermion are dirichlet?\n");
        //}

        _mul(res, eta_mu);
        _add(result, res);
    }

    _mul(pResultData[uiSiteIndex], f2am);
    //if (__cuCabsSqf(result.m_ve[0]) + __cuCabsSqf(result.m_ve[1]) + __cuCabsSqf(result.m_ve[2]) < 1.0e-10)
    //{
    //    printf("zero result!\n");
    //}

    if (bDDagger)
    {
        _sub(pResultData[uiSiteIndex], result);
    }
    else
    {
        _add(pResultData[uiSiteIndex], result);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(pResultData[uiSiteIndex], fCoeff);
        break;
    case EOCT_Complex:
        _mul(pResultData[uiSiteIndex], cCoeff);
        break;
    }
}

/**
 * Calculate Force
 */
 template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_DT(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceVector* const* __restrict__ pFermionPointers,
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
            const deviceVector* phi_i = pFermionPointers[uiR];
            const deviceVector* phi_id = pFermionPointers[uiR + uiRational];
            //const deviceGauge gaugelink = _deviceGetGaugeBCSU3Dir(pGauge, uiBigIndex, idir);
            deviceVector toContract = _mulVec(pGauge[linkIndex], phi_i[x_p_mu_Fermion.m_uiSiteIndex]);
            deviceGauge thisTerm = _makeContract<deviceGauge, deviceVector>(phi_id[uiSiteIndex], toContract);

            toContract = _mulVec(pGauge[linkIndex], phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
            _add(thisTerm, _makeContract<deviceGauge, deviceVector>(toContract, phi_i[uiSiteIndex]));

            if (x_p_mu_Fermion.NeedToOpposite())
            {
                _mul(thisTerm, eta_mu * pNumerators[uiR] * F(-1.0));
            }
            else
            {
                _mul(thisTerm, eta_mu * pNumerators[uiR]);
            }

            _ta(thisTerm);
            _sub(pForce[linkIndex], thisTerm);
        }
    }

}


#pragma endregion

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTD<deviceVector, deviceGauge, vectorN>::DOperatorKS(void* pTargetBuffer, const void* pBuffer,
    const void* pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceVector* pTarget = (deviceVector*)pTargetBuffer;
    const deviceVector* pSource = (const deviceVector*)pBuffer;
    const deviceGauge* pGauge = (const deviceGauge*)pGaugeBuffer;

    preparethread;
    _kernelDFermionKS_DT << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[this->m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pMoveCache[this->m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        f2am,
        this->m_byFieldId,
        byGaugeFieldId,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

/**
 * partial D_{st0} / partial omega
 * Make sure m_pMDNumerator and m_pRationalFieldPointers are filled
 */
template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTD<deviceVector, deviceGauge, vectorN>::DerivateD0(
    void* pForce,
    const void* pGaugeBuffer,
    BYTE byGaugeFieldId) const
{
    preparethread;
    _kernelDFermionKSForce_DT << <block, threads >> > (
        (const deviceGauge*)pGaugeBuffer,
        (deviceGauge*)pForce,
        appGetLattice()->m_pIndexCache->m_pMoveCache[this->m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        this->m_pRationalFieldPointers,
        this->m_pMDNumerator,
        this->m_rMD.m_uiDegree,
        this->m_byFieldId);
}

#pragma endregion

#pragma region Other Kernel

template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_FixBoundaryT(
    deviceVector* pDeviceData,
    BYTE byFieldId)
{
    intokernalInt4;
    const UINT uiBigIndex = __bi(sSite4);
    if (__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIndex].IsDirichlet())
    {
        //printf("has dirichlet\n");
        pDeviceData[uiSiteIndex] = _makeZero<deviceVector>();
    }
}

#pragma endregion

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTD<deviceVector, deviceGauge, vectorN>::FixBoundary()
{
    preparethread;
    _kernelDFermionKS_FixBoundaryT << <block, threads >> > (this->m_pDeviceData, this->m_byFieldId);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTD<deviceVector, deviceGauge, vectorN>::PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson)
{
    this->InitialField(EFIT_RandomGaussian);
    FixBoundary();
    this->D_MC(gaugeNum, bosonNum, gaugeFields, pBoson);
    FixBoundary();

    if (NULL != appGetFermionSolver(this->m_byFieldId) && !appGetFermionSolver(this->m_byFieldId)->IsAbsoluteAccuracy())
    {
        this->m_fLength = Dot(this).x;
    }

}

template<typename deviceVector, typename deviceGauge, INT vectorN>
CCString CFieldFermionKSTD<deviceVector, deviceGauge, vectorN>::GetInfos(const CCString& tab) const
{
    CCString sRet = CFieldFermionKST<deviceVector, deviceGauge, vectorN>::GetInfos(tab);

    SSmallInt4 boundary = appGetLattice()->m_pIndex->GetBoudanryCondition()->GetFieldBC(this->m_byFieldId);
    sRet = sRet + tab +appToString(boundary) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================