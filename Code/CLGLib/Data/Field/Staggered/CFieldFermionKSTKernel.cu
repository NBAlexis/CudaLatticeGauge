//=============================================================================
// FILENAME : CFieldFermionKSTKernel.cu
// 
// DESCRIPTION:
// This is the device implementations of Wilson fermion
//
// This implementation assumes SU3 and square lattice
//
// REVISION:
//  [07/21/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "Data/Field/Gauge/CFieldGaugeLink.h"
#include "CFieldFermionKSTKernel.h"
#include "CFieldFermionKST.h"

__BEGIN_NAMESPACE

template<typename deviceVector, typename deviceGauge, INT vectorN>
UINT CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::TestAntiHermitianS(BYTE byFieldId, const CFieldGauge* pGauge)
{
    const UINT uiVolume = _HC_Volume;
    const UINT uiRealVolume = vectorN * uiVolume;
    CLGComplex* matrixElement = (CLGComplex*)malloc(sizeof(CLGComplex) * uiRealVolume * uiRealVolume);
    deviceVector* hostData = (deviceVector*)malloc(sizeof(deviceVector) * uiVolume);
    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* v = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(byFieldId));

    for (UINT i = 0; i < uiVolume; ++i)
    {
        const SSmallInt4 point = __hostSiteIndexToInt4(i);
        for (UINT j = 0; j < vectorN; ++j)
        {
            SFermionBosonSource source;
            source.m_byColorIndex = static_cast<BYTE>(j);
            source.m_eSourceType = EFS_Point;
            source.m_sSourcePoint = point;
            v->InitialAsSource(source);
            v->D0S(pGauge);

            checkCudaErrors(cudaMemcpy(hostData, v->m_pDeviceData, sizeof(deviceVector) * uiVolume, cudaMemcpyDeviceToHost));

            const UINT x = i * vectorN + j;
            for (UINT k = 0; k < uiVolume; ++k)
            {
                for (UINT kcolor = 0; kcolor < vectorN; ++kcolor)
                {
                    matrixElement[(vectorN * k + kcolor) * uiRealVolume + x] = _make_cuComplex(_element(hostData[k], 2 * kcolor), _element(hostData[k], 2 * kcolor + 1));
                }
            }
            appGeneral(_T("%d / %d have been done\n"), x, uiRealVolume);
        }
    }

    UINT uiE = 0;
    UINT uiWrong = 0;
    //List all results
    for (UINT i = 0; i < uiRealVolume * uiRealVolume; ++i)
    {
        const UINT x = i / uiRealVolume;
        const UINT y = i % uiRealVolume;
        const SSmallInt4 xSite = __hostSiteIndexToInt4(x / vectorN);
        const SSmallInt4 ySite = __hostSiteIndexToInt4(y / vectorN);
        const UINT daggerIdx = y * uiRealVolume + x;
        const BYTE cx = x % vectorN;
        const BYTE cy = y % vectorN;

        if (_cuCabsf(matrixElement[i]) > F(0.0000001))
        {
            ++uiE;
            if (appAbs(matrixElement[i].x + matrixElement[daggerIdx].x) > F(0.0000001)
                || appAbs(matrixElement[i].y - matrixElement[daggerIdx].y) > F(0.0000001))
            {
                ++uiWrong;
                appGeneral(_T("[(%d, %d, %d, %d)_(%d)-(%d, %d, %d, %d)_(%d)]: D = %f + %f I   Ddagger = %f + %f I\n"),
                    xSite.x, xSite.y, xSite.z, xSite.w, cx,
                    ySite.x, ySite.y, ySite.z, ySite.w, cy,
                    matrixElement[i].x, matrixElement[i].y,
                    matrixElement[daggerIdx].x, matrixElement[daggerIdx].y);
            }
        }
    }
    v->Return();
    appSafeFree(matrixElement);
    appSafeFree(hostData);
    appGeneral(_T("%d none zero element checked, %d wrong found...\n"), uiE, uiWrong);
    return uiWrong;
}

#pragma region DOperator

#pragma region kernel

/**
* Dks = 2am + \sum _{\mu} \eta_{\mu} (n) (U_{\mu}(n) \delta _{n,n+\mu} -U^+_{\mu}(n-\mu) \delta _{n,n-\mu})
* U act on su3
* gamma act on spinor
*
* If bDagger, it is just \eta 5(n) Dks \eta _5(n)
* A easier way is D_{ks}(m=0) is anti-Hermitian, so D^+ = - D_{ks,m=0} + 2am
*
*/
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKST(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    Real f2am,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernaldir;

    deviceVector result = _makeZero<deviceVector>();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        const deviceGauge& x_Gauge_element = pGauge[linkIndex];
        deviceGauge x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            _dagger(x_m_mu_Gauge_element);
        }

        //U(x,mu) phi(x+ mu)
        deviceVector u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            _mul(u_phi_x_p_m, F(-1.0));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceVector u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            _add(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        else
        {
            _sub(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        _mul(u_phi_x_p_m, eta_mu);
        _add(result, u_phi_x_p_m);
    }

    _mul(pResultData[uiSiteIndex], f2am);
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
 * For some strange boundary condition
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSPlusEtaT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    Real f2am,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernaldir;

    deviceVector result = _makeZero<deviceVector>();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //This is in fact, -1 * eta(n + mu)
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        const Real eta_mu2 = (1 == ((pEtaTable[x_m_mu_Fermion.m_uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        const deviceGauge& x_Gauge_element = pGauge[linkIndex];
        deviceGauge x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            _dagger(x_m_mu_Gauge_element);
        }

        //U(x,mu) phi(x+ mu)
        deviceVector u_phi_x_p_m = _mulVec(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            _mul(u_phi_x_p_m, F(-1.0) * eta_mu);
        }
        else
        {
            _mul(u_phi_x_p_m, eta_mu);
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceVector u_dagger_phi_x_m_m = _mulVec(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        _mul(u_dagger_phi_x_m_m, eta_mu2);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            _add(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        else
        {
            _sub(u_phi_x_p_m, u_dagger_phi_x_m_m);
        }
        //_mul(u_phi_x_p_m, eta_mu);
        _add(result, u_phi_x_p_m);
    }

    _mul(pResultData[uiSiteIndex], f2am);
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
_kernelDFermionKSForceT(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId)
{
    intokernaldir;

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];

        for (UINT uiR = 0; uiR < uiRational; ++uiR)
        {
            const deviceVector* phi_i = pFermionPointers[uiR];
            const deviceVector* phi_id = pFermionPointers[uiR + uiRational];

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
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKS(UBOOL bEachSiteEta, deviceVector* pTarget, const deviceVector* pSource,
    const deviceGauge* pGauge, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    if (bEachSiteEta)
    {
        _kernelDFermionKSPlusEtaT << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            pTarget,
            f2am,
            byFieldId,
            bDagger,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }
    else
    {
        _kernelDFermionKST << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            pTarget,
            f2am,
            byFieldId,
            bDagger,
            eOCT,
            fRealCoeff,
            cCmpCoeff);
    }
}

/**
 * partial D_{st0} / partial omega
 * Make sure m_pMDNumerator and m_pRationalFieldPointers are filled
 */
template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DerivateD0(
    const deviceVector* pFermion,
    BYTE byFieldId,
    deviceGauge* pForce,
    const deviceGauge* pGaugeBuffer,
    BYTE byGaugeFieldId, 
    const deviceVector* const* pRationalFields, 
    const Real* pNumerator, 
    UINT uiRationApproxOrder)
{
    preparethread;
    _kernelDFermionKSForceT << <block, threads >> > (
        pGaugeBuffer,
        pForce,
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pRationalFields,
        pNumerator,
        uiRationApproxOrder,
        byFieldId);
}

#pragma endregion

#pragma region DOperator Dirichlet

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
            const deviceGauge x_m_mu_Gauge_element = _deviceGetGaugeBCT(byGaugeFieldId, pGauge, x_m_mu_Gauge);
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
    BYTE byFieldId,
    BYTE byGaugeFieldId)
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
        if (__idx->_deviceIsBondOnSurface(uiBigIndex, byGaugeFieldId, static_cast<BYTE>(idir)))
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
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKS_D(UBOOL bEachSiteEta, deviceVector* pTarget, const deviceVector* pSource,
    const deviceGauge* pGauge, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    _kernelDFermionKS_DT << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        f2am,
        byFieldId,
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
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DerivateD0_D(
    const deviceVector* pFermion,
    BYTE byFieldId,
    deviceGauge* pForce,
    const deviceGauge* pGaugeBuffer,
    BYTE byGaugeFieldId,
    const deviceVector* const* pRationalFields,
    const Real* pNumerator,
    UINT uiRationApproxOrder)
{
    preparethread;
    _kernelDFermionKSForce_DT << <block, threads >> > (
        pGaugeBuffer,
        pForce,
        appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pRationalFields,
        pNumerator,
        uiRationApproxOrder,
        byFieldId,
        byGaugeFieldId);
}

#pragma endregion


#pragma region Helper functions to implement higher orders

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLinkT(
    const deviceGauge* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fCoefficient,
    BYTE byEtaIndex,
    const INT* __restrict__ path,
    BYTE pathLength)
{
    intokernalInt4;
    INT pathLeft[_kLinkMaxLength];
    INT pathRight[_kLinkMaxLength];
    for (BYTE iSeperation = 0; iSeperation <= pathLength; ++iSeperation)
    {
        BYTE LLength = 0;
        BYTE RLength = 0;

        _deviceSeperate(path, iSeperation, pathLength, pathLeft, pathRight, LLength, RLength);

        const UBOOL bHasLeft = (LLength > 0) && (pathLeft[0] > 0);
        const UBOOL bHasRight = (RLength > 0) && (pathRight[0] > 0);

        if (bHasLeft || bHasRight)
        {
            //=================================
            // 1. Find n1, n2
            const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, pathLeft, LLength);
            const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, pathRight, RLength);
            const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
            const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];
            INT iEtaMu1 = (pEtaTable[sn1.m_uiSiteIndex] >> byEtaIndex);
            if (sn1.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            if (sn2.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            //=================================
            // 2. Find V(n,n1), V(n,n2)
            const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
            const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

            for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
            {
                const deviceVector* phi_i = pFermionPointers[rfieldId];
                const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                //=================================
                // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                deviceVector phi1 = _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]);
                deviceVector phi2 = _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]);

                deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                //This Add is required by partial(D^+D)
                phi1 = _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]);
                phi2 = _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]);
                _add(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
                _ta(res);
                _mul(res, fCoefficient * pNumerators[rfieldId]);

                if (bHasLeft)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        _sub(pForce[linkIndex], res);
                    }
                    else
                    {
                        _add(pForce[linkIndex], res);
                    }
                }

                if (bHasRight)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
                    }
                }
            }
        }
    }
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_OneLinkT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fCoefficient,
    const INT* __restrict__ path,
    BYTE pathLength,
    BYTE byEtaIdx,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;
    INT pathBuffer[_kLinkMaxLength];
    deviceVector result = _makeZero<deviceVector>();

    SSmallInt4 siten = _deviceSmallInt4OffsetC(sSite4, path, pathLength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten)];
    deviceGauge vn = _deviceLinkT(pGauge, sSite4, pathLength, byGaugeFieldId, path);
    INT etamu = (pEtaTable[uiSiteIndex] >> byEtaIdx) & 1;
    if (sn1.NeedToOpposite())
    {
        etamu = etamu + 1;
    }
    if (etamu & 1)
    {
        _sub(result, _mulVec(vn, pDeviceData[sn1.m_uiSiteIndex]));
    }
    else
    {
        _add(result, _mulVec(vn, pDeviceData[sn1.m_uiSiteIndex]));
    }

    _devicePathDagger(path, pathBuffer, pathLength);
    siten = _deviceSmallInt4OffsetC(sSite4, pathBuffer, pathLength);
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten)];
    vn = _deviceLinkT(pGauge, sSite4, pathLength, byGaugeFieldId, pathBuffer);
    etamu = (pEtaTable[sn2.m_uiSiteIndex] >> byEtaIdx) & 1;
    if (sn2.NeedToOpposite())
    {
        etamu = etamu + 1;
    }
    if (etamu & 1)
    {
        _add(result, _mulVec(vn, pDeviceData[sn2.m_uiSiteIndex]));
    }
    else
    {
        _sub(result, _mulVec(vn, pDeviceData[sn2.m_uiSiteIndex]));
    }

    if (bDDagger)
    {
        fCoefficient = fCoefficient * F(-1.0);
    }
    _mul(result, fCoefficient);

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    }

    _add(pResultData[uiSiteIndex], result);
}


template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_OnlyMassT(
    const deviceVector* __restrict__ pDeviceData,
    deviceVector* pResultData,
    Real f2am,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernal;
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];
    _mul(pResultData[uiSiteIndex], f2am);

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


template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::OnlyMass(const deviceVector* pSource, deviceVector* pTarget, Real fm, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    _kernelDFermionKS_OnlyMassT << <block, threads >> > (
        pSource,
        pTarget,
        fm,
        eOCT,
        fRealCoeff,
        cCmpCoeff
        );
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::OneLinkS(
    const deviceVector* pSource, 
    BYTE byFieldId,
    const deviceGauge* pGauge,
    BYTE byGaugeFieldId,
    deviceVector* pTarget,
    Real fCoefficient,
    const INT* pDevicePath,
    BYTE pathLength,
    BYTE byEtaIdx,
    UBOOL bDagger,
    EOperatorCoefficientType eOCT,
    Real fRealCoeff,
    const CLGComplex& cCmpCoeff)
{
    assert(pathLength <= _kLinkMaxLength);
    preparethread;
    _kernelDFermionKS_OneLinkT << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fCoefficient,
        pDevicePath,
        pathLength,
        byEtaIdx,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff
        );
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::OneLinkForceS(
    const deviceVector* pFermion, 
    BYTE byFieldId,
    const deviceGauge* pGauge,
    BYTE byGaugeFieldId,
    deviceGauge* pForce,
    Real fCoefficient,
    const INT* pDevicePath,
    BYTE pathLength,
    BYTE byEtaIdx, 
    const deviceVector* const* pRationalFields, 
    const Real* pNumerator, 
    UINT uiRationApproxOrder)
{
    assert(pathLength <= _kLinkMaxLength);
    preparethread;
    _kernelDFermionKSForce_WithLinkT << <block, threads >> > (
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pForce,
        pRationalFields,
        pNumerator,
        uiRationApproxOrder,
        byFieldId,
        byGaugeFieldId,
        fCoefficient,
        byEtaIdx,
        pDevicePath,
        pathLength
        );
}

#pragma endregion

#pragma region Gamma for KS

#pragma region gamma kernels

/**
 * gamma_i ^dagger = gamma_i for i = 1,2,3,4,5, so 'bDDagger' is used for only coefficient
 *
 * 1/2a eta_mu(x) sum_{mu=+-} psibar(x) U_mu psi(x+mu)
 * the 1/2a is absorbed (be sure to add when measuring)
 */
    template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma1234T(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bDDagger,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byDir)
{
    intokernal;

    const Real eta_mu = ((pEtaTable[uiSiteIndex] >> byDir) & 1) ? F(-1.0) : F(1.0);
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byDir);
    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    const deviceGauge& x_Gauge_element = pGauge[linkIndex];
    deviceGauge x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, byDir)];
    if (x_m_mu_Gauge.NeedToDagger())
    {
        _dagger(x_m_mu_Gauge_element);
    }

    deviceVector result = _mulVec(x_Gauge_element, pOther[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        _mul(result, F(-1.0));
    }
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        _sub(result, _mulVec(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _add(result, _mulVec(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    //result.MulReal(eta_mu);

    if (bDDagger)
    {
        fGammCoefficient = -fGammCoefficient;
    }
    _mul(result, _make_cuComplex(F(0.0), fGammCoefficient * eta_mu));

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

/**
* This function applys i Gamma_i on the field
*/
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGammaEta1234T(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bDDagger,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byDir)
{
    intokernal;

    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byDir);
    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    BYTE eta_mu = pEtaTable[uiSiteIndex] >> byDir;
    BYTE eta_mu2 = pEtaTable[x_m_mu_Fermion.m_uiSiteIndex] >> byDir;

    const deviceGauge& x_Gauge_element = pGauge[linkIndex];
    deviceGauge x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, x_m_mu_Gauge.m_byDir)];
    if (x_m_mu_Gauge.NeedToDagger())
    {
        _dagger(x_m_mu_Gauge_element);
    }

    deviceVector result = _mulVec(x_Gauge_element, pOther[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        eta_mu = eta_mu + 1;
    }

    if (eta_mu & 1)
    {
        _mul(result, F(-1.0));
    }

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        eta_mu2 = eta_mu2 + 1;
    }
    if (eta_mu2 & 1)
    {
        _sub(result, _mulVec(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        _add(result, _mulVec(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }

    if (bDDagger)
    {
        fGammCoefficient = -fGammCoefficient;
    }
    _mul(result, _make_cuComplex(F(0.0), fGammCoefficient));

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelFermionKSForceGamma1234T(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE idir,
    Real fCoeff,
    BYTE byFieldId)
{
    intokernal;

    const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> idir) & 1)) ? F(-1.0) : F(1.0);
    //x, mu
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];

    for (UINT uiR = 0; uiR < uiRational; ++uiR)
    {
        const deviceVector* phi_i = pFermionPointers[uiR];
        const deviceVector* phi_id = pFermionPointers[uiR + uiRational];

        deviceVector toContract = _mulVec(pGauge[linkIndex], phi_i[x_p_mu_Fermion.m_uiSiteIndex]);
        deviceGauge thisTerm = _makeContract<deviceGauge, deviceVector>(phi_id[uiSiteIndex], toContract);

        toContract = _mulVec(pGauge[linkIndex], phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
        _sub(thisTerm, _makeContract<deviceGauge, deviceVector>(toContract, phi_i[uiSiteIndex]));

        if (x_p_mu_Fermion.NeedToOpposite())
        {
            _mul(thisTerm, eta_mu * pNumerators[uiR] * F(-1.0));
        }
        else
        {
            _mul(thisTerm, eta_mu * pNumerators[uiR]);
        }

        _mul(thisTerm, _make_cuComplex(F(0.0), fCoeff));
        _ta(thisTerm);
        _sub(pForce[linkIndex], thisTerm);
    }
}


template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGammaSigmaIJT(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    UBOOL bDDagger,
    Real fGammCoefficient,
    BYTE byEtaShift,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    SBYTE byDir1,
    SBYTE byDir2,
    BYTE byFieldId,
    BYTE byGaugeFieldId)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

#pragma unroll
    for (UINT idx = 0; idx < 4; ++idx)
    {
        const UBOOL bPlus12[2] = { (0 != (idx & 1)), (0 != (idx & 2)) };

        SSmallInt4 sOffset = sSite4;
        SBYTE dim12[2] = {
            bPlus12[0] ? static_cast<SBYTE>(byDir1 + 1) : static_cast<SBYTE>(-byDir1 - 1),
            bPlus12[1] ? static_cast<SBYTE>(byDir2 + 1) : static_cast<SBYTE>(-byDir2 - 1)
        };
        sOffset.m_byData4[byDir1] = sOffset.m_byData4[byDir1] + (bPlus12[0] ? 1 : -1);
        sOffset.m_byData4[byDir2] = sOffset.m_byData4[byDir2] + (bPlus12[1] ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceVector right = _mulVec(_devicePlaneDiagonalT(pGauge, sSite4, byGaugeFieldId, dim12[0], dim12[1]), pOther[sTargetBigIndex.m_uiSiteIndex]);

        //eta12 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta2 = _deviceEta2(pEtaTable[uiSiteIndex], byDir1, byDir2);
        if (1 == byEtaShift || 2 == byEtaShift)
        {
            //type 1, it is not the corner, and if it cross the Y-boundary
            if (!bPlus12[byEtaShift - 1])
            {
                eta2 = _deviceEta2(pEtaTable[sTargetBigIndex.m_uiSiteIndex], byDir1, byDir2) + 1;
            }
        }

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta2 = eta2 + 1;
        }

        if (eta2 & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, -F(0.5) * fGammCoefficient);
    }
    else
    {
        _mul(result, F(0.5) * fGammCoefficient);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

/**
 * gamma5i corresponds to diagonal links of cubic in the other 3 dimensions.
 * For each diagonal link, there are 6 different ways to add the gauge field
 * We simply use the average of all 6 links
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma51234T(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
    UBOOL bDDagger,
    BYTE byMissingDir,
    SBYTE byEtaShift,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byFieldId,
    BYTE byGaugeFieldId)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

#pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlus123[3] = { (0 != (idx & 1)), (0 != (idx & 2)), (0 != (idx & 4)) };

        SSmallInt4 sOffset = sSite4;
        SBYTE dim123[3] = { 0, 0, 0 };
        BYTE byDimIndex = 0;
        for (SBYTE byCubeDir = 0; byCubeDir < 4; ++byCubeDir)
        {
            if (byCubeDir != static_cast<SBYTE>(byMissingDir))
            {
                sOffset.m_byData4[byCubeDir] = sOffset.m_byData4[byCubeDir] + (bPlus123[byDimIndex] ? 1 : -1);
                dim123[byDimIndex] = bPlus123[byDimIndex] ? (byCubeDir + 1) : (-byCubeDir - 1);
                byDimIndex++;
            }
        }

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceVector right = _mulVec(_deviceCubicDiagonalT(pGauge, sSite4, byGaugeFieldId, dim123[0], dim123[1], dim123[2]), pOther[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta124 of site is almost always -target, so use left or right is same
        //The only exception is on the boundary
        INT eta3 = _deviceEta3(sSite4, byMissingDir);
        if (byEtaShift >= 0 && byEtaShift <= 2 && !bPlus123[byEtaShift])
        {
            eta3 = _deviceEta3(site_target, byMissingDir) + 1;
        }

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta3 = eta3 + 1;
        }

        if (eta3 & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, -F(0.25) * fGammCoefficient);
    }
    else
    {
        _mul(result, F(0.25) * fGammCoefficient);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

/**
 *
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGamma5T(
    deviceVector* pMe,
    const deviceVector* __restrict__ pOther,
    const deviceGauge* __restrict__ pGauge,
    UBOOL bDDagger,
    UBOOL bEtaShift,
    Real fGammCoefficient,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff,
    BYTE byFieldId,
    BYTE byGaugeFieldId)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

#pragma unroll
    for (UINT idx = 0; idx < 16; ++idx)
    {
        const UBOOL bPlus1234[4] = { (0 != (idx & 1)), (0 != (idx & 2)), (0 != (idx & 4)), (0 != (idx & 8)) };

        SSmallInt4 sOffset = sSite4;
        SBYTE dim1234[4] =
        {
            bPlus1234[0] ? static_cast<SBYTE>(1) : static_cast<SBYTE>(-1),
            bPlus1234[1] ? static_cast<SBYTE>(2) : static_cast<SBYTE>(-2),
            bPlus1234[2] ? static_cast<SBYTE>(3) : static_cast<SBYTE>(-3),
            bPlus1234[3] ? static_cast<SBYTE>(4) : static_cast<SBYTE>(-4)
        };
        sOffset.m_byData4[0] = sOffset.m_byData4[0] + (bPlus1234[0] ? 1 : -1);
        sOffset.m_byData4[1] = sOffset.m_byData4[1] + (bPlus1234[1] ? 1 : -1);
        sOffset.m_byData4[2] = sOffset.m_byData4[2] + (bPlus1234[2] ? 1 : -1);
        sOffset.m_byData4[3] = sOffset.m_byData4[3] + (bPlus1234[3] ? 1 : -1);

        //We have anti-periodic boundary, so we need to use index out of lattice to get the correct sign
        const SIndex& sTargetBigIndex = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(sOffset)];

        const deviceVector right = _mulVec(_deviceHyperCubicDiagonalT(pGauge, sSite4, byGaugeFieldId, dim1234[0], dim1234[1], dim1234[2], dim1234[3]), pOther[sTargetBigIndex.m_uiSiteIndex]);
        const SSmallInt4 site_target = __deviceSiteIndexToInt4(sTargetBigIndex.m_uiSiteIndex);

        //eta51 is gamma5 (x+z)
        INT eta4 = _deviceEta3(sSite4, 0);
        if (bEtaShift && !bPlus1234[3])
        {
            //target is almost always site4 except for boundaries
            eta4 = _deviceEta3(site_target, 0);
        }

        if (sTargetBigIndex.NeedToOpposite())
        {
            eta4 = eta4 + 1;
        }

        if (eta4 & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, _make_cuComplex(F(0.0), -F(0.125) * fGammCoefficient));
    }
    else
    {
        _mul(result, _make_cuComplex(F(0.0), F(0.125) * fGammCoefficient));
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    }

    _add(pMe[uiSiteIndex], result);
}

/**
 * similar as _kernelDFermionKSForce_WithLink.
 * but _kernelDFermionKSForce_WithLink use gamma mu as the gamma matrix
 * we use sigma12
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLink_SigmaIJT(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    const BYTE* __restrict__ pEtaTable,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fGammCoefficient,
    BYTE bDir1,
    BYTE bDir2,
    const INT* __restrict__ path)
{
    intokernalInt4;
    INT pathLeft[2];
    INT pathRight[2];
    for (BYTE iSeperation = 0; iSeperation <= 2; ++iSeperation)
    {
        BYTE LLength = 0;
        BYTE RLength = 0;

        _deviceSeperate(path, iSeperation, 2, pathLeft, pathRight, LLength, RLength);

        const UBOOL bHasLeft = (LLength > 0) && (pathLeft[0] > 0);
        const UBOOL bHasRight = (RLength > 0) && (pathRight[0] > 0);

        if (bHasLeft || bHasRight)
        {
            //=================================
            // 1. Find n1, n2
            const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, pathLeft, LLength);
            const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, pathRight, RLength);
            const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
            const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];
            INT iEtaMu1 = _deviceEta2(pEtaTable[sn1.m_uiSiteIndex], bDir1, bDir2);

            if (sn1.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            if (sn2.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            //=================================
            // 2. Find V(n,n1), V(n,n2)
            const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
            const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

            for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
            {
                const deviceVector* phi_i = pFermionPointers[rfieldId];
                const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                //=================================
                // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                deviceVector phi1 = _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]);
                deviceVector phi2 = _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]);
                //deviceVector phi3 = _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]);
                //deviceVector phi4 = _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]);

                deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                //This Add is required by partial(D^+D)
                phi1 = _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]);
                phi2 = _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]);
                _add(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
                //res.Add(deviceGauge::makeSU3ContractV(phi4, phi3));
                _ta(res);
                _mul(res, fGammCoefficient * pNumerators[rfieldId]);

                if (bHasLeft)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        _sub(pForce[linkIndex], res);
                    }
                    else
                    {
                        _add(pForce[linkIndex], res);
                    }
                }

                if (bHasRight)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
                    }
                }
            }
        }
    }
}

/**
 * similar as _kernelDFermionKSForce_WithLink.
 * but _kernelDFermionKSForce_WithLink use gamma mu as the gamma matrix
 * we use gamma5i
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLink_Gamma51234T(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fGammCoefficient,
    BYTE byMissingDir,
    const INT* __restrict__ path)
{
    intokernalInt4;
    INT pathLeft[3];
    INT pathRight[3];
    for (BYTE iSeperation = 0; iSeperation <= 3; ++iSeperation)
    {
        BYTE LLength = 0;
        BYTE RLength = 0;

        _deviceSeperate(path, iSeperation, 3, pathLeft, pathRight, LLength, RLength);

        const UBOOL bHasLeft = (LLength > 0) && (pathLeft[0] > 0);
        const UBOOL bHasRight = (RLength > 0) && (pathRight[0] > 0);

        if (bHasLeft || bHasRight)
        {
            //=================================
            // 1. Find n1, n2
            const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, pathLeft, LLength);
            const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, pathRight, RLength);
            const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
            const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];
            INT iEtaMu1 = _deviceEta3(__deviceSiteIndexToInt4(sn1.m_uiSiteIndex), byMissingDir);
            //INT iEtaMu2 = _deviceEta3(__deviceSiteIndexToInt4(sn2.m_uiSiteIndex), byMissingDir);

            if (sn1.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            if (sn2.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            //=================================
            // 2. Find V(n,n1), V(n,n2)
            const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
            const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

            for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
            {
                const deviceVector* phi_i = pFermionPointers[rfieldId];
                const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                //=================================
                // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                deviceVector phi1 = _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]);
                deviceVector phi2 = _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]);
                //deviceVector phi3 = _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]);
                //deviceVector phi4 = _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]);

                deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                //This Add is required by partial(D^+D)
                phi1 = _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]);
                phi2 = _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]);
                _add(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
                //res.Add(deviceGauge::makeSU3ContractV(phi4, phi3));
                _ta(res);
                _mul(res, fGammCoefficient * pNumerators[rfieldId]);

                if (bHasLeft)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        _sub(pForce[linkIndex], res);
                    }
                    else
                    {
                        _add(pForce[linkIndex], res);
                    }
                }

                if (bHasRight)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
                    }
                }
            }
        }
    }
}

//similar as above
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLink_Gamma5T(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    Real fGammCoefficient,
    const INT* __restrict__ path)
{
    intokernalInt4;
    INT pathLeft[4];
    INT pathRight[4];

    //if (0 == uiSiteIndex)
    //{
    //    printf("%f, %d %d %d %d\n", fGammCoefficient, path[0], path[1], path[2], path[3]);
    //}

    for (BYTE iSeperation = 0; iSeperation <= 4; ++iSeperation)
    {
        BYTE LLength = 0;
        BYTE RLength = 0;

        _deviceSeperate(path, iSeperation, 4, pathLeft, pathRight, LLength, RLength);

        const UBOOL bHasLeft = (LLength > 0) && (pathLeft[0] > 0);
        const UBOOL bHasRight = (RLength > 0) && (pathRight[0] > 0);

        if (bHasLeft || bHasRight)
        {
            //=================================
            // 1. Find n1, n2
            const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, pathLeft, LLength);
            const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, pathRight, RLength);
            const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
            const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];
            INT iEtaMu1 = _deviceEta3(__deviceSiteIndexToInt4(sn1.m_uiSiteIndex), 0);
            //INT iEtaMu2 = _deviceEta3(__deviceSiteIndexToInt4(sn2.m_uiSiteIndex), byMissingDir);

            if (sn1.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            if (sn2.NeedToOpposite())
            {
                iEtaMu1 = iEtaMu1 + 1;
            }
            //=================================
            // 2. Find V(n,n1), V(n,n2)
            const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
            const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

            for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
            {
                const deviceVector* phi_i = pFermionPointers[rfieldId];
                const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

                //=================================
                // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                deviceVector phi1 = _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]);
                deviceVector phi2 = _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]);
                //deviceVector phi3 = _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]);
                //deviceVector phi4 = _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]);

                deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
                //This Add is required by partial(D^+D)
                phi1 = _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]);
                phi2 = _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]);
                _sub(res, _makeContract<deviceGauge, deviceVector>(phi1, phi2));
                //res.Add(deviceGauge::makeSU3ContractV(phi4, phi3));

                _mul(res, _make_cuComplex(F(0.0), fGammCoefficient * pNumerators[rfieldId]));
                //_mul(res, fGammCoefficient * pNumerators[rfieldId]);
                _ta(res);

                if (bHasLeft)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        _sub(pForce[linkIndex], res);
                    }
                    else
                    {
                        _add(pForce[linkIndex], res);
                    }
                }

                if (bHasRight)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        _add(pForce[linkIndex], res);
                    }
                    else
                    {
                        _sub(pForce[linkIndex], res);
                    }
                }
            }
        }
    }
}


#pragma endregion


template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::appApplyGammaKS(
    void* pTargetBuffer,
    const void* pBuffer,
    const void* pGaugeBuffer,
    EGammaMatrix eGamma,
    UBOOL bShiftCenter,
    UBOOL bDagger,
    Real fGammaCoeff,
    EOperatorCoefficientType eOCT,
    Real fRealCoeff,
    CLGComplex cCmpCoeff,
    BYTE byFieldID,
    BYTE byGaugeFieldID)
{
    deviceVector* pTarget = (deviceVector*)pTargetBuffer;
    const deviceVector* pSource = (const deviceVector*)pBuffer;
    const deviceGauge* pGauge = (const deviceGauge*)pGaugeBuffer;
    preparethread;

    switch (eGamma)
    {
    case GAMMA1:
    case GAMMA2:
    case GAMMA3:
    case GAMMA4:
    {
        INT iDir = static_cast<INT>(eGamma) - 1;

        if (bShiftCenter)
        {
            _kernelKSApplyGammaEta1234T << <block, threads >> > (
                pTarget,
                pSource,
                pGauge,
                appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldID],
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bDagger,
                fGammaCoeff,
                eOCT,
                fRealCoeff,
                cCmpCoeff,
                static_cast<BYTE>(iDir));
        }
        else
        {
            _kernelKSApplyGamma1234T << <block, threads >> > (
                pTarget,
                pSource,
                pGauge,
                appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[byFieldID],
                appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bDagger,
                fGammaCoeff,
                eOCT,
                fRealCoeff,
                cCmpCoeff,
                static_cast<BYTE>(iDir));
        }
    }
    break;
    case SIGMA12:
        if (bShiftCenter)
        {
            appCrucial(_T("Sigma 12 in projective plane boundary condition is not supported!\n"));
        }
        else
        {
            _kernelKSApplyGammaSigmaIJT << <block, threads >> > (
                pTarget,
                pSource,
                pGauge,
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                bDagger,
                fGammaCoeff,
                bShiftCenter ? 3 : 0,
                eOCT,
                fRealCoeff,
                cCmpCoeff,
                0,
                1,
                byFieldID,
                byGaugeFieldID);
        }
        break;
    case SIGMA31:
        //this is sigma 13
        _kernelKSApplyGammaSigmaIJT << <block, threads >> > (
            pTarget,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fGammaCoeff,
            bShiftCenter ? 1 : 0,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            0,
            2,
            byFieldID,
            byGaugeFieldID);
        break;
    case SIGMA41:
        //this is sigma 14
        _kernelKSApplyGammaSigmaIJT << <block, threads >> > (
            pTarget,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fGammaCoeff,
            bShiftCenter ? 1 : 0,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            0,
            3,
            byFieldID,
            byGaugeFieldID);
        break;
    case SIGMA23:
        _kernelKSApplyGammaSigmaIJT << <block, threads >> > (
            pTarget,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fGammaCoeff,
            bShiftCenter ? 2 : 0,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            1,
            2,
            byFieldID,
            byGaugeFieldID);
        break;
    case SIGMA42:
        //this is sigma 24
        _kernelKSApplyGammaSigmaIJT << <block, threads >> > (
            pTarget,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fGammaCoeff,
            bShiftCenter ? 2 : 0,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            1,
            3,
            byFieldID,
            byGaugeFieldID);
        break;
    case SIGMA43:
        //this is sigma 34
        _kernelKSApplyGammaSigmaIJT << <block, threads >> > (
            pTarget,
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            bDagger,
            fGammaCoeff,
            0,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            2,
            3,
            byFieldID,
            byGaugeFieldID);
        break;
    case GAMMA51:
    case GAMMA52:
    case GAMMA53:
    case GAMMA54:
    {
        //eta shift is:
        //x->y
        //y->x
        //z->t
        //t->z

        const BYTE byMissingDir = static_cast<BYTE>(eGamma - GAMMA51);
        SBYTE etaShift = -1;
        if (bShiftCenter)
        {
            if (byMissingDir < 2)
            {
                etaShift = 0;
            }
            else
            {
                etaShift = 2;
            }
        }

        _kernelKSApplyGamma51234T << <block, threads >> > (
            pTarget,
            pSource,
            pGauge,
            bDagger,
            byMissingDir,
            etaShift,
            fGammaCoeff,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            byFieldID,
            byGaugeFieldID);
    }
    break;
    case GAMMA5:
        _kernelKSApplyGamma5T << <block, threads >> > (
            pTarget,
            pSource,
            pGauge,
            bDagger,
            bShiftCenter,
            fGammaCoeff,
            eOCT,
            fRealCoeff,
            cCmpCoeff,
            byFieldID,
            byGaugeFieldID);
        break;
    default:
        appGeneral(_T("not implimented!\n"));
        break;

    }
}


template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::GammaKSForce(
    void* pForce,
    const void* pGaugeBuffer,
    const deviceVector* const* pRationalFields,
    const Real* pRationalNumerator,
    UINT uiRationalDegree,
    Real fCoeff,
    EGammaMatrix eGamma,
    INT* devicePathBuffer,
    BYTE byFieldID,
    BYTE byGaugeFieldID)
{
    preparethread;

    switch (eGamma)
    {
    case GAMMA1:
    case GAMMA2:
    case GAMMA3:
    case GAMMA4:
    {
        BYTE byDir = static_cast<BYTE>(eGamma) - 1;
        _kernelFermionKSForceGamma1234T << <block, threads >> > (
            (const deviceGauge*)pGaugeBuffer,
            (deviceGauge*)pForce,
            appGetLattice()->m_pIndexCache->m_pMoveCache[byFieldID],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            pRationalFields,
            pRationalNumerator,
            uiRationalDegree,
            byDir,
            fCoeff,
            byFieldID);
    }
    break;
    case SIGMA12:
    case SIGMA31:
    case SIGMA41:
    case SIGMA23:
    case SIGMA42:
    case SIGMA43:
    {
        SBYTE byDirs[2] = { 0, 1 };
        if (SIGMA31 == eGamma)
        {
            byDirs[0] = 0; byDirs[1] = 2;
        }
        else if (SIGMA41 == eGamma)
        {
            byDirs[0] = 0; byDirs[1] = 3;
        }
        else if (SIGMA23 == eGamma)
        {
            byDirs[0] = 1; byDirs[1] = 2;
        }
        else if (SIGMA42 == eGamma)
        {
            byDirs[0] = 1; byDirs[1] = 3;
        }
        else if (SIGMA43 == eGamma)
        {
            byDirs[0] = 2; byDirs[1] = 3;
        }
        for (INT idx = 0; idx < 2; ++idx)
        {
            UBOOL bPlus2 = (0 == (idx & 1));

            INT dimordered[2];
            for (INT order = 0; order < 2; ++order)
            {
                dimordered[0] = byDirs[0] + 1;
                dimordered[1] = bPlus2 ? static_cast<SBYTE>(byDirs[1] + 1) : static_cast<SBYTE>(-byDirs[1] - 1);
                if (1 == order)
                {
                    dimordered[0] = bPlus2 ? static_cast<SBYTE>(byDirs[1] + 1) : static_cast<SBYTE>(-byDirs[1] - 1);
                    dimordered[1] = byDirs[0] + 1;
                }
                checkCudaErrors(cudaMemcpy(devicePathBuffer, dimordered, sizeof(INT) * 2, cudaMemcpyHostToDevice));

                _kernelDFermionKSForce_WithLink_SigmaIJT << <block, threads >> > (
                    (const deviceGauge*)pGaugeBuffer,
                    (deviceGauge*)pForce,
                    pRationalFields,
                    pRationalNumerator,
                    uiRationalDegree,
                    appGetLattice()->m_pIndexCache->m_pEtaMu,
                    byFieldID,
                    byGaugeFieldID,
                    fCoeff * F(0.25),
                    byDirs[0],
                    byDirs[1],
                    devicePathBuffer
                    );
            }
        }
    }
    break;
    case GAMMA51:
    case GAMMA52:
    case GAMMA53:
    case GAMMA54:
    {
        const BYTE byMissingDir = static_cast<BYTE>(eGamma - GAMMA51);
        for (INT idx = 0; idx < 4; ++idx)
        {
            UBOOL bPlus123[3] = { (0 == (idx & 1)), (0 == (idx & 2)), TRUE };
            if (byMissingDir < 2)
            {
                bPlus123[0] = TRUE;
                bPlus123[1] = (0 == (idx & 1));
                bPlus123[2] = (0 == (idx & 2));
            }

            INT dim123[3];
            BYTE byDimIndex = 0;
            for (INT byCubeDir = 0; byCubeDir < 4; ++byCubeDir)
            {
                if (byCubeDir != static_cast<INT>(byMissingDir))
                {
                    dim123[byDimIndex] = bPlus123[byDimIndex] ? (byCubeDir + 1) : (-byCubeDir - 1);
                    byDimIndex++;
                }
            }

            INT dimordered[3];
            for (INT order = 0; order < 6; ++order)
            {
                switch (order)
                {
                case 1:
                    dimordered[0] = dim123[0];
                    dimordered[1] = dim123[2];
                    dimordered[2] = dim123[1];
                    break;
                case 2:
                    dimordered[0] = dim123[1];
                    dimordered[1] = dim123[0];
                    dimordered[2] = dim123[2];
                    break;
                case 3:
                    dimordered[0] = dim123[1];
                    dimordered[1] = dim123[2];
                    dimordered[2] = dim123[0];
                    break;
                case 4:
                    dimordered[0] = dim123[2];
                    dimordered[1] = dim123[0];
                    dimordered[2] = dim123[1];
                    break;
                case 5:
                    dimordered[0] = dim123[2];
                    dimordered[1] = dim123[1];
                    dimordered[2] = dim123[0];
                    break;
                default:
                    dimordered[0] = dim123[0];
                    dimordered[1] = dim123[1];
                    dimordered[2] = dim123[2];
                    break;
                }

                checkCudaErrors(cudaMemcpy(devicePathBuffer, dimordered, sizeof(INT) * 3, cudaMemcpyHostToDevice));

                _kernelDFermionKSForce_WithLink_Gamma51234T << <block, threads >> > (
                    (const deviceGauge*)pGaugeBuffer,
                    (deviceGauge*)pForce,
                    pRationalFields,
                    pRationalNumerator,
                    uiRationalDegree,
                    byFieldID,
                    byGaugeFieldID,
                    fCoeff * OneOver24,
                    byMissingDir,
                    devicePathBuffer
                    );
            }
        }
    }
    break;
    case GAMMA5:
    {
        for (INT idx = 0; idx < 8; ++idx)
        {
            const UBOOL bPlus1234[4] =
            {
                (0 == (idx & 1)),
                (0 == (idx & 2)),
                (0 == (idx & 4)),
                TRUE
            };

            const INT dim1234[4] =
            {
                bPlus1234[0] ? 1 : -1,
                bPlus1234[1] ? 2 : -2,
                bPlus1234[2] ? 3 : -3,
                4
            };

            INT dimordered[4];
            INT dim234[3];
            for (BYTE k = 0; k < 4; ++k)
            {
                dimordered[0] = dim1234[k];
                for (BYTE k2 = 0; k2 < 3; ++k2)
                {
                    BYTE idx2 = k2 + 1 + k;
                    idx2 = idx2 > 3 ? (idx2 - 4) : idx2;
                    dim234[k2] = dim1234[idx2];
                }

                for (BYTE order2 = 0; order2 < 6; ++order2)
                {
                    switch (order2)
                    {
                    case 1:
                        dimordered[1] = dim234[0];
                        dimordered[2] = dim234[2];
                        dimordered[3] = dim234[1];
                        break;
                    case 2:
                        dimordered[1] = dim234[1];
                        dimordered[2] = dim234[0];
                        dimordered[3] = dim234[2];
                        break;
                    case 3:
                        dimordered[1] = dim234[1];
                        dimordered[2] = dim234[2];
                        dimordered[3] = dim234[0];
                        break;
                    case 4:
                        dimordered[1] = dim234[2];
                        dimordered[2] = dim234[0];
                        dimordered[3] = dim234[1];
                        break;
                    case 5:
                        dimordered[1] = dim234[2];
                        dimordered[2] = dim234[1];
                        dimordered[3] = dim234[0];
                        break;
                    default:
                        dimordered[1] = dim234[0];
                        dimordered[2] = dim234[1];
                        dimordered[3] = dim234[2];
                        break;
                    }
                    //appGeneral(_T("dimordered=%d %d %d %d\n"), dimordered[0], dimordered[1], dimordered[2], dimordered[3]);
                    checkCudaErrors(cudaMemcpy(devicePathBuffer, dimordered, sizeof(INT) * 4, cudaMemcpyHostToDevice));
                    _kernelDFermionKSForce_WithLink_Gamma5T << <block, threads >> > (
                        (const deviceGauge*)pGaugeBuffer,
                        (deviceGauge*)pForce,
                        pRationalFields,
                        pRationalNumerator,
                        uiRationalDegree,
                        byFieldID,
                        byGaugeFieldID,
                        fCoeff * OneOver192,
                        devicePathBuffer
                        );
                }
            }
        }
    }
    break;
    default:
        appGeneral(_T("not implimented!\n"));
        break;

    }
}

#pragma endregion

#pragma region Matrix Operator

/**
* Assume m >= k
* V=(vk[0], ... , vk[k - 1])
* W=(vk[0], ... , vk[k - 1], vmk[0], ..., vmk[m-k-1])
*
* V(v1,v2,...,vk) = W(w1,w2,...,wm) (m11, ..., m1k)
*                                   (..., ..., ...)
*                                   (mm1, ..., mmk)
* I think this is expansive... the FLOP of Ax is about 100n, but this has m x k x n
*/
template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelMatrixMultiplyKST(
    INT vectorN,
    deviceVector** pRes,
    deviceVector** pLeft,
    const CLGComplex* __restrict__ pMatrix,
    UINT uiDimX, UINT uiDimY) //x=m,y=k
{
    intokernalE(vectorN);

    CLGComplex result[CFieldMatrixOperation::_kFieldMatrixMaxDim];

    for (UINT i = 0; i < uiDimY; ++i)
    {
        result[i] = _make_cuComplex(F(0.0), F(0.0));
        for (UINT j = 0; j < uiDimX; ++j)
        {
            const CLGComplex a = _make_cuComplex(_element(pRes[j][uiSiteIndex], 2 * elementIdx), _element(pRes[j][uiSiteIndex], 2 * elementIdx + 1));
            const CLGComplex b = _make_cuComplex(_element(pLeft[j - uiDimY][uiSiteIndex], 2 * elementIdx), _element(pLeft[j - uiDimY][uiSiteIndex], 2 * elementIdx + 1));
            result[i] = _cuCaddf(result[i], _cuCmulf(j < uiDimY ? a : b, pMatrix[j * uiDimY + i]));
        }
    }

    for (UINT i = 0; i < uiDimY; ++i)
    {
        _setelement(pRes[i][uiSiteIndex], 2 * elementIdx, result[i].x);
        _setelement(pRes[i][uiSiteIndex], 2 * elementIdx + 1, result[i].y);
    }
}

/**
* Assume m >= k
* V=(vk[0], ... , vk[k - 1])
* W=(vk[0], ... , vk[k - 1], vmk[0], ..., vmk[m-k-1])
*
* v1 = (m11, ..., m1m)  w1
* ..   (..., ..., ...)  ..
* vk   (mk1, ..., mkm)  wk
*                       wk+1
*                       ...
*                       wm
*/
template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::VectorMultiplyMatrix(deviceVector** hostResBuffer, deviceVector** hostLeftBuffer, deviceVector** resBuffer, deviceVector** leftBuffer,
    TArray<CField*>& res, const TArray<CField*>& left, const CLGComplex* deviceMatrix, UINT uiDimX, UINT uiDimY)
{
    for (UINT i = 0; i < uiDimY; ++i)
    {
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pF = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(res[i]);
        if (NULL == pF)
        {
            appCrucial(_T("CFieldMatrixOperationKSSU3 only work with CFieldFermionKST<deviceVector, deviceGauge, vectorN>!\n"));
            return;
        }
        hostResBuffer[i] = pF->m_pDeviceData;
    }

    for (UINT i = 0; i < uiDimX - uiDimY; ++i)
    {
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pF = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(left[i]);
        if (NULL == pF)
        {
            appCrucial(_T("CFieldMatrixOperationKSSU3 only work with CFieldFermionKST<deviceVector, deviceGauge, vectorN>!\n"));
            return;
        }
        hostLeftBuffer[i] = pF->m_pDeviceData;
    }

    checkCudaErrors(cudaMemcpy(resBuffer, hostResBuffer, sizeof(deviceVector*) * uiDimY, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(leftBuffer, hostLeftBuffer, sizeof(deviceVector*) * (uiDimX - uiDimY), cudaMemcpyHostToDevice));

    preparethreadE(vectorN);
    _kernelMatrixMultiplyKST << <block, threads >> > (vectorN, resBuffer, leftBuffer, deviceMatrix, uiDimX, uiDimY);
}

#pragma endregion

#pragma endregion

#pragma region DOperator Rotation

#pragma region kernel

/**
* When link n and n+mu, the coordinate is stick with n
* When link n and n-mu, the coordinate is stick with n-mu
* Irrelavent with tau
* Optimization: bXorY removed, block.x *= 2
*/
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTermT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    SSmallInt4 sCenter,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();
    //const INT eta_tau = ((pEtaTable[uiSiteIndex] >> 3) & 1);
    const INT eta_tau = pEtaTable[uiSiteIndex] >> 3;

#pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusMu = idx & 2;
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

        deviceVector right = _mulVec(_deviceVXXTauOptimizedT(pGauge, sSite4, byGaugeFieldId, bXorY, bPlusMu, bPlusTau),
            pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        //when bXorY = 1, it is y partial _x, so is [1]
        //when bXorY = 0, it is x partial _y, so is [0]
        _mul(right, sMidSite.m_byData4[bXorY] - sCenter.m_byData4[bXorY] + F(0.5));

        if (!bPlusMu)
        {
            //for -2x, -2y terms, there is another minus sign
            this_eta_tau = this_eta_tau + 1;
        }

        if (this_eta_tau & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, F(-0.25) * fOmega);
    }
    else
    {
        _mul(result, F(0.25) * fOmega);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    }

    _add(pResultData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTermRealT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    SSmallInt4 sCenter,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();
    //const INT eta_tau = ((pEtaTable[uiSiteIndex] >> 3) & 1);
    const INT eta_tau = pEtaTable[uiSiteIndex] >> 3;

#pragma unroll
    for (UINT idx = 0; idx < 8; ++idx)
    {
        const UBOOL bPlusMu = idx & 2;
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

        deviceVector right = _mulVec(_deviceVXXTauOptimizedT(pGauge, sSite4, byGaugeFieldId, bXorY, bPlusMu, bPlusTau),
            pDeviceData[sTargetBigIndex.m_uiSiteIndex]);

        //when bXorY = 1, it is y partial _x, so is [1]
        //when bXorY = 0, it is x partial _y, so is [0]
        _mul(right, sMidSite.m_byData4[bXorY] - sCenter.m_byData4[bXorY] + F(0.5));

        if (!bPlusMu)
        {
            //for -2x, -2y terms, there is another minus sign
            this_eta_tau = this_eta_tau + 1;
        }

        if (this_eta_tau & 1)
        {
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    //if (bDDagger)
    //{
    //    _mul(result, F(-0.25) * fOmega);
    //}
    //else
    //{
        //_mul(result, F(0.25) * fOmega);
    //}
    _mul(result, _make_cuComplex(F(0.0), F(-0.25) * fOmega));

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    }

    _add(pResultData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTau_TermT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

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

        const deviceVector right = _mulVec(_deviceVXYTOptimizedT(pGauge, sSite4, byGaugeFieldId, bPlusX, bPlusY, bPlusT), pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
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
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    if (bDDagger)
    {
        _mul(result, -F(0.125) * fOmega);
    }
    else
    {
        _mul(result, F(0.125) * fOmega);
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    }

    _add(pResultData[uiSiteIndex], result);
}

template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_PR_XYTau_TermRealT(
    const deviceVector* __restrict__ pDeviceData,
    const deviceGauge* __restrict__ pGauge,
    deviceVector* pResultData,
    BYTE byFieldId,
    BYTE byGaugeFieldId,
    DOUBLE fOmega,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernalInt4;

    deviceVector result = _makeZero<deviceVector>();

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

        const deviceVector right = _mulVec(_deviceVXYTOptimizedT(pGauge, sSite4, byGaugeFieldId, bPlusX, bPlusY, bPlusT), pDeviceData[sTargetBigIndex.m_uiSiteIndex]);
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
            _sub(result, right);
        }
        else
        {
            _add(result, right);
        }
    }

    //if (bDDagger)
    //{
    //    _mul(result, -F(0.125) * fOmega);
    //}
    //else
    //{
    //    _mul(result, F(0.125) * fOmega);
    //}
    _mul(result, _make_cuComplex(F(0.0), F(-0.125) * fOmega));

    switch (eCoeff)
    {
    case EOCT_Real:
        _mul(result, fCoeff);
        break;
    case EOCT_Complex:
        _mul(result, cCoeff);
        break;
    }

    _add(pResultData[uiSiteIndex], result);
}

#pragma endregion

#pragma region Derivate

/**
 * Have n, n->n1, n->n2,
 * 1. we need to obtain V_(n, n1) , V_(n, n2)
 * 2. we need phi(n1), phi(n2), phid(n1), phid(n2)
 *
 * byContribution: 0 for mu, 1 for tau, 2 for both mu and tau
 *
 * iTau = 1 for +t, -1 for -t
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PR_XYTermT(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const BYTE* __restrict__ pEtaTable,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    DOUBLE fOmega,
    BYTE byMu, INT iTau,
    INT pathLdir1, INT pathLdir2, INT pathLdir3, BYTE Llength,
    INT pathRdir1, INT pathRdir2, INT pathRdir3, BYTE Rlength,
    BYTE byContribution)
{
    intokernalInt4;
    //const UINT uiBigIdx = __bi(sSite4);

    //=================================
    // 1. Find n1, n2
    INT Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    INT Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    SSmallInt4 site_n1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(_deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength))];
    //const SSmallInt4 middleSite = _deviceSmallInt4OffsetC(site_n1, byMu + 1);
    //From now on, site_n1 is smiddle
    site_n1 = _deviceSmallInt4OffsetC(site_n1, byMu + 1);
    const SIndex& smiddle = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(site_n1)];

    site_n1 = __deviceSiteIndexToInt4(smiddle.m_uiSiteIndex);
    //y Dx and -x Dy
    const Real fNv = (0 == byMu)
        ? static_cast<Real>(site_n1.y - _DC_Centery + F(0.5))
        : static_cast<Real>(_DC_Centerx - site_n1.x - F(0.5));

    //=================================
    // 2. Find V(n,n1), V(n,n2)
    const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, Llength, 1, Ldirs);
    const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, Rlength, 1, Rdirs);

    for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
    {
        const deviceVector* phi_i = pFermionPointers[rfieldId];
        const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];
        //=================================
        // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
        deviceVector phi1 = _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]);
        deviceVector phi2 = _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]);
        deviceVector phi3 = _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]);
        deviceVector phi4 = _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]);
        if (sn1.NeedToOpposite())
        {
            _mul(phi1, F(-1.0));
            _mul(phi3, F(-1.0));
        }
        if (sn2.NeedToOpposite())
        {
            _mul(phi2, F(-1.0));
            _mul(phi4, F(-1.0));
        }
        deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
        _add(res, _makeContract<deviceGauge, deviceVector>(phi4, phi3));
        _ta(res);
        const Real eta_tau = (iTau > 0 ?
            ((pEtaTable[sn1.m_uiSiteIndex] >> 3) & 1)
            : ((pEtaTable[sn2.m_uiSiteIndex] >> 3) & 1))
            ? F(-1.0) : F(1.0);
        _mul(res, OneOver12 * fOmega * fNv * pNumerators[rfieldId] * eta_tau);

        //For mu
        if (0 == byContribution || 2 == byContribution)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byMu);
            _sub(pForce[linkIndex], res);
        }

        //For tau
        if (1 == byContribution || 2 == byContribution)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, 3);
            if (iTau > 0)
            {
                _sub(pForce[linkIndex], res);
            }
            else
            {
                _add(pForce[linkIndex], res);
            }
        }
    }
}

/**
 *
 */
template<typename deviceVector, typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_PR_XYTau_TermT(
    const deviceGauge* __restrict__ pGauge,
    deviceGauge* pForce,
    const deviceVector* const* __restrict__ pFermionPointers,
    const Real* __restrict__ pNumerators,
    UINT uiRational,
    BYTE byFieldId,
    DOUBLE fOmega,
    INT pathLdir1, INT pathLdir2, INT pathLdir3, BYTE Llength,
    INT pathRdir1, INT pathRdir2, INT pathRdir3, BYTE Rlength)
{
    intokernalInt4;
    //const UINT uiBigIdx = __bi(sSite4);

    //=================================
    // 1. Find n1, n2
    INT Ldirs[3] = { pathLdir1, pathLdir2, pathLdir3 };
    INT Rdirs[3] = { pathRdir1, pathRdir2, pathRdir3 };
    const SSmallInt4 siten1 = _deviceSmallInt4OffsetC(sSite4, Ldirs, Llength);
    const SSmallInt4 siten2 = _deviceSmallInt4OffsetC(sSite4, Rdirs, Rlength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten1)];
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten2)];

    //Why use sn2? shouldn't it be sn1?
    const Real eta124 = _deviceEta124(__deviceSiteIndexToInt4(sn1.m_uiSiteIndex));
    //=================================
    // 2. Find V(n,n1), V(n,n2)
    const deviceGauge vnn1 = _deviceLinkT(pGauge, sSite4, Llength, 1, Ldirs);
    const deviceGauge vnn2 = _deviceLinkT(pGauge, sSite4, Rlength, 1, Rdirs);

    for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
    {
        const deviceVector* phi_i = pFermionPointers[rfieldId];
        const deviceVector* phi_id = pFermionPointers[rfieldId + uiRational];

        //=================================
        // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
        deviceVector phi1 = _mulVec(vnn1, phi_id[sn1.m_uiSiteIndex]);
        deviceVector phi2 = _mulVec(vnn2, phi_i[sn2.m_uiSiteIndex]);
        deviceVector phi3 = _mulVec(vnn1, phi_i[sn1.m_uiSiteIndex]);
        deviceVector phi4 = _mulVec(vnn2, phi_id[sn2.m_uiSiteIndex]);
        if (sn1.NeedToOpposite())
        {
            _mul(phi1, F(-1.0));
            _mul(phi3, F(-1.0));
        }
        if (sn2.NeedToOpposite())
        {
            _mul(phi2, F(-1.0));
            _mul(phi4, F(-1.0));
        }
        deviceGauge res = _makeContract<deviceGauge, deviceVector>(phi1, phi2);
        //This was phi2 phi1+ * eta124(n1) - phi3 phi4+ * eta124(n2)
        //The sign of the second term is because of 'dagger'
        //However, eta124(n1) = -eta124(n2), so use Add directly.
        _add(res, _makeContract<deviceGauge, deviceVector>(phi4, phi3));
        _ta(res);
        _mul(res, OneOver48 * static_cast<Real>(fOmega) * pNumerators[rfieldId] * eta124);

        //Use eta124 of n2 so Add left Sub right
        //Change to use eta124 of n1, Sub left and Add right
        if (pathLdir1 > 0)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLdir1 - 1);
            _add(pForce[linkIndex], res);
        }

        if (pathRdir1 > 0)
        {
            const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRdir1 - 1);
            _sub(pForce[linkIndex], res);
        }
    }

}

#pragma endregion

#pragma region D and derivate

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKS_R_RealRotation(
    DOUBLE fOmega, 
    UBOOL bEachSiteEta, deviceVector* pTarget, const deviceVector* pSource,
    const deviceGauge* pGauge, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    //CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DOperatorKS(pTargetBuffer, pBuffer, pGaugeBuffer, byGaugeFieldId, f2am, bDagger, eOCT, fRealCoeff, cCmpCoeff);

    preparethread;
    if (bDagger)
    {
        appCrucial(_T("D dagger is not supported for real rotation!\n"));
    }

    _kernelDFermionKS_PR_XYTermRealT << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        _HC_Center,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    _kernelDFermionKS_PR_XYTau_TermRealT << <block, threads >> > (
        pSource,
        pGauge,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DOperatorKS_R_ImaginaryRotation(
    DOUBLE fOmega,
    UBOOL bEachSiteEta, deviceVector* pTarget, const deviceVector* pSource,
    const deviceGauge* pGauge, BYTE byFieldId, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    _kernelDFermionKS_PR_XYTermT << <block, threads >> > (
        pSource,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        _HC_Center,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);

    _kernelDFermionKS_PR_XYTau_TermT << <block, threads >> > (
        pSource,
        pGauge,
        pTarget,
        byFieldId,
        byGaugeFieldId,
        fOmega,
        bDagger,
        eOCT,
        fRealCoeff,
        cCmpCoeff);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKSTKernel<deviceVector, deviceGauge, vectorN>::DerivateD0_R(
    DOUBLE fOmega,
    const deviceVector* pFermion,
    BYTE byFieldId,
    deviceGauge* pForce,
    const deviceGauge* pGaugeBuffer,
    BYTE byGaugeFieldId,
    const deviceVector* const* pRationalFields,
    const Real* pNumerator,
    UINT uiRationApproxOrder)
{


    preparethread;

#pragma region X Y Term
    INT mu[2] = { 0, 1 };
    for (INT imu = 0; imu < 2; ++imu)
    {
        INT dirs[6][3] =
        {
            {4, mu[imu] + 1, mu[imu] + 1},
            {mu[imu] + 1, 4, mu[imu] + 1},
            {mu[imu] + 1, mu[imu] + 1, 4},
            //{4, -mu[imu] - 1, -mu[imu] - 1},
            //{-mu[imu] - 1, 4, -mu[imu] - 1},
            //{-mu[imu] - 1, -mu[imu] - 1, 4},
            {mu[imu] + 1, mu[imu] + 1, -4},
            {mu[imu] + 1, -4, mu[imu] + 1},
            {-4, mu[imu] + 1, mu[imu] + 1},
        };

        INT iTau[6] = { 1, 1, 1, -1, -1, -1 };
        BYTE contributionOf[6][4] =
        {
            {1, 0, 0, 3},
            {0, 1, 0, 3},
            {0, 0, 1, 3},
            //{1, 3, 0, 0},
            //{3, 2, 3, 0},
            //{3, 0, 2, 3},
            {0, 0, 3, 1},
            {0, 3, 2, 3},
            {3, 2, 0, 3},
        };

        for (INT pathidx = 0; pathidx < 6; ++pathidx)
        {
            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                if (3 == contributionOf[pathidx][iSeperation])
                {
                    continue;
                }

                INT L[3] = { 0, 0, 0 };
                INT R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(dirs[pathidx], iSeperation, L, R, LLength, RLength);

                _kernelDFermionKSForce_PR_XYTermT << <block, threads >> > (
                    pGaugeBuffer,
                    pForce,
                    appGetLattice()->m_pIndexCache->m_pEtaMu,
                    pRationalFields,
                    pNumerator,
                    uiRationApproxOrder,
                    byFieldId,
                    fOmega,
                    static_cast<BYTE>(imu), iTau[pathidx],
                    L[0], L[1], L[2], LLength,
                    R[0], R[1], R[2], RLength,
                    contributionOf[pathidx][iSeperation]
                    );
            }
        }
    }

#pragma endregion

#pragma region Polarization term

    //===========================
    //polarization terms
    //ilinkType is +-x +-y +t,
    //INT linkTypes[4][3] =
    //{
    //    {1, 2, 4},
    //    {1, 2, -4},
    //    {-1, 2, 4},
    //    {-1, 2, -4}
    //};
    INT linkTypes[4][3] =
    {
        {1, 2, 4},
        {1, -2, 4},
        {-1, 2, 4},
        {-1, -2, 4}
    };

    for (INT ilinkType = 0; ilinkType < 4; ++ilinkType)
    {
        INT sixlinks[6][3] =
        {
            {linkTypes[ilinkType][0], linkTypes[ilinkType][1], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][0], linkTypes[ilinkType][2], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][0], linkTypes[ilinkType][2]},
            {linkTypes[ilinkType][1], linkTypes[ilinkType][2], linkTypes[ilinkType][0]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][0], linkTypes[ilinkType][1]},
            {linkTypes[ilinkType][2], linkTypes[ilinkType][1], linkTypes[ilinkType][0]}
        };

        for (INT isixtype = 0; isixtype < 6; ++isixtype)
        {
            //bearly no change of time, because force calculation is not frequent
            /*
            _giveupkernelDFermionKSForce_PR_XYTau_Term2 << <block, threads >> > (
                (const deviceGauge*)pGaugeBuffer,
                (deviceGauge*)pForce,
                m_pRationalFieldPointers,
                m_pMDNumerator,
                m_rMD.m_uiDegree,
                m_byFieldId,
                CCommonData::m_fOmega,
                sixlinks[isixtype][0], sixlinks[isixtype][1], sixlinks[isixtype][2]
                );
            */

            for (INT iSeperation = 0; iSeperation < 4; ++iSeperation)
            {
                INT L[3] = { 0, 0, 0 };
                INT R[3] = { 0, 0, 0 };
                BYTE LLength = 0;
                BYTE RLength = 0;

                Seperate(sixlinks[isixtype], iSeperation, L, R, LLength, RLength);

                const UBOOL bHasLeft = (LLength > 0) && (L[0] > 0);
                const UBOOL bHasRight = (RLength > 0) && (R[0] > 0);

                if (bHasLeft || bHasRight)
                {
                    _kernelDFermionKSForce_PR_XYTau_TermT << <block, threads >> > (
                        pGaugeBuffer,
                        pForce,
                        pRationalFields,
                        pNumerator,
                        uiRationApproxOrder,
                        byFieldId,
                        fOmega,
                        L[0], L[1], L[2], LLength,
                        R[0], R[1], R[2], RLength
                        );
                }
            }
        }
    }

#pragma endregion
}

#pragma endregion

#pragma endregion


template class CFieldFermionKSTKernel<CLGComplex, CLGComplex, 1>;
template class CFieldFermionKSTKernel<deviceSU2Vector, deviceSU2, 2>;
template class CFieldFermionKSTKernel<deviceSU4Vector, deviceSU4, 4>;
template class CFieldFermionKSTKernel<deviceSU5Vector, deviceSU5, 5>;
template class CFieldFermionKSTKernel<deviceSU6Vector, deviceSU6, 6>;
template class CFieldFermionKSTKernel<deviceSU7Vector, deviceSU7, 7>;
template class CFieldFermionKSTKernel<deviceSU8Vector, deviceSU8, 8>;

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================