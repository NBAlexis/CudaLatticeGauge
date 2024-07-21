//=============================================================================
// FILENAME : CFieldFermionKST<deviceVector, deviceGauge, vectorN>.cu
// 
// DESCRIPTION:
// This is the device implementations of Wilson fermion
//
// This implementation assumes SU3 and square lattice
//
// REVISION:
//  [12/08/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "Data/Field/Gauge/CFieldGaugeLink.h"
#include "CFieldFermionKST.h"
//#include "CFieldFermionKSSU3Gamma.h"

__BEGIN_NAMESPACE

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
    const deviceVector * __restrict__ pDeviceData,
    const deviceGauge * __restrict__ pGauge,
    const SIndex * __restrict__ pGaugeMove,
    const SIndex * __restrict__ pFermionMove,
    const BYTE * __restrict__ pEtaTable,
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
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DOperatorKS(void* pTargetBuffer, const void * pBuffer,
    const void * pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceVector* pTarget = (deviceVector*)pTargetBuffer;
    const deviceVector* pSource = (const deviceVector*)pBuffer;
    const deviceGauge* pGauge = (const deviceGauge*)pGaugeBuffer;

    preparethread;
    if (m_bEachSiteEta)
    {
        _kernelDFermionKSPlusEtaT << <block, threads >> > (
            pSource,
            pGauge,
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            pTarget,
            f2am,
            m_byFieldId,
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
            appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
            appGetLattice()->m_pIndexCache->m_pEtaMu,
            pTarget,
            f2am,
            m_byFieldId,
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
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DerivateD0(
    void* pForce, 
    const void* pGaugeBuffer,
    BYTE byGaugeFieldId) const
{
    preparethread;
    _kernelDFermionKSForceT << <block, threads >> > (
        (const deviceGauge*)pGaugeBuffer,
        (deviceGauge*)pForce,
        appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        m_pRationalFieldPointers,
        m_pMDNumerator,
        m_rMD.m_uiDegree,
        m_byFieldId);
}

#pragma endregion

#pragma region Staggered

/**
* Initialize
*/
template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialFermionKST(deviceVector* pDevicePtr, BYTE byFieldId, EFieldInitialType eInitialType)
{
    intokernalInt4;

    switch (eInitialType)
    {
    case EFIT_Zero:
    {
        pDevicePtr[uiSiteIndex] = _makeZero<deviceVector>();
    }
    break;
    case EFIT_Identity:
    {
        pDevicePtr[uiSiteIndex] = _makeId<deviceVector>();
    }
    break;
    case EFIT_RandomGaussian:
    {
        pDevicePtr[uiSiteIndex] = _makeGaussian<deviceVector>(_deviceGetFatIndex(uiSiteIndex, 0));
    }
    break;
    case EFIT_RandomZ4:
    {
        pDevicePtr[uiSiteIndex] = _makeZ4<deviceVector>(_deviceGetFatIndex(uiSiteIndex, 0));
    }
    break;
    default:
    {
        printf("Wilson Fermion Field cannot be initialized with this type!");
    }
    break;
    }
}

/**
* generate phi by gaussian random.
* phi = (D^+D)^{1/4} phi
*/
template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson)
{
    preparethread;
    _kernelInitialFermionKST << <block, threads >> > (
        m_pDeviceData,
        m_byFieldId,
        EFIT_RandomGaussian);

    D_MC(gaugeNum, bosonNum, gaugeFields, pBoson);

    if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
    {
        m_fLength = Dot(this).x;
    }

    //For KS, we generally use shifted solver, so do NOT cache the last result
}

/**
 * Use \sqrt{a} and b of rational A^{-1/2}
 */
template<typename deviceVector, typename deviceGauge, INT vectorN>
UBOOL CFieldFermionKST<deviceVector, deviceGauge, vectorN>::CalculateForceS(
    const CFieldGauge* pGauge,
    CFieldGauge* pForce,
    ESolverPhase ePhase) const
{
    if (NULL == pGauge || vectorN != pGauge->MatrixN())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
        return FALSE;
    }
    if (NULL == pForce || vectorN != pForce->MatrixN())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
        return FALSE;
    }

    TArray<CField*> phii;
    TArray<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*> phiid;
    for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
    {
        CField* pPhi_i = dynamic_cast<CField*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
        phii.AddItem(pPhi_i);
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPhi_id = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
        phiid.AddItem(pPhi_id);
    }

    CMultiShiftSolver* solver = appGetMultiShiftSolver(m_byFieldId);
    if (NULL == solver)
    {
        appCrucial(_T("No multi solver found!"));
        _FAIL_EXIT;
    }

    TArray<CLGComplex> shifts;
    for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
    {
        shifts.AddItem(_make_cuComplex(m_rMD.m_lstB[i], F(0.0)));
    }
    TArray<const CFieldGauge*> gaues;
    gaues.AddItem(pGauge);
    solver->Solve(phii, shifts, this, 1, 0, gaues.GetData(), NULL, EFO_F_DDdagger);

    const UINT uiBufferSize = sizeof(deviceVector*) * 2 * m_rMD.m_uiDegree;
    deviceVector** hostPointers = (deviceVector**)appAlloca(uiBufferSize);
    for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
    {
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* phi_ks = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(phii[i]);
        phi_ks->CopyTo(phiid[i]);
        if (m_bDiagonalMass)
        {
            phiid[i]->DS(pGauge);
        }
        else
        {
            phiid[i]->D0S(pGauge);
        }

        hostPointers[i] = phi_ks->m_pDeviceData;
        hostPointers[i + m_rMD.m_uiDegree] = phiid[i]->m_pDeviceData;
    }
    checkCudaErrors(cudaMemcpy(m_pRationalFieldPointers, hostPointers, uiBufferSize, cudaMemcpyHostToDevice));

    const CFieldGaugeLink<deviceGauge, vectorN>* pGaugeSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
    CFieldGaugeLink<deviceGauge, vectorN>* pForceSU3 = dynamic_cast<CFieldGaugeLink<deviceGauge, vectorN>*>(pForce);

    DerivateD0(pForceSU3->m_pDeviceData, pGaugeSU3->m_pDeviceData, pGaugeSU3->m_byFieldId);
    
    for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
    {
        phii[i]->Return();
        phiid[i]->Return();
    }

    return TRUE;
}

#pragma endregion

#pragma region Kernel

template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusFermionKST(deviceVector* pMe, const deviceVector* __restrict__ pOther)
{
    intokernal;
    _add(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusFermionKST(deviceVector* pMe, const deviceVector* __restrict__ pOther)
{
    intokernal;
    _sub(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyComplexFermionKST(deviceVector* pMe, const deviceVector* __restrict__ pOther, CLGComplex a)
{
    intokernal;
    _add(pMe[uiSiteIndex], _mulC(pOther[uiSiteIndex], a));
}

template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelMulFermionKST(deviceVector* pMe, const deviceVector* __restrict__ pOther, UBOOL bConj)
{
    intokernal;
    if (bConj)
    {
        _dagger(pMe[uiSiteIndex]);
    }
    _mul(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyRealFermionKST(deviceVector* pMe, const deviceVector* __restrict__ pOther, Real a)
{
    intokernal;
    _add(pMe[uiSiteIndex], _mulC(pOther[uiSiteIndex], a));
}

template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelDotFermionKST(const deviceVector* __restrict__ pMe, const deviceVector* __restrict__ pOther, cuDoubleComplex* result
)
{
    intokernal;
    result[uiSiteIndex] = _cToDouble(_dot(pMe[uiSiteIndex], pOther[uiSiteIndex]));
}

template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyComplexKST(deviceVector* pMe, CLGComplex a)
{
    intokernal;
    _mul(pMe[uiSiteIndex], a);
}

template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyRealKST(deviceVector* pMe, Real a)
{
    intokernal;
    _mul(pMe[uiSiteIndex], a);
}

template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelFermionKSConjugateT(deviceVector* pDeviceData)
{
    intokernal;
    _dagger(pDeviceData[uiSiteIndex]);
}


template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelMakePointSourceKST(deviceVector* pDeviceData, UINT uiDesiredSite, BYTE byColor)
{
    intokernal;
    if (uiSiteIndex == uiDesiredSite)
    {
        pDeviceData[uiSiteIndex] = _makeColorVector<deviceVector>(byColor);
    }
    else
    {
        pDeviceData[uiSiteIndex] = _makeZero<deviceVector>();
    }
}

template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelMakePointSourceKSOneT(deviceVector* pDeviceData, UINT uiDesiredSite)
{
    intokernal;
    if (uiSiteIndex == uiDesiredSite)
    {
        pDeviceData[uiSiteIndex] = _makeId<deviceVector>();
    }
    else
    {
        pDeviceData[uiSiteIndex] = _makeZero<deviceVector>();
    }
}

template<typename deviceVector>
__global__ void _CLG_LAUNCH_BOUND
_kernelMakeWallSourceKST(deviceVector* pDeviceData, 
    INT uiDesiredT, UINT uiShift, BYTE color, BYTE byFieldID)
{
    intokernalOnlyInt4;

    //pDeviceData[uiSiteIndex] = _makeZero<deviceVector>_makeZero();
    //We shall not set zero here! [2024/7/21/ why not: because offset site will be set to non-zero!!]

    if ( (0 == (sSite4.x & 1))
      && (0 == (sSite4.y & 1))
      && (0 == (sSite4.z & 1))
      && (uiDesiredT < 0 || uiDesiredT == sSite4.w))
    {
        //sSite4 is no longer used
        sSite4.x = sSite4.x + static_cast<SBYTE>(uiShift & 1);
        sSite4.y = sSite4.y + static_cast<SBYTE>((uiShift >> 1) & 1);
        sSite4.z = sSite4.z + static_cast<SBYTE>((uiShift >> 2) & 1);
        const SIndex& sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldID][__bi(sSite4)];
        if (!sIdx.IsDirichlet())
        {
            pDeviceData[sIdx.m_uiSiteIndex] = _makeColorVector<deviceVector>(color);
        }
    }
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
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::OnlyMass(void* pTarget, Real fm, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    _kernelDFermionKS_OnlyMassT << <block, threads >> > (
        m_pDeviceData,
        (deviceVector*)pTarget,
        fm,
        eOCT,
        fRealCoeff,
        cCmpCoeff
        );
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::OneLinkS(
    const void* pGuage, 
    BYTE byGaugeFieldId, 
    void* pTarget,
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
        m_pDeviceData,
        (const deviceGauge*)pGuage,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        (deviceVector*)pTarget,
        m_byFieldId,
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
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::OneLinkForceS(
    const void* pGuage, 
    BYTE byGaugeFieldId, 
    void* pForce,
    Real fCoefficient,
    const INT* pDevicePath, 
    BYTE pathLength, 
    BYTE byEtaIdx) const
{
    assert(pathLength <= _kLinkMaxLength);
    preparethread;
    _kernelDFermionKSForce_WithLinkT << <block, threads >> > (
        (const deviceGauge*)pGuage,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        (deviceGauge*)pForce,
        m_pRationalFieldPointers,
        m_pMDNumerator,
        m_rMD.m_uiDegree,
        m_byFieldId,
        byGaugeFieldId,
        fCoefficient,
        byEtaIdx,
        pDevicePath,
        pathLength
        );
}

#pragma endregion

template<typename deviceVector, typename deviceGauge, INT vectorN>
CFieldFermionKST<deviceVector, deviceGauge, vectorN>::CFieldFermionKST()
    : CFieldFermionKS()
    , m_pRationalFieldPointers(NULL)
{
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount));
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
CFieldFermionKST<deviceVector, deviceGauge, vectorN>::~CFieldFermionKST()
{
    checkCudaErrors(__cudaFree(m_pDeviceData));
    if (NULL != m_pRationalFieldPointers)
    {
        checkCudaErrors(cudaFree(m_pRationalFieldPointers));
    }
    if (NULL != m_pMDNumerator)
    {
        checkCudaErrors(cudaFree(m_pMDNumerator));
    }
}

/**
*
*/
template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialField(EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialFermionKST << <block, threads >> > (m_pDeviceData, m_byFieldId, eInitialType);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::Dagger()
{
    preparethread;
    _kernelFermionKSConjugateT << <block, threads >> > (m_pDeviceData);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFieldType)
{
    if (eFieldType != EFFT_CLGBin)
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialFieldWithFile: Only support CLG Bin File\n"));
        return;
    }

    UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * vectorN * m_uiSiteCount);
    BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
    if (NULL == data)
    {
        appCrucial(_T("File not found: %s\n"), sFileName.c_str());
        _FAIL_EXIT;
    }
    if (uiSize != static_cast<UINT>(sizeof(Real) * 2 * vectorN * m_uiSiteCount))
    {
        appCrucial(_T("File size not correct: expecting: %d, found: %d\n"), static_cast<UINT>(sizeof(Real) * 2 * vectorN * m_uiSiteCount), uiSize);
        _FAIL_EXIT;
    }
    InitialWithByte(data);
    free(data);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialWithByte(BYTE* byData)
{
    deviceVector* readData = (deviceVector*)malloc(sizeof(deviceVector) * m_uiSiteCount);
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        Real thisSite[2 * vectorN];
        memcpy(thisSite, byData + i * sizeof(Real) * 2 * vectorN, sizeof(Real) * 2 * vectorN);
        for (UINT k = 0; k < 2 * vectorN; ++k)
        {
            _setelement(readData[i], k, thisSite[k]);
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyHostToDevice));
    free(readData);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialOtherParameters(CParameters& params)
{
    CFieldFermionKS::InitialOtherParameters(params);

    if (NULL != m_pRationalFieldPointers)
    {
        checkCudaErrors(cudaFree(m_pRationalFieldPointers));
    }
    checkCudaErrors(cudaMalloc((void**)&m_pRationalFieldPointers, sizeof(deviceVector*) * 2 * m_rMD.m_uiDegree));
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DebugPrintMe() const
{
    deviceVector* toprint = (deviceVector*)malloc(sizeof(deviceVector) * m_uiSiteCount);
    checkCudaErrors(cudaMemcpy(toprint, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT uiSite = 0; uiSite < m_uiSiteCount; ++uiSite)
    {
        const SSmallInt4 site4 = __hostSiteIndexToInt4(uiSite);
        appGeneral(_T(" --- %d,%d,%d,%d --- %s\n"),
            site4.x, site4.y, site4.z, site4.w, appToString(toprint[uiSite]).c_str()
        );
    }

    appSafeFree(toprint);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::CopyTo(CField* U) const
{
    if (NULL == U || GetFieldType() != U->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
        return;
    }

    CFieldFermionKS::CopyTo(U);

    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(U);
    checkCudaErrors(cudaMemcpy(pField->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    if (NULL != pField->m_pMDNumerator)
    {
        checkCudaErrors(cudaFree(pField->m_pMDNumerator));
    }
    if (NULL != pField->m_pRationalFieldPointers)
    {
        checkCudaErrors(cudaFree(pField->m_pRationalFieldPointers));
    }

    checkCudaErrors(cudaMalloc((void**)&pField->m_pRationalFieldPointers, sizeof(deviceVector*) * 2 * m_rMD.m_uiDegree));
    checkCudaErrors(cudaMalloc((void**)&pField->m_pMDNumerator, sizeof(Real) * m_rMD.m_uiDegree));
    checkCudaErrors(cudaMemcpy(pField->m_pMDNumerator, m_pMDNumerator, sizeof(Real) * m_rMD.m_uiDegree, cudaMemcpyDeviceToDevice));
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::AxpyPlus(const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
        return;
    }
    const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);

    preparethread;
    _kernelAxpyPlusFermionKST << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::AxpyMinus(const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
        return;
    }
    const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);

    preparethread;
    _kernelAxpyMinusFermionKST << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::Axpy(Real a, const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
        return;
    }
    const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);

    preparethread;
    _kernelAxpyRealFermionKST << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
        return;
    }
    const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);

    preparethread;
    _kernelAxpyComplexFermionKST << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::Mul(const CField* other, UBOOL bDagger)
{
    if (NULL == other || GetFieldType() != other->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
        return;
    }
    const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(other);

    preparethread;
    _kernelMulFermionKST << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, bDagger);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
cuDoubleComplex CFieldFermionKST<deviceVector, deviceGauge, vectorN>::Dot(const CField* x) const
{
    if (NULL == x || GetFieldType() != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only copy to CFieldFermionKST<deviceVector, deviceGauge, vectorN>!"));
        return make_cuDoubleComplex(0, 0);
    }
    const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pField = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(x);
    preparethread;
    _kernelDotFermionKST << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, _D_ComplexThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::ScalarMultply(const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplyComplexKST << <block, threads >> > (m_pDeviceData, a);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplyRealKST << <block, threads >> > (m_pDeviceData, a);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::ApplyGamma(EGammaMatrix eGamma)
{
    appCrucial(_T("Not implemented yet...\n"));
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::ApplyGammaKSS(const CFieldGauge* pGauge, EGammaMatrix eGamma)
{
    if (NULL == pGauge || vectorN != pGauge->MatrixN())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
    InitialField(EFIT_Zero);

    //If it was gamma_mu or gamma_5 or sigmaij, it is i gamma mu and i gamma 5, therefore multiply -i
    UBOOL bImag = (GAMMA1 == eGamma) || (GAMMA2 == eGamma) || (GAMMA3 == eGamma) || (GAMMA4 == eGamma) || (GAMMA5 == eGamma)
        || (SIGMA12 == eGamma) || (SIGMA31 == eGamma) || (SIGMA41 == eGamma) || (SIGMA23 == eGamma) || (SIGMA42 == eGamma) || (SIGMA43 == eGamma);

    appApplyGammaKS(
        m_pDeviceData,
        pPooled->m_pDeviceData,
        pFieldSU3->m_pDeviceData,
        eGamma,
        m_bEachSiteEta,
        FALSE,
        F(0.5),
        bImag ? EOCT_Complex : EOCT_None,
        F(1.0),
        bImag ? _make_cuComplex(F(0.0), -F(1.0)) : _onec,
        m_byFieldId,
        pGauge->m_byFieldId
    );

    pPooled->Return();
}


//Kai should be part of D operator
template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
        return;
    }
    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
        return;
    }
    if (m_bDiagonalMass)
    {
        appCrucial(_T("In the cass mass is not a number, should not in here!\n"));
    }

    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::D0S(const CField* pGauge)
{
    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
        return;
    }
    if (m_bDiagonalMass)
    {
        appCrucial(_T("In the cass mass is not a number, should not in here except for check anti-hermiticity!\n"));
    }

    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, F(0.0),
        FALSE, EOCT_None, F(1.0), _onec);

    pPooled->Return();
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
UINT CFieldFermionKST<deviceVector, deviceGauge, vectorN>::TestAntiHermitianS(const CFieldGauge* pGauge) const
{
    const UINT uiVolume = _HC_Volume;
    const UINT uiRealVolume = vectorN * uiVolume;
    CLGComplex* matrixElement = (CLGComplex*)malloc(sizeof(CLGComplex) * uiRealVolume * uiRealVolume);
    deviceVector* hostData = (deviceVector*)malloc(sizeof(deviceVector) * uiVolume);
    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* v = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

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

//Kai should be part of D operator
template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
        return;
    }
    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
        TRUE, eCoeffType, fRealCoeff, cCompCoeff);


    pPooled->Return();
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DdaggerWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
        return;
    }
    if (m_bDiagonalMass)
    {
        appCrucial(_T("In the cass mass is not a number, should not in here!\n"));
    }

    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
        TRUE, eCoeffType, fRealCoeff, cCompCoeff);


    pPooled->Return();
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DDdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
        return;
    }

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }
    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
        TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //why only apply coeff in the next step?
    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DDS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
        return;
    }

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }
    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
        FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //why only apply coeff in the next step?
    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DDdaggerWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
        return;
    }
    if (m_bDiagonalMass)
    {
        appCrucial(_T("In the cass mass is not a number, should not in here!\n"));
    }

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }
    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    DOperatorKS(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
        TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //why only apply coeff in the next step?
    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::DDWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    const CFieldGaugeLink<deviceGauge, vectorN>* pFieldSU3 = dynamic_cast<const CFieldGaugeLink<deviceGauge, vectorN>*>(pGauge);
    if (NULL == pFieldSU3 || vectorN != pFieldSU3->MatrixN())
    {
        appCrucial(_T("CFieldFermionKST<deviceVector, deviceGauge, vectorN> can only play with gauge SU3!"));
        return;
    }

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }
    CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pPooled = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    DOperatorKS(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
        FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //why only apply coeff in the next step?
    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, pFieldSU3->m_byFieldId, fMass,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::InitialAsSource(const SFermionBosonSource& sourceData)
{
    const UINT uiSiteIndex = _hostGetSiteIndex(sourceData.m_sSourcePoint);
    switch (sourceData.m_eSourceType)
    {
    case EFS_Point:
    {
        preparethread;
        if (sourceData.m_byColorIndex < vectorN)
        {
            _kernelMakePointSourceKST << <block, threads >> > (m_pDeviceData, uiSiteIndex, sourceData.m_byColorIndex);
        }
        else
        {
            _kernelMakePointSourceKSOneT << <block, threads >> > (m_pDeviceData, uiSiteIndex);
        }
    }
    break;
    case EFS_Wall:
    {
        preparethread;
        _kernelInitialFermionKST << <block, threads >> > (m_pDeviceData, m_byFieldId,EFIT_Zero);
        _kernelMakeWallSourceKST << <block, threads >> > (
            m_pDeviceData,
            static_cast<INT>(sourceData.m_sSourcePoint.w),
            sourceData.m_bySpinIndex,
            sourceData.m_byColorIndex,
            m_byFieldId);
    }
    break;
    default:
        appCrucial(_T("The source type %s not implemented yet!\n"), __ENUM_TO_STRING(EFermionBosonSource, sourceData.m_eSourceType).c_str());
        break;
    }
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
BYTE* CFieldFermionKST<deviceVector, deviceGauge, vectorN>::CopyDataOut(UINT& uiSize) const
{
    deviceVector* toSave = (deviceVector*)malloc(sizeof(deviceVector) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(Real) * m_uiSiteCount * 2 * vectorN);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        Real oneSite[2 * vectorN];
        for (UINT k = 0; k < 2 * vectorN; ++k)
        {
            oneSite[k] = _element(toSave[i], k);
        }
        memcpy(saveData + sizeof(Real) * i * 2 * vectorN, oneSite, sizeof(Real) * 2 * vectorN);
    }

    //appGetFileSystem()->WriteAllBytes(fileName.c_str(), saveData, uiSize);
    //free(saveData);
    free(toSave);
    return saveData;
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
BYTE* CFieldFermionKST<deviceVector, deviceGauge, vectorN>::CopyDataOutFloat(UINT& uiSize) const
{
    deviceVector* toSave = (deviceVector*)malloc(sizeof(deviceVector) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(FLOAT) * m_uiSiteCount * 2 * vectorN);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        FLOAT oneSite[2 * vectorN];
        for (UINT k = 0; k < 2 * vectorN; ++k)
        {
            oneSite[k] = static_cast<FLOAT>(_element(toSave[i], k));
        }
        memcpy(saveData + sizeof(FLOAT) * i * 2 * vectorN, oneSite, sizeof(FLOAT) * 2 * vectorN);
    }

    free(toSave);
    return saveData;
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
BYTE* CFieldFermionKST<deviceVector, deviceGauge, vectorN>::CopyDataOutDouble(UINT& uiSize) const
{
    deviceVector* toSave = (deviceVector*)malloc(sizeof(deviceVector) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(DOUBLE) * m_uiSiteCount * 2 * vectorN);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceVector) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        DOUBLE oneSite[2 * vectorN];
        for (UINT k = 0; k < 2 * vectorN; ++k)
        {
            oneSite[k] = static_cast<DOUBLE>(_element(toSave[i], k));
        }
        memcpy(saveData + sizeof(DOUBLE) * i * 2 * vectorN, oneSite, sizeof(DOUBLE) * 2 * vectorN);
    }

    free(toSave);
    return saveData;
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
TArray<CFieldFermion*> CFieldFermionKST<deviceVector, deviceGauge, vectorN>::GetSourcesAtSiteFromPool(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const SSmallInt4& site) const
{
    TArray<CFieldFermion*> ret;
    for (UINT j = 0; j < vectorN; ++j)
    {
        ret.AddItem(dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId)));
        if (NULL == ret[j])
        {
            appCrucial(_T("GetSourcesAtSiteFromPool failed!\n"));
            _FAIL_EXIT;
        }
    }

    for (BYTE c = 0; c < vectorN; ++c)
    {
        SFermionBosonSource sourceData;
        sourceData.m_eSourceType = EFS_Point;
        sourceData.m_sSourcePoint = site;
        sourceData.m_byColorIndex = c;
        sourceData.m_bySpinIndex = 0;

        ret[c]->InitialAsSource(sourceData);

        if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
        {
            ret[c]->m_fLength = ret[c]->Dot(ret[c]).x;
        }

        ret[c]->InverseD(gaugeNum, bosonNum, gaugeFields, pBoson);
    }
    return ret;
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::PrepareForHMCOnlyRandomize()
{
    preparethread;
    _kernelInitialFermionKST << <block, threads >> > (
        m_pDeviceData,
        m_byFieldId,
        EFIT_RandomGaussian);
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::PrepareForHMCNotRandomize(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson)
{
    D_MC(gaugeNum, bosonNum, gaugeFields, pBoson);
}

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
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::appApplyGammaKS(
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
void CFieldFermionKST<deviceVector, deviceGauge, vectorN>::GammaKSForce(
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


#pragma region Field Matrix Operation

template<typename deviceVector, typename deviceGauge, INT vectorN>
CFieldMatrixOperationKST<deviceVector, deviceGauge, vectorN>::CFieldMatrixOperationKST()
{
    m_pHostResBuffer = (deviceVector**)malloc(sizeof(deviceVector*) * _kFieldMatrixMaxDim);
    m_pHostLeftBuffer = (deviceVector**)malloc(sizeof(deviceVector*) * _kFieldMatrixMaxDim);
    checkCudaErrors(cudaMalloc((void**)&m_pResBuffer, sizeof(deviceVector*) * _kFieldMatrixMaxDim));
    checkCudaErrors(cudaMalloc((void**)&m_pLeftBuffer, sizeof(deviceVector*) * _kFieldMatrixMaxDim));
}

template<typename deviceVector, typename deviceGauge, INT vectorN>
CFieldMatrixOperationKST<deviceVector, deviceGauge, vectorN>::~CFieldMatrixOperationKST()
{
    free(m_pHostResBuffer);
    free(m_pHostLeftBuffer);

    checkCudaErrors(cudaFree(m_pResBuffer));
    checkCudaErrors(cudaFree(m_pLeftBuffer));
}

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
void CFieldMatrixOperationKST<deviceVector, deviceGauge, vectorN>::VectorMultiplyMatrix(TArray<CField*>& res, const TArray<CField*>& left, const CLGComplex* deviceMatrix, UINT uiDimX, UINT uiDimY)
{
    for (UINT i = 0; i < uiDimY; ++i)
    {
        CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pF = dynamic_cast<CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(res[i]);
        if (NULL == pF)
        {
            appCrucial(_T("CFieldMatrixOperationKSSU3 only work with CFieldFermionKST<deviceVector, deviceGauge, vectorN>!\n"));
            return;
        }
        m_pHostResBuffer[i] = pF->m_pDeviceData;
    }

    for (UINT i = 0; i < uiDimX - uiDimY; ++i)
    {
        const CFieldFermionKST<deviceVector, deviceGauge, vectorN>* pF = dynamic_cast<const CFieldFermionKST<deviceVector, deviceGauge, vectorN>*>(left[i]);
        if (NULL == pF)
        {
            appCrucial(_T("CFieldMatrixOperationKSSU3 only work with CFieldFermionKST<deviceVector, deviceGauge, vectorN>!\n"));
            return;
        }
        m_pHostLeftBuffer[i] = pF->m_pDeviceData;
    }

    checkCudaErrors(cudaMemcpy(m_pResBuffer, m_pHostResBuffer, sizeof(deviceVector*) * uiDimY, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(m_pLeftBuffer, m_pHostLeftBuffer, sizeof(deviceVector*) * (uiDimX - uiDimY), cudaMemcpyHostToDevice));

    preparethreadE(vectorN);
    _kernelMatrixMultiplyKST << <block, threads >> > (vectorN, m_pResBuffer, m_pLeftBuffer, deviceMatrix, uiDimX, uiDimY);
}

#pragma endregion

__CLG_FORCETEMPLATE_CONSTRUCTOR(CFieldMatrixOperationKST, U1, CLGComplex, CLGComplex, 1)

__CLGIMPLEMENT_CLASS(CFieldFermionKSU1)

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================