//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3.cu
// 
// DESCRIPTION:
// This is the device implementations of Wilson fermion
//
// This implementation assumes SU3 and square lattice
//
// REVISION:
//  [10/01/2021 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionKSU1)

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
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSU1(
    const CLGComplex * __restrict__ pDeviceData,
    const CLGComplex* __restrict__ pGauge,
    const SIndex * __restrict__ pGaugeMove,
    const SIndex * __restrict__ pFermionMove,
    const BYTE * __restrict__ pEtaTable,
    CLGComplex* pResultData,
    Real f2am,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernaldir;

    CLGComplex result = _zeroc;
    //pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];

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
        const CLGComplex& x_Gauge_element = pGauge[linkIndex];
        CLGComplex x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.y = -x_m_mu_Gauge_element.y;
        }

        //U(x,mu) phi(x+ mu)
        CLGComplex u_phi_x_p_m = _cuCmulf(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.x = -u_phi_x_p_m.x;
            u_phi_x_p_m.y = -u_phi_x_p_m.y;
        }

        //U^{dagger}(x-mu) phi(x-mu)
        CLGComplex u_dagger_phi_x_m_m = _cuCmulf(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.x += u_dagger_phi_x_m_m.x;
            u_phi_x_p_m.y += u_dagger_phi_x_m_m.y;
        }
        else
        {
            u_phi_x_p_m.x -= u_dagger_phi_x_m_m.x;
            u_phi_x_p_m.y -= u_dagger_phi_x_m_m.y;
        }
        u_phi_x_p_m.x *= eta_mu;
        u_phi_x_p_m.y *= eta_mu;
        result.x += u_phi_x_p_m.x;
        result.y += u_phi_x_p_m.y;
    }

    pResultData[uiSiteIndex].x = pDeviceData[uiSiteIndex].x * f2am;
    pResultData[uiSiteIndex].y = pDeviceData[uiSiteIndex].y * f2am;
    if (bDDagger)
    {
        pResultData[uiSiteIndex].x -= result.x;
        pResultData[uiSiteIndex].y -= result.y;
    }
    else
    {
        pResultData[uiSiteIndex].x += result.x;
        pResultData[uiSiteIndex].y += result.y;
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        pResultData[uiSiteIndex].x *= fCoeff;
        pResultData[uiSiteIndex].y *= fCoeff;
        break;
    case EOCT_Complex:
        pResultData[uiSiteIndex] = _cuCmulf(pResultData[uiSiteIndex], cCoeff);
        break;
    }
}

/**
 * For some strange boundary condition
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSPlusEtaU1(
    const CLGComplex* __restrict__ pDeviceData,
    const CLGComplex* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    CLGComplex* pResultData,
    Real f2am,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernaldir;

    CLGComplex result = _zeroc;
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
        const CLGComplex& x_Gauge_element = pGauge[linkIndex];
        CLGComplex x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.y = -x_m_mu_Gauge_element.y;
        }

        //U(x,mu) phi(x+ mu)
        CLGComplex u_phi_x_p_m = _cuCmulf(x_Gauge_element, pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.x *= F(-1.0) * eta_mu;
            u_phi_x_p_m.y *= F(-1.0) * eta_mu;
        }
        else
        {
            u_phi_x_p_m.x *= eta_mu;
            u_phi_x_p_m.y *= eta_mu;
        }

        //U^{dagger}(x-mu) phi(x-mu)
        CLGComplex u_dagger_phi_x_m_m = _cuCmulf(x_m_mu_Gauge_element, pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        u_dagger_phi_x_m_m.x *= eta_mu2;
        u_dagger_phi_x_m_m.y *= eta_mu2;
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.x += u_dagger_phi_x_m_m.x;
            u_phi_x_p_m.y += u_dagger_phi_x_m_m.y;
        }
        else
        {
            u_phi_x_p_m.x -= u_dagger_phi_x_m_m.x;
            u_phi_x_p_m.y -= u_dagger_phi_x_m_m.y;
        }
        result.x += u_phi_x_p_m.x;
        result.y += u_phi_x_p_m.y;
    }

    pResultData[uiSiteIndex].x *= f2am;
    pResultData[uiSiteIndex].y *= f2am;
    if (bDDagger)
    {
        pResultData[uiSiteIndex].x -= result.x;
        pResultData[uiSiteIndex].y -= result.y;
    }
    else
    {
        pResultData[uiSiteIndex].x += result.x;
        pResultData[uiSiteIndex].y += result.y;
    }

    switch (eCoeff)
    {
    case EOCT_Real:
        pResultData[uiSiteIndex].x *= fCoeff;
        pResultData[uiSiteIndex].y *= fCoeff;
        break;
    case EOCT_Complex:
        pResultData[uiSiteIndex] = _cuCmulf(pResultData[uiSiteIndex], cCoeff);
        break;
    }
}

/**
 * Calculate Force
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForceU1(
    const CLGComplex* __restrict__ pGauge,
    CLGComplex* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const CLGComplex* const* __restrict__ pFermionPointers,
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
            const CLGComplex* phi_i = pFermionPointers[uiR];
            const CLGComplex* phi_id = pFermionPointers[uiR + uiRational];

            CLGComplex toContract = _cuCmulf(pGauge[linkIndex], phi_i[x_p_mu_Fermion.m_uiSiteIndex]);
            CLGComplex thisTerm = _cuCmulf(_cuConjf(phi_id[uiSiteIndex]), toContract);

            toContract = _cuCmulf(pGauge[linkIndex], phi_id[x_p_mu_Fermion.m_uiSiteIndex]);
            toContract.y = -toContract.y;
            thisTerm = _cuCaddf(thisTerm, _cuCmulf(toContract, phi_i[uiSiteIndex]));

            //if (x_p_mu_Fermion.NeedToOpposite())
            //{
            //    thisTerm = cuCmulf_cr(thisTerm, eta_mu * pNumerators[uiR] * F(-1.0));
            //}
            //else
            //{
            //    thisTerm = cuCmulf_cr(thisTerm, eta_mu * pNumerators[uiR]);
            //}

            pForce[linkIndex].y -= thisTerm.y * (x_p_mu_Fermion.NeedToOpposite() ? 
                (eta_mu * pNumerators[uiR] * F(-1.0)) : 
                (eta_mu * pNumerators[uiR]));
        }
    }
}

#pragma endregion

void CFieldFermionKSU1::DOperatorKS(void* pTargetBuffer, const void * pBuffer,
    const void * pGaugeBuffer, BYTE byGaugeFieldId, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    CLGComplex* pTarget = (CLGComplex*)pTargetBuffer;
    const CLGComplex* pSource = (const CLGComplex*)pBuffer;
    const CLGComplex* pGauge = (const CLGComplex*)pGaugeBuffer;
    //appGeneral(_T("Shift: %d, am: %f\n"), m_bEachSiteEta, f2am);

    preparethread;
    if (m_bEachSiteEta)
    {
        _kernelDFermionKSPlusEtaU1 << <block, threads >> > (
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
        _kernelDFermionKSU1 << <block, threads >> > (
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
void CFieldFermionKSU1::DerivateD0(
    void* pForce, 
    const void* pGaugeBuffer, BYTE byGaugeFieldId) const
{
    preparethread;
    _kernelDFermionKSForceU1 << <block, threads >> > (
        (const CLGComplex*)pGaugeBuffer,
        (CLGComplex*)pForce,
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
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialFermionKSU1(CLGComplex *pDevicePtr, BYTE byFieldId, EFieldInitialType eInitialType)
{
    intokernalInt4;

    switch (eInitialType)
    {
    case EFIT_Zero:
    {
        pDevicePtr[uiSiteIndex] = _zeroc;
    }
    break;
    case EFIT_Identity:
    {
        pDevicePtr[uiSiteIndex] = _onec;
    }
    break;
    case EFIT_RandomGaussian:
    {
        pDevicePtr[uiSiteIndex] = _deviceRandomGaussC(_deviceGetFatIndex(uiSiteIndex, 0));
    }
    break;
    case EFIT_RandomZ4:
    {
        pDevicePtr[uiSiteIndex] = _deviceRandomZ4(_deviceGetFatIndex(uiSiteIndex, 0));
    }
    break;
    default:
    {
        printf("Wilson Fermion Field cannot be initialized with this type!");
    }
    break;
    }
}

//void CFieldFermionKSSU3::D_EN(const CField* pGauge)
//{
//    RationalApproximation(EFO_F_DDdagger, pGauge, &m_rEN);
//}

/**
* generate phi by gaussian random.
* phi = (D^+D)^{1/4} phi
*/
void CFieldFermionKSU1::PrepareForHMC(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson)
{
    preparethread;
    _kernelInitialFermionKSU1 << <block, threads >> > (
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
UBOOL CFieldFermionKSU1::CalculateForceS(
    const CFieldGauge* pGauge,
    CFieldGauge* pForce,
    ESolverPhase ePhase) const
{
    if (NULL == pGauge || EFT_GaugeU1 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    if (NULL == pForce || EFT_GaugeU1 != pForce->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }

    TArray<CField*> phii;
    TArray<CFieldFermionKSU1*> phiid;
    for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
    {
        CField* pPhi_i = dynamic_cast<CField*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
        phii.AddItem(pPhi_i);
        CFieldFermionKSU1* pPhi_id = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
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
    TArray<const CFieldGauge*> gauge;
    gauge.AddItem(pGauge);
    solver->Solve(phii, shifts, this, 1, 0, gauge.GetData(), NULL, EFO_F_DDdagger);

    const UINT uiBufferSize = sizeof(CLGComplex*) * 2 * m_rMD.m_uiDegree;
    CLGComplex** hostPointers = (CLGComplex**)appAlloca(uiBufferSize);
    for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
    {
        CFieldFermionKSU1* phi_ks = dynamic_cast<CFieldFermionKSU1*>(phii[i]);
        phi_ks->CopyTo(phiid[i]);
        phiid[i]->D0S(pGauge);

        hostPointers[i] = phi_ks->m_pDeviceData;
        hostPointers[i + m_rMD.m_uiDegree] = phiid[i]->m_pDeviceData;
    }
    checkCudaErrors(cudaMemcpy(m_pRationalFieldPointers, hostPointers, uiBufferSize, cudaMemcpyHostToDevice));

    const CFieldGaugeU1* pGaugeU1 = dynamic_cast<const CFieldGaugeU1*>(pGauge);
    CFieldGaugeU1* pForceU1 = dynamic_cast<CFieldGaugeU1*>(pForce);

    DerivateD0(pForceU1->m_pDeviceData, pGaugeU1->m_pDeviceData, pGaugeU1->m_byFieldId);
    //preparethread;
    //_kernelDFermionKSForce << <block, threads >> > (
    //    pGaugeSU3->m_pDeviceData,
    //    pForceSU3->m_pDeviceData,
    //    appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
    //    appGetLattice()->m_pIndexCache->m_pEtaMu,
    //    m_pRationalFieldPointers,
    //    m_pMDNumerator,
    //    m_rMD.m_uiDegree,
    //    m_byFieldId);

    
    for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
    {
        phii[i]->Return();
        phiid[i]->Return();
    }

    return TRUE;
}

#pragma endregion

#pragma region Kernel

__global__ void _CLG_LAUNCH_BOUND
_kernelPrintFermionKSU1(const CLGComplex * __restrict__ pData)
{
    intokernal;

    printf("%d=(%1.2f %1.2fi)\n",
        uiSiteIndex,
        pData[uiSiteIndex].x, pData[uiSiteIndex].y);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusFermionKSU1(CLGComplex* pMe, const CLGComplex* __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex].x += pOther[uiSiteIndex].x;
    pMe[uiSiteIndex].y += pOther[uiSiteIndex].y;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusFermionKSU1(CLGComplex* pMe, const CLGComplex* __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex].x -= pOther[uiSiteIndex].x;
    pMe[uiSiteIndex].y -= pOther[uiSiteIndex].y;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyComplexFermionKSU1(CLGComplex* pMe, const CLGComplex* __restrict__ pOther, CLGComplex a)
{
    intokernal;
    pMe[uiSiteIndex] = _cuCaddf(pMe[uiSiteIndex], _cuCmulf(pOther[uiSiteIndex], a));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMulFermionKSU1(CLGComplex* pMe, const CLGComplex* __restrict__ pOther, UBOOL bConj)
{
    intokernal;
    if (bConj)
    {
        pMe[uiSiteIndex].y = -pMe[uiSiteIndex].y;
    }
    pMe[uiSiteIndex] = _cuCmulf(pMe[uiSiteIndex], pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyRealFermionKSU1(CLGComplex* pMe, const CLGComplex* __restrict__ pOther, Real a)
{
    intokernal;
    pMe[uiSiteIndex].x += pOther[uiSiteIndex].x * a;
    pMe[uiSiteIndex].y += pOther[uiSiteIndex].y * a;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDotFermionKSU1(const CLGComplex* __restrict__ pMe, const CLGComplex* __restrict__ pOther,
#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex* result
#else
    CLGComplex* result
#endif
)
{
    intokernal;
#if !_CLG_DOUBLEFLOAT
    result[uiSiteIndex] = _cToDouble(_cuCmulf(_cuConjf(pMe[uiSiteIndex]), pOther[uiSiteIndex]));
#else
    result[uiSiteIndex] = _cuCmulf(_cuConjf(pMe[uiSiteIndex]), pOther[uiSiteIndex]);
#endif
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyComplexKSU1(CLGComplex* pMe, CLGComplex a)
{
    intokernal;
    pMe[uiSiteIndex] = _cuCmulf(pMe[uiSiteIndex], a);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyRealKSU1(CLGComplex* pMe, Real a)
{
    intokernal;
    pMe[uiSiteIndex].x *= a;
    pMe[uiSiteIndex].y *= a;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFermionKSConjugateU1(CLGComplex* pDeviceData)
{
    intokernal;
    pDeviceData[uiSiteIndex].y = -pDeviceData[uiSiteIndex].y;
}


__global__ void _CLG_LAUNCH_BOUND
_kernelMakePointSourceKSU1(CLGComplex* pDeviceData, UINT uiDesiredSite)
{
    intokernal;
    if (uiSiteIndex == uiDesiredSite)
    {
        pDeviceData[uiSiteIndex] = _onec;
    }
    else
    {
        pDeviceData[uiSiteIndex] = _zeroc;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMakeWallSourceKSU1(CLGComplex* pDeviceData,
    INT uiDesiredT, UINT uiShift, BYTE byFieldID)
{
    intokernalOnlyInt4;

    //pDeviceData[uiSiteIndex] = deviceSU3Vector::makeZeroSU3Vector();
    //We shall not set zero here!

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
            pDeviceData[sIdx.m_uiSiteIndex] = _onec;
        }
    }
}


__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGammaU1(
    CLGComplex* pMe,
    const CLGComplex* __restrict__ pOther,
    const CLGComplex* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    BYTE byDir)
{
    intokernal;
    const Real eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> byDir) & 1)) ? F(-1.0) : F(1.0);
    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byDir);
    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    const CLGComplex& x_Gauge_element = pGauge[linkIndex];
    CLGComplex x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, byDir)];
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.y = -x_m_mu_Gauge_element.y;
    }

    pMe[uiSiteIndex] = _cuCmulf(x_Gauge_element, pOther[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        pMe[uiSiteIndex].x = pMe[uiSiteIndex].x * F(-1.0);
        pMe[uiSiteIndex].y = pMe[uiSiteIndex].y * F(-1.0);
    }
    if (x_m_mu_Fermion.NeedToOpposite())
    {
        pMe[uiSiteIndex] = _cuCsubf(pMe[uiSiteIndex], _cuCmulf(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {
        pMe[uiSiteIndex] = _cuCaddf(pMe[uiSiteIndex], _cuCmulf(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    pMe[uiSiteIndex] = cuCmulf_cr(pMe[uiSiteIndex], F(0.5) * eta_mu);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelKSApplyGammaEtaU1(
    CLGComplex* pMe,
    const CLGComplex* __restrict__ pOther,
    const CLGComplex* __restrict__ pGauge,
    const SIndex* __restrict__ pGaugeMove,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    BYTE byDir)
{
    intokernal;

    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, byDir);
    const SIndex& x_m_mu_Gauge = pGaugeMove[linkIndex];
    const SIndex& x_p_mu_Fermion = pFermionMove[2 * linkIndex];
    const SIndex& x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

    BYTE eta_mu = (1 == ((pEtaTable[uiSiteIndex] >> byDir) & 1));
    BYTE eta_mu2 = (1 == ((pEtaTable[x_m_mu_Fermion.m_uiSiteIndex] >> byDir) & 1));

    const CLGComplex& x_Gauge_element = pGauge[linkIndex];
    CLGComplex x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, byDir)];
    if (x_m_mu_Gauge.NeedToDagger())
    {
        x_m_mu_Gauge_element.y = -x_m_mu_Gauge_element.y;
    }

    pMe[uiSiteIndex] = _cuCmulf(x_Gauge_element, pOther[x_p_mu_Fermion.m_uiSiteIndex]);
    if (x_p_mu_Fermion.NeedToOpposite())
    {
        eta_mu = eta_mu + 1;
    }

    if (eta_mu & 1)
    {
        pMe[uiSiteIndex] = cuCmulf_cr(pMe[uiSiteIndex], F(-1.0));
    }

    if (x_m_mu_Fermion.NeedToOpposite())
    {
        eta_mu2 = eta_mu2 + 1;
    }
    if (eta_mu2 & 1)
    {
        pMe[uiSiteIndex] = _cuCsubf(pMe[uiSiteIndex], _cuCmulf(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    else
    {

        pMe[uiSiteIndex] = _cuCaddf(pMe[uiSiteIndex], _cuCmulf(x_m_mu_Gauge_element, pOther[x_m_mu_Fermion.m_uiSiteIndex]));
    }
    pMe[uiSiteIndex] = cuCmulf_cr(pMe[uiSiteIndex], F(0.5));
}


#pragma endregion


#pragma region Helper functions to implement higher orders

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSForce_WithLinkU1(
    const CLGComplex* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    CLGComplex* pForce,
    const CLGComplex* const* __restrict__ pFermionPointers,
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
    INT pathLeft[CFieldFermionKSSU3::_kKSLinkLength];
    INT pathRight[CFieldFermionKSSU3::_kKSLinkLength];
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
            const CLGComplex vnn1 = _deviceLinkU1(pGauge, sSite4, LLength, byGaugeFieldId, pathLeft);
            const CLGComplex vnn2 = _deviceLinkU1(pGauge, sSite4, RLength, byGaugeFieldId, pathRight);

            for (BYTE rfieldId = 0; rfieldId < uiRational; ++rfieldId)
            {
                const CLGComplex* phi_i = pFermionPointers[rfieldId];
                const CLGComplex* phi_id = pFermionPointers[rfieldId + uiRational];

                //=================================
                // 3. Find phi_{1,2,3,4}(n1), phi_i(n2)
                CLGComplex phi1 = _cuCmulf(vnn1, phi_id[sn1.m_uiSiteIndex]);
                CLGComplex phi2 = _cuCmulf(vnn2, phi_i[sn2.m_uiSiteIndex]);
                //deviceSU3Vector phi3 = vnn1.MulVector(phi_i[sn1.m_uiSiteIndex]);
                //deviceSU3Vector phi4 = vnn2.MulVector(phi_id[sn2.m_uiSiteIndex]);

                CLGComplex res = _cuCmulf(_cuConjf(phi1), phi2);
                //This Add is required by partial(D^+D)
                phi1 = _cuCmulf(vnn2, phi_id[sn2.m_uiSiteIndex]);
                phi2 = _cuCmulf(vnn1, phi_i[sn1.m_uiSiteIndex]);
                res = _cuCaddf(res, _cuCmulf(_cuConjf(phi1), phi2));
                //res.Add(deviceSU3::makeSU3ContractV(phi4, phi3));
                res.x = 0;
                res.y = res.y * fCoefficient * pNumerators[rfieldId];

                if (bHasLeft)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathLeft[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        pForce[linkIndex].y = pForce[linkIndex].y - res.y;
                    }
                    else
                    {
                        pForce[linkIndex].y = pForce[linkIndex].y + res.y;
                    }
                }

                if (bHasRight)
                {
                    const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, pathRight[0] - 1);
                    if (iEtaMu1 & 1)
                    {
                        pForce[linkIndex].y = pForce[linkIndex].y + res.y;
                    }
                    else
                    {
                        pForce[linkIndex].y = pForce[linkIndex].y - res.y;
                    }
                }
            }
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_OneLinkU1(
    const CLGComplex* __restrict__ pDeviceData,
    const CLGComplex* __restrict__ pGauge,
    const BYTE* __restrict__ pEtaTable,
    CLGComplex* pResultData,
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
    INT pathBuffer[CFieldFermionKSSU3::_kKSLinkLength];
    CLGComplex result = _zeroc;

    SSmallInt4 siten = _deviceSmallInt4OffsetC(sSite4, path, pathLength);
    const SIndex& sn1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten)];
    CLGComplex vn = _deviceLinkU1(pGauge, sSite4, pathLength, byGaugeFieldId, path);
    INT etamu = (pEtaTable[uiSiteIndex] >> byEtaIdx) & 1;
    if (sn1.NeedToOpposite())
    {
        etamu = etamu + 1;
    }
    if (etamu & 1)
    {
        result = _cuCsubf(result, _cuCmulf(vn, pDeviceData[sn1.m_uiSiteIndex]));
    }
    else
    {
        result = _cuCaddf(result, _cuCmulf(vn, pDeviceData[sn1.m_uiSiteIndex]));
    }

    _devicePathDagger(path, pathBuffer, pathLength);
    siten = _deviceSmallInt4OffsetC(sSite4, pathBuffer, pathLength);
    const SIndex& sn2 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(siten)];
    vn = _deviceLinkU1(pGauge, sSite4, pathLength, byGaugeFieldId, pathBuffer);
    etamu = (pEtaTable[sn2.m_uiSiteIndex] >> byEtaIdx) & 1;
    if (sn2.NeedToOpposite())
    {
        etamu = etamu + 1;
    }
    if (etamu & 1)
    {
        result = _cuCaddf(result, _cuCmulf(vn, pDeviceData[sn2.m_uiSiteIndex]));
    }
    else
    {
        result = _cuCsubf(result, _cuCmulf(vn, pDeviceData[sn2.m_uiSiteIndex]));
    }

    if (bDDagger)
    {
        fCoefficient = fCoefficient * F(-1.0);
    }
    result = cuCmulf_cr(result, fCoefficient);

    switch (eCoeff)
    {
    case EOCT_Real:
        result = cuCmulf_cr(result, fCoeff);
        break;
    case EOCT_Complex:
        result = _cuCmulf(result, cCoeff);
        break;
    }

    pResultData[uiSiteIndex] = _cuCaddf(pResultData[uiSiteIndex], result);
}


__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKS_OnlyMassU1(
    const CLGComplex* __restrict__ pDeviceData,
    CLGComplex* pResultData,
    Real f2am,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernal;
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];
    pResultData[uiSiteIndex] = cuCmulf_cr(pResultData[uiSiteIndex], f2am);

    switch (eCoeff)
    {
    case EOCT_Real:
        pResultData[uiSiteIndex] = cuCmulf_cr(pResultData[uiSiteIndex], fCoeff);
        break;
    case EOCT_Complex:
        pResultData[uiSiteIndex] = _cuCmulf(pResultData[uiSiteIndex], cCoeff);
        break;
    }
}

void CFieldFermionKSU1::OnlyMass(void* pTarget, Real fm, EOperatorCoefficientType eOCT, Real fRealCoeff, const CLGComplex& cCmpCoeff)
{
    preparethread;
    _kernelDFermionKS_OnlyMassU1 << <block, threads >> > (
        m_pDeviceData,
        (CLGComplex*)pTarget,
        fm,
        eOCT,
        fRealCoeff,
        cCmpCoeff
        );
}

void CFieldFermionKSU1::OneLinkS(
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
    assert(pathLength <= CFieldFermionKS::_kKSLinkLength);
    preparethread;
    _kernelDFermionKS_OneLinkU1 << <block, threads >> > (
        m_pDeviceData,
        (const CLGComplex*)pGuage,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        (CLGComplex*)pTarget,
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

void CFieldFermionKSU1::OneLinkForceS(
    const void* pGuage,
    BYTE byGaugeFieldId, 
    void* pForce,
    Real fCoefficient,
    const INT* pDevicePath, 
    BYTE pathLength, 
    BYTE byEtaIdx) const
{
    assert(pathLength <= CFieldFermionKS::_kKSLinkLength);
    preparethread;
    _kernelDFermionKSForce_WithLinkU1 << <block, threads >> > (
        (const CLGComplex*)pGuage,
        appGetLattice()->m_pIndexCache->m_pEtaMu,
        (CLGComplex*)pForce,
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

CFieldFermionKSU1::CFieldFermionKSU1()
    : CFieldFermionKS()
    , m_pRationalFieldPointers(NULL)
{
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount));
}

CFieldFermionKSU1::~CFieldFermionKSU1()
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
void CFieldFermionKSU1::InitialField(EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialFermionKSU1 << <block, threads >> > (m_pDeviceData, m_byFieldId, eInitialType);
}

void CFieldFermionKSU1::Dagger()
{
    preparethread;
    _kernelFermionKSConjugateU1 << <block, threads >> > (m_pDeviceData);
}

void CFieldFermionKSU1::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFieldType)
{
    if (eFieldType != EFFT_CLGBin)
    {
        appCrucial(_T("CFieldFermionKSSU3::InitialFieldWithFile: Only support CLG Bin File\n"));
        return;
    }

    UINT uiSize = static_cast<UINT>(sizeof(Real) * 2 * m_uiSiteCount);
    BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
    InitialWithByte(data);
    free(data);
}

void CFieldFermionKSU1::InitialWithByte(BYTE* byData)
{
    CLGComplex* readData = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiSiteCount);
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        Real thisSite[2];
        memcpy(thisSite, byData + i * sizeof(Real) * 2, sizeof(Real) * 2);
        readData[i] = _make_cuComplex(
            thisSite[0],
            thisSite[1]);
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyHostToDevice));
    free(readData);
}

void CFieldFermionKSU1::InitialOtherParameters(CParameters& params)
{
    CFieldFermionKS::InitialOtherParameters(params);
    if (NULL != m_pRationalFieldPointers)
    {
        checkCudaErrors(cudaFree(m_pRationalFieldPointers));
    }
    checkCudaErrors(cudaMalloc((void**)&m_pRationalFieldPointers, sizeof(CLGComplex*) * 2 * m_rMD.m_uiDegree));
}

void CFieldFermionKSU1::DebugPrintMe() const
{
    CLGComplex* toprint = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiSiteCount);
    checkCudaErrors(cudaMemcpy(toprint, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT uiSite = 0; uiSite < m_uiSiteCount; ++uiSite)
    {
        const SSmallInt4 site4 = __hostSiteIndexToInt4(uiSite);
        appGeneral(_T(" --- %d,%d,%d,%d --- %f %s %f I\n"),
            site4.x, site4.y, site4.z, site4.w,
            toprint[uiSite].x,
            toprint[uiSite].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(toprint[uiSite].y)
        );
    }

    appSafeFree(toprint);
}

void CFieldFermionKSU1::CopyTo(CField* U) const
{
    if (NULL == U || EFT_FermionStaggeredU1 != U->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only copy to CFieldFermionKSU1!"));
        return;
    }

    CFieldFermionKS::CopyTo(U);

    CFieldFermionKSU1* pField = dynamic_cast<CFieldFermionKSU1*>(U);

    if (NULL != pField->m_pMDNumerator)
    {
        checkCudaErrors(cudaFree(pField->m_pMDNumerator));
    }
    if (NULL != pField->m_pRationalFieldPointers)
    {
        checkCudaErrors(cudaFree(pField->m_pRationalFieldPointers));
    }

    checkCudaErrors(cudaMalloc((void**)&pField->m_pRationalFieldPointers, sizeof(CLGComplex*) * 2 * m_rMD.m_uiDegree));
    checkCudaErrors(cudaMalloc((void**)&pField->m_pMDNumerator, sizeof(Real) * m_rMD.m_uiDegree));

    checkCudaErrors(cudaMemcpy(pField->m_pDeviceData, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
    checkCudaErrors(cudaMemcpy(pField->m_pMDNumerator, m_pMDNumerator, sizeof(Real) * m_rMD.m_uiDegree, cudaMemcpyDeviceToDevice));
}

void CFieldFermionKSU1::AxpyPlus(const CField* x)
{
    if (NULL == x || EFT_FermionStaggeredU1 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only copy to CFieldFermionKSU1!"));
        return;
    }
    const CFieldFermionKSU1* pField = dynamic_cast<const CFieldFermionKSU1*>(x);

    preparethread;
    _kernelAxpyPlusFermionKSU1 << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldFermionKSU1::AxpyMinus(const CField* x)
{
    if (NULL == x || EFT_FermionStaggeredU1 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only copy to CFieldFermionKSU1!"));
        return;
    }
    const CFieldFermionKSU1* pField = dynamic_cast<const CFieldFermionKSU1*>(x);

    preparethread;
    _kernelAxpyMinusFermionKSU1 << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldFermionKSU1::Axpy(Real a, const CField* x)
{
    if (NULL == x || EFT_FermionStaggeredU1 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only copy to CFieldFermionKSU1!"));
        return;
    }
    const CFieldFermionKSU1* pField = dynamic_cast<const CFieldFermionKSU1*>(x);

    preparethread;
    _kernelAxpyRealFermionKSU1 << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

void CFieldFermionKSU1::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || EFT_FermionStaggeredU1 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only copy to CFieldFermionKSU1!"));
        return;
    }
    const CFieldFermionKSU1* pField = dynamic_cast<const CFieldFermionKSU1*>(x);

    preparethread;
    _kernelAxpyComplexFermionKSU1 << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

void CFieldFermionKSU1::Mul(const CField* other, UBOOL bDagger)
{
    if (NULL == other || EFT_FermionStaggeredU1 != other->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only copy to CFieldFermionKSU1!"));
        return;
    }
    const CFieldFermionKSU1* pField = dynamic_cast<const CFieldFermionKSU1*>(other);

    preparethread;
    _kernelMulFermionKSU1 << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, bDagger);
}

cuDoubleComplex CFieldFermionKSU1::Dot(const CField* x) const
{
    if (NULL == x || EFT_FermionStaggeredU1 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only copy to CFieldFermionKSU1!"));
        return make_cuDoubleComplex(0, 0);
    }
    const CFieldFermionKSU1* pField = dynamic_cast<const CFieldFermionKSU1*>(x);
    preparethread;
    _kernelDotFermionKSU1 << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, _D_ComplexThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

void CFieldFermionKSU1::ScalarMultply(const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplyComplexKSU1 << <block, threads >> > (m_pDeviceData, a);
}

void CFieldFermionKSU1::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplyRealKSU1 << <block, threads >> > (m_pDeviceData, a);
}

void CFieldFermionKSU1::ApplyGamma(EGammaMatrix eGamma)
{
    appCrucial(_T("Not implemented yet...\n"));
}

void CFieldFermionKSU1::ApplyGammaKSS(const CFieldGauge* pGauge, EGammaMatrix eGamma)
{
    INT iDir = -1;
    switch (eGamma)
    {
    case UNITY:
        {
            return;
        }
    case GAMMA1:
        {
            iDir = 0;
        }
        break;
    case GAMMA2:
        {
            iDir = 1;
        }
        break;
    case GAMMA3:
        {
            iDir = 2;
        }
        break;
    case GAMMA4:
        {
            iDir = 3;
        }
        break;
    }

    if (iDir >= 0)
    {
        if (NULL == pGauge || EFT_GaugeU1 != pGauge->GetFieldType())
        {
            appCrucial(_T("CFieldFermionKSU1 can only play with gauge U1!"));
            return;
        }
        const CFieldGaugeU1* pFieldU1 = dynamic_cast<const CFieldGaugeU1*>(pGauge);
        CFieldFermionKSU1* pPooled = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
        checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
        preparethread;
        if (m_bEachSiteEta)
        {
            _kernelKSApplyGammaEtaU1 << <block, threads >> > (
                m_pDeviceData,
                pPooled->m_pDeviceData,
                pFieldU1->m_pDeviceData,
                appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                static_cast<BYTE>(iDir));
        }
        else
        {
            _kernelKSApplyGammaU1 << <block, threads >> > (
                m_pDeviceData,
                pPooled->m_pDeviceData,
                pFieldU1->m_pDeviceData,
                appGetLattice()->m_pIndexCache->m_pGaugeMoveCache[m_byFieldId],
                appGetLattice()->m_pIndexCache->m_pMoveCache[m_byFieldId],
                appGetLattice()->m_pIndexCache->m_pEtaMu,
                static_cast<BYTE>(iDir));
        }

        pPooled->Return();
        return;
    }

    appCrucial(_T("Not implemented yet...\n"));

}


//Kai should be part of D operator
void CFieldFermionKSU1::DS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeU1 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only play with gauge U1!"));
        return;
    }
    const CFieldGaugeU1* pFieldU1 = dynamic_cast<const CFieldGaugeU1*>(pGauge);
    CFieldFermionKSU1* pPooled = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldU1->m_pDeviceData, pFieldU1->m_byFieldId,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

void CFieldFermionKSU1::DWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeU1 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only play with gauge U1!"));
        return;
    }
    const CFieldGaugeU1* pFieldU1 = dynamic_cast<const CFieldGaugeU1*>(pGauge);
    CFieldFermionKSU1* pPooled = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldU1->m_pDeviceData, pFieldU1->m_byFieldId, fMass,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

void CFieldFermionKSU1::D0S(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeU1 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only play with gauge U1!"));
        return;
    }
    const CFieldGaugeU1* pFieldU1 = dynamic_cast<const CFieldGaugeU1*>(pGauge);
    CFieldFermionKSU1* pPooled = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldU1->m_pDeviceData, pFieldU1->m_byFieldId, F(0.0),
        FALSE, EOCT_None, F(1.0), _onec);

    pPooled->Return();
}

UINT CFieldFermionKSU1::TestAntiHermitianS(const CFieldGauge* pGauge) const
{
    const UINT uiVolume = _HC_Volume;
    CLGComplex* matrixElement = (CLGComplex*)malloc(sizeof(CLGComplex) * uiVolume * uiVolume);
    CLGComplex* hostData = (CLGComplex*)malloc(sizeof(CLGComplex) * uiVolume);
    CFieldFermionKSU1* v = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    for (UINT x = 0; x < uiVolume; ++x)
    {
        const SSmallInt4 point = __hostSiteIndexToInt4(x);
        SFermionSource source;
        source.m_byColorIndex = 0;
        source.m_eSourceType = EFS_Point;
        source.m_sSourcePoint = point;
        v->InitialAsSource(source);
        v->D0S(pGauge);

        checkCudaErrors(cudaMemcpy(hostData, v->m_pDeviceData, sizeof(CLGComplex) * uiVolume, cudaMemcpyDeviceToHost));

        for (UINT y = 0; y < uiVolume; ++y)
        {
            matrixElement[y * uiVolume + x] = hostData[y];
        }
        appGeneral(_T("%d / %d have been done\n"), x, uiVolume);
    }

    UINT uiE = 0;
    UINT uiWrong = 0;
    //List all results
    for (UINT i = 0; i < uiVolume * uiVolume; ++i)
    {
        const UINT x = i / uiVolume;
        const UINT y = i % uiVolume;
        const SSmallInt4 xSite = __hostSiteIndexToInt4(x);
        const SSmallInt4 ySite = __hostSiteIndexToInt4(y);
        const UINT daggerIdx = y * uiVolume + x;

        if (_cuCabsf(matrixElement[i]) > F(0.0000001))
        {
            ++uiE;
            if (appAbs(matrixElement[i].x + matrixElement[daggerIdx].x) > F(0.0000001)
             || appAbs(matrixElement[i].y - matrixElement[daggerIdx].y) > F(0.0000001))
            {
                ++uiWrong;
                appGeneral(_T("[(%d, %d, %d, %d)-(%d, %d, %d, %d)]: D = %f + %f I   Ddagger = %f + %f I\n"),
                    xSite.x, xSite.y, xSite.z, xSite.w, 
                    ySite.x, ySite.y, ySite.z, ySite.w, 
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
void CFieldFermionKSU1::DdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeU1 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only play with gauge U1!"));
        return;
    }
    const CFieldGaugeU1* pFieldU1 = dynamic_cast<const CFieldGaugeU1*>(pGauge);
    CFieldFermionKSU1* pPooled = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldU1->m_pDeviceData, pFieldU1->m_byFieldId,
        TRUE, eCoeffType, fRealCoeff, cCompCoeff);


    pPooled->Return();
}

void CFieldFermionKSU1::DdaggerWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeU1 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only play with gauge U1!"));
        return;
    }
    const CFieldGaugeU1* pFieldU1 = dynamic_cast<const CFieldGaugeU1*>(pGauge);
    CFieldFermionKSU1* pPooled = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldU1->m_pDeviceData, pFieldU1->m_byFieldId, fMass,
        TRUE, eCoeffType, fRealCoeff, cCompCoeff);


    pPooled->Return();
}

void CFieldFermionKSU1::DDdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeU1 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only play with gauge U1!"));
        return;
    }
    const CFieldGaugeU1* pFieldU1 = dynamic_cast<const CFieldGaugeU1*>(pGauge);

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }
    CFieldFermionKSU1* pPooled = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldU1->m_pDeviceData, pFieldU1->m_byFieldId,
        TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //why only apply coeff in the next step?
    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldU1->m_pDeviceData, pFieldU1->m_byFieldId,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

void CFieldFermionKSU1::DDS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeU1 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only play with gauge U1!"));
        return;
    }
    const CFieldGaugeU1* pFieldU1 = dynamic_cast<const CFieldGaugeU1*>(pGauge);

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }
    CFieldFermionKSU1* pPooled = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldU1->m_pDeviceData, pFieldU1->m_byFieldId,
        FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //why only apply coeff in the next step?
    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldU1->m_pDeviceData, pFieldU1->m_byFieldId,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

void CFieldFermionKSU1::DDdaggerWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeU1 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only play with gauge U1!"));
        return;
    }
    const CFieldGaugeU1* pFieldU1 = dynamic_cast<const CFieldGaugeU1*>(pGauge);

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }
    CFieldFermionKSU1* pPooled = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    DOperatorKS(pPooled->m_pDeviceData, m_pDeviceData, pFieldU1->m_pDeviceData, pFieldU1->m_byFieldId, fMass,
        TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //why only apply coeff in the next step?
    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldU1->m_pDeviceData, pFieldU1->m_byFieldId, fMass,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

void CFieldFermionKSU1::DDWithMassS(const CField* pGauge, Real fMass, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeU1 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSU1 can only play with gauge U1!"));
        return;
    }
    const CFieldGaugeU1* pFieldU1 = dynamic_cast<const CFieldGaugeU1*>(pGauge);

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }
    CFieldFermionKSU1* pPooled = dynamic_cast<CFieldFermionKSU1*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    DOperatorKS(pPooled->m_pDeviceData, m_pDeviceData, pFieldU1->m_pDeviceData, pFieldU1->m_byFieldId, fMass,
        FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //why only apply coeff in the next step?
    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldU1->m_pDeviceData, pFieldU1->m_byFieldId, fMass,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

void CFieldFermionKSU1::InitialAsSource(const SFermionSource& sourceData)
{
    const UINT uiSiteIndex = _hostGetSiteIndex(sourceData.m_sSourcePoint);
    switch (sourceData.m_eSourceType)
    {
    case EFS_Point:
    {
        preparethread;
        _kernelMakePointSourceKSU1 << <block, threads >> > (m_pDeviceData, uiSiteIndex);
    }
    break;
    case EFS_Wall:
    {
        preparethread;
        _kernelInitialFermionKSU1 << <block, threads >> > (m_pDeviceData, m_byFieldId, EFIT_Zero);
        _kernelMakeWallSourceKSU1 << <block, threads >> > (
            m_pDeviceData,
            static_cast<INT>(sourceData.m_sSourcePoint.w),
            sourceData.m_bySpinIndex,
            m_byFieldId);
    }
    break;
    default:
        appCrucial(_T("The source type %s not implemented yet!\n"), __ENUM_TO_STRING(EFermionSource, sourceData.m_eSourceType).c_str());
        break;
    }
}

BYTE* CFieldFermionKSU1::CopyDataOut(UINT& uiSize) const
{
    CLGComplex* toSave = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(Real) * m_uiSiteCount * 2);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        Real oneSite[2];
        oneSite[0] = static_cast<Real>(toSave[i].x);
        oneSite[1] = static_cast<Real>(toSave[i].y);
        memcpy(saveData + sizeof(Real) * i * 2, oneSite, sizeof(Real) * 2);
    }

    //appGetFileSystem()->WriteAllBytes(fileName.c_str(), saveData, uiSize);
    //free(saveData);
    free(toSave);
    return saveData;
}

BYTE* CFieldFermionKSU1::CopyDataOutFloat(UINT& uiSize) const
{
    CLGComplex* toSave = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(FLOAT) * m_uiSiteCount * 2);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        FLOAT oneSite[2];
        oneSite[0] = static_cast<FLOAT>(toSave[i].x);
        oneSite[1] = static_cast<FLOAT>(toSave[i].y);
        memcpy(saveData + sizeof(FLOAT) * i * 2, oneSite, sizeof(FLOAT) * 2);
    }

    //appGetFileSystem()->WriteAllBytes(fileName.c_str(), saveData, uiSize);
    //free(saveData);
    free(toSave);
    return saveData;
}

BYTE* CFieldFermionKSU1::CopyDataOutDouble(UINT& uiSize) const
{
    CLGComplex* toSave = (CLGComplex*)malloc(sizeof(CLGComplex) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(DOUBLE) * m_uiSiteCount * 2);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(CLGComplex) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        DOUBLE oneSite[2];
        oneSite[0] = static_cast<DOUBLE>(toSave[i].x);
        oneSite[1] = static_cast<DOUBLE>(toSave[i].y);
        memcpy(saveData + sizeof(DOUBLE) * i * 2, oneSite, sizeof(DOUBLE) * 2);
    }

    //appGetFileSystem()->WriteAllBytes(fileName.c_str(), saveData, uiSize);
    //free(saveData);
    free(toSave);
    return saveData;
}

TArray<CFieldFermion*> CFieldFermionKSU1::GetSourcesAtSiteFromPool(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const SSmallInt4& site) const
{
    TArray<CFieldFermion*> ret;
    ret.AddItem(dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId)));

    SFermionSource sourceData;
    sourceData.m_eSourceType = EFS_Point;
    sourceData.m_sSourcePoint = site;
    sourceData.m_bySpinIndex = 0;

    ret[0]->InitialAsSource(sourceData);

    if (NULL != appGetFermionSolver(m_byFieldId) 
            && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
    {
        ret[0]->m_fLength = ret[0]->Dot(ret[0]).x;
    }

    ret[0]->InverseD(gaugeNum, bosonNum, gaugeFields, pBoson);
    return ret;
}

void CFieldFermionKSU1::PrepareForHMCOnlyRandomize()
{
    preparethread;
    _kernelInitialFermionKSU1 << <block, threads >> > (
        m_pDeviceData,
        m_byFieldId,
        EFIT_RandomGaussian);
}

void CFieldFermionKSU1::PrepareForHMCNotRandomize(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson)
{
    D_MC(gaugeNum, bosonNum, gaugeFields, pBoson);
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================