//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3.cu
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

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionKSSU3)

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
_kernelDFermionKS(
    const deviceSU3Vector * __restrict__ pDeviceData,
    const deviceSU3 * __restrict__ pGauge,
    const SIndex * __restrict__ pGaugeMove,
    const SIndex * __restrict__ pFermionMove,
    const BYTE * __restrict__ pEtaTable,
    deviceSU3Vector* pResultData,
    Real f2am,
    BYTE byFieldId,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    CLGComplex cCoeff)
{
    intokernaldir;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
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
        const deviceSU3& x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

        //U(x,mu) phi(x+ mu)
        deviceSU3Vector u_phi_x_p_m = x_Gauge_element.MulVector(pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.MulReal(F(-1.0));
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceSU3Vector u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulVector(pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.Add(u_dagger_phi_x_m_m);
        }
        else
        {
            u_phi_x_p_m.Sub(u_dagger_phi_x_m_m);
        }
        u_phi_x_p_m.MulReal(eta_mu);
        result.Add(u_phi_x_p_m);
    }

    pResultData[uiSiteIndex].MulReal(f2am);
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
 * For some strange boundary condition
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionKSPlusEta(
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
    intokernaldir;

    deviceSU3Vector result = deviceSU3Vector::makeZeroSU3Vector();
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
        const deviceSU3& x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

        //U(x,mu) phi(x+ mu)
        deviceSU3Vector u_phi_x_p_m = x_Gauge_element.MulVector(pDeviceData[x_p_mu_Fermion.m_uiSiteIndex]);
        if (x_p_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.MulReal(F(-1.0) * eta_mu);
        }
        else
        {
            u_phi_x_p_m.MulReal(eta_mu);
        }

        //U^{dagger}(x-mu) phi(x-mu)
        deviceSU3Vector u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulVector(pDeviceData[x_m_mu_Fermion.m_uiSiteIndex]);
        u_dagger_phi_x_m_m.MulReal(eta_mu2);
        if (x_m_mu_Fermion.NeedToOpposite())
        {
            u_phi_x_p_m.Add(u_dagger_phi_x_m_m);
        }
        else
        {
            u_phi_x_p_m.Sub(u_dagger_phi_x_m_m);
        }
        //u_phi_x_p_m.MulReal(eta_mu);
        result.Add(u_phi_x_p_m);
    }

    pResultData[uiSiteIndex].MulReal(f2am);
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
_kernelDFermionKSForce(
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    const SIndex* __restrict__ pFermionMove,
    const BYTE* __restrict__ pEtaTable,
    const deviceSU3Vector* const* __restrict__ pFermionPointers,
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
            const deviceSU3Vector* phi_i = pFermionPointers[uiR];
            const deviceSU3Vector* phi_id = pFermionPointers[uiR + uiRational];

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

void CFieldFermionKSSU3::DOperatorKS(void* pTargetBuffer, const void * pBuffer,
    const void * pGaugeBuffer, Real f2am,
    UBOOL bDagger, EOperatorCoefficientType eOCT,
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceSU3Vector* pTarget = (deviceSU3Vector*)pTargetBuffer;
    const deviceSU3Vector* pSource = (const deviceSU3Vector*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    preparethread;
    if (m_bEachSiteEta)
    {
        _kernelDFermionKSPlusEta << <block, threads >> > (
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
    else
    {
        _kernelDFermionKS << <block, threads >> > (
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
}

void CFieldFermionKSSU3::DerivateDOperator(void* pForce, const void * pDphi, const void * pDDphi, const void * pGaugeBuffer) const
{
    appCrucial(_T("Do not call KS DerivateDOperator!"));
    _FAIL_EXIT;
}

/**
 * partial D_{st0} / partial omega
 * Make sure m_pMDNumerator and m_pRationalFieldPointers are filled
 */
void CFieldFermionKSSU3::DerivateD0(
    void* pForce, 
    const void* pGaugeBuffer) const
{
    preparethread;
    _kernelDFermionKSForce << <block, threads >> > (
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

#pragma region Staggered

/**
* Initialize
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialFermionKS(deviceSU3Vector* pDevicePtr, BYTE byFieldId, EFieldInitialType eInitialType)
{
    intokernalInt4;

    switch (eInitialType)
    {
    case EFIT_Zero:
    {
        pDevicePtr[uiSiteIndex] = deviceSU3Vector::makeZeroSU3Vector();
    }
    break;
    case EFIT_Identity:
    {
        pDevicePtr[uiSiteIndex] = deviceSU3Vector::makeOneSU3Vector();
    }
    break;
    case EFIT_RandomGaussian:
    {
        pDevicePtr[uiSiteIndex] = deviceSU3Vector::makeRandomGaussian(_deviceGetFatIndex(uiSiteIndex, 0));
    }
    break;
    case EFIT_RandomZ4:
    {
        pDevicePtr[uiSiteIndex] = deviceSU3Vector::makeRandomZ4(_deviceGetFatIndex(uiSiteIndex, 0));
    }
    break;
    default:
    {
        printf("Wilson Fermion Field cannot be initialized with this type!");
    }
    break;
    }
}


void CFieldFermionKSSU3::D_MD(const CField* pGauge)
{
    RationalApproximation(EFO_F_DDdagger, pGauge, &m_rMD);
}

void CFieldFermionKSSU3::D_MC(const CField* pGauge)
{
    RationalApproximation(EFO_F_DDdagger, pGauge, &m_rMC);
}

//void CFieldFermionKSSU3::D_EN(const CField* pGauge)
//{
//    RationalApproximation(EFO_F_DDdagger, pGauge, &m_rEN);
//}

/**
* generate phi by gaussian random.
* phi = (D^+D)^{1/4} phi
*/
void CFieldFermionKSSU3::PrepareForHMC(const CFieldGauge* pGauge)
{
    preparethread;
    _kernelInitialFermionKS << <block, threads >> > (
        m_pDeviceData,
        m_byFieldId,
        EFIT_RandomGaussian);

    D_MC(pGauge);

    if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
    {
        m_fLength = Dot(this).x;
    }

    //For KS, we generally use shifted solver, so do NOT cache the last result
}

/**
 * Use \sqrt{a} and b of rational A^{-1/2}
 */
UBOOL CFieldFermionKSSU3::CalculateForce(
    const CFieldGauge* pGauge,
    CFieldGauge* pForce,
    ESolverPhase ePhase) const
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    if (NULL == pForce || EFT_GaugeSU3 != pForce->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }

    TArray<CField*> phii;
    TArray<CFieldFermionKSSU3*> phiid;
    for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
    {
        CField* pPhi_i = dynamic_cast<CField*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
        phii.AddItem(pPhi_i);
        CFieldFermionKSSU3* pPhi_id = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
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
    solver->Solve(phii, shifts, this, pGauge, EFO_F_DDdagger);

    const UINT uiBufferSize = sizeof(deviceSU3Vector*) * 2 * m_rMD.m_uiDegree;
    deviceSU3Vector** hostPointers = (deviceSU3Vector**)appAlloca(uiBufferSize);
    for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
    {
        CFieldFermionKSSU3* phi_ks = dynamic_cast<CFieldFermionKSSU3*>(phii[i]);
        phi_ks->CopyTo(phiid[i]);
        phiid[i]->D0(pGauge);

        hostPointers[i] = phi_ks->m_pDeviceData;
        hostPointers[i + m_rMD.m_uiDegree] = phiid[i]->m_pDeviceData;
    }
    checkCudaErrors(cudaMemcpy(m_pRationalFieldPointers, hostPointers, uiBufferSize, cudaMemcpyHostToDevice));

    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3* pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);

    DerivateD0(pForceSU3->m_pDeviceData, pGaugeSU3->m_pDeviceData);
    //preparethread;
    //_kernelDFermionKSForce << <block, threads >> > (
    //    pGaugeSU3->m_pDeviceData,
    //    pForceSU3->m_pDeviceData,
    //    appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
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
_kernelPrintFermionKS(const deviceSU3Vector * __restrict__ pData)
{
    intokernal;

    printf("%d=(%1.2f %1.2fi, %1.2f %1.2fi, %1.2f %1.2fi))\n",
        uiSiteIndex,
        pData[uiSiteIndex].m_ve[0].x, pData[uiSiteIndex].m_ve[0].y,
        pData[uiSiteIndex].m_ve[1].x, pData[uiSiteIndex].m_ve[1].y,
        pData[uiSiteIndex].m_ve[2].x, pData[uiSiteIndex].m_ve[2].y
    );
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusFermionKS(deviceSU3Vector* pMe, const deviceSU3Vector* __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusFermionKS(deviceSU3Vector* pMe, const deviceSU3Vector* __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex].Sub(pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyComplexFermionKS(deviceSU3Vector* pMe, const deviceSU3Vector* __restrict__ pOther, CLGComplex a)
{
    intokernal;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex].MulCompC(a));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyRealFermionKS(deviceSU3Vector* pMe, const deviceSU3Vector* __restrict__ pOther, Real a)
{
    intokernal;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex].MulRealC(a));
}

__global__ void _CLG_LAUNCH_BOUND
_kernelDotFermionKS(const deviceSU3Vector* __restrict__ pMe, const deviceSU3Vector* __restrict__ pOther, CLGComplex* result)
{
    intokernal;
    result[uiSiteIndex] = pMe[uiSiteIndex].ConjugateDotC(pOther[uiSiteIndex]);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyComplexKS(deviceSU3Vector* pMe, CLGComplex a)
{
    intokernal;
    pMe[uiSiteIndex].MulComp(a);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyRealKS(deviceSU3Vector* pMe, Real a)
{
    intokernal;
    pMe[uiSiteIndex].MulReal(a);
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFermionKSConjugate(deviceSU3Vector* pDeviceData)
{
    intokernal;
    pDeviceData[uiSiteIndex].Conjugate();
}


__global__ void _CLG_LAUNCH_BOUND
_kernelMakePointSourceKS(deviceSU3Vector* pDeviceData, UINT uiDesiredSite, BYTE byColor)
{
    intokernal;
    if (uiSiteIndex == uiDesiredSite)
    {
        pDeviceData[uiSiteIndex] = deviceSU3Vector::makeOneSU3VectorColor(byColor);
    }
    else
    {
        pDeviceData[uiSiteIndex] = deviceSU3Vector::makeZeroSU3Vector();
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMakeWallSourceKS(deviceSU3Vector* pDeviceData, UINT uiDesiredT, BYTE byColor, BYTE byFieldID, UINT uiVolumn)
{
    intokernalInt4;
    const Real fDenominator = F(1.0) / uiVolumn;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldID][uiBigIdx];

    if (uiDesiredT == sSite4.w && !sIdx.IsDirichlet())
    {
        pDeviceData[uiSiteIndex] = deviceSU3Vector::makeOneSU3VectorColor(byColor);
        pDeviceData[uiSiteIndex].MulReal(fDenominator);
    }
    else
    {
        pDeviceData[uiSiteIndex] = deviceSU3Vector::makeZeroSU3Vector();
    }
}

#pragma endregion

CFieldFermionKSSU3::CFieldFermionKSSU3()
    : CFieldFermion()
    , m_bEachSiteEta(FALSE)
    , m_f2am(F(0.01))
    , m_pRationalFieldPointers(NULL)
    , m_pMDNumerator(NULL)
{
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount));
}

CFieldFermionKSSU3::~CFieldFermionKSSU3()
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
void CFieldFermionKSSU3::InitialField(EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialFermionKS << <block, threads >> > (m_pDeviceData, m_byFieldId, eInitialType);
}

void CFieldFermionKSSU3::Dagger()
{
    preparethread;
    _kernelFermionKSConjugate << <block, threads >> > (m_pDeviceData);
}

void CFieldFermionKSSU3::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFieldType)
{
    if (eFieldType != EFFT_CLGBin)
    {
        appCrucial(_T("CFieldFermionKSSU3::InitialFieldWithFile: Only support CLG Bin File\n"));
        return;
    }

    UINT uiSize = static_cast<UINT>(sizeof(Real) * 6 * m_uiSiteCount);
    BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
    InitialWithByte(data);
    free(data);
}

void CFieldFermionKSSU3::InitialWithByte(BYTE* byData)
{
    deviceSU3Vector* readData = (deviceSU3Vector*)malloc(sizeof(deviceSU3Vector) * m_uiSiteCount);
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        Real thisSite[24];
        memcpy(thisSite, byData + i * sizeof(Real) * 24, sizeof(Real) * 24);
        for (UINT k = 0; k < 3; ++k)
        {
            readData[i].m_ve[k] = _make_cuComplex(
                thisSite[2 * k],
                thisSite[2 * k + 1]);
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyHostToDevice));
    free(readData);
}

void CFieldFermionKSSU3::InitialOtherParameters(CParameters& params)
{
    params.FetchValueReal(_T("Mass"), m_f2am);
    if (m_f2am < F(0.00000001))
    {
        appCrucial(_T("CFieldFermionKSSU3: Mass is nearly 0!\n"));
    }

    INT iEachEta = 0;
    params.FetchValueINT(_T("EachSiteEta"), iEachEta);
    m_bEachSiteEta = (0 != iEachEta);

    TArray<Real> coeffs;
    params.FetchValueArrayReal(_T("MD"), coeffs);
    m_rMD.Initial(coeffs);

    params.FetchValueArrayReal(_T("MC"), coeffs);
    m_rMC.Initial(coeffs);

    //params.FetchValueArrayReal(_T("EN"), coeffs);
    //m_rEN.Initial(coeffs);

    if (NULL != m_pRationalFieldPointers)
    {
        checkCudaErrors(cudaFree(m_pRationalFieldPointers));
    }
    if (NULL != m_pMDNumerator)
    {
        checkCudaErrors(cudaFree(m_pMDNumerator));
    }
    checkCudaErrors(cudaMalloc((void**)&m_pRationalFieldPointers, sizeof(deviceSU3Vector*) * 2 * m_rMD.m_uiDegree));
    checkCudaErrors(cudaMalloc((void**)&m_pMDNumerator, sizeof(Real) * m_rMD.m_uiDegree));
    Real* hostNumerator = (Real*)appAlloca(sizeof(Real) * m_rMD.m_uiDegree);
    for (UINT i = 0; i < m_rMD.m_uiDegree; ++i)
    {
        hostNumerator[i] = m_rMD.m_lstA[i];
    }
    checkCudaErrors(cudaMemcpy(m_pMDNumerator, hostNumerator, sizeof(Real) * m_rMD.m_uiDegree, cudaMemcpyHostToDevice));

}

void CFieldFermionKSSU3::DebugPrintMe() const
{
    deviceSU3Vector* toprint = (deviceSU3Vector*)malloc(sizeof(deviceSU3Vector) * m_uiSiteCount);
    checkCudaErrors(cudaMemcpy(toprint, m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT uiSite = 0; uiSite < m_uiSiteCount; ++uiSite)
    {
        appGeneral(_T(" --- %d --- {%f %s %f I, %f %s %f I, %f %s %f I}\n"),
            uiSite,
            toprint[uiSite].m_ve[0].x,
            toprint[uiSite].m_ve[0].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(toprint[uiSite].m_ve[0].y),

            toprint[uiSite].m_ve[1].x,
            toprint[uiSite].m_ve[1].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(toprint[uiSite].m_ve[1].y),

            toprint[uiSite].m_ve[2].x,
            toprint[uiSite].m_ve[2].y > F(0.0) ? _T("+") : _T("-"),
            appAbs(toprint[uiSite].m_ve[2].y)
        );
    }

    appSafeFree(toprint);
}

void CFieldFermionKSSU3::CopyTo(CField* U) const
{
    if (NULL == U || EFT_FermionStaggeredSU3 != U->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSSU3 can only copy to CFieldFermionKSSU3!"));
        return;
    }

    CField::CopyTo(U);

    CFieldFermionKSSU3* pField = dynamic_cast<CFieldFermionKSSU3*>(U);
    checkCudaErrors(cudaMemcpy(pField->m_pDeviceData, m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
    pField->m_byFieldId = m_byFieldId;
    pField->m_f2am = m_f2am;
    pField->m_rMC = m_rMC;
    pField->m_rMD = m_rMD;
    pField->m_bEachSiteEta = m_bEachSiteEta;
    //pField->m_rEN = m_rEN;

    if (NULL != pField->m_pMDNumerator)
    {
        checkCudaErrors(cudaFree(pField->m_pMDNumerator));
    }
    if (NULL != pField->m_pRationalFieldPointers)
    {
        checkCudaErrors(cudaFree(pField->m_pRationalFieldPointers));
    }

    checkCudaErrors(cudaMalloc((void**)&pField->m_pRationalFieldPointers, sizeof(deviceSU3Vector*) * 2 * m_rMD.m_uiDegree));
    checkCudaErrors(cudaMalloc((void**)&pField->m_pMDNumerator, sizeof(Real) * m_rMD.m_uiDegree));
    checkCudaErrors(cudaMemcpy(pField->m_pMDNumerator, m_pMDNumerator, sizeof(Real) * m_rMD.m_uiDegree, cudaMemcpyDeviceToDevice));
}

void CFieldFermionKSSU3::AxpyPlus(const CField* x)
{
    if (NULL == x || EFT_FermionStaggeredSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSSU3 can only copy to CFieldFermionKSSU3!"));
        return;
    }
    const CFieldFermionKSSU3* pField = dynamic_cast<const CFieldFermionKSSU3*>(x);

    preparethread;
    _kernelAxpyPlusFermionKS << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldFermionKSSU3::AxpyMinus(const CField* x)
{
    if (NULL == x || EFT_FermionStaggeredSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSSU3 can only copy to CFieldFermionKSSU3!"));
        return;
    }
    const CFieldFermionKSSU3* pField = dynamic_cast<const CFieldFermionKSSU3*>(x);

    preparethread;
    _kernelAxpyMinusFermionKS << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldFermionKSSU3::Axpy(Real a, const CField* x)
{
    if (NULL == x || EFT_FermionStaggeredSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSSU3 can only copy to CFieldFermionKSSU3!"));
        return;
    }
    const CFieldFermionKSSU3* pField = dynamic_cast<const CFieldFermionKSSU3*>(x);

    preparethread;
    _kernelAxpyRealFermionKS << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

void CFieldFermionKSSU3::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || EFT_FermionStaggeredSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSSU3 can only copy to CFieldFermionKSSU3!"));
        return;
    }
    const CFieldFermionKSSU3* pField = dynamic_cast<const CFieldFermionKSSU3*>(x);

    preparethread;
    _kernelAxpyComplexFermionKS << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

CLGComplex CFieldFermionKSSU3::Dot(const CField* x) const
{
    if (NULL == x || EFT_FermionStaggeredSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSSU3 can only copy to CFieldFermionKSSU3!"));
        return _make_cuComplex(0, 0);
    }
    const CFieldFermionKSSU3* pField = dynamic_cast<const CFieldFermionKSSU3*>(x);

    preparethread;
    _kernelDotFermionKS << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, _D_ComplexThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

void CFieldFermionKSSU3::ScalarMultply(const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplyComplexKS << <block, threads >> > (m_pDeviceData, a);
}

void CFieldFermionKSSU3::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplyRealKS << <block, threads >> > (m_pDeviceData, a);
}

void CFieldFermionKSSU3::ApplyGamma(EGammaMatrix eGamma)
{
    appCrucial(_T("Not implemented yet...\n"));
}


//Kai should be part of D operator
void CFieldFermionKSSU3::D(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldFermionKSSU3* pPooled = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

void CFieldFermionKSSU3::D0(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldFermionKSSU3* pPooled = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    DOperatorKS(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, F(0.0),
        FALSE, EOCT_None, F(1.0), _onec);

    pPooled->Return();
}

UINT CFieldFermionKSSU3::TestAntiHermitian(const CFieldGauge* pGauge) const
{
    const UINT uiVolume = _HC_Volume;
    const UINT uiRealVolume = 3 * uiVolume;
    CLGComplex* matrixElement = (CLGComplex*)malloc(sizeof(CLGComplex) * uiRealVolume * uiRealVolume);
    deviceSU3Vector* hostData = (deviceSU3Vector*)malloc(sizeof(deviceSU3Vector) * uiVolume);
    CFieldFermionKSSU3* v = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    for (UINT i = 0; i < uiVolume; ++i)
    {
        const SSmallInt4 point = __hostSiteIndexToInt4(i);
        for (UINT j = 0; j < 3; ++j)
        {
            SFermionSource source;
            source.m_byColorIndex = static_cast<BYTE>(j);
            source.m_eSourceType = EFS_Point;
            source.m_sSourcePoint = point;
            v->InitialAsSource(source);
            v->D0(pGauge);

            checkCudaErrors(cudaMemcpy(hostData, v->m_pDeviceData, sizeof(deviceSU3Vector) * uiVolume, cudaMemcpyDeviceToHost));

            const UINT x = i * 3 + j;
            for (UINT k = 0; k < uiVolume; ++k)
            {
                matrixElement[(3 * k + 0) * uiRealVolume + x] = hostData[k].m_ve[0];
                matrixElement[(3 * k + 1) * uiRealVolume + x] = hostData[k].m_ve[1];
                matrixElement[(3 * k + 2) * uiRealVolume + x] = hostData[k].m_ve[2];
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
        const SSmallInt4 xSite = __hostSiteIndexToInt4(x / 3);
        const SSmallInt4 ySite = __hostSiteIndexToInt4(y / 3);
        const UINT daggerIdx = y * uiRealVolume + x;
        const BYTE cx = x % 3;
        const BYTE cy = y % 3;

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
void CFieldFermionKSSU3::Ddagger(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldFermionKSSU3* pPooled = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData,
        TRUE, eCoeffType, fRealCoeff, cCompCoeff);


    pPooled->Return();
}

void CFieldFermionKSSU3::DDdagger(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionKSSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }
    CFieldFermionKSSU3* pPooled = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData,
        TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //why only apply coeff in the next step?
    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

UBOOL CFieldFermionKSSU3::InverseD(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    //Find a solver to solve me.
    return appGetFermionSolver(m_byFieldId)->Solve(this, /*this is const*/this, pFieldSU3, EFO_F_D);
}

UBOOL CFieldFermionKSSU3::InverseDdagger(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    //Find a solver to solve me.
    return appGetFermionSolver(m_byFieldId)->Solve(this, /*this is const*/this, pFieldSU3, EFO_F_Ddagger);
}

UBOOL CFieldFermionKSSU3::InverseDDdagger(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    //Find a solver to solve me.
    return appGetFermionSolver(m_byFieldId)->Solve(this, /*this is const*/this, pFieldSU3, EFO_F_DDdagger);
}


void CFieldFermionKSSU3::InitialAsSource(const SFermionSource& sourceData)
{
    const UINT uiSiteIndex = _hostGetSiteIndex(sourceData.m_sSourcePoint);
    switch (sourceData.m_eSourceType)
    {
    case EFS_Point:
    {
        preparethread;
        _kernelMakePointSourceKS << <block, threads >> > (m_pDeviceData, uiSiteIndex, sourceData.m_byColorIndex);
    }
    break;
    case EFS_Wall:
    {
        preparethread;
        _kernelMakeWallSourceKS << <block, threads >> > (
            m_pDeviceData,
            sourceData.m_sSourcePoint.w,
            sourceData.m_byColorIndex,
            m_byFieldId,
            appGetLattice()->m_pIndexCache->m_uiSiteXYZ);
    }
    break;
    default:
        appCrucial(_T("The source type %s not implemented yet!\n"), __ENUM_TO_STRING(EFermionSource, sourceData.m_eSourceType).c_str());
        break;
    }
}

void CFieldFermionKSSU3::SetMass(Real f2am)
{
    m_f2am = f2am;
}

void CFieldFermionKSSU3::SaveToFile(const CCString& fileName) const
{
    UINT uiSize = 0;
    BYTE* saveData = CopyDataOut(uiSize);
    appGetFileSystem()->WriteAllBytes(fileName.c_str(), saveData, uiSize);
    free(saveData);
}

BYTE* CFieldFermionKSSU3::CopyDataOut(UINT& uiSize) const
{
    deviceSU3Vector* toSave = (deviceSU3Vector*)malloc(sizeof(deviceSU3Vector) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(Real) * m_uiSiteCount * 6);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceSU3Vector) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        Real oneSite[6];
        for (UINT k = 0; k < 3; ++k)
        {
            oneSite[2 * k] = static_cast<Real>(toSave[i].m_ve[k].x);
            oneSite[2 * k + 1] = static_cast<Real>(toSave[i].m_ve[k].y);
        }
        memcpy(saveData + sizeof(Real) * i * 6, oneSite, sizeof(Real) * 6);
    }

    //appGetFileSystem()->WriteAllBytes(fileName.c_str(), saveData, uiSize);
    //free(saveData);
    free(toSave);
    return saveData;
}

TArray<CFieldFermion*> CFieldFermionKSSU3::GetSourcesAtSiteFromPool(const class CFieldGauge* pGauge, const SSmallInt4& site) const
{
    TArray<CFieldFermion*> ret;
    for (UINT j = 0; j < 3; ++j)
    {
        ret.AddItem(dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId)));
        if (NULL == ret[j])
        {
            appCrucial(_T("GetSourcesAtSiteFromPool failed!\n"));
            _FAIL_EXIT;
        }
    }

    for (BYTE c = 0; c < 3; ++c)
    {
        SFermionSource sourceData;
        sourceData.m_eSourceType = EFS_Point;
        sourceData.m_sSourcePoint = site;
        sourceData.m_byColorIndex = c;
        sourceData.m_bySpinIndex = 0;

        ret[c]->InitialAsSource(sourceData);

        if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
        {
            ret[c]->m_fLength = ret[c]->Dot(ret[c]).x;
        }
        ret[c]->InverseD(pGauge);
    }
    return ret;
}

CCString CFieldFermionKSSU3::GetInfos(const CCString& tab) const
{
    CCString sRet = tab + _T("Name : CFieldFermionKSSU3\n");
    sRet = sRet + tab + _T("Mass (2am) : ") + appFloatToString(m_f2am) + _T("\n");
    sRet = sRet + tab + _T("MD Rational (c) : ") + appFloatToString(m_rMD.m_fC) + _T("\n");
    sRet = sRet + tab + _T("MC Rational (c) : ") + appFloatToString(m_rMC.m_fC) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================