//=============================================================================
// FILENAME : CFieldFermionWilsonSquareSU3.cu
// 
// DESCRIPTION:
// This is the device implementations of Wilson fermion
//
// This implementation assumes SU3 and square lattice
//
// REVISION:
//  [12/27/2018 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CFieldFermionWilsonSquareSU3)

#pragma region DOperator

#pragma region kernel

/**
* Dw phi(x) = phi(x) - kai sum _mu (1-gamma _mu) U(x,mu) phi(x+ mu) + (1+gamma _mu) U^{dagger}(x-mu) phi(x-mu)
* U act on su3
* gamma act on spinor
*
* If bDagger, it is gamma5, D, gamma5
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDFermionWilsonSquareSU3(
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
    intokernaldir;

    //const SSmallInt4 test = __deviceSiteIndexToInt4(uiSiteIndex);

    const gammaMatrix & gamma5 = __chiralGamma[GAMMA5];
    deviceWilsonVectorSU3 result = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    pResultData[uiSiteIndex] = pDeviceData[uiSiteIndex];
    if (bDDagger)
    {
        pResultData[uiSiteIndex] = gamma5.MulWilsonC(pResultData[uiSiteIndex]);
    }

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        const gammaMatrix & gammaMu = __chiralGamma[GAMMA1 + idir];

        //x, mu
        const UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex & x_m_mu_Gauge = pGaugeMove[linkIndex];

        const SIndex & x_p_mu_Fermion = pFermionMove[2 * linkIndex];
        const SIndex & x_m_mu_Fermion = pFermionMove[2 * linkIndex + 1];

        //=====================
        //check links
        //const SSmallInt4 site1 = __deviceSiteIndexToInt4(uiSiteIndex);
        //const SSmallInt4 site2 = __deviceSiteIndexToInt4(x_p_mu_Fermion.m_uiSiteIndex);
        //const SSmallInt4 site3 = __deviceSiteIndexToInt4(x_m_mu_Fermion.m_uiSiteIndex);

        //if (site1.x == 0 && site1.y == 3 && site1.z == 0 && site1.w == 0
        //    && site2.x == 3 && site2.y == 0 && site2.z == 0 && site2.w == 0)
        //{
        //    printf("link1 = %d\n", linkIndex);
        //}

        //if (site3.x == 0 && site3.y == 3 && site3.z == 0 && site3.w == 0
        //    && site1.x == 3 && site1.y == 0 && site1.z == 0 && site1.w == 0)
        //{
        //    printf("link2 = %d\n", _deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir));
        //}

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        const deviceSU3 & x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        if (x_m_mu_Gauge.NeedToDagger())
        {
            x_m_mu_Gauge_element.Dagger();
        }

        deviceWilsonVectorSU3 x_p_mu_Fermion_element = pDeviceData[x_p_mu_Fermion.m_uiSiteIndex];
        deviceWilsonVectorSU3 x_m_mu_Fermion_element = pDeviceData[x_m_mu_Fermion.m_uiSiteIndex];

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
    pResultData[uiSiteIndex].Sub(result);

    if (bDDagger)
    {
        pResultData[uiSiteIndex] = gamma5.MulWilsonC(pResultData[uiSiteIndex]);
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
* The output is on a gauge field
* Therefor cannot make together with _kernelDWilson
*
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelDWilsonForceSU3(
    const deviceWilsonVectorSU3* __restrict__ pInverseD,
    const deviceWilsonVectorSU3* __restrict__ pInverseDDdagger,
    const deviceSU3* __restrict__ pGauge,
    const SIndex* __restrict__ pFermionMove,
    deviceSU3* pForce,
    Real fKai,
    BYTE byFieldId)
{
    intokernaldir;

    const deviceWilsonVectorSU3 & x_Left  = pInverseDDdagger[uiSiteIndex];
    const deviceWilsonVectorSU3 & x_Right = pInverseD[uiSiteIndex];

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        const gammaMatrix & gammaMu = __chiralGamma[GAMMA1 + idir];

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        const SIndex& x_p_mu_Fermion = pFermionMove[linkIndex * 2]; // __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, (idir + 1));

        const deviceWilsonVectorSU3& x_p_mu_Right = pInverseD[x_p_mu_Fermion.m_uiSiteIndex];
        const deviceWilsonVectorSU3& x_p_mu_Left = pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex];

        const deviceSU3& x_Gauge_element = pGauge[linkIndex];

        //if (961 == linkIndex)
        //{
        //    const SSmallInt4 s0 = __deviceSiteIndexToInt4(240);
        //    const SSmallInt4 s1 = __deviceSiteIndexToInt4(uiSiteIndex);
        //    const SSmallInt4 s2 = __deviceSiteIndexToInt4(x_p_mu_Fermion.m_uiSiteIndex);
        //    printf("s0 = %d, %d, %d, %d, s1 = %d, %d, %d, %d, s2 = %d, %d, %d, %d;\n", 
        //        s0.x, s0.y, s0.z, s0.w, s1.x, s1.y, s1.z, s1.w, s2.x, s2.y, s2.z, s2.w);
        //}

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

#pragma endregion

void CFieldFermionWilsonSquareSU3::DOperator(void* pTargetBuffer, const void* pBuffer, 
    const void* pGaugeBuffer, 
    UBOOL bDagger, EOperatorCoefficientType eOCT, 
    Real fRealCoeff, const CLGComplex& cCmpCoeff) const
{
    deviceWilsonVectorSU3* pTarget = (deviceWilsonVectorSU3*)pTargetBuffer;
    const deviceWilsonVectorSU3* pSource = (deviceWilsonVectorSU3*)pBuffer;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;

    preparethread;
    _kernelDFermionWilsonSquareSU3 << <block, threads >> > (
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

void CFieldFermionWilsonSquareSU3::DerivateDOperator(void* pForce, const void* pDphi, const void* pDDphi, const void* pGaugeBuffer) const
{
    deviceSU3* pForceSU3 = (deviceSU3*)pForce;
    const deviceSU3* pGauge = (const deviceSU3*)pGaugeBuffer;
    const deviceWilsonVectorSU3* pDphiBuffer = (deviceWilsonVectorSU3*)pDphi;
    const deviceWilsonVectorSU3* pDDphiBuffer = (deviceWilsonVectorSU3*)pDDphi;

    preparethread;
    _kernelDWilsonForceSU3 << <block, threads >> > (
        pDphiBuffer,
        pDDphiBuffer,
        pGauge,
        appGetLattice()->m_pIndexCache->m_pFermionMoveCache[m_byFieldId],
        pForceSU3,
        m_fKai, m_byFieldId);
}

#pragma endregion

#pragma region Kernel

__global__ void _CLG_LAUNCH_BOUND
_kernelPrintFermionWilsonSquareSU3(const deviceWilsonVectorSU3 * __restrict__ pData)
{
    intokernal;

    printf("%d=((%1.2f %1.2fi, %1.2f %1.2fi, %1.2f %1.2fi),(%1.2f %1.2fi, %1.2f %1.2fi, %1.2f %1.2fi),(%1.2f %1.2fi, %1.2f %1.2fi, %1.2f %1.2fi),(%1.2f %1.2fi, %1.2f %1.2fi, %1.2f %1.2fi))\n",
        uiSiteIndex,
        pData[uiSiteIndex].m_d[0].m_ve[0].x, pData[uiSiteIndex].m_d[0].m_ve[0].y,
        pData[uiSiteIndex].m_d[0].m_ve[1].x, pData[uiSiteIndex].m_d[0].m_ve[1].y,
        pData[uiSiteIndex].m_d[0].m_ve[2].x, pData[uiSiteIndex].m_d[0].m_ve[2].y,

        pData[uiSiteIndex].m_d[1].m_ve[0].x, pData[uiSiteIndex].m_d[1].m_ve[0].y,
        pData[uiSiteIndex].m_d[1].m_ve[1].x, pData[uiSiteIndex].m_d[1].m_ve[1].y,
        pData[uiSiteIndex].m_d[1].m_ve[2].x, pData[uiSiteIndex].m_d[1].m_ve[2].y,

        pData[uiSiteIndex].m_d[2].m_ve[0].x, pData[uiSiteIndex].m_d[2].m_ve[0].y,
        pData[uiSiteIndex].m_d[2].m_ve[1].x, pData[uiSiteIndex].m_d[2].m_ve[1].y,
        pData[uiSiteIndex].m_d[2].m_ve[2].x, pData[uiSiteIndex].m_d[2].m_ve[2].y,

        pData[uiSiteIndex].m_d[3].m_ve[0].x, pData[uiSiteIndex].m_d[3].m_ve[0].y,
        pData[uiSiteIndex].m_d[3].m_ve[1].x, pData[uiSiteIndex].m_d[3].m_ve[1].y,
        pData[uiSiteIndex].m_d[3].m_ve[2].x, pData[uiSiteIndex].m_d[3].m_ve[2].y
    );
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyPlusFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 * pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex]);
}

//__global__ void _CLG_LAUNCH_BOUND_(_QUICK_AXPY_BLOCK)
//_kernelAxpyPlusQuick(
//    deviceWilsonVectorSU3 * pMe,
//    const deviceWilsonVectorSU3 * __restrict__ pOther)
//{
//    intokernalE(24);
//    pMe[uiSiteIndex].m_rme[elementIdx] += pOther[uiSiteIndex].m_rme[elementIdx];
//}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyMinusFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 * pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex].Sub(pOther[uiSiteIndex]);
}

//__global__ void _CLG_LAUNCH_BOUND_(_QUICK_AXPY_BLOCK)
//_kernelAxpyMinusQuick(
//    deviceWilsonVectorSU3 * pMe,
//    const deviceWilsonVectorSU3 * __restrict__ pOther)
//{
//    intokernalE(24);
//    pMe[uiSiteIndex].m_rme[elementIdx] -= pOther[uiSiteIndex].m_rme[elementIdx];
//}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyComplexFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 * pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther, CLGComplex a)
{
    intokernal;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex].MulCompC(a));
}

//__global__ void _CLG_LAUNCH_BOUND_(_QUICK_AXPY_BLOCK)
//_kernelAxpyComplexQuick(
//    deviceWilsonVectorSU3 * pMe,
//    const deviceWilsonVectorSU3 * __restrict__ pOther, CLGComplex a)
//{
//    intokernalE(12);
//    pMe[uiSiteIndex].m_me[elementIdx] = 
//        _cuCaddf(pMe[uiSiteIndex].m_me[elementIdx], 
//            _cuCmulf(a, pOther[uiSiteIndex].m_me[elementIdx]));
//}

__global__ void _CLG_LAUNCH_BOUND
_kernelAxpyRealFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 * pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther, Real a)
{
    intokernal;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex].MulRealC(a));
}

//__global__ void _CLG_LAUNCH_BOUND_(_QUICK_AXPY_BLOCK)
//_kernelAxpyRealQuick(
//    deviceWilsonVectorSU3 * pMe,
//    const deviceWilsonVectorSU3 * __restrict__ pOther, Real a)
//{
//    intokernalE(24);
//    pMe[uiSiteIndex].m_rme[elementIdx] += pOther[uiSiteIndex].m_rme[elementIdx] * a;
//}

__global__ void _CLG_LAUNCH_BOUND
_kernelDotFermionWilsonSquareSU3(
    const deviceWilsonVectorSU3 * __restrict__ pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther,
#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex* result
#else
    CLGComplex * result
#endif
)
{
    intokernal;
#if !_CLG_DOUBLEFLOAT
    result[uiSiteIndex] = _cToDouble(pMe[uiSiteIndex].ConjugateDotC(pOther[uiSiteIndex]));
#else
    result[uiSiteIndex] = pMe[uiSiteIndex].ConjugateDotC(pOther[uiSiteIndex]);
#endif
}

//__global__ void _CLG_LAUNCH_BOUND_(_QUICK_AXPY_BLOCK)
//_kernelDotQuick(
//    const deviceWilsonVectorSU3 * __restrict__ pMe,
//    const deviceWilsonVectorSU3 * __restrict__ pOther, Real* x, Real* y)
//{
//    intokernalE(12);
//    CLGComplex toAdd = _cuCmulf(_cuConjf(pMe[uiSiteIndex].m_me[elementIdx]), pOther[uiSiteIndex].m_me[elementIdx]);
//    atomicAdd(x, toAdd.x);
//    atomicAdd(y, toAdd.y);
//}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyComplex(
    deviceWilsonVectorSU3 * pMe, 
    CLGComplex a)
{
    intokernal;
    pMe[uiSiteIndex].MulComp(a);
}

//__global__ void _CLG_LAUNCH_BOUND_(_QUICK_AXPY_BLOCK)
//_kernelScalarQuick(
//    deviceWilsonVectorSU3 * pMe,
//    CLGComplex a)
//{
//    intokernalE(12);
//    pMe[uiSiteIndex].m_me[elementIdx] = _cuCmulf(pMe[uiSiteIndex].m_me[elementIdx], a);
//}

__global__ void _CLG_LAUNCH_BOUND
_kernelScalarMultiplyReal(
    deviceWilsonVectorSU3 * pMe, 
    Real a)
{
    intokernal;
    pMe[uiSiteIndex].MulReal(a);
}

//__global__ void _CLG_LAUNCH_BOUND_(_QUICK_AXPY_BLOCK)
//_kernelScalarQuickReal(
//    deviceWilsonVectorSU3 * pMe,
//    Real a)
//{
//    intokernalE(24);
//    pMe[uiSiteIndex].m_rme[elementIdx] *= a;
//}

/**
* Initialize
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 *pDevicePtr, 
    BYTE byFieldId,
    EFieldInitialType eInitialType)
{
    intokernalInt4;

    switch (eInitialType)
    {
    case EFIT_Zero:
    {
        pDevicePtr[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    }
    break;
    case EFIT_Identity:
    {
        pDevicePtr[uiSiteIndex] = deviceWilsonVectorSU3::makeOneWilsonVectorSU3();
    }
    break;
    case EFIT_RandomGaussian:
    {
        pDevicePtr[uiSiteIndex] = deviceWilsonVectorSU3::makeRandomGaussian(_deviceGetFatIndex(uiSiteIndex, 0));
    }
    break;
    case EFIT_RandomZ4:
    {
        pDevicePtr[uiSiteIndex] = deviceWilsonVectorSU3::makeRandomZ4(_deviceGetFatIndex(uiSiteIndex, 0));
    }
    break;
    default:
    {
        printf("Wilson Fermion Field cannot be initialized with this type!");
    }
    break;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelFermionWilsonSquareSU3Conjugate(deviceWilsonVectorSU3* pDevicePtr)
{
    intokernal;
    pDevicePtr[uiSiteIndex].Conjugate();
}

__global__ void _CLG_LAUNCH_BOUND
_kernelApplyGammaSU3(deviceWilsonVectorSU3* pDeviceData, UINT uiGamma)
{
    intokernal;
    pDeviceData[uiSiteIndex] = __chiralGamma[uiGamma].MulWilsonC(pDeviceData[uiSiteIndex]);
}


__global__ void _CLG_LAUNCH_BOUND
_kernelMakePointSource(deviceWilsonVectorSU3* pDeviceData, UINT uiDesiredSite, BYTE bySpin, BYTE byColor)
{
    intokernal;
    if (uiSiteIndex == uiDesiredSite)
    {
        pDeviceData[uiSiteIndex] = deviceWilsonVectorSU3::makeOneWilsonVectorSU3SpinColor(bySpin, byColor);
    }
    else
    {
        pDeviceData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMakePointSourceOne(deviceWilsonVectorSU3* pDeviceData, UINT uiDesiredSite)
{
    intokernal;
    if (uiSiteIndex == uiDesiredSite)
    {
        pDeviceData[uiSiteIndex] = deviceWilsonVectorSU3::makeOneWilsonVectorSU3();
    }
    else
    {
        pDeviceData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelMakeWallSource(deviceWilsonVectorSU3* pDeviceData, UINT uiDesiredT, BYTE bySpin, BYTE byColor, BYTE byFieldID, UINT uiVolumn)
{
    intokernalInt4;
    const Real fDenominator = F(1.0) / uiVolumn;
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex sIdx = __idx->m_pDeviceIndexPositionToSIndex[byFieldID][uiBigIdx];

    if (uiDesiredT == sSite4.w && !sIdx.IsDirichlet())
    {
        pDeviceData[uiSiteIndex] = deviceWilsonVectorSU3::makeOneWilsonVectorSU3SpinColor(bySpin, byColor);
        pDeviceData[uiSiteIndex].MulReal(fDenominator);
    }
    else
    {
        pDeviceData[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    }
}

#pragma endregion

CFieldFermionWilsonSquareSU3::CFieldFermionWilsonSquareSU3() 
    : CFieldFermion()
    , m_fKai(F(0.125))
    , m_tmpBuffer2(NULL)
{
    checkCudaErrors(__cudaMalloc((void**)&m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount));
    checkCudaErrors(cudaMalloc((void**)&m_tmpBuffer2, sizeof(Real) * 2));
}

CFieldFermionWilsonSquareSU3::~CFieldFermionWilsonSquareSU3()
{
    checkCudaErrors(__cudaFree(m_pDeviceData));
    checkCudaErrors(cudaFree(m_tmpBuffer2));
}

/**
*
*/
void CFieldFermionWilsonSquareSU3::InitialField(EFieldInitialType eInitialType)
{
    preparethread;

    _kernelInitialFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceData, m_byFieldId, eInitialType);
}

void CFieldFermionWilsonSquareSU3::Dagger()
{
    preparethread;

    _kernelFermionWilsonSquareSU3Conjugate << <block, threads >> > (m_pDeviceData);
}

void CFieldFermionWilsonSquareSU3::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFieldType)
{
    if (eFieldType != EFFT_CLGBin)
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3::InitialFieldWithFile: Only support CLG Bin File\n"));
        return;
    }

    UINT uiSize = static_cast<UINT>(sizeof(Real) * 24 * m_uiSiteCount);
    BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
    InitialWithByte(data);
    free(data);
}

void CFieldFermionWilsonSquareSU3::InitialWithByte(BYTE* byData)
{
    deviceWilsonVectorSU3* readData = (deviceWilsonVectorSU3*)malloc(sizeof(deviceWilsonVectorSU3) * m_uiSiteCount);
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        Real thisSite[24];
        memcpy(thisSite, byData + i * sizeof(Real) * 24, sizeof(Real) * 24);
        for (UINT j = 0; j < 4; ++j)
        {
            for (UINT k = 0; k < 3; ++k)
            {
                readData[i].m_d[j].m_ve[k] = _make_cuComplex(
                    thisSite[2 * (j * 3 + k)],
                    thisSite[2 * (j * 3 + k) + 1]);
            }
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount, cudaMemcpyHostToDevice));
    free(readData);
}

void CFieldFermionWilsonSquareSU3::InitialOtherParameters(CParameters& params)
{
    CFieldFermion::InitialOtherParameters(params);

    params.FetchValueReal(_T("Hopping"), m_fKai);
    if (m_fKai < F(0.00000001))
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3: Kai is nearly 0, such that Dphi \approx phi! This will cause problem!\n"));
    }
    CCommonData::m_fKai = m_fKai;

    //INT iEvenFieldId = -1;
    //params.FetchValueINT(_T("EvenFieldId"), iEvenFieldId);
    //if (iEvenFieldId > 0)
    //{
    //    m_byEvenFieldId = static_cast<SBYTE>(iEvenFieldId);
    //}
}

void CFieldFermionWilsonSquareSU3::DebugPrintMe() const
{
    preparethread;
    _kernelPrintFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceData);
}

void CFieldFermionWilsonSquareSU3::CopyTo(CField* U) const
{
    if (NULL == U || EFT_FermionWilsonSquareSU3 != U->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }

    CField::CopyTo(U);

    CFieldFermionWilsonSquareSU3 * pField = dynamic_cast<CFieldFermionWilsonSquareSU3*>(U);
    checkCudaErrors(cudaMemcpy(pField->m_pDeviceData, m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
    pField->m_byFieldId = m_byFieldId;
    pField->m_fKai = m_fKai;
}

#if 0
#pragma region Quick Axpys (Although this is faster in debug mode, this is not fast in release Mode)

void CFieldFermionWilsonSquareSU3::AxpyPlus1(const CField* x)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3 * pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethreadE(24);
    _kernelAxpyPlusQuick << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldFermionWilsonSquareSU3::AxpyMinus1(const CField* x)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3 * pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethreadE(24);
    _kernelAxpyMinusQuick << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldFermionWilsonSquareSU3::Axpy1(Real a, const CField* x)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3 * pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethreadE(24);
    _kernelAxpyRealQuick << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

void CFieldFermionWilsonSquareSU3::Axpy1(const CLGComplex& a, const CField* x)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3 * pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethreadE(12);
    _kernelAxpyComplexQuick << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

void CFieldFermionWilsonSquareSU3::ScalarMultply1(const CLGComplex& a)
{
    preparethreadE(12);
    _kernelScalarQuick << <block, threads >> >(m_pDeviceData, a);
}

void CFieldFermionWilsonSquareSU3::ScalarMultply1(Real a)
{
    preparethreadE(24);
    _kernelScalarQuickReal << <block, threads >> >(m_pDeviceData, a);
}

CLGComplex CFieldFermionWilsonSquareSU3::Dot1(const CField* x) const
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return _make_cuComplex(0, 0);
    }
    const CFieldFermionWilsonSquareSU3 * pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    Real tmp[] = { F(0.0), F(0.0) };
    checkCudaErrors(cudaMemcpy(m_tmpBuffer2, tmp, sizeof(Real) * 2, cudaMemcpyHostToDevice));
    preparethreadE(12);
    _kernelDotQuick << <block, threads >> > 
        (m_pDeviceData, pField->m_pDeviceData, m_tmpBuffer2, m_tmpBuffer2 + 1);
    checkCudaErrors(cudaMemcpy(tmp, m_tmpBuffer2, sizeof(Real) * 2, cudaMemcpyDeviceToHost));

    return _make_cuComplex(tmp[0], tmp[1]);
}

#pragma endregion
#endif

void CFieldFermionWilsonSquareSU3::AxpyPlus(const CField* x)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3 * pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethread;
    _kernelAxpyPlusFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldFermionWilsonSquareSU3::AxpyMinus(const CField* x)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3 * pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethread;
    _kernelAxpyMinusFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData);
}

void CFieldFermionWilsonSquareSU3::Axpy(Real a, const CField* x)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3 * pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethread;
    _kernelAxpyRealFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

void CFieldFermionWilsonSquareSU3::Axpy(const CLGComplex& a, const CField* x)
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return;
    }
    const CFieldFermionWilsonSquareSU3 * pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethread;
    _kernelAxpyComplexFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, a);
}

cuDoubleComplex CFieldFermionWilsonSquareSU3::Dot(const CField* x) const
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return make_cuDoubleComplex(0,0);
    }
    const CFieldFermionWilsonSquareSU3 * pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);
    preparethread;
    _kernelDotFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, _D_ComplexThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

void CFieldFermionWilsonSquareSU3::ScalarMultply(const CLGComplex& a)
{
    preparethread;
    _kernelScalarMultiplyComplex << <block, threads >> >(m_pDeviceData, a);
}

void CFieldFermionWilsonSquareSU3::ScalarMultply(Real a)
{
    preparethread;
    _kernelScalarMultiplyReal << <block, threads >> >(m_pDeviceData, a);
}

void CFieldFermionWilsonSquareSU3::ApplyGamma(EGammaMatrix eGamma)
{
    preparethread;
    _kernelApplyGammaSU3 << <block, threads >> >(m_pDeviceData, static_cast<UINT>(eGamma));
}
/**
* generate phi by gaussian random.
* phi = D phi
*/
void CFieldFermionWilsonSquareSU3::PrepareForHMCS(const CFieldGauge* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldFermionWilsonSquareSU3* pPooled = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    preparethread;
    _kernelInitialFermionWilsonSquareSU3 << <block, threads >> > (
        pPooled->m_pDeviceData,
        m_byFieldId,
        EFIT_RandomGaussian);

    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData, 
        FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));

    pPooled->Return();

    if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
    {
        m_fLength = Dot(this).x;
    }
    //cache a inverse DDdagger field
    if (CCommonData::m_bStoreLastSolution)
    {
        CFieldCache* pCache = appGetLattice()->m_pFieldCache;
        CField* pField = pCache->GetCachedField(CFieldCache::CachedInverseDDdaggerField);
        if (NULL == pField)
        {
            pField = GetCopy();
            pCache->CacheField(CFieldCache::CachedInverseDDdaggerField, pField);
        }
        else
        {
            CopyTo(pField);
        }
    }
}

//Kai should be part of D operator
void CFieldFermionWilsonSquareSU3::DS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldFermionWilsonSquareSU3* pPooled = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

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

//Kai should be part of D operator
void CFieldFermionWilsonSquareSU3::DdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldFermionWilsonSquareSU3* pPooled = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    checkCudaErrors(cudaMemcpy(pPooled->m_pDeviceData, m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

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

void CFieldFermionWilsonSquareSU3::DDdaggerS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    Real fRealCoeff = fCoeffReal;
    const CLGComplex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }
    CFieldFermionWilsonSquareSU3* pPooled = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData,
        TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //why only apply coeff in the next step?
    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

void CFieldFermionWilsonSquareSU3::DDS(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
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
    CFieldFermionWilsonSquareSU3* pPooled = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    DOperator(pPooled->m_pDeviceData, m_pDeviceData, pFieldSU3->m_pDeviceData,
        FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //why only apply coeff in the next step?
    DOperator(m_pDeviceData, pPooled->m_pDeviceData, pFieldSU3->m_pDeviceData,
        FALSE, eCoeffType, fRealCoeff, cCompCoeff);

    pPooled->Return();
}

UBOOL CFieldFermionWilsonSquareSU3::CalculateForceS(
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

    const CFieldGaugeSU3 * pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3 * pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);

    CField * pDDaggerPhi = appGetLattice()->GetPooledFieldById(m_byFieldId);
    CField * pDPhi = appGetLattice()->GetPooledFieldById(m_byFieldId);
    CField * pCachedField = CCommonData::m_bStoreLastSolution ? 
        appGetLattice()->m_pFieldCache->GetCachedField(CFieldCache::CachedInverseDDdaggerField)
       : NULL;

    if (NULL == pDDaggerPhi || EFT_FermionWilsonSquareSU3 != pDDaggerPhi->GetFieldType()
     || NULL == pDPhi || EFT_FermionWilsonSquareSU3 != pDPhi->GetFieldType())
    {
        appCrucial(_T("Pooled field not found!\n"));
        if (NULL != pDDaggerPhi)
        {
            pDDaggerPhi->Return();
        }
        if (NULL != pDPhi)
        {
            pDPhi->Return();
        }
        return FALSE;
    }
    CFieldFermionWilsonSquareSU3* pDDaggerPhiWilson = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pDDaggerPhi);
    CFieldFermionWilsonSquareSU3* pDPhiWilson = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pDPhi);
    //if (!pDDaggerPhiWilson->InverseDDdagger(pGaugeSU3))

    //if (m_byEvenFieldId > 0)
    //{
    //    CopyTo(pDDaggerPhiWilson);
    //    pDDaggerPhiWilson->InverseDDdagger(pGaugeSU3);
    //}
    //else
    //{
    TArray<const CFieldGauge*> gauge;
    gauge.AddItem(pGaugeSU3);
        if (!appGetFermionSolver(m_byFieldId)->Solve(
            pDDaggerPhiWilson, this, 1, 0, gauge.GetData(), NULL,
            EFO_F_DDdagger, ePhase, pCachedField))
        {
            appCrucial(_T("Sparse Linear Solver failed...\n"));
            pDDaggerPhi->Return();
            pDPhi->Return();
            return FALSE;
        }
    //}

    //phi 2 = D^{-1}phi = D+ (DD+)^{-1} phi
    //It is faster to calcuate D+ phi2 then D^{-1} phi
    pDDaggerPhiWilson->CopyTo(pDPhiWilson);
    if (NULL != pCachedField)
    {
        //The gauge field is changing slowly, and D depends only on gauge, also change slowly
        //Use the last solution as start point will accelerate the solver, so we cache it
        pDDaggerPhiWilson->CopyTo(pCachedField);
    }
    pDPhiWilson->DdaggerS(pGaugeSU3);
    DerivateDOperator(
        pForceSU3->m_pDeviceData, 
        pDPhiWilson->m_pDeviceData, 
        pDDaggerPhiWilson->m_pDeviceData, 
        pGaugeSU3->m_pDeviceData);

    pDDaggerPhi->Return();
    pDPhi->Return();

    return TRUE;
}

void CFieldFermionWilsonSquareSU3::InitialAsSource(const SFermionSource& sourceData)
{
    const UINT uiSiteIndex = _hostGetSiteIndex(sourceData.m_sSourcePoint);
    switch (sourceData.m_eSourceType)
    {
    case EFS_Point:
    {
        preparethread;
        if (sourceData.m_byColorIndex >= 3)
        {
            _kernelMakePointSourceOne << <block, threads >> > (m_pDeviceData, uiSiteIndex);
        }
        else
        {
            _kernelMakePointSource << <block, threads >> > (m_pDeviceData, uiSiteIndex, sourceData.m_bySpinIndex, sourceData.m_byColorIndex);
        }
    }
    break;
    case EFS_Wall:
    {
        preparethread;
        _kernelMakeWallSource << <block, threads >> > (
            m_pDeviceData, 
            sourceData.m_sSourcePoint.w, 
            sourceData.m_bySpinIndex, 
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

void CFieldFermionWilsonSquareSU3::SetKai(Real fKai)
{
    m_fKai = fKai;
    CCommonData::m_fKai = fKai;
}

BYTE* CFieldFermionWilsonSquareSU3::CopyDataOut(UINT &uiSize) const
{
    deviceWilsonVectorSU3* toSave = (deviceWilsonVectorSU3*)malloc(sizeof(deviceWilsonVectorSU3) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(Real) * m_uiSiteCount * 24);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        Real oneSite[24];
        for (UINT j = 0; j < 4; ++j)
        {
            for (UINT k = 0; k < 3; ++k)
            {
                oneSite[2 * (j * 3 + k)] = static_cast<Real>(toSave[i].m_d[j].m_ve[k].x);
                oneSite[2 * (j * 3 + k) + 1] = static_cast<Real>(toSave[i].m_d[j].m_ve[k].y);
            }
        }
        memcpy(saveData + sizeof(Real) * i * 24, oneSite, sizeof(Real) * 24);
    }
    
    //appGetFileSystem()->WriteAllBytes(fileName.c_str(), saveData, uiSize);
    //free(saveData);
    free(toSave);
    return saveData;
}

BYTE* CFieldFermionWilsonSquareSU3::CopyDataOutFloat(UINT& uiSize) const
{
    deviceWilsonVectorSU3* toSave = (deviceWilsonVectorSU3*)malloc(sizeof(deviceWilsonVectorSU3) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(FLOAT) * m_uiSiteCount * 24);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        FLOAT oneSite[24];
        for (UINT j = 0; j < 4; ++j)
        {
            for (UINT k = 0; k < 3; ++k)
            {
                oneSite[2 * (j * 3 + k)] = static_cast<FLOAT>(toSave[i].m_d[j].m_ve[k].x);
                oneSite[2 * (j * 3 + k) + 1] = static_cast<FLOAT>(toSave[i].m_d[j].m_ve[k].y);
            }
        }
        memcpy(saveData + sizeof(FLOAT) * i * 24, oneSite, sizeof(FLOAT) * 24);
    }

    //appGetFileSystem()->WriteAllBytes(fileName.c_str(), saveData, uiSize);
    //free(saveData);
    free(toSave);
    return saveData;
}

BYTE* CFieldFermionWilsonSquareSU3::CopyDataOutDouble(UINT& uiSize) const
{
    deviceWilsonVectorSU3* toSave = (deviceWilsonVectorSU3*)malloc(sizeof(deviceWilsonVectorSU3) * m_uiSiteCount);
    uiSize = static_cast<UINT>(sizeof(DOUBLE) * m_uiSiteCount * 24);
    BYTE* saveData = (BYTE*)malloc(static_cast<size_t>(uiSize));
    checkCudaErrors(cudaMemcpy(toSave, m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount, cudaMemcpyDeviceToHost));
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        DOUBLE oneSite[24];
        for (UINT j = 0; j < 4; ++j)
        {
            for (UINT k = 0; k < 3; ++k)
            {
                oneSite[2 * (j * 3 + k)] = static_cast<DOUBLE>(toSave[i].m_d[j].m_ve[k].x);
                oneSite[2 * (j * 3 + k) + 1] = static_cast<DOUBLE>(toSave[i].m_d[j].m_ve[k].y);
            }
        }
        memcpy(saveData + sizeof(DOUBLE) * i * 24, oneSite, sizeof(DOUBLE) * 24);
    }

    //appGetFileSystem()->WriteAllBytes(fileName.c_str(), saveData, uiSize);
    //free(saveData);
    free(toSave);
    return saveData;
}

TArray<CFieldFermion*> CFieldFermionWilsonSquareSU3::GetSourcesAtSiteFromPool(INT gaugeNum, INT bosonNum, const CFieldGauge* const* gaugeFields, const CFieldBoson* const* pBoson, const SSmallInt4& site) const
{
    TArray<CFieldFermion*> ret;
    for (UINT j = 0; j < 12; ++j)
    {
        ret.AddItem(dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId)));
        if (NULL == ret[j])
        {
            appCrucial(_T("GetSourcesAtSiteFromPool failed!\n"));
            _FAIL_EXIT;
        }
    }

    for (BYTE s = 0; s < 4; ++s)
    {
        for (BYTE c = 0; c < 3; ++c)
        {
            SFermionSource sourceData;
            sourceData.m_eSourceType = EFS_Point;
            sourceData.m_sSourcePoint = site;
            sourceData.m_byColorIndex = c;
            sourceData.m_bySpinIndex = s;

            ret[s * 3 + c]->InitialAsSource(sourceData);

            if (NULL != appGetFermionSolver(m_byFieldId) && !appGetFermionSolver(m_byFieldId)->IsAbsoluteAccuracy())
            {
                ret[s * 3 + c]->m_fLength = ret[s * 3 + c]->Dot(ret[s * 3 + c]).x;
            }
            ret[s * 3 + c]->InverseD(gaugeNum, bosonNum, gaugeFields, pBoson);
        }
    }
    return ret;
}

CCString CFieldFermionWilsonSquareSU3::GetInfos(const CCString &tab) const
{
    CCString sRet = CFieldFermion::GetInfos(tab);
    sRet = sRet + tab + _T("Hopping : ") + appToString(CCommonData::m_fKai) + _T("\n");
    return sRet;
}

#pragma region Even Odd Preconditioner

#if Discard_Even_odd_decomposition

/**
 * Calculate phi_even
 * Solve D z_even = phi_even
 * Calculate z_odd = D_oo^{-1} (phi_odd + D_oe z_even)
 */
UBOOL CFieldFermionWilsonSquareSU3::InverseD_eo(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    CFieldFermionWilsonSU3DEven* pEvenField = dynamic_cast<CFieldFermionWilsonSU3DEven*>(appGetLattice()->GetPooledFieldById(m_byEvenFieldId));
    if (NULL == pEvenField)
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 NO even field!"));
        return FALSE;
    }

    pEvenField->WriteEvenSites(this, pFieldSU3, FALSE);

    //Find a solver to solve me.
    if (!appGetFermionSolver(m_byEvenFieldId)->Solve(pEvenField, /*this is const*/pEvenField, pFieldSU3, EFO_F_D))
    {
        if (NULL != pEvenField)
        {
            pEvenField->Return();
        }
        return FALSE;
    }

    //z_even obtained, we calculate phi_odd + D_oe z_even
    pEvenField->WriteBackEvenSites(this, pFieldSU3, FALSE);
    pEvenField->Return();
    return TRUE;
}

UBOOL CFieldFermionWilsonSquareSU3::InverseDdagger_eo(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    const CFieldGaugeSU3* pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    CFieldFermionWilsonSU3DEven* pEvenField = dynamic_cast<CFieldFermionWilsonSU3DEven*>(appGetLattice()->GetPooledFieldById(m_byEvenFieldId));
    if (NULL == pEvenField)
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 NO even field!"));
        return FALSE;
    }

    pEvenField->WriteEvenSites(this, pFieldSU3, TRUE);

    //Find a solver to solve me.
    if (!appGetFermionSolver(m_byEvenFieldId)->Solve(pEvenField, /*this is const*/pEvenField, pFieldSU3, EFO_F_Ddagger))
    {
        if (NULL != pEvenField)
        {
            pEvenField->Return();
        }
        return FALSE;
    }

    //z_even obtained, we calculate phi_odd + D_oe z_even
    pEvenField->WriteBackEvenSites(this, pFieldSU3, TRUE);
    pEvenField->Return();

    return TRUE;
}

/**
 * Our solver can not solve this field (or, we need more than one solver, so we solve DDdagger twice
 * Maybe change to multi-solver
 * (D D+)^{-1} = D+^{-1} D^{-1}
 */
UBOOL CFieldFermionWilsonSquareSU3::InverseDDdagger_eo(const CField* pGauge)
{
    if (!InverseD_eo(pGauge))
    {
        return FALSE;
    }
    return InverseDdagger_eo(pGauge);
}

#endif

#pragma endregion

#pragma region Test Gamma5 Hermitian

UINT CFieldFermionWilsonSquareSU3::TestGamma5Hermitian(const CFieldGauge* pGauge, UBOOL bTestGamma5) const
{
    const UINT uiVolume = _HC_Volume;
    const UINT uiRealVolume = 12 * uiVolume;
    CLGComplex* matrixElement = (CLGComplex*)malloc(sizeof(CLGComplex) * uiRealVolume * uiRealVolume);
    CLGComplex* matrixElementd = NULL;
    if (bTestGamma5)
    {
        matrixElementd = (CLGComplex*)malloc(sizeof(CLGComplex) * uiRealVolume * uiRealVolume);
    }
    deviceWilsonVectorSU3* hostData = (deviceWilsonVectorSU3*)malloc(sizeof(deviceWilsonVectorSU3) * uiVolume);
    CFieldFermionWilsonSquareSU3* v = dynamic_cast<CFieldFermionWilsonSquareSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

    for (UINT i = 0; i < uiVolume; ++i)
    {
        const SSmallInt4 point = __hostSiteIndexToInt4(i);
        for (UINT j = 0; j < 12; ++j)
        {
            SFermionSource source;
            source.m_byColorIndex = static_cast<BYTE>(j % 3);
            source.m_bySpinIndex = static_cast<BYTE>(j / 3);
            source.m_eSourceType = EFS_Point;
            source.m_sSourcePoint = point;
            v->InitialAsSource(source);
            v->DS(pGauge);
            if (!bTestGamma5)
            {
                v->ApplyGamma(GAMMA5);
            }

            checkCudaErrors(cudaMemcpy(hostData, v->m_pDeviceData, sizeof(deviceWilsonVectorSU3) * uiVolume, cudaMemcpyDeviceToHost));

            const UINT x = i * 12 + j;
            for (UINT k = 0; k < uiVolume; ++k)
            {
                matrixElement[(12 * k +  0) * uiRealVolume + x] = hostData[k].m_d[0].m_ve[0];
                matrixElement[(12 * k +  1) * uiRealVolume + x] = hostData[k].m_d[0].m_ve[1];
                matrixElement[(12 * k +  2) * uiRealVolume + x] = hostData[k].m_d[0].m_ve[2];

                matrixElement[(12 * k +  3) * uiRealVolume + x] = hostData[k].m_d[1].m_ve[0];
                matrixElement[(12 * k +  4) * uiRealVolume + x] = hostData[k].m_d[1].m_ve[1];
                matrixElement[(12 * k +  5) * uiRealVolume + x] = hostData[k].m_d[1].m_ve[2];

                matrixElement[(12 * k +  6) * uiRealVolume + x] = hostData[k].m_d[2].m_ve[0];
                matrixElement[(12 * k +  7) * uiRealVolume + x] = hostData[k].m_d[2].m_ve[1];
                matrixElement[(12 * k +  8) * uiRealVolume + x] = hostData[k].m_d[2].m_ve[2];

                matrixElement[(12 * k +  9) * uiRealVolume + x] = hostData[k].m_d[3].m_ve[0];
                matrixElement[(12 * k + 10) * uiRealVolume + x] = hostData[k].m_d[3].m_ve[1];
                matrixElement[(12 * k + 11) * uiRealVolume + x] = hostData[k].m_d[3].m_ve[2];
            }

            if (bTestGamma5)
            {
                v->InitialAsSource(source);
                v->DdaggerS(pGauge);

                checkCudaErrors(cudaMemcpy(hostData, v->m_pDeviceData, sizeof(deviceWilsonVectorSU3) * uiVolume, cudaMemcpyDeviceToHost));

                for (UINT k = 0; k < uiVolume; ++k)
                {
                    matrixElementd[(12 * k +  0) * uiRealVolume + x] = hostData[k].m_d[0].m_ve[0];
                    matrixElementd[(12 * k +  1) * uiRealVolume + x] = hostData[k].m_d[0].m_ve[1];
                    matrixElementd[(12 * k +  2) * uiRealVolume + x] = hostData[k].m_d[0].m_ve[2];

                    matrixElementd[(12 * k +  3) * uiRealVolume + x] = hostData[k].m_d[1].m_ve[0];
                    matrixElementd[(12 * k +  4) * uiRealVolume + x] = hostData[k].m_d[1].m_ve[1];
                    matrixElementd[(12 * k +  5) * uiRealVolume + x] = hostData[k].m_d[1].m_ve[2];

                    matrixElementd[(12 * k +  6) * uiRealVolume + x] = hostData[k].m_d[2].m_ve[0];
                    matrixElementd[(12 * k +  7) * uiRealVolume + x] = hostData[k].m_d[2].m_ve[1];
                    matrixElementd[(12 * k +  8) * uiRealVolume + x] = hostData[k].m_d[2].m_ve[2];

                    matrixElementd[(12 * k +  9) * uiRealVolume + x] = hostData[k].m_d[3].m_ve[0];
                    matrixElementd[(12 * k + 10) * uiRealVolume + x] = hostData[k].m_d[3].m_ve[1];
                    matrixElementd[(12 * k + 11) * uiRealVolume + x] = hostData[k].m_d[3].m_ve[2];
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
        const SSmallInt4 xSite = __hostSiteIndexToInt4(x / 12);
        const SSmallInt4 ySite = __hostSiteIndexToInt4(y / 12);
        const UINT daggerIdx = y * uiRealVolume + x;
        const BYTE cx = (x % 12) % 3;
        const BYTE sx = (x % 12) / 3;
        const BYTE cy = (y % 12) % 3;
        const BYTE sy = (y % 12) / 3;

        if (_cuCabsf(matrixElement[i]) > F(0.0000001))
        {
            ++uiE;
            if (bTestGamma5)
            {
                if (appAbs(matrixElement[i].x - matrixElementd[daggerIdx].x) > F(0.0000001)
                 || appAbs(matrixElement[i].y + matrixElementd[daggerIdx].y) > F(0.0000001))
                {
                    ++uiWrong;
                    appGeneral(_T("[(%d, %d, %d, %d)_(%d %d)-(%d, %d, %d, %d)_(%d %d)]: D = %f + %f I   Ddagger = %f + %f I\n"),
                        xSite.x, xSite.y, xSite.z, xSite.w, cx, sx,
                        ySite.x, ySite.y, ySite.z, ySite.w, cy, sy,
                        matrixElement[i].x, matrixElement[i].y,
                        matrixElement[daggerIdx].x, matrixElement[daggerIdx].y);
                }
            }
            else
            {
                if (appAbs(matrixElement[i].x - matrixElement[daggerIdx].x) > F(0.0000001)
                    || appAbs(matrixElement[i].y + matrixElement[daggerIdx].y) > F(0.0000001))
                {
                    ++uiWrong;
                    appGeneral(_T("[(%d, %d, %d, %d)_(%d %d)-(%d, %d, %d, %d)_(%d %d)]: D = %f + %f I   Ddagger = %f + %f I\n"),
                        xSite.x, xSite.y, xSite.z, xSite.w, cx, sx,
                        ySite.x, ySite.y, ySite.z, ySite.w, cy, sy,
                        matrixElement[i].x, matrixElement[i].y,
                        matrixElement[daggerIdx].x, matrixElement[daggerIdx].y);
                }
            }
        }
    }
    v->Return();
    appSafeFree(matrixElement);
    if (bTestGamma5)
    {
        appSafeFree(matrixElementd);
    }
    appSafeFree(hostData);
    appGeneral(_T("%d none zero element checked, %d wrong found...\n"), uiE, uiWrong);
    return uiWrong;
}

#pragma endregion

#pragma region Field Matrix Operation

CFieldMatrixOperationWilsonSquareSU3::CFieldMatrixOperationWilsonSquareSU3()
{
    m_pHostResBuffer = (deviceWilsonVectorSU3**)malloc(sizeof(deviceWilsonVectorSU3*) * _kFieldMatrixMaxDim);
    m_pHostLeftBuffer = (deviceWilsonVectorSU3**)malloc(sizeof(deviceWilsonVectorSU3*) * _kFieldMatrixMaxDim);
    checkCudaErrors(cudaMalloc((void**)&m_pResBuffer, sizeof(deviceWilsonVectorSU3*) * _kFieldMatrixMaxDim));
    checkCudaErrors(cudaMalloc((void**)&m_pLeftBuffer, sizeof(deviceWilsonVectorSU3*) * _kFieldMatrixMaxDim));
}

CFieldMatrixOperationWilsonSquareSU3::~CFieldMatrixOperationWilsonSquareSU3()
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
__global__ void _CLG_LAUNCH_BOUND
_kernelMatrixMultiply(
    deviceWilsonVectorSU3 ** pRes,
    deviceWilsonVectorSU3 ** pLeft,
    const CLGComplex * __restrict__ pMatrix,
    UINT uiDimX, UINT uiDimY) //x=m,y=k
{
    intokernalE(12);

    CLGComplex result[CFieldMatrixOperation::_kFieldMatrixMaxDim];

    for (UINT i = 0; i < uiDimY; ++i)
    {
        result[i] = _make_cuComplex(F(0.0), F(0.0));
        for (UINT j = 0; j < uiDimX; ++j)
        {
            result[i] = _cuCaddf(result[i], _cuCmulf(
                j < uiDimY ? 
                  pRes[j][uiSiteIndex].m_me[elementIdx]
                : pLeft[j - uiDimY][uiSiteIndex].m_me[elementIdx],
                pMatrix[j * uiDimY + i]
            ));
        }
    }

    for (UINT i = 0; i < uiDimY; ++i)
    {
        pRes[i][uiSiteIndex].m_me[elementIdx] = result[i];
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

void CFieldMatrixOperationWilsonSquareSU3::VectorMultiplyMatrix(TArray<CField*>& res, const TArray<CField*>& left, const CLGComplex* deviceMatrix, UINT uiDimX, UINT uiDimY)
{
    for (UINT i = 0; i < uiDimY; ++i)
    {
        CFieldFermionWilsonSquareSU3* pF = dynamic_cast<CFieldFermionWilsonSquareSU3*>(res[i]);
        if (NULL == pF)
        {
            appCrucial(_T("CFieldMatrixOperationWilsonSquareSU3 only work with CFieldFermionWilsonSquareSU3!\n"));
            return;
        }
        m_pHostResBuffer[i] = pF->m_pDeviceData;
    }

    for (UINT i = 0; i < uiDimX - uiDimY; ++i)
    {
        const CFieldFermionWilsonSquareSU3* pF = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(left[i]);
        if (NULL == pF)
        {
            appCrucial(_T("CFieldMatrixOperationWilsonSquareSU3 only work with CFieldFermionWilsonSquareSU3!\n"));
            return;
        }
        m_pHostLeftBuffer[i] = pF->m_pDeviceData;
    }

    checkCudaErrors(cudaMemcpy(m_pResBuffer, m_pHostResBuffer, sizeof(deviceWilsonVectorSU3*) * uiDimY, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(m_pLeftBuffer, m_pHostLeftBuffer, sizeof(deviceWilsonVectorSU3*) * (uiDimX - uiDimY), cudaMemcpyHostToDevice));

    preparethreadE(12);
    _kernelMatrixMultiply << <block, threads >> > (m_pResBuffer, m_pLeftBuffer, deviceMatrix, uiDimX, uiDimY);
}

#pragma endregion

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================