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

#pragma region Kernel

__global__ void _kernelPrintFermionWilsonSquareSU3(const deviceWilsonVectorSU3 * __restrict__ pData)
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

__global__ void _kernelAxpyPlusFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 * pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex]);
}

__global__ void _kernelAxpyMinusFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 * pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther)
{
    intokernal;
    pMe[uiSiteIndex].Sub(pOther[uiSiteIndex]);
}

__global__ void _kernelAxpyComplexFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 * pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther, _Complex a)
{
    intokernal;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex].MulCompC(a));
}

__global__ void _kernelAxpyRealFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 * pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther, Real a)
{
    intokernal;
    pMe[uiSiteIndex].Add(pOther[uiSiteIndex].MulRealC(a));
}

__global__ void _kernelDotFermionWilsonSquareSU3(
    const deviceWilsonVectorSU3 * __restrict__ pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther, 
    _Complex * result)
{
    intokernal;
    result[uiSiteIndex] = pMe[uiSiteIndex].ConjugateDotC(pOther[uiSiteIndex]);
}

__global__ void _kernelScalarMultiplyComplex(
    deviceWilsonVectorSU3 * pMe, 
    _Complex a)
{
    intokernal;
    pMe[uiSiteIndex].MulComp(a);
}

__global__ void _kernelScalarMultiplyReal(
    deviceWilsonVectorSU3 * pMe, 
    Real a)
{
    intokernal;
    pMe[uiSiteIndex].MulReal(a);
}

/**
*
*/
__global__ void _kernelInitialFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 *pDevicePtr, 
    EFieldInitialType eInitialType)
{
    intokernal;

    switch (eInitialType)
    {
    case EFIT_Zero:
    {
        pDevicePtr[uiSiteIndex] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
    }
    break;
    case EFIT_RandomGaussian:
    {
        pDevicePtr[uiSiteIndex] = deviceWilsonVectorSU3::makeRandomGaussian(_deviceGetFatIndex(uiSiteIndex, 0));
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
* Dw phi(x) = phi(x) - kai sum _mu (1-gamma _mu) U(x,mu) phi(x+ mu) + (1+gamma _mu) U^{dagger}(x-mu) phi(x-mu)
* U act on su3
* gamma act on spinor
*
* If bDagger, it is gamma5, D, gamma5
*
*/
__global__ void _kernelDFermionWilsonSquareSU3(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    deviceWilsonVectorSU3* pResultData,
    Real kai,
    BYTE byFieldId,
    UBOOL bDiracChiralGamma,
    UBOOL bDDagger,
    EOperatorCoefficientType eCoeff,
    Real fCoeff,
    _Complex cCoeff)
{
    intokernaldir;

    gammaMatrix gamma5 = bDiracChiralGamma ? __diracGamma[GAMMA5] : __chiralGamma[GAMMA5];

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
        gammaMatrix gammaMu = bDiracChiralGamma ? __diracGamma[GAMMA1 + idir] : __chiralGamma[GAMMA1 + idir];

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        SIndex x_m_mu_Gauge = __idx->_deviceGaugeIndexWalk(uiSiteIndex, -(idir + 1));
        SIndex x_p_mu_Fermion = __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, (idir + 1));
        SIndex x_m_mu_Fermion = __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, -(idir + 1));

        //Assuming periodic
        //get U(x,mu), U^{dagger}(x-mu), 
        deviceSU3 x_Gauge_element = pGauge[linkIndex];
        deviceSU3 x_m_mu_Gauge_element = pGauge[_deviceGetLinkIndex(x_m_mu_Gauge.m_uiSiteIndex, idir)];
        x_m_mu_Gauge_element.Dagger();

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
        result.Add(u_phi_x_p_m);

        //- gammamu U(x,mu) phi(x+ mu)
        result.Sub(gammaMu.MulWilsonC(u_phi_x_p_m));

        //U^{dagger}(x-mu) phi(x-mu)
        deviceWilsonVectorSU3 u_dagger_phi_x_m_m = x_m_mu_Gauge_element.MulWilsonVector(x_m_mu_Fermion_element);
        result.Add(u_dagger_phi_x_m_m);

        //gammamu U^{dagger}(x-mu) phi(x-mu)
        result.Add(gammaMu.MulWilsonC(u_dagger_phi_x_m_m));
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
__global__ void _kernelDWilsonForceSU3(
    const deviceWilsonVectorSU3* __restrict__ pInverseD,
    const deviceWilsonVectorSU3* __restrict__ pInverseDDdagger,
    const deviceSU3* __restrict__ pGauge,
    deviceSU3* pForce,
    deviceSU3* pCachedForce,
    Real fKai,
    BYTE byFieldId,
    UBOOL bDiracChiralGamma)
{
    intokernaldir;

    deviceWilsonVectorSU3 x_Left(pInverseDDdagger[uiSiteIndex]);
    deviceWilsonVectorSU3 x_Right(pInverseD[uiSiteIndex]);

    //idir = mu
    for (UINT idir = 0; idir < uiDir; ++idir)
    {
        //Get Gamma mu
        gammaMatrix gammaMu = bDiracChiralGamma ? __diracGamma[GAMMA1 + idir] : __chiralGamma[GAMMA1 + idir];

        //x, mu
        UINT linkIndex = _deviceGetLinkIndex(uiSiteIndex, idir);

        SIndex x_m_mu_Gauge = __idx->_deviceGaugeIndexWalk(uiSiteIndex, -(idir + 1));
        SIndex x_p_mu_Fermion = __idx->_deviceFermionIndexWalk(byFieldId, uiSiteIndex, (idir + 1));

        deviceWilsonVectorSU3 x_p_mu_Right(pInverseD[x_p_mu_Fermion.m_uiSiteIndex]);
        deviceWilsonVectorSU3 x_p_mu_Left(pInverseDDdagger[x_p_mu_Fermion.m_uiSiteIndex]);
        deviceSU3 x_Gauge_element = pGauge[linkIndex];

        deviceWilsonVectorSU3 right1(x_p_mu_Right);
        right1.Sub(gammaMu.MulWilsonC(right1));
        deviceSU3 mid = deviceSU3::makeSU3Contract(x_Left, right1);

        deviceWilsonVectorSU3 right2(x_Right);
        right2.Add(gammaMu.MulWilsonC(right2));
        mid.Add(deviceSU3::makeSU3Contract(right2, x_p_mu_Left));

        deviceSU3 forceOfThisLink = x_Gauge_element.MulC(mid);
        forceOfThisLink.Ta();
        forceOfThisLink.MulReal(fKai);

        pForce[linkIndex].Add(forceOfThisLink);
        if (NULL != pCachedForce)
        {
            pCachedForce[linkIndex] = forceOfThisLink;
        }
    }
}

__global__ void _kernelApplyGammaSU3(deviceWilsonVectorSU3* pDeviceData, UINT uiGamma, UBOOL bDiracChiralGamma)
{
    intokernal;
    pDeviceData[uiSiteIndex] = (bDiracChiralGamma ? __diracGamma[uiGamma] : __chiralGamma[uiGamma]).MulWilsonC(pDeviceData[uiSiteIndex]);
}

#pragma endregion

CFieldFermionWilsonSquareSU3::CFieldFermionWilsonSquareSU3() : CFieldFermion(), m_fKai(F(0.125))
{
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDataCopy, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount));
}

CFieldFermionWilsonSquareSU3::~CFieldFermionWilsonSquareSU3()
{
    checkCudaErrors(cudaFree(m_pDeviceData));
    checkCudaErrors(cudaFree(m_pDeviceDataCopy));
}

/**
*
*/
void CFieldFermionWilsonSquareSU3::InitialField(EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceData, eInitialType);
}

void CFieldFermionWilsonSquareSU3::InitialFieldWithFile(const CCString& sFileName, EFieldFileType eFieldType)
{
    if (eFieldType != EFFT_CLGBin)
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3::InitialFieldWithFile: Only support CLG Bin File\n"));
        return;
    }

    UINT uiSize = static_cast<UINT>(sizeof(FLOAT) * 24 * m_uiSiteCount);
    deviceWilsonVectorSU3* readData = (deviceWilsonVectorSU3*)malloc(sizeof(deviceWilsonVectorSU3) * m_uiSiteCount);
    BYTE* data = appGetFileSystem()->ReadAllBytes(sFileName.c_str(), uiSize);
    for (UINT i = 0; i < m_uiSiteCount; ++i)
    {
        FLOAT thisSite[24];
        memcpy(thisSite, data + i * sizeof(FLOAT) * 24, sizeof(FLOAT) * 24);
        for (UINT j = 0; j < 4; ++j)
        {
            for (UINT k = 0; k < 3; ++k)
            {
                readData[i].m_d[j].m_ve[k] = _make_cuComplex(
                    static_cast<Real>(thisSite[2 * (j * 3 + k)]), 
                    static_cast<Real>(thisSite[2 * (j * 3 + k) + 1]));
            }
        }
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceData, readData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount, cudaMemcpyHostToDevice));
    free(data);
    free(readData);
}

void CFieldFermionWilsonSquareSU3::InitialOtherParameters(CParameters& params)
{
    params.FetchValueReal(_T("Hopping"), m_fKai);
    if (m_fKai < F(0.00000001))
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3: Kai is nearly 0, such that Dphi \approx phi! This will cause problem!\n"));
    }
    CCommonData::m_fKai = m_fKai;
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

void CFieldFermionWilsonSquareSU3::Axpy(const _Complex& a, const CField* x)
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

_Complex CFieldFermionWilsonSquareSU3::Dot(const CField* x) const
{
    if (NULL == x || EFT_FermionWilsonSquareSU3 != x->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only copy to CFieldFermionWilsonSquareSU3!"));
        return _make_cuComplex(0,0);
    }
    const CFieldFermionWilsonSquareSU3 * pField = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(x);

    preparethread;
    _kernelDotFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceData, pField->m_pDeviceData, _D_ComplexThreadBuffer);

    return appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
}

void CFieldFermionWilsonSquareSU3::ScalarMultply(const _Complex& a)
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
    _kernelApplyGammaSU3 << <block, threads >> >(m_pDeviceData, static_cast<UINT>(eGamma), TRUE);
}
/**
* generate phi by gaussian random.
* phi = D phi
*/
void CFieldFermionWilsonSquareSU3::PrepareForHMC(const CFieldGauge* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    preparethread;
    _kernelInitialFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceDataCopy, EFIT_RandomGaussian);
    _kernelDFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceDataCopy, pFieldSU3->m_pDeviceData, m_pDeviceData, m_fKai, m_byFieldId, TRUE, FALSE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));

    //cache a inverse DDdagger field
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
    CFieldFermionWilsonSquareSU3* pCachedSU3 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pField);
    pCachedSU3->InverseDDdagger(pGauge);
}

//Kai should be part of D operator
void CFieldFermionWilsonSquareSU3::D(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    checkCudaErrors(cudaMemcpy(m_pDeviceDataCopy, m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    _Complex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    preparethread;
    _kernelDFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceDataCopy, pFieldSU3->m_pDeviceData, m_pDeviceData, m_fKai, m_byFieldId, TRUE, FALSE, eCoeffType, fRealCoeff, cCompCoeff);
}

//Kai should be part of D operator
void CFieldFermionWilsonSquareSU3::Ddagger(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    checkCudaErrors(cudaMemcpy(m_pDeviceDataCopy, m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    Real fRealCoeff = fCoeffReal;
    _Complex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    preparethread;
    _kernelDFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceDataCopy, pFieldSU3->m_pDeviceData, m_pDeviceData, m_fKai, m_byFieldId, TRUE, TRUE, eCoeffType, fRealCoeff, cCompCoeff);
}

void CFieldFermionWilsonSquareSU3::DDdagger(const CField* pGauge, EOperatorCoefficientType eCoeffType, Real fCoeffReal, Real fCoeffImg)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    Real fRealCoeff = fCoeffReal;
    _Complex cCompCoeff = _make_cuComplex(fCoeffReal, fCoeffImg);
    if (EOCT_Minus == eCoeffType)
    {
        eCoeffType = EOCT_Real;
        fRealCoeff = F(-1.0);
    }

    preparethread;
    //Ddagger first, m_pDeviceDataCopy = D+ m_pDeviceData
    _kernelDFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceData, pFieldSU3->m_pDeviceData, m_pDeviceDataCopy, m_fKai, m_byFieldId, TRUE, TRUE, EOCT_None, F(1.0), _make_cuComplex(F(1.0), F(0.0)));
    //Then D, m_pDeviceData = D m_pDeviceDataCopy
    _kernelDFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceDataCopy, pFieldSU3->m_pDeviceData, m_pDeviceData, m_fKai, m_byFieldId, TRUE, FALSE, eCoeffType, fRealCoeff, cCompCoeff);
}

UBOOL CFieldFermionWilsonSquareSU3::InverseD(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    //Find a solver to solve me.
    return appGetFermionSolver()->Solve(this, /*this is const*/this, pFieldSU3, EFO_F_D);
}

UBOOL CFieldFermionWilsonSquareSU3::InverseDdagger(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    //Find a solver to solve me.
    return appGetFermionSolver()->Solve(this, /*this is const*/this, pFieldSU3, EFO_F_Ddagger);
}

UBOOL CFieldFermionWilsonSquareSU3::InverseDDdagger(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    //Find a solver to solve me.
    return appGetFermionSolver()->Solve(this, /*this is const*/this, pFieldSU3, EFO_F_DDdagger);
}

UBOOL CFieldFermionWilsonSquareSU3::CalculateForce(const CFieldGauge* pGauge, CFieldGauge* pForce, CFieldGauge* pCachedForce) const
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
    if (NULL != pCachedForce && EFT_GaugeSU3 != pCachedForce->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return FALSE;
    }

    const CFieldGaugeSU3 * pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);
    CFieldGaugeSU3 * pForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pForce);
    CFieldGaugeSU3 * pCacheForceForceSU3 = NULL;
    if (NULL != pCachedForce)
    {
        pCacheForceForceSU3 = dynamic_cast<CFieldGaugeSU3*>(pCachedForce);
    }

    CField * pDDaggerPhi = appGetLattice()->GetPooledFieldById(m_byFieldId);
    CField * pDPhi = appGetLattice()->GetPooledFieldById(m_byFieldId);
    CField * pCachedField = appGetLattice()->m_pFieldCache->GetCachedField(CFieldCache::CachedInverseDDdaggerField);
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
    if (!appGetFermionSolver()->Solve(pDDaggerPhiWilson, this, pGaugeSU3, EFO_F_DDdagger, pCachedField))
    {
        appCrucial(_T("Sparse Linear Solver failed...\n"));
        pDDaggerPhi->Return();
        pDPhi->Return();
        return FALSE;
    }
    //phi 2 = D^{-1}phi = D+ (DD+)^{-1} phi
    //It is faster to calcuate D+ phi2 then D^{-1} phi
    pDDaggerPhiWilson->CopyTo(pDPhiWilson);
    if (NULL != pCachedField)
    {
        //The gauge field is changing slowly, and D depends only on gauge, also change slowly
        //Use the last solution as start point will accelerate the solver, so we cache it
        pDDaggerPhiWilson->CopyTo(pCachedField);
    }
    pDPhiWilson->Ddagger(pGaugeSU3);

    preparethread;
    _kernelDWilsonForceSU3 << <block, threads >> > (
        pDPhiWilson->m_pDeviceData,
        pDDaggerPhiWilson->m_pDeviceData,
        pGaugeSU3->m_pDeviceData,
        pForceSU3->m_pDeviceData,
        NULL == pCacheForceForceSU3 ? NULL : pCacheForceForceSU3->m_pDeviceData,
        m_fKai, m_byFieldId, TRUE);

    pDDaggerPhi->Return();
    pDPhi->Return();

    return TRUE;
}

void CFieldFermionWilsonSquareSU3::SetKai(Real fKai)
{
    m_fKai = fKai;
    CCommonData::m_fKai = fKai;
}

void CFieldFermionWilsonSquareSU3::SaveToFile(const CCString &fileName) const
{
    deviceWilsonVectorSU3* toSave = (deviceWilsonVectorSU3*)malloc(sizeof(deviceWilsonVectorSU3) * m_uiSiteCount);
    UINT uiSize = static_cast<UINT>(sizeof(FLOAT) * m_uiSiteCount * 24);
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
    appGetFileSystem()->WriteAllBytes(fileName.c_str(), saveData, uiSize);
    free(saveData);
    free(toSave);
}

CCString CFieldFermionWilsonSquareSU3::GetInfos(const CCString &tab) const
{
    CCString sRet;
    sRet = tab + _T("Name : CFieldFermionWilsonSquareSU3\n");
    sRet = sRet + tab + _T("Hopping : ") + appFloatToString(CCommonData::m_fKai) + _T("\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================