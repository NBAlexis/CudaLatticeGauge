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

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        UINT siteIndexX = _deviceGetSiteIndex(coord);

        printf("%d,%d,%d,%d=((%1.2f %1.2fi, %1.2f %1.2fi, %1.2f %1.2fi),(%1.2f %1.2fi, %1.2f %1.2fi, %1.2f %1.2fi),(%1.2f %1.2fi, %1.2f %1.2fi, %1.2f %1.2fi),(%1.2f %1.2fi, %1.2f %1.2fi, %1.2f %1.2fi))\n", 
            coord[0], coord[1], coord[2], coord[3],
            pData[siteIndexX].m_d[0].m_ve[0].x, pData[siteIndexX].m_d[0].m_ve[0].y,
            pData[siteIndexX].m_d[0].m_ve[1].x, pData[siteIndexX].m_d[0].m_ve[1].y,
            pData[siteIndexX].m_d[0].m_ve[2].x, pData[siteIndexX].m_d[0].m_ve[2].y,

            pData[siteIndexX].m_d[1].m_ve[0].x, pData[siteIndexX].m_d[1].m_ve[0].y,
            pData[siteIndexX].m_d[1].m_ve[1].x, pData[siteIndexX].m_d[1].m_ve[1].y,
            pData[siteIndexX].m_d[1].m_ve[2].x, pData[siteIndexX].m_d[1].m_ve[2].y,

            pData[siteIndexX].m_d[2].m_ve[0].x, pData[siteIndexX].m_d[2].m_ve[0].y,
            pData[siteIndexX].m_d[2].m_ve[1].x, pData[siteIndexX].m_d[2].m_ve[1].y,
            pData[siteIndexX].m_d[2].m_ve[2].x, pData[siteIndexX].m_d[2].m_ve[2].y,

            pData[siteIndexX].m_d[3].m_ve[0].x, pData[siteIndexX].m_d[3].m_ve[0].y,
            pData[siteIndexX].m_d[3].m_ve[1].x, pData[siteIndexX].m_d[3].m_ve[1].y,
            pData[siteIndexX].m_d[3].m_ve[2].x, pData[siteIndexX].m_d[3].m_ve[2].y
            );
    }
}

__global__ void _kernelAxpyPlusFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 * pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther)
{
    intokernal;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        UINT siteIndex = _deviceGetSiteIndex(coord);
        pMe[siteIndex].Add(pOther[siteIndex]);
    }
}

__global__ void _kernelAxpyMinusFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 * pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther)
{
    intokernal;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        UINT siteIndex = _deviceGetSiteIndex(coord);
        pMe[siteIndex].Sub(pOther[siteIndex]);
    }
}

__global__ void _kernelAxpyComplexFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 * pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther, _Complex a)
{
    intokernal;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        UINT siteIndex = _deviceGetSiteIndex(coord);
        pMe[siteIndex].Add(pOther[siteIndex].MulCompC(a));
    }
}

__global__ void _kernelAxpyRealFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 * pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther, Real a)
{
    intokernal;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        UINT siteIndex = _deviceGetSiteIndex(coord);
        pMe[siteIndex].Add(pOther[siteIndex].MulRealC(a));
    }
}

__global__ void _kernelDotFermionWilsonSquareSU3(
    const deviceWilsonVectorSU3 * __restrict__ pMe, 
    const deviceWilsonVectorSU3 * __restrict__ pOther, 
    _Complex * result)
{
    intokernal;
    _Complex res = _make_cuComplex(F(0.0), F(0.0));
    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        UINT siteIndex = _deviceGetSiteIndex(coord);
        res = _cuCaddf(res, pMe[siteIndex].ConjugateDotC(pOther[siteIndex]));
    }
    result[__thread_id] = res;
}

__global__ void _kernelScalarMultiplyComplex(
    deviceWilsonVectorSU3 * pMe, 
    _Complex a)
{
    intokernal;
    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        UINT siteIndex = _deviceGetSiteIndex(coord);
        pMe[siteIndex].MulComp(a);
    }
}

__global__ void _kernelScalarMultiplyReal(
    deviceWilsonVectorSU3 * pMe, 
    Real a)
{
    intokernal;
    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        UINT siteIndex = _deviceGetSiteIndex(coord);
        pMe[siteIndex].MulReal(a);
    }
}

/**
* phi dagger, phi
*/
__global__ void _kernel_This_IsNot_Dot_FermionWilsonSquareSU3(
    const deviceWilsonVectorSU3 * __restrict__ pLeft,
    const deviceWilsonVectorSU3 * __restrict__ pRight,
    deviceSU3* result)
{
    intokernaldir;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;

        for (int idir = 0; idir < uiDir; ++idir)
        {
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);

            deviceSU3 resultThisLink = deviceSU3::makeSU3Zero();
            for (int i = 0; i < 8; ++i)
            {
                _Complex omega = pLeft[linkIndex * 8 + i].ConjugateDotC(pRight[linkIndex * 8 + i]);
                resultThisLink.Add(__SU3Generators[i]->MulCompC(omega));
            }
            result[linkIndex] = resultThisLink;
        }
    }
}

/**
*
*/
__global__ void _kernelInitialFermionWilsonSquareSU3(
    deviceWilsonVectorSU3 *pDevicePtr, 
    EFieldInitialType eInitialType)
{
    intokernal;

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        UINT siteIndexX = _deviceGetSiteIndex(coord);
        UINT fatIndex = _deviceGetFatIndex(siteIndexX, 0);

        switch (eInitialType)
        {
            case EFIT_Zero:
                {
                    pDevicePtr[siteIndexX] = deviceWilsonVectorSU3::makeZeroWilsonVectorSU3();
                }
                break;
            case EFIT_RandomGaussian:
                {
                    pDevicePtr[siteIndexX] = deviceWilsonVectorSU3::makeRandomGaussian(fatIndex);
                }
                break;
            default:
                {
                    printf("Wilson Fermion Field cannot be initialized with this type!");
                }
            break;
        }
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
    UBOOL bDDagger)
{
    intokernaldir;

    gammaMatrix gamma5 = bDiracChiralGamma ? __diracGamma->m_gm[gammaMatrixSet::GAMMA5] : __chiralGamma->m_gm[gammaMatrixSet::GAMMA5];

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        //x
        UINT siteIndexX = _deviceGetSiteIndex(coord);
        deviceWilsonVectorSU3 result;
        pResultData[siteIndexX] = pDeviceData[siteIndexX];
        if (bDDagger)
        {
            pResultData[siteIndexX] = gamma5.MulWilsonC(pResultData[siteIndexX]);
        }

        //idir = mu
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            //Get Gamma mu
            gammaMatrix gammaMu = bDiracChiralGamma ? 
                  __diracGamma->m_gm[gammaMatrixSet::GAMMA1 + idir]
                : __chiralGamma->m_gm[gammaMatrixSet::GAMMA1 + idir];

            //x, mu
            UINT linkIndex = _deviceGetLinkIndex(siteIndexX, idir);

            SIndex x_m_mu_Gauge = __idx->_deviceGaugeIndexWalk(siteIndexX, -(idir + 1));
            SIndex x_p_mu_Fermion = __idx->_deviceFermionIndexWalk(byFieldId, siteIndexX, (idir + 1));
            SIndex x_m_mu_Fermion = __idx->_deviceFermionIndexWalk(byFieldId, siteIndexX, -(idir + 1));          

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
            for (UINT iSpinor = 0; iSpinor < 4; ++iSpinor) //Wilson fermion is 4-spinor
            {
                //U(x,mu) phi(x+ mu)
                result.m_d[iSpinor].Add(x_Gauge_element.MulVector(x_p_mu_Fermion_element.m_d[iSpinor]));

                //- gammamu U(x,mu) phi(x+ mu)
                result.m_d[iSpinor].Sub(x_Gauge_element.MulVector(gammaMu.MulVectorC(x_p_mu_Fermion_element, iSpinor)));

                //U^{dagger}(x-mu) phi(x-mu)
                result.m_d[iSpinor].Add(x_m_mu_Gauge_element.MulVector(x_m_mu_Fermion_element.m_d[iSpinor]));

                //gammamu U^{dagger}(x-mu) phi(x-mu)
                result.m_d[iSpinor].Add(x_m_mu_Gauge_element.MulVector(gammaMu.MulVectorC(x_m_mu_Fermion_element, iSpinor)));
            }
        }

        //result = phi(x) - kai sum _mu result
        result.MulReal(kai);
        pResultData[siteIndexX].Sub(result);
        if (bDDagger)
        {
            pResultData[siteIndexX] = gamma5.MulWilsonC(pResultData[siteIndexX]);
        }
    }
}

/**
* The output is on a gauge field
* Therefor cannot make together with _kernelDWilson
*
*/
__global__ void _kernelDWilsonMuSU3(
    const deviceWilsonVectorSU3* __restrict__ pDeviceData,
    const deviceSU3* __restrict__ pGauge,
    deviceWilsonVectorSU3* pResultDataArray,
    Real kai,
    BYTE byFieldId,
    UBOOL bDiracChiralGamma,
    UBOOL bDDagger,
    UBOOL bPartialOmega)
{
    intokernaldir;

    gammaMatrix gamma5 = bDiracChiralGamma ? __diracGamma->m_gm[gammaMatrixSet::GAMMA5] : __chiralGamma->m_gm[gammaMatrixSet::GAMMA5];

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        //x
        UINT siteIndexX = _deviceGetSiteIndex(coord);
        deviceWilsonVectorSU3 x_Fermion_element = pDeviceData[siteIndexX];
        if (bDDagger)
        {
            x_Fermion_element = gamma5.MulWilsonC(x_Fermion_element);
        }

        //idir = mu
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
            deviceWilsonVectorSU3 result[8];

            //Get Gamma mu
            gammaMatrix gammaMu = bDiracChiralGamma ?
                  __diracGamma->m_gm[gammaMatrixSet::GAMMA1 + idir]
                : __chiralGamma->m_gm[gammaMatrixSet::GAMMA1 + idir];

            //x, mu
            UINT linkIndex = _deviceGetLinkIndex(coord, idir);

            SIndex x_m_mu_Gauge = __idx->_deviceGaugeIndexWalk(siteIndexX, -(idir + 1));
            SIndex x_p_mu_Fermion = __idx->_deviceFermionIndexWalk(byFieldId, siteIndexX, (idir + 1));
            SIndex x_m_mu_Fermion = __idx->_deviceFermionIndexWalk(byFieldId, siteIndexX, -(idir + 1));

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
            for (UINT iSpinor = 0; iSpinor < 4; ++iSpinor) //Wilson fermion is 4-spinor
            {
                for (int i = 0; i < 8; ++i)
                {
                    if (!bPartialOmega)
                    {
                        //U(x,mu) phi(x+ mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].AddC(x_Gauge_element.MulVector(x_p_mu_Fermion_element.m_d[iSpinor]));

                        //- gammamu U(x,mu) phi(x+ mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].SubC(x_Gauge_element.MulVector(gammaMu.MulVectorC(x_p_mu_Fermion_element, iSpinor)));

                        //U^{dagger}(x-mu) phi(x-mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].AddC(x_m_mu_Gauge_element.MulVector(x_m_mu_Fermion_element.m_d[iSpinor]));

                        //gammamu U^{dagger}(x-mu) phi(x-mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].AddC(x_m_mu_Gauge_element.MulVector(gammaMu.MulVectorC(x_m_mu_Fermion_element, iSpinor)));
                    }
                    else
                    {
                        //U(x,mu) phi(x+ mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].AddC(__SU3Generators[i]->MulC(x_Gauge_element).MulVector(x_p_mu_Fermion_element.m_d[iSpinor]));

                        //- gammamu U(x,mu) phi(x+ mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].SubC(__SU3Generators[i]->MulC(x_Gauge_element).MulVector(gammaMu.MulVectorC(x_p_mu_Fermion_element, iSpinor)));

                        //U^{dagger}(x-mu) phi(x-mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].AddC(__SU3Generators[i]->MulC(x_m_mu_Gauge_element).MulVector(x_m_mu_Fermion_element.m_d[iSpinor]));

                        //gammamu U^{dagger}(x-mu) phi(x-mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].AddC(__SU3Generators[i]->MulC(x_m_mu_Gauge_element).MulVector(gammaMu.MulVectorC(x_m_mu_Fermion_element, iSpinor)));
                    }
                }

            }

            for (int i = 0; i < 8; ++i)
            {
                if (!bPartialOmega)
                {
                    //result = phi(x) - kai sum _mu result
                    result[i].MulReal(kai);
                    pResultDataArray[linkIndex * 8 + i] = x_Fermion_element.SubC(result[i]);
                    if (bDDagger)
                    {
                        pResultDataArray[linkIndex * 8 + i] = gamma5.MulWilsonC(pResultDataArray[linkIndex * 8 + i]);
                    }
                }
                else
                {
                    result[i].MulComp(_make_cuComplex(0, -kai));
                    pResultDataArray[linkIndex * 8 + i] = result[i];
                    if (bDDagger)
                    {
                        pResultDataArray[linkIndex * 8 + i] = gamma5.MulWilsonC(pResultDataArray[linkIndex * 8 + i]);
                    }
                }
            }
        }
    }
}

#pragma endregion

CFieldFermionWilsonSquareSU3::CFieldFermionWilsonSquareSU3() : CFieldFermion(), m_fKai(F(1.0))
{
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDataCopy, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount));

    //checkCudaErrors(cudaMalloc((void**)&m_pForceRightVector, sizeof(deviceWilsonVectorSU3) * m_uiLinkeCount * 8));
    //checkCudaErrors(cudaMalloc((void**)&m_pForceRightVectorCopy, sizeof(deviceWilsonVectorSU3) * m_uiLinkeCount * 8));
    //checkCudaErrors(cudaMalloc((void**)&m_pForceLeftVector, sizeof(deviceWilsonVectorSU3) * m_uiLinkeCount * 8));
    //checkCudaErrors(cudaMalloc((void**)&m_pForceLeftVectorCopy, sizeof(deviceWilsonVectorSU3) * m_uiLinkeCount * 8));
}

CFieldFermionWilsonSquareSU3::~CFieldFermionWilsonSquareSU3()
{
    checkCudaErrors(cudaFree(m_pDeviceData));
    checkCudaErrors(cudaFree(m_pDeviceDataCopy));

    //checkCudaErrors(cudaFree(m_pForceRightVector));
    //checkCudaErrors(cudaFree(m_pForceRightVectorCopy));
    //checkCudaErrors(cudaFree(m_pForceLeftVector));
    //checkCudaErrors(cudaFree(m_pForceLeftVectorCopy));
}

/**
*
*/
void CFieldFermionWilsonSquareSU3::InitialField(EFieldInitialType eInitialType)
{
    preparethread;
    _kernelInitialFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceData, eInitialType);
}

void CFieldFermionWilsonSquareSU3::InitialOtherParameters(CParameters& params)
{
    params.FetchValueReal(_T("Hopping"), m_fKai);
    if (appAbs(m_fKai) < F(0.00000001))
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3: Kai is nearly 0, such that Dphi \approx phi! This will cause problem!\n"));
    }
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
    CFieldFermionWilsonSquareSU3 * pField = dynamic_cast<CFieldFermionWilsonSquareSU3*>(U);
    checkCudaErrors(cudaMemcpy(pField->m_pDeviceData, m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount, cudaMemcpyDeviceToDevice));
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
    _kernelDFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceDataCopy, pFieldSU3->m_pDeviceData, m_pDeviceData, m_fKai, m_byFieldId, TRUE, FALSE);

}

//Kai should be part of D operator
void CFieldFermionWilsonSquareSU3::D(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    checkCudaErrors(cudaMemcpy(m_pDeviceDataCopy, m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    preparethread;
    _kernelDFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceDataCopy, pFieldSU3->m_pDeviceData, m_pDeviceData, m_fKai, m_byFieldId, TRUE, FALSE);
}

//Kai should be part of D operator
void CFieldFermionWilsonSquareSU3::Ddagger(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    checkCudaErrors(cudaMemcpy(m_pDeviceDataCopy, m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount, cudaMemcpyDeviceToDevice));

    preparethread;
    _kernelDFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceDataCopy, pFieldSU3->m_pDeviceData, m_pDeviceData, m_fKai, m_byFieldId, TRUE, TRUE);
}

void CFieldFermionWilsonSquareSU3::DDdagger(const CField* pGauge)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CFieldFermionWilsonSquareSU3 can only play with gauge SU3!"));
        return;
    }
    const CFieldGaugeSU3 * pFieldSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    preparethread;
    _kernelDFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceData, pFieldSU3->m_pDeviceData, m_pDeviceDataCopy, m_fKai, m_byFieldId, TRUE, TRUE);
    _kernelDFermionWilsonSquareSU3 << <block, threads >> > (m_pDeviceDataCopy, pFieldSU3->m_pDeviceData, m_pDeviceData, m_fKai, m_byFieldId, TRUE, FALSE);
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

void CFieldFermionWilsonSquareSU3::CalculateForce(const CFieldGauge* pGauge, CFieldGauge* pForce)
{

}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================