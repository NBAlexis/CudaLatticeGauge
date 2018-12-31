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

        printf("xyzt:%d,%d,%d,%d = ((%f+%f i, %f+%f i, %f+%f i),(%f+%f i, %f+%f i, %f+%f i),(%f+%f i, %f+%f i, %f+%f i),(%f+%f i, %f+%f i, %f+%f i))\n", 
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

/**
* phi dagger, phi
*/
__global__ void _kernelDotFermionWilsonSquareSU3(const deviceWilsonVectorSU3 * __restrict__ pLeft,
                                           const deviceWilsonVectorSU3 * __restrict__ pRight,
                                           deviceSU3* result)
{
    intokernal;

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
                resultThisLink.Add(__SU3Generators[i]->Mulc(omega));
            }
            result[linkIndex] = resultThisLink;
        }
    }
}

/**
*
*/
__global__ void _kernelInitialFermionWilsonSquareSU3(deviceWilsonVectorSU3 *pDevicePtr, EFieldInitialType eInitialType)
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
                    pDevicePtr[siteIndexX].MakeZero();
                }
                break;
            case EFIT_RandomGaussian:
                {
                    pDevicePtr[siteIndexX].MakeRandomGaussian(fatIndex);
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
__global__ void _kernelDFermionWilsonSquareSU3(const deviceWilsonVectorSU3* __restrict__ pDeviceData,
                                  const deviceSU3* __restrict__ pGauge,
                                  deviceWilsonVectorSU3* pResultData,
                                  Real kai,
                                  BYTE byFieldId,
                                  UBOOL bDiracChiralGamma,
                                  UBOOL bDDagger)
{
    intokernal;

    gammaMatrix gamma5 = bDiracChiralGamma ? __diracGamma->m_gm[gammaMatrixSet::GAMMA5] : __chiralGamma->m_gm[gammaMatrixSet::GAMMA5];

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        //x
        UINT siteIndexX = _deviceGetSiteIndex(coord);
        deviceWilsonVectorSU3 result;
        deviceWilsonVectorSU3 x_Fermion_element = pDeviceData[siteIndexX];
        if (bDDagger)
        {
            x_Fermion_element = gamma5.MulC(x_Fermion_element);
        }

        //idir = mu
        for (UINT idir = 0; idir < uiDir; ++idir)
        {
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
                x_p_mu_Fermion_element = gamma5.MulC(x_p_mu_Fermion_element);
                x_m_mu_Fermion_element = gamma5.MulC(x_m_mu_Fermion_element);
            }

            //hopping terms
            for (UINT iSpinor = 0; iSpinor < 4; ++iSpinor) //Wilson fermion is 4-spinor
            {
                //U(x,mu) phi(x+ mu)
                result.m_d[iSpinor] = result.m_d[iSpinor].AddC(x_Gauge_element.Mul(x_p_mu_Fermion_element.m_d[iSpinor]));

                //- gammamu U(x,mu) phi(x+ mu)
                result.m_d[iSpinor] = result.m_d[iSpinor].SubC(x_Gauge_element.Mul(gammaMu.MulC(x_p_mu_Fermion_element, iSpinor)));

                //U^{dagger}(x-mu) phi(x-mu)
                result.m_d[iSpinor] = result.m_d[iSpinor].AddC(x_m_mu_Gauge_element.Mul(x_m_mu_Fermion_element.m_d[iSpinor]));

                //gammamu U^{dagger}(x-mu) phi(x-mu)
                result.m_d[iSpinor] = result.m_d[iSpinor].AddC(x_m_mu_Gauge_element.Mul(gammaMu.MulC(x_m_mu_Fermion_element, iSpinor)));
            }
        }

        //result = phi(x) - kai sum _mu result
        result.Mul(_make_cuComplex(kai, 0));
        pResultData[siteIndexX] = x_Fermion_element.SubC(result);
        if (bDDagger)
        {
            pResultData[siteIndexX] = gamma5.MulC(pResultData[siteIndexX]);
        }
    }
}

/**
* The output is on a gauge field
* Therefor cannot make together with _kernelDWilson
*
*/
__global__ void _kernelDWilsonMuSU3(const deviceWilsonVectorSU3* __restrict__ pDeviceData,
                                    const deviceSU3* __restrict__ pGauge,
                                    deviceWilsonVectorSU3* pResultDataArray,
                                    Real kai,
                                    BYTE byFieldId,
                                    UBOOL bDiracChiralGamma,
                                    UBOOL bDDagger,
                                    UBOOL bPartialOmega)
{
    intokernal;

    gammaMatrix gamma5 = bDiracChiralGamma ? __diracGamma->m_gm[gammaMatrixSet::GAMMA5] : __chiralGamma->m_gm[gammaMatrixSet::GAMMA5];

    for (UINT it = 0; it < uiTLength; ++it)
    {
        coord[3] = it;
        //x
        UINT siteIndexX = _deviceGetSiteIndex(coord);
        deviceWilsonVectorSU3 x_Fermion_element = pDeviceData[siteIndexX];
        if (bDDagger)
        {
            x_Fermion_element = gamma5.MulC(x_Fermion_element);
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
                x_p_mu_Fermion_element = gamma5.MulC(x_p_mu_Fermion_element);
                x_m_mu_Fermion_element = gamma5.MulC(x_m_mu_Fermion_element);
            }

            //hopping terms
            for (UINT iSpinor = 0; iSpinor < 4; ++iSpinor) //Wilson fermion is 4-spinor
            {
                for (int i = 0; i < 8; ++i)
                {
                    if (!bPartialOmega)
                    {
                        //U(x,mu) phi(x+ mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].AddC(x_Gauge_element.Mul(x_p_mu_Fermion_element.m_d[iSpinor]));

                        //- gammamu U(x,mu) phi(x+ mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].SubC(x_Gauge_element.Mul(gammaMu.MulC(x_p_mu_Fermion_element, iSpinor)));

                        //U^{dagger}(x-mu) phi(x-mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].AddC(x_m_mu_Gauge_element.Mul(x_m_mu_Fermion_element.m_d[iSpinor]));

                        //gammamu U^{dagger}(x-mu) phi(x-mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].AddC(x_m_mu_Gauge_element.Mul(gammaMu.MulC(x_m_mu_Fermion_element, iSpinor)));
                    }
                    else
                    {
                        //U(x,mu) phi(x+ mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].AddC(__SU3Generators[i]->Mulc(x_Gauge_element).Mul(x_p_mu_Fermion_element.m_d[iSpinor]));

                        //- gammamu U(x,mu) phi(x+ mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].SubC(__SU3Generators[i]->Mulc(x_Gauge_element).Mul(gammaMu.MulC(x_p_mu_Fermion_element, iSpinor)));

                        //U^{dagger}(x-mu) phi(x-mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].AddC(__SU3Generators[i]->Mulc(x_m_mu_Gauge_element).Mul(x_m_mu_Fermion_element.m_d[iSpinor]));

                        //gammamu U^{dagger}(x-mu) phi(x-mu)
                        result[i].m_d[iSpinor] = result[i].m_d[iSpinor].AddC(__SU3Generators[i]->Mulc(x_m_mu_Gauge_element).Mul(gammaMu.MulC(x_m_mu_Fermion_element, iSpinor)));
                    }
                }

            }

            for (int i = 0; i < 8; ++i)
            {
                if (!bPartialOmega)
                {
                    //result = phi(x) - kai sum _mu result
                    result[i].Mul(_make_cuComplex(kai, 0));
                    pResultDataArray[linkIndex * 8 + i] = x_Fermion_element.SubC(result[i]);
                    if (bDDagger)
                    {
                        pResultDataArray[linkIndex * 8 + i] = gamma5.MulC(pResultDataArray[linkIndex * 8 + i]);
                    }
                }
                else
                {
                    result[i].Mul(_make_cuComplex(0, -kai));
                    pResultDataArray[linkIndex * 8 + i] = result[i];
                    if (bDDagger)
                    {
                        pResultDataArray[linkIndex * 8 + i] = gamma5.MulC(pResultDataArray[linkIndex * 8 + i]);
                    }
                }
            }
        }
    }
}

#pragma endregion

extern "C"
{
    void _cInitialFermionWilsonSquareSU3(deviceWilsonVectorSU3 *pDevicePtr, EFieldInitialType eInitialType)
    {
        preparethread;
        _kernelInitialFermionWilsonSquareSU3 << <block, threads >> > (pDevicePtr, eInitialType);
    }

    void _cDFermionWilsonSquareSU3(const deviceWilsonVectorSU3* __restrict__ pDeviceData,
                      const deviceSU3* __restrict__ pGauge,
                      deviceWilsonVectorSU3* pResultData,
                      Real kai,
                      BYTE byFieldId,
                      UBOOL bDiracChiralGamma,
                      UBOOL bDDagger)
    {
        preparethread;
        _kernelDFermionWilsonSquareSU3 << <block, threads >> > (pDeviceData, pGauge, pResultData, kai, byFieldId, bDiracChiralGamma, bDDagger);
    }

    void _cDWilsonMuSU3(const deviceWilsonVectorSU3* __restrict__ pDeviceData,
                     const deviceSU3* __restrict__ pGauge,
                     deviceWilsonVectorSU3* pResultDataArray,
                     Real kai,
                     BYTE byFieldId,
                     UBOOL bDiracChiralGamma,
                     UBOOL bDDagger,
                     UBOOL bPartialOmega)
    {
        preparethread;
        _kernelDWilsonMuSU3 << <block, threads >> > (pDeviceData, pGauge, pResultDataArray, kai, byFieldId, bDiracChiralGamma, bDiracChiralGamma, bDDagger, bPartialOmega);
    }


    void _cDotFermionWilsonSquareSU3(const deviceWilsonVectorSU3 * __restrict__ pLeft,
                  const deviceWilsonVectorSU3 * __restrict__ pRight,
                  _Complex* result)
    {
        preparethread;
        _kernelDotFermionWilsonSquareSU3 << <block, threads >> > (pLeft, pRight, result);
    }

    void _cPrintFermionWilsonSquareSU3(const deviceWilsonVectorSU3 * __restrict__ pData)
    {
        preparethread;
        _kernelPrintFermionWilsonSquareSU3 << <block, threads >> > (pData);
    }
}

CFieldFermionWilsonSquareSU3::CFieldFermionWilsonSquareSU3() : CFieldFermion()
{
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceData, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDataCopy, sizeof(deviceWilsonVectorSU3) * m_uiSiteCount));

    checkCudaErrors(cudaMalloc((void**)&m_pForceRightVector, sizeof(deviceWilsonVectorSU3) * m_uiLinkeCount * 8));
    checkCudaErrors(cudaMalloc((void**)&m_pForceRightVectorCopy, sizeof(deviceWilsonVectorSU3) * m_uiLinkeCount * 8));
    checkCudaErrors(cudaMalloc((void**)&m_pForceLeftVector, sizeof(deviceWilsonVectorSU3) * m_uiLinkeCount * 8));
    checkCudaErrors(cudaMalloc((void**)&m_pForceLeftVectorCopy, sizeof(deviceWilsonVectorSU3) * m_uiLinkeCount * 8));
}

CFieldFermionWilsonSquareSU3::~CFieldFermionWilsonSquareSU3()
{
    checkCudaErrors(cudaFree(m_pDeviceData));
    checkCudaErrors(cudaFree(m_pDeviceDataCopy));
    checkCudaErrors(cudaFree(m_pForceRightVector));
    checkCudaErrors(cudaFree(m_pForceRightVectorCopy));
    checkCudaErrors(cudaFree(m_pForceLeftVector));
    checkCudaErrors(cudaFree(m_pForceLeftVectorCopy));
}

/**
*
*/
void CFieldFermionWilsonSquareSU3::InitialField(EFieldInitialType eInitialType)
{
    _cInitialFermionWilsonSquareSU3(m_pDeviceData, eInitialType);
}

void CFieldFermionWilsonSquareSU3::DebugPrintMe() const
{
    _cPrintFermionWilsonSquareSU3(m_pDeviceData);
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

    _cInitialFermionWilsonSquareSU3(m_pDeviceDataCopy, EFIT_RandomGaussian);
    _cDFermionWilsonSquareSU3(m_pDeviceDataCopy, pFieldSU3->m_pDeviceData, m_pDeviceData, m_fKai, m_byFieldId, TRUE, FALSE);
}

void CFieldFermionWilsonSquareSU3::CalculateForce(const CFieldGauge* pGauge, CFieldGauge* pForce)
{

}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================