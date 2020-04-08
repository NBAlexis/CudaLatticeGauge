//=============================================================================
// FILENAME : CMeasureChiralCondensate.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [06/13/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureChiralCondensate)

#pragma region kernels

/**
 * Psi^bar Psi
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDotAndGatherXYChiral(
    const deviceWilsonVectorSU3 * __restrict__ pMe,
    const deviceWilsonVectorSU3 * __restrict__ pOther,
    CLGComplex* resultXYPlaneChiral,
    CLGComplex * result)
{
    intokernal;

    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
    result[uiSiteIndex] = pMe[uiSiteIndex].ConjugateDotC(pOther[uiSiteIndex]);
    atomicAdd(&resultXYPlaneChiral[uiXY].x, result[uiSiteIndex].x);
    atomicAdd(&resultXYPlaneChiral[uiXY].y, result[uiSiteIndex].y);
}

/**
 * Psi^bar gamma5 Psi
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDotAndGatherXYPion(
    const deviceWilsonVectorSU3* __restrict__ pMe,
    const deviceWilsonVectorSU3* __restrict__ pOther,
    CLGComplex* resultXYPlanePion,
    CLGComplex* result)
{
    intokernal;

    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
    result[uiSiteIndex] = pMe[uiSiteIndex].ConjugateDotC(__chiralGamma[GAMMA5].MulWilsonC(pOther[uiSiteIndex]));
    atomicAdd(&resultXYPlanePion[uiXY].x, result[uiSiteIndex].x);
    atomicAdd(&resultXYPlanePion[uiXY].y, result[uiSiteIndex].y);
}

/**
 * Psi^bar gamma_4 Psi
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDotAndGatherXYRhon(
    const deviceWilsonVectorSU3* __restrict__ pMe,
    const deviceWilsonVectorSU3* __restrict__ pOther,
    CLGComplex* resultXYPlaneRhon,
    CLGComplex* result)
{
    intokernal;

    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
    result[uiSiteIndex] = pMe[uiSiteIndex].ConjugateDotC(__chiralGamma[GAMMA4].MulWilsonC(pOther[uiSiteIndex]));
    atomicAdd(&resultXYPlaneRhon[uiXY].x, result[uiSiteIndex].x);
    atomicAdd(&resultXYPlaneRhon[uiXY].y, result[uiSiteIndex].y);
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelChiralCondensateInitialDist(UINT* pCount, CLGComplex* pChiral, CLGComplex* pPion, CLGComplex* pRhon)
{
    pCount[threadIdx.x] = 0;
    pChiral[threadIdx.x] = _zeroc;
    pPion[threadIdx.x] = _zeroc;
    pRhon[threadIdx.x] = _zeroc;
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelChiralCondensateMeasureDist(
    const CLGComplex* __restrict__ chiralXY,
    const CLGComplex* __restrict__ pionXY,
    const CLGComplex* __restrict__ rhonXY,
    SSmallInt4 sCenter, UINT uiMax, BYTE byFieldId,
    UINT* counter, 
    CLGComplex* chiral,
    CLGComplex* pion,
    CLGComplex* rhon
)
{
    UINT uiXY = (threadIdx.x + blockIdx.x * blockDim.x);
    INT uiX = static_cast<SBYTE>(uiXY / _DC_Ly);
    INT uiY = static_cast<SBYTE>(uiXY % _DC_Ly);
    UINT uiC = (sCenter.x - uiX) * (sCenter.x - uiX)
        + (sCenter.y - uiY) * (sCenter.y - uiY);

    SSmallInt4 sSite4;
    sSite4.z = sCenter.z;
    sSite4.w = sCenter.w;
    sSite4.x = static_cast<SBYTE>(uiX);
    sSite4.y = static_cast<SBYTE>(uiY);
    if (uiC <= uiMax && !__idx->_deviceGetMappingIndex(sSite4, byFieldId).IsDirichlet())
    {
        atomicAdd(&counter[uiC], 1);
        atomicAdd(&chiral[uiC].x, chiralXY[uiXY].x);
        atomicAdd(&pion[uiC].x, pionXY[uiXY].x);
        atomicAdd(&rhon[uiC].x, rhonXY[uiXY].x);
        atomicAdd(&chiral[uiC].y, chiralXY[uiXY].y);
        atomicAdd(&pion[uiC].y, pionXY[uiXY].y);
        atomicAdd(&rhon[uiC].y, rhonXY[uiXY].y);
    }
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelChiralAverageDist(UINT* pCount, CLGComplex* pChiral, CLGComplex* pPion, CLGComplex* pRhon)
{
    const UINT uiIdx = threadIdx.x;
    if (pCount[uiIdx] > 0)
    {
        pChiral[uiIdx].x = pChiral[uiIdx].x / static_cast<Real>(pCount[uiIdx]);
        pChiral[uiIdx].y = pChiral[uiIdx].y / static_cast<Real>(pCount[uiIdx]);
        pPion[uiIdx].x = pPion[uiIdx].x / static_cast<Real>(pCount[uiIdx]);
        pPion[uiIdx].y = pPion[uiIdx].y / static_cast<Real>(pCount[uiIdx]);
        pRhon[uiIdx].x = pRhon[uiIdx].x / static_cast<Real>(pCount[uiIdx]);
        pRhon[uiIdx].y = pRhon[uiIdx].y / static_cast<Real>(pCount[uiIdx]);
    }
}

#pragma endregion

CMeasureChiralCondensate::~CMeasureChiralCondensate()
{
    if (NULL != m_pDeviceXYBufferChiral)
    {
        checkCudaErrors(cudaFree(m_pDeviceXYBufferChiral));
        checkCudaErrors(cudaFree(m_pDeviceXYBufferPion));
        checkCudaErrors(cudaFree(m_pDeviceXYBufferRhon));
        free(m_pHostXYBuffer);
    }

    if (NULL != m_pDistributionR)
    {
        checkCudaErrors(cudaFree(m_pDistributionR));
    }

    if (NULL != m_pDistributionChiral)
    {
        checkCudaErrors(cudaFree(m_pDistributionChiral));
        checkCudaErrors(cudaFree(m_pDistributionPion));
        checkCudaErrors(cudaFree(m_pDistributionRhon));
    }

    if (NULL != m_pHostDistributionR)
    {
        free(m_pHostDistributionR);
    }

    if (NULL != m_pHostDistributionChiral)
    {
        free(m_pHostDistributionChiral);
        free(m_pHostDistributionPion);
        free(m_pHostDistributionRhon);
    }
}

void CMeasureChiralCondensate::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasureStochastic::Initial(pOwner, pLatticeData, param, byId);

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceXYBufferChiral, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceXYBufferPion, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceXYBufferRhon, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    m_pHostXYBuffer = (CLGComplex*)malloc(sizeof(CLGComplex) * _HC_Lx * _HC_Ly);

    Reset();

    INT iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;

    iValue = 1;
    param.FetchValueINT(_T("MeasureDist"), iValue);
    m_bMeasureDistribution = iValue != 0;

    if (m_bMeasureDistribution)
    {
        //assuming the center is really at center
        m_uiMaxR = ((_HC_Lx + 1) / 2 ) * ((_HC_Lx + 1) / 2 )
            + ((_HC_Ly + 1) / 2 ) * ((_HC_Ly + 1) / 2 );

        checkCudaErrors(cudaMalloc((void**)&m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1)));
        checkCudaErrors(cudaMalloc((void**)&m_pDistributionChiral, sizeof(CLGComplex) * (m_uiMaxR + 1)));
        checkCudaErrors(cudaMalloc((void**)&m_pDistributionPion, sizeof(CLGComplex) * (m_uiMaxR + 1)));
        checkCudaErrors(cudaMalloc((void**)&m_pDistributionRhon, sizeof(CLGComplex) * (m_uiMaxR + 1)));

        m_pHostDistributionR = (UINT*)malloc(sizeof(UINT) * (m_uiMaxR + 1));
        m_pHostDistributionChiral = (CLGComplex*)malloc(sizeof(CLGComplex) * (m_uiMaxR + 1));
        m_pHostDistributionPion = (CLGComplex*)malloc(sizeof(CLGComplex) * (m_uiMaxR + 1));
        m_pHostDistributionRhon = (CLGComplex*)malloc(sizeof(CLGComplex) * (m_uiMaxR + 1));
    }
}

void CMeasureChiralCondensate::OnConfigurationAcceptedZ4(
    const class CFieldGauge* pAcceptGauge, 
    const class CFieldGauge* pCorrespondingStaple, 
    const class CFieldFermion* pZ4, 
    const class CFieldFermion* pInverseZ4, 
    UBOOL bStart, 
    UBOOL bEnd)
{
    if (bStart)
    {
        _ZeroXYPlaneC(m_pDeviceXYBufferChiral);
        _ZeroXYPlaneC(m_pDeviceXYBufferPion);
        _ZeroXYPlaneC(m_pDeviceXYBufferRhon);
        m_cTmpSumChiral = _zeroc;
        m_cTmpSumPion = _zeroc;
        m_cTmpSumRhon = _zeroc;
    }

    const Real oneOuiVolume = F(1.0) / appGetLattice()->m_pIndexCache->m_uiSiteNumber[m_byFieldId];
    const CFieldFermionWilsonSquareSU3 * pF1W = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(pInverseZ4);
    const CFieldFermionWilsonSquareSU3 * pF2W = dynamic_cast<const CFieldFermionWilsonSquareSU3*>(pZ4);   

    
#pragma region Dot

    // The results are Atomic Add to m_pDeviceXYBuffer

    preparethread;
    _kernelDotAndGatherXYChiral << <block, threads >> > (
        pF1W->m_pDeviceData,
        pF2W->m_pDeviceData,
        m_pDeviceXYBufferChiral,
        _D_ComplexThreadBuffer);

    const CLGComplex thisSumChiral = appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);

    m_cTmpSumChiral = _cuCaddf(m_cTmpSumChiral, cuCmulf_cr(thisSumChiral, oneOuiVolume));

    _kernelDotAndGatherXYPion << <block, threads >> > (
        pF1W->m_pDeviceData,
        pF2W->m_pDeviceData,
        m_pDeviceXYBufferPion,
        _D_ComplexThreadBuffer);

    const CLGComplex thisSumPion = appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);

    m_cTmpSumPion = _cuCaddf(m_cTmpSumPion, cuCmulf_cr(thisSumPion, oneOuiVolume));

    _kernelDotAndGatherXYRhon << <block, threads >> > (
        pF1W->m_pDeviceData,
        pF2W->m_pDeviceData,
        m_pDeviceXYBufferRhon,
        _D_ComplexThreadBuffer);

    const CLGComplex thisSumRhon = appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);

    m_cTmpSumRhon = _cuCaddf(m_cTmpSumRhon, cuCmulf_cr(thisSumRhon, oneOuiVolume));

#pragma endregion

    if (bEnd)
    {
        if (m_bMeasureDistribution)
        {
            dim3 block2(_HC_DecompX, 1, 1);
            dim3 threads2(_HC_DecompLx, 1, 1);
            dim3 block3(m_uiMaxR + 1, 1, 1);
            dim3 threads3(m_uiMaxR + 1, 1, 1);

            _kernelChiralCondensateInitialDist << <block3, threads3 >> >(m_pDistributionR, 
                m_pDistributionChiral,
                m_pDistributionPion,
                m_pDistributionRhon);

            _kernelChiralCondensateMeasureDist << <block2, threads2 >> >(
                m_pDeviceXYBufferChiral,
                m_pDeviceXYBufferPion,
                m_pDeviceXYBufferRhon,
                CCommonData::m_sCenter,
                m_uiMaxR,
                m_byFieldId,
                m_pDistributionR,
                m_pDistributionChiral,
                m_pDistributionPion,
                m_pDistributionRhon
                );

            _kernelChiralAverageDist << <block3, threads3 >> >(m_pDistributionR, 
                m_pDistributionChiral,
                m_pDistributionPion,
                m_pDistributionRhon);

            //extract res
            checkCudaErrors(cudaMemcpy(m_pHostDistributionR, m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
            checkCudaErrors(cudaMemcpy(m_pHostDistributionChiral, m_pDistributionChiral, sizeof(CLGComplex) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
            checkCudaErrors(cudaMemcpy(m_pHostDistributionPion, m_pDistributionPion, sizeof(CLGComplex) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
            checkCudaErrors(cudaMemcpy(m_pHostDistributionRhon, m_pDistributionRhon, sizeof(CLGComplex) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));

            const Real fDivider = F(1.0) / (m_uiFieldCount * _HC_Lz * _HC_Lt);
            if (0 == m_uiConfigurationCount)
            {
                assert(0 == m_lstR.Num());
                assert(0 == m_lstChiral.Num());
                assert(0 == m_lstPion.Num());
                assert(0 == m_lstRhon.Num());
                for (UINT uiL = 0; uiL <= m_uiMaxR; ++uiL)
                {
                    if (m_pHostDistributionR[uiL] > 0)
                    {
                        m_lstR.AddItem(uiL);
                        m_lstChiral.AddItem(cuCmulf_cr(m_pHostDistributionChiral[uiL], fDivider));
                        m_lstPion.AddItem(cuCmulf_cr(m_pHostDistributionPion[uiL], fDivider));
                        m_lstRhon.AddItem(cuCmulf_cr(m_pHostDistributionRhon[uiL], fDivider));

                        if (m_bShowResult)
                        {
                            appDetailed(_T("C(%f)=%f, %f, %f\n"),
                                _hostsqrt(static_cast<Real>(uiL)),
                                m_pHostDistributionChiral[uiL].x,
                                m_pHostDistributionPion[uiL].x,
                                m_pHostDistributionRhon[uiL].x
                            );
                        }
                    }
                }
            }
            else
            {
                for (INT i = 0; i < m_lstR.Num(); ++i)
                {
                    assert(m_pHostDistributionR[m_lstR[i]] > 0);
                    m_lstChiral.AddItem(cuCmulf_cr(m_pHostDistributionChiral[m_lstR[i]], fDivider));
                    m_lstPion.AddItem(cuCmulf_cr(m_pHostDistributionPion[m_lstR[i]], fDivider));
                    m_lstRhon.AddItem(cuCmulf_cr(m_pHostDistributionRhon[m_lstR[i]], fDivider));

                    if (m_bShowResult)
                    {
                        appDetailed(_T("C(%f)=%f, %f, %f\n"),
                            _hostsqrt(static_cast<Real>(m_lstR[i])),
                            m_pHostDistributionChiral[m_lstR[i]].x,
                            m_pHostDistributionPion[m_lstR[i]].x,
                            m_pHostDistributionRhon[m_lstR[i]].x
                        );
                    }
                }
            }
        }

        //we in fact don't care about XY distribution now...
        checkCudaErrors(cudaMemcpy(m_pHostXYBuffer, m_pDeviceXYBufferChiral, sizeof(CLGComplex) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
        if (m_bShowResult)
        {
            appDetailed(_T("\n ------ Densisty -----\n"));
        }

        for (UINT i = static_cast<UINT>(CCommonData::m_sCenter.x); i < _HC_Lx; ++i)
        {
            CLGComplex cvalue = m_pHostXYBuffer[i * _HC_Ly + CCommonData::m_sCenter.y];
            cvalue.x = cvalue.x / (m_uiFieldCount * _HC_Lz * _HC_Lt);
            cvalue.y = cvalue.y / (m_uiFieldCount * _HC_Lz * _HC_Lt);
            m_lstCondensateDensity.AddItem(cvalue);
            if (m_bShowResult)
            {
                appDetailed(_T("(%d,%d)=%1.6f %s %1.6f I   "), i, CCommonData::m_sCenter.y,
                    cvalue.x,
                    cvalue.y < F(0.0) ? _T("") : _T("+"),
                    appAbs(cvalue.y));
            }
        }
        if (m_bShowResult)
        {
            appDetailed(_T("\n ------ Densisty -----\n"));
        }

        const Real fDiv2 = F(1.0) / m_uiFieldCount;
        m_cTmpSumChiral = cuCmulf_cr(m_cTmpSumChiral, fDiv2);
        m_cTmpSumPion = cuCmulf_cr(m_cTmpSumPion, fDiv2);
        m_cTmpSumRhon = cuCmulf_cr(m_cTmpSumRhon, fDiv2);

        appDetailed(_T("\nChiral Condensate = %2.12f + %2.12f\n"), m_cTmpSumChiral.x, m_cTmpSumChiral.y);
        appDetailed(_T("\nPion Condensate = %2.12f + %2.12f\n"), m_cTmpSumPion.x, m_cTmpSumPion.y);
        appDetailed(_T("\nRhon Condensate = %2.12f + %2.12f\n"), m_cTmpSumRhon.x, m_cTmpSumRhon.y);
        ++m_uiConfigurationCount;
        m_lstChiralAll.AddItem(m_cTmpSumChiral);
        m_lstPionAll.AddItem(m_cTmpSumPion);
        m_lstRhonAll.AddItem(m_cTmpSumRhon);
    }
}

void CMeasureChiralCondensate::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{

}

void CMeasureChiralCondensate::Average(UINT )
{
    //nothing to do
}

void CMeasureChiralCondensate::Report()
{
    assert(m_uiConfigurationCount == static_cast<UINT>(m_lstChiralAll.Num()));
    assert(m_uiConfigurationCount == static_cast<UINT>(m_lstPionAll.Num()));
    assert(m_uiConfigurationCount == static_cast<UINT>(m_lstRhonAll.Num()));

    assert(static_cast<UINT>(m_uiConfigurationCount * CCommonData::m_sCenter.x)
        == static_cast<UINT>(m_lstCondensateDensity.Num()));

    appSetLogDate(FALSE);
    CLGComplex tmpChargeSum = _make_cuComplex(F(0.0), F(0.0));
    m_lstAverageCondensateDensity.RemoveAll();

    appGeneral(_T("\n==========================================================================\n"));
    appGeneral(_T("==================== Chiral Condensate (%d con)============================\n"), m_uiConfigurationCount);

    if (m_uiConfigurationCount > 1)
    {
        appGeneral(_T("\n ----------- each configuration ------------- \n"));
        appGeneral(_T("{"));

        for (UINT i = 0; i < m_uiConfigurationCount; ++i)
        {
            tmpChargeSum.x += m_lstChiralAll[i].x;
            tmpChargeSum.y += m_lstChiralAll[i].y;
            LogGeneralComplex(m_lstChiralAll[i]);
        }
        appGeneral(_T("}\n"));

        tmpChargeSum.x = tmpChargeSum.x / m_uiConfigurationCount;
        tmpChargeSum.y = tmpChargeSum.y / m_uiConfigurationCount;
        appGeneral(_T("\n ----------- average condensate = %2.12f + %2.12f ------------- \n"),
            tmpChargeSum.x, tmpChargeSum.y);

        m_cAverageCondensate = tmpChargeSum;
    }
    else
    {
        appGeneral(_T("\n ----------- average condensate = %2.12f + %2.12f ------------- \n"),
            m_lstChiralAll[0].x,
            m_lstChiralAll[0].y);

        m_cAverageCondensate = m_lstChiralAll[0];
    }

    appGeneral(_T("\n ----------- condensate density------------- \n"));
    appGeneral(_T("{\n"));
    for (UINT k = 0; k < m_uiConfigurationCount; ++k)
    {
        appGeneral(_T("{"));
        for (UINT i = 0; i < static_cast<UINT>(CCommonData::m_sCenter.x); ++i)
        {
            LogGeneralComplex(m_lstCondensateDensity[k * CCommonData::m_sCenter.x + i]);

            if (0 == k)
            {
                m_lstAverageCondensateDensity.AddItem(m_lstCondensateDensity[k * CCommonData::m_sCenter.x + i]);
            }
            else
            {
                m_lstAverageCondensateDensity[i] = _cuCaddf(m_lstAverageCondensateDensity[i], m_lstCondensateDensity[k * CCommonData::m_sCenter.x + i]);
            }

            if (k == m_uiConfigurationCount - 1)
            {
                m_lstAverageCondensateDensity[i].x = m_lstAverageCondensateDensity[i].x / m_uiConfigurationCount;
                m_lstAverageCondensateDensity[i].y = m_lstAverageCondensateDensity[i].y / m_uiConfigurationCount;
            }
        }
        appGeneral(_T("},\n"));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("==========================================================================\n"));
    appSetLogDate(TRUE);
}

void CMeasureChiralCondensate::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstChiralAll.RemoveAll();
    m_lstPionAll.RemoveAll();
    m_lstRhonAll.RemoveAll();

    m_lstCondensateDensity.RemoveAll();

    m_lstR.RemoveAll();
    m_lstChiral.RemoveAll();
    m_lstPion.RemoveAll();
    m_lstRhon.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================