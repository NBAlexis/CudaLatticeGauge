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

__global__ void _CLG_LAUNCH_BOUND
_kernelDotAndGatherXY(
    const deviceWilsonVectorSU3 * __restrict__ pMe,
    const deviceWilsonVectorSU3 * __restrict__ pOther,
    CLGComplex * resultXYPlane,
    CLGComplex * result)
{
    intokernal;

    UINT uiXY = threadIdx.x + blockIdx.x * blockDim.x;
    result[uiSiteIndex] = pMe[uiSiteIndex].ConjugateDotC(pOther[uiSiteIndex]);
    atomicAdd(&resultXYPlane[uiXY].x, result[uiSiteIndex].x);
    atomicAdd(&resultXYPlane[uiXY].y, result[uiSiteIndex].y);
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelChiralCondensateInitialDist(UINT* pCount, Real* pValue)
{
    pCount[threadIdx.x] = 0;
    pValue[threadIdx.x] = F(0.0);
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelChiralCondensateMeasureDist(
    const CLGComplex* __restrict__ chiralXY,
    SSmallInt4 sCenter, UINT uiMax, BYTE byFieldId,
    UINT* counter, Real* correlator)
{
    UINT uiXY = (threadIdx.x + blockIdx.x * blockDim.x);
    SBYTE uiX = static_cast<SBYTE>(uiXY / _DC_Ly);
    SBYTE uiY = static_cast<SBYTE>(uiXY % _DC_Ly);
    UINT uiC = (sCenter.x - uiX) * (sCenter.x - uiX)
        + (sCenter.y - uiY) * (sCenter.y - uiY);

    SSmallInt4 sSite4;
    sSite4.z = sCenter.z;
    sSite4.w = sCenter.w;
    sSite4.x = uiX;
    sSite4.y = uiY;
    if (uiC <= uiMax && !__idx->_deviceGetMappingIndex(sSite4, byFieldId).IsDirichlet())
    {
        atomicAdd(&counter[uiC], 1);
        atomicAdd(&correlator[uiC], chiralXY[uiXY].x);
    }
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelChiralAverageDist(UINT* pCount, Real* pValue)
{
    UINT uiIdx = threadIdx.x;
    if (pCount[uiIdx] > 0)
    {
        pValue[uiIdx] = pValue[uiIdx] / static_cast<Real>(pCount[uiIdx]);
    }
}

#pragma endregion

CMeasureChiralCondensate::~CMeasureChiralCondensate()
{
    if (NULL != m_pDeviceXYBuffer)
    {
        checkCudaErrors(cudaFree(m_pDeviceXYBuffer));
        free(m_pHostXYBuffer);
    }

    if (NULL != m_pDistributionR)
    {
        checkCudaErrors(cudaFree(m_pDistributionR));
    }

    if (NULL != m_pDistributionC)
    {
        checkCudaErrors(cudaFree(m_pDistributionC));
    }

    if (NULL != m_pHostDistributionR)
    {
        free(m_pHostDistributionR);
    }

    if (NULL != m_pHostDistributionC)
    {
        free(m_pHostDistributionC);
    }
}

void CMeasureChiralCondensate::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    checkCudaErrors(cudaMalloc((void**)&m_pDeviceXYBuffer, sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    m_pHostXYBuffer = (CLGComplex*)malloc(sizeof(CLGComplex) * _HC_Lx * _HC_Ly);

    Reset();

    INT iValue = 2;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("FieldCount"), iValue);
    m_uiFieldCount = static_cast<UINT>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;

    iValue = 1;
    param.FetchValueINT(_T("MeasureDist"), iValue);
    m_bMeasureDistribution = iValue != 0;

    if (m_bMeasureDistribution)
    {
        //assuming the center is really at center
        m_uiMaxR = ((_HC_Lx + 1) / 2 - 1) * ((_HC_Lx + 1) / 2 - 1)
            + ((_HC_Ly + 1) / 2 - 1) * ((_HC_Ly + 1) / 2 - 1);

        checkCudaErrors(cudaMalloc((void**)&m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1)));
        checkCudaErrors(cudaMalloc((void**)&m_pDistributionC, sizeof(Real) * (m_uiMaxR + 1)));

        m_pHostDistributionR = (UINT*)malloc(sizeof(UINT) * (m_uiMaxR + 1));
        m_pHostDistributionC = (Real*)malloc(sizeof(Real) * (m_uiMaxR + 1));
    }
}

void CMeasureChiralCondensate::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{
    CLGComplex res = _make_cuComplex(F(0.0), F(0.0));

    CFieldFermion* pF1 = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    CFieldFermion* pF2 = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    UINT uiVolume = appGetLattice()->m_pIndexCache->m_uiSiteNumber[m_byFieldId];
    CFieldFermionWilsonSquareSU3 * pF1W = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pF1);
    CFieldFermionWilsonSquareSU3 * pF2W = dynamic_cast<CFieldFermionWilsonSquareSU3*>(pF2);

    _ZeroXYPlaneC(m_pDeviceXYBuffer);
    for (UINT i = 0; i < m_uiFieldCount; ++i)
    {
        pF1->InitialField(EFIT_RandomZ4);
        pF1->FixBoundary();
        //pF1->DebugPrintMe();
        //pGauge->DebugPrintMe();
        pF1->CopyTo(pF2);
        pF1->InverseD(pGauge);       

#pragma region Dot

        preparethread;
        _kernelDotAndGatherXY << <block, threads >> > (
            pF1W->m_pDeviceData,
            pF2W->m_pDeviceData,
            m_pDeviceXYBuffer,
            _D_ComplexThreadBuffer);

        CLGComplex thisSum = appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
#pragma endregion

        res.x = res.x + thisSum.x / uiVolume;
        res.y = res.y + thisSum.y / uiVolume;

    }

    if (m_bMeasureDistribution)
    {
        dim3 block2(_HC_DecompX, 1, 1);
        dim3 threads2(_HC_DecompLx, 1, 1);
        dim3 block3(m_uiMaxR + 1, 1, 1);
        dim3 threads3(m_uiMaxR + 1, 1, 1);

        _kernelChiralCondensateInitialDist << <block3, threads3 >> >(m_pDistributionR, m_pDistributionC);

        _kernelChiralCondensateMeasureDist << <block2, threads2 >> >(
            m_pDeviceXYBuffer,
            CCommonData::m_sCenter,
            m_uiMaxR,
            m_byFieldId,
            m_pDistributionR,
            m_pDistributionC
            );

        _kernelChiralAverageDist << <block3, threads3 >> >(m_pDistributionR, m_pDistributionC);

        //extract res
        checkCudaErrors(cudaMemcpy(m_pHostDistributionR, m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(m_pHostDistributionC, m_pDistributionC, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));

        if (0 == m_uiConfigurationCount)
        {
            assert(0 == m_lstR.Num());
            assert(0 == m_lstC.Num());

            for (UINT uiL = 0; uiL <= m_uiMaxR; ++uiL)
            {
                if (m_pHostDistributionR[uiL] > 0)
                {
                    m_lstR.AddItem(uiL);
                    m_lstC.AddItem(m_pHostDistributionC[uiL] / (m_uiFieldCount * _HC_Lz * _HC_Lt));

                    if (m_bShowResult)
                    {
                        appDetailed(_T("C(%f)=%f, \n"),
                            _hostsqrt(static_cast<Real>(uiL)),
                            m_pHostDistributionC[uiL]);
                    }
                }
            }
        }
        else
        {
            for (INT i = 0; i < m_lstR.Num(); ++i)
            {
                assert(m_pHostDistributionR[m_lstR[i]] > 0);
                m_lstC.AddItem(m_pHostDistributionC[m_lstR[i]] / (m_uiFieldCount * _HC_Lz * _HC_Lt));

                if (m_bShowResult)
                {
                    appDetailed(_T("C(%f)=%f, \n"),
                        _hostsqrt(static_cast<Real>(m_lstR[i])),
                        m_pHostDistributionC[m_lstR[i]]);
                }
            }
        }
    }

    checkCudaErrors(cudaMemcpy(m_pHostXYBuffer, m_pDeviceXYBuffer, sizeof(CLGComplex) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
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

    pF1->Return();
    pF2->Return();
    res.x = res.x / m_uiFieldCount;
    res.y = res.y / m_uiFieldCount;
    appDetailed(_T("\nChiral Condensate = %2.12f + %2.12f\n"), res.x, res.y);
    ++m_uiConfigurationCount;
    m_lstCondensate.AddItem(res);
}

void CMeasureChiralCondensate::Average(UINT )
{
    //nothing to do
}

void CMeasureChiralCondensate::Report()
{
    assert(m_uiConfigurationCount == static_cast<UINT>(m_lstCondensate.Num()));
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
            tmpChargeSum.x += m_lstCondensate[i].x;
            tmpChargeSum.y += m_lstCondensate[i].y;
            LogGeneralComplex(m_lstCondensate[i]);
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
            m_lstCondensate[0].x,
            m_lstCondensate[0].y);

        m_cAverageCondensate = m_lstCondensate[0];
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
    m_lstCondensate.RemoveAll();
    m_lstCondensateDensity.RemoveAll();

    m_lstR.RemoveAll();
    m_lstC.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================