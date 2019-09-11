//=============================================================================
// FILENAME : CMeasureAMomemtumJG.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [05/21/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureAMomentumJG)

#pragma region kernles

/**
* Initial as zero
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialZero_XYPlane(Real* pBuffer)
{
    intokernalOnlyInt4;

    if (0 == sSite4.z && 0 == sSite4.w)
    {
        pBuffer[sSite4.x * _DC_Ly + sSite4.y] = F(0.0);
    }
}

/**
* Average over z and t
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelAverageOverZT_XYPlane(Real* pBuffer)
{
    intokernalOnlyInt4;

    if (0 == sSite4.z && 0 == sSite4.w)
    {
        const UINT uiIdx = sSite4.x * _DC_Ly + sSite4.y;
        pBuffer[uiIdx] = pBuffer[uiIdx] / (_DC_Lz * _DC_Lt);
    }
}

/**
* calculate momentum, and sum over y and z
*/
__global__ void _CLG_LAUNCH_BOUND
_kernelCalculateAngularMomentumJG(
    const deviceSU3* __restrict__ pDeviceData,
    Real* pBuffer, 
    SSmallInt4 sCenter, 
    Real betaOverN,
    BYTE byFieldId)
{
    intokernalOnlyInt4;

    const UINT uiN = __idx->_deviceGetBigIndex(sSite4);
    Real fRes = F(0.0);
    if (!__idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiN].IsDirichlet())
    {

        //======================================================
        //4-chair terms except for the last one
        betaOverN = F(0.125) * betaOverN;

        const Real fX = (sSite4.x - sCenter.x);

        //===============
        //+x Omega V412
        const Real fV412 = fX * _deviceChairTerm(pDeviceData, 3, 0, 1, uiN);

        //===============
        //+x Omega V432
        const Real fV432 = fX * _deviceChairTerm(pDeviceData, 3, 2, 1, uiN);

        const Real fY = -(sSite4.y - sCenter.y);

        //===============
        //-y Omega V421
        const Real fV421 = fY * _deviceChairTerm(pDeviceData, 3, 1, 0, uiN);

        //===============
        //-y Omega V431
        const Real fV431 = fY * _deviceChairTerm(pDeviceData, 3, 2, 0, uiN);

        fRes = (fV412 + fV432 + fV421 + fV431) * betaOverN;
    }

    atomicAdd(&pBuffer[sSite4.x * _DC_Ly + sSite4.y], fRes);
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelMomentumJGInitialDist(UINT* pCount, Real* pValue)
{
    pCount[threadIdx.x] = 0;
    pValue[threadIdx.x] = F(0.0);
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelMomentumJGMeasureDist(
    const Real* __restrict__ jgXY,
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
        atomicAdd(&correlator[uiC], jgXY[uiXY]);
    }
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelMomentumJGAverageDist(UINT* pCount, Real* pValue)
{
    const UINT uiIdx = threadIdx.x;
    if (pCount[uiIdx] > 0)
    {
        pValue[uiIdx] = pValue[uiIdx] / static_cast<Real>(pCount[uiIdx]);
    }
}

#pragma endregion

/**
* array[x, y] = array[x, y] / (lz * lt)
*/
CLGAPI void _AverageXYPlane(Real* pDeviceRes)
{
    preparethread;
    _kernelAverageOverZT_XYPlane << <block, threads >> > (pDeviceRes);
}

/**
* array[x, y] = 0
*/
CLGAPI void _ZeroXYPlane(Real* pDeviceRes)
{
    preparethread;
    _kernelInitialZero_XYPlane << <block, threads >> > (pDeviceRes);
}

CMeasureAMomentumJG::~CMeasureAMomentumJG()
{
    if (NULL != m_pHostDataBuffer)
    {
        free(m_pHostDataBuffer);
    }
    if (NULL != m_pDeviceDataBufferOneConfig)
    {
        checkCudaErrors(cudaFree(m_pDeviceDataBufferOneConfig));
    }

    if (NULL != m_pDistributionR)
    {
        checkCudaErrors(cudaFree(m_pDistributionR));
    }

    if (NULL != m_pDistributionJG)
    {
        checkCudaErrors(cudaFree(m_pDistributionJG));
    }

    if (NULL != m_pHostDistributionR)
    {
        free(m_pHostDistributionR);
    }

    if (NULL != m_pHostDistributionJG)
    {
        free(m_pHostDistributionJG);
    }
}

void CMeasureAMomentumJG::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    m_pHostDataBuffer = (Real*)malloc(sizeof(Real) * _HC_Lx * _HC_Ly);
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDataBufferOneConfig, sizeof(Real) * _HC_Lx * _HC_Ly));
    Reset();

    INT iValue = 1;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;

    iValue = 1;
    param.FetchValueINT(_T("MeasureDist"), iValue);
    m_bMeasureDistribution = iValue != 0;

    if (m_bMeasureDistribution)
    {
        //assuming the center is really at center
        m_uiMaxR = ((_HC_Lx + 1) / 2) * ((_HC_Lx + 1) / 2)
                 + ((_HC_Ly + 1) / 2) * ((_HC_Ly + 1) / 2);

        m_uiEdgeR = ((_HC_Lx + 1) / 2 - 1) * ((_HC_Lx + 1) / 2 - 1);

        checkCudaErrors(cudaMalloc((void**)&m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1)));
        checkCudaErrors(cudaMalloc((void**)&m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1)));

        m_pHostDistributionR = (UINT*)malloc(sizeof(UINT) * (m_uiMaxR + 1));
        m_pHostDistributionJG = (Real*)malloc(sizeof(Real) * (m_uiMaxR + 1));
    }
}

void CMeasureAMomentumJG::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pGauge);

    const Real fBetaOverN = CCommonData::m_fBeta / static_cast<Real>(_HC_SUN);

    _ZeroXYPlane(m_pDeviceDataBufferOneConfig);

    preparethread;

    _kernelCalculateAngularMomentumJG << <block, threads >> > (
        pGaugeSU3->m_pDeviceData, 
        m_pDeviceDataBufferOneConfig,
        CCommonData::m_sCenter,
        fBetaOverN,
        m_byFieldId);

    _AverageXYPlane(m_pDeviceDataBufferOneConfig);

    if (m_bMeasureDistribution)
    {
        dim3 block2(_HC_DecompX, 1, 1);
        dim3 threads2(_HC_DecompLx, 1, 1);
        dim3 block3(m_uiMaxR + 1, 1, 1);
        dim3 threads3(m_uiMaxR + 1, 1, 1);

        _kernelMomentumJGInitialDist << <block3, threads3 >> >(m_pDistributionR, m_pDistributionJG);

        _kernelMomentumJGMeasureDist << <block2, threads2 >> >(
            m_pDeviceDataBufferOneConfig,
            CCommonData::m_sCenter,
            m_uiMaxR,
            m_byFieldId,
            m_pDistributionR,
            m_pDistributionJG
            );

        _kernelMomentumJGAverageDist << <block3, threads3 >> >(m_pDistributionR, m_pDistributionJG);

        //extract res
        checkCudaErrors(cudaMemcpy(m_pHostDistributionR, m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));
        checkCudaErrors(cudaMemcpy(m_pHostDistributionJG, m_pDistributionJG, sizeof(Real) * (m_uiMaxR + 1), cudaMemcpyDeviceToHost));

        Real fAverageJGInner = F(0.0);
        UINT uiInnerPointsInner = 0;
        Real fAverageJGAll = F(0.0);
        UINT uiInnerPointsAll = 0;
        if (0 == m_uiConfigurationCount)
        {
            assert(0 == m_lstR.Num());
            assert(0 == m_lstJG.Num());

            for (UINT uiL = 0; uiL <= m_uiMaxR; ++uiL)
            {
                if (m_pHostDistributionR[uiL] > 0)
                {
                    m_lstR.AddItem(uiL);
                    m_lstJG.AddItem(m_pHostDistributionJG[uiL]);

                    uiInnerPointsAll += m_pHostDistributionR[uiL];
                    fAverageJGAll += m_pHostDistributionR[uiL] * m_pHostDistributionJG[uiL];
                    if (uiL < m_uiEdgeR)
                    {
                        uiInnerPointsInner += m_pHostDistributionR[uiL];
                        fAverageJGInner += m_pHostDistributionR[uiL] * m_pHostDistributionJG[uiL];
                    }

                    if (m_bShowResult)
                    {
                        appDetailed(_T("JG(%f)=%f, \n"),
                            _hostsqrt(static_cast<Real>(uiL)),
                            m_pHostDistributionJG[uiL]);
                    }
                }
            }
        }
        else
        {
            for (INT i = 0; i < m_lstR.Num(); ++i)
            {
                assert(m_pHostDistributionR[m_lstR[i]] > 0);
                m_lstJG.AddItem(m_pHostDistributionJG[m_lstR[i]]);

                uiInnerPointsAll += m_pHostDistributionR[m_lstR[i]];
                fAverageJGAll += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJG[m_lstR[i]];
                if (m_lstR[i] < m_uiEdgeR)
                {
                    uiInnerPointsInner += m_pHostDistributionR[m_lstR[i]];
                    fAverageJGInner += m_pHostDistributionR[m_lstR[i]] * m_pHostDistributionJG[m_lstR[i]];
                }

                if (m_bShowResult)
                {
                    appDetailed(_T("JG(%f)=%f, \n"),
                        _hostsqrt(static_cast<Real>(m_lstR[i])),
                        m_pHostDistributionJG[m_lstR[i]]);
                }
            }
        }

        if (uiInnerPointsAll > 0)
        {
            fAverageJGAll = fAverageJGAll / uiInnerPointsAll;
        }
        if (uiInnerPointsInner > 0)
        {
            fAverageJGInner = fAverageJGInner / uiInnerPointsInner;
        }
        m_lstJGAll.AddItem(fAverageJGAll);
        m_lstJGInner.AddItem(fAverageJGInner);
    }

    ++m_uiConfigurationCount;
    checkCudaErrors(cudaMemcpy(m_pHostDataBuffer, m_pDeviceDataBufferOneConfig, sizeof(Real) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));

    if (m_bShowResult)
    {
        appDetailed(_T(" === Angular Momentum JG of site y=%d ======\n"), CCommonData::m_sCenter.y);
    }
    for (UINT i = 1; i < _HC_Ly; ++i)
    {
        for (UINT j = 1; j < _HC_Lx; ++j)
        {
            m_lstRes.AddItem(m_pHostDataBuffer[j * _HC_Ly + i]);
        }
    }

    for (UINT i = 1; i < _HC_Lx; ++i)
    {
        if (m_bShowResult)
        {
            appDetailed(_T("(%d,%d)=%1.6f  "), 
                i, 
                CCommonData::m_sCenter.y, 
                m_pHostDataBuffer[i * _HC_Ly + CCommonData::m_sCenter.y]);
        }
    }
    if (m_bShowResult)
    {
        appDetailed(_T("\n"));
    }

}

void CMeasureAMomentumJG::Average(UINT )
{
    //nothing to do
}

void CMeasureAMomentumJG::Report()
{
    appSetLogDate(FALSE);

    assert(m_uiConfigurationCount * (_HC_Lx - 1) * (_HC_Ly - 1) 
        == static_cast<UINT>(m_lstRes.Num()));
    appGeneral(_T("\n===================================================\n"));
    appGeneral(_T("=========== Angular Momentum JG of sites ==========\n"), CCommonData::m_sCenter.x);
    appGeneral(_T("===================================================\n"));

    TArray<Real> tmp;
    appGeneral(_T("{\n"));
    for (UINT k = 0; k < m_uiConfigurationCount; ++k)
    {
        appGeneral(_T("{"));
        for (UINT i = 0; i < _HC_Ly - 1; ++i)
        {
            appGeneral(_T("{"));
            for (UINT j = 0; j < _HC_Lx - 1; ++j)
            {
                const UINT idx = k * (_HC_Lx - 1) * (_HC_Ly - 1) + i * (_HC_Lx - 1) + j;
                appGeneral(_T("%1.6f, "), m_lstRes[idx]);
                if (0 == k)
                {
                    tmp.AddItem(m_lstRes[idx]);
                }
                else
                {
                    tmp[i * (_HC_Lx - 1) + j] += m_lstRes[idx];
                }
            }
            appGeneral(_T("}, "));
        }
        appGeneral(_T("}\n"));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n -------------------- Average -------------------------\n\n"));
    
    for (UINT i = 0; i < _HC_Ly - 1; ++i)
    {
        for (UINT j = 0; j < _HC_Lx - 1; ++j)
        {
            appGeneral(_T("(x=%d,y=%d)%2.8f,   "), j + 1, i + 1, tmp[i * (_HC_Lx - 1) + j] / m_uiConfigurationCount);
        }
        appGeneral(_T("\n"));
    }

    appGeneral(_T("===================================================\n"));
    appGeneral(_T("===================================================\n"));

    appSetLogDate(TRUE);
}

void CMeasureAMomentumJG::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstRes.RemoveAll();

    m_lstR.RemoveAll();
    m_lstJG.RemoveAll();

    m_lstJGAll.RemoveAll();
    m_lstJGInner.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================