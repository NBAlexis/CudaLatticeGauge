//=============================================================================
// FILENAME : CMeasurePolyakov.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [07/07/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasurePolyakov)

#pragma region kernles 

__global__ void
_CLG_LAUNCH_BOUND
_kernelPolyakovTraceOfSite(
    const deviceSU3* __restrict__ resXYZ,
    CLGComplex* traceXYZ,
    CLGComplex* resSum)
{
    UINT uiXYZ = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    CLGComplex trres = resXYZ[uiXYZ].Tr();
    traceXYZ[uiXYZ] = trres;
    atomicAdd(&resSum[0].x, trres.x);
    atomicAdd(&resSum[0].y, trres.y);
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelInitialStaticPotentialCorrelatorOfSite(UINT* counter, CLGComplex* correlator)
{
    counter[threadIdx.x] = 0;
    correlator[threadIdx.x] = _make_cuComplex(F(0.0), F(0.0));
}

/**
* C(n)=P(n)P(center)^+
* c = x^2+y^2+z^2
* counter[c] = counter[c] + 1
* correlator[c] = correlator[c] + C(n)
*/
__global__ void
_CLG_LAUNCH_BOUND
_kernelStaticPotentialCorrelatorOfSite(
    const CLGComplex* __restrict__ traceXYZ,
    SSmallInt4 sCenter, UINT uiMax,
    UINT* counter,  CLGComplex* correlator)
{
    UINT uiXY = (threadIdx.x + blockIdx.x * blockDim.x);
    INT uiX = static_cast<INT>(uiXY / _DC_Ly);
    INT uiY = static_cast<INT>(uiXY % _DC_Ly);
    INT uiZ = static_cast<INT>(threadIdx.y + blockIdx.y * blockDim.y);
    UINT uiXYZ = uiXY * _DC_Lz + uiZ;
    const UINT uiCenter = (sCenter.x * _DC_Ly + sCenter.y) * _DC_Lz + sCenter.z;
    INT uiC = (sCenter.x - uiX) * (sCenter.x - uiX) 
             + (sCenter.y - uiY) * (sCenter.y - uiY) 
             + (sCenter.z - uiZ) * (sCenter.z - uiZ);

    if (uiC < uiMax)
    {
        CLGComplex correlatorres = _cuCmulf(traceXYZ[uiXYZ], _cuConjf(traceXYZ[uiCenter]));
        atomicAdd(&counter[uiC], 1);
        atomicAdd(&correlator[uiC].x, correlatorres.x);
        atomicAdd(&correlator[uiC].y, correlatorres.y);
    }
}

/**
 * Block.x = x * Ly + y, Thread.x = Lz
 * Block.y = shift xy, Thread.z = shift z
 */
__global__ void
_CLG_LAUNCH_BOUND
_kernelStaticPotentialCorrelatorOfSite2(
    const CLGComplex* __restrict__ traceXYZ, UINT uiMax,
    UINT* counter, CLGComplex* correlator)
{
    INT uiX = static_cast<INT>(blockIdx.x / _DC_Ly);
    INT uiY = static_cast<INT>(blockIdx.x % _DC_Ly);
    INT uiZ = static_cast<INT>(threadIdx.x);

    INT uiShiftX = static_cast<INT>(blockIdx.y / _DC_Ly) - static_cast<INT>(_DC_Lx >> 1);
    INT uiShiftY = static_cast<INT>(blockIdx.y % _DC_Ly) - static_cast<INT>(_DC_Ly >> 1);
    INT uiShiftZ = static_cast<INT>(threadIdx.y);

    INT uiX2 = (uiX + uiShiftX + static_cast<INT>(_DC_Lx)) % static_cast<INT>(_DC_Lx);
    INT uiY2 = (uiY + uiShiftY + static_cast<INT>(_DC_Ly)) % static_cast<INT>(_DC_Ly);
    INT uiZ2 = (uiZ + uiShiftZ + static_cast<INT>(_DC_Lz)) % static_cast<INT>(_DC_Lz);

    INT uiXYZ1 = blockIdx.x * _DC_Lz + uiZ;
    INT uiXYZ2 = ((uiX2 * _DC_Ly) + uiY2) * _DC_Lz + uiZ2;

    INT uiC = (uiX2 - uiX) * (uiX2 - uiX)
            + (uiY2 - uiY) * (uiY2 - uiY)
            + (uiZ2 - uiZ) * (uiZ2 - uiZ);

    if (uiC < uiMax)
    {
        CLGComplex correlatorres = _cuCmulf(traceXYZ[uiXYZ1], _cuConjf(traceXYZ[uiXYZ2]));
        atomicAdd(&counter[uiC], 1);
        atomicAdd(&correlator[uiC].x, correlatorres.x);
        atomicAdd(&correlator[uiC].y, correlatorres.y);
    }
}

__global__ void
_CLG_LAUNCH_BOUND
_kernelAverageStaticPotentialCorrelatorOfSite(UINT* counter, CLGComplex* correlator)
{
    const UINT uiIdx = threadIdx.x + blockIdx.x * blockDim.x;
    if (counter[uiIdx] > 0)
    {
        correlator[uiIdx].x =
            correlator[uiIdx].x / static_cast<Real>(counter[uiIdx]);
        correlator[uiIdx].y =
            correlator[uiIdx].y / static_cast<Real>(counter[uiIdx]);
    }
}

#pragma endregion

CMeasurePolyakov::~CMeasurePolyakov()
{
    if (NULL != m_pHostCorrelator)
    {
        free(m_pHostCorrelator);
    }

    if (NULL != m_pHostCorrelatorCounter)
    {
        free(m_pHostCorrelatorCounter);
    }

    if (NULL != m_pTmpDeviceSum)
    {
        checkCudaErrors(cudaFree(m_pTmpDeviceSum));
    }

    if (NULL != m_pTraceRes)
    {
        checkCudaErrors(cudaFree(m_pTraceRes));
    }

    if (NULL != m_pTmpLoop)
    {
        checkCudaErrors(cudaFree(m_pTmpLoop));
    }

    if (NULL != m_pCorrelator)
    {
        checkCudaErrors(cudaFree(m_pCorrelator));
    }

    if (NULL != m_pCorrelatorCounter)
    {
        checkCudaErrors(cudaFree(m_pCorrelatorCounter));
    }
}

void CMeasurePolyakov::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    checkCudaErrors(cudaMalloc((void**)&m_pTmpDeviceSum, sizeof(CLGComplex)));
    checkCudaErrors(cudaMalloc((void**)&m_pTmpLoop, sizeof(deviceSU3) * _HC_Lx * _HC_Ly * _HC_Lz));
    checkCudaErrors(cudaMalloc((void**)&m_pTraceRes, sizeof(CLGComplex) * _HC_Lx * _HC_Ly * _HC_Lz));

    //We assume the center is really at center
    m_uiMaxLengthSq = ((_HC_Lx + 1) / 2 - 1) * ((_HC_Lx + 1) / 2 - 1)
                    + ((_HC_Ly + 1) / 2 - 1) * ((_HC_Ly + 1) / 2 - 1)
                    + ((_HC_Lz + 1) / 2 - 1) * ((_HC_Lz + 1) / 2 - 1)
                    + 1;

    if (m_uiMaxLengthSq > _HC_ThreadConstraint)
    {
        m_uiMaxLengthSq = _HC_ThreadConstraint;
    }
    if (m_uiMaxLengthSq > _HC_ThreadConstraintX)
    {
        m_uiMaxLengthSq = _HC_ThreadConstraintX;
    }

    checkCudaErrors(cudaMalloc((void**)&m_pCorrelator, sizeof(CLGComplex) * m_uiMaxLengthSq));
    checkCudaErrors(cudaMalloc((void**)&m_pCorrelatorCounter, sizeof(UINT) * m_uiMaxLengthSq));
    
    m_pHostCorrelator = (CLGComplex*)(malloc(sizeof(CLGComplex) * m_uiMaxLengthSq));
    m_pHostCorrelatorCounter = (UINT*)(malloc(sizeof(UINT) * m_uiMaxLengthSq));

    Reset();

    INT iValue = 1;
    param.FetchValueINT(_T("FieldId"), iValue);
    m_byFieldId = static_cast<BYTE>(iValue);

    iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;
}

void CMeasurePolyakov::OnConfigurationAccepted(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pAcceptGauge || EFT_GaugeSU3 != pAcceptGauge->GetFieldType())
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SU3!\n"));
        return;
    }
    const CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<const CFieldGaugeSU3*>(pAcceptGauge);

    //_HC_DecompX * _HC_DecompLx =  Lx * Ly
    //_HC_DecompY * _HC_DecompLy = Lz
    dim3 block1(_HC_DecompX, _HC_DecompY, 1); 
    dim3 threads1(_HC_DecompLx, _HC_DecompLy, 1);

    dim3 block2(1, 1, 1);
    dim3 threads2(m_uiMaxLengthSq, 1, 1);

    CLGComplex res[1];
    res[0] = _make_cuComplex(F(0.0), F(0.0));
    checkCudaErrors(cudaMemcpy(m_pTmpDeviceSum, res, sizeof(CLGComplex), cudaMemcpyHostToDevice));

    _PolyakovAtSite(pGaugeSU3->m_pDeviceData, m_pTmpLoop);
    _kernelPolyakovTraceOfSite << <block1, threads1 >> > (m_pTmpLoop, m_pTraceRes, m_pTmpDeviceSum);

    _kernelInitialStaticPotentialCorrelatorOfSite << <block2, threads2 >> > (m_pCorrelatorCounter, m_pCorrelator);

    //========= It is impossible to calculate every pairs of sites
    //_kernelStaticPotentialCorrelatorOfSite << <block1, threads1 >> > (
    //    m_pTraceRes,
    //    CCommonData::m_sCenter,
    //    m_uiMaxLengthSq,
    //    m_pCorrelatorCounter,
    //    m_pCorrelator
    //);
    
    dim3 block3(_HC_Lx * _HC_Ly, _HC_Lx * _HC_Ly, 1);
    dim3 threads3(_HC_Lz, _HC_Lz / 2, 1);
    _kernelStaticPotentialCorrelatorOfSite2 << <block3, threads3 >> > (
        m_pTraceRes,
        m_uiMaxLengthSq,
        m_pCorrelatorCounter,
        m_pCorrelator
        );

    _kernelAverageStaticPotentialCorrelatorOfSite << <block2, threads2 >> > (m_pCorrelatorCounter, m_pCorrelator);


    //extract res
    
    checkCudaErrors(cudaMemcpy(res, m_pTmpDeviceSum, sizeof(CLGComplex), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostCorrelatorCounter, m_pCorrelatorCounter, sizeof(UINT) * m_uiMaxLengthSq, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostCorrelator, m_pCorrelator, sizeof(CLGComplex) * m_uiMaxLengthSq, cudaMemcpyDeviceToHost));

    const UINT uiSiteNumber = appGetLattice()->m_pIndexCache->m_uiSiteXYZ;
    res[0].x = res[0].x / uiSiteNumber;
    res[0].y = res[0].y / uiSiteNumber;
    m_lstAverageLoop.AddItem(res[0]);

    if (m_bShowResult)
    {
        appSetLogDate(FALSE);
        appDetailed(_T("==================== Polyakov Loop ============================\n"));
        appDetailed(_T("=== <P> = %f + %f I\n"), res[0].x, res[0].y);
    }

    if (0 == m_uiConfigurationCount)
    {
        assert(0 == m_lstR.Num());
        assert(0 == m_lstC.Num());

        for (UINT uiL = 1; uiL < m_uiMaxLengthSq; ++uiL)
        {
            if (m_pHostCorrelatorCounter[uiL] > 0)
            {
                m_lstR.AddItem(uiL);
                m_lstC.AddItem(m_pHostCorrelator[uiL]);

                if (m_bShowResult)
                {
                    appDetailed(_T("C(%f)=%f + %f I\n"), 
                        _hostsqrt(static_cast<Real>(uiL)),
                        m_pHostCorrelator[uiL].x,
                        m_pHostCorrelator[uiL].y);
                }
            }
        }
    }
    else
    {
        for (INT i = 0; i < m_lstR.Num(); ++i)
        {
            assert(m_pHostCorrelatorCounter[m_lstR[i]] > 0);
            m_lstC.AddItem(m_pHostCorrelator[m_lstR[i]]);

            if (m_bShowResult)
            {
                appDetailed(_T("C(%f)=%f + %f I\n"),
                    _hostsqrt(static_cast<Real>(m_lstR[i])),
                    m_pHostCorrelator[m_lstR[i]].x,
                    m_pHostCorrelator[m_lstR[i]].y);
            }
        }
    }

    if (m_bShowResult)
    {
        appSetLogDate(TRUE);
    }
    ++m_uiConfigurationCount;
}

void CMeasurePolyakov::Average(UINT )
{
    //nothing to do
}

void CMeasurePolyakov::Report()
{
    assert(m_uiConfigurationCount == static_cast<UINT>(m_lstAverageLoop.Num()));
    assert(static_cast<UINT>(m_uiConfigurationCount * m_lstR.Num())
        == static_cast<UINT>(m_lstC.Num()));

    appSetLogDate(FALSE);
    CLGComplex tmpChargeSum = _make_cuComplex(F(0.0), F(0.0));


    TArray<CLGComplex> tmpLoop;
    TArray<CLGComplex> tmpCorrelator;

    appGeneral(_T("\n\n==========================================================================\n"));
    appGeneral(_T("==================== Polyakov Loop (%d con)============================\n"), m_uiConfigurationCount);

    appGeneral(_T("\n ----------- Loop ------------- \n"));

    appGeneral(_T("{"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        tmpChargeSum.x += m_lstAverageLoop[i].x;
        tmpChargeSum.y += m_lstAverageLoop[i].y;
        LogGeneralComplex(m_lstAverageLoop[i]);
    }
    appGeneral(_T("}\n"));

    tmpChargeSum.x = tmpChargeSum.x / m_uiConfigurationCount;
    tmpChargeSum.y = tmpChargeSum.y / m_uiConfigurationCount;
    m_cAverageLoop = tmpChargeSum;
    appGeneral(_T("\n ----------- average Loop |<P>| = %2.12f arg(P) = %2.12f ------------- \n"), _cuCabsf(tmpChargeSum), __cuCargf(tmpChargeSum));

    m_lstAverageC.RemoveAll();

    appGeneral(_T("\n ----------- Static Quark Correlator ------------- \n"));

    appGeneral(_T("lengthR = {"));
    for (INT i = 0; i < m_lstR.Num(); ++i)
    {
        appGeneral(_T("%2.12f%s "), 
            _hostsqrt(static_cast<Real>(m_lstR[i])),
            (i == m_lstR.Num() - 1) ? _T("") : _T(","));
    }
    appGeneral(_T("}\n"));

    if (m_bShowResult)
    {
        appGeneral(_T("correlator = {\n"));
    }
    
    for (INT k = 0; k < m_lstR.Num(); ++k)
    {
        if (m_bShowResult)
        {
            appGeneral(_T("{"));
        }
        CLGComplex averageOfC_R = _make_cuComplex(F(0.0), F(0.0));
        for (UINT i = 0; i < m_uiConfigurationCount; ++i)
        {
            if (m_bShowResult)
            {
                LogGeneralComplex(m_lstC[i * m_lstR.Num() + k]);
            }
            averageOfC_R = _cuCaddf(averageOfC_R, m_lstC[i * m_lstR.Num() + k]);
        }
        if (m_bShowResult)
        {
            appGeneral(_T("},\n"));
        }

        averageOfC_R.x = averageOfC_R.x / m_uiConfigurationCount;
        averageOfC_R.y = averageOfC_R.y / m_uiConfigurationCount;
        m_lstAverageC.AddItem(averageOfC_R);
    }
    if (m_bShowResult)
    {
        appGeneral(_T("}\n"));
    }

    appGeneral(_T("averagecorrelator = {\n"));
    for (INT k = 0; k < m_lstR.Num(); ++k)
    {
        LogGeneralComplex(m_lstAverageC[k]);
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n==========================================================================\n"));
    appGeneral(_T("==========================================================================\n\n"));
    appSetLogDate(TRUE);
}

void CMeasurePolyakov::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstR.RemoveAll();
    m_lstC.RemoveAll();
    m_lstAverageLoop.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================