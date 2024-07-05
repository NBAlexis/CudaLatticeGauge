//=============================================================================
// FILENAME : CMeasurePolyakovSUN.cu
// 
// DESCRIPTION:
//
//
// REVISION:
//  [07/07/2019 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "CMeasurePolyakovSUN.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasurePolyakov2)

#pragma region kernles 


__global__ void
_CLG_LAUNCH_BOUND
_kernelInitialStaticPotentialCorrelatorOfSiteSUN(UINT* counter, CLGComplex* correlator)
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
__global__ void _CLG_LAUNCH_BOUND
_kernelStaticPotentialCorrelatorOfSiteSUN(
    const cuDoubleComplex* __restrict__ traceXYZ,
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
        CLGComplex correlatorres = _cuCmulf(_cToRealC(traceXYZ[uiXYZ]), _cuConjf(_cToRealC(traceXYZ[uiCenter])));
        atomicAdd(&counter[uiC], 1);
        atomicAdd(&correlator[uiC].x, correlatorres.x);
        atomicAdd(&correlator[uiC].y, correlatorres.y);
    }
}

/**
 * Block.x = x * Ly + y, Thread.x = Lz
 * Block.y = shift xy, Thread.z = shift z
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelStaticPotentialCorrelatorOfSite2SUN(
    const cuDoubleComplex* __restrict__ traceXYZ, 
    UINT uiMax,
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
        CLGComplex correlatorres = _cuCmulf(_cToRealC(traceXYZ[uiXYZ1]), _cuConjf(_cToRealC(traceXYZ[uiXYZ2])));
        atomicAdd(&counter[uiC], 1);
        atomicAdd(&correlator[uiC].x, correlatorres.x);
        atomicAdd(&correlator[uiC].y, correlatorres.y);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelAverageStaticPotentialCorrelatorOfSiteSUN(UINT* counter, CLGComplex* correlator)
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

CMeasurePolyakov2::~CMeasurePolyakov2()
{
    if (NULL != m_pHostCorrelator)
    {
        free(m_pHostCorrelator);
    }

    if (NULL != m_pHostCorrelatorCounter)
    {
        free(m_pHostCorrelatorCounter);
    }

    if (NULL != m_pTraceRes)
    {
        checkCudaErrors(cudaFree(m_pTraceRes));
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

void CMeasurePolyakov2::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    checkCudaErrors(cudaMalloc((void**)&m_pTraceRes, sizeof(cuDoubleComplex) * _HC_Lx * _HC_Ly * _HC_Lz));

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
}

void CMeasurePolyakov2::OnConfigurationAcceptedSingleField(const class CFieldGauge* pAcceptGauge, const class CFieldGauge* pCorrespondingStaple)
{
    if (NULL == pAcceptGauge)
    {
        appCrucial(_T("CMeasureMesonCorrelator only implemented with gauge SUN!\n"));
        return;
    }

    pAcceptGauge->PolyakovOnSpatialSite(m_pTraceRes);

    dim3 block2(1, 1, 1);
    dim3 threads2(m_uiMaxLengthSq, 1, 1);
    _kernelInitialStaticPotentialCorrelatorOfSiteSUN << <block2, threads2 >> > (m_pCorrelatorCounter, m_pCorrelator);
    
    dim3 block3(_HC_Lx * _HC_Ly, _HC_Lx * _HC_Ly, 1);
    dim3 threads3(_HC_Lz, _HC_Lz / 2, 1);
    _kernelStaticPotentialCorrelatorOfSite2SUN << <block3, threads3 >> > (
        m_pTraceRes,
        m_uiMaxLengthSq,
        m_pCorrelatorCounter,
        m_pCorrelator
        );

    _kernelAverageStaticPotentialCorrelatorOfSiteSUN << <block2, threads2 >> > (m_pCorrelatorCounter, m_pCorrelator);


    //extract res
    checkCudaErrors(cudaMemcpy(m_pHostCorrelatorCounter, m_pCorrelatorCounter, sizeof(UINT) * m_uiMaxLengthSq, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(m_pHostCorrelator, m_pCorrelator, sizeof(CLGComplex) * m_uiMaxLengthSq, cudaMemcpyDeviceToHost));

    CLGComplex res = _cToRealC(appGetCudaHelper()->ReduceComplex(m_pTraceRes, _HC_Lx * _HC_Ly * _HC_Lz));
    const UINT uiSiteNumber = appGetLattice()->m_pIndexCache->m_uiSiteXYZ;
    res.x = res.x / uiSiteNumber;
    res.y = res.y / uiSiteNumber;
    UpdateComplexResult(res, FALSE);

    if (m_bShowResult)
    {
        appPushLogDate(FALSE);
        appDetailed(_T("==================== Polyakov Loop ============================\n"));
        appDetailed(_T("=== <P> = %f + %f I\n"), res.x, res.y);
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
        appPopLogDate();
    }
    ++m_uiConfigurationCount;
}

void CMeasurePolyakov2::Report()
{
    Average();
    assert(static_cast<UINT>(m_uiConfigurationCount * m_lstR.Num())
        == static_cast<UINT>(m_lstC.Num()));

    appPushLogDate(FALSE);

    TArray<CLGComplex> tmpLoop;
    TArray<CLGComplex> tmpCorrelator;

    appGeneral(_T("\n\n==========================================================================\n"));
    appGeneral(_T("==================== Polyakov Loop (%d con)============================\n"), m_uiConfigurationCount);

    appGeneral(_T("\n ----------- Loop ------------- \n"));

    appGeneral(_T("{"));
    for (UINT i = 0; i < m_uiConfigurationCount; ++i)
    {
        LogGeneralComplex(CmpResAtI(i));
    }
    appGeneral(_T("}\n"));

    appGeneral(_T("\n ----------- average Loop |<P>| = %2.12f arg(P) = %2.12f ------------- \n"), _cuCabsf(GetAverageCmpRes()), __cuCargf(GetAverageCmpRes()));

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
    appPopLogDate();
}

void CMeasurePolyakov2::Reset()
{
    CMeasure::Reset();
    m_lstR.RemoveAll();
    m_lstC.RemoveAll();
}


__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================