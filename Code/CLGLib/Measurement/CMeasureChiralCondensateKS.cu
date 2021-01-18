//=============================================================================
// FILENAME : CMeasureChiralCondensateKS.cpp
// 
// DESCRIPTION:
// almost copy from CMeasureChiralCondensate.cpp, but with Wilson SU3 vector to SU3 vector
//
// REVISION:
//  [10/01/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureChiralCondensateKS)

#pragma region kernels

/**
 * 
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelDotMeasureAllKS(
    const deviceSU3Vector* __restrict__ pZ4,
    const deviceSU3Vector* __restrict__ pApplied,
    CLGComplex* resultXYPlan,
#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex* result
#else
    CLGComplex* result
#endif
)
{
    intokernalInt4;

#if !_CLG_DOUBLEFLOAT
    result[uiSiteIndex] = _cToDouble(pZ4[uiSiteIndex].ConjugateDotC(pApplied[uiSiteIndex]));
    atomicAdd(&resultXYPlan[_ixy].x, static_cast<Real>(result[uiSiteIndex].x));
    atomicAdd(&resultXYPlan[_ixy].y, static_cast<Real>(result[uiSiteIndex].y));
#else
    result[uiSiteIndex] = pZ4[uiSiteIndex].ConjugateDotC(pApplied[uiSiteIndex]);
    atomicAdd(&resultXYPlan[_ixy].x, result[uiSiteIndex].x);
    atomicAdd(&resultXYPlan[_ixy].y, result[uiSiteIndex].y);
#endif
}

__global__ void _CLG_LAUNCH_BOUND
_kernelChiralCondensateInitialDistRKS(UINT* pCount)
{
    pCount[threadIdx.x] = 0;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelChiralCondensateInitialDistCondKS(CLGComplex* pCond)
{
    pCond[threadIdx.x] = _zeroc;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelChiralCondensateMeasureDistKS(
    const CLGComplex* __restrict__ CondXY,
    SSmallInt4 sCenter, UINT uiMax, BYTE byFieldId, UBOOL bCalcR,
    UINT* counter, 
    CLGComplex* CondR
)
{
    UINT uiXY = (threadIdx.x + blockIdx.x * blockDim.x);
    INT uiX = static_cast<INT>(uiXY / _DC_Ly);
    INT uiY = static_cast<INT>(uiXY % _DC_Ly);
    UINT uiC = (sCenter.x - uiX) * (sCenter.x - uiX)
        + (sCenter.y - uiY) * (sCenter.y - uiY);

    SSmallInt4 sSite4;
    sSite4.z = sCenter.z;
    sSite4.w = sCenter.w;
    sSite4.x = static_cast<SBYTE>(uiX);
    sSite4.y = static_cast<SBYTE>(uiY);
    if (uiC <= uiMax && !__idx->_deviceGetMappingIndex(sSite4, byFieldId).IsDirichlet())
    {
        if (bCalcR)
        {
            atomicAdd(&counter[uiC], 1);
        }
        
        atomicAdd(&CondR[uiC].x, CondXY[uiXY].x);
        atomicAdd(&CondR[uiC].y, CondXY[uiXY].y);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelChiralAverageDistKS(UINT* pCount, CLGComplex* pCond)
{
    const UINT uiIdx = threadIdx.x;
    if (pCount[uiIdx] > 0)
    {
        pCond[uiIdx].x = pCond[uiIdx].x / static_cast<Real>(pCount[uiIdx]);
        pCond[uiIdx].y = pCond[uiIdx].y / static_cast<Real>(pCount[uiIdx]);
    }
}

#pragma endregion

void CMeasureChiralCondensateKS::KSTraceEndZ4(
    UINT uiMaxR,
    BYTE byFieldId,
    UINT uiFieldCount,
    UINT uiMeasureCount,
    UINT uiCurrentMeasureCount,
    const CLGComplex* const* pXYBuffers,
    UINT* pCountBuffer,
    CLGComplex* pValueBuffer,
    UINT* pHostCountBuffer,
    CLGComplex* pHostValueBuffer,
    TArray<UINT>& lstR,
    TArray<CLGComplex>* lstValues,
    UBOOL bLogResult)
{
    dim3 block2(_HC_DecompX, 1, 1);
    dim3 threads2(_HC_DecompLx, 1, 1);
    dim3 block3(1, 1, 1);
    dim3 threads3(uiMaxR + 1, 1, 1);

    const Real fDivider = F(1.0) / (uiFieldCount * _HC_Lz * _HC_Lt);
    _kernelChiralCondensateInitialDistRKS << <block3, threads3 >> > (pCountBuffer);
    for (UINT i = 0; i < uiMeasureCount; ++i)
    {
        _kernelChiralCondensateInitialDistCondKS << <block3, threads3 >> > (pValueBuffer);

        _kernelChiralCondensateMeasureDistKS << <block2, threads2 >> > (
            pXYBuffers[i],
            CCommonData::m_sCenter,
            uiMaxR,
            byFieldId,
            0 == i,
            pCountBuffer,
            pValueBuffer
            );

        _kernelChiralAverageDistKS << <block3, threads3 >> > (pCountBuffer, pValueBuffer);

        if (0 == i)
        {
            checkCudaErrors(cudaMemcpy(pHostCountBuffer, pCountBuffer, sizeof(UINT) * (uiMaxR + 1), cudaMemcpyDeviceToHost));
        }

        checkCudaErrors(cudaMemcpy(pHostValueBuffer, pValueBuffer, sizeof(CLGComplex) * (uiMaxR + 1), cudaMemcpyDeviceToHost));

        if (0 == uiCurrentMeasureCount)
        {
            if (0 == i)
            {
                assert(0 == lstR.Num());
            }
            assert(0 == lstValues[i].Num());
            for (UINT uiL = 0; uiL <= uiMaxR; ++uiL)
            {
                if (pHostCountBuffer[uiL] > 0)
                {
                    if (0 == i)
                    {
                        lstR.AddItem(uiL);
                    }
                    lstValues[i].AddItem(cuCmulf_cr(pHostValueBuffer[uiL], fDivider));

                    if (bLogResult)
                    {
                        appDetailed(_T("Cond %d (r = %f)= %f + %f i\n"),
                            i,
                            _hostsqrt(static_cast<Real>(uiL)),
                            pHostValueBuffer[uiL].x,
                            pHostValueBuffer[uiL].y
                        );
                    }
                }
            }
        }
        else
        {
            for (INT j = 0; j < lstR.Num(); ++j)
            {
                assert(pHostCountBuffer[lstR[j]] > 0);
                lstValues[i].AddItem(cuCmulf_cr(pHostValueBuffer[lstR[j]], fDivider));

                if (bLogResult)
                {
                    appDetailed(_T("Cond %d (r = %f)=%f + %f i\n"),
                        i,
                        _hostsqrt(static_cast<Real>(lstR[j])),
                        pHostValueBuffer[lstR[j]].x,
                        pHostValueBuffer[lstR[j]].y
                    );
                }
            }
        }
    }
}

CMeasureChiralCondensateKS::~CMeasureChiralCondensateKS()
{
    if (NULL != m_pDeviceXYBuffer[0])
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            checkCudaErrors(cudaFree(m_pDeviceXYBuffer[i]));
        }
    }

    if (NULL != m_pHostXYBuffer)
    {
        free(m_pHostXYBuffer);
    }

    if (NULL != m_pDistributionR)
    {
        checkCudaErrors(cudaFree(m_pDistributionR));
    }

    if (NULL != m_pDistribution)
    {
        checkCudaErrors(cudaFree(m_pDistribution));
    }

    if (NULL != m_pHostDistributionR)
    {
        free(m_pHostDistributionR);
    }

    if (NULL != m_pHostDistribution)
    {
        free(m_pHostDistribution);
    }
}

void CMeasureChiralCondensateKS::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasureStochastic::Initial(pOwner, pLatticeData, param, byId);

    for (UINT i = 0; i < ChiralKSMax; ++i)
    {
        checkCudaErrors(cudaMalloc((void**)&m_pDeviceXYBuffer[i], sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    }    
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
        checkCudaErrors(cudaMalloc((void**)&m_pDistribution, sizeof(CLGComplex) * (m_uiMaxR + 1)));

        m_pHostDistributionR = (UINT*)malloc(sizeof(UINT) * (m_uiMaxR + 1));
        m_pHostDistribution = (CLGComplex*)malloc(sizeof(CLGComplex) * (m_uiMaxR + 1));
    }
}

void CMeasureChiralCondensateKS::OnConfigurationAcceptedZ4(
    const class CFieldGauge* pAcceptGauge, 
    const class CFieldGauge* pCorrespondingStaple, 
    const class CFieldFermion* pZ4, 
    const class CFieldFermion* pInverseZ4, 
    UBOOL bStart, 
    UBOOL bEnd)
{
    if (bStart)
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            _ZeroXYPlaneC(m_pDeviceXYBuffer[i]);
            m_cTmpSum[i] = _zeroc;
        }
    }

    const Real oneOuiVolume = F(1.0) / appGetLattice()->m_pIndexCache->m_uiSiteNumber[m_byFieldId];
    const CFieldFermionKSSU3 * pF1W = dynamic_cast<const CFieldFermionKSSU3*>(pZ4);
    const CFieldFermionKSSU3* pF2W = dynamic_cast<const CFieldFermionKSSU3*>(pInverseZ4);
    CFieldFermionKSSU3* pAfterApplied = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));

#pragma region Dot

    // The results are Atomic Add to m_pDeviceXYBuffer
    preparethread;
    for (BYTE i = 0; i < ChiralKSMax; ++i)
    {
        switch ((EChiralMeasureTypeKS)i)
        {
        case ChiralKS:
            {
                pF2W->CopyTo(pAfterApplied);
            }
            break;
        case ConnectSusp:
            {
                pF2W->CopyTo(pAfterApplied);
                pAfterApplied->InverseD(pAcceptGauge);
            }
            break;
        }

        _kernelDotMeasureAllKS << <block, threads >> > (
            pF1W->m_pDeviceData,
            pF2W->m_pDeviceData,
            m_pDeviceXYBuffer[i],
            _D_ComplexThreadBuffer
            );

#if !_CLG_DOUBLEFLOAT
        const CLGComplex thisSum = _cToFloat(appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer));
#else
        const CLGComplex thisSum = appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
#endif
        m_cTmpSum[i] = _cuCaddf(m_cTmpSum[i], cuCmulf_cr(thisSum, oneOuiVolume));
    }
    pAfterApplied->Return();

#pragma endregion

    if (bEnd)
    {
        if (m_bMeasureDistribution)
        {
            KSTraceEndZ4(
                m_uiMaxR,
                m_byFieldId,
                m_uiFieldCount,
                ChiralKSMax,
                m_uiConfigurationCount,
                m_pDeviceXYBuffer,
                m_pDistributionR,
                m_pDistribution,
                m_pHostDistributionR,
                m_pHostDistribution,
                m_lstR,
                m_lstCond,
                m_bShowResult);
        }

        const Real fDiv2 = F(1.0) / m_uiFieldCount;
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            m_cTmpSum[i] = cuCmulf_cr(m_cTmpSum[i], fDiv2);
            appDetailed(_T("\n Condensate %d = %2.12f + %2.12f\n"), i, m_cTmpSum[i].x, m_cTmpSum[i].y);
            m_lstCondAll[i].AddItem(m_cTmpSum[i]);
        }

        ++m_uiConfigurationCount;
    }
}

void CMeasureChiralCondensateKS::OnConfigurationAccepted(const CFieldGauge* pGauge, const CFieldGauge* pCorrespondingStaple)
{

}

void CMeasureChiralCondensateKS::Average(UINT )
{
    //nothing to do
}

void CMeasureChiralCondensateKS::Report()
{
    for (UINT i = 0; i < ChiralKSMax; ++i)
    {
        assert(m_uiConfigurationCount == static_cast<UINT>(m_lstCondAll[i].Num()));

        appGeneral(_T("\n==========================================================================\n"));
        appGeneral(_T("==================== Condensate No %d (%d con)============================\n"), i, m_uiConfigurationCount);
        CLGComplex tmpChargeSum = _zeroc;
        if (m_uiConfigurationCount > 1)
        {
            appGeneral(_T("\n ----------- each configuration ------------- \n"));
            appGeneral(_T("{"));

            for (UINT j = 0; j < m_uiConfigurationCount; ++j)
            {
                tmpChargeSum.x += m_lstCondAll[i][j].x;
                tmpChargeSum.y += m_lstCondAll[i][j].y;
                LogGeneralComplex(m_lstCondAll[i][j]);
            }
            appGeneral(_T("}\n"));

            tmpChargeSum.x = tmpChargeSum.x / m_uiConfigurationCount;
            tmpChargeSum.y = tmpChargeSum.y / m_uiConfigurationCount;
            appGeneral(_T("\n ----------- average condensate = %2.12f + %2.12f ------------- \n"),
                tmpChargeSum.x, tmpChargeSum.y);

            //m_cAverageCondensate = tmpChargeSum;
        }
        else
        {
            appGeneral(_T("\n ----------- average condensate = %2.12f + %2.12f ------------- \n"),
                m_lstCondAll[i][0].x,
                m_lstCondAll[i][0].y);

            //m_cAverageCondensate = m_lstCondAll[i][0];
        }
    }

    appGeneral(_T("==========================================================================\n"));
    appSetLogDate(TRUE);
}

void CMeasureChiralCondensateKS::Reset()
{
    m_uiConfigurationCount = 0;
    for (UINT i = 0; i < ChiralKSMax; ++i)
    {
        m_lstCondAll[i].RemoveAll();
        m_lstCond[i].RemoveAll();
    }
    m_lstR.RemoveAll();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================