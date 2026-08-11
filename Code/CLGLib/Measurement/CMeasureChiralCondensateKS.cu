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
#include "Data/Field/Staggered/CFieldFermionKST.h"
#include "CMeasureChiralCondensateKS.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CMeasureChiralCondensateKS)

#pragma region kernels

/**
 * 
 */
template<class vectorData>
__global__ void _CLG_LAUNCH_BOUND
_kernelDotMeasureAllKS(
    const vectorData* __restrict__ pZ4,
    const vectorData* __restrict__ pApplied,
    CLGComplex* resultXYPlan,
    CLGComplex* resX,
    CLGComplex* resY,
    CLGComplex* resZ,
    CLGComplex* resT,
    cuDoubleComplex* result
)
{
    intokernalInt4;
    const UINT _ixy = uiSiteIndex / _DC_GridDimZT;

#if !_CLG_DOUBLEFLOAT
    result[uiSiteIndex] = _cToDouble(_dot(pZ4[uiSiteIndex], pApplied[uiSiteIndex]));
    atomicAdd(&resultXYPlan[_ixy].x, static_cast<Real>(result[uiSiteIndex].x));
    atomicAdd(&resultXYPlan[_ixy].y, static_cast<Real>(result[uiSiteIndex].y));
    if (NULL != resX)
    {
        atomicAdd(&resX[sSite4.x].x, static_cast<Real>(result[uiSiteIndex].x));
        atomicAdd(&resX[sSite4.x].y, static_cast<Real>(result[uiSiteIndex].y));
    }
    if (NULL != resY)
    {
        atomicAdd(&resY[sSite4.y].x, static_cast<Real>(result[uiSiteIndex].x));
        atomicAdd(&resY[sSite4.y].y, static_cast<Real>(result[uiSiteIndex].y));
    }
    if (NULL != resZ)
    {
        atomicAdd(&resZ[sSite4.z].x, static_cast<Real>(result[uiSiteIndex].x));
        atomicAdd(&resZ[sSite4.z].y, static_cast<Real>(result[uiSiteIndex].y));
}
    if (NULL != resT)
    {
        atomicAdd(&resT[sSite4.w].x, static_cast<Real>(result[uiSiteIndex].x));
        atomicAdd(&resT[sSite4.w].y, static_cast<Real>(result[uiSiteIndex].y));
    }
#else
    result[uiSiteIndex] = _dot(pZ4[uiSiteIndex], pApplied[uiSiteIndex]);
    atomicAdd(&resultXYPlan[_ixy].x, result[uiSiteIndex].x);
    atomicAdd(&resultXYPlan[_ixy].y, result[uiSiteIndex].y);
    if (NULL != resX)
    {
        atomicAdd(&resX[sSite4.x].x, result[uiSiteIndex].x);
        atomicAdd(&resX[sSite4.x].y, result[uiSiteIndex].y);
    }
    if (NULL != resY)
    {
        atomicAdd(&resY[sSite4.y].x, result[uiSiteIndex].x);
        atomicAdd(&resY[sSite4.y].y, result[uiSiteIndex].y);
    }
    if (NULL != resZ)
    {
        atomicAdd(&resZ[sSite4.z].x, result[uiSiteIndex].x);
        atomicAdd(&resZ[sSite4.z].y, result[uiSiteIndex].y);
    }
    if (NULL != resT)
    {
        atomicAdd(&resT[sSite4.w].x, result[uiSiteIndex].x);
        atomicAdd(&resT[sSite4.w].y, result[uiSiteIndex].y);
    }
#endif
}

//__global__ void
//_CLG_LAUNCH_BOUND
//_kernelFillZSlice(
//    const CLGComplex* __restrict__ res,
//    CLGComplex** resZ)
//{
//    UINT uiXY= (threadIdx.x + blockIdx.x * blockDim.x);
//    UINT uiT = (threadIdx.z + blockIdx.z * blockDim.z);
//    UINT uiZ = threadIdx.y + blockIdx.y * blockDim.y;
//    resZ[uiZ][uiXY * _DC_Lt + uiT]
//    = res[uiXY * _DC_GridDimZT + uiZ * _DC_Lt + uiT];
//}

//__global__ void
//_CLG_LAUNCH_BOUND
//_kernelInitialZSliceChiralKS(CLGComplex* resZ, UINT uiMax)
//{
//    const UINT idx = threadIdx.x + blockIdx.x * blockDim.x;
//    if (idx < uiMax)
//    {
//        resZ[idx] = _zeroc;
//    }
//}

#pragma endregion


CMeasureChiralCondensateKS::~CMeasureChiralCondensateKS()
{
    if (NULL != m_pDeviceXYBuffer[0])
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            checkCudaErrors(__cudaFree(m_pDeviceXYBuffer[i]));
        }
    }

    if (NULL != m_pDeviceXBuffer[0])
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            checkCudaErrors(__cudaFree(m_pDeviceXBuffer[i]));
        }
    }
    if (NULL != m_pDeviceYBuffer[0])
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            checkCudaErrors(__cudaFree(m_pDeviceYBuffer[i]));
        }
    }
    if (NULL != m_pDeviceZBuffer[0])
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            checkCudaErrors(__cudaFree(m_pDeviceZBuffer[i]));
        }
    }
    if (NULL != m_pDeviceTBuffer[0])
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            checkCudaErrors(__cudaFree(m_pDeviceTBuffer[i]));
        }
    }

    if (NULL != m_pHostXYBuffer)
    {
        free(m_pHostXYBuffer);
    }

    if (NULL != m_pHostXBuffer)
    {
        free(m_pHostXBuffer);
    }
    if (NULL != m_pHostYBuffer)
    {
        free(m_pHostYBuffer);
    }
    if (NULL != m_pHostZBuffer)
    {
        free(m_pHostZBuffer);
    }
    if (NULL != m_pHostTBuffer)
    {
        free(m_pHostTBuffer);
    }

    if (NULL != m_pDistributionR)
    {
        checkCudaErrors(__cudaFree(m_pDistributionR));
    }

    if (NULL != m_pDistribution)
    {
        checkCudaErrors(__cudaFree(m_pDistribution));
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
        checkCudaErrors(__cudaMalloc((void**)&m_pDeviceXYBuffer[i], sizeof(CLGComplex) * _HC_Lx * _HC_Ly));
    }    
    m_pHostXYBuffer = (CLGComplex*)malloc(sizeof(CLGComplex) * _HC_Lx * _HC_Ly);

    Reset();

    INT iValue = 0;
    param.FetchValueINT(_T("ShiftCenter"), iValue);
    m_bShiftCenter = iValue != 0;

    iValue = 0;
    param.FetchValueINT(_T("MeasureSigma12"), iValue);
    m_bMeasureSigma12 = iValue != 0;

    iValue = 0;
    param.FetchValueINT(_T("MeasureConnect"), iValue);
    m_bMeasureConnect = iValue != 0;

    iValue = 0;
    param.FetchValueINT(_T("XSlice"), iValue);
    m_bMeasureXSlice = iValue != 0;
    iValue = 0;
    param.FetchValueINT(_T("YSlice"), iValue);
    m_bMeasureYSlice = iValue != 0;
    iValue = 0;
    param.FetchValueINT(_T("ZSlice"), iValue);
    m_bMeasureZSlice = iValue != 0;
    iValue = 0;
    param.FetchValueINT(_T("TSlice"), iValue);
    m_bMeasureTSlice = iValue != 0;

    if (m_bMeasureXSlice)
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            checkCudaErrors(__cudaMalloc((void**)&m_pDeviceXBuffer[i], sizeof(CLGComplex) * _HC_Lx));
        }
        m_pHostXBuffer = (CLGComplex*)malloc(sizeof(CLGComplex) * _HC_Lx);
    }
    else
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            m_pDeviceXBuffer[i] = NULL;
        }
    }
    if (m_bMeasureYSlice)
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            checkCudaErrors(__cudaMalloc((void**)&m_pDeviceYBuffer[i], sizeof(CLGComplex) * _HC_Ly));
        }
        m_pHostYBuffer = (CLGComplex*)malloc(sizeof(CLGComplex) * _HC_Ly);
    }
    else
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            m_pDeviceZBuffer[i] = NULL;
        }
    }
    if (m_bMeasureZSlice)
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            checkCudaErrors(__cudaMalloc((void**)&m_pDeviceZBuffer[i], sizeof(CLGComplex) * _HC_Lz));
        }
        m_pHostZBuffer = (CLGComplex*)malloc(sizeof(CLGComplex) * _HC_Lz);
    }
    else
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            m_pDeviceZBuffer[i] = NULL;
        }
    }
    if (m_bMeasureTSlice)
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            checkCudaErrors(__cudaMalloc((void**)&m_pDeviceTBuffer[i], sizeof(CLGComplex) * _HC_Lt));
        }
        m_pHostTBuffer = (CLGComplex*)malloc(sizeof(CLGComplex) * _HC_Lt);
    }
    else
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            m_pDeviceTBuffer[i] = NULL;
        }
    }

    //assuming the center is really at center
    SetMaxAndEdge(&m_uiMaxR, &m_uiEdge, m_bShiftCenter);

    checkCudaErrors(__cudaMalloc((void**)&m_pDistributionR, sizeof(UINT) * (m_uiMaxR + 1)));
    checkCudaErrors(__cudaMalloc((void**)&m_pDistribution, sizeof(CLGComplex) * (m_uiMaxR + 1)));

    m_pHostDistributionR = (UINT*)malloc(sizeof(UINT) * (m_uiMaxR + 1));
    m_pHostDistribution = (CLGComplex*)malloc(sizeof(CLGComplex) * (m_uiMaxR + 1));
}

void CMeasureChiralCondensateKS::OnConfigurationAcceptedZ4(
    INT gaugeNum,
    INT bosonNum,
    INT tensor2Num,
    const class CFieldGauge* const* pAcceptGauge, 
    const class CFieldBoson* const* pAcceptBoson,
    const class CFieldTensor2* const* tensor2Fields,
    const class CFieldGauge* const* pCorrespondingStaple, 
    const class CFieldFermion* pZ4, 
    const class CFieldFermion* pInverseZ4, 
    UBOOL bStart, 
    UBOOL bEnd)
{
    if (bStart)
    {
        for (UINT i = 0; i < ChiralKSMax; ++i)
        {
            _ZeroXYPlane(m_pDeviceXYBuffer[i]);
            //m_cTmpSum[i] = _zeroc;
            if (m_bDebugDivation)
            {
                m_lstDebugData[i].RemoveAll();
            }

            if (m_bMeasureXSlice)
            {
                _ZeroSlice(m_pDeviceXBuffer[i], 0);
            }

            if (m_bMeasureYSlice)
            {
                _ZeroSlice(m_pDeviceYBuffer[i], 1);
            }

            if (m_bMeasureZSlice)
            {
                _ZeroSlice(m_pDeviceZBuffer[i], 2);
            }

            if (m_bMeasureTSlice)
            {
                _ZeroSlice(m_pDeviceTBuffer[i], 3);
            }
        }
    }

    // oneOuiVolume is only used in debug deviation, "-1" is for <qbar M q> = -tr[MD^{-1}]
    //P4-3.4: global volume for the debug per-site normalization.
    DOUBLE fCondVolume = static_cast<DOUBLE>(appGetLattice()->m_pIndexCache->m_uiSiteNumber[GetFermionFieldId()]);
#if _CLG_MULTI_GPU
    if (NULL != appGetComm())
    {
        const UINT* pG = appGetComm()->GlobalLattice();
        fCondVolume = static_cast<DOUBLE>(pG[0]) * pG[1] * pG[2] * pG[3];
    }
#endif
    const Real oneOuiVolume = F(-1.0) / fCondVolume;
    const CFieldFermionKS * pF1W = dynamic_cast<const CFieldFermionKS*>(pZ4);
    const CFieldFermionKS* pF2W = dynamic_cast<const CFieldFermionKS*>(pInverseZ4);
    CFieldFermionKS* pAfterApplied = dynamic_cast<CFieldFermionKS*>(appGetLattice()->GetPooledFieldById(GetFermionFieldId(), _T(__FILE__), __LINE__));

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
        case CMTKSGamma1:
        case CMTKSGamma2:
        case CMTKSGamma3:
        case CMTKSGamma4:
        case CMTKSGamma5:
        case CMTKSGamma51:
        case CMTKSGamma52:
        case CMTKSGamma53:
        case CMTKSGamma54:
            {
                pF2W->CopyTo(pAfterApplied);
                pAfterApplied->ApplyGamma(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson, (EGammaMatrix)(i - 1));
            }
            break;
        case CMTKSSigma12:
            {
                pF2W->CopyTo(pAfterApplied);
                if (m_bMeasureSigma12)
                {
                    pAfterApplied->ApplyGamma(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson, SIGMA12);
                }
            }
            break;
        case CMTKSSigma13:
            {
                pF2W->CopyTo(pAfterApplied);
                pAfterApplied->ApplyGamma(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson, SIGMA31);
            }
            break;
        case CMTKSSigma14:
            {
                pF2W->CopyTo(pAfterApplied);
                pAfterApplied->ApplyGamma(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson, SIGMA41);
            }
            break;
        case CMTKSSigma23:
            {
                pF2W->CopyTo(pAfterApplied);
                pAfterApplied->ApplyGamma(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson, SIGMA23);
            }
            break;
        case CMTKSSigma24:
            {
                pF2W->CopyTo(pAfterApplied);
                pAfterApplied->ApplyGamma(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson, SIGMA42);
            }
            break;
        case CMTKSSigma34:
            {
                pF2W->CopyTo(pAfterApplied);
                pAfterApplied->ApplyGamma(gaugeNum, bosonNum, pAcceptGauge, pAcceptBoson, SIGMA43);
            }
            break;
        case ConnectSusp:
            {
                pF2W->CopyTo(pAfterApplied);
                if (m_bMeasureConnect)
                {
                    pAfterApplied->InverseD(gaugeNum, bosonNum, tensor2Num, pAcceptGauge, pAcceptBoson, tensor2Fields);
                }
            }
            break;
        default:
            appCrucial(_T("This is impossible!\n"));
            break;
        }

        switch (pF1W->GetFieldType())
        {
        case EFT_FermionStaggeredSU3:
            {
                const CFieldFermionKSSU3* pF1WSU3 = dynamic_cast<const CFieldFermionKSSU3*>(pF1W);
                const CFieldFermionKSSU3* pAfterSU3 = dynamic_cast<const CFieldFermionKSSU3*>(pAfterApplied);
                _LAUNCH_KERNEL(_kernelDotMeasureAllKS<deviceSU3Vector>, block, threads, 
                    pF1WSU3->m_pDeviceData,
                    pAfterSU3->m_pDeviceData,
                    m_pDeviceXYBuffer[i],
                    m_bMeasureXSlice ? m_pDeviceXBuffer[i] : NULL,
                    m_bMeasureYSlice ? m_pDeviceYBuffer[i] : NULL,
                    m_bMeasureZSlice ? m_pDeviceZBuffer[i] : NULL,
                    m_bMeasureTSlice ? m_pDeviceTBuffer[i] : NULL,
                    _D_ComplexThreadBuffer
                    );
            }
            break;
        case EFT_FermionStaggeredU1:
            {
                const CFieldFermionKSU1* pF1WU1 = dynamic_cast<const CFieldFermionKSU1*>(pF1W);
                const CFieldFermionKSU1* pAfterU1 = dynamic_cast<const CFieldFermionKSU1*>(pAfterApplied);
                _LAUNCH_KERNEL(_kernelDotMeasureAllKS<CLGComplex>, block, threads,
                    pF1WU1->m_pDeviceData,
                    pAfterU1->m_pDeviceData,
                    m_pDeviceXYBuffer[i],
                    m_bMeasureXSlice ? m_pDeviceXBuffer[i] : NULL,
                    m_bMeasureYSlice ? m_pDeviceYBuffer[i] : NULL,
                    m_bMeasureZSlice ? m_pDeviceZBuffer[i] : NULL,
                    m_bMeasureTSlice ? m_pDeviceTBuffer[i] : NULL,
                    _D_ComplexThreadBuffer
                    );
            }
            break;
        default:
            {
                appCrucial(_T("CMeasureChiralCondensateKS unsupported field type!\n"));
            }
            break;
        }


#if !_CLG_DOUBLEFLOAT
        CLGComplex thisSum = _cToFloat(appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer));
#else
        CLGComplex thisSum = appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
#endif
        //P4-3.4: ThreadBufferSum is the LOCAL partial sum; reduce across ranks.
        GlobalSumComplex(thisSum);
        //m_cTmpSum[i] = _cuCaddf(m_cTmpSum[i], cuCmulf_cr(thisSum, oneOuiVolume));
        if (m_bDebugDivation)
        {
            m_lstDebugData[i].AddItem(cuCmulf_cr(thisSum, oneOuiVolume));
        }
    }
    pAfterApplied->Return();

#pragma endregion

    if (bEnd)
    {
        if (m_bDebugDivation)
        {
            appGeneral(_T("Debug data:\n"));
            for (BYTE i = 0; i < ChiralKSMax; ++i)
            {
                appGeneral(_T("{"));
                for (INT j = 0; j < m_lstDebugData[i].Num(); ++j)
                {
                    LogGeneralComplex(m_lstDebugData[i][j]);
                }
                appGeneral(_T("}\n"));
            }
        }

#if _CLG_MULTI_GPU
        //P4-3.4: the XY distribution accumulates per-rank partial sums over the
        //local z/t extent; sum across ranks so the R-distribution is global.
        //Requires x/y NOT split (same (x,y) index set per rank).
        if (NULL != appGetComm())
        {
            if (appGetComm()->GpuGrid()[0] > 1 || appGetComm()->GpuGrid()[1] > 1)
            {
                appCrucial(_T("CMeasureChiralCondensateKS: the XY distribution is not supported on multi-GPU with a split x/y direction. Rejected.\n"));
                return;
            }
            for (INT i = 0; i < static_cast<INT>(ChiralKSMax); ++i)
            {
                checkCudaErrors(cudaMemcpy(m_pHostXYBuffer, m_pDeviceXYBuffer[i], sizeof(CLGComplex) * _HC_Lx * _HC_Ly, cudaMemcpyDeviceToHost));
                GlobalSumComplexArray(m_pHostXYBuffer, _HC_Lx * _HC_Ly);
                checkCudaErrors(cudaMemcpy(m_pDeviceXYBuffer[i], m_pHostXYBuffer, sizeof(CLGComplex) * _HC_Lx * _HC_Ly, cudaMemcpyHostToDevice));
            }
        }
#endif

        TransformFromXYDataToRData(
            TRUE,
            m_bShiftCenter,
            m_uiMaxR,
            m_uiEdge,
            GetFermionFieldId(),
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
            m_lstCondAll,
            m_lstCondIn
        );

        if (m_bMeasureXSlice)
        {
            //P4-3.4: global normalization; per-rank slice partials summed.
            const Real fDemon = F(-1.0) / static_cast<Real> (m_uiFieldCount * GlobalL(1) * GlobalL(2) * GlobalL(3));
            for (INT i = 0; i < static_cast<INT>(ChiralKSMax); ++i)
            {
                checkCudaErrors(cudaMemcpy(m_pHostXBuffer, m_pDeviceXBuffer[i], sizeof(CLGComplex) * _HC_Lx, cudaMemcpyDeviceToHost));
                GlobalSumComplexArray(m_pHostXBuffer, _HC_Lx);
                for (UINT j = 0; j < _HC_Lx; ++j)
                {
                    m_lstCondXSlice[i].AddItem(cuCmulf_cr(m_pHostXBuffer[j], fDemon));
                }
            }
        }
        if (m_bMeasureYSlice)
        {
            //P4-3.4: global normalization; per-rank slice partials summed.
            const Real fDemon = F(-1.0) / static_cast<Real> (m_uiFieldCount * GlobalL(0) * GlobalL(2) * GlobalL(3));
            for (INT i = 0; i < static_cast<INT>(ChiralKSMax); ++i)
            {
                checkCudaErrors(cudaMemcpy(m_pHostYBuffer, m_pDeviceYBuffer[i], sizeof(CLGComplex) * _HC_Ly, cudaMemcpyDeviceToHost));
                GlobalSumComplexArray(m_pHostYBuffer, _HC_Ly);
                for (UINT j = 0; j < _HC_Ly; ++j)
                {
                    m_lstCondYSlice[i].AddItem(cuCmulf_cr(m_pHostYBuffer[j], fDemon));
                }
            }
        }
        if (m_bMeasureZSlice)
        {
#if _CLG_MULTI_GPU
            //P4-3.4: with a split z direction each rank holds a different
            //global-z slice; assembling the global profile is a gather, not a
            //sum. Reject explicitly instead of a silently-wrong partial profile.
            if (NULL != appGetComm() && appGetComm()->GpuGrid()[2] > 1)
            {
                appCrucial(_T("CMeasureChiralCondensateKS: m_bMeasureZSlice is not supported on multi-GPU with a split z direction. Rejected.\n"));
                return;
            }
#endif
            // "-1" comes from <qbar M q> = -tr[M D^{-1}]
            const Real fDemon = F(-1.0) / static_cast<Real> (m_uiFieldCount * GlobalL(0) * GlobalL(1) * GlobalL(3));
            for (INT i = 0; i < static_cast<INT>(ChiralKSMax); ++i)
            {
                checkCudaErrors(cudaMemcpy(m_pHostZBuffer, m_pDeviceZBuffer[i], sizeof(CLGComplex) * _HC_Lz, cudaMemcpyDeviceToHost));
                GlobalSumComplexArray(m_pHostZBuffer, _HC_Lz);
                for (UINT j = 0; j < _HC_Lz; ++j)
                {
                    m_lstCondZSlice[i].AddItem(cuCmulf_cr(m_pHostZBuffer[j], fDemon));
                }
            }
        }
        if (m_bMeasureTSlice)
        {
            //P4-3.4: global normalization; per-rank slice partials summed.
            const Real fDemon = F(-1.0) / static_cast<Real> (m_uiFieldCount * GlobalL(0) * GlobalL(1) * GlobalL(2));
            for (INT i = 0; i < static_cast<INT>(ChiralKSMax); ++i)
            {
                checkCudaErrors(cudaMemcpy(m_pHostTBuffer, m_pDeviceTBuffer[i], sizeof(CLGComplex) * _HC_Lt, cudaMemcpyDeviceToHost));
                GlobalSumComplexArray(m_pHostTBuffer, _HC_Lt);
                for (UINT j = 0; j < _HC_Lt; ++j)
                {
                    m_lstCondTSlice[i].AddItem(cuCmulf_cr(m_pHostTBuffer[j], fDemon));
                }
            }
        }
        UpdateRealResult(m_lstCondAll[0][m_uiConfigurationCount].x, FALSE);
        UpdateComplexResult(m_lstCondAll[0][m_uiConfigurationCount], FALSE);

        if (NULL != m_pOwner)
        {
            static const CCString sChannelNames[ChiralKSMax] =
            {
                _T("ChiralKS"), _T("ConnectSusp"),
                _T("Gamma1"), _T("Gamma2"), _T("Gamma3"), _T("Gamma4"),
                _T("Gamma5"), _T("Gamma51"), _T("Gamma52"), _T("Gamma53"), _T("Gamma54"),
                _T("Sigma12"), _T("Sigma13"), _T("Sigma14"),
                _T("Sigma23"), _T("Sigma24"), _T("Sigma34")
            };

            for (INT i = 0; i < static_cast<INT>(ChiralKSMax); ++i)
            {
                m_pOwner->AddOneConfigurationResult(this, _T("CondAll_") + sChannelNames[i], m_lstCondAll[i][m_uiConfigurationCount]);
                m_pOwner->AddOneConfigurationResult(this, _T("CondIn_") + sChannelNames[i], m_lstCondIn[i][m_uiConfigurationCount]);

                if (m_lstR.Num() > 0)
                {
                    const INT iRadialStart = m_lstCond[i].Num() - m_lstR.Num();
                    TArray<CLGComplex> radial;
                    for (INT r = 0; r < m_lstR.Num(); ++r)
                    {
                        radial.AddItem(m_lstCond[i][iRadialStart + r]);
                    }
                    m_pOwner->AddOneConfigurationResult(this, _T("CondRadial_") + sChannelNames[i], radial);
                }
            }
        }

        ++m_uiConfigurationCount;
    }
}

void CMeasureChiralCondensateKS::Report()
{
    appPushLogDate(FALSE);
    for (UINT i = 0; i < ChiralKSMax; ++i)
    {
        appAssert(m_uiConfigurationCount == static_cast<UINT>(m_lstCondAll[i].Num()));

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
    appPopLogDate();
}

void CMeasureChiralCondensateKS::Reset()
{
    CMeasureStochastic::Reset();
    for (UINT i = 0; i < ChiralKSMax; ++i)
    {
        m_lstCondAll[i].RemoveAll();
        m_lstCondIn[i].RemoveAll();
        m_lstCond[i].RemoveAll();
        m_lstCondXSlice[i].RemoveAll();
        m_lstCondYSlice[i].RemoveAll();
        m_lstCondZSlice[i].RemoveAll();
        m_lstCondTSlice[i].RemoveAll();
    }
    m_lstR.RemoveAll();
}

TArray<TArray<CLGComplex>> CMeasureChiralCondensateKS::ExportDiagnal(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* gaugeFields, const class CFieldBoson* const* bosonFields, const class CFieldTensor2* const* tensor2Fields, class CFieldFermion* pooled1, class CFieldFermion* pooled2)
{
    TArray<TArray<CLGComplex>> ret;
    CFieldFermionKS* pF1 = dynamic_cast<CFieldFermionKS*>(pooled1);
    CFieldFermionKS* pF2 = dynamic_cast<CFieldFermionKS*>(pooled2);
    if (NULL == pF1 || NULL == pF2)
    {
        appCrucial(_T("CMeasureChiralCondensateKS only work with CFieldFermionKS"));
        return ret;
    }
    
    UINT uiSiteCount = pF1->GetSiteCount();
    TArray<CLGComplex> rets[ChiralKSMax];
    for (UINT x = 0; x < uiSiteCount; ++x)
    {
        BYTE maxC = 0;
        switch (pF1->GetFieldType())
        {
        case EFT_FermionStaggeredSU3:
            maxC = 3;
            break;
        case EFT_FermionStaggeredU1:
            maxC = 1;
            break;
        default:
            appCrucial(_T("This is impossible!\n"));
            break;
        }

        for (BYTE c = 0; c < maxC; ++c)
        {
            SFermionBosonSource source;
            source.m_eSourceType = EFS_Point;
            source.m_byColorIndex = c;
            source.m_sSourcePoint = __hostSiteIndexToInt4(x);
            pF1->InitialAsSource(source);
            pF1->InverseD(gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields);

            for (BYTE i = 0; i < ChiralKSMax; ++i)
            {
                switch ((EChiralMeasureTypeKS)i)
                {
                case CMTKSGamma1:
                case CMTKSGamma2:
                case CMTKSGamma3:
                case CMTKSGamma4:
                case CMTKSGamma5:
                case CMTKSGamma51:
                case CMTKSGamma52:
                case CMTKSGamma53:
                case CMTKSGamma54:
                    {
                        pF1->CopyTo(pF2);
                        pF2->ApplyGamma(gaugeNum, bosonNum, gaugeFields, bosonFields, (EGammaMatrix)(i - 1));
                    }
                    break;
                case CMTKSSigma12:
                    {
                        pF1->CopyTo(pF2);
                        if (m_bMeasureSigma12)
                        {
                            pF2->ApplyGamma(gaugeNum, bosonNum, gaugeFields, bosonFields, SIGMA12);
                        }
                    }
                    break;
                case CMTKSSigma13:
                    {
                        pF1->CopyTo(pF2);
                        pF2->ApplyGamma(gaugeNum, bosonNum, gaugeFields, bosonFields, SIGMA31);
                    }
                    break;
                case CMTKSSigma14:
                    {
                        pF1->CopyTo(pF2);
                        pF2->ApplyGamma(gaugeNum, bosonNum, gaugeFields, bosonFields, SIGMA41);
                    }
                    break;
                case CMTKSSigma23:
                    {
                        pF1->CopyTo(pF2);
                        pF2->ApplyGamma(gaugeNum, bosonNum, gaugeFields, bosonFields, SIGMA23);
                    }
                    break;
                case CMTKSSigma24:
                    {
                        pF1->CopyTo(pF2);
                        pF2->ApplyGamma(gaugeNum, bosonNum, gaugeFields, bosonFields, SIGMA42);
                    }
                    break;
                case CMTKSSigma34:
                    {
                        pF1->CopyTo(pF2);
                        pF2->ApplyGamma(gaugeNum, bosonNum, gaugeFields, bosonFields, SIGMA43);
                    }
                    break;
                case ConnectSusp:
                    {
                        pF1->CopyTo(pF2);
                        if (m_bMeasureConnect)
                        {
                            pF2->InverseD(gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields);
                        }
                    }
                    break;
                default:
                    pF1->CopyTo(pF2);
                    break;
                }

                switch (pF1->GetFieldType())
                {
                case EFT_FermionStaggeredSU3:
                    {
                        CFieldFermionKSSU3* pF2SU3 = dynamic_cast<CFieldFermionKSSU3*>(pF2);
                        deviceSU3Vector hostv[1];
                        checkCudaErrors(cudaMemcpy(hostv, pF2SU3->m_pDeviceData + x, sizeof(deviceSU3Vector), cudaMemcpyDeviceToHost));
                        rets[i].AddItem(hostv->m_ve[c]);
                    }
                    break;
                case EFT_FermionStaggeredU1:
                    {
                        CFieldFermionKSU1* pF2U1 = dynamic_cast<CFieldFermionKSU1*>(pF2);
                        CLGComplex hostv[1];
                        checkCudaErrors(cudaMemcpy(hostv, pF2U1->m_pDeviceData + x, sizeof(CLGComplex), cudaMemcpyDeviceToHost));
                        rets[i].AddItem(hostv[0]);
                    }
                    break;
                default:
                    appCrucial(_T("not implemented yet!\n"));
                    break;
                }
            }
        }
    }

    for (BYTE i = 0; i < ChiralKSMax; ++i)
    {
        ret.AddItem(rets[i]);
    }
    return ret;
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================