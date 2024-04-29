//=============================================================================
// FILENAME : CMeasureMesonCorrelatorStaggered.cpp
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION:
//  [09/28/2020 nbale]
//=============================================================================

#include "CLGLib_Private.h"

__BEGIN_NAMESPACE

#pragma region kernels

__global__ void _CLG_LAUNCH_BOUND
_kernelPickPropagatorsSimple(
    const deviceSU3Vector* __restrict__ propagator,
    Real* res)
{
    intokernalInt4;
    if (0 == sSite4.w)
    {
        res[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple::_kMesonCorrelatorTypeSimple + 0] = F(0.0);
        res[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple::_kMesonCorrelatorTypeSimple + 1] = F(0.0);
        res[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple::_kMesonCorrelatorTypeSimple + 2] = F(0.0);
        res[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple::_kMesonCorrelatorTypeSimple + 3] = F(0.0);
        return;
    }
    //pick propagator with a phase
    Real fV = F(0.0);
    #pragma unroll
    for (BYTE byC = 0; byC < 3; ++byC)
    {
        fV += __cuCabsSqf(propagator[uiSiteIndex].m_ve[0]);
        fV += __cuCabsSqf(propagator[uiSiteIndex].m_ve[1]);
        fV += __cuCabsSqf(propagator[uiSiteIndex].m_ve[2]);
    }

    #pragma unroll
    for (BYTE byType = 0; byType < CMeasureMesonCorrelatorStaggeredSimple::_kMesonCorrelatorTypeSimple; ++byType)
    {
        const Real fPhase = static_cast<Real>(_deviceStaggeredFermionSimplePhase(sSite4, byType));
        res[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple::_kMesonCorrelatorTypeSimple + byType] = fV * fPhase;
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPickEveryTimeSliceSimple(
    const Real* __restrict__ pAll,
    BYTE byT, BYTE byType,
#if !_CLG_DOUBLEFLOAT
    DOUBLE* res
#else
    Real* res
#endif
)
{
    const UINT uiVolumnIdx = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    const UINT uiSiteIndex = uiVolumnIdx * _DC_Lt + byT;
    res[uiVolumnIdx] = pAll[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple::_kMesonCorrelatorTypeSimple + byType];
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CMeasureMesonCorrelatorStaggeredSimple)


CMeasureMesonCorrelatorStaggeredSimple::~CMeasureMesonCorrelatorStaggeredSimple()
{
    checkCudaErrors(cudaFree(m_pDevicePropogators));

    appSafeFree(m_pResPropogators);
}

void CMeasureMesonCorrelatorStaggeredSimple::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);

    checkCudaErrors(cudaMalloc((void**)&m_pDevicePropogators, sizeof(Real) * _HC_Volume * _kMesonCorrelatorTypeSimple));
#if !_CLG_DOUBLEFLOAT
    m_pResPropogators = (DOUBLE*)malloc(sizeof(DOUBLE) * _kMesonCorrelatorTypeSimple * (_HC_Lt - 1));
#else
    m_pResPropogators = (Real*)malloc(sizeof(Real) * _kMesonCorrelatorTypeSimple * (_HC_Lt - 1));
#endif
}

void CMeasureMesonCorrelatorStaggeredSimple::OnConfigurationAcceptedSingleField(const CFieldGauge* pGaugeField, const CFieldGauge* pStapleField)
{
    CFieldFermionKSSU3* pFermion = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(GetFermionFieldId()));
    assert(NULL != pFermion);
    SFermionSource pointSource;
    pointSource.m_byColorIndex = 4;
    pointSource.m_eSourceType = EFS_Point;
    pointSource.m_sSourcePoint = SSmallInt4(0, 0, 0, 0);
    pFermion->InitialAsSource(pointSource);
    pFermion->InverseD(pGaugeField);

    preparethread;
    dim3 block1 = dim3(block.x, block.y, 1);
    dim3 thread1 = dim3(threads.x, threads.y, 1);
    _kernelPickPropagatorsSimple << <block, threads >> > (pFermion->m_pDeviceData, m_pDevicePropogators);

    for (BYTE byType = 0; byType < _kMesonCorrelatorTypeSimple; ++byType)
    {
        for (BYTE byT = 1; byT < _HC_Lt; ++byT)
        {
            _kernelPickEveryTimeSliceSimple << <block1, thread1 >> > (
                m_pDevicePropogators, byT, byType, _D_RealThreadBuffer);
#if !_CLG_DOUBLEFLOAT
            const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                _D_RealThreadBuffer, _HC_Volume_xyz);
#else
            const Real sum = appGetCudaHelper()->ReduceReal(
                _D_RealThreadBuffer, _HC_Volume_xyz);
#endif

            m_pResPropogators[byType * (_HC_Lt - 1) + byT - 1] = sum;
        }
    }
    pFermion->Return();

    //========== extract result ===========
    if (m_bShowResult)
    {
        appGeneral(_T("==================== correlators ===============\n"));
    }
#if !_CLG_DOUBLEFLOAT
    TArray<TArray<DOUBLE>> thisConf;
#else
    TArray<TArray<Real>> thisConf;
#endif
    for (INT i = 0; i < _kMesonCorrelatorTypeSimple; ++i)
    {
        if (m_bShowResult)
        {
            appGeneral(_T("Type%d:"), i);
        }
#if !_CLG_DOUBLEFLOAT
        TArray<DOUBLE> thisType;
#else
        TArray<Real> thisType;
#endif
        for (INT j = 0; j < _HC_Lti - 1; ++j)
        {
#if !_CLG_DOUBLEFLOAT
            const DOUBLE res = m_pResPropogators[i * (_HC_Lt - 1) + j];
#else
            const Real res = m_pResPropogators[i * (_HC_Lt - 1) + j];
#endif
            if (m_bShowResult)
            {
                appGeneral(_T("%2.12f, "), res);
            }
            thisType.AddItem(res);
        }
        thisConf.AddItem(thisType);
        if (m_bShowResult)
        {
            appGeneral(_T("\n"));
        }
    }
    m_lstResults.AddItem(thisConf);

    ++m_uiConfigurationCount;
}

void CMeasureMesonCorrelatorStaggeredSimple::Report()
{
    appGeneral(_T(" =====================================================\n"));
    appGeneral(_T(" =================== Staggered Meson =================\n"));
    appGeneral(_T(" =====================================================\n\n"));
    m_lstAverageResults.RemoveAll();
    for (INT ty = 0; ty < _kMesonCorrelatorTypeSimple; ++ty)
    {
        appGeneral(_T("(* ======================= Type:%d=================*)\ntabres%d={\n"), ty, ty);
#if !_CLG_DOUBLEFLOAT
        TArray<DOUBLE> thisType;
#else
        TArray<Real> thisType;
#endif
        for (INT conf = 0; conf < m_lstResults.Num(); ++conf)
        {
            appGeneral(_T("{"));
            for (INT t = 0; t < _HC_Lti - 1; ++t)
            {
                appGeneral(_T("%2.12f%s"), m_lstResults[conf][ty][t], (t != (_HC_Lti - 2)) ? _T(",") : _T(""));
                if (0 == conf)
                {
                    thisType.AddItem(m_lstResults[conf][ty][t]);
                }
                else
                {
                    thisType[t] = thisType[t] + m_lstResults[conf][ty][t];
                }
            }
            appGeneral(_T("}%s"), (conf == m_lstResults.Num() - 1) ? _T("\n};\n") : _T(",\n"));
        }

        for (INT t = 0; t < _HC_Lti - 1; ++t)
        {
            thisType[t] = thisType[t] / m_lstResults.Num();
        }
        m_lstAverageResults.AddItem(thisType);
    }


    appGeneral(_T("(* ======================= All Type averages =================*)\navr={\n"));
    for (INT ty = 0; ty < _kMesonCorrelatorTypeSimple; ++ty)
    {
        appGeneral(_T("{"));
        for (INT t = 0; t < _HC_Lti - 1; ++t)
        {
            appGeneral(_T("%2.12f%s"), m_lstAverageResults[ty][t], (t != (_HC_Lti - 2)) ? _T(",") : _T(""));
        }
        appGeneral(_T("}%s"), ty == 19 ? _T("\n};\n") : _T(",\n"));
    }
}

void CMeasureMesonCorrelatorStaggeredSimple::Reset()
{
    CMeasure::Reset();
    m_lstResults.RemoveAll();
    m_lstAverageResults.RemoveAll();
}




__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================