//=============================================================================
// FILENAME : CMeasureMesonCorrelatorStaggeredSimple2.cpp
// 
// DESCRIPTION:
// This is the class for one measurement
//
// REVISION: [dd-mm-yy]
//  [08/12/2022 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Data/Field/Staggered/CFieldFermionKSSU3.h"
#include "CMeasureMesonCorrelatorStaggeredSimple2.h"

__BEGIN_NAMESPACE

#pragma region kernels

/**
 * sum _Ac1c2 phi(A) p1_c1c2 p2_c1c2*
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelPickPropagatorsSimple2(
    const deviceSU3Vector* __restrict__ propagatorf1r,
    const deviceSU3Vector* __restrict__ propagatorf1g,
    const deviceSU3Vector* __restrict__ propagatorf1b,
    const deviceSU3Vector* __restrict__ propagatorf2r,
    const deviceSU3Vector* __restrict__ propagatorf2g,
    const deviceSU3Vector* __restrict__ propagatorf2b,
    Real* res)
{
    intokernalInt4;
    if (0 == sSite4.w)
    {
        for (INT i = 0; i < static_cast<INT>(CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2); ++i)
        {
            res[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2 + i] = F(0.0);
        }

        return;
    }

    const deviceSU3Vector* pf1[3] = { propagatorf1r , propagatorf1g, propagatorf1b };
    const deviceSU3Vector* pf2[3] = { propagatorf2r , propagatorf2g, propagatorf2b };

    #pragma unroll
    for (BYTE byType = 0; byType < CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2; ++byType)
    {
        //pick propagator with a phase
        CLGComplex sumovercolor = _zeroc;
        #pragma unroll
        for (BYTE byC1 = 0; byC1 < 3; ++byC1)
        {
            #pragma unroll
            for (BYTE byC2 = 0; byC2 < 3; ++byC2)
            {
                sumovercolor = _cuCaddf(sumovercolor, _cuCmulf(pf1[byC1][uiSiteIndex].m_ve[byC2], _cuConjf(pf2[byC1][uiSiteIndex].m_ve[byC2])));
            }
        }

        const Real fPhase = static_cast<Real>(_deviceStaggeredFermionSimplePhase(sSite4, byType));
        res[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2 + byType] = _cuCabsf(sumovercolor) * fPhase;
    }
}

/**
 * sum over phase
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelPickEveryTimeSliceSimple2(
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
    res[uiVolumnIdx] = pAll[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2 + byType];
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CMeasureMesonCorrelatorStaggeredSimple2)


CMeasureMesonCorrelatorStaggeredSimple2::~CMeasureMesonCorrelatorStaggeredSimple2()
{
    checkCudaErrors(cudaFree(m_pDevicePropogators));

    appSafeFree(m_pResPropogators);
}

void CMeasureMesonCorrelatorStaggeredSimple2::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);
    INT iValue = 0;
    param.FetchValueINT(_T("FieldId2"), iValue);
    m_byFieldID2 = static_cast<BYTE>(iValue);

    checkCudaErrors(cudaMalloc((void**)&m_pDevicePropogators, sizeof(Real) * _HC_Volume * _kMesonCorrelatorTypeSimple2));
    if (HasOtherField())
    {
#if !_CLG_DOUBLEFLOAT
        m_pResPropogators = (DOUBLE*)malloc(sizeof(DOUBLE) * _kMesonCorrelatorTypeSimple2 * (_HC_Lt - 1) * 4);
#else
        m_pResPropogators = (Real*)malloc(sizeof(Real) * _kMesonCorrelatorTypeSimple2 * (_HC_Lt - 1) * 4);
#endif
    }
    else
    {
#if !_CLG_DOUBLEFLOAT
        m_pResPropogators = (DOUBLE*)malloc(sizeof(DOUBLE) * _kMesonCorrelatorTypeSimple2 * (_HC_Lt - 1));
#else
        m_pResPropogators = (Real*)malloc(sizeof(Real) * _kMesonCorrelatorTypeSimple2 * (_HC_Lt - 1));
#endif        
    }

    iValue = 1;
    param.FetchValueINT(_T("WallSource"), iValue);
    m_bWallSource = iValue != 0;
}

void CMeasureMesonCorrelatorStaggeredSimple2::OnConfigurationAccepted(INT gn, INT bn, const CFieldGauge* const* gs, const CFieldBoson* const* bs, const CFieldGauge* const* pStapleField)
{
    BuildSource();
    IniverseSource(gn, bn, gs, bs);

    preparethread;
    if (HasOtherField())
    {
        //Nf=1+1
        //uu
        _kernelPickPropagatorsSimple2 << <block, threads >> > (
            m_pSources[0]->m_pDeviceData,
            m_pSources[1]->m_pDeviceData,
            m_pSources[2]->m_pDeviceData,
            m_pSources[6]->m_pDeviceData,
            m_pSources[7]->m_pDeviceData,
            m_pSources[8]->m_pDeviceData,
            m_pDevicePropogators);

        dim3 block1 = dim3(block.x, block.y, 1);
        dim3 thread1 = dim3(threads.x, threads.y, 1);
        for (BYTE byType = 0; byType < _kMesonCorrelatorTypeSimple2; ++byType)
        {
            for (BYTE byT = 1; byT < _HC_Lt; ++byT)
            {
                _kernelPickEveryTimeSliceSimple2 << <block1, thread1 >> > (
                    m_pDevicePropogators, byT, byType, _D_RealThreadBuffer);
#if !_CLG_DOUBLEFLOAT
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyz);
#else
                const Real sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyz);
#endif

                m_pResPropogators[(byType * 4 + 0) * (_HC_Lt - 1) + byT - 1] = sum;
            }
        }

        //ud
        _kernelPickPropagatorsSimple2 << <block, threads >> > (
            m_pSources[0]->m_pDeviceData,
            m_pSources[1]->m_pDeviceData,
            m_pSources[2]->m_pDeviceData,
            m_pSources[9]->m_pDeviceData,
            m_pSources[10]->m_pDeviceData,
            m_pSources[11]->m_pDeviceData,
            m_pDevicePropogators);

        for (BYTE byType = 0; byType < _kMesonCorrelatorTypeSimple2; ++byType)
        {
            for (BYTE byT = 1; byT < _HC_Lt; ++byT)
            {
                _kernelPickEveryTimeSliceSimple2 << <block1, thread1 >> > (
                    m_pDevicePropogators, byT, byType, _D_RealThreadBuffer);
#if !_CLG_DOUBLEFLOAT
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyz);
#else
                const Real sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyz);
#endif

                m_pResPropogators[(byType * 4 + 1) * (_HC_Lt - 1) + byT - 1] = sum;
            }
        }

        //du
        _kernelPickPropagatorsSimple2 << <block, threads >> > (
            m_pSources[3]->m_pDeviceData,
            m_pSources[4]->m_pDeviceData,
            m_pSources[5]->m_pDeviceData,
            m_pSources[6]->m_pDeviceData,
            m_pSources[7]->m_pDeviceData,
            m_pSources[8]->m_pDeviceData,
            m_pDevicePropogators);

        for (BYTE byType = 0; byType < _kMesonCorrelatorTypeSimple2; ++byType)
        {
            for (BYTE byT = 1; byT < _HC_Lt; ++byT)
            {
                _kernelPickEveryTimeSliceSimple2 << <block1, thread1 >> > (
                    m_pDevicePropogators, byT, byType, _D_RealThreadBuffer);
#if !_CLG_DOUBLEFLOAT
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyz);
#else
                const Real sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyz);
#endif

                m_pResPropogators[(byType * 4 + 2) * (_HC_Lt - 1) + byT - 1] = sum;
            }
        }

        //dd
        _kernelPickPropagatorsSimple2 << <block, threads >> > (
            m_pSources[3]->m_pDeviceData,
            m_pSources[4]->m_pDeviceData,
            m_pSources[5]->m_pDeviceData,
            m_pSources[9]->m_pDeviceData,
            m_pSources[10]->m_pDeviceData,
            m_pSources[11]->m_pDeviceData,
            m_pDevicePropogators);

        for (BYTE byType = 0; byType < _kMesonCorrelatorTypeSimple2; ++byType)
        {
            for (BYTE byT = 1; byT < _HC_Lt; ++byT)
            {
                _kernelPickEveryTimeSliceSimple2 << <block1, thread1 >> > (
                    m_pDevicePropogators, byT, byType, _D_RealThreadBuffer);
#if !_CLG_DOUBLEFLOAT
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyz);
#else
                const Real sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyz);
#endif

                m_pResPropogators[(byType * 4 + 3) * (_HC_Lt - 1) + byT - 1] = sum;
            }
        }
    }
    else
    {
        //Nf=2
        _kernelPickPropagatorsSimple2 << <block, threads >> > (
            m_pSources[0]->m_pDeviceData,
            m_pSources[1]->m_pDeviceData,
            m_pSources[2]->m_pDeviceData,
            m_pSources[3]->m_pDeviceData,
            m_pSources[4]->m_pDeviceData,
            m_pSources[5]->m_pDeviceData,
            m_pDevicePropogators);

        dim3 block1 = dim3(block.x, block.y, 1);
        dim3 thread1 = dim3(threads.x, threads.y, 1);
        for (BYTE byType = 0; byType < _kMesonCorrelatorTypeSimple2; ++byType)
        {
            for (BYTE byT = 1; byT < _HC_Lt; ++byT)
            {
                _kernelPickEveryTimeSliceSimple2 << <block1, thread1 >> > (
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
    }

    ReleaseSource();

    //========== extract result ===========
    if (m_bShowResult)
    {
        appPushLogDate(FALSE);
        appGeneral(_T("==================== correlators ===============\n"));
    }
#if !_CLG_DOUBLEFLOAT
    TArray<TArray<DOUBLE>> thisConf;
#else
    TArray<TArray<Real>> thisConf;
#endif
    const INT totalType = HasOtherField() ? (_kMesonCorrelatorTypeSimple2 * 4) : _kMesonCorrelatorTypeSimple2;
    for (INT i = 0; i < totalType; ++i)
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
    if (m_bShowResult)
    {
        appPopLogDate();
    }

    ++m_uiConfigurationCount;
}

void CMeasureMesonCorrelatorStaggeredSimple2::BuildSource()
{
    m_pSources.RemoveAll();

    for (BYTE color = 0; color < 3; ++color)
    {
        CFieldFermionKSSU3* pFermion = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(GetFermionFieldId()));
        assert(NULL != pFermion);
        SFermionBosonSource source;
        source.m_byColorIndex = color;
        source.m_bySpinIndex = 0;
        source.m_eSourceType = m_bWallSource ? EFS_Wall : EFS_Point;
        source.m_sSourcePoint = SSmallInt4(0, 0, 0, 0);
        pFermion->InitialAsSource(source);

        m_pSources.AddItem(pFermion);
    }

    if (HasOtherField())
    {
        for (BYTE color = 0; color < 3; ++color)
        {
            CFieldFermionKSSU3* pFermion = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldID2));
            assert(NULL != pFermion);
            SFermionBosonSource source;
            source.m_byColorIndex = color;
            source.m_bySpinIndex = 0;
            source.m_eSourceType = m_bWallSource ? EFS_Wall : EFS_Point;
            source.m_sSourcePoint = SSmallInt4(0, 0, 0, 0);
            pFermion->InitialAsSource(source);

            m_pSources.AddItem(pFermion);
        }
    }
}

/**
 * if nf=2, the first three is D^{-1}s, the next three is Ddagger^{-1}s
 * if nf=1+1, the first 6 is D^{-1}(u,d), the next three is Ddagger^{-1}(u,d)s
 */
void CMeasureMesonCorrelatorStaggeredSimple2::IniverseSource(INT gn, INT bn, const CFieldGauge* const* gs, const CFieldBoson* const* bs)
{
    INT totalNumOfSource = m_pSources.Num();
    for (INT color = 0; color < totalNumOfSource; ++color)
    {
        m_pSources[color]->InverseDDdagger(gn, bn, gs, bs);
        //appGeneral(_T("byfield id: %d\n"), m_pSources[color]->m_byFieldId);
    }

    for (BYTE color = 0; color < totalNumOfSource; ++color)
    {
        CFieldFermionKSSU3* pFermion = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_pSources[color]->m_byFieldId));
        m_pSources[color]->CopyTo(pFermion);
        pFermion->D(gn, bn, gs, bs);
        m_pSources.AddItem(pFermion);
    }

    for (INT color = 0; color < totalNumOfSource; ++color)
    {
        m_pSources[color]->Ddagger(gn, bn, gs, bs);
    }
}

void CMeasureMesonCorrelatorStaggeredSimple2::ReleaseSource()
{
    for (INT i = 0; i < m_pSources.Num(); ++i)
    {
        m_pSources[i]->Return();
    }
    m_pSources.RemoveAll();
}

void CMeasureMesonCorrelatorStaggeredSimple2::Report()
{
    appPushLogDate(FALSE);
    appGeneral(_T(" =====================================================\n"));
    appGeneral(_T(" =================== Staggered Meson =================\n"));
    appGeneral(_T(" =====================================================\n\n"));
    m_lstAverageResults.RemoveAll();
    for (INT ty = 0; ty < _kMesonCorrelatorTypeSimple2; ++ty)
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
    for (INT ty = 0; ty < _kMesonCorrelatorTypeSimple2; ++ty)
    {
        appGeneral(_T("{"));
        for (INT t = 0; t < _HC_Lti - 1; ++t)
        {
            appGeneral(_T("%2.12f%s"), m_lstAverageResults[ty][t], (t != (_HC_Lti - 2)) ? _T(",") : _T(""));
        }
        appGeneral(_T("}%s"), ty == 19 ? _T("\n};\n") : _T(",\n"));
    }
    appPopLogDate();
}

void CMeasureMesonCorrelatorStaggeredSimple2::Reset()
{
    CMeasure::Reset();
    m_lstResults.RemoveAll();
    m_lstAverageResults.RemoveAll();
    ReleaseSource();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================