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
#include "Data/Field/Staggered/CFieldFermionKST.h"
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
    //if (0 == sSite4.w)
    //{
    //    for (INT i = 0; i < static_cast<INT>(CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2); ++i)
    //    {
    //        res[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2 + i] = F(0.0);
    //    }

    //    return;
    //}

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
    SCHAR uiT, BYTE byType,
    DOUBLE* res
)
{
    intokernalInt4_S(uiT);
    res[uiSiteIndex3D] = pAll[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2 + byType];
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPickEveryXSliceSimple2(
    const Real* __restrict__ pAll,
    SCHAR uiX, BYTE byType,
    DOUBLE* res
)
{
    intokernalInt4_Syzt(uiX);
    res[uiSiteIndex3DYZT] = pAll[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2 + byType];
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPickEveryYSliceSimple2(
    const Real* __restrict__ pAll,
    SCHAR uiY, BYTE byType,
    DOUBLE* res
)
{
    intokernalInt4_Sxzt(uiY);
    res[uiSiteIndex3DXZT] = pAll[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2 + byType];
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPickEveryZSliceSimple2(
    const Real* __restrict__ pAll,
    SCHAR uiZ, BYTE byType,
    DOUBLE* res
)
{
    intokernalInt4_Sxyt(uiZ);
    res[uiSiteIndex3DXYT] = pAll[uiSiteIndex * CMeasureMesonCorrelatorStaggeredSimple2::_kMesonCorrelatorTypeSimple2 + byType];
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CMeasureMesonCorrelatorStaggeredSimple2)


CMeasureMesonCorrelatorStaggeredSimple2::~CMeasureMesonCorrelatorStaggeredSimple2()
{
    checkCudaErrors(__cudaFree(m_pDevicePropogators));

    appSafeFree(m_pResPropogators);
    appSafeFree(m_pResPropogatorsX);
    appSafeFree(m_pResPropogatorsY);
    appSafeFree(m_pResPropogatorsZ);
}

void CMeasureMesonCorrelatorStaggeredSimple2::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);
    INT iValue = 0;
    param.FetchValueINT(_T("FieldId2"), iValue);
    m_byFieldID2 = static_cast<BYTE>(iValue);

    checkCudaErrors(__cudaMalloc((void**)&m_pDevicePropogators, sizeof(Real) * _HC_Volume * _kMesonCorrelatorTypeSimple2));
    if (HasOtherField())
    {
        m_pResPropogators = (DOUBLE*)malloc(sizeof(DOUBLE) * _kMesonCorrelatorTypeSimple2 * _HC_Lt * 4);
        m_pResPropogatorsX = (DOUBLE*)malloc(sizeof(DOUBLE) * _kMesonCorrelatorTypeSimple2 * _HC_Lx * 4);
        m_pResPropogatorsY = (DOUBLE*)malloc(sizeof(DOUBLE) * _kMesonCorrelatorTypeSimple2 * _HC_Ly * 4);
        m_pResPropogatorsZ = (DOUBLE*)malloc(sizeof(DOUBLE) * _kMesonCorrelatorTypeSimple2 * _HC_Lz * 4);
    }
    else
    {
        m_pResPropogators = (DOUBLE*)malloc(sizeof(DOUBLE) * _kMesonCorrelatorTypeSimple2 * _HC_Lt);
        m_pResPropogatorsX = (DOUBLE*)malloc(sizeof(DOUBLE) * _kMesonCorrelatorTypeSimple2 * _HC_Lx);
        m_pResPropogatorsY = (DOUBLE*)malloc(sizeof(DOUBLE) * _kMesonCorrelatorTypeSimple2 * _HC_Ly);
        m_pResPropogatorsZ = (DOUBLE*)malloc(sizeof(DOUBLE) * _kMesonCorrelatorTypeSimple2 * _HC_Lz);
    }

    CCString sValue = _T("EFS_Point");
    param.FetchStringValue(_T("Source"), sValue);
    m_eSource = __STRING_TO_ENUM(EFermionBosonSource, sValue);
}

void CMeasureMesonCorrelatorStaggeredSimple2::OnConfigurationAccepted(INT gn, INT bn, INT tensor2Num, const CFieldGauge* const* gs, const CFieldBoson* const* bs, const CFieldTensor2* const* tensor2Fields, const CFieldGauge* const* pStapleField)
{
    BuildSource();
    IniverseSource(gn, bn, gs, bs);

    #pragma region Time slice

    preparethread;
    preparethread_S;
    preparethread_Sxyt;
    preparethread_Sxzt;
    preparethread_Syzt;

    if (HasOtherField())
    {
        //Nf=1+1
        //uu
        _LAUNCH_KERNEL(_kernelPickPropagatorsSimple2, block, threads, 
            m_pSources[0]->m_pDeviceData,
            m_pSources[1]->m_pDeviceData,
            m_pSources[2]->m_pDeviceData,
            m_pSources[6]->m_pDeviceData,
            m_pSources[7]->m_pDeviceData,
            m_pSources[8]->m_pDeviceData,
            m_pDevicePropogators);

        for (BYTE byType = 0; byType < _kMesonCorrelatorTypeSimple2; ++byType)
        {
            for (SCHAR byT = 0; byT < _HC_Lti; ++byT)
            {
                _LAUNCH_KERNEL(_kernelPickEveryTimeSliceSimple2, block3d, threads3d, 
                    m_pDevicePropogators, byT, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyz);
                m_pResPropogators[(byType * 4 + 0) * _HC_Lt + byT] = sum;
            }

            for (SCHAR byX = 0; byX < _HC_Lxi; ++byX)
            {
                _LAUNCH_KERNEL(_kernelPickEveryXSliceSimple2, block3dyzt, threads3dyzt, 
                    m_pDevicePropogators, byX, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_yzt);
                m_pResPropogatorsX[(byType * 4 + 0) * _HC_Lx + byX] = sum;
            }

            for (SCHAR byY = 0; byY < _HC_Lyi; ++byY)
            {
                _LAUNCH_KERNEL(_kernelPickEveryYSliceSimple2, block3dxzt, threads3dxzt, 
                    m_pDevicePropogators, byY, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xzt);
                m_pResPropogatorsY[(byType * 4 + 0) * _HC_Ly + byY] = sum;
            }

            for (SCHAR byZ = 0; byZ < _HC_Lzi; ++byZ)
            {
                _LAUNCH_KERNEL(_kernelPickEveryZSliceSimple2, block3dxyt, threads3dxyt, 
                    m_pDevicePropogators, byZ, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyt);
                m_pResPropogatorsZ[(byType * 4 + 0) * _HC_Lz + byZ] = sum;
            }
        }

        //ud
        _LAUNCH_KERNEL(_kernelPickPropagatorsSimple2, block, threads, 
            m_pSources[0]->m_pDeviceData,
            m_pSources[1]->m_pDeviceData,
            m_pSources[2]->m_pDeviceData,
            m_pSources[9]->m_pDeviceData,
            m_pSources[10]->m_pDeviceData,
            m_pSources[11]->m_pDeviceData,
            m_pDevicePropogators);

        for (BYTE byType = 0; byType < _kMesonCorrelatorTypeSimple2; ++byType)
        {
            for (SCHAR byT = 0; byT < _HC_Lti; ++byT)
            {
                _LAUNCH_KERNEL(_kernelPickEveryTimeSliceSimple2, block3d, threads3d, 
                    m_pDevicePropogators, byT, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyz);
                m_pResPropogators[(byType * 4 + 1) * _HC_Lt + byT] = sum;
            }

            for (SCHAR byX = 0; byX < _HC_Lxi; ++byX)
            {
                _LAUNCH_KERNEL(_kernelPickEveryXSliceSimple2, block3dyzt, threads3dyzt, 
                    m_pDevicePropogators, byX, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_yzt);
                m_pResPropogatorsX[(byType * 4 + 1) * _HC_Lx + byX] = sum;
            }

            for (SCHAR byY = 0; byY < _HC_Lyi; ++byY)
            {
                _LAUNCH_KERNEL(_kernelPickEveryYSliceSimple2, block3dxzt, threads3dxzt, 
                    m_pDevicePropogators, byY, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xzt);
                m_pResPropogatorsY[(byType * 4 + 1) * _HC_Ly + byY] = sum;
            }

            for (SCHAR byZ = 0; byZ < _HC_Lzi; ++byZ)
            {
                _LAUNCH_KERNEL(_kernelPickEveryZSliceSimple2, block3dxyt, threads3dxyt, 
                    m_pDevicePropogators, byZ, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyt);
                m_pResPropogatorsZ[(byType * 4 + 1) * _HC_Lz + byZ] = sum;
            }
        }

        //du
        _LAUNCH_KERNEL(_kernelPickPropagatorsSimple2, block, threads, 
            m_pSources[3]->m_pDeviceData,
            m_pSources[4]->m_pDeviceData,
            m_pSources[5]->m_pDeviceData,
            m_pSources[6]->m_pDeviceData,
            m_pSources[7]->m_pDeviceData,
            m_pSources[8]->m_pDeviceData,
            m_pDevicePropogators);

        for (BYTE byType = 0; byType < _kMesonCorrelatorTypeSimple2; ++byType)
        {
            for (SCHAR byT = 0; byT < _HC_Lti; ++byT)
            {
                _LAUNCH_KERNEL(_kernelPickEveryTimeSliceSimple2, block3d, threads3d, 
                    m_pDevicePropogators, byT, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyz);
                m_pResPropogators[(byType * 4 + 2) * _HC_Lt + byT] = sum;
            }

            for (SCHAR byX = 0; byX < _HC_Lxi; ++byX)
            {
                _LAUNCH_KERNEL(_kernelPickEveryXSliceSimple2, block3dyzt, threads3dyzt, 
                    m_pDevicePropogators, byX, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_yzt);
                m_pResPropogatorsX[(byType * 4 + 2) * _HC_Lx + byX] = sum;
            }

            for (SCHAR byY = 0; byY < _HC_Lyi; ++byY)
            {
                _LAUNCH_KERNEL(_kernelPickEveryYSliceSimple2, block3dxzt, threads3dxzt, 
                    m_pDevicePropogators, byY, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xzt);
                m_pResPropogatorsY[(byType * 4 + 2) * _HC_Ly + byY] = sum;
            }

            for (SCHAR byZ = 0; byZ < _HC_Lzi; ++byZ)
            {
                _LAUNCH_KERNEL(_kernelPickEveryZSliceSimple2, block3dxyt, threads3dxyt, 
                    m_pDevicePropogators, byZ, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyt);
                m_pResPropogatorsZ[(byType * 4 + 2) * _HC_Lz + byZ] = sum;
            }
        }

        //dd
        _LAUNCH_KERNEL(_kernelPickPropagatorsSimple2, block, threads, 
            m_pSources[3]->m_pDeviceData,
            m_pSources[4]->m_pDeviceData,
            m_pSources[5]->m_pDeviceData,
            m_pSources[9]->m_pDeviceData,
            m_pSources[10]->m_pDeviceData,
            m_pSources[11]->m_pDeviceData,
            m_pDevicePropogators);

        for (BYTE byType = 0; byType < _kMesonCorrelatorTypeSimple2; ++byType)
        {
            for (SCHAR byT = 0; byT < _HC_Lti; ++byT)
            {
                _LAUNCH_KERNEL(_kernelPickEveryTimeSliceSimple2, block3d, threads3d, 
                    m_pDevicePropogators, byT, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyz);
                m_pResPropogators[(byType * 4 + 3) * _HC_Lt + byT] = sum;
            }

            for (SCHAR byX = 0; byX < _HC_Lxi; ++byX)
            {
                _LAUNCH_KERNEL(_kernelPickEveryXSliceSimple2, block3dyzt, threads3dyzt, 
                    m_pDevicePropogators, byX, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_yzt);
                m_pResPropogatorsX[(byType * 4 + 3) * _HC_Lx + byX] = sum;
            }

            for (SCHAR byY = 0; byY < _HC_Lyi; ++byY)
            {
                _LAUNCH_KERNEL(_kernelPickEveryYSliceSimple2, block3dxzt, threads3dxzt, 
                    m_pDevicePropogators, byY, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xzt);
                m_pResPropogatorsY[(byType * 4 + 3) * _HC_Ly + byY] = sum;
            }

            for (SCHAR byZ = 0; byZ < _HC_Lzi; ++byZ)
            {
                _LAUNCH_KERNEL(_kernelPickEveryZSliceSimple2, block3dxyt, threads3dxyt, 
                    m_pDevicePropogators, byZ, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyt);
                m_pResPropogatorsZ[(byType * 4 + 3) * _HC_Lz + byZ] = sum;
            }
        }
    }
    else
    {
        //Nf=2
        _LAUNCH_KERNEL(_kernelPickPropagatorsSimple2, block, threads, 
            m_pSources[0]->m_pDeviceData,
            m_pSources[1]->m_pDeviceData,
            m_pSources[2]->m_pDeviceData,
            m_pSources[3]->m_pDeviceData,
            m_pSources[4]->m_pDeviceData,
            m_pSources[5]->m_pDeviceData,
            m_pDevicePropogators);

        for (BYTE byType = 0; byType < _kMesonCorrelatorTypeSimple2; ++byType)
        {
            for (SCHAR byT = 0; byT < _HC_Lti; ++byT)
            {
                _LAUNCH_KERNEL(_kernelPickEveryTimeSliceSimple2, block3d, threads3d, 
                    m_pDevicePropogators, byT, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyz);
                m_pResPropogators[byType * _HC_Lt + byT] = sum;
            }

            for (SCHAR byX = 0; byX < _HC_Lxi; ++byX)
            {
                _LAUNCH_KERNEL(_kernelPickEveryXSliceSimple2, block3dyzt, threads3dyzt, 
                    m_pDevicePropogators, byX, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_yzt);
                m_pResPropogatorsX[byType * _HC_Lx + byX] = sum;
            }

            for (SCHAR byY = 0; byY < _HC_Lyi; ++byY)
            {
                _LAUNCH_KERNEL(_kernelPickEveryYSliceSimple2, block3dxzt, threads3dxzt, 
                    m_pDevicePropogators, byY, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xzt);
                m_pResPropogatorsY[byType * _HC_Ly + byY] = sum;
            }

            for (SCHAR byZ = 0; byZ < _HC_Lzi; ++byZ)
            {
                _LAUNCH_KERNEL(_kernelPickEveryZSliceSimple2, block3dxyt, threads3dxyt, 
                    m_pDevicePropogators, byZ, byType, _D_RealThreadBuffer);
                const DOUBLE sum = appGetCudaHelper()->ReduceReal(
                    _D_RealThreadBuffer, _HC_Volume_xyt);
                m_pResPropogatorsZ[byType * _HC_Lz + byZ] = sum;
            }
        }
    }

    #pragma endregion

    ReleaseSource();

#if _CLG_MULTI_GPU
    //P4-3.6: every per-t / per-x / per-y slice value is the LOCAL spatial
    //partial sum (ReduceReal over this rank's volume); sum across ranks so the
    //reported profile is global. Requires t/x/y NOT split (each rank then
    //holds the same slice index set); a split t/x/y would need a global-index
    //gather. The z-profile with a split z direction is rejected below.
    if (NULL != appGetComm())
    {
        if (appGetComm()->GpuGrid()[0] > 1 || appGetComm()->GpuGrid()[1] > 1 || appGetComm()->GpuGrid()[3] > 1)
        {
            appCrucial(_T("CMeasureMesonCorrelatorStaggeredSimple2: the t/x/y slice profiles are not supported on multi-GPU with a split t/x/y direction. Rejected.\n"));
            return;
        }
        const UINT uiFactor = HasOtherField() ? 4U : 1U;
        const UINT uiCountT = static_cast<UINT>(_kMesonCorrelatorTypeSimple2) * uiFactor * _HC_Lt;
        GlobalSumRealArray(m_pResPropogators, uiCountT);
        const UINT uiCountX = static_cast<UINT>(_kMesonCorrelatorTypeSimple2) * uiFactor * _HC_Lx;
        GlobalSumRealArray(m_pResPropogatorsX, uiCountX);
        const UINT uiCountY = static_cast<UINT>(_kMesonCorrelatorTypeSimple2) * uiFactor * _HC_Ly;
        GlobalSumRealArray(m_pResPropogatorsY, uiCountY);
        if (appGetComm()->GpuGrid()[2] > 1)
        {
            //z-profile needs a global-z gather; skip the reduction (keeps the
            //local partial, documented limitation) so the t/x/y profiles, which
            //are the primary correlator outputs, still reduce correctly.
            appGeneral(_T("CMeasureMesonCorrelatorStaggeredSimple2: z slice profile is local-only on multi-GPU with a split z direction.\n"));
        }
        else
        {
            const UINT uiCountZ = static_cast<UINT>(_kMesonCorrelatorTypeSimple2) * uiFactor * _HC_Lz;
            GlobalSumRealArray(m_pResPropogatorsZ, uiCountZ);
        }
    }
#endif

    //========== extract result ===========
    if (m_bShowResult)
    {
        appPushLogDate(FALSE);
        appGeneral(_T("==================== correlators ===============\n"));
    }
    TArray<TArray<DOUBLE>> thisConf;
    TArray<TArray<DOUBLE>> thisConfX;
    TArray<TArray<DOUBLE>> thisConfY;
    TArray<TArray<DOUBLE>> thisConfZ;
    const INT totalType = HasOtherField() ? (_kMesonCorrelatorTypeSimple2 * 4) : _kMesonCorrelatorTypeSimple2;
    for (INT i = 0; i < totalType; ++i)
    {
        if (m_bShowResult)
        {
            appGeneral(_T("Type%d:"), i);
        }
        TArray<DOUBLE> thisType;
        TArray<DOUBLE> thisTypeX;
        TArray<DOUBLE> thisTypeY;
        TArray<DOUBLE> thisTypeZ;

        for (INT j = 0; j < _HC_Lti; ++j)
        {
            const DOUBLE res = m_pResPropogators[i * _HC_Lt + j];
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

        //========= screen mass ==================
        for (INT j = 0; j < _HC_Lxi; ++j)
        {
            const DOUBLE res = m_pResPropogatorsX[i * _HC_Lx + j];
            if (m_bShowResult)
            {
                appGeneral(_T("%2.12f, "), res);
            }
            thisTypeX.AddItem(res);
        }
        thisConfX.AddItem(thisTypeX);
        if (m_bShowResult)
        {
            appGeneral(_T("\n"));
        }

        //========= screen mass ==================
        for (INT j = 0; j < _HC_Lyi; ++j)
        {
            const DOUBLE res = m_pResPropogatorsZ[i * _HC_Ly + j];
            if (m_bShowResult)
            {
                appGeneral(_T("%2.12f, "), res);
            }
            thisTypeY.AddItem(res);
        }
        thisConfY.AddItem(thisTypeY);
        if (m_bShowResult)
        {
            appGeneral(_T("\n"));
        }

        //========= screen mass ==================
        for (INT j = 0; j < _HC_Lzi; ++j)
        {
            const DOUBLE res = m_pResPropogatorsZ[i * _HC_Lz + j];
            if (m_bShowResult)
            {
                appGeneral(_T("%2.12f, "), res);
            }
            thisTypeZ.AddItem(res);
        }
        thisConfZ.AddItem(thisTypeZ);
        if (m_bShowResult)
        {
            appGeneral(_T("\n"));
        }
    }

    m_lstResults.AddItem(thisConf);
    m_lstResultsX.AddItem(thisConfX);
    m_lstResultsY.AddItem(thisConfY);
    m_lstResultsZ.AddItem(thisConfZ);
    if (NULL != m_pOwner)
    {
        m_pOwner->AddOneConfigurationResult(this, _T("CorrelatorT"), thisConf);
        m_pOwner->AddOneConfigurationResult(this, _T("CorrelatorX"), thisConfX);
        m_pOwner->AddOneConfigurationResult(this, _T("CorrelatorY"), thisConfY);
        m_pOwner->AddOneConfigurationResult(this, _T("CorrelatorZ"), thisConfZ);
    }
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
        CFieldFermionKSSU3* pFermion = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(GetFermionFieldId(), _T(__FILE__), __LINE__));
        appAssert(NULL != pFermion);
        SFermionBosonSource source;
        source.m_byColorIndex = color;
        source.m_bySpinIndex = 0;
        source.m_eSourceType = m_eSource;
        source.m_sSourcePoint = SSmallInt4(0, 0, 0, 0);
        pFermion->InitialAsSource(source);

        m_pSources.AddItem(pFermion);
    }

    if (HasOtherField())
    {
        for (BYTE color = 0; color < 3; ++color)
        {
            CFieldFermionKSSU3* pFermion = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldID2, _T(__FILE__), __LINE__));
            appAssert(NULL != pFermion);
            SFermionBosonSource source;
            source.m_byColorIndex = color;
            source.m_bySpinIndex = 0;
            source.m_eSourceType = m_eSource;
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
    appAssert(HasOtherField() ? (6 == totalNumOfSource) : (3 == totalNumOfSource));
    for (INT color = 0; color < totalNumOfSource; ++color)
    {
        m_pSources[color]->InverseDDdagger(gn, bn, 0, gs, bs, NULL);
        //appGeneral(_T("byfield id: %d\n"), m_pSources[color]->m_byFieldId);
    }

    //The 6-11(or 3-5 if only one field) was InverseDdagger
    for (BYTE color = 0; color < totalNumOfSource; ++color)
    {
        CFieldFermionKSSU3* pFermion = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_pSources[color]->m_byFieldId, _T(__FILE__), __LINE__));
        m_pSources[color]->CopyTo(pFermion);
        pFermion->D(gn, bn, 0, gs, bs, NULL);
        m_pSources.AddItem(pFermion);
    }

    //The 0-5(or 0-2 if only one field) was InverseD
    for (INT color = 0; color < totalNumOfSource; ++color)
    {
        m_pSources[color]->Ddagger(gn, bn, 0, gs, bs, NULL);
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
    m_lstAverageResultsX.RemoveAll();
    m_lstAverageResultsY.RemoveAll();
    m_lstAverageResultsZ.RemoveAll();

    for (INT ty = 0; ty < _kMesonCorrelatorTypeSimple2; ++ty)
    {
        appGeneral(_T("(* ======================= Type:%d=================*)\ntabres%d={\n"), ty, ty);
        TArray<DOUBLE> thisType;
        for (INT conf = 0; conf < m_lstResults.Num(); ++conf)
        {
            appGeneral(_T("{"));
            for (INT t = 0; t < _HC_Lti; ++t)
            {
                appGeneral(_T("%2.12f%s"), m_lstResults[conf][ty][t], (t != (_HC_Lti - 1)) ? _T(",") : _T(""));
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
        for (INT t = 0; t < _HC_Lti; ++t)
        {
            thisType[t] = thisType[t] / m_lstResults.Num();
        }
        m_lstAverageResults.AddItem(thisType);

        //========================= screen mass =============================
        thisType.RemoveAll();
        appGeneral(_T("(* ======================= X Type:%d=================*)\ntabresx%d={\n"), ty, ty);
        for (INT conf = 0; conf < m_lstResultsX.Num(); ++conf)
        {
            appGeneral(_T("{"));
            for (INT x = 0; x < _HC_Lxi; ++x)
            {
                appGeneral(_T("%2.12f%s"), m_lstResultsX[conf][ty][x], (x != (_HC_Lxi - 1)) ? _T(",") : _T(""));
                if (0 == conf)
                {
                    thisType.AddItem(m_lstResultsX[conf][ty][x]);
                }
                else
                {
                    thisType[x] = thisType[x] + m_lstResultsX[conf][ty][x];
                }
            }
            appGeneral(_T("}%s"), (conf == m_lstResultsX.Num() - 1) ? _T("\n};\n") : _T(",\n"));
        }
        for (INT x = 0; x < _HC_Lxi; ++x)
        {
            thisType[x] = thisType[x] / m_lstResultsX.Num();
        }
        m_lstAverageResultsX.AddItem(thisType);

        //========================= screen mass =============================
        thisType.RemoveAll();
        appGeneral(_T("(* ======================= Y Type:%d=================*)\ntabresy%d={\n"), ty, ty);
        for (INT conf = 0; conf < m_lstResultsY.Num(); ++conf)
        {
            appGeneral(_T("{"));
            for (INT y = 0; y < _HC_Lyi; ++y)
            {
                appGeneral(_T("%2.12f%s"), m_lstResultsY[conf][ty][y], (y != (_HC_Lyi - 1)) ? _T(",") : _T(""));
                if (0 == conf)
                {
                    thisType.AddItem(m_lstResultsY[conf][ty][y]);
                }
                else
                {
                    thisType[y] = thisType[y] + m_lstResultsY[conf][ty][y];
                }
            }
            appGeneral(_T("}%s"), (conf == m_lstResultsY.Num() - 1) ? _T("\n};\n") : _T(",\n"));
        }
        for (INT y = 0; y < _HC_Lyi; ++y)
        {
            thisType[y] = thisType[y] / m_lstResultsY.Num();
        }
        m_lstAverageResultsY.AddItem(thisType);

        //========================= screen mass =============================
        thisType.RemoveAll();
        appGeneral(_T("(* ======================= Z Type:%d=================*)\ntabresz%d={\n"), ty, ty);
        for (INT conf = 0; conf < m_lstResultsZ.Num(); ++conf)
        {
            appGeneral(_T("{"));
            for (INT z = 0; z < _HC_Lzi; ++z)
            {
                appGeneral(_T("%2.12f%s"), m_lstResultsZ[conf][ty][z], (z != (_HC_Lzi - 1)) ? _T(",") : _T(""));
                if (0 == conf)
                {
                    thisType.AddItem(m_lstResultsZ[conf][ty][z]);
                }
                else
                {
                    thisType[z] = thisType[z] + m_lstResultsZ[conf][ty][z];
                }
            }
            appGeneral(_T("}%s"), (conf == m_lstResultsZ.Num() - 1) ? _T("\n};\n") : _T(",\n"));
        }
        for (INT z = 0; z < _HC_Lzi; ++z)
        {
            thisType[z] = thisType[z] / m_lstResultsZ.Num();
        }
        m_lstAverageResultsZ.AddItem(thisType);
    }

    appGeneral(_T("(* ======================= All Type averages =================*)\navr={\n"));
    for (INT ty = 0; ty < _kMesonCorrelatorTypeSimple2; ++ty)
    {
        appGeneral(_T("{"));
        for (INT t = 0; t < _HC_Lti; ++t)
        {
            appGeneral(_T("%2.12f%s"), m_lstAverageResults[ty][t], (t != (_HC_Lti - 1)) ? _T(",") : _T(""));
        }
        appGeneral(_T("}%s"), ty == (static_cast<INT>(_kMesonCorrelatorTypeSimple2) - 1) ? _T("\n};\n") : _T(",\n"));
    }

    appGeneral(_T("(* ======================= All Type averages X =================*)\navrx={\n"));
    for (INT ty = 0; ty < _kMesonCorrelatorTypeSimple2; ++ty)
    {
        appGeneral(_T("{"));
        for (INT x = 0; x < _HC_Lxi; ++x)
        {
            appGeneral(_T("%2.12f%s"), m_lstAverageResultsX[ty][x], (x != (_HC_Lxi - 1)) ? _T(",") : _T(""));
        }
        appGeneral(_T("}%s"), ty == (static_cast<INT>(_kMesonCorrelatorTypeSimple2) - 1) ? _T("\n};\n") : _T(",\n"));
    }

    appGeneral(_T("(* ======================= All Type averages Y =================*)\navry={\n"));
    for (INT ty = 0; ty < _kMesonCorrelatorTypeSimple2; ++ty)
    {
        appGeneral(_T("{"));
        for (INT y = 0; y < _HC_Lyi; ++y)
        {
            appGeneral(_T("%2.12f%s"), m_lstAverageResultsY[ty][y], (y != (_HC_Lyi - 1)) ? _T(",") : _T(""));
        }
        appGeneral(_T("}%s"), ty == (static_cast<INT>(_kMesonCorrelatorTypeSimple2) - 1) ? _T("\n};\n") : _T(",\n"));
    }

    appGeneral(_T("(* ======================= All Type averages Z =================*)\navrz={\n"));
    for (INT ty = 0; ty < _kMesonCorrelatorTypeSimple2; ++ty)
    {
        appGeneral(_T("{"));
        for (INT z = 0; z < _HC_Lzi; ++z)
        {
            appGeneral(_T("%2.12f%s"), m_lstAverageResultsZ[ty][z], (z != (_HC_Lzi - 1)) ? _T(",") : _T(""));
        }
        appGeneral(_T("}%s"), ty == (static_cast<INT>(_kMesonCorrelatorTypeSimple2) - 1) ? _T("\n};\n") : _T(",\n"));
    }

    appPopLogDate();
}

void CMeasureMesonCorrelatorStaggeredSimple2::Reset()
{
    CMeasure::Reset();
    m_lstResults.RemoveAll();
    m_lstResultsX.RemoveAll();
    m_lstResultsY.RemoveAll();
    m_lstResultsZ.RemoveAll();
    m_lstAverageResults.RemoveAll();
    m_lstAverageResultsX.RemoveAll();
    m_lstAverageResultsY.RemoveAll();
    m_lstAverageResultsZ.RemoveAll();
    ReleaseSource();
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================