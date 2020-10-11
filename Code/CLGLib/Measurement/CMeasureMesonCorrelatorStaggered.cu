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

/**
 * block.x = Lt -  2
 * thread.x = 20
 */
//__global__ void _CLG_LAUNCH_BOUND
//_kernelInitialPropagatorsKS(CLGComplex* res)
//{
//    res[threadIdx.x * (_DC_Lt - 1) + blockIdx.x] = _zeroc;
//}

/**
 * block.x = 24
 * thread.x = 24
 * block.y = Lt -  1
 * block.z = 20
 */
//__global__ void _CLG_LAUNCH_BOUND
//_kernelPickPropagators(
//    const CLGComplex* __restrict__ propagators,
//    const BYTE* __restrict__ signTable,
//    const BYTE* __restrict__ deltaTable,
//    CLGComplex* res)
//{
//    //run for all A, B, C1, C2
//    const UINT uiA = blockIdx.x;
//    const UINT uiB = threadIdx.x;
//    const UINT uiAv = blockIdx.x / 3;
//    const UINT uiBv = threadIdx.x / 3;
//    //UINT byC1 = blockIdx.x % 3;
//    //UINT byC2 = threadIdx.x % 3;
//
//    const UINT byT = blockIdx.y;
//    const UINT uiType = blockIdx.z;
//    const UINT uiDelta = deltaTable[uiType];
//    const UINT byA_p_d = (uiAv ^ uiDelta) * 3 + (uiA % 3);
//    const UINT byB_p_d = (uiBv ^ uiDelta) * 3 + (uiB % 3);
//
//    //the coord with t should times 24
//    //B,0 to A,t this is propagator 2-3
//    const UINT uiPropagator1 = byT * 1152 + 576 + uiA * 24 + uiB;
//    //B_p_d,0 tp A_p_d,t this is propagator 1-4
//    const UINT uiPropagator2 = byT * 1152 + byA_p_d * 24 + byB_p_d;
//    const CLGComplex propagator = _cuCmulf(_cuConjf(propagators[uiPropagator1]), propagators[uiPropagator2]);
//    const INT sign = signTable[uiType * 8 + uiAv] + signTable[uiType * 8 + uiBv];
//    const UINT resIdx = uiType * (_DC_Lt - 1) + byT;
//    //if (5 == uiType && 1 == byT)
//    //{
//    //    printf("A=%d,B=%d,d=%d,sign=%d\n", uiAv, uiBv, uiDelta, sign);
//    //}
//    if (sign & 1)
//    {
//        //minus sign, but we have a total "-"
//        atomicAdd(&res[resIdx].x, -propagator.x);
//        atomicAdd(&res[resIdx].y, -propagator.y);
//    }
//    else
//    {
//        //plus sign, but we have a total "-"
//        atomicAdd(&res[resIdx].x, propagator.x);
//        atomicAdd(&res[resIdx].y, propagator.y);
//    }
//}

__global__ void _CLG_LAUNCH_BOUND
_kernelPickPropagators(
    const deviceSU3Vector* __restrict__ const * __restrict__ w1,
    const deviceSU3Vector* __restrict__ const * __restrict__ w2,
    const BYTE* __restrict__ signTable,
    const BYTE* __restrict__ deltaTable,
    BYTE byType, BYTE byFieldId,
    CLGComplex* res)
{
    intokernalInt4;

    res[uiSiteIndex] = _zeroc;
    if (0 == sSite4.w)
    {
        return;
    }

    if (0 != (sSite4.x & 1)
     || 0 != (sSite4.y & 1)
     || 0 != (sSite4.z & 1))
    {
        return;
    }

    const SBYTE byDelta = static_cast<SBYTE>(deltaTable[byType]);
    CLGComplex thisSiteSum = _zeroc;
    //we are picking for t = sSite4.w component
    //sum over x is sSite4.xyz
    #pragma unroll
    for (SBYTE byA1 = 0; byA1 < 8; ++byA1)
    {
        const SBYTE byA1_p_delta = byA1 ^ byDelta;

        const SBYTE a1x = byA1 & 1;
        const SBYTE a1y = (byA1 >> 1) & 1;
        const SBYTE a1z = (byA1 >> 2) & 1;

        const SBYTE a1pdx = byA1_p_delta & 1;
        const SBYTE a1pdy = (byA1_p_delta >> 1) & 1;
        const SBYTE a1pdz = (byA1_p_delta >> 2) & 1;

        #pragma unroll
        for (SBYTE byA2 = 0; byA2 < 8; ++byA2)
        {
            const SBYTE byA2v = byA2 * 3;
            const SBYTE byA2v_p_delta = (byA2 ^ byDelta) * 3;
            CLGComplex respropagator = _zeroc;

            #pragma unroll
            for (SBYTE byShift = 0; byShift < 8; ++byShift)
            {
                //shift.xyz = +-1
                const SBYTE shiftx = ((byShift & 1) << 1) - 1;
                const SBYTE shifty = (((byShift >> 1) & 1) << 1) - 1;
                const SBYTE shiftz = (((byShift >> 2) & 1) << 1) - 1;

                const SBYTE shifta1x = shiftx * a1x;
                const SBYTE shifta1y = shifty * a1y;
                const SBYTE shifta1z = shiftz * a1z;

                const SBYTE shifta1pdx = shiftx * a1pdx;
                const SBYTE shifta1pdy = shifty * a1pdy;
                const SBYTE shifta1pdz = shiftz * a1pdz;

                //Discard the terms outside of the unit cube
                //So we also do not need to consider the sites outside the lattice
                if (shifta1x < 0 || shifta1y < 0 || shifta1z < 0
                 || shifta1pdx < 0 || shifta1pdy < 0 || shifta1pdz < 0)
                {
                    continue;
                }
                //real shift
                SSmallInt4 shifted_A1 = sSite4;
                shifted_A1.x = shifted_A1.x + shifta1x;
                shifted_A1.y = shifted_A1.y + shifta1y;
                shifted_A1.z = shifted_A1.z + shifta1z;
                SSmallInt4 shifted_A1_p_delta = sSite4;
                shifted_A1_p_delta.x = shifted_A1_p_delta.x + shifta1pdx;
                shifted_A1_p_delta.y = shifted_A1_p_delta.y + shifta1pdy;
                shifted_A1_p_delta.z = shifted_A1_p_delta.z + shifta1pdz;

                //const UINT uiSiteShiftedA1 = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(shifted_A1)].m_uiSiteIndex;
                //const UINT uiSiteShiftedA1_p_delta = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__bi(shifted_A1_p_delta)].m_uiSiteIndex;;
                const UINT uiSiteShiftedA1 = _deviceGetSiteIndex(shifted_A1);
                const UINT uiSiteShiftedA1_p_delta = _deviceGetSiteIndex(shifted_A1_p_delta);

                #pragma unroll
                for (BYTE c = 0; c < 9; ++c)
                {
                    //sink(DA).w1(A2_p_d)
                    //sink(DA_p_d).w2(A2)
                    const BYTE c1 = c / 3;
                    const BYTE c2 = c % 3;
                    respropagator = _cuCaddf(respropagator,
                        _cuCmulf(
                            w1[byA2v_p_delta + c1][uiSiteShiftedA1_p_delta].m_ve[c2]
                            ,
                            _cuConjf(
                              w2[byA2v + c1][uiSiteShiftedA1].m_ve[c2]
                            )
                            )
                        );              
                }
            }

            const BYTE sign = signTable[byType * 8 + byA1] + signTable[byType * 8 + byA2];
            if (sign & 1)
            {
                thisSiteSum = _cuCaddf(thisSiteSum, respropagator);
            }
            else
            {
                thisSiteSum = _cuCsubf(thisSiteSum, respropagator);
            }
        }
    }
    //printf("just kankan %f + %f\n", thisSiteSum.x, thisSiteSum.y);
    res[uiSiteIndex] = thisSiteSum;
}

__global__ void _CLG_LAUNCH_BOUND
_kernelPickEveryTimeSlice(
    const CLGComplex* __restrict__ pAll,
    BYTE byT,
#if !_CLG_DOUBLEFLOAT
    cuDoubleComplex* res
#else
    CLGComplex* res
#endif
)
{
    const UINT uiVolumnIdx = (threadIdx.x + blockIdx.x * blockDim.x) * _DC_Lz + (threadIdx.y + blockIdx.y * blockDim.y);
    const UINT uiSiteIndex = uiVolumnIdx * _DC_Lt + byT;
#if !_CLG_DOUBLEFLOAT
    res[uiVolumnIdx] = _cToDouble(pAll[uiSiteIndex]);
#else
    res[uiVolumnIdx] = pAll[uiSiteIndex];
#endif
}

#pragma endregion

__CLGIMPLEMENT_CLASS(CMeasureMesonCorrelatorStaggered)

#pragma region Help functions

static inline INT __eta(INT x, INT y, INT z, INT i)
{
    switch (i)
    {
    case 1:
        return 0;
    case 2:
        return x;
    case 3:
        return x + y;
    case 4:
        return x + y + z;
    default:
        break;
    }
    return x + y + z;
}

static inline INT __xi(INT x, INT y, INT z, INT i)
{
    switch (i)
    {
    case 1:
        return y + z;
    case 2:
        return z;
    case 3:
        return 1;
    case 4:
        return 1;
    default:
        break;
    }
    return x + y + z;
}

void CMeasureMesonCorrelatorStaggered::InitialSignTable()
{
    for (INT i = 0; i < 20; ++i)
    {
        for (INT j = 0; j < 8; ++j)
        {
            const INT x = j & 1;
            const INT y = (j >> 1) & 1;
            const INT z = (j >> 2) & 1;
            const INT idx = i * 8 + j;
            //eta_1=xi_3=xi_4=0

            switch (i)
            {
            case 0:
                //1
                m_pSignTable[idx] = 0;
                break;
            case 1:
                //eta_4
                m_pSignTable[idx] = (__eta(x, y, z, 4)) & 1;
                break;
            case 2:
                //eta_2
                m_pSignTable[idx] = (__eta(x, y, z, 2)) & 1;
                break;
            case 3:
                //xi_1
                m_pSignTable[idx] = (__xi(x, y, z, 1)) & 1;
                break;
            case 4:
                //1
                m_pSignTable[idx] = 0;
                break;
            case 5:
                //eta_4
                m_pSignTable[idx] = (__eta(x, y, z, 4)) & 1;
                break;
            case 6:
                //eta_2
                m_pSignTable[idx] = (__eta(x, y, z, 2)) & 1;
                break;
            case 7:
                //xi_1
                m_pSignTable[idx] = (__xi(x, y, z, 1)) & 1;
            case 8:
                //1
                m_pSignTable[idx] = 0;
                break;
            case 9:
                //eta_4
                m_pSignTable[idx] = (__eta(x, y, z, 4)) & 1;
                break;
            case 10:
                //eta_2
                m_pSignTable[idx] = (__eta(x, y, z, 2)) & 1;
                break;
            case 11:
                //xi_1
                m_pSignTable[idx] = (__xi(x, y, z, 1)) & 1;
                break;
            case 12:
                //xi_1 xi_2
                m_pSignTable[idx] = (__xi(x, y, z, 1) + __xi(x, y, z, 2)) & 1;
                break;
            case 13:
                //eta_2 xi_2
                m_pSignTable[idx] = (__eta(x, y, z, 2) + __xi(x, y, z, 2)) & 1;
                break;
            case 14:
                //eta_4
                m_pSignTable[idx] = (__eta(x, y, z, 4)) & 1;
                break;
            case 15:
                //1
                m_pSignTable[idx] = 0;
                break;
            case 16:
                //eta_2 eta_3
                m_pSignTable[idx] = (__eta(x, y, z, 2) + __eta(x, y, z, 3)) & 1;
                break;
            case 17:
                //eta_2 xi_2
                m_pSignTable[idx] = (__eta(x, y, z, 2) + __xi(x, y, z, 2)) & 1;
                break;
            case 18:
                //eta_3
                m_pSignTable[idx] = (__eta(x, y, z, 3)) & 1;
                break;
            case 19:
                //xi_2
                m_pSignTable[idx] = (__xi(x, y, z, 2)) & 1;
                break;
            default:
                break;
            }
        }

        switch (i)
        {
        case 0:
            m_pDeltaTable[i] = 0;
            break;
        case 1:
            m_pDeltaTable[i] = 0;
            break;
        case 2:
            m_pDeltaTable[i] = 0;
            break;
        case 3:
            m_pDeltaTable[i] = 0;
            break;
        case 4:
            m_pDeltaTable[i] = 1;
            break;
        case 5:
            m_pDeltaTable[i] = 1;
            break;
        case 6:
            m_pDeltaTable[i] = 1;
            break;
        case 7:
            m_pDeltaTable[i] = 1;
            break;
        case 8:
            m_pDeltaTable[i] = 2;
            break;
        case 9:
            m_pDeltaTable[i] = 2;
            break;
        case 10:
            m_pDeltaTable[i] = 3;
            break;
        case 11:
            m_pDeltaTable[i] = 3;
            break;
        case 12:
            m_pDeltaTable[i] = 3;
            break;
        case 13:
            m_pDeltaTable[i] = 3;
            break;
        case 14:
            m_pDeltaTable[i] = 3;
            break;
        case 15:
            m_pDeltaTable[i] = 3;
            break;
        case 16:
            m_pDeltaTable[i] = 7;
            break;
        case 17:
            m_pDeltaTable[i] = 7;
            break;
        case 18:
            m_pDeltaTable[i] = 7;
            break;
        case 19:
            m_pDeltaTable[i] = 7;
            break;
        default:
            break;
        }
    }

    checkCudaErrors(cudaMemcpy(m_pDeviceSignTable, m_pSignTable, sizeof(BYTE) * 160, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(m_pDeviceDeltaTable, m_pDeltaTable, sizeof(BYTE) * 20, cudaMemcpyHostToDevice));
}

void CMeasureMesonCorrelatorStaggered::CalculateSources(const CFieldGauge* pGauge)
{
    for (BYTE shift = 0; shift < 8; ++shift)
    {
        for (BYTE c = 0; c < 3; ++c)
        {
            const INT idx = shift * 3 + c;
            SFermionSource source;
            source.m_bySpinIndex = shift;
            source.m_byColorIndex = c;
            source.m_eSourceType = EFS_Wall;
            source.m_sSourcePoint = SSmallInt4(0, 0, 0, 0);
            m_pW1[idx]->InitialAsSource(source);
            appParanoiac(_T("generating source(%d)...\n"), idx);
            m_pW1[idx]->InverseDDdagger(pGauge);
            m_pW1[idx]->CopyTo(m_pW2[idx]);

            m_pW1[idx]->Ddagger(pGauge);
            m_pW2[idx]->D(pGauge);
        }
    }
}

void CMeasureMesonCorrelatorStaggered::CalculatePropogators()
{
    deviceSU3Vector* w1[24];
    deviceSU3Vector* w2[24];
    for (BYTE shift = 0; shift < 24; ++shift)
    {
        w1[shift] = m_pW1[shift]->m_pDeviceData;
        w2[shift] = m_pW2[shift]->m_pDeviceData;
    }
    checkCudaErrors(cudaMemcpy(m_pDeviceW1, w1, sizeof(deviceSU3Vector*) * 24, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(m_pDeviceW2, w2, sizeof(deviceSU3Vector*) * 24, cudaMemcpyHostToDevice));

    preparethread;
    dim3 block1 = dim3(block.x, block.y, 1);
    dim3 thread1 = dim3(threads.x, threads.y, 1);
    for (INT i = 0; i < _kMesonCorrelatorType; ++i)
    {
        _kernelPickPropagators << <block, threads >> > (
            m_pDeviceW1,
            m_pDeviceW2,
            m_pDeviceSignTable,
            m_pDeviceDeltaTable,
            static_cast<BYTE>(i),
            m_byFieldId,
            m_pDevicePropogators);


        for (INT t = 1; t < _HC_Lti; ++t)
        {
            _kernelPickEveryTimeSlice << <block1, thread1 >> > (
                m_pDevicePropogators, 
                static_cast<BYTE>(t), 
                m_pDevicePropogatorsEveryTimeSlice);
#if !_CLG_DOUBLEFLOAT
            const cuDoubleComplex sum = appGetCudaHelper()->ReduceComplex(
                m_pDevicePropogatorsEveryTimeSlice, _HC_Volume_xyz);
#else
            const CLGComplex sum = appGetCudaHelper()->ReduceComplex(
                m_pDevicePropogatorsEveryTimeSlice, _HC_Volume_xyz);
#endif
            m_pResPropogators[i * (_HC_Lti - 1) + t - 1] = sum;
            //appParanoiac(_T("Type=%d t=%d, res = %2.18f + %2.18f\n"), i, t, sum.x, sum.y);
        }
    }
}

void CMeasureMesonCorrelatorStaggered::InitialBuffers()
{
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceSignTable, sizeof(BYTE) * _kMesonCorrelatorType * 8));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDeltaTable, sizeof(BYTE) * _kMesonCorrelatorType));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceW1, sizeof(deviceSU3Vector*) * 24));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceW2, sizeof(deviceSU3Vector*) * 24));
    checkCudaErrors(cudaMalloc((void**)&m_pDevicePropogators, sizeof(CLGComplex) * _HC_Volume));
#if !_CLG_DOUBLEFLOAT
    checkCudaErrors(cudaMalloc((void**)&m_pDevicePropogatorsEveryTimeSlice, sizeof(cuDoubleComplex) * _HC_Volume_xyz));
    m_pResPropogators = (cuDoubleComplex*)(malloc(sizeof(cuDoubleComplex) * _kMesonCorrelatorType * (_HC_Lt - 1)));
#else
    checkCudaErrors(cudaMalloc((void**)&m_pDevicePropogatorsEveryTimeSlice, sizeof(CLGComplex) * _HC_Volume_xyz));
    m_pResPropogators = (CLGComplex*)(malloc(sizeof(CLGComplex) * _kMesonCorrelatorType * (_HC_Lt - 1)));
#endif
    
}

#pragma endregion

CMeasureMesonCorrelatorStaggered::~CMeasureMesonCorrelatorStaggered()
{
    checkCudaErrors(cudaFree(m_pDeviceW1));
    checkCudaErrors(cudaFree(m_pDeviceW2));
    checkCudaErrors(cudaFree(m_pDeviceSignTable));
    checkCudaErrors(cudaFree(m_pDeviceDeltaTable));
    checkCudaErrors(cudaFree(m_pDevicePropogators));
    checkCudaErrors(cudaFree(m_pDevicePropogatorsEveryTimeSlice));

    appSafeFree(m_pResPropogators);
    appSafeDelete(m_pGaugeFixing);
}

void CMeasureMesonCorrelatorStaggered::Initial(CMeasurementManager* pOwner, CLatticeData* pLatticeData, const CParameters& param, BYTE byId)
{
    CMeasure::Initial(pOwner, pLatticeData, param, byId);
    INT iValue = 1;
    param.FetchValueINT(_T("ShowResult"), iValue);
    m_bShowResult = iValue != 0;

    iValue = 1;
    param.FetchValueINT(_T("GaugeFixing"), iValue);
    m_bGaugeFixing = iValue != 0;

    InitialBuffers();
    InitialSignTable();
}

void CMeasureMesonCorrelatorStaggered::OnConfigurationAccepted(const CFieldGauge* pGaugeField, const CFieldGauge* pStapleField)
{
    for (BYTE i = 0; i < 24; ++i)
    {
        m_pW1[i] = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
        m_pW2[i] = dynamic_cast<CFieldFermionKSSU3*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    }

    if (m_bGaugeFixing)
    {
        if (NULL == m_pGaugeFixing)
        {
            m_pGaugeFixing = dynamic_cast<CFieldGauge*>(pGaugeField->GetCopy());
        }
        else
        {
            pGaugeField->CopyTo(m_pGaugeFixing);
        }
        if (NULL != appGetLattice()->m_pGaugeFixing)
        {
            appGetLattice()->m_pGaugeFixing->GaugeFixing(m_pGaugeFixing);
        }
        CalculateSources(m_pGaugeFixing);
        CalculatePropogators();
    }
    else
    {
        CalculateSources(pGaugeField);
        CalculatePropogators();
    }

    //========== extract result ===========
    if (m_bShowResult)
    {
        appGeneral(_T("==================== correlators ===============\n"));
    }
#if !_CLG_DOUBLEFLOAT
    TArray<TArray<cuDoubleComplex>> thisConf;
#else
    TArray<TArray<CLGComplex>> thisConf;
#endif
    for (INT i = 0; i < _kMesonCorrelatorType; ++i)
    {
        if (m_bShowResult)
        {
            appGeneral(_T("Type%d:"), i);
        }
#if !_CLG_DOUBLEFLOAT
        TArray<cuDoubleComplex> thisType;
#else
        TArray<CLGComplex> thisType;
#endif
        for (INT j = 0; j < _HC_Lti - 1; ++j)
        {
#if !_CLG_DOUBLEFLOAT
            const cuDoubleComplex& res = m_pResPropogators[i * (_HC_Lt - 1) + j];
#else
            const CLGComplex& res = m_pResPropogators[i * (_HC_Lt - 1) + j];
#endif
            if (m_bShowResult)
            {
                LogGeneralComplex(res);
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

    for (BYTE i = 0; i < 24; ++i)
    {
        m_pW1[i]->Return();
        m_pW2[i]->Return();
    }
    ++m_uiConfigurationCount;
}

void CMeasureMesonCorrelatorStaggered::Average(UINT)
{

}

void CMeasureMesonCorrelatorStaggered::Report()
{
    appGeneral(_T(" =====================================================\n"));
    appGeneral(_T(" =================== Staggered Meson =================\n"));
    appGeneral(_T(" =====================================================\n\n"));
    m_lstAverageResults.RemoveAll();
    for (INT ty = 0; ty < _kMesonCorrelatorType; ++ty)
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
                appGeneral(_T("%2.12f%s"), m_lstResults[conf][ty][t].x, (t != (_HC_Lti - 2)) ? _T(",") : _T(""));
                if (0 == conf)
                {
                    thisType.AddItem(m_lstResults[conf][ty][t].x);
                }
                else
                {
                    thisType[t] = thisType[t] + m_lstResults[conf][ty][t].x;
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
    for (INT ty = 0; ty < _kMesonCorrelatorType; ++ty)
    {
        appGeneral(_T("{"));
        for (INT t = 0; t < _HC_Lti - 1; ++t)
        {
            appGeneral(_T("%2.12f%s"), m_lstAverageResults[ty][t], (t != (_HC_Lti - 2)) ? _T(",") : _T(""));
        }
        appGeneral(_T("}%s"), ty == 19 ? _T("\n};\n") : _T(",\n"));
    }
}

void CMeasureMesonCorrelatorStaggered::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstResults.RemoveAll();
    m_lstAverageResults.RemoveAll();
}




__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================