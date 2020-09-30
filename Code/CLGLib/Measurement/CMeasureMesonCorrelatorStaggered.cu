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
__global__ void _CLG_LAUNCH_BOUND
_kernelInitialPropagatorsKS(CLGComplex* res)
{
    res[threadIdx.x * (_DC_Lt - 1) + blockIdx.x] = _zeroc;
}

/**
 * block.x = 24
 * thread.x = 24
 * block.y = Lt -  2
 * block.z = 20
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelPickPropagators(
    const CLGComplex* __restrict__ propagators,
    const BYTE* __restrict__ signTable,
    const BYTE* __restrict__ deltaTable,
    CLGComplex* res)
{
    //run for all A, B, C1, C2
    const UINT uiA = blockIdx.x;
    const UINT uiB = threadIdx.x;
    const UINT uiAv = blockIdx.x / 3;
    const UINT uiBv = threadIdx.x / 3;
    //UINT byC1 = blockIdx.x % 3;
    //UINT byC2 = threadIdx.x % 3;

    const UINT byT = blockIdx.y;
    const UINT uiType = blockIdx.z;
    const UINT uiDelta = deltaTable[uiType];
    const UINT byA_p_d = (uiAv ^ uiDelta) * 3 + (uiA % 3);
    const UINT byB_p_d = (uiBv ^ uiDelta) * 3 + (uiB % 3);

    //the coord with t should times 24
    //B,0 to A,t this is propagator 2-3
    const UINT uiPropagator1 = byT * 1152 + 576 + uiA * 24 + uiB;
    //B_p_d,0 tp A_p_d,t this is propagator 1-4
    const UINT uiPropagator2 = byT * 1152 + byA_p_d * 24 + byB_p_d;
    const CLGComplex propagator = _cuCmulf(_cuConjf(propagators[uiPropagator1]), propagators[uiPropagator2]);
    const INT sign = signTable[uiType * 8 + uiAv] + signTable[uiType * 8 + uiBv];
    const UINT resIdx = uiType * (_DC_Lt - 1) + byT;

    if (sign & 1)
    {
        //minus sign, but we have a total "-"
        atomicAdd(&res[resIdx].x, propagator.x);
        atomicAdd(&res[resIdx].y, propagator.y);
    }
    else
    {
        //plus sign, but we have a total "-"
        atomicAdd(&res[resIdx].x, -propagator.x);
        atomicAdd(&res[resIdx].y, -propagator.y);
    }
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
            source.m_eSourceType = EFS_Wall;
            source.m_byColorIndex = c;
            source.m_sSourcePoint = SSmallInt4(0, 0, 0, 0);
            m_pSources[idx]->InitialAsSource(source);
            appParanoiac(_T("generating source(%d)...\n"), idx);
            m_pSources[idx]->InverseDDdagger(pGauge);
        }
    }
}

void CMeasureMesonCorrelatorStaggered::CalculatePropogators(const CFieldGauge* pGauge)
{
    m_pContract = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    m_pToBeContract = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
    for (INT i = 1; i < _HC_Lti; ++i)
    {
        appParanoiac(_T("calculating propogator for lt=(%d)...\n"), i);
        //To optimize: contract(shif1,c1).source(shift2,c2)=contract(shif2,c1).source(shift1,c2)
        //Can we do this optimization?
        for (BYTE shift = 0; shift < 8; ++shift)
        {
            for (BYTE c = 0; c < 3; ++c)
            {
                const INT idx = shift * 3 + c;
                SFermionSource source;
                source.m_bySpinIndex = shift;
                source.m_eSourceType = EFS_Wall;
                source.m_byColorIndex = c;
                source.m_sSourcePoint = SSmallInt4(0, 0, 0, static_cast<SBYTE>(i));
                m_pContract->InitialAsSource(source);

                for (BYTE targetIdx = 0; targetIdx < 24; ++targetIdx)
                //for (BYTE targetShift = 0; targetShift < 8; ++targetShift)
                {
                    //for (BYTE targetC = 0; targetC < 3; ++targetC)
                    //{
                        //const BYTE targetIdx = targetShift * 3 + targetC;
                        m_pSources[targetIdx]->CopyTo(m_pToBeContract);
                        m_pToBeContract->Ddagger(pGauge);
                        const CLGComplex prop14 = m_pContract->Dot(m_pToBeContract);
                        m_pSources[targetIdx]->CopyTo(m_pToBeContract);
                        m_pToBeContract->D(pGauge);
                        const CLGComplex prop23 = m_pContract->Dot(m_pToBeContract);

                        const INT totalIdx1 = idx * 24 + targetIdx;
                        m_pPropogators[(i - 1) * 1152 + totalIdx1] = prop14;
                        m_pPropogators[(i - 1) * 1152 + totalIdx1 + 576] = prop23;
                    //}
                }
            }
        }
    }
    m_pContract->Return();
    m_pToBeContract->Return();

    //appGeneral(_T("{\n"));
    //for (INT i = 0; i < 24; ++i)
    //{
    //    appGeneral(_T("{"));
    //    for (INT j = 0; j < 24; ++j)
    //    {
    //        LogGeneralComplex(m_pPropogators[1152 + i * 24 + j], 23 != j);
    //    }
    //    if (23 == i)
    //    {
    //        appGeneral(_T("}\n}\n"));
    //    }
    //    else
    //    {
    //        appGeneral(_T("},\n"));
    //    }
    //}
    //appGeneral(_T("{\n"));
    //for (INT i = 0; i < 24; ++i)
    //{
    //    appGeneral(_T("{"));
    //    for (INT j = 0; j < 24; ++j)
    //    {
    //        LogGeneralComplex(m_pPropogators[1152 + 576 + i * 24 + j], 23 != j);
    //    }
    //    if (23 == i)
    //    {
    //        appGeneral(_T("}\n}\n"));
    //    }
    //    else
    //    {
    //        appGeneral(_T("},\n"));
    //    }
    //}
}

void CMeasureMesonCorrelatorStaggered::PickPropagatorAndSign()
{
    dim3 block0(_HC_Lt - 1, 1, 1);
    dim3 threads0(_kMesonCorrelatorType, 1, 1);
    _kernelInitialPropagatorsKS << <block0, threads0 >> > (m_pDeviceResultPropogators);

    checkCudaErrors(cudaMemcpy(m_pDevicePropogators, m_pPropogators, sizeof(CLGComplex) * 1152 * (_HC_Lt - 1), cudaMemcpyHostToDevice));
    dim3 block(24, _HC_Lt - 1, _kMesonCorrelatorType);
    dim3 threads(24, 1, 1);

    _kernelPickPropagators << <block, threads >> > (m_pDevicePropogators, m_pDeviceSignTable, m_pDeviceDeltaTable, m_pDeviceResultPropogators);

    checkCudaErrors(cudaMemcpy(m_pFinalPropagators, m_pDeviceResultPropogators, sizeof(CLGComplex) * _kMesonCorrelatorType * (_HC_Lt - 1), cudaMemcpyDeviceToHost));
}

void CMeasureMesonCorrelatorStaggered::InitialBuffers()
{
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceSignTable, sizeof(BYTE) * _kMesonCorrelatorType * 8));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceDeltaTable, sizeof(BYTE) * _kMesonCorrelatorType));
    checkCudaErrors(cudaMalloc((void**)&m_pDevicePropogators, sizeof(CLGComplex) * 1152 * (_HC_Lt - 1)));
    checkCudaErrors(cudaMalloc((void**)&m_pDeviceResultPropogators, sizeof(CLGComplex) * _kMesonCorrelatorType * (_HC_Lt - 1)));

    m_pPropogators = (CLGComplex*)(malloc(sizeof(CLGComplex) * 1152 * (_HC_Lt - 1)));
    m_pFinalPropagators = (CLGComplex*)(malloc(sizeof(CLGComplex) * _kMesonCorrelatorType * (_HC_Lt - 1)));
}

#pragma endregion

CMeasureMesonCorrelatorStaggered::~CMeasureMesonCorrelatorStaggered()
{
    checkCudaErrors(cudaFree(m_pDeviceSignTable));
    checkCudaErrors(cudaFree(m_pDeviceDeltaTable));
    checkCudaErrors(cudaFree(m_pDevicePropogators));
    checkCudaErrors(cudaFree(m_pDeviceResultPropogators));

    appSafeFree(m_pPropogators);
    appSafeFree(m_pFinalPropagators);
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
        m_pSources[i] = dynamic_cast<CFieldFermion*>(appGetLattice()->GetPooledFieldById(m_byFieldId));
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
        CalculatePropogators(m_pGaugeFixing);
    }
    else
    {
        CalculateSources(pGaugeField);
        CalculatePropogators(pGaugeField);
    }

    PickPropagatorAndSign();

    //========== extract result ===========
    if (m_bShowResult)
    {
        appGeneral(_T("==================== correlators ===============\n"));
    }
    TArray<TArray<CLGComplex>> thisConf;
    for (INT i = 0; i < _kMesonCorrelatorType; ++i)
    {
        if (m_bShowResult)
        {
            appGeneral(_T("Type%d:"), i);
        }
        TArray<CLGComplex> thisType;
        for (INT j = 0; j < _HC_Lti - 1; ++j)
        {
            const CLGComplex& res = m_pFinalPropagators[i * (_HC_Lt - 1) + j];
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
        m_pSources[i]->Return();
    }
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
        appGeneral(_T(" ======================= Type:%d=================\n{\n"), ty);
        TArray<Real> thisType;
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
            appGeneral(_T("}%s"), (conf == m_lstResults.Num() - 1) ? _T("\n}\n") : _T(",\n"));
        }

        for (INT t = 0; t < _HC_Lti - 1; ++t)
        {
            thisType[t] = thisType[t] / m_lstResults.Num();
        }
        m_lstAverageResults.AddItem(thisType);
    }


    appGeneral(_T(" ======================= All Type averages =================\n{\n"));
    for (INT ty = 0; ty < _kMesonCorrelatorType; ++ty)
    {
        appGeneral(_T("{"));
        for (INT t = 0; t < _HC_Lti - 1; ++t)
        {
            appGeneral(_T("%2.12f%s"), m_lstAverageResults[ty][t], (t != (_HC_Lti - 2)) ? _T(",") : _T(""));
        }
        appGeneral(_T("}%s"), ty != 19 ? _T("\n}\n") : _T(",\n"));
    }
}

void CMeasureMesonCorrelatorStaggered::Reset()
{
    m_uiConfigurationCount = 0;
    m_lstResults.RemoveAll();
}




__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================