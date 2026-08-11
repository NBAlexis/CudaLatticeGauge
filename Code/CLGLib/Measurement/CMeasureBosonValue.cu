//=============================================================================
// FILENAME : CMeasureBosonValue.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [mm/dd/yy]
//  [11/07/2024 nbale]
//=============================================================================

#include "CLGLib_Private.h"
#include "Tools/Math/DeviceInlineTemplate.h"
#include "Data/Field/Boson/CFieldBosonVN.h"
#include "CMeasureBosonValue.h"

__BEGIN_NAMESPACE

#pragma region kernels

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelBosonValueComplex(
    const deviceDataBoson* __restrict__ pDeviceData,
    cuDoubleComplex* pResultData)
{
    intokernalInt4;

    const UINT uiXY = sSite4.x * _DC_Ly + sSite4.y;
    const UINT uiZ = sSite4.z;
    const UINT uiT = sSite4.w;
    //const UINT uiSiteIndex = uiXY * _DC_GridDimZT + uiZ * _DC_Lt + uiT;
    const UINT uiWriteIdx = uiZ * _DC_Lx * _DC_Ly * _DC_Lt + uiXY * _DC_Lt + uiT;

    pResultData[uiWriteIdx] = _cToDouble(_detv(pDeviceData[uiSiteIndex]));
}

template<typename deviceDataBoson>
__global__ void _CLG_LAUNCH_BOUND
_kernelBosonValueReal(
    const deviceDataBoson* __restrict__ pDeviceData,
    DOUBLE* pResultData)
{
    intokernalInt4;

    const UINT uiXY = sSite4.x * _DC_Ly + sSite4.y;
    const UINT uiZ = sSite4.z;
    const UINT uiT = sSite4.w;
    //const UINT uiSiteIndex = uiXY * _DC_GridDimZT + uiZ * _DC_Lt + uiT;
    const UINT uiWriteIdx = uiZ * _DC_Lx * _DC_Ly * _DC_Lt + uiXY * _DC_Lt + uiT;

    pResultData[uiWriteIdx] = static_cast<DOUBLE>(_absv(pDeviceData[uiSiteIndex]));
}

template<>
__global__ void _CLG_LAUNCH_BOUND
_kernelBosonValueReal<Real>(
    const Real* __restrict__ pDeviceData,
    DOUBLE* pResultData)
{
    intokernalInt4;

    const UINT uiXY = sSite4.x * _DC_Ly + sSite4.y;
    const UINT uiZ = sSite4.z;
    const UINT uiT = sSite4.w;
    //const UINT uiSiteIndex = uiXY * _DC_GridDimZT + uiZ * _DC_Lt + uiT;
    const UINT uiWriteIdx = uiZ * _DC_Lx * _DC_Ly * _DC_Lt + uiXY * _DC_Lt + uiT;

    pResultData[uiWriteIdx] = static_cast<DOUBLE>(pDeviceData[uiSiteIndex]);
}

#pragma endregion

template<typename deviceDataBoson, typename deviceDataGauge>
void TMeasureBosonValue<deviceDataBoson, deviceDataGauge>::OnConfigurationAccepted(INT gaugeNum, INT bosonNum, INT tensor2Num, const class CFieldGauge* const* pAcceptGauge, const class CFieldBoson* const* pAcceptBoson, const class CFieldTensor2* const* tensor2Fields, const class CFieldGauge* const* pCorrespondingStaple)
{
    if (m_lstBosonFieldIds.Num() < 1)
    {
        appCrucial(_T("TMeasureBosonValue, BosonFields not properly set\n"));
        return;
    }
    INT ibosonidx = CLatticeData::GetBosonFieldIndexById(bosonNum, pAcceptBoson, m_lstBosonFieldIds[0]);
    if (ibosonidx < 0 || ibosonidx >= bosonNum)
    {
        appCrucial(_T("TMeasureBosonValue, BosonFields not properly set\n"));
        return;
    }

    const CFieldBosonVN<deviceDataBoson, deviceDataGauge>* bosonfield = dynamic_cast<const CFieldBosonVN<deviceDataBoson, deviceDataGauge>*>(pAcceptBoson[ibosonidx]);
    if (NULL == bosonfield)
    {
        appCrucial(_T("TMeasureBosonValue, Only work with CFieldBosonVN\n"));
        return;
    }

    const deviceDataBoson* pBuffer = bosonfield->m_pDeviceData;

    const UINT xyt = _HC_Lx * _HC_Ly * _HC_Lt;
    preparethread;
    _LAUNCH_KERNEL(_kernelBosonValueComplex<deviceDataBoson>, block, threads, pBuffer, _D_ComplexThreadBuffer);
    _LAUNCH_KERNEL(_kernelBosonValueReal<deviceDataBoson>, block, threads, pBuffer, _D_RealThreadBuffer);

    //reduce sum will mass up the buffer
    cuDoubleComplex sumc = appGetCudaHelper()->ThreadBufferSum(_D_ComplexThreadBuffer);
    sumc.x = sumc.x / _HC_Volume;
    sumc.y = sumc.y / _HC_Volume;
    DOUBLE sumd = appGetCudaHelper()->ThreadBufferSum(_D_RealThreadBuffer);
    sumd = sumd / _HC_Volume;

    m_lstEveryConfigurationC.AddItem(sumc);
    m_lstEveryConfigurationR.AddItem(sumd);

    _LAUNCH_KERNEL(_kernelBosonValueComplex<deviceDataBoson>, block, threads, pBuffer, _D_ComplexThreadBuffer);
    _LAUNCH_KERNEL(_kernelBosonValueReal<deviceDataBoson>, block, threads, pBuffer, _D_RealThreadBuffer);

    TArray<cuDoubleComplex> zsliceC;
    TArray<DOUBLE> zsliceR;
    for (UINT z = 0; z < _HC_Lz; ++z)
    {
        cuDoubleComplex sumc_zslice = CCudaHelper::ReduceComplex(_D_ComplexThreadBuffer + xyt * z, xyt);
        DOUBLE sumd_zslice = CCudaHelper::ReduceReal(_D_RealThreadBuffer + xyt * z, xyt);

        sumc_zslice.x = sumc_zslice.x / xyt;
        sumc_zslice.y = sumc_zslice.y / xyt;
        sumd_zslice = sumd_zslice / xyt;

        zsliceC.AddItem(sumc_zslice);
        zsliceR.AddItem(sumd_zslice);
    }
    m_lstEveryConfigurationZsliceC.AddItem(zsliceC);
    m_lstEveryConfigurationZsliceR.AddItem(zsliceR);

    if (NULL != m_pOwner)
    {
        m_pOwner->AddOneConfigurationResult(this, _T("BosonValueComplex"), _cToRealC(sumc));
        m_pOwner->AddOneConfigurationResult(this, _T("BosonValueReal"), sumd);

        TArray<CLGComplex> zsliceCClg;
        for (INT i = 0; i < zsliceC.Num(); ++i)
        {
            zsliceCClg.AddItem(_cToRealC(zsliceC[i]));
        }
        m_pOwner->AddOneConfigurationResult(this, _T("BosonValueZSliceComplex"), zsliceCClg);
        m_pOwner->AddOneConfigurationResult(this, _T("BosonValueZSliceReal"), zsliceR);
    }

    ++m_uiConfigurationCount;
}

template<typename deviceDataBoson, typename deviceDataGauge>
void TMeasureBosonValue<deviceDataBoson, deviceDataGauge>::Reset()
{
    CMeasure::Reset();
    m_lstEveryConfigurationC.RemoveAll();
    m_lstEveryConfigurationR.RemoveAll();
    m_lstEveryConfigurationZsliceC.RemoveAll();
    m_lstEveryConfigurationZsliceR.RemoveAll();
}

__CLGIMPLEMENT_CLASS(CMeasureBosonValueReal)
__CLGIMPLEMENT_CLASS(CMeasureBosonValueU1)
__CLGIMPLEMENT_CLASS(CMeasureBosonValueSU2)
__CLGIMPLEMENT_CLASS(CMeasureBosonValueSU3)

#if _CLG_SU4_BOSON
__CLGIMPLEMENT_CLASS(CMeasureBosonValueSU4)
#endif
#if _CLG_SU5_BOSON
__CLGIMPLEMENT_CLASS(CMeasureBosonValueSU5)
#endif
#if _CLG_SU6_BOSON
__CLGIMPLEMENT_CLASS(CMeasureBosonValueSU6)
#endif
#if _CLG_SU7_BOSON
__CLGIMPLEMENT_CLASS(CMeasureBosonValueSU7)
#endif
#if _CLG_SU8_BOSON
__CLGIMPLEMENT_CLASS(CMeasureBosonValueSU8)
#endif

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================