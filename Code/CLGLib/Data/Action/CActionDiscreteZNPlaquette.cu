//=============================================================================
// FILENAME : CActionDiscreteZNPlaquette.cu
//
// DESCRIPTION:
// Z_N Wilson plaquette action -- device functions and class implementation
//
// REVISION:
//  [05/16/2026 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionDiscreteZNPlaquette.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionDiscreteZNPlaquette)

template<INT N>
__device__ DOUBLE _deviceZNPlaquetteLinkActionWeight(const void* __restrict__ pGauge, const void* __restrict__ pStaple, UINT linkIndex, UINT k, DOUBLE beta)
{
    UN_USE(pGauge);
    const deviceZN<N>* stapleField = static_cast<const deviceZN<N>*>(pStaple);
    CLGComplex staple = stapleField[linkIndex].m_me;
    DOUBLE theta = 2.0 * PI * k / N;
    DOUBLE re = staple.x * _cos(theta) - staple.y * _sin(theta);
    return beta * re;
}

// __device__ __constant__ array of function pointers -- initialized at
// compile time with device addresses.  cudaMemcpyFromSymbol can copy from
// this array to host because the array itself is a proper device symbol.
// Pattern from CMeasureChargeAndCurrents.cu.
__device__ __constant__ DeviceDiscreteLinkActionFunc _cZNFuncs[5] =
{
    _deviceZNPlaquetteLinkActionWeight<2>,
    _deviceZNPlaquetteLinkActionWeight<3>,
    _deviceZNPlaquetteLinkActionWeight<4>,
    _deviceZNPlaquetteLinkActionWeight<5>,
    _deviceZNPlaquetteLinkActionWeight<6>
};

CActionDiscreteZNPlaquette::CActionDiscreteZNPlaquette()
    : CActionDiscreteGauge()
{
}

DeviceDiscreteLinkActionFunc CActionDiscreteZNPlaquette::GetDeviceFunc() const
{
    if (m_byGaugeFieldIds.Num() < 1)
    {
        appCrucial(_T("CActionDiscreteZNPlaquette: no gauge field configured!\n"));
        _FAIL_EXIT;
    }

    const CFieldGauge* pGauge = dynamic_cast<const CFieldGauge*>(m_pOwner->GetFieldById(GetPrimaryGaugeFieldId()));
    if (NULL == pGauge)
    {
        appCrucial(_T("CActionDiscreteZNPlaquette: failed to get gauge field!\n"));
        _FAIL_EXIT;
    }

    INT iIdx = -1;
    switch (pGauge->GetFieldType())
    {
    case EFT_GaugeZ2: iIdx = 0; break;
    case EFT_GaugeZ3: iIdx = 1; break;
    case EFT_GaugeZ4: iIdx = 2; break;
    case EFT_GaugeZ5: iIdx = 3; break;
    case EFT_GaugeZ6: iIdx = 4; break;
    default:
        appCrucial(_T("CActionDiscreteZNPlaquette: unsupported gauge field type!\n"));
        _FAIL_EXIT;
    }

    DeviceDiscreteLinkActionFunc hFuncs[5];
    checkCudaErrors(cudaMemcpyFromSymbol(hFuncs, _cZNFuncs, sizeof(DeviceDiscreteLinkActionFunc) * 5));
    return hFuncs[iIdx];
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
