//=============================================================================
// FILENAME : CActionDiscreteDNPlaquette.cu
//
// DESCRIPTION:
// D_N Wilson plaquette action -- device functions and class implementation.
// Action weight: W(V_k) = -beta * Re[Tr(V_k^dag * S)]
// where V_k iterates over the 2N elements of D_N (k = 0..2N-1).
//
// REVISION:
//  [05/16/2026 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "CActionDiscreteDNPlaquette.h"

__BEGIN_NAMESPACE

__CLGIMPLEMENT_CLASS(CActionDiscreteDNPlaquette)

template<INT N>
__device__ DOUBLE _deviceDNPlaquetteLinkActionWeight(const void* __restrict__ pGauge, const void* __restrict__ pStaple, UINT linkIndex, UINT k, DOUBLE beta)
{
    UN_USE(pGauge);
    const deviceDN<N>* stapleField = static_cast<const deviceDN<N>*>(pStaple);
    const deviceDN<N>& staple = stapleField[linkIndex];

    DOUBLE theta = 2.0 * PI / N;
    INT parity = k / N;
    INT kk = k % N;
    DOUBLE thk = theta * kk;
    DOUBLE cosk = cos(thk);
    DOUBLE sink = sin(thk);
    DOUBLE trRe;
    if (0 == parity)
    {
        // Rotation r^{kk}: diag(e^{i*theta*kk}, e^{-i*theta*kk})
        trRe = cosk * (staple.m_me[0].x + staple.m_me[3].x)
             + sink * (staple.m_me[3].y - staple.m_me[0].y);
    }
    else
    {
        // Reflection r^{kk} s: anti-diag(e^{i*theta*kk}, e^{-i*theta*kk})
        trRe = cosk * (staple.m_me[2].x + staple.m_me[1].x)
             + sink * (staple.m_me[1].y - staple.m_me[2].y);
    }
    return beta * trRe;
}

__device__ __constant__ DeviceDiscreteLinkActionFunc _cDNFuncs[3] =
{
    _deviceDNPlaquetteLinkActionWeight<3>,
    _deviceDNPlaquetteLinkActionWeight<4>,
    _deviceDNPlaquetteLinkActionWeight<8>
};

CActionDiscreteDNPlaquette::CActionDiscreteDNPlaquette()
    : CActionDiscreteGauge()
{
}

DeviceDiscreteLinkActionFunc CActionDiscreteDNPlaquette::GetDeviceFunc() const
{
    if (m_byGaugeFieldIds.Num() < 1)
    {
        appCrucial(_T("CActionDiscreteDNPlaquette: no gauge field configured!\n"));
        _FAIL_EXIT;
    }

    const CFieldGauge* pGauge = dynamic_cast<const CFieldGauge*>(m_pOwner->GetFieldById(GetPrimaryGaugeFieldId()));
    if (NULL == pGauge)
    {
        appCrucial(_T("CActionDiscreteDNPlaquette: failed to get gauge field!\n"));
        _FAIL_EXIT;
    }

    INT iIdx = -1;
    switch (pGauge->GetFieldType())
    {
#if _CLG_D3_GAUGE
    case EFT_GaugeD3: iIdx = 0; break;
#endif
#if _CLG_D4_GAUGE
    case EFT_GaugeD4: iIdx = 1; break;
#endif
#if _CLG_D8_GAUGE
    case EFT_GaugeD8: iIdx = 2; break;
#endif
    default:
        appCrucial(_T("CActionDiscreteDNPlaquette: unsupported gauge field type!\n"));
        _FAIL_EXIT;
    }

    DeviceDiscreteLinkActionFunc hFuncs[3];
    checkCudaErrors(cudaMemcpyFromSymbol(hFuncs, _cDNFuncs, sizeof(DeviceDiscreteLinkActionFunc) * 3));
    return hFuncs[iIdx];
}

__END_NAMESPACE

//=============================================================================
// END OF FILE
//=============================================================================
