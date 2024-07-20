//=============================================================================
// FILENAME : CGaugeFixingRandom.cpp
// 
// DESCRIPTION:
//
//
// REVISION:
//  [09/25/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3.h"
#include "Data/Field/Staggered/CFieldFermionKSSU3.h"
#include "CGaugeFixingRandom.h"

__BEGIN_NAMESPACE

#pragma region kernels

__global__ void _CLG_LAUNCH_BOUND
_kernelRandomGauge(deviceSU3* pGx)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    if (site.IsDirichlet())
    {
        pGx[uiSiteIndex] = deviceSU3::makeSU3Id();
    }
    else
    {
        pGx[uiSiteIndex] = deviceSU3::makeSU3Random(_deviceGetFatIndex(uiSiteIndex, 0));
    }
}

/**
 * g(n) U_mu(n) g(n+mu)^dagger
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformRandom(
    BYTE byFieldId,
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pGauge)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const deviceSU3 left(pGx[uiSiteIndex]);

    for (BYTE dir = 0; dir < _DC_Dir; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
            deviceSU3 res(pGauge[uiLinkDir]);
            SSmallInt4 sWalking = _deviceSmallInt4OffsetC(sSite4, dir + 1);
            const SIndex site_p_mu = __idx->m_pDeviceIndexPositionToSIndex[1][__idx->_deviceGetBigIndex(sWalking)];
            if (!site_p_mu.IsDirichlet())
            {
                res.MulDagger(pGx[site_p_mu.m_uiSiteIndex]);
            }

            pGauge[uiLinkDir] = left.MulC(res);
        }
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformFermionWilsonSU3(
    const deviceSU3* __restrict__ pGx,
    deviceWilsonVectorSU3* pFermion)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    if (!site.IsDirichlet())
    {
        pFermion[uiSiteIndex] = pGx[uiSiteIndex].MulWilsonVector(pFermion[uiSiteIndex]);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformFermionKSSU3(
    const deviceSU3* __restrict__ pGx,
    deviceSU3Vector* pFermion)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    if (!site.IsDirichlet())
    {
        pFermion[uiSiteIndex] = pGx[uiSiteIndex].MulVector(pFermion[uiSiteIndex]); 
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformAPhys(
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pAphys)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[1][uiBigIdx];

    if (!site.IsDirichlet())
    {
        for (BYTE dir = 0; dir < uiDir; ++dir)
        {
            const UINT uiLinkIndex = _deviceGetLinkIndex(uiSiteIndex, dir);
            pAphys[uiLinkIndex] = pGx[uiSiteIndex].MulC(pAphys[uiLinkIndex]);
            pAphys[uiLinkIndex].MulDagger(pGx[uiSiteIndex]);
        }
    }
}

#pragma endregion


#pragma endregion

__CLGIMPLEMENT_CLASS(CGaugeFixingRandom)

void CGaugeFixingRandom::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;
    checkCudaErrors(cudaMalloc((void**)& m_pG, _HC_Volume * sizeof(deviceSU3)));
}

void CGaugeFixingRandom::GaugeFixing(CFieldGauge* pResGauge)
{
    if (NULL == pResGauge || EFT_GaugeSU3 != pResGauge->GetFieldType())
    {
        appCrucial(_T("CGaugeFixingLandauCornell only implemented with gauge SU3!\n"));
        return;
    }
    CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pResGauge);

    preparethread;
    _kernelRandomGauge << <block, threads >> > (m_pG);
    _kernelGaugeTransformRandom << <block, threads >> > (pGaugeSU3->m_byFieldId, m_pG, pGaugeSU3->m_pDeviceData);
}

void CGaugeFixingRandom::AlsoFixingFermion(CFieldFermion* pFermion) const
{
    if (EFT_FermionWilsonSquareSU3 == pFermion->GetFieldType())
    {
        AlsoFixingFermionWilsonSU3(dynamic_cast<CFieldFermionWilsonSquareSU3*>(pFermion));
    }
    else if (EFT_FermionStaggeredSU3 == pFermion->GetFieldType())
    {
        AlsoFixingFermionKSSU3(dynamic_cast<CFieldFermionKSSU3*>(pFermion));
    }
}

void CGaugeFixingRandom::AlsoFixingFermionWilsonSU3(CFieldFermionWilsonSquareSU3* pFermion) const
{
    preparethread;
    _kernelGaugeTransformFermionWilsonSU3 << <block, threads >> > (m_pG, pFermion->m_pDeviceData);
}

void CGaugeFixingRandom::AlsoFixingFermionKSSU3(CFieldFermionKSSU3* pFermion) const
{
    preparethread;
    _kernelGaugeTransformFermionKSSU3 << <block, threads >> > (m_pG, pFermion->m_pDeviceData);
}

void CGaugeFixingRandom::AlsoFixingAphys(CFieldGauge* pGauge) const
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CGaugeFixingLandauCornell only implemented with gauge SU3!\n"));
        return;
    }
    CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pGauge);

    preparethread;
    _kernelGaugeTransformAPhys << <block, threads >> > (m_pG, pGaugeSU3->m_pDeviceData);
}

CCString CGaugeFixingRandom::GetInfos(const CCString& tab) const
{
    CCString sRet;
    sRet = sRet + tab + _T("Name : CGaugeFixingRandom\n");
    return sRet;
}

__END_NAMESPACE


//=============================================================================
// END OF FILE
//=============================================================================