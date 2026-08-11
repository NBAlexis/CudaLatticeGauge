//=============================================================================
// FILENAME : CGaugeFixingRandom.cpp
//
// DESCRIPTION:
// Random gauge transform for testing gauge invariance.
// Supports SU(2) and SU(3) gauge fields.
//
// REVISION:
//  [09/25/2019 nbale]
//=============================================================================
#include "CLGLib_Private.h"
#include "Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3.h"
#include "Data/Field/Staggered/CFieldFermionKST.h"
#include "CGaugeFixingRandom.h"

__BEGIN_NAMESPACE

#pragma region SU2 kernels

__global__ void _CLG_LAUNCH_BOUND
_kernelRandomGaugeSU2(deviceSU2* pGx, BYTE byFieldId)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (site.IsDirichlet())
    {
        pGx[uiSiteIndex] = deviceSU2::makeSU2Id();
    }
    else
    {
        pGx[uiSiteIndex] = deviceSU2::makeSU2Random(_deviceGetLinkIndex(uiSiteIndex, 0));
    }
}

/**
 * g(n) U_mu(n) g(n+mu)^dagger
 */
__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformRandomSU2(
    BYTE byFieldId,
    const deviceSU2* __restrict__ pGx,
    deviceSU2* pGauge)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const deviceSU2 left(pGx[uiSiteIndex]);

    for (BYTE dir = 0; dir < _DC_Dir; ++dir)
    {
        if (!__idx->_deviceIsBondOnSurface(uiBigIdx, byFieldId, dir))
        {
            UINT uiLinkDir = _deviceGetLinkIndex(uiSiteIndex, dir);
            deviceSU2 res(pGauge[uiLinkDir]);
            SSmallInt4 sWalking = _deviceSmallInt4OffsetC(sSite4, dir + 1);
            const SIndex site_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__idx->_deviceGetBigIndex(sWalking)];
            if (!site_p_mu.IsDirichlet())
            {
                res.MulDagger(pGx[site_p_mu.m_uiSiteIndex]);
            }

            pGauge[uiLinkDir] = left.MulC(res);
        }
    }
}

#pragma endregion

#pragma region SU3 kernels

__global__ void _CLG_LAUNCH_BOUND
_kernelRandomGauge(deviceSU3* pGx, BYTE byFieldId)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

    if (site.IsDirichlet())
    {
        pGx[uiSiteIndex] = deviceSU3::makeSU3Id();
    }
    else
    {
        pGx[uiSiteIndex] = deviceSU3::makeSU3Random(_deviceGetLinkIndex(uiSiteIndex, 0));
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
            const SIndex site_p_mu = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][__idx->_deviceGetBigIndex(sWalking)];
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
    deviceWilsonVectorSU3* pFermion,
    BYTE byFermionFieldId)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFermionFieldId][uiBigIdx];

    if (!site.IsDirichlet())
    {
        pFermion[uiSiteIndex] = pGx[uiSiteIndex].MulWilsonVector(pFermion[uiSiteIndex]);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformFermionKSSU3(
    const deviceSU3* __restrict__ pGx,
    deviceSU3Vector* pFermion,
    BYTE byFermionFieldId)
{
    intokernalInt4;

    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFermionFieldId][uiBigIdx];

    if (!site.IsDirichlet())
    {
        pFermion[uiSiteIndex] = pGx[uiSiteIndex].MulVector(pFermion[uiSiteIndex]);
    }
}

__global__ void _CLG_LAUNCH_BOUND
_kernelGaugeTransformAPhys(
    const deviceSU3* __restrict__ pGx,
    deviceSU3* pAphys,
    BYTE byFieldId)
{
    intokernalInt4;

    const BYTE uiDir = static_cast<BYTE>(_DC_Dir);
    const UINT uiBigIdx = __idx->_deviceGetBigIndex(sSite4);
    const SIndex site = __idx->m_pDeviceIndexPositionToSIndex[byFieldId][uiBigIdx];

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

__CLGIMPLEMENT_CLASS(CGaugeFixingRandom)

void CGaugeFixingRandom::Initial(class CLatticeData* pOwner, const CParameters& params)
{
    m_pOwner = pOwner;
    checkCudaErrors(__cudaMalloc((void**)& m_pGSU2, _HC_Volume * sizeof(deviceSU2)));
    checkCudaErrors(__cudaMalloc((void**)& m_pG, _HC_Volume * sizeof(deviceSU3)));
}

void CGaugeFixingRandom::GaugeFixing(CFieldGauge* pResGauge)
{
    if (NULL == pResGauge)
    {
        appCrucial(_T("CGaugeFixingRandom::GaugeFixing: pResGauge is NULL!\n"));
        return;
    }

#if _CLG_MULTI_GPU
    if (NULL != appGetComm() && appGetComm()->Size() > 1)
    {
        //P4-2.5 fix: gather -> rank0 global-context transform -> scatter, same as
        //the other fixers. The random transform draws from _deviceRandomF, whose
        //curand states are sized to the LOCAL lattice volume (Random.cu); under
        //the GLOBAL context the site index would reach the global volume, so the
        //states are rebuilt via CRandom::EnterGlobalContext() (global volume,
        //globally-consistent reseed) for the transform and restored afterwards.
        const EFieldType eType = pResGauge->GetFieldType();
        CFieldGaugeSU2* pGaugeSU2 = (EFT_GaugeSU2 == eType) ? dynamic_cast<CFieldGaugeSU2*>(pResGauge) : NULL;
        CFieldGaugeSU3* pGaugeSU3 = (EFT_GaugeSU3 == eType) ? dynamic_cast<CFieldGaugeSU3*>(pResGauge) : NULL;
        if (NULL == pGaugeSU2 && NULL == pGaugeSU3)
        {
            appCrucial(_T("CGaugeFixingRandom: unsupported field type %d!\n"), pResGauge->GetFieldType());
            return;
        }

        const UINT uiElemSize = (NULL != pGaugeSU2) ? static_cast<UINT>(sizeof(deviceSU2)) : static_cast<UINT>(sizeof(deviceSU3));
        const UINT uiBytesPerSite = uiElemSize * static_cast<UINT>(_HC_Dir);
        const UINT uiLocalBytes = uiElemSize * static_cast<UINT>(_HC_LinkCount);
        const BYTE* pLocalData = (NULL != pGaugeSU2) ? (const BYTE*)pGaugeSU2->m_pDeviceData : (const BYTE*)pGaugeSU3->m_pDeviceData;

        BYTE* pHostLocal = (BYTE*)malloc(uiLocalBytes);
        checkCudaErrors(cudaMemcpy(pHostLocal, pLocalData, uiLocalBytes, cudaMemcpyDeviceToHost));

        UINT uiGlobalBytes = 0;
        BYTE* pGlobal = appGetComm()->GatherFieldToRoot(pHostLocal, uiBytesPerSite, uiGlobalBytes);
        free(pHostLocal);

        if (appGetComm()->IsRoot())
        {
            BYTE* pDevGlobal = NULL;
            checkCudaErrors(__cudaMalloc((void**)&pDevGlobal, uiGlobalBytes));
            checkCudaErrors(cudaMemcpy(pDevGlobal, pGlobal, uiGlobalBytes, cudaMemcpyHostToDevice));
            free(pGlobal);
            pGlobal = NULL;

            MGEnterGlobalFixerContext();
            ResizeBuffersToGlobal();

            CRandom* pRNG = (NULL != appGetLattice()) ? appGetLattice()->m_pRandom : NULL;
            if (NULL != pRNG)
            {
                pRNG->EnterGlobalContext();
            }

            preparethread;
            if (NULL != pGaugeSU2)
            {
                _LAUNCH_KERNEL(_kernelRandomGaugeSU2, block, threads, m_pGSU2, pGaugeSU2->m_byFieldId);
                _LAUNCH_KERNEL(_kernelGaugeTransformRandomSU2, block, threads, pGaugeSU2->m_byFieldId, m_pGSU2, reinterpret_cast<deviceSU2*>(pDevGlobal));
            }
            else
            {
                _LAUNCH_KERNEL(_kernelRandomGauge, block, threads, m_pG, pGaugeSU3->m_byFieldId);
                _LAUNCH_KERNEL(_kernelGaugeTransformRandom, block, threads, pGaugeSU3->m_byFieldId, m_pG, reinterpret_cast<deviceSU3*>(pDevGlobal));
            }

            if (NULL != pRNG)
            {
                pRNG->ExitGlobalContext();
            }
            RestoreLocalBuffers();
            MGExitGlobalFixerContext();

            pGlobal = (BYTE*)malloc(uiGlobalBytes);
            checkCudaErrors(cudaMemcpy(pGlobal, pDevGlobal, uiGlobalBytes, cudaMemcpyDeviceToHost));
            checkCudaErrors(__cudaFree(pDevGlobal));
        }

        BYTE* pLocalOut = (BYTE*)malloc(uiLocalBytes);
        appGetComm()->ScatterFieldFromRoot(pGlobal, uiBytesPerSite, pLocalOut);
        checkCudaErrors(cudaMemcpy((void*)(NULL != pGaugeSU2 ? (void*)pGaugeSU2->m_pDeviceData : (void*)pGaugeSU3->m_pDeviceData), pLocalOut, uiLocalBytes, cudaMemcpyHostToDevice));
        free(pLocalOut);
        if (NULL != pGlobal)
        {
            free(pGlobal);
        }

        //Every local link was rewritten by the scatter memcpy, which bypasses
        //guarded launches; bump the owner handle so the next reader's guard
        //refills the halo (no-op for an unbound handle).
        pResGauge->NotifyWritten();
        return;
    }
#endif

    preparethread;

    switch (pResGauge->GetFieldType())
    {
    case EFT_GaugeSU2:
    {
        CFieldGaugeSU2* pGaugeSU2 = dynamic_cast<CFieldGaugeSU2*>(pResGauge);
        _LAUNCH_KERNEL(_kernelRandomGaugeSU2, block, threads, m_pGSU2, pGaugeSU2->m_byFieldId);
        _LAUNCH_KERNEL(_kernelGaugeTransformRandomSU2, block, threads, pGaugeSU2->m_byFieldId, m_pGSU2, pGaugeSU2->m_pDeviceData);
        break;
    }
    case EFT_GaugeSU3:
    {
        CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pResGauge);
        _LAUNCH_KERNEL(_kernelRandomGauge, block, threads, m_pG, pGaugeSU3->m_byFieldId);
        _LAUNCH_KERNEL(_kernelGaugeTransformRandom, block, threads, pGaugeSU3->m_byFieldId, m_pG, pGaugeSU3->m_pDeviceData);
        break;
    }
    default:
        appCrucial(_T("CGaugeFixingRandom: unsupported field type %d!\n"), pResGauge->GetFieldType());
        break;
    }
}

#if _CLG_MULTI_GPU
void CGaugeFixingRandom::ResizeBuffersToGlobal()
{
    //Called under the temporary GLOBAL lattice context (rank 0 only), where
    //_HC_Volume is the GLOBAL value: save the local pointers and re-allocate the
    //transform buffers to the global volume.
    m_pSavedGSU2 = m_pGSU2;
    m_pSavedG = m_pG;
    const UINT uiVol = _HC_Volume;
    checkCudaErrors(__cudaMalloc((void**)&m_pGSU2, uiVol * sizeof(deviceSU2)));
    checkCudaErrors(__cudaMalloc((void**)&m_pG, uiVol * sizeof(deviceSU3)));
}

void CGaugeFixingRandom::RestoreLocalBuffers()
{
    //Free the temporary global buffers and put the local ones back. The
    //destructor frees whatever is current, so it must always see the local
    //pointers afterwards.
    checkCudaErrors(__cudaFree(m_pGSU2));
    checkCudaErrors(__cudaFree(m_pG));
    m_pGSU2 = m_pSavedGSU2;
    m_pG = m_pSavedG;
    m_pSavedGSU2 = NULL;
    m_pSavedG = NULL;
}
#endif

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
    _LAUNCH_KERNEL(_kernelGaugeTransformFermionWilsonSU3, block, threads, m_pG, pFermion->m_pDeviceData, pFermion->m_byFieldId);
}

void CGaugeFixingRandom::AlsoFixingFermionKSSU3(CFieldFermionKSSU3* pFermion) const
{
    preparethread;
    _LAUNCH_KERNEL(_kernelGaugeTransformFermionKSSU3, block, threads, m_pG, pFermion->m_pDeviceData, pFermion->m_byFieldId);
}

void CGaugeFixingRandom::AlsoFixingAphys(CFieldGauge* pGauge) const
{
    if (NULL == pGauge || EFT_GaugeSU3 != pGauge->GetFieldType())
    {
        appCrucial(_T("CGaugeFixingRandom::AlsoFixingAphys only implemented with gauge SU3!\n"));
        return;
    }
    CFieldGaugeSU3* pGaugeSU3 = dynamic_cast<CFieldGaugeSU3*>(pGauge);

    preparethread;
    _LAUNCH_KERNEL(_kernelGaugeTransformAPhys, block, threads, m_pG, pGaugeSU3->m_pDeviceData, pGaugeSU3->m_byFieldId);
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
