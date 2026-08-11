//=============================================================================
// FILENAME : TestImprove6.cpp
//
// DESCRIPTION:
// Improve-1 I6 (multi-GPU-improve1.md section 3.4, work item "在 storage
// owner 层注册所有受管理 buffer"): every managed device buffer owned OUTSIDE
// a field class is registered with the halo registry by its storage owner, so
// the I5 automatic launch guard discovers it through the ORIGINAL public
// launch macros alone -- no MG parameter and no new macro at any call site.
//
// Covered owners (appendix A.5 dispositions):
//   - CFieldFermionWilsonSquareCloverBuffer: the PhidPhi singleton is
//     HaloCapable (site-major after I6, per-site halo tail appended);
//   - CStapleCacheT: staple/plaquette device pointer arrays are buffer sets
//     over the pooled member fields; the flat Fmunu / rotation-link buffers
//     are LocalOnly (consumed strictly site-locally);
//   - CRationalFieldPointer: the rational field pointer array is a buffer set
//     re-bound to the current pooled membership after every fill;
//   - CFieldGaugeU1Real::m_pFmunu: lazy LocalOnly extent.
//
// REVISION:
//  [08/09/2026 Improve-1 I6 nbale]
//=============================================================================

#include "CLGTest.h"

#include "TestImproveHaloCommon.h"

/**
 * I6 gate 1: the clover PhidPhi singleton is a registered HaloCapable extent,
 * discovered by pointer-range lookup (base, interior, halo tail), re-bound on
 * reallocation, Ensured by the pre-launch guard, and its per-site halo refill
 * is content-correct against the deterministic global-site pattern (this also
 * proves the I6 site-major relayout of the two clover kernels end to end).
 */
UINT TestMGOwnerCloverBuffer(CParameters&)
{
    UINT uiErrors = 0;
    const UINT uiVolume = _HC_Volume;
    const UINT uiHaloSites = static_cast<UINT>(_HC_HaloSiteCount());
    const UINT uiBytesPerSite = 6U * static_cast<UINT>(sizeof(deviceSU3));
    const UINT uiElemPerSite = uiBytesPerSite / static_cast<UINT>(sizeof(Real));

    //--- Part 1: registration metadata of the singleton extent ---
    CFieldFermionWilsonSquareCloverBuffer* pSingleton = CFieldFermionWilsonSquareCloverBuffer::GetInstance();
    deviceSU3* pBuf = pSingleton->GetPhidPhiBuffer<deviceSU3>(uiVolume * static_cast<ULONGLONG>(uiBytesPerSite));
    CHaloBufferHandle* pHandle = pSingleton->GetHaloBufferHandle();
    if (NULL == pBuf || !pHandle->IsBound())
    {
        appGeneral(_T("clover buffer: allocation or handle bind missing.\n"));
        ++uiErrors;
    }
    else
    {
        const SHaloBufferInfo& info = pHandle->Info();
        if (!info.m_bHaloCapable)
        {
            appGeneral(_T("clover buffer: must be HaloCapable.\n"));
            ++uiErrors;
        }
        if (info.m_uiBytesPerSite != uiBytesPerSite || info.m_uiLocalSiteCount != uiVolume
            || info.m_uiHaloSiteCount != uiHaloSites
            || info.m_uiCapacityBytes != static_cast<size_t>(uiVolume + uiHaloSites) * uiBytesPerSite)
        {
            appGeneral(_T("clover buffer: extent mismatch (bps %u/%u, sites %u/%u+%u, capacity %llu).\n"),
                info.m_uiBytesPerSite, uiBytesPerSite, info.m_uiLocalSiteCount, uiVolume, uiHaloSites,
                static_cast<ULONGLONG>(info.m_uiCapacityBytes));
            ++uiErrors;
        }
        if (CLG_PseudoFieldIdCloverPhidPhi != info.m_byFieldId)
        {
            appGeneral(_T("clover buffer: pseudo tag id %u != %u.\n"), info.m_byFieldId, CLG_PseudoFieldIdCloverPhidPhi);
            ++uiErrors;
        }
        //Range lookup: base, one-site-in, and the LAST byte of the halo tail.
        BYTE* pBytes = reinterpret_cast<BYTE*>(pBuf);
        if (appGetHaloManager()->LookupByPointer(pBytes) != pHandle
            || appGetHaloManager()->LookupByPointer(pBytes + uiBytesPerSite) != pHandle
            || appGetHaloManager()->LookupByPointer(pBytes + info.m_uiCapacityBytes - 1) != pHandle)
        {
            appGeneral(_T("clover buffer: base/interior/halo-tail lookup miss.\n"));
            ++uiErrors;
        }
    }

    //--- Part 2: reallocation re-binds (LOCAL instance, never the singleton) ---
    {
        const UINT uiCount0 = appGetHaloManager()->RegisteredHandleCount();
        {
            CFieldFermionWilsonSquareCloverBuffer localBuf;
            localBuf.GetPhidPhiBuffer<deviceSU3>(uiVolume * static_cast<ULONGLONG>(uiBytesPerSite));
            if (!localBuf.GetHaloBufferHandle()->IsBound()
                || appGetHaloManager()->RegisteredHandleCount() != uiCount0 + 1)
            {
                appGeneral(_T("clover local: first bind failed (count %u -> %u).\n"),
                    uiCount0, appGetHaloManager()->RegisteredHandleCount());
                ++uiErrors;
            }
            //Grow by exactly one site: forces free + alloc + re-Bind.
            localBuf.GetPhidPhiBuffer<deviceSU3>((uiVolume + 1) * static_cast<ULONGLONG>(uiBytesPerSite));
            if (!localBuf.GetHaloBufferHandle()->IsBound()
                || appGetHaloManager()->RegisteredHandleCount() != uiCount0 + 1)
            {
                appGeneral(_T("clover local: re-bind after grow failed (count %u).\n"),
                    appGetHaloManager()->RegisteredHandleCount());
                ++uiErrors;
            }
        }
        if (appGetHaloManager()->RegisteredHandleCount() != uiCount0)
        {
            appGeneral(_T("clover local: destructor did not unbind (count %u != %u).\n"),
                appGetHaloManager()->RegisteredHandleCount(), uiCount0);
            ++uiErrors;
        }
    }

    //--- Part 3: guard integration through the original public argument ---
    pHandle->NotifyWritten();
    const ULONGLONG ullV0 = pHandle->InteriorVersion();
    CHaloLaunchGuard guard = appGetHaloManager()->BeginFromArguments(pBuf, 42, static_cast<const BYTE*>(NULL));
    if (1 != guard.RecordedCount())
    {
        appGeneral(_T("clover guard: RecordedCount %d != 1.\n"), guard.RecordedCount());
        ++uiErrors;
    }
#if _CLG_MULTI_GPU
    if (uiHaloSites > 0)
    {
        if (pHandle->HaloVersion() != ullV0 || static_cast<UINT>(_HC_HaloWidth) != pHandle->ValidWidth())
        {
            appGeneral(_T("clover guard: pre-launch Ensure did not sync the halo.\n"));
            ++uiErrors;
        }
    }
#endif
    guard.Commit();
    guard.Commit();
    if (pHandle->HaloVersion() == pHandle->InteriorVersion())
    {
        appGeneral(_T("clover guard: Commit did not invalidate the halo.\n"));
        ++uiErrors;
    }

    //--- Part 4: halo content (site-major stride end to end, decomposed only) ---
    SImprove1HaloCtx ctx;
    if (Improve1BuildCtx(ctx))
    {
        Improve1FillLocalDeterministic(pBuf, uiBytesPerSite, uiElemPerSite, ctx, F(0.0));
        //The host fill is invisible to the guard; mark and Ensure via the handle.
        pHandle->NotifyWritten();
        appGetHaloManager()->Ensure(*pHandle, static_cast<UINT>(_HC_HaloWidth));
        uiErrors += Improve1CheckHalo(pBuf, uiBytesPerSite, uiElemPerSite, ctx, F(0.0), _T("CloverPhidPhi"));
        //Leave a dirty state so the next real consumer refills its own data.
        pHandle->NotifyWritten();
    }
    return uiErrors;
}
___REGIST_TEST(TestMGOwnerCloverBuffer, MG, TestMGOwnerCloverBuffer, MGOwnerCloverBuffer, _TEST_MULTIGPU);

/**
 * I6 gate 2: the staple cache registers its two device pointer arrays as
 * buffer sets over the pooled member gauge fields, and its flat Fmunu /
 * rotation-link buffers as LocalOnly extents.
 */
UINT TestMGOwnerStapleCacheSets(CParameters&)
{
    UINT uiErrors = 0;
    CField* pGauge = appGetLattice()->m_pGaugeField[0];
    CStapleCacheSU3* pCache = new CStapleCacheSU3();
    pCache->m_bCacheStaple = TRUE;
    pCache->m_bCachePlaqutte = TRUE;
    pCache->m_bCacheFmunu = TRUE;
    pCache->m_bCacheRotationLink = TRUE;
    pCache->InitialBuffers(pGauge->m_byFieldId);

    const CHaloBufferSetHandle* pStapleSet = appGetHaloManager()->LookupSetByPointer(pCache->m_pDeviceStaplePtr);
    const CHaloBufferSetHandle* pPlaqSet = appGetHaloManager()->LookupSetByPointer(pCache->m_pDevicePlaqPtr);
    if (NULL == pStapleSet || 6 != pStapleSet->MemberCount())
    {
        appGeneral(_T("staple cache: staple array set missing / %d members != 6.\n"),
            (NULL == pStapleSet) ? -1 : static_cast<INT>(pStapleSet->MemberCount()));
        ++uiErrors;
    }
    if (NULL == pPlaqSet || 6 != pPlaqSet->MemberCount())
    {
        appGeneral(_T("staple cache: plaquette array set missing / %d members != 6.\n"),
            (NULL == pPlaqSet) ? -1 : static_cast<INT>(pPlaqSet->MemberCount()));
        ++uiErrors;
    }
    if (NULL != pStapleSet)
    {
        for (UINT i = 0; i < pStapleSet->MemberCount(); ++i)
        {
            if (NULL == pStapleSet->Member(i) || !pStapleSet->Member(i)->IsBound())
            {
                appGeneral(_T("staple cache: member %u unbound.\n"), i);
                ++uiErrors;
            }
        }
    }

    //The flat buffers are LocalOnly; interior pointers resolve to the same handle.
    CHaloBufferHandle* pFmunu = appGetHaloManager()->LookupByPointer(pCache->m_pDeviceFmunuPtr);
    if (NULL == pFmunu || pFmunu->Info().m_bHaloCapable || 0 != pFmunu->Info().m_uiHaloSiteCount)
    {
        appGeneral(_T("staple cache: Fmunu buffer is not a registered LocalOnly extent.\n"));
        ++uiErrors;
    }
    CHaloBufferHandle* pRotation = appGetHaloManager()->LookupByPointer(pCache->m_pDeviceGaugeRotationLinksPtr);
    if (NULL == pRotation || pRotation->Info().m_bHaloCapable)
    {
        appGeneral(_T("staple cache: rotation-link buffer is not a registered LocalOnly extent.\n"));
        ++uiErrors;
    }
    else if (appGetHaloManager()->LookupByPointer(pCache->m_pDeviceGaugeRotationLinksPtr + _HC_Volume * 8) != pRotation)
    {
        //The cached-gauge consumers pass pCachedGauge + V * 8 for the tau half.
        appGeneral(_T("staple cache: rotation-link interior (tau half) lookup miss.\n"));
        ++uiErrors;
    }

    //Guard expansion: the array argument alone records ALL 6 members.
    CHaloLaunchGuard guard = appGetHaloManager()->BeginFromArguments(pCache->m_pDeviceStaplePtr);
    if (6 != guard.RecordedCount())
    {
        appGeneral(_T("staple cache: guard recorded %d handles != 6 members.\n"), guard.RecordedCount());
        ++uiErrors;
    }
    guard.Commit();

    //Save the outer extent address before the owner dies; after destruction
    //the set registry must no longer resolve it.
    const void* pStapleArray = pCache->m_pDeviceStaplePtr;
    appSafeDelete(pCache);
    if (NULL != appGetHaloManager()->LookupSetByPointer(pStapleArray))
    {
        appGeneral(_T("staple cache: set survives its owner.\n"));
        ++uiErrors;
    }
    return uiErrors;
}
___REGIST_TEST(TestMGOwnerStapleCacheSets, MG, TestMGOwnerStapleCacheSets, MGOwnerStapleCacheSets, _TEST_MULTIGPU);

/**
 * I6 gate 3: the rational field pointer array is a buffer set re-bound to the
 * CURRENT pooled membership (membership changes on every force evaluation).
 */
UINT TestMGOwnerRationalSet(CParameters&)
{
    UINT uiErrors = 0;
    CField* pGauge = appGetLattice()->m_pGaugeField[0];
    CField* pPooled = appGetLattice()->GetPooledFieldById(pGauge->m_byFieldId, _T(__FILE__), __LINE__);

    CRationalFieldPointer* pRFP = CRationalFieldPointer::GetInstance();
    deviceSU3** pArray = pRFP->GetRationPoint<deviceSU3>(2);

    TArray<const CHaloBufferHandle*> members;
    members.AddItem(pGauge->GetHaloBufferHandle());
    members.AddItem(pPooled->GetHaloBufferHandle());
    pRFP->RebindMemberHandles(members);

    const CHaloBufferSetHandle* pSet = appGetHaloManager()->LookupSetByPointer(pArray);
    if (NULL == pSet || 2 != pSet->MemberCount())
    {
        appGeneral(_T("rational set: missing / %d members != 2.\n"),
            (NULL == pSet) ? -1 : static_cast<INT>(pSet->MemberCount()));
        ++uiErrors;
    }
    //An interior element pointer resolves to the same set.
    if (appGetHaloManager()->LookupSetByPointer(reinterpret_cast<BYTE*>(pArray) + sizeof(void*)) != pSet)
    {
        appGeneral(_T("rational set: interior element lookup miss.\n"));
        ++uiErrors;
    }
    //Guard expansion records exactly the two member handles.
    CHaloLaunchGuard guard = appGetHaloManager()->BeginFromArguments(static_cast<deviceSU3* const*>(pArray));
    if (2 != guard.RecordedCount())
    {
        appGeneral(_T("rational set: guard recorded %d handles != 2 members.\n"), guard.RecordedCount());
        ++uiErrors;
    }
    guard.Commit();

    //Re-bind with empty membership detaches the set (membership-change rule).
    pRFP->RebindMemberHandles(TArray<const CHaloBufferHandle*>());
    if (NULL != appGetHaloManager()->LookupSetByPointer(pArray))
    {
        appGeneral(_T("rational set: empty re-bind did not detach.\n"));
        ++uiErrors;
    }

    pPooled->Return();
    return uiErrors;
}
___REGIST_TEST(TestMGOwnerRationalSet, MG, TestMGOwnerRationalSet, MGOwnerRationalSet, _TEST_MULTIGPU);

/**
 * I6 gate 4: the lazy U1Real Fmunu buffer is registered LocalOnly by its
 * owner when it materializes (the CopyBufferTo propagation path).
 */
UINT TestMGOwnerU1RealFmunu(CParameters&)
{
    UINT uiErrors = 0;
    CFieldGaugeU1Real* pSource = new CFieldGaugeU1Real();
    CFieldGaugeU1Real* pTarget = new CFieldGaugeU1Real();

    //Materialize the lazy cache on the source; CopyBufferTo then propagates it
    //(and the owner registration) to the target.
    checkCudaErrors(__cudaMalloc((void**)&pSource->m_pFmunu, sizeof(Real) * _HC_Volume * 6));
    pSource->m_bCacheFmunu = TRUE;
    pSource->CopyBufferTo(pTarget);

    if (!pTarget->m_FmunuHaloBuffer.IsBound())
    {
        appGeneral(_T("U1Real Fmunu: target handle not bound after propagation.\n"));
        ++uiErrors;
    }
    else
    {
        const SHaloBufferInfo& info = pTarget->m_FmunuHaloBuffer.Info();
        if (info.m_bHaloCapable || 0 != info.m_uiHaloSiteCount
            || info.m_uiBytesPerSite != 6 * static_cast<UINT>(sizeof(Real))
            || info.m_uiCapacityBytes != sizeof(Real) * _HC_Volume * 6)
        {
            appGeneral(_T("U1Real Fmunu: not a LocalOnly 6*V extent.\n"));
            ++uiErrors;
        }
        if (appGetHaloManager()->LookupByPointer(pTarget->m_pFmunu) != &pTarget->m_FmunuHaloBuffer)
        {
            appGeneral(_T("U1Real Fmunu: pointer lookup miss.\n"));
            ++uiErrors;
        }
    }

    appSafeDelete(pSource);
    appSafeDelete(pTarget);
    return uiErrors;
}
___REGIST_TEST(TestMGOwnerU1RealFmunu, MG, TestMGOwnerU1RealFmunu, MGOwnerU1RealFmunu, _TEST_MULTIGPU);
