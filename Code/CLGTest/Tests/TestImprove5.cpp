//=============================================================================
// FILENAME : TestImprove5.cpp
//
// DESCRIPTION:
// Improve-1 I5 (multi-GPU-improve1.md section 3.3/3.4.4): the automatic
// launch guard. Every public _LAUNCH_KERNEL* macro inspects its arguments
// exactly once before the backend launch (CHaloManager::BeginFromArguments):
// a pointer hitting a registered HaloCapable buffer Ensures its halo, every
// managed hit is recorded once, and the successful launch conservatively
// NotifyWritten-s every recorded handle (Commit). The rank-local variant
// never Ensures and fail-fasts on a HaloCapable hit.
//
// Covered gates (I5):
//   - registry range lookup: base AND interior pointers of one extent hit the
//     same handle and are de-duplicated; scalars / NULL / unmanaged pointers
//     stay ordinary arguments;
//   - ordering: the Ensure has already synced the halo to the pre-launch
//     interior version BEFORE Commit, and Commit only invalidates;
//   - Commit is idempotent;
//   - LocalOnly buffers are recorded WITHOUT an Ensure (no MPI);
//   - rank-local: LocalOnly recorded, HaloCapable hit rejected (fail-fast);
//   - a real launch through the public macro (InitialField) shows the whole
//     sequence on the MG build, and byte-for-byte historical behaviour
//     (no guard effects at all) on the single-GPU build.
//
// REVISION:
//  [08/09/2026 Improve-1 I5 nbale]
//=============================================================================

#include "CLGTest.h"

/**
 * I5 gate: guard traversal -- range lookup, de-dup, ordering, idempotence,
 * and the treatment of unmanaged / LocalOnly arguments.
 */
UINT TestMGLaunchGuardScan(CParameters&)
{
    UINT uiErrors = 0;
    CField* pGauge = appGetLattice()->m_pGaugeField[0];
    CHaloBufferHandle* pHandle = pGauge->GetHaloBufferHandle();
    if (NULL == pHandle || !pHandle->IsBound())
    {
        appGeneral(_T("guard scan: gauge handle missing.\n"));
        return 1;
    }
    const SHaloBufferInfo& info = pHandle->Info();
    BYTE* pData = const_cast<BYTE*>(info.m_pDeviceData);

    //Start from a dirty halo so the pre-launch Ensure is observable.
    pHandle->NotifyWritten();
    const ULONGLONG ullV0 = pHandle->InteriorVersion();

    //Two pointers into the SAME extent (base + one site in), a scalar, a NULL.
    CHaloLaunchGuard guard = appGetHaloManager()->BeginFromArguments(
        pData, pData + info.m_uiBytesPerSite, 42, static_cast<const BYTE*>(NULL));

    if (1 != guard.RecordedCount())
    {
        appGeneral(_T("guard scan: RecordedCount %d != 1 (dedup / range lookup wrong).\n"), guard.RecordedCount());
        ++uiErrors;
    }
    if (pHandle->InteriorVersion() != ullV0)
    {
        appGeneral(_T("guard scan: interior version changed BEFORE Commit.\n"));
        ++uiErrors;
    }
#if _CLG_MULTI_GPU
    if (_HC_HaloSiteCount() > 0)
    {
        //The pre-launch Ensure must already have synced the halo to the
        //CURRENT (pre-Commit) interior version.
        if (pHandle->HaloVersion() != ullV0)
        {
            appGeneral(_T("guard scan: halo not Ensured pre-launch (halo %llu != interior %llu).\n"),
                pHandle->HaloVersion(), ullV0);
            ++uiErrors;
        }
        if (static_cast<UINT>(_HC_HaloWidth) != pHandle->ValidWidth())
        {
            appGeneral(_T("guard scan: valid width %d != configured %d after pre-launch Ensure.\n"),
                pHandle->ValidWidth(), static_cast<UINT>(_HC_HaloWidth));
            ++uiErrors;
        }
    }
#endif

    //Commit invalidates exactly once.
    guard.Commit();
    if (pHandle->InteriorVersion() == ullV0)
    {
        appGeneral(_T("guard scan: Commit did not bump the interior version.\n"));
        ++uiErrors;
    }
    if (0 != pHandle->ValidWidth())
    {
        appGeneral(_T("guard scan: Commit left valid width %d.\n"), pHandle->ValidWidth());
        ++uiErrors;
    }
#if _CLG_MULTI_GPU
    if (_HC_HaloSiteCount() > 0 && pHandle->HaloVersion() != ullV0)
    {
        appGeneral(_T("guard scan: Commit disturbed the pre-launch halo sync (halo %llu, want %llu).\n"),
            pHandle->HaloVersion(), ullV0);
        ++uiErrors;
    }
#endif
    const ULONGLONG ullV1 = pHandle->InteriorVersion();
    guard.Commit();
    if (pHandle->InteriorVersion() != ullV1)
    {
        appGeneral(_T("guard scan: second Commit not a no-op.\n"));
        ++uiErrors;
    }

    //Unmanaged pointer + scalar: nothing recorded, nothing touched.
    BYTE* pFake = NULL;
    checkCudaErrors(__cudaMalloc((void**)&pFake, 4 * info.m_uiBytesPerSite));
    if (NULL == pFake)
    {
        appGeneral(_T("guard scan: fake allocation failed.\n"));
        return uiErrors + 1;
    }
    {
        CHaloLaunchGuard gFree = appGetHaloManager()->BeginFromArguments(pFake, 17);
        if (0 != gFree.RecordedCount())
        {
            appGeneral(_T("guard scan: unmanaged pointer recorded.\n"));
            ++uiErrors;
        }
        gFree.Commit();
        if (pHandle->InteriorVersion() != ullV1)
        {
            appGeneral(_T("guard scan: unmanaged Commit touched the gauge handle.\n"));
            ++uiErrors;
        }
    }

    //LocalOnly buffer: recorded WITHOUT an Ensure (no MPI, no width change).
    {
        SHaloBufferInfo local = info;
        local.m_pDeviceData = pFake;
        local.m_uiLocalSiteCount = 4;
        local.m_uiHaloSiteCount = 0;
        local.m_uiCapacityBytes = 4 * info.m_uiBytesPerSite;
        local.m_bHaloCapable = FALSE;
        CHaloBufferHandle hLocal;
        if (!hLocal.Bind(local))
        {
            appGeneral(_T("guard scan: LocalOnly fake Bind rejected.\n"));
            ++uiErrors;
        }
        else
        {
            const ULONGLONG ullL0 = hLocal.InteriorVersion();
            CHaloLaunchGuard gLocal = appGetHaloManager()->BeginFromArguments(pFake);
            if (1 != gLocal.RecordedCount())
            {
                appGeneral(_T("guard scan: LocalOnly hit not recorded.\n"));
                ++uiErrors;
            }
            if (hLocal.InteriorVersion() != ullL0 || 0 != hLocal.ValidWidth())
            {
                appGeneral(_T("guard scan: LocalOnly hit triggered an Ensure / version change.\n"));
                ++uiErrors;
            }
            gLocal.Commit();
            if (hLocal.InteriorVersion() == ullL0)
            {
                appGeneral(_T("guard scan: LocalOnly Commit did not invalidate.\n"));
                ++uiErrors;
            }
            hLocal.Unbind();
        }
    }
    checkCudaErrors(__cudaFree(pFake));
    return uiErrors;
}
___REGIST_TEST(TestMGLaunchGuardScan, MG, TestMGLaunchGuardScan, MGLaunchGuardScan, _TEST_MULTIGPU);

/**
 * I5 gate (3.4.4): the rank-local variant. A LocalOnly scratch buffer is
 * recorded for post-launch invalidation without any Ensure; a HaloCapable
 * hit is a caller bug -- the manager fail-fasts (a CRUCIAL line in the log is
 * EXPECTED from this test) and does NOT record it.
 */
UINT TestMGLaunchGuardRankLocal(CParameters&)
{
    UINT uiErrors = 0;
    CField* pGauge = appGetLattice()->m_pGaugeField[0];
    CHaloBufferHandle* pHandle = pGauge->GetHaloBufferHandle();
    if (NULL == pHandle || !pHandle->IsBound())
    {
        appGeneral(_T("rank-local: gauge handle missing.\n"));
        return 1;
    }
    const SHaloBufferInfo& info = pHandle->Info();

    BYTE* pFake = NULL;
    checkCudaErrors(__cudaMalloc((void**)&pFake, 4 * info.m_uiBytesPerSite));
    if (NULL == pFake)
    {
        appGeneral(_T("rank-local: fake allocation failed.\n"));
        return 1;
    }

    //LocalOnly scratch: recorded, never Ensured.
    SHaloBufferInfo local = info;
    local.m_pDeviceData = pFake;
    local.m_uiLocalSiteCount = 4;
    local.m_uiHaloSiteCount = 0;
    local.m_uiCapacityBytes = 4 * info.m_uiBytesPerSite;
    local.m_bHaloCapable = FALSE;
    CHaloBufferHandle hLocal;
    if (!hLocal.Bind(local))
    {
        appGeneral(_T("rank-local: LocalOnly fake Bind rejected.\n"));
        checkCudaErrors(__cudaFree(pFake));
        return uiErrors + 1;
    }
    const ULONGLONG ullL0 = hLocal.InteriorVersion();
    CHaloLaunchGuard gLocal = appGetHaloManager()->BeginFromArgumentsRankLocal(pFake);
    if (1 != gLocal.RecordedCount())
    {
        appGeneral(_T("rank-local: LocalOnly scratch not recorded.\n"));
        ++uiErrors;
    }
    if (hLocal.InteriorVersion() != ullL0 || 0 != hLocal.ValidWidth())
    {
        appGeneral(_T("rank-local: LocalOnly scratch was Ensured (must never happen).\n"));
        ++uiErrors;
    }
    gLocal.Commit();
    if (hLocal.InteriorVersion() == ullL0)
    {
        appGeneral(_T("rank-local: Commit did not invalidate the scratch.\n"));
        ++uiErrors;
    }
    hLocal.Unbind();
    checkCudaErrors(__cudaFree(pFake));

    //Negative: a HaloCapable argument is a caller bug. The manager logs a
    //CRUCIAL (expected!) and does not record the hit.
    const ULONGLONG ullV0 = pHandle->InteriorVersion();
    CHaloLaunchGuard gBad = appGetHaloManager()->BeginFromArgumentsRankLocal(
        const_cast<BYTE*>(info.m_pDeviceData));
    if (0 != gBad.RecordedCount())
    {
        appGeneral(_T("rank-local: HaloCapable hit was recorded (fail-fast broken).\n"));
        ++uiErrors;
    }
    gBad.Commit();
    if (pHandle->InteriorVersion() != ullV0)
    {
        appGeneral(_T("rank-local: rejected hit still invalidated the gauge handle.\n"));
        ++uiErrors;
    }
    return uiErrors;
}
___REGIST_TEST(TestMGLaunchGuardRankLocal, MG, TestMGLaunchGuardRankLocal, MGLaunchGuardRankLocal, _TEST_MULTIGPU);

/**
 * I5 gate: a REAL launch through the public macro. InitialField(EFIT_Zero)
 * launches _kernelInitialLink via _LAUNCH_KERNEL with the managed gauge data
 * pointer as an argument.
 *
 * MG build: the guard Ensures the halo to the pre-launch interior version
 * before the kernel runs and invalidates after it; the field's own explicit
 * NotifyWritten then bumps once more. Post-state proof: halo version == the
 * PRE-launch interior version (the Ensure ran), interior version != it (the
 * invalidations ran), width 0.
 *
 * Single-GPU build: the macro is the raw backend AND the field's explicit
 * NotifyWritten is _CLG_MULTI_GPU-only, so NOTHING may change at all --
 * byte-for-byte historical behaviour.
 */
UINT TestMGLaunchGuardMacro(CParameters&)
{
    UINT uiErrors = 0;
    CField* pGauge = appGetLattice()->m_pGaugeField[0];
    CHaloBufferHandle* pHandle = pGauge->GetHaloBufferHandle();
    if (NULL == pHandle || !pHandle->IsBound())
    {
        appGeneral(_T("guard macro: gauge handle missing.\n"));
        return 1;
    }

    //Dirty halo, known pre-launch versions.
    pHandle->NotifyWritten();
    const ULONGLONG ullV0 = pHandle->InteriorVersion();
    const ULONGLONG ullH0 = pHandle->HaloVersion();

    pGauge->InitialField(EFIT_Zero);

#if _CLG_MULTI_GPU
    if (pHandle->InteriorVersion() == ullV0)
    {
        appGeneral(_T("guard macro: launch did not invalidate the gauge handle.\n"));
        ++uiErrors;
    }
    if (0 != pHandle->ValidWidth())
    {
        appGeneral(_T("guard macro: launch left valid width %d.\n"), pHandle->ValidWidth());
        ++uiErrors;
    }
    if (_HC_HaloSiteCount() > 0 && pHandle->HaloVersion() != ullV0)
    {
        appGeneral(_T("guard macro: halo %llu != pre-launch interior %llu (pre-launch Ensure missing).\n"),
            pHandle->HaloVersion(), ullV0);
        ++uiErrors;
    }
#else
    //Raw backend, no MG write notification: no version may move at all.
    if (pHandle->InteriorVersion() != ullV0 || pHandle->HaloVersion() != ullH0)
    {
        appGeneral(_T("guard macro: single-GPU build moved a version (guard/notify leaked).\n"));
        ++uiErrors;
    }
#endif
    return uiErrors;
}
___REGIST_TEST(TestMGLaunchGuardMacro, MG, TestMGLaunchGuardMacro, MGLaunchGuardMacro, _TEST_MULTIGPU);

/**
 * I5 gate (3.4.1): the buffer-set handle. A device pointer array's outer
 * extent is bound to a fixed member set (the gauge field's handle plus a
 * pooled copy's handle -- two REAL HaloCapable handles); a launch argument
 * hitting the outer extent expands the set: both members Ensured pre-Commit,
 * both invalidated post-Commit. Also covers overlap rejection, double-Bind
 * rejection, dedup between a direct member hit and the set expansion, and
 * Unbind restoring the unmanaged state.
 */
UINT TestMGLaunchGuardSet(CParameters&)
{
    UINT uiErrors = 0;
    CField* pGauge = appGetLattice()->m_pGaugeField[0];
    CField* pCopy = appGetLattice()->GetPooledCopy(pGauge, __FILE__, __LINE__);
    if (NULL == pCopy)
    {
        appGeneral(_T("guard set: pool copy failed.\n"));
        return 1;
    }
    CHaloBufferHandle* pH1 = pGauge->GetHaloBufferHandle();
    CHaloBufferHandle* pH2 = pCopy->GetHaloBufferHandle();

    BYTE* pOuter = NULL;
    checkCudaErrors(__cudaMalloc((void**)&pOuter, 4 * sizeof(void*)));
    if (NULL == pOuter)
    {
        appGeneral(_T("guard set: outer allocation failed.\n"));
        pCopy->Return();
        return 1;
    }

    //Negative: outer extent overlapping a real handle extent must be rejected.
    {
        TArray<const CHaloBufferHandle*> members;
        members.AddItem(pH2);
        CHaloBufferSetHandle bad;
        if (bad.Bind(const_cast<BYTE*>(pH1->Info().m_pDeviceData), 2 * sizeof(void*), members))
        {
            appGeneral(_T("guard set: outer extent overlapping a handle accepted.\n"));
            bad.Unbind();
            ++uiErrors;
        }
    }

    TArray<const CHaloBufferHandle*> members;
    members.AddItem(pH1);
    members.AddItem(pH2);
    CHaloBufferSetHandle set;
    if (!set.Bind(pOuter, 4 * sizeof(void*), members))
    {
        appGeneral(_T("guard set: Bind rejected.\n"));
        ++uiErrors;
    }
    else
    {
        //Double Bind rejected.
        if (set.Bind(pOuter, 4 * sizeof(void*), members))
        {
            appGeneral(_T("guard set: double Bind accepted.\n"));
            ++uiErrors;
        }

        //Dirty both members, then hit the set with TWO outer pointers (the
        //second re-expansion must dedup to the same two members).
        pH1->NotifyWritten();
        pH2->NotifyWritten();
        const ULONGLONG ullV1 = pH1->InteriorVersion();
        const ULONGLONG ullV2 = pH2->InteriorVersion();

        CHaloLaunchGuard guard = appGetHaloManager()->BeginFromArguments(
            pOuter, pOuter + sizeof(void*));
        if (2 != guard.RecordedCount())
        {
            appGeneral(_T("guard set: RecordedCount %d != 2 after expansion.\n"), guard.RecordedCount());
            ++uiErrors;
        }
        if (pH1->InteriorVersion() != ullV1 || pH2->InteriorVersion() != ullV2)
        {
            appGeneral(_T("guard set: a member version changed BEFORE Commit.\n"));
            ++uiErrors;
        }
#if _CLG_MULTI_GPU
        if (_HC_HaloSiteCount() > 0
            && (pH1->HaloVersion() != ullV1 || pH2->HaloVersion() != ullV2))
        {
            appGeneral(_T("guard set: a member halo not Ensured pre-launch.\n"));
            ++uiErrors;
        }
#endif
        guard.Commit();
        if (pH1->InteriorVersion() == ullV1 || pH2->InteriorVersion() == ullV2)
        {
            appGeneral(_T("guard set: Commit did not invalidate both members.\n"));
            ++uiErrors;
        }

        //A direct member pointer plus the outer extent in one launch: the
        //member must still be processed only once.
        CHaloLaunchGuard g2 = appGetHaloManager()->BeginFromArguments(
            pOuter, const_cast<BYTE*>(pH1->Info().m_pDeviceData));
        if (2 != g2.RecordedCount())
        {
            appGeneral(_T("guard set: direct+set dedup gave RecordedCount %d != 2.\n"), g2.RecordedCount());
            ++uiErrors;
        }
        g2.Commit();

        set.Unbind();
        CHaloLaunchGuard g3 = appGetHaloManager()->BeginFromArguments(pOuter);
        if (0 != g3.RecordedCount())
        {
            appGeneral(_T("guard set: outer pointer still managed after Unbind.\n"));
            ++uiErrors;
        }
        g3.Commit();
    }
    checkCudaErrors(__cudaFree(pOuter));
    pCopy->Return();
    return uiErrors;
}
___REGIST_TEST(TestMGLaunchGuardSet, MG, TestMGLaunchGuardSet, MGLaunchGuardSet, _TEST_MULTIGPU);

/**
 * I5 gate (3.3): single evaluation. Drives appLaunchGuardSingleEvalProbe,
 * which launches a probe kernel through the PUBLIC macro with side-effecting
 * block / thread / argument expressions. Every counter must be exactly 1 and
 * the device marker must read 7 + 8. On the raw (<<<>>>) build the same test
 * binary exercises the lambda path -- run manually there per the I5 gate.
 * Single-GPU build: the probe is MG-only; the historical macro evaluates
 * each expression once by construction, so the test trivially passes.
 */
UINT TestMGLaunchGuardSingleEval(CParameters&)
{
    UINT uiErrors = 0;
#if _CLG_MULTI_GPU
    UINT uiCounts[4] = { 0, 0, 0, 0 };
    UINT* pMarker = NULL;
    checkCudaErrors(__cudaMalloc((void**)&pMarker, sizeof(UINT)));
    checkCudaErrors(cudaMemset(pMarker, 0, sizeof(UINT)));

    appLaunchGuardSingleEvalProbe(uiCounts, pMarker);

    for (UINT i = 0; i < 4; ++i)
    {
        if (1 != uiCounts[i])
        {
            appGeneral(_T("single-eval: expression %d evaluated %d times (want 1).\n"), i, uiCounts[i]);
            ++uiErrors;
        }
    }
    UINT uiMarker = 0;
    checkCudaErrors(cudaMemcpy(&uiMarker, pMarker, sizeof(UINT), cudaMemcpyDeviceToHost));
    if (15 != uiMarker)
    {
        appGeneral(_T("single-eval: device marker %d != 15 (kernel did not run with 7 and 8).\n"), uiMarker);
        ++uiErrors;
    }
    checkCudaErrors(__cudaFree(pMarker));
#else
    appGeneral(_T("single-eval: single-GPU build, probe is MG-only (historical macro is single-eval by construction).\n"));
#endif
    return uiErrors;
}
___REGIST_TEST(TestMGLaunchGuardSingleEval, MG, TestMGLaunchGuardSingleEval, MGLaunchGuardSingleEval, _TEST_MULTIGPU);

//=============================================================================
// END OF FILE
//=============================================================================
