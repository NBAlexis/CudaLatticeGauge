//=============================================================================
// FILENAME : TestImprove2.cpp
//
// DESCRIPTION:
// Improve-1 I2 (multi-GPU-improve1.md section 3.1/3.2): CHaloBufferHandle and
// the explicit per-field halo descriptor. These tests belong to the normal MG
// suite (they must PASS once I2 lands) -- unlike the TestImprove1.cpp group
// which reproduces still-open defects.
//
// Covered gates (I2):
//   - fields expose a bound handle whose pointer/stride/capacity match the
//     real allocation (metadata, not guessed);
//   - same field id but a different buffer extent (pool copy) -> different
//     handle;
//   - version semantics: NotifyWritten bumps interior and drops valid width,
//     NotifyHaloSynchronized sets halo==interior and records the width;
//   - the manager registry rejects overlapping/duplicate extents and
//     capacity-short Binds; a freed (Unbound) extent can be re-Bound.
//
// REVISION:
//  [08/09/2026 Improve-1 I2 nbale]
//=============================================================================

#include "CLGTest.h"

/**
 * Validate one field's bound handle against its real allocation. Returns the
 * number of failed checks (prints each mismatch).
 */
static UINT Improve2CheckFieldHandle(const CField* pField, const TCHAR* sWhat, UBOOL bGauge)
{
    UINT uiErrors = 0;
    if (NULL == pField)
    {
        appGeneral(_T("%s: field missing, skipped.\n"), sWhat);
        return 0;
    }
    const CHaloBufferHandle* pHandle = pField->GetHaloBufferHandle();
    if (NULL == pHandle)
    {
        appGeneral(_T("%s: GetHaloBufferHandle() is NULL, field not wired.\n"), sWhat);
        return 1;
    }
    if (!pHandle->IsBound())
    {
        appGeneral(_T("%s: handle not bound.\n"), sWhat);
        return 1;
    }
    const SHaloBufferInfo& info = pHandle->Info();

    if (info.m_pDeviceData != static_cast<const BYTE*>(pField->GetData()))
    {
        appGeneral(_T("%s: info pointer %p != field GetData %p.\n"), sWhat, info.m_pDeviceData, pField->GetData());
        ++uiErrors;
    }
    if (0 == info.m_uiBytesPerSite)
    {
        appGeneral(_T("%s: bytesPerSite is 0.\n"), sWhat);
        ++uiErrors;
    }
    if (info.m_uiLocalSiteCount != static_cast<UINT>(_HC_Volume))
    {
        appGeneral(_T("%s: localSiteCount %d != volume %d.\n"), sWhat, info.m_uiLocalSiteCount, _HC_Volume);
        ++uiErrors;
    }
    if (info.m_uiHaloSiteCount != _HC_HaloSiteCount())
    {
        appGeneral(_T("%s: haloSiteCount %d != _HC_HaloSiteCount %d.\n"), sWhat, info.m_uiHaloSiteCount, _HC_HaloSiteCount());
        ++uiErrors;
    }
    const size_t uiNeed = (static_cast<size_t>(info.m_uiLocalSiteCount) + info.m_uiHaloSiteCount) * info.m_uiBytesPerSite;
    if (info.m_uiCapacityBytes != uiNeed)
    {
        appGeneral(_T("%s: capacity %llu != (local+halo)*bytesPerSite %llu.\n"),
            sWhat, static_cast<ULONGLONG>(info.m_uiCapacityBytes), static_cast<ULONGLONG>(uiNeed));
        ++uiErrors;
    }
    if (bGauge)
    {
        //Link field: Dir matrices of MatrixN x MatrixN complex per site.
        const CFieldGauge* pGauge = dynamic_cast<const CFieldGauge*>(pField);
        const UINT uiExpect = _HC_Dir * pGauge->MatrixN() * pGauge->MatrixN() * 2 * sizeof(Real);
        if (info.m_uiBytesPerSite != uiExpect)
        {
            appGeneral(_T("%s: gauge bytesPerSite %d != Dir*N*N*2*Real %d.\n"), sWhat, info.m_uiBytesPerSite, uiExpect);
            ++uiErrors;
        }
    }
    if (0 == info.m_ullLayoutGeneration || info.m_ullLayoutGeneration != appGetLayoutGeneration())
    {
        appGeneral(_T("%s: layout generation %llu stale (current %llu).\n"),
            sWhat, info.m_ullLayoutGeneration, appGetLayoutGeneration());
        ++uiErrors;
    }
    if (0 == pHandle->InteriorVersion())
    {
        appGeneral(_T("%s: interior version is 0 after Bind.\n"), sWhat);
        ++uiErrors;
    }
    return uiErrors;
}

/**
 * I2 gate: every wired field family exposes correct pointer/stride/capacity.
 * Single GPU: halo counts are 0 and the capacity must equal the pre-halo
 * allocation (behaviour unchanged).
 */
UINT TestMGHaloHandleMeta(CParameters&)
{
    UINT uiErrors = 0;
    uiErrors += Improve2CheckFieldHandle(appGetLattice()->m_pGaugeField[0], _T("gauge0"), TRUE);
    if (appGetLattice()->m_pFermionField.Num() > 0)
    {
        uiErrors += Improve2CheckFieldHandle(appGetLattice()->m_pFermionField[0], _T("fermion0"), FALSE);
    }
    if (appGetLattice()->m_pBosonField.Num() > 0)
    {
        uiErrors += Improve2CheckFieldHandle(appGetLattice()->m_pBosonField[0], _T("boson0"), FALSE);
    }
    return uiErrors;
}
___REGIST_TEST(TestMGHaloHandleMeta, MG, TestMGHaloHandleMeta, MGHaloHandleMeta, _TEST_MULTIGPU);

/**
 * I2 gate: a pool copy shares the field ID but owns a different buffer
 * extent, so it must carry a DIFFERENT bound handle.
 */
UINT TestMGHaloHandlePool(CParameters&)
{
    UINT uiErrors = 0;
    CField* pGauge = appGetLattice()->m_pGaugeField[0];
    CField* pCopy = appGetLattice()->GetPooledCopy(pGauge, __FILE__, __LINE__);
    if (NULL == pCopy)
    {
        appGeneral(_T("pool copy failed.\n"));
        return 1;
    }

    const CHaloBufferHandle* pH1 = pGauge->GetHaloBufferHandle();
    const CHaloBufferHandle* pH2 = pCopy->GetHaloBufferHandle();
    if (NULL == pH2 || !pH2->IsBound())
    {
        appGeneral(_T("pool copy: handle missing or unbound.\n"));
        ++uiErrors;
    }
    else
    {
        if (pH1 == pH2)
        {
            appGeneral(_T("pool copy: shares the OWNER's handle (version state would clash).\n"));
            ++uiErrors;
        }
        if (pH2->Info().m_pDeviceData == pH1->Info().m_pDeviceData)
        {
            appGeneral(_T("pool copy: same device pointer as owner.\n"));
            ++uiErrors;
        }
        if (pH2->Info().m_byFieldId != pH1->Info().m_byFieldId)
        {
            appGeneral(_T("pool copy: field id not propagated (%d vs %d).\n"),
                pH2->Info().m_byFieldId, pH1->Info().m_byFieldId);
            ++uiErrors;
        }
        if (pH2->Info().m_uiCapacityBytes != pH1->Info().m_uiCapacityBytes
            || pH2->Info().m_uiBytesPerSite != pH1->Info().m_uiBytesPerSite)
        {
            appGeneral(_T("pool copy: extent metadata differs from owner.\n"));
            ++uiErrors;
        }
    }
    pCopy->Return();
    return uiErrors;
}
___REGIST_TEST(TestMGHaloHandlePool, MG, TestMGHaloHandlePool, MGHaloHandlePool, _TEST_MULTIGPU);

/**
 * I2 version semantics: Bind -> interior non-zero, halo invalid;
 * NotifyWritten -> interior changes, width 0; NotifyHaloSynchronized ->
 * halo == interior, width recorded. Drives the handle directly; leaves the
 * field with a dirty halo state afterwards (NotifyWritten last).
 */
UINT TestMGHaloHandleVersion(CParameters&)
{
    UINT uiErrors = 0;
    CHaloBufferHandle* pHandle = appGetLattice()->m_pGaugeField[0]->GetHaloBufferHandle();
    if (NULL == pHandle || !pHandle->IsBound())
    {
        appGeneral(_T("version test: gauge handle missing.\n"));
        return 1;
    }

    const ULONGLONG ullV0 = pHandle->InteriorVersion();
    pHandle->NotifyWritten();
    if (pHandle->InteriorVersion() == ullV0)
    {
        appGeneral(_T("version test: NotifyWritten did not change interior version.\n"));
        ++uiErrors;
    }
    if (0 != pHandle->ValidWidth())
    {
        appGeneral(_T("version test: NotifyWritten left valid width %d.\n"), pHandle->ValidWidth());
        ++uiErrors;
    }

    pHandle->NotifyHaloSynchronized(1);
    if (pHandle->HaloVersion() != pHandle->InteriorVersion())
    {
        appGeneral(_T("version test: halo version != interior after sync.\n"));
        ++uiErrors;
    }
    if (1 != pHandle->ValidWidth())
    {
        appGeneral(_T("version test: valid width %d != 1 after sync.\n"), pHandle->ValidWidth());
        ++uiErrors;
    }

    //A second sync without an intermediate write must NOT change any version.
    const ULONGLONG ullV1 = pHandle->InteriorVersion();
    const ULONGLONG ullH1 = pHandle->HaloVersion();
    pHandle->NotifyHaloSynchronized(1);
    if (pHandle->InteriorVersion() != ullV1 || pHandle->HaloVersion() != ullH1)
    {
        appGeneral(_T("version test: repeated sync changed a version.\n"));
        ++uiErrors;
    }

    //Leave the halo stale (honest state after a version-test drive).
    pHandle->NotifyWritten();
    return uiErrors;
}
___REGIST_TEST(TestMGHaloHandleVersion, MG, TestMGHaloHandleVersion, MGHaloHandleVersion, _TEST_MULTIGPU);

/**
 * I2 registry: the manager rejects a second handle over the same or an
 * overlapping extent (exact duplicate and sub-range of a live field buffer),
 * rejects a capacity-short Bind, and accepts a disjoint extent; Unbind makes
 * the extent available again.
 */
UINT TestMGHaloHandleRegistry(CParameters&)
{
    UINT uiErrors = 0;
    const CHaloBufferHandle* pField = appGetLattice()->m_pGaugeField[0]->GetHaloBufferHandle();
    if (NULL == pField || !pField->IsBound())
    {
        appGeneral(_T("registry test: gauge handle missing.\n"));
        return 1;
    }
    const SHaloBufferInfo& info = pField->Info();

    //Exact duplicate extent -> must be rejected.
    {
        CHaloBufferHandle dup;
        if (dup.Bind(info))
        {
            appGeneral(_T("registry test: duplicate extent accepted.\n"));
            dup.Unbind();
            ++uiErrors;
        }
    }

    //Sub-range INSIDE the live extent -> overlap, must be rejected. Counts and
    //capacity stay layout-consistent (OnHandleBound validates them now), so the
    //rejection comes from the overlap check, not the count check.
    {
        SHaloBufferInfo sub = info;
        sub.m_pDeviceData = info.m_pDeviceData + info.m_uiBytesPerSite;
        CHaloBufferHandle hSub;
        if (hSub.Bind(sub))
        {
            appGeneral(_T("registry test: overlapping sub-range accepted.\n"));
            hSub.Unbind();
            ++uiErrors;
        }
    }

    //Disjoint fresh extent -> accepted; duplicate over it -> rejected; after
    //Unbind the same extent binds again (new generation, no inherited state).
    //Counts must stay consistent with the current layout (OnHandleBound rejects
    //mismatched local/halo counts), so the probe uses full-field counts even
    //though it only exercises registry semantics.
    const UINT uiProbeLocal = static_cast<UINT>(_HC_Volume);
    const UINT uiProbeHalo = static_cast<UINT>(_HC_HaloSiteCount());
    const UINT uiProbeBytes = (uiProbeLocal + uiProbeHalo) * info.m_uiBytesPerSite;
    BYTE* pFake = NULL;
    checkCudaErrors(__cudaMalloc((void**)&pFake, uiProbeBytes));
    if (NULL == pFake)
    {
        appGeneral(_T("registry test: fake allocation failed.\n"));
        return uiErrors + 1;
    }
    {
        SHaloBufferInfo fake = info;
        fake.m_pDeviceData = pFake;
        fake.m_uiLocalSiteCount = uiProbeLocal;
        fake.m_uiHaloSiteCount = uiProbeHalo;
        fake.m_uiCapacityBytes = uiProbeBytes;

        //Capacity one byte short -> rejected before any registry decision.
        SHaloBufferInfo shortCap = fake;
        shortCap.m_uiCapacityBytes = fake.m_uiCapacityBytes - 1;
        CHaloBufferHandle hShort;
        if (hShort.Bind(shortCap))
        {
            appGeneral(_T("registry test: capacity-short Bind accepted.\n"));
            hShort.Unbind();
            ++uiErrors;
        }

        CHaloBufferHandle h1;
        if (!h1.Bind(fake))
        {
            appGeneral(_T("registry test: disjoint extent rejected.\n"));
            ++uiErrors;
        }
        CHaloBufferHandle h2;
        if (h2.Bind(fake))
        {
            appGeneral(_T("registry test: second handle over same disjoint extent accepted.\n"));
            h2.Unbind();
            ++uiErrors;
        }
        h1.Unbind();
        if (!h2.Bind(fake))
        {
            appGeneral(_T("registry test: re-Bind after Unbind rejected.\n"));
            ++uiErrors;
        }
        if (h2.IsBound() && 0 == h2.InteriorVersion())
        {
            appGeneral(_T("registry test: re-Bind inherited a zero/old interior version.\n"));
            ++uiErrors;
        }
        h2.Unbind();
    }
    checkCudaErrors(__cudaFree(pFake));
    return uiErrors;
}
___REGIST_TEST(TestMGHaloHandleRegistry, MG, TestMGHaloHandleRegistry, MGHaloHandleRegistry, _TEST_MULTIGPU);

//=============================================================================
// END OF FILE
//=============================================================================
