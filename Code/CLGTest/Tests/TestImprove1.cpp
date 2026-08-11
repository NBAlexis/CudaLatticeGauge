//=============================================================================
// FILENAME : TestImprove1.cpp
//
// DESCRIPTION:
// Improve-1 (multi-GPU-improve1.md) halo defect reproduction tests, originally
// the opt-in category "Improve1KnownFailures". Every defect they exposed is
// fixed as of I4 (generic handle refill I3; handle-version Ensure + uniform
// host-write NotifyWritten I4), so the group was dissolved: all five tests
// below now live in the normal MG suite and MUST stay green.
//
// All halo-content checks use a deterministic global-site-index pattern
// (value(g) = g + bias, exactly representable in both precisions), never RNG:
// every halo slot is checked against the GLOBAL site it represents, per
// element. The slot -> global-site decode below mirrors TestMGHaloSelfExchange
// (FACE/EDGE/CORNER per CLGHaloLayout.h); I8 replaces it with the baked
// 32-bit global-coordinate table (which also covers codimension-4).
//
// Pre-fix failure commands (Release_MG, from Bin/Ubuntu, historical record;
// current driver: `runtest.sh <TestName> 1` from the repo root resolves the
// rank count and grid from the test metadata):
//   mpiexec -n 2 ./CLGTest TestMGHaloSU3_12 --mg-worker --gpu-grid 1,1,1,2 --device-per-node 1
//   mpiexec -n 2 ./CLGTest TestMGHaloU1Real --mg-worker --gpu-grid 1,1,1,2 --device-per-node 1
//   mpiexec -n 2 ./CLGTest TestMGHaloBoson --mg-worker --gpu-grid 1,1,1,2 --device-per-node 1
//   mpiexec -n 2 ./CLGTest TestMGHaloBufferAlias --mg-worker --gpu-grid 1,1,1,2 --device-per-node 1
//   mpiexec -n 2 ./CLGTest TestMGHaloWriteVersion --mg-worker --gpu-grid 1,1,1,2 --device-per-node 1
//   mpiexec -n 2 ./CLGTest TestMGBosonConsistency --mg-worker --gpu-grid 1,1,1,2 --device-per-node 1
// (-n 1 keeps the trivial pass: halo count 0 / nothing split.)
//
// REVISION:
//  [08/08/2026 Improve-1 I1 nbale]
//=============================================================================

#include "CLGTest.h"

#include "TestImproveHaloCommon.h"

/**
 * I1-2a: compact gauge (SU3_12, 12 reals/link) halo capacity + stride.
 *
 * Defect exposed: CFieldGaugeSU3_12 allocates LOCAL links only
 * (m_uiLinkeCount * sizeof(deviceSU3_12), no halo tail), and the manager's
 * gauge stride guess (_HC_Dir * MatrixN^2 * sizeof(CLGComplex) = 72 reals/site
 * for N=3) does not match the real storage (_HC_Dir * 12 = 48 reals/site).
 * The capacity assertion fails first and is deterministic; the functional
 * refill check runs only when capacity exists (avoids an out-of-bounds write
 * against the current allocation) and then proves the stride via halo bytes.
 */
UINT TestMGHaloSU3_12(CParameters&)
{
    UINT uiError = 0;
    CFieldGaugeSU3_12* pField = new CFieldGaugeSU3_12();

    const UINT uiElemPerSite = _HC_Dir * 12; //deviceSU3_12 = 6 CLGComplex
    const UINT uiBytesPerSite = uiElemPerSite * static_cast<UINT>(sizeof(Real));

    //1. Capacity: a halo-capable link field must own halo link slots.
    const UINT uiNeedHaloLinks = _HC_HaloLinkCount();
    if (pField->GetHaloLinkCount() < uiNeedHaloLinks)
    {
        ++uiError;
        CCString sProblem;
        sProblem.Format(_T("SU3_12 halo link capacity %u < required %u: halo allocation missing"),
            pField->GetHaloLinkCount(), uiNeedHaloLinks);
        appGeneral(_T("%s\n"), sProblem.c_str());
        LastProbem(sProblem);
    }

    //2. Functional: deterministic fill + refill + per-slot content check.
    SImprove1HaloCtx ctx;
    if (0 == uiError && Improve1BuildCtx(ctx))
    {
        Improve1FillLocalDeterministic(pField->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(0.0));
        Improve1TriggerRefill(pField);
        uiError += Improve1CheckHalo(pField->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(0.0), _T("SU3_12"));
    }
    else if (0 == uiError)
    {
        appGeneral(_T("Halo count is 0 (nothing split), functional part trivially passes.\n"));
    }

    appSafeDelete(pField);
    return uiError;
}
___REGIST_TEST(TestMGHaloSU3_12, MG, TestMGHaloSU3_12, HaloSU3_12, _TEST_MULTIGPU);

/**
 * I1-2b: background gauge (U1Real, 1 Real/link) halo capacity + stride.
 *
 * Defect exposed: CFieldGaugeU1Real allocates LOCAL links only, yet its
 * CalculateForceAndStaple / CalculatePlaqutteEnergy read the per-field staple
 * and plaquette caches, which redirect split-direction neighbours to halo
 * slots under a decomposed grid -- out-of-bounds on the current allocation.
 * The manager's gauge stride guess (Dir*1*1*sizeof(CLGComplex) = 8 reals/site)
 * is also wrong for the real Dir*sizeof(Real) = 4 reals/site storage.
 */
UINT TestMGHaloU1Real(CParameters&)
{
    UINT uiError = 0;
    CFieldGaugeU1Real* pField = new CFieldGaugeU1Real();

    const UINT uiElemPerSite = _HC_Dir; //1 Real per link
    const UINT uiBytesPerSite = uiElemPerSite * static_cast<UINT>(sizeof(Real));

    const UINT uiNeedHaloLinks = _HC_HaloLinkCount();
    if (pField->GetHaloLinkCount() < uiNeedHaloLinks)
    {
        ++uiError;
        CCString sProblem;
        sProblem.Format(_T("U1Real halo link capacity %u < required %u: halo allocation missing"),
            pField->GetHaloLinkCount(), uiNeedHaloLinks);
        appGeneral(_T("%s\n"), sProblem.c_str());
        LastProbem(sProblem);
    }

    SImprove1HaloCtx ctx;
    if (0 == uiError && Improve1BuildCtx(ctx))
    {
        Improve1FillLocalDeterministic(pField->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(0.0));
        Improve1TriggerRefill(pField);
        uiError += Improve1CheckHalo(pField->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(0.0), _T("U1Real"));
    }
    else if (0 == uiError)
    {
        appGeneral(_T("Halo count is 0 (nothing split), functional part trivially passes.\n"));
    }

    appSafeDelete(pField);
    return uiError;
}
___REGIST_TEST(TestMGHaloU1Real, MG, TestMGHaloU1Real, HaloU1Real, _TEST_MULTIGPU);

/**
 * I1-1 companion: boson generic-refill gap (deterministic, in-process).
 *
 * Defect exposed: CHaloManager::RefillHalo resolves the field by id and only
 * dynamic_casts to CFieldGauge / CFieldFermion; a boson id silently returns,
 * and Ensure() then records a valid width for a halo that was NEVER filled.
 * The boson kernels (CFieldBosonVNKernel::DFromSource / ForceOnGauge) call
 * exactly this Ensure+RefillHalo pair before their stencil reads, so on a
 * split grid the D operator reads the halo tail uninitialised.
 *
 * The check uses CFieldBosonSU3 (deviceSU3Vector = 3 complex = 6 reals/site,
 * matching the yaml block). End-to-end 1-vs-N coverage: TestMGBosonConsistency.
 */
UINT TestMGHaloBoson(CParameters&)
{
    UINT uiError = 0;
    if (appGetLattice()->m_pBosonField.Num() <= 0)
    {
        LastProbem(_T("Need a boson field (yaml BosonFieldCount 1)"));
        return 1;
    }
    CFieldBosonSU3* pBoson = dynamic_cast<CFieldBosonSU3*>(appGetLattice()->m_pBosonField[0]);
    if (NULL == pBoson)
    {
        LastProbem(_T("TestMGHaloBoson requires CFieldBosonSU3 (6 reals/site)"));
        return 1;
    }

    const UINT uiElemPerSite = 6; //deviceSU3Vector = 3 CLGComplex
    const UINT uiBytesPerSite = uiElemPerSite * static_cast<UINT>(sizeof(Real));

    SImprove1HaloCtx ctx;
    if (!Improve1BuildCtx(ctx))
    {
        appGeneral(_T("Halo count is 0 (nothing split), test trivially passes.\n"));
        return 0;
    }

    Improve1FillLocalDeterministic(pBoson->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(0.0));

    //Exactly what the production boson stencil does before reading neighbours.
    appGetHaloManager()->Ensure(*pBoson->GetHaloBufferHandle(), ctx.uiHaloWidth);
    appGetHaloManager()->RefillHalo(*pBoson->GetHaloBufferHandle());

    uiError += Improve1CheckHalo(pBoson->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(0.0), _T("Boson"));
    if (0 == uiError)
    {
        appGeneral(_T("Boson halo verified: all %u slots hold the expected global site data.\n"), ctx.uiHaloSites);
    }
    return uiError;
}
___REGIST_TEST(TestMGHaloBoson, MG, TestMGHaloBoson, HaloBoson, _TEST_MULTIGPU);

/**
 * I1-3: same field id, different buffer -- the per-id halo state cross-talk.
 *
 * The registered gauge field A (id 1) and a GetCopy() B (shares id 1, own
 * buffer) are filled with DIFFERENT deterministic patterns. Ensure(id) refills
 * only A and records m_uiValidWidth[1]; a later Ensure for the buffer a
 * stencil actually reads (B -- the HMC UPrime/pool-copy pattern) is a no-op,
 * so B's halo is never refilled. Post-fix (handle identity, 3.1) A and B have
 * independent states and B's Ensure refills B.
 */
UINT TestMGHaloBufferAlias(CParameters&)
{
    UINT uiError = 0;
    CFieldGauge* pA = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField[0]);
    if (NULL == pA)
    {
        LastProbem(_T("No gauge field"));
        return 1;
    }
    CFieldGauge* pB = dynamic_cast<CFieldGauge*>(pA->GetCopy());
    if (NULL == pB)
    {
        LastProbem(_T("GetCopy did not return a gauge field"));
        return 1;
    }

    const UINT uiMatrixN = pA->MatrixN();
    const UINT uiElemPerSite = _HC_Dir * uiMatrixN * uiMatrixN * 2;
    const UINT uiBytesPerSite = uiElemPerSite * static_cast<UINT>(sizeof(Real));

    SImprove1HaloCtx ctx;
    if (!Improve1BuildCtx(ctx))
    {
        appGeneral(_T("Halo count is 0 (nothing split), test trivially passes.\n"));
        appSafeDelete(pB);
        return 0;
    }

    Improve1FillLocalDeterministic(pA->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(0.0));
    Improve1FillLocalDeterministic(pB->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(500000.0));

    //Ensure the REGISTERED field: refills A's halo and marks its handle valid.
    appGetHaloManager()->Ensure(*pA->GetHaloBufferHandle(), ctx.uiHaloWidth);
    uiError += Improve1CheckHalo(pA->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(0.0), _T("Alias-A"));
    if (0 != uiError)
    {
        //A itself is broken: this is a plain halo defect, not the alias state
        //bug this test targets; report and stop.
        LastProbem(_T("Alias sanity check on the registered field failed"));
        appSafeDelete(pB);
        return uiError;
    }

    //The stencil is about to read the COPY B. I4: Ensure B through its OWN
    //handle (3.6) -- independent of A's state, so B's halo is refilled from
    //B's interior even though A was just Ensured under the same field id.
    appGetHaloManager()->Ensure(*pB->GetHaloBufferHandle(), ctx.uiHaloWidth);
    uiError += Improve1CheckHalo(pB->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(500000.0), _T("Alias-B"));
    if (0 == uiError)
    {
        appGeneral(_T("Alias verified: two buffers sharing a field id hold independent halos.\n"));
    }
    appSafeDelete(pB);
    return uiError;
}
___REGIST_TEST(TestMGHaloBufferAlias, MG, TestMGHaloBufferAlias, HaloBufferAlias, _TEST_MULTIGPU);

/**
 * I1-4: Ensure -> BLAS write -> stencil read; the stale-halo gap.
 *
 * After Ensure() the gauge interior is overwritten by a BLAS call (AxpyPlus,
 * a kernel launch that does NOT MarkDirty -- multi-GPU-improve1.md problem 4).
 * The next Ensure() must refill; today it early-returns because the per-id
 * valid width is still recorded, so the stencil reads the pre-write halo.
 * AxpyPlus covers local links only (m_uiLinkeCount), so it cannot touch the
 * halo tail itself.
 */
UINT TestMGHaloWriteVersion(CParameters&)
{
    UINT uiError = 0;
    CFieldGauge* pA = dynamic_cast<CFieldGauge*>(appGetLattice()->m_pGaugeField[0]);
    if (NULL == pA)
    {
        LastProbem(_T("No gauge field"));
        return 1;
    }
    CFieldGauge* pAdd = dynamic_cast<CFieldGauge*>(pA->GetCopy());
    if (NULL == pAdd)
    {
        LastProbem(_T("GetCopy did not return a gauge field"));
        return 1;
    }

    const UINT uiMatrixN = pA->MatrixN();
    const UINT uiElemPerSite = _HC_Dir * uiMatrixN * uiMatrixN * 2;
    const UINT uiBytesPerSite = uiElemPerSite * static_cast<UINT>(sizeof(Real));

    SImprove1HaloCtx ctx;
    if (!Improve1BuildCtx(ctx))
    {
        appGeneral(_T("Halo count is 0 (nothing split), test trivially passes.\n"));
        appSafeDelete(pAdd);
        return 0;
    }

    //V1 interior, halo refilled (valid), sanity check.
    Improve1FillLocalDeterministic(pA->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(0.0));
    appGetHaloManager()->Ensure(*pA->GetHaloBufferHandle(), ctx.uiHaloWidth);
    uiError += Improve1CheckHalo(pA->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(0.0), _T("WriteVersion-V1"));
    if (0 != uiError)
    {
        LastProbem(_T("WriteVersion sanity check failed"));
        appSafeDelete(pAdd);
        return uiError;
    }

    //BLAS write of the interior: A = A + 1 (elementwise, local links only).
    Improve1FillLocalConstant(pAdd->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(1.0));
    pA->AxpyPlus(pAdd);

    //The next stencil's Ensure must see the write and refill. The expected
    //halo content is now (globalSite + 1).
    appGetHaloManager()->Ensure(*pA->GetHaloBufferHandle(), ctx.uiHaloWidth);
    uiError += Improve1CheckHalo(pA->GetData(), uiBytesPerSite, uiElemPerSite, ctx, F(1.0), _T("WriteVersion-V2"));
    if (0 == uiError)
    {
        appGeneral(_T("Write version verified: BLAS write invalidated the halo.\n"));
    }
    appSafeDelete(pAdd);
    return uiError;
}
___REGIST_TEST(TestMGHaloWriteVersion, MG, TestMGHaloWriteVersion, HaloWriteVersion, _TEST_MULTIGPU);

//=============================================================================
// END OF FILE
//=============================================================================
